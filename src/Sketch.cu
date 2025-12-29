
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "JEM.h"
#include "timers.h"
#include <fstream>
#include <chrono>
#include <iostream>

#include <cstdlib>

#include <unordered_map>
#include <unordered_set>


#include <numeric>
#include <cstdint>



#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

using std::ifstream; using std::ostringstream;

extern int rank, size;
extern int num_threads;
extern int no_trials;
extern int read_length;
extern int w_size;
extern std::string primeFileName;
extern std::string AFileName;
extern std::string BFileName;

int Ax[150];
int Bx[150];
int Px[150];


extern std::vector<std::vector<kmer_t>> kmer_sets; 

void read_array()
{
    std::ifstream input(primeFileName);

    for (int i = 0; i < 150; i++) {
        input >> Px[i];
    }
    std::ifstream input1(AFileName);

    for (int i = 0; i < 150; i++) {
        input1 >> Ax[i];
    }
    std::ifstream input2(BFileName);

    for (int i = 0; i < 150; i++) {
        input2 >> Bx[i];
    }
}


/* Get minimizer */
void recalculate_min_kmer (std::string ptr, kmer_t *m_kmer, int *fact, int *pos)
{
    kmer_t min_kmer=3074457345618258602, tmp_kmer=0;

    kmer_t kmer=0;

    // int i;

    int tracker = 0;

    int special_char= 0; //tracking N
    for(int i=0; i<KMER_LENGTH-1; i++) {

        if(special_char> 0)
        {
           special_char--;
        }
        if (ptr[i] == 'N' || ptr[i] == 'Y' || ptr[i] == 'S' || ptr[i] == 'R' || ptr[i] == 'I' || ptr[i] == 'E' || ptr[i] == 'K')
        {
            special_char= KMER_LENGTH;
        }

        kmer = kmer_shift(kmer, char_to_el(ptr[i]));
 
    }
    int start = 0;
    for(int i=KMER_LENGTH-1; i < w_size; i++) {
    ;
            if(special_char> 0)
            {
               special_char--;
            }
            if (ptr[i] == 'N' || ptr[i] == 'Y' || ptr[i] == 'S' || ptr[i] == 'R' || ptr[i] == 'I' || ptr[i] == 'E' || ptr[i] == 'K')
            {
                special_char= KMER_LENGTH;
            }

            kmer = kmer_shift(kmer, char_to_el(ptr[i]));
            
            if(special_char<= 0)
            {
                tracker = 1;
                
                if (min_kmer > kmer)
                {
                    min_kmer = kmer;
                    start = i - KMER_LENGTH;
                }
            }
            
    } 
    *m_kmer = min_kmer;
    *fact = tracker;
    *pos = start;
}


char convert_to_char (char given_char, int len, int id)
{
     if (given_char == 'A')
         return 'T';
     else if (given_char == 'C')
         return 'G';
     else if (given_char == 'G')
         return 'C';
     else if (given_char == 'T')
         return 'A';
     else if (given_char == 'N')
         return 'N';
     else if (given_char == 'Y')
         return 'Y';
     else if (given_char == 'S')
         return 'S';
     else if (given_char == 'R')
         return 'R';
     else if (given_char == 'I')
         return 'I';
     else if (given_char == 'E')
         return 'E';
     else if (given_char == 'K')
         return 'K';
     else {        

            return given_char;
     }
}

using Clock = std::chrono::high_resolution_clock;
using sec   = std::chrono::duration<double>;



// Cuda Kernels


// Convert each bp into 2 bit
#ifdef USE_CUDA
// cudaSetDevice(1);
__device__ inline ElType char_to_el_dev(char ch) {

    switch (ch) {
        case 'A': return (ElType)0; 
        case 'C': return (ElType)1;
        case 'G': return (ElType)2;
        case 'T': return (ElType)3;
        default:  return (ElType)0;
    }
}

__device__ inline kmer_t kmer_shift_dev(kmer_t kmer_in, ElType el) {
    return (kmer_t)(((kmer_in << 2) | (kmer_t)el) & (kmer_t)KMER_MASK);
}

__device__ inline kmer_t smer_shift_dev(kmer_t kmer_in, ElType el) {
    return (kmer_t)(((kmer_in << 2) | (kmer_t)el) & (kmer_t)KMER_MASK);
}
#endif

// Kernel for Minimizer per window
// Each thread handles one window of size w_size
#ifdef USE_CUDA
// #define smer_len 11   // s-mer length (fixed)
// #define KMER_LENGTH 15
__global__ void compute_minimizer_per_window_kernel(
    const char* seq, int seq_len, int w_size,
    kmer_t* out_min_kmer, int* out_tracker, int* out_pos)
{
    // Cuda thread Index; one thread handles one window
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    int smer_len = 11;
    int kmer_len = 15;


    int max_windows = seq_len - w_size + 1; // max windows = 10000-200 + 1
    if (idx < 0 || idx >= max_windows) return;

    const char* window_ptr = seq + idx; 
    // compute minimizer across this window
    kmer_t min_kmer = (kmer_t)3074457345618258602ULL;

    kmer_t min_smer = (kmer_t)3074457345618258602ULL;
    kmer_t kmer = 0;
    kmer_t smer = 0;

    int tracker = 0;
    int start_pos = 0;
    int special_char = 0;
    
    for (int i = 0; i < smer_len - 1; ++i) {
        char ch = window_ptr[i];

        if (special_char > 0) special_char--;
        if (ch == 'N' || ch == 'Y' || ch == 'S' ||
            ch == 'R' || ch == 'I' || ch == 'E' || ch == 'K') {
            special_char = smer_len;
        }

        ElType el = char_to_el_dev(ch);
        kmer = kmer_shift_dev(kmer, el);
        smer = smer_shift_dev(smer, el);
    }

    for (int i = smer_len - 1; i < kmer_len; ++i) {
        char ch = window_ptr[i];

        if (special_char > 0) special_char--;
        if (ch == 'N' || ch == 'Y' || ch == 'S' ||
            ch == 'R' || ch == 'I' || ch == 'E' || ch == 'K') {
            special_char = smer_len;
        }

        ElType el = char_to_el_dev(ch);
        kmer = kmer_shift_dev(kmer, el);
        smer = smer_shift_dev(smer, el);

        if (special_char <= 0) {
            tracker = 1;
            if (smer < min_smer) {
                min_smer = smer;
                start_pos = i - (smer_len - 1);  // i - 10
            }
        }
    }

    if (start_pos == 0 || start_pos == (kmer_len - smer_len)) {
        out_min_kmer[idx] = kmer;
        out_tracker[idx] = tracker;
        out_pos[idx] = start_pos;
    } else {
        out_min_kmer[idx] = 0;
        out_tracker[idx] = 0;
        out_pos[idx] = start_pos;
    }
}


//     out_min_kmer[idx] = min_kmer;
//     out_tracker[idx] = tracker;
//     out_pos[idx] = start_pos;
// }
#endif 

#ifdef USE_CUDA

static void compute_minimizers_on_gpu(const std::string &seq,
                                      int read_len, int w_size,
                                      std::vector<kmer_t> &forward_set,
                                      std::vector<int> &pos_set,
                                      std::vector<int> &forward_set_tracker)
{
    if (read_len < w_size) {
        return;
    }

    int windows = read_len - w_size + 1;

    char* d_seq = nullptr;
    kmer_t* d_min = nullptr;
    int* d_tracker = nullptr;
    int* d_pos = nullptr;

    cudaError_t cerr;
    cerr = cudaMalloc((void**)&d_seq, (size_t)read_len * sizeof(char));
    if (cerr != cudaSuccess) {
        fprintf(stderr, "cudaMalloc d_seq failed: %s\n", cudaGetErrorString(cerr));
        return;
    }
    cerr = cudaMalloc((void**)&d_min, (size_t)windows * sizeof(kmer_t));
    if (cerr != cudaSuccess) { cudaFree(d_seq); return; }
    cerr = cudaMalloc((void**)&d_tracker, (size_t)windows * sizeof(int));
    if (cerr != cudaSuccess) { cudaFree(d_seq); cudaFree(d_min); return; }
    cerr = cudaMalloc((void**)&d_pos, (size_t)windows * sizeof(int));
    if (cerr != cudaSuccess) { cudaFree(d_seq); cudaFree(d_min); cudaFree(d_tracker); return; }

    cerr = cudaMemcpy(d_seq, seq.data(), (size_t)read_len * sizeof(char), cudaMemcpyHostToDevice);
    if (cerr != cudaSuccess) { 
        fprintf(stderr, "cudaMemcpy d_seq failed: %s\n", cudaGetErrorString(cerr));
        cudaFree(d_seq); cudaFree(d_min); cudaFree(d_tracker); cudaFree(d_pos); 
        return;
    }

    int threads = 128;
    int blocks = (windows + threads - 1) / threads;

    compute_minimizer_per_window_kernel<<<blocks, threads>>>(d_seq, read_len, w_size, d_min, d_tracker, d_pos);


    cerr = cudaGetLastError();
    if (cerr != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(cerr));
        cudaFree(d_seq); cudaFree(d_min); cudaFree(d_tracker); cudaFree(d_pos);
        return;
    }

    std::vector<kmer_t> h_min(windows);
    std::vector<int> h_tracker(windows);
    std::vector<int> h_pos(windows);

    cerr = cudaMemcpy(h_min.data(), d_min, windows * sizeof(kmer_t), cudaMemcpyDeviceToHost);
    if (cerr != cudaSuccess) { fprintf(stderr, "cudaMemcpy back h_min failed\n"); }

    cerr = cudaMemcpy(h_tracker.data(), d_tracker, windows * sizeof(int), cudaMemcpyDeviceToHost);
    if (cerr != cudaSuccess) { fprintf(stderr, "cudaMemcpy back h_tracker failed\n"); }

    cerr = cudaMemcpy(h_pos.data(), d_pos, windows * sizeof(int), cudaMemcpyDeviceToHost);
    if (cerr != cudaSuccess) { fprintf(stderr, "cudaMemcpy back h_pos failed\n"); }

    // free device memory
    cudaFree(d_seq);
    cudaFree(d_min);
    cudaFree(d_tracker);
    cudaFree(d_pos);

    // populate forward_set, pos_set and tracker vectors to match original host logic:
    forward_set.clear(); forward_set.reserve(windows);
    pos_set.clear(); pos_set.reserve(windows);
    forward_set_tracker.clear(); forward_set_tracker.reserve(windows);

    for (int i = 0; i < windows; ++i) {
        forward_set.push_back(h_min[i]);
        pos_set.push_back(h_pos[i]);
        forward_set_tracker.push_back(h_tracker[i]);
    }
}
#endif // USE_CUDA

// Cuda Kernels for MinHash

// #define CUDA_CHECK(call) \
// do { \
//   cudaError_t err = (call); \
//   if (err != cudaSuccess) { \
//     fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
//     exit(1); \
//   } \
// } while(0)


// __global__
// void compute_min_per_trial_kernel(const kmer_t* d_kmers, long long k_begin, long long k_end,
//                                   const kmer_t* d_Ax, const kmer_t* d_Bx, const kmer_t* d_Px,
//                                   int no_trials,
//                                   kmer_t* d_min_val_out, kmer_t* d_corr_kmer_out)
// {
//     int trial = blockIdx.x * blockDim.x + threadIdx.x;
//     if (trial >= no_trials) return;

//     kmer_t best_val = (kmer_t)3074457345618258602ULL;
//     kmer_t best_kmer = 0;

//     kmer_t A = d_Ax[trial];
//     kmer_t B = d_Bx[trial];
//     kmer_t P = d_Px[trial];

//     for (long long kk = k_begin; kk < k_end; ++kk) {
//         kmer_t kmer = d_kmers[kk];
//         unsigned __int128 prod = (unsigned __int128)A * (unsigned __int128)kmer;
//         kmer_t val = (kmer_t)((prod + B) % P);
//         if (val < best_val) {
//             best_val = val;
//             best_kmer = kmer;
//         }
//     }

//     d_min_val_out[trial] = best_val;
//     d_corr_kmer_out[trial] = best_kmer;
// }



#include <cuda_runtime.h>
#include <limits>

struct Window {
    int start;
    int end;
};

__global__ void window_min_kernel(
    const kmer_t* kmers,
    const int* pos,
    const Window* windows_hash,
    int num_windows,
    int no_trials,
    const int* Ax,
    const int* Bx,
    const int* Px,
    kmer_t* out_min_vals
) {
    int w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= num_windows) return;

    int start = windows_hash[w].start;
    int end   = windows_hash[w].end;

    for (int trial = 0; trial < no_trials; ++trial) {
        kmer_t min_val = (kmer_t)3074457345618258602;

        for (int i = start; i < end; ++i) {
            kmer_t val = (Ax[trial] * kmers[i] + Bx[trial]) % Px[trial];
            if (val < min_val) min_val = val;
        }

        out_min_vals[w * no_trials + trial] = min_val;
    }
}





void print_trial_map_stats(
    const std::vector<std::unordered_map<kmer_t, std::vector<int>>>& trial_maps)
{
    for (size_t t = 0; t < trial_maps.size(); ++t)
    {
        const auto &mp = trial_maps[t];

        size_t total_keys = mp.size();
        size_t min_vec = std::numeric_limits<size_t>::max();
        size_t max_vec = 0;
        size_t total_subject_count = 0;

        for (const auto &kv : mp)
        {
            const std::vector<int> &v = kv.second;
            size_t sz = v.size();

            total_subject_count += sz;

            if (sz < min_vec) min_vec = sz;
            if (sz > max_vec) max_vec = sz;
        }

        if (total_keys == 0)   // handle empty trial
        {
            min_vec = 0;
            max_vec = 0;
        }

        // std::cout << "Trial " << t << " stats:\n";
        std::cout << "Trial:" << t << " Keys: " << total_keys <<  "Total subjects" << total_subject_count <<"\n";
        // std::cout << "  Min subject list size: " << min_vec << "\n";
        // std::cout << "  Max subject list size: " << max_vec;
        // std::cout << "  Total subjects across all keys: " << total_subject_count << "\n\n";
    }
}

size_t estimate_map_memory(const std::unordered_map<kmer_t, std::vector<int>>& map) {
    size_t total = sizeof(map); // unordered_map object itself

    for (const auto& kv : map) {
        total += sizeof(kv.first);                  // key
        total += sizeof(kv.second);                 // vector object (3 pointers)
        total += sizeof(void*);                     // pointer to next node (approximation)
        total += kv.second.capacity() * sizeof(int); // vector data
    }

    return total;
}

// Estimate memory for the vector of maps
size_t memory_usage_trial_maps(const std::vector<std::unordered_map<kmer_t, std::vector<int>>>& trial_maps) {
    size_t total = sizeof(trial_maps); // vector object itself

    for (const auto& map : trial_maps) {
        total += estimate_map_memory(map);

        // std::cout << "Estimated memory for one trial map: " << total << " bytes\n";
    }

    return total;
}




struct BatchedRead {
    int subject_id;
    int offset_kmer;
    int offset_pos;
    int count;
};

struct Window_batch {
    int start;
    int end;
    int read_id;
};

__global__
void window_min_kernel_batch(
    const kmer_t* kmers,
    const int* pos,
    const Window_batch* windows,
    int num_windows,
    int no_trials,
    const int* Ax,
    const int* Bx,
    const int* Px,
    kmer_t* out
) {
    int wid = blockIdx.x * blockDim.x + threadIdx.x;
    if (wid >= num_windows) return;

    Window_batch w = windows[wid];

    for (int t = 0; t < no_trials; t++) {
        kmer_t minv = 3074457345618258602ULL;

        for (int i = w.start; i < w.end; i++) {
            kmer_t h = (Ax[t] * kmers[i] + Bx[t]) % Px[t];
            if (h < minv) minv = h;
        }
        out[wid * no_trials + t] = minv;
    }
}


void run_batched_min_hash(
    const std::vector<BatchedRead>& batch_reads,
    const std::vector<kmer_t>& batch_kmers,
    const std::vector<int>& batch_pos,
    int no_trials,
    int* Ax, int* Bx, int* Px,
    int read_length,
    std::vector<std::unordered_map<kmer_t, std::vector<int>>>& trial_maps
    ) {
    // --------------------
    // BUILD WINDOWS
    // --------------------
    std::vector<Window_batch> windows;

    for (int r = 0; r < (int)batch_reads.size(); r++) {
        const auto& br = batch_reads[r];
        int j = 0, k = 0;

        while (j < br.count && k < br.count) {
            if (batch_pos[br.offset_pos + k] -
                batch_pos[br.offset_pos + j] <= read_length) {
                k++;
            } else {
                windows.push_back({
                    br.offset_kmer + j,
                    br.offset_kmer + k,
                    r
                });
                j++;
            }
        }
        if (j < k) {
            windows.push_back({
                br.offset_kmer + j,
                br.offset_kmer + k,
                r
            });
        }
    }

    int num_windows = windows.size();
    if (num_windows == 0) return;

    // --------------------
    // DEVICE ALLOCATIONS
    // --------------------
    Window_batch* d_windows;
    kmer_t* d_kmers;
    int* d_pos;
    int* d_Ax;
    int* d_Bx;
    int* d_Px;
    kmer_t* d_out;

    cudaMalloc(&d_windows, num_windows * sizeof(Window_batch));
    cudaMalloc(&d_kmers, batch_kmers.size() * sizeof(kmer_t));
    cudaMalloc(&d_pos, batch_pos.size() * sizeof(int));
    cudaMalloc(&d_Ax, no_trials * sizeof(int));
    cudaMalloc(&d_Bx, no_trials * sizeof(int));
    cudaMalloc(&d_Px, no_trials * sizeof(int));
    cudaMalloc(&d_out, num_windows * no_trials * sizeof(kmer_t));

    cudaMemcpy(d_windows, windows.data(),
               num_windows * sizeof(Window_batch), cudaMemcpyHostToDevice);
    cudaMemcpy(d_kmers, batch_kmers.data(),
               batch_kmers.size() * sizeof(kmer_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos, batch_pos.data(),
               batch_pos.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Ax, Ax, no_trials * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Bx, Bx, no_trials * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Px, Px, no_trials * sizeof(int), cudaMemcpyHostToDevice);

    // --------------------
    // KERNEL LAUNCH
    // --------------------
    int threads = 256;
    int blocks = (num_windows + threads - 1) / threads;

    // std::cout << "Windows: " << num_windows  << " Blocks: " << blocks << "\n";


    window_min_kernel_batch<<<blocks, threads>>>(
        d_kmers,
        d_pos,
        d_windows,
        num_windows,
        no_trials,
        d_Ax, d_Bx, d_Px,
        d_out
    );

    cudaDeviceSynchronize();

    std::vector<kmer_t> h_out(num_windows * no_trials);
    cudaMemcpy(h_out.data(), d_out,
               h_out.size() * sizeof(kmer_t),
               cudaMemcpyDeviceToHost);

    
    std::vector<std::vector<std::unordered_set<kmer_t>>> per_read_trial_keys(
        batch_reads.size(),
        std::vector<std::unordered_set<kmer_t>>(no_trials)
    );
    for (int w = 0; w < num_windows; ++w) {
    int read_id = windows[w].read_id;

    for (int t = 0; t < no_trials; ++t) {
        kmer_t key = h_out[w * no_trials + t];
        per_read_trial_keys[read_id][t].insert(key);
    }
}
for (int r = 0; r < (int)batch_reads.size(); ++r) {
    int subject_id = batch_reads[r].subject_id;

    for (int t = 0; t < no_trials; ++t) {
        for (const auto& key : per_read_trial_keys[r][t]) {
            trial_maps[t][key].push_back(subject_id);
        }
    }
}


    // for (int w = 0; w < num_windows; w++) {
    //     int read_id = windows[w].read_id;
    //     int subject_id = batch_reads[read_id].subject_id;

    //     for (int t = 0; t < no_trials; t++) {
    //         kmer_t val = h_out[w * no_trials + t];
    //         trial_maps[t][val].push_back(subject_id);
    //     }
    // }

    cudaFree(d_windows);
    cudaFree(d_kmers);
    cudaFree(d_pos);
    cudaFree(d_Ax);
    cudaFree(d_Bx);
    cudaFree(d_Px);
    cudaFree(d_out);
}


//  Sliding_window that uses GPU
void Sliding_window (char *ptr, size_t length, int *M_for_individual_process, int *num_subjects,
                     std::vector<MinHashPairs> &initial_sets, int s_index)
{

    size_t p=0;
    int max_set_length = 0;
    std::vector<kmer_t> forward_set; // set of kmers from forward strand (from GPU)
    std::vector<int> pos_set; // position of kmers (from GPU)
    std::vector<kmer_t> kmer_set; 
    std::vector<kmer_t> kmer_set_pos;
    std::vector<kmer_t> set_of_distinct_kmers;
    std::vector<int> set_of_distinct_pos;
    std::vector<kmer_t> set_of_distinct_kmers_rev;
    std::vector<int> set_of_distinct_pos_rev;
    std::vector<kmer_t> set_of_dist_kmers;
    std::vector<MinHashPairs> set_of_dist_minhash_pairs;
    std::vector<int> forward_set_tracker;
    int total_subjects = 0;
  
    for(; ptr[p]!='>' && p<length; p++) { }

    // kmer_t kmer = 0; 
    kmer_t min_kmer = 0;
    kmer_t rev_min_kmer = 0;
    // kmer_t min_lmer_freq = 0; 
    int min_pos = 0;
    int min_tracker = 0;

    int total_k_mer = 0;
    int total_min_hashes = 0;

    using Clock = std::chrono::high_resolution_clock;
    using sec   = std::chrono::duration<double>;
    sec forward_time{0}, reverse_time{0}, distinct_time{0}, minhash_time{0}, check_time{0};
    sec calculation_time{0}, hash_check_time{0}, hash_populate_time{0};

    int count = 0;

    std::vector<std::unordered_map<kmer_t, std::vector<int>>> trial_maps;
    trial_maps.resize(no_trials);

    int max_distinct_kmers = 0;
    int count_window = 0;

    constexpr int BATCH_SIZE = 100;

    std::vector<kmer_t> batch_kmers;
    std::vector<int>    batch_pos;
    std::vector<BatchedRead> batch_reads;


    while(p<length) {

        // Iterate across the entries

        count += 1;


        if (p > 0 && ptr[p-1] == '>') {
            --p; // step back so assert sees >
        }
        assert(ptr[p] == '>');
        // Skip to next line


        // skip header name
        for(; p<length && ptr[p]!='\n'; p++) { } 
        size_t seq_start = ++p;
        size_t seq_len = 0;
        
        while (seq_start + seq_len < length && ptr[seq_start + seq_len] != '>') seq_len++;
        p++; 
        total_subjects++;

        if (p + w_size > length) break;

        std::string s; 
        int i;

        // s.reserve(seq_len);
        std::string str; 
        // str.reserve(seq_len);
        int read_len = 0;
        size_t q = seq_start;

        while (q < seq_start + seq_len && !isspace(ptr[q])) {
            char orig = ptr[q];
            str.push_back(orig);               // original sequence used for reverse computation
            s.push_back(convert_to_char(orig, read_len, total_subjects)); // forward complement version
            q++;
            read_len++;

            // std::cout << s << "\n";

        }

        p = seq_start + read_len;

        // GPU: compute forward_set, pos_set, forward_set_tracker for the full read 'str'
#ifdef USE_CUDA
        auto gstart = Clock::now();
        compute_minimizers_on_gpu(str, read_len, w_size, forward_set, pos_set, forward_set_tracker);
        auto gend = Clock::now();
        forward_time += gend - gstart;
#else

        auto t1 = Clock::now();
        for(i=0; !isspace(ptr[p]) && i<w_size-1; i++) {
      
            s.push_back(convert_to_char(ptr[p], read_len, total_subjects));
            str.push_back(ptr[p]);
            p++;
            read_len++;
            
            
        }
        while(p<length && !isspace(ptr[p])) {  
            s.push_back(convert_to_char(ptr[p], read_len, total_subjects));
            str.push_back(ptr[p]);
            //s.push_back(ptr[p]);
            p++;
            read_len++;
            //rev_set.push_back(min_kmer);
            recalculate_min_kmer(str.substr(read_len - w_size, w_size), &min_kmer, &min_tracker, &min_pos);
            forward_set.push_back(min_kmer);
            pos_set.push_back(min_pos);
            if (min_tracker > 0)
            {
                forward_set_tracker.push_back(1);
            }
            else
            {
                forward_set_tracker.push_back(0);
            }
            
        } 

        auto t2 = Clock::now();
        forward_time += t2 - t1;
#endif

        auto t_rev_start = Clock::now();
        // reverse(s) etc. (same as original):
        reverse(s.begin(), s.end());

        #ifdef USE_CUDA

        std::string s_rev = s;

        std::vector<kmer_t> reverse_set;
        std::vector<int> reverse_pos;
        std::vector<int> reverse_tracker;

        compute_minimizers_on_gpu(s_rev, read_len, w_size, reverse_set, reverse_pos, reverse_tracker);

        int windows = (read_len >= w_size) ? (read_len - w_size + 1) : 0;

        // std::cout << "Reverse set " << reverse_set.size() << " Window " << windows << w_size << "\n";


        if ((int)reverse_set.size() != windows || (int)forward_set.size() != windows) {

            std::cout<<"Yes GPU failed";
            int length_tracker_local = 0;
            int itr_local = 0;
            for(int i=0; i<w_size-1; i++) { length_tracker_local++; }
            while(length_tracker_local < read_len) {
                length_tracker_local++;
                recalculate_min_kmer(s.substr(length_tracker_local - w_size, w_size), &rev_min_kmer, &min_tracker, &min_pos);

                if (rev_min_kmer <= forward_set[read_len - w_size - itr_local]) {
                    if (forward_set_tracker[read_len - w_size - itr_local] == 1) {
                        if (rev_min_kmer != 3074457345618258602ULL) {
                            kmer_set.push_back(rev_min_kmer);
                            kmer_set_pos.push_back(read_len - w_size - itr_local + (w_size - 1) - (min_pos + KMER_LENGTH - 1));
                        }
                    }
                } else {
                    if (forward_set_tracker[read_len - w_size - itr_local] == 1) {
                        kmer_set.push_back(forward_set[read_len - w_size - itr_local]);
                        kmer_set_pos.push_back(pos_set[read_len - w_size - itr_local] + read_len - w_size - itr_local);
                    }
                }
                itr_local++;
            }
        } else {

            for (int itr_idx = 0; itr_idx < windows; ++itr_idx) {
                int fidx = read_len - w_size - itr_idx; // index into forward_set / pos_set / forward_set_tracker
                int ridx = itr_idx;                      // index into reverse_set / reverse_pos / reverse_tracker

                kmer_t rev_k = reverse_set[ridx];
                int f_tracker = forward_set_tracker[fidx];
                if (rev_k <= forward_set[fidx]) {
                    if (f_tracker == 1) {
                        if (rev_k != (kmer_t)3074457345618258602ULL) {

                            int minpos = reverse_pos[ridx];
                            int computed_pos = read_len - w_size - itr_idx + (w_size - 1) - (minpos + (int)KMER_LENGTH - 1);
                            kmer_set.push_back(rev_k);
                            kmer_set_pos.push_back((kmer_t)computed_pos);
                        }
                    }
                } else {
                    if (f_tracker == 1) {
                        int computed_pos = pos_set[fidx] + read_len - w_size - itr_idx;
                        kmer_set.push_back(forward_set[fidx]);
                        kmer_set_pos.push_back((kmer_t)computed_pos);
                    }
                }
            }
        }
#else
        // If CUDA not available,  back to original host-side reverse loop

        int length_tracker = 0;
        int itr = 0;
        for(int i=0; i<w_size-1; i++) { length_tracker++; }

        while(length_tracker < read_len) {
            length_tracker++;
            recalculate_min_kmer(s.substr(length_tracker - w_size, w_size), &rev_min_kmer, &min_tracker, &min_pos);

            if (rev_min_kmer <= forward_set[read_len - w_size - itr]) {
                if (forward_set_tracker[read_len - w_size - itr] == 1) {
                    if (rev_min_kmer != 3074457345618258602ULL) {
                        kmer_set.push_back(rev_min_kmer);
                        kmer_set_pos.push_back(read_len - w_size - itr + (w_size - 1) - (min_pos + KMER_LENGTH - 1));
                    }
                }
            } else {
                if (forward_set_tracker[read_len - w_size - itr] == 1) {
                    kmer_set.push_back(forward_set[read_len - w_size - itr]);
                    kmer_set_pos.push_back(pos_set[read_len - w_size - itr] + read_len - w_size - itr);
                }
            }
            itr++;
        }
#endif
        auto t_rev_end = Clock::now();
        reverse_time += t_rev_end - t_rev_start;


        auto t_dist_start = Clock::now();
        // std::cout <<"Syncmer size: " << kmer_set.size() << "\n";



        if (read_len >= 1000 && kmer_set.size() > 0) {
            kmer_t prev = kmer_set[0];
            int i;
            for (i = 1; i < (int)kmer_set.size(); i++) {
                if (kmer_set[i] != prev) {
                    set_of_distinct_kmers.push_back(prev);
                    prev = kmer_set[i];
                    set_of_distinct_pos.push_back(kmer_set_pos[i-1]);
                }
            }
            set_of_distinct_kmers.push_back(prev);
            set_of_distinct_pos.push_back(kmer_set_pos[i-1]);
            set_of_distinct_kmers.push_back(prev);
            set_of_distinct_pos.push_back(kmer_set_pos[i-1]);
            //printf("Wsize Subject =%d\n", set_of_distinct_kmers.size());
            kmer_set.clear();
            kmer_set.shrink_to_fit();
            kmer_set_pos.clear();
            kmer_set_pos.shrink_to_fit();
        }


        for (int reverse = (int)set_of_distinct_pos.size() - 1; reverse >= 0; --reverse) {
            set_of_distinct_pos_rev.push_back(set_of_distinct_pos[reverse]);
            set_of_distinct_kmers_rev.push_back(set_of_distinct_kmers[reverse]);
        }

        // for (int reverse = (int)kmer_set_pos.size() - 1; reverse >= 0; --reverse) {
        //     set_of_distinct_pos_rev.push_back(kmer_set_pos[reverse]);
        //     set_of_distinct_kmers_rev.push_back(kmer_set[reverse]);
        // }


        kmer_set.clear(); kmer_set.shrink_to_fit();
        kmer_set_pos.clear(); kmer_set_pos.shrink_to_fit();

        // std::cout << "Distinct size: " << set_of_distinct_kmers_rev.size() << "\n";

        // break;

        // size_t total_mem =
        //     sizeof(set_of_distinct_kmers_rev) +
        //     set_of_distinct_kmers_rev.capacity() * sizeof(kmer_t);

        // total_k_mer = total_k_mer + set_of_distinct_kmers.size();

        // if (max_distinct_kmers < (int)set_of_distinct_kmers.size())
        // {
        //     max_distinct_kmers = set_of_distinct_kmers.size();
        // }


        // total_k_mer

        BatchedRead br;
        br.subject_id   = total_subjects + s_index;
        br.offset_kmer  = batch_kmers.size();
        br.offset_pos   = batch_pos.size();
        br.count        = set_of_distinct_kmers_rev.size();

        batch_kmers.insert(batch_kmers.end(),
                        set_of_distinct_kmers_rev.begin(),
                        set_of_distinct_kmers_rev.end());

        batch_pos.insert(batch_pos.end(),
                        set_of_distinct_pos_rev.begin(),
                        set_of_distinct_pos_rev.end());

        batch_reads.push_back(br);

        // clear per-read buffers (as before)
        set_of_distinct_kmers.clear();
        set_of_distinct_pos.clear();
        set_of_distinct_kmers_rev.clear();
        set_of_distinct_pos_rev.clear();

        // ====================
        // LAUNCH WHEN BATCH FULL
        // ====================

        if ((int)batch_reads.size() == BATCH_SIZE) {
           auto t1 = Clock::now();

            run_batched_min_hash(
                batch_reads,
                batch_kmers,
                batch_pos,
                no_trials,
                Ax, Bx, Px,
                read_length,
                trial_maps
            );

            batch_reads.clear();
            batch_kmers.clear();
            batch_pos.clear();
        
        
            auto t2 = Clock::now();
            minhash_time += t2 - t1;
        }
        
        
        set_of_distinct_kmers.clear();
        set_of_distinct_kmers.shrink_to_fit();
        // window_min_sets.clear();
        // window_min_sets.shrink_to_fit();

        
        set_of_distinct_pos.clear();
        set_of_distinct_pos.shrink_to_fit();
        set_of_distinct_kmers_rev.clear();
        set_of_distinct_kmers_rev.shrink_to_fit();
        
        set_of_distinct_pos_rev.clear();
        set_of_distinct_pos_rev.shrink_to_fit();
        set_of_dist_kmers.clear();
        set_of_dist_kmers.shrink_to_fit();
        p++; 
        p++;

    //    auto t2 = Clock::now();
        // minhash_time += t2 - t1;
    }

    print_trial_map_stats(trial_maps);




    size_t mem = memory_usage_trial_maps(trial_maps);
    std::cout << "Memory used by trial_maps = " << mem << " bytes\n";


        
    std::cout << "\nTiming \n";
    // std::cout << count_window;


    std::cout << "Total minimizers " << total_k_mer << " s\n";
    std::cout << "Total minhashes " << count_window << " \n";
    std::cout << "Maximum  minimizers across read" << max_distinct_kmers << " s\n";
    std::cout << "Maximum  minhashes across read" << max_distinct_kmers << " s\n";



    std::cout << "Forward minimizers (GPU/host): " << forward_time.count() << " s\n";
    std::cout << "Reverse minimizers: " << reverse_time.count() << " s\n";
    // std::cout << "Distinct minimizers: " << distinct_time.count() << " s\n";
    std::cout << "MinHash: " << minhash_time.count()  << " s\n\n";

    std::cout << "MinHash population: " << hash_populate_time.count()  << " s\n\n";
    std::cout << "MinHash window + minval calculation: " << calculation_time.count()  << " s\n\n";

}



void generate_set_of_subjects (char *read_data, size_t length, int s_index, char *r_data, size_t r_length, int start_index, int total_q, int *M_final, int *num_subjects)
{
    double time_l1 = omp_get_wtime();
    read_array();

    std::vector< MinHashPairs > minhash_from_set_of_subjects;

    int M_for_individual_processes = 0;
    int n_subjects;
    // int total_subjects;

    // Sliding_window - compute MinHash sketches from minimizer k-mers across subjects, with strand direction 
    auto t1 = Clock::now();

    Sliding_window (read_data, length, &M_for_individual_processes, &n_subjects, minhash_from_set_of_subjects, s_index);

}
