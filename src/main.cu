// Main file
// Sketch-Arch

// main.cpp (single-process, non-MPI version)

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../includes/JEM.h"
#include "../includes/timers.h"   // optional: keep if your project uses it

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::ostringstream;

long int MAX_KMER_COUNT = 0;
int rank = 0, size = 1; // single-process defaults (kept for compatibility)
int coverage = 0;

// Files
std::string inputFileName;
std::string queryFileName;
std::string primeFileName;
std::string AFileName;
std::string BFileName;
int read_length = 0;
int w_size = 0;
int no_trials = 0;

int num_batch_transfers = 0;

std::vector<KmerPairs> kmer_proc_buf;
std::vector<std::vector<kmer_t>> kmer_sets;

// std::vector<kmer_t> global_kmer_frequency(15,0);

// Hash seeds
uint64_t hasha = 68111;
uint64_t hashb = 105929;

void parseCommandLine(const int argc, char * const argv[]);

std::string readFileIntoString(const std::string& path) {
    ifstream input_file(path, std::ios::in | std::ios::binary);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '" << path << "'" << endl;
        exit(EXIT_FAILURE);
    }
    return string((std::istreambuf_iterator<char>(input_file)), std::istreambuf_iterator<char>());
}

// Safer get_file_size
long long get_file_size(const std::string &filename) {
    FILE *p_file = fopen(filename.c_str(), "rb");
    if (!p_file) {
        perror(("fopen: " + filename).c_str());
        return -1;
    }
    if (fseek(p_file, 0, SEEK_END) != 0) {
        perror(("fseek: " + filename).c_str());
        fclose(p_file);
        return -1;
    }
    long long size = ftell(p_file);
    fclose(p_file);
    return size;
}

int main(int argc, char **argv) {

    parseCommandLine(argc, argv);
    double time_l1 = omp_get_wtime();

    cout << "Subject" << endl;

    input_read_data sdata = perform_input_reading(inputFileName, 5000000);
    int M = 0;
    int total_subjects = 0;

   //  cout << "Query" << endl;
   //  input_read_data qdata = perform_input_reading(queryFileName, read_length);
   //  std::cout << "L (read_length): " << read_length << std::endl;

    generate_set_of_subjects (sdata.read_data, sdata.read_data_size, sdata.start_index, sdata.read_data, sdata.read_data_size, sdata.start_index, sdata.total, &M, &total_subjects);

    // int total_number_of_subs_in_p = kmer_sets.size();
    // kmer_t** Hash_table = new kmer_t*[total_subjects];

    // int total_hash_functions = 150;
    // for (int i = 0; i < total_subjects; i++) {
    //     Hash_table[i] = new kmer_t[total_hash_functions];
    // }
    // genereate_hash_table(M, total_subjects, Hash_table);


    double time_l2 = omp_get_wtime();
    std::cout << "Elapsed time (s): " << (time_l2 - time_l1) << std::endl;

    // free buffers returned from perform_input_reading if needed
    if (sdata.read_data) {
        free(sdata.read_data);
        sdata.read_data = nullptr;
    }
   //  if (qdata.read_data) {
   //      free(qdata.read_data);
   //      qdata.read_data = nullptr;
   //  }

    return 0;
}



void parseCommandLine(const int argc, char * const argv[])
{
  int ret;

  while ((ret = getopt(argc, argv, "s:q:p:a:b:t:w:l:")) != -1) {
    switch (ret) {
    case 's':
       inputFileName.assign(optarg);
       std::cout << inputFileName << std::endl;
       break;
    case 'q':
       queryFileName.assign(optarg);
       std::cout << inputFileName << std::endl;
       break;
    case 'p':
       primeFileName.assign(optarg);
       std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 'a':
       AFileName.assign(optarg);
       std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 'b':
       BFileName.assign(optarg);
       std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 't':
       no_trials = atoi(optarg);
       std::cout << no_trials << std::endl;
       break;
    case 'w':
       w_size = atoi(optarg);
       std::cout << no_trials << std::endl;
       std::cout << "W (window size): " << w_size << std::endl;
       break;
    case 'l':
       read_length = atoi(optarg);
       std::cout << read_length << std::endl;
       break;
    
    
    default:
       assert(0 && "Should not reach here!!");
       break;
    }
  }


} 

