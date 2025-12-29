// JEM.h — MPI removed, single-process read API
// JEM-Mapper: A C++ implementation for Jaccard Estimate MinHash-based sequence-to-sequence mapping
// (authors & license header omitted for brevity — keep your original header text if desired)

#ifndef JEM_H
#define JEM_H

#include <stdio.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <queue>          // std::priority_queue
#include <string>
#include <assert.h>
#include <cstring>
#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <atomic>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <cstdint>

// #define KMER_LENGTH                    (WINDW_SIZE+1)
#define KMER_LENGTH                    (WINDW_SIZE)

#define LMER_SIZE                      (pow(4, LMER_LENGTH))
#define MN_LENGTH                      (KMER_LENGTH-1)

#define MOD 2147483647
#define HL 31

typedef uint64_t kmer_t;

using BasePair = uint8_t;
typedef BasePair ElType;
typedef uint64_t lmer_t;

#define KMER_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*KMER_LENGTH))
#define LMER_MASK           ((~0UL)) >> ((sizeof(lmer_t)*8) - (2*KMER_LENGTH))
#define MN_MASK             ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*MN_LENGTH))
#define SUCC_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*(KMER_LENGTH-1)))

extern int rank, size;  // keep for compatibility; set these in your single-process main (rank=0,size=1)

template<typename T>
inline bool
operator == (const std::vector<T>& v1,
             const std::vector<T>& v2) {
  if(v1.size() != v2.size()) {
    return false;
  }
  return std::equal(v1.begin(), v1.end(), v2.begin());
}

inline char el_to_char(unsigned x) {
      static char symbols[] = {'A', 'C', 'T', 'G', '\0', '\0', 'N'};
        return symbols[(ElType)x];
}

inline kmer_t mnmer_shift(kmer_t kmer_in, 
                        ElType el) {
  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)MN_MASK;
}

inline kmer_t mn_extract_pred(kmer_t kmer_in, size_t dist) {
    kmer_t mask = ((kmer_t)1 << (dist*2))-1;
    size_t rem = MN_LENGTH - dist;
    kmer_t pmask = mask << (rem*2);
    return ((pmask & kmer_in) >> (rem*2));
}

inline kmer_t mn_extract_succ(kmer_t kmer_in, size_t dist) {
    kmer_t mask = ((kmer_t)1 << (dist*2))-1;
    return (mask & kmer_in);
}

class Comp_rev{
    const std::vector<std::pair<int, int>> & _v;
  public:
    Comp_rev(const std::vector<std::pair<int, int>> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return ((_v[i].second > _v[j].second) ||
                 ((_v[i].second == _v[j].second) &&
                  ( _v[i].first > _v[j].first))
                );
   }
};

class Comp_rev1{
    const std::vector<std::pair<kmer_t, int>> & _v;
  public:
    Comp_rev1(const std::vector<std::pair<kmer_t, int>> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return ((_v[i].first < _v[j].first) ||
                 ((_v[i].first == _v[j].first) &&
                  ( _v[i].second > _v[j].second))
                );
   }
};

struct WireInfo {
  int suffix_id;
  int offset_in_suffix;
  int count;
  WireInfo()=default;
};

struct PrefixInfo {
  int prefix_pos;
  int num_wires;
  PrefixInfo()=default;
};

typedef struct begin_kmer_id
{
    kmer_t node;
    int terminal_prefix_id;
} BeginMN;

struct SubjectMapping{
    int s_id;
    int q_id;
    int extn;
};

typedef struct __attribute__ ((__packed__)) kmer_pairs
{ 
  kmer_t seq;
  int k_count;
} KmerPairs;
static_assert(sizeof(KmerPairs) == (sizeof(kmer_t)+sizeof(int)), "struct KmerPairs shouldn't be padded");

typedef struct __attribute__ ((__packed__)) MinHash_pairs
{ 
  int trial;
  kmer_t seq;
  int start_ind;
  int subject_id;
} MinHashPairs;
static_assert(sizeof(MinHashPairs) == sizeof(int)+(sizeof(kmer_t)+sizeof(int)+ sizeof(int)), "struct MinHashPairs shouldn't be padded");

typedef struct __attribute__ ((__packed__)) Top_Hit
{ 
  int sub;
  int score;
} TopHit;
static_assert(sizeof(MinHashPairs) == sizeof(int)+(sizeof(kmer_t)+sizeof(int)+ sizeof(int)), "struct MinHashPairs shouldn't be padded");

typedef struct Rd_Sequence
{
    // Read data
    char *read_data;
    size_t read_data_size;
    int start_index;
    int local_count;
    int total;
} input_read_data;

typedef struct each_contig_entry
{
       long int my_contig_id;
       std::string contig_data;
} contig_entry;

typedef struct contig_series
{
       int contig_count;
       std::vector < contig_entry> c_series;
} contig_thrd_list;

class Comp{
    const std::vector<kmer_t> & _v;
  public:
    Comp(const std::vector<kmer_t> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return _v[i] < _v[j];
   }
};

inline ElType kmerel(kmer_t kmer, unsigned pos) {
  assert(pos < KMER_LENGTH);
  return ElType(((kmer) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline ElType lmerel(lmer_t kmer, unsigned pos) {
  assert(pos < LMER_LENGTH);
  return ElType(((kmer) >> ((LMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline lmer_t kmer_to_lmer(kmer_t kmer_in, unsigned pos, lmer_t kmer_out) {
  assert(pos < KMER_LENGTH);
  return (lmer_t)((kmer_out<<2) | (lmer_t)(ElType(((kmer_in) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3)))) & (lmer_t)LMER_MASK;
}

inline ElType char_to_el(char ch) {
    return (ElType)((((ElType)ch)>>1) & 0x7);
}

inline kmer_t kmer_cons(kmer_t kmer_in, 
                        unsigned pos,
                        ElType el) {
  assert(pos < KMER_LENGTH);
  return kmer_in | ((kmer_t)el << ((KMER_LENGTH-pos-1)*2));
}

inline kmer_t kmer_shift(kmer_t kmer_in, 
                        ElType el) {
  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)KMER_MASK;
}

inline kmer_t smer_shift(kmer_t kmer_in, 
                        ElType el) {
  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)KMER_MASK;
}

inline kmer_t mod_kmer(kmer_t kmer_in, 
                        kmer_t kmer_out) {
  return (kmer_t)((kmer_in<<30) | (kmer_t) kmer_out) & (kmer_t)KMER_MASK;
}

inline lmer_t lmer_shift(lmer_t kmer_in, 
                        ElType el) {
  return (lmer_t)((kmer_in<<2) | (lmer_t)el) & (lmer_t)LMER_MASK;
}

inline kmer_t tokmer(const char *kmer_str,
                     int kmer_len) {
  int i;
  assert(kmer_len == KMER_LENGTH);
  assert(kmer_str != NULL);
  assert(kmer_str[kmer_len] == '\0');

  kmer_t km = 0;
  for(i=0; i<kmer_len; i++) {
    ElType el = char_to_el(kmer_str[i]);
    km = kmer_cons(km, i, el);
  }  
  return km;
}

template <typename T> 
inline long uhash31( uint64_t a, uint64_t b, T x)
{
  T result;
  long lresult;  
  result=(a * x) + b;
  result = ((result >> HL) + result) & MOD;
  lresult=(long) result; 
  return(lresult);
}

inline long hash31(long long a, long long b, long long x)
{
  long long result;
  long lresult;  
  result=(a * x) + b;
  result = ((result >> HL) + result) & MOD;
  lresult=(long) result; 
  return(lresult);
}

/* ------------------------
   Function declarations
   ------------------------ */

/* hashing, utilities */
int hash_fcn (const char* word, unsigned M);
void parse_alphabet_file (FILE * fp);
void free_reads (int num_reads);
bool is_sigma (char c);
void change_to_num (char *pred_kmer, int *num);
void change_to_char (char *pred_kmer, int num);
int find_max (int array[]);
int find_min (int array[], int *min);
int compare (const void * a, const void * b);
int FindIndex( const int a[], int size, int value);


char* DivideReads(const std::string &filename,
                  uint64_t *nlines,
                  size_t *data_size);

input_read_data perform_input_reading (const std::string &fileName, int read_length);
// input_read_data perform_input_reading(const std::string &fileName, int read_length)
/* legacy / alternative variant — kept but converted to single-process signature */
input_read_data perform_input_readingC (const std::string &fileName, int read_length);

/* other previously-declared functions — left intact */
void Sliding_window_l (const char *ptr, size_t length);
void Sliding_window (char *ptr, size_t length, int *M_for_individual_process, int *num_subjects,
                     std::vector<MinHashPairs> &initial_sets, std::vector<int> &subject_size, int s_index);
                     
void process_remaining_kmers(
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts); 

int convert_hash_fcn (const char* word);
unsigned int hash_str(const char *str);

int convert_to_int (char given_char);
std::string readFileIntoString(const std::string& path);
long long get_file_size(std::string filename);
void read_array();
void get_hash_value(kmer_t **A1, int M, kmer_t **Prefix);
void generate_set_of_subjects (char *ptr, size_t length, int s_index,char *read_data, size_t r_length, int start_index, int total_q, int *M_final, int *num_subjects);
void genereate_hash_table(int M, int total_subjects, kmer_t **Ag_Hash_Table);
void get_hash_value_queires(std::vector<std::vector<kmer_t>> &modified_sets, kmer_t **A1);
void Sliding_window_queires (char *ptr, size_t length, int *num_queries,
                     std::vector<std::unordered_map<kmer_t, std::unordered_map<int, std::vector<std::pair<int, int>>>>> Tl, int start, int total_subjects, int total_q, std::vector< int > subject_size, std::vector<SubjectMapping> &subjects);
void generate_set_of_queries (const char *read_data, size_t length, int start_index, int total_subjects, int M, kmer_t **Ag_Hash_Table);
void generate_modified_set_queries(int M, std::vector<std::vector<kmer_t>> &previous_sets);
char convert_to_char (char given_char);
void prefix_val(kmer_t **A1, int M);

#endif // JEM_H
