// single_process_io.cpp
// Replacement for MPI-based input reader — single-process / non-MPI version.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <chrono>
#include <string>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include "../includes/JEM.h"   // keeps your input_read_data type and other declarations

// preserved global timing variables (now measured locally)
double file_open_time=0.0, global_file_open_time=0.0;
double file_read_time=0.0, global_file_read_time=0.0;
double file_copy_time=0.0, global_file_copy_time=0.0;
double file_close_time=0.0, global_file_close_time=0.0;
double input_process_time=0.0, global_input_process_time=0.0;

/*
 * DivideReads - single-process variant
 * Reads entire file into memory, counts records starting with '>'
 * Returns a malloc'd char* which the caller must free when appropriate.
 *
 * Parameters:
 *   filename  - name of input file
 *   nlines    - output: number of reads (lines starting with '>')
 *   data_size - output: size in bytes of returned buffer
 */
char* DivideReads(const std::string &filename, uint64_t *nlines, size_t *data_size)
{
    using clock = std::chrono::steady_clock;
    auto t_start = clock::now();

    // open file
    FILE *f = fopen(filename.c_str(), "rb");
    if (!f) {
        fprintf(stderr, "Couldn't open input FASTA file: %s\n", filename.c_str());
        return nullptr;
    }
    auto t_open = clock::now();

    // determine file size
    if (fseek(f, 0, SEEK_END) != 0) {
        fprintf(stderr, "fseek failed for file: %s\n", filename.c_str());
        fclose(f);
        return nullptr;
    }
    long filesize = ftell(f);
    if (filesize < 0) {
        fprintf(stderr, "ftell failed for file: %s\n", filename.c_str());
        fclose(f);
        return nullptr;
    }
    rewind(f);

    // allocate buffer (add one extra byte for null-terminator)
    size_t localsize = (size_t)filesize;
    char *data = (char*) malloc(localsize + 1);
    if (!data) {
        fprintf(stderr, "malloc failed for size %zu\n", localsize + 1);
        fclose(f);
        return nullptr;
    }

    // read entire file
    size_t read_total = 0;
    while (read_total < localsize) {
        size_t r = fread(data + read_total, 1, localsize - read_total, f);
        if (r == 0) {
            if (feof(f)) break;
            if (ferror(f)) {
                fprintf(stderr, "Error while reading file: %s\n", filename.c_str());
                free(data);
                fclose(f);
                return nullptr;
            }
        }
        read_total += r;
    }
    data[read_total] = '\0'; // null-terminate
    auto t_read = clock::now();

    // count number of reads (lines starting with '>')
    uint64_t lines = 0;
    for (size_t i = 0; i < read_total; ++i) {
        if (data[i] == '>') ++lines;
    }

    // close file
    fclose(f);
    auto t_close = clock::now();

    // set outputs
    *nlines = lines;
    *data_size = read_total;

    // fill timing globals (single-process; global_* mirror local values)
    file_open_time = std::chrono::duration<double>(t_open - t_start).count();
    file_read_time = std::chrono::duration<double>(t_read - t_open).count();
    file_copy_time = 0.0; // not applicable in this single read
    file_close_time = std::chrono::duration<double>(t_close - t_read).count();
    input_process_time = std::chrono::duration<double>(t_close - t_start).count();

    global_file_open_time = file_open_time;
    global_file_read_time = file_read_time;
    global_file_copy_time = file_copy_time;
    global_file_close_time = file_close_time;
    global_input_process_time = input_process_time;

    // print timings (similar messages as MPI version; single-process)
    fprintf(stderr, "File open time (secs): %f\n", file_open_time);
    fprintf(stderr, "File read time (secs): %f\n", file_read_time);
    fprintf(stderr, "File close time (secs): %f\n", file_close_time);
    fprintf(stderr, "Total input processing time (secs): %f\n", input_process_time);

    return data;
}

// input_read_data perform_input_reading (const std::string &fileName, int read_length)


// {
//     input_read_data input_rdata;
//     input_rdata.read_data = nullptr;
//     input_rdata.read_data_size = 0;
//     input_rdata.start_index = 0;
//     input_rdata.local_count = 0;
//     input_rdata.total = 0;

//     double start_t = 0.0;

//     // read file into memory
//     uint64_t nlines = 0;
//     size_t data_size = 0;

//     char *buffer = DivideReads(fileName, &nlines, &data_size);
//     if (!buffer) {
//         fprintf(stderr, "Failed to read file: %s\n", fileName.c_str());
//         return input_rdata;
//     }

//     // Populate input_rdata fields
//     input_rdata.read_data = buffer;
//     input_rdata.read_data_size = data_size;
//     input_rdata.start_index = 0;     // single-process: start at 0
//     input_rdata.local_count = nlines;
//     input_rdata.total = nlines;

//     fprintf(stderr, "Total number of reads: %lu\n", (unsigned long) nlines);
//     fprintf(stderr, "Completed Reading the input dataset\n");

//     return input_rdata;
// }


input_read_data perform_input_reading (const std::string &fileName, int read_length)
{
    input_read_data input_rdata;
    input_rdata.read_data = nullptr;
    input_rdata.read_data_size = 0;
    input_rdata.start_index = 0;
    input_rdata.local_count = 0;
    input_rdata.total = 0;

    // read file into memory (whole file)
    uint64_t nlines = 0;
    size_t data_size = 0;
    char *buffer = DivideReads(fileName, &nlines, &data_size);
    if (!buffer) {
        fprintf(stderr, "Failed to read file: %s\n", fileName.c_str());
        return input_rdata;
    }

    // --- Build vector of offsets where each read starts (positions of '>') ---
    std::vector<size_t> starts;
    starts.reserve((size_t)nlines);
    for (size_t i = 0; i < data_size; ++i) {
        if (buffer[i] == '>') starts.push_back(i);
    }
    // sanity check
    if (starts.empty()) {
        // nothing to return; keep original buffer ownership to caller (or free and return empty)
        free(buffer);
        fprintf(stderr, "No reads (no '>' found) in file: %s\n", fileName.c_str());
        return input_rdata;
    }

    // --- Partitioning parameters ---
    const int parts = 1;          // divide into 4 parts; change if you want other granularity
    const int selected_part = 0;  // which part to return: 0..parts-1 (0==first quarter)
                                  // change this to 1,2,3 to get other quarters

    // clamp selected_part to valid range
    int sel = selected_part;
    if (sel < 0) sel = 0;
    if (sel >= parts) sel = parts - 1;

    // compute read indices for the selected partition
    uint64_t reads_per_part = (uint64_t) (nlines / parts);
    uint64_t start_read_idx = (uint64_t) sel * reads_per_part;
    uint64_t end_read_idx = (sel == parts - 1) ? nlines : start_read_idx + reads_per_part;
    if (end_read_idx > nlines) end_read_idx = nlines;

    // determine byte offsets in the buffer for copying
    size_t start_byte = starts[(size_t)start_read_idx];
    size_t end_byte;
    if (end_read_idx < starts.size()) {
        end_byte = starts[(size_t)end_read_idx];
    } else {
        // last partition goes to end of file
        end_byte = data_size;
    }

    size_t part_size = (end_byte > start_byte) ? (end_byte - start_byte) : 0;

    // Copy selected partition into a new buffer that we'll return
    char *part_buf = nullptr;
    if (part_size > 0) {
        part_buf = (char*) malloc(part_size + 1); // +1 for null-terminator
        if (!part_buf) {
            fprintf(stderr, "malloc failed for partition size %zu\n", part_size);
            free(buffer);
            return input_rdata;
        }
        // measure copy time if you want (optional)
        memcpy(part_buf, buffer + start_byte, part_size);
        part_buf[part_size] = '\0';
    } else {
        // nothing in this partition (shouldn't usually happen)
        part_buf = (char*) malloc(1);
        part_buf[0] = '\0';
    }

    // free full-file buffer (caller will free the returned partition buffer)
    free(buffer);

    // populate input_rdata with partition info
    input_rdata.read_data = part_buf;
    input_rdata.read_data_size = part_size;
    input_rdata.start_index = 0; // start of this returned buffer
    input_rdata.local_count = (int)(end_read_idx - start_read_idx);
    input_rdata.total = (int)nlines; // total reads in original file

    fprintf(stderr, "Total original reads: %lu\n", (unsigned long)nlines);
    fprintf(stderr, "Returning partition %d/%d : reads %lu .. %lu (count %d)\n",
            sel+1, parts, (unsigned long)start_read_idx, (unsigned long)(end_read_idx-1),
            input_rdata.local_count);
    fprintf(stderr, "Partition byte-range: %zu .. %zu (size %zu)\n", start_byte, end_byte, part_size);
    fprintf(stderr, "Completed Reading the input dataset (partition returned)\n");

    return input_rdata;
}
