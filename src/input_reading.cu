
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
#include "../includes/JEM.h"  


double file_open_time=0.0, global_file_open_time=0.0;
double file_read_time=0.0, global_file_read_time=0.0;
double file_copy_time=0.0, global_file_copy_time=0.0;
double file_close_time=0.0, global_file_close_time=0.0;
double input_process_time=0.0, global_input_process_time=0.0;

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

    uint64_t lines = 0;
    for (size_t i = 0; i < read_total; ++i) {
        if (data[i] == '>') ++lines;
    }

    fclose(f);
    auto t_close = clock::now();

    *nlines = lines;
    *data_size = read_total;

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

    fprintf(stderr, "File open time (secs): %f\n", file_open_time);
    fprintf(stderr, "File read time (secs): %f\n", file_read_time);
    fprintf(stderr, "File close time (secs): %f\n", file_close_time);
    fprintf(stderr, "Total input processing time (secs): %f\n", input_process_time);

    return data;
}

input_read_data perform_input_reading (const std::string &fileName, int read_length)
{
    input_read_data input_rdata;
    input_rdata.read_data = nullptr;
    input_rdata.read_data_size = 0;
    input_rdata.start_index = 0;
    input_rdata.local_count = 0;
    input_rdata.total = 0;

    uint64_t nlines = 0;
    size_t data_size = 0;
    char *buffer = DivideReads(fileName, &nlines, &data_size);
    if (!buffer) {
        fprintf(stderr, "Failed to read file: %s\n", fileName.c_str());
        return input_rdata;
    }

    std::vector<size_t> starts;
    starts.reserve((size_t)nlines);
    for (size_t i = 0; i < data_size; ++i) {
        if (buffer[i] == '>') starts.push_back(i);
    }
    // sanity check
    if (starts.empty()) {

        free(buffer);
        fprintf(stderr, "No reads (no '>' found) in file: %s\n", fileName.c_str());
        return input_rdata;
    }

    const int parts = 1;          // divide into n parts
    const int selected_part = 0; 

    int sel = selected_part;
    if (sel < 0) sel = 0;
    if (sel >= parts) sel = parts - 1;

    uint64_t reads_per_part = (uint64_t) (nlines / parts);
    uint64_t start_read_idx = (uint64_t) sel * reads_per_part;
    uint64_t end_read_idx = (sel == parts - 1) ? nlines : start_read_idx + reads_per_part;
    if (end_read_idx > nlines) end_read_idx = nlines;

    size_t start_byte = starts[(size_t)start_read_idx];
    size_t end_byte;
    if (end_read_idx < starts.size()) {
        end_byte = starts[(size_t)end_read_idx];
    } else {
        end_byte = data_size;
    }

    size_t part_size = (end_byte > start_byte) ? (end_byte - start_byte) : 0;

    char *part_buf = nullptr;
    if (part_size > 0) {
        part_buf = (char*) malloc(part_size + 1); // +1 for null-terminator
        if (!part_buf) {
            fprintf(stderr, "malloc failed for partition size %zu\n", part_size);
            free(buffer);
            return input_rdata;
        }
        // measure copy time
        memcpy(part_buf, buffer + start_byte, part_size);
        part_buf[part_size] = '\0';
    } else {
        part_buf = (char*) malloc(1);
        part_buf[0] = '\0';
    }

    free(buffer);

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
