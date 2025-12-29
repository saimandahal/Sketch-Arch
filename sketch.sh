#!/bin/bash

set -e


# Long reads

# Minimizers
./sketch_gpu -s data/Ecoli_sub.fasta -q data/Ecoli_query.fa -l 2000 -w 100 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/PA_subject.fasta -q ../Sketch/data/PA_query.fa -l 2000 -w 100 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/Human8_subject.fasta -q ../Sketch/data/Human8_query.fa -l 2000 -w 100 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/Sativa_subject.fasta -q ../Sketch/data/Sativa_query.fa -l 2000 -w 100 -a A.txt -b B.txt -t 30 -p Prime.txt


# Short reads

# ./sketch_gpu -s ../Sketch/data/Ecoli_short_read.fa -q ../Sketch/data/Ecoli_short_read.fa -l 80 -w 40 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/C_elegans_pg_sr_100x.fa -q ../Sketch/data/C_elegans_pg_sr_100x.fa -l 80 -w 40 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/HM7_short_reads.fa -q ../Sketch/data/HM7_short_reads.fa -l 80 -w 40 -a A.txt -b B.txt -t 30 -p Prime.txt


# Syncmers

# ./sketch_gpu -s data/Ecoli_sub.fasta -q data/Ecoli_query.fa -l 1000 -w 15 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/PA_subject.fasta -q ..1/Sketch/data/PA_query.fa -l 1000 -w 15 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/Human8_subject.fasta -q ../Sketch/data/Human8_query.fa -l 1000 -w 15 -a A.txt -b B.txt -t 30 -p Prime.txt
# ./sketch_gpu -s ../Sketch/data/Sativa_subject.fasta -q ../Sketch/data/Sativa_query.fa -l 1000 -w 15 -a A.txt -b B.txt -t 30 -p Prime.txt
