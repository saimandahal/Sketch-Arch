#!/bin/bash

set -e
# Minimizers
./sketch_gpu -s data/Ecoli_sub.fasta -q data/Ecoli_query.fa -l 2000 -w 100 -a A.txt -b B.txt -t 30 -p Prime.txt

# Syncmers
# ./sketch_gpu -s data/Ecoli_sub.fasta -q data/Ecoli_query.fa -l 1000 -w 15 -a A.txt -b B.txt -t 30 -p Prime.txt
