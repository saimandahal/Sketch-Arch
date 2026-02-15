#!/bin/bash

set -e

# chmod +x sketch.sh
# ./sketch_gpu_run.sh syncmer

# deafult - minimizer
MODE=${1:-minimizer}

if [[ "$MODE" == "minimizer" ]]; then
    echo "Running sketch_arch "
    ./sketch_arch -s data/Ecoli_sub.fasta -q data/Ecoli_query.fa -l 2000 -w 100 -a A.txt -b B.txt -t 30 -p Prime.txt

elif [[ "$MODE" == "syncmer" ]]; then
    echo "Running sketch_arch"
    ./sketch_arch -s data/Ecoli_sub.fasta -q data/Ecoli_query.fa -l 1000 -w 15 -a A.txt -b B.txt -t 30 -p Prime.txt
else
    echo "Unknown mode: $MODE"
    echo "Usage: $0 [minimizer|syncmer]"
    exit 1
fi
