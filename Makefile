

# CC = nvcc
# CXX = g++
# CFLAGS = -O3 -arch=sm_70 -std=c++14 -DLMER_LENGTH=8 
# CFLAGS += -DWINDW_SIZE=31
# CXXFLAGS = -O3 -std=c++14
# TARGET = sketch_gpu
# SOURCES = src/main.cu src/input_reading.cu

# $(TARGET): $(SOURCES)
# 	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCES)

# run: $(TARGET)
# 	./$(TARGET) -s subject.fasta -q query.fasta -l 100

# clean:
# 	rm -f $(TARGET)

# .PHONY: run clean


CC = nvcc
CXXFLAGS = -O3 -std=c++14
# CFLAGS = -O3 -arch=sm_70 -std=c++14 -DWINDW_SIZE=16 -DLMER_LENGTH=8 -Xcompiler -fopenmp

# CFLAGS = -O3 -arch=sm_70 -std=c++14 -DWINDW_SIZE=16 -DLMER_LENGTH=8 -Xcompiler -fopenmp -DUSE_CUDA

# Syncmers
CFLAGS = -O3 -arch=sm_70 -std=c++14 -DWINDW_SIZE=11 -DLMER_LENGTH=15 -Xcompiler -fopenmp -DUSE_CUDA
# CFLAGS = -O3 -arch=sm_70 -std=c++14 -DWINDW_SIZE=11 -DLMER_LENGTH=8 -Xcompiler -fopenmp 


TARGET = sketch_gpu

SRC = $(wildcard src/*.cu)
HDR = $(wildcard includes/*.h)

$(TARGET): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -Iincludes -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET) -s subject.fasta -q query.fasta -l 100

clean:
	rm -f $(TARGET)

.PHONY: run clean
