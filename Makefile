# Compiler
CC = nvcc
CXXFLAGS = -O3 -std=c++14

# make mode=minimizer

# User can specify mode: minimizer or syncmer
MODE ?= minimizer

# Default flags for each mode
ifeq ($(MODE),minimizer)
	CFLAGS = -O3 -arch=sm_70 -std=c++14 -DWINDW_SIZE=16 -DLMER_LENGTH=8 -Xcompiler -fopenmp -DUSE_CUDA
else ifeq ($(MODE),syncmer)
	CFLAGS = -O3 -arch=sm_70 -std=c++14 -DWINDW_SIZE=11 -DLMER_LENGTH=15 -Xcompiler -fopenmp -DUSE_CUDA
else
	$(error Unknown mode "$(MODE)", must be "minimizer" or "syncmer")
endif

TARGET = sketch_arch

SRC = $(wildcard src/*.cu)
HDR = $(wildcard includes/*.h)

# Build target
$(TARGET): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -Iincludes -o $(TARGET) $(SRC)

# Run target: passes mode as first argument to executable
run: $(TARGET)
	./$(TARGET) $(MODE)

# Clean
clean:
	rm -f $(TARGET)

.PHONY: run clean
