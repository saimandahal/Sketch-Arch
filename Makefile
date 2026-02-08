CC := nvcc
CXXFLAGS := -O3 -std=c++14

# default values for device and mode
device ?= gpu        
mode   ?= minimizer  

COMMON_FLAGS := -O3 -std=c++14 -Xcompiler -fopenmp

# mode
ifeq ($(mode),minimizer)
    WINDOW_SIZE := 16
    LMER_LENGTH := 8
else ifeq ($(mode),syncmer)
    WINDOW_SIZE := 11
    LMER_LENGTH := 15
else
    $(error Invalid mode: $(mode))
endif

# device
ifeq ($(device),gpu)
    DEVICE_FLAGS := -arch=sm_70 -DUSE_CUDA
else ifeq ($(device),cpu)
    DEVICE_FLAGS :=
else
    $(error Invalid device: $(device))
endif

CFLAGS := $(COMMON_FLAGS) \
          -DWINDW_SIZE=$(WINDOW_SIZE) \
          -DLMER_LENGTH=$(LMER_LENGTH) \
          $(DEVICE_FLAGS)



TARGET = sketch_arch

SRC = $(wildcard src/*.cu)
HDR = $(wildcard includes/*.h)

$(TARGET): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -Iincludes -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET) -s subject.fasta -q query.fasta -l 100

clean:
	rm -f $(TARGET)

.PHONY: run clean
