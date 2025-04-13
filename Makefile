# === Global compiler settings ===
CXX := g++
NVCC := $(HOME)/cuda-11.8/bin/nvcc
CXXFLAGS := -Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable -O3 -std=c++14
NVCCFLAGS := -g -Xcompiler="-Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable" -O3 -std=c++14

INCLUDES := -I/usr/include
CUDA_INC := -I$(HOME)/cuda-11.8/include
CUDA_LIB := -L$(HOME)/cuda-11.8/lib64

# === Directory settings ===
CPU_LIB_DIR := etoheom-cpu/lib
CPU_SRC_DIR := etoheom-cpu/src
GPU_LIB_DIR := etoheom-gpu/lib
GPU_SRC_DIR := etoheom-gpu/src

# === File lists ===

## CPU
CPU_LIB_SRCS := $(wildcard $(CPU_LIB_DIR)/*.cpp)
CPU_LIB_OBJS := $(CPU_LIB_SRCS:.cpp=.o)
CPU_LIB := $(CPU_LIB_DIR)/libcpu2des.a

CPU_SRC_SRCS := $(wildcard $(CPU_SRC_DIR)/*.cpp)
CPU_SRC_OBJS := $(CPU_SRC_SRCS:.cpp=.o)
CPU_BIN := $(CPU_SRC_DIR)/CPU_2DES $(CPU_SRC_DIR)/CPU_dynamics

## GPU
GPU_LIB_SRCS := $(wildcard $(GPU_LIB_DIR)/*.cu)
GPU_LIB_OBJS := $(GPU_LIB_SRCS:.cu=.o)
GPU_LIB := $(GPU_LIB_DIR)/libgpu2des.a

GPU_SRC_SRCS := $(wildcard $(GPU_SRC_DIR)/*.cu)
GPU_SRC_OBJS := $(GPU_SRC_SRCS:.cu=.o)
GPU_BIN := $(GPU_SRC_DIR)/GPU_2DES $(GPU_SRC_DIR)/GPU_dynamics

# === Default target ===
all: cpu gpu

# === CPU targets ===
cpu: $(CPU_LIB) $(CPU_BIN)

$(CPU_LIB): $(CPU_LIB_OBJS)
	ar rcs $@ $^

$(CPU_LIB_DIR)/%.o: $(CPU_LIB_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(CPU_SRC_DIR)/CPU_2DES: $(CPU_SRC_DIR)/CPU_2DES.o $(CPU_LIB)
	$(CXX) $(CXXFLAGS) -I$(CPU_LIB_DIR) -o $@ $< -L$(CPU_LIB_DIR) -lcpu2des -lgsl -lgslcblas -lm

$(CPU_SRC_DIR)/CPU_dynamics: $(CPU_SRC_DIR)/CPU_dynamics.o $(CPU_LIB)
	$(CXX) $(CXXFLAGS) -I$(CPU_LIB_DIR) -o $@ $< -L$(CPU_LIB_DIR) -lcpu2des -lgsl -lgslcblas -lm

$(CPU_SRC_DIR)/%.o: $(CPU_SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -I$(CPU_LIB_DIR) -c $< -o $@

# === GPU targets ===
gpu: $(GPU_LIB) $(GPU_BIN)

$(GPU_LIB): $(GPU_LIB_OBJS)
	ar rcs $@ $^

$(GPU_LIB_DIR)/%.o: $(GPU_LIB_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $(CUDA_INC) $(INCLUDES) -dc $< -o $@

$(GPU_SRC_DIR)/GPU_2DES: $(GPU_SRC_DIR)/GPU_2DES.o $(GPU_LIB)
	$(NVCC) $(NVCCFLAGS) $(CUDA_INC) -o $@ $< -L$(GPU_LIB_DIR) -lgpu2des $(CUDA_LIB) -lcublas

$(GPU_SRC_DIR)/GPU_dynamics: $(GPU_SRC_DIR)/GPU_dynamics.o $(GPU_LIB)
	$(NVCC) $(NVCCFLAGS) $(CUDA_INC) -o $@ $< -L$(GPU_LIB_DIR) -lgpu2des $(CUDA_LIB) -lcublas

$(GPU_SRC_DIR)/%.o: $(GPU_SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $(CUDA_INC) -I$(GPU_LIB_DIR) -dc $< -o $@

# === Clean ===
clean:
	rm -f $(CPU_LIB_OBJS) $(CPU_SRC_OBJS) $(CPU_BIN) $(CPU_LIB)
	rm -f $(GPU_LIB_OBJS) $(GPU_SRC_OBJS) $(GPU_BIN) $(GPU_LIB)

.PHONY: all clean cpu gpu
