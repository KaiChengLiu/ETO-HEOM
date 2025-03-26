# Define CUDA 11.8 path
CUDA_HOME = $(HOME)/cuda-11.8

# Compiler and flags
CXX = $(CUDA_HOME)/bin/nvcc
CXXFLAGS = -g -Xcompiler="-Wall -Wno-sign-compare -Wno-unused-variable" -O3 -std=c++14

# Include directories for CUDA and lib folder
INCLUDES = -I/usr/include -I$(CUDA_HOME)/include -I../lib

# Library path and linked libraries
LIB_PATH = ../lib
LDFLAGS = -L$(CUDA_HOME)/lib64 -L$(LIB_PATH) -lgpu2des -lcublas

# Source files in src/
SRCS = GPU_2DES.cu GPU_dynamics.cu

# Object files
OBJS = $(SRCS:.cu=.o)

# Executable names
EXECUTABLES = GPU_2DES GPU_dynamics

# Default target (build all executables)
all: $(EXECUTABLES)

# Compile CUDA object files
%.o: %.cu
	$(CXX) $(CXXFLAGS) $(INCLUDES) -dc $< -o $@

# Link each .o file into its corresponding executable
GPU_2DES: GPU_2DES.o
	$(CXX) $(CXXFLAGS) GPU_2DES.o -o GPU_2DES $(LDFLAGS)

GPU_dynamics: GPU_dynamics.o
	$(CXX) $(CXXFLAGS) GPU_dynamics.o -o GPU_dynamics $(LDFLAGS)

# Clean up build files
clean:
	rm -f $(OBJS) $(EXECUTABLES)

# Phony targets
.PHONY: all clean
