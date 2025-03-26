# CUDA 11.8 installation path
CUDA_HOME = $(HOME)/cuda-11.8

# Compiler and flags
CXX = $(CUDA_HOME)/bin/nvcc
CXXFLAGS = -g -Xcompiler="-Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable" -O3 -std=c++14

# Include directories for CUDA and other headers
INCLUDES = -I/usr/include -I$(CUDA_HOME)/include

# Library directories
LDFLAGS = -L$(CUDA_HOME)/lib64 -lcublas

# Source files
SRCS = HEOM_constant.cu HEOM_dynamics.cu HEOM_polar.cu HEOM_TD_hamiltonian.cu HEOM_param.cu HEOM_utilize.cu HEOM_kernel_func.cu

# Header files
HDRS = HEOM_constant.h HEOM_dynamics.h HEOM_polar.h HEOM_TD_hamiltonian.h HEOM_param.h HEOM_utilize.h HEOM_kernel_func.h

# Object files
OBJS = $(SRCS:.cu=.o)

# Static library name
LIB = libgpu2des.a

# Default target
all: $(LIB)

# Compile object files
%.o: %.cu $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -dc $< -o $@

# Create static library
$(LIB): $(OBJS)
	ar rcs $(LIB) $(OBJS)

# Clean up the build directory
clean:
	rm -f $(OBJS) $(LIB) *.a

# Phony targets
.PHONY: all clean
