# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable -O3 -std=c++14 -fopenmp

# Include directories for GSL and other headers
INCLUDES = -I/usr/include

# Source files
SRCS = constant.cpp dynamics.cpp polar.cpp TD_hamiltonian.cpp param.cpp utilize.cpp

# Header files
HDRS = constant.h dynamics.h polar.h TD_hamiltonian.h param.h utilize.h

# Object files
OBJS = $(SRCS:.cpp=.o)

# Static library name
LIB = libcpu2des.a

# Default target
all: $(LIB)

# Compile object files
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create static library
$(LIB): $(OBJS)
	ar rcs $(LIB) $(OBJS)

# Clean up the build directory
clean:
	rm -f $(OBJS) $(LIB) *.a

# Phony targets
.PHONY: all clean