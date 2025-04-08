# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable -O3 -std=c++14 -fopenmp

# Include directories for GSL, other headers, and the library directory
INCLUDES = -I/usr/include -I../lib

# Libraries to link against
LIBS = -L../lib -lgsl -lgslcblas -lm -lcpu2des

# Source files
SRCS = CPU_2DES.cpp CPU_dynamics.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable names
EXEC_2DES = CPU_2DES
EXEC_DYNAMICS = CPU_dynamics

# Default target (build both executables)
all: $(EXEC_2DES) $(EXEC_DYNAMICS)

# Link the CPU_2DES executable
$(EXEC_2DES): CPU_2DES.o
	$(CXX) $(CXXFLAGS) -o $(EXEC_2DES) CPU_2DES.o $(LIBS)

# Link the CPU_dynamics executable
$(EXEC_DYNAMICS): CPU_dynamics.o
	$(CXX) $(CXXFLAGS) -o $(EXEC_DYNAMICS) CPU_dynamics.o $(LIBS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up the build directory
clean:
	rm -f $(OBJS) $(EXEC_2DES) $(EXEC_DYNAMICS)

# Phony targets
.PHONY: all clean
