#ifndef HEOM_UTILIZE_H
#define HEOM_UTILIZE_H

#include "HEOM_kernel_func.h"
#include "HEOM_constant.h"
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include "cublas_v2.h"
#include "cusolverDn.h"
#include "cuComplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
#include "HEOM_param.h"

using namespace std;

#define cudaError(ans, file, line) { gpuAssert((ans), (file), (line)); }
#define cublasError(ans, file, line) { cublasAssert((ans), (file), (line)); }

inline void gpuAssert(cudaError_t code, const char* file, int line) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		exit(EXIT_FAILURE);
	}
}

inline void cublasAssert(cublasStatus_t code, const char* file, int line) {
	if (code != CUBLAS_STATUS_SUCCESS) {
		fprintf(stderr, "cuBLASassert: %d %s %d\n", code, file, line);
		exit(EXIT_FAILURE);
	}
}

//extern vector<string> key;

double generateRandomGaussian(double mean, double stddev);

data_type* To_Device(const vector<data_type>& A);

data_type* To_Device(const data_type& A);

void To_Host(const data_type* d_A, vector<data_type>& A);

void To_Host(const data_type* d_A, data_type& A, const int offset);


void cusolverError(cusolverStatus_t status, string file, int line);


void print_matrix(const vector<data_type> m, const int size);
void print_matrix_real(const std::vector<data_type> m, const int size);
void print_matrix(const data_type* m, const int size);
void print_Hal(const vector<data_type> m, const int size);

float cuComplex_abs(const data_type c);

#endif