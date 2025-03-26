#include "HEOM_utilize.h"


/*
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
*/

double generateRandomGaussian(double mean, double stddev) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<double> dist(mean, stddev);

	return dist(gen);
}

data_type* To_Device(const std::vector<data_type>& A) {
	data_type* d_A = nullptr;
	cudaError(cudaMalloc((void**)&d_A, A.size() * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_A, A.data(), A.size() * sizeof(data_type), cudaMemcpyHostToDevice), __FILE__, __LINE__);
	return d_A;
}

data_type* To_Device(const data_type& A) {
	data_type* d_A = nullptr;
	cudaError(cudaMalloc((void**)&d_A, sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_A, &A, sizeof(data_type), cudaMemcpyHostToDevice), __FILE__, __LINE__);

	return d_A;
}

void To_Host(const data_type* d_A, std::vector<data_type>& A){
	cudaError(cudaMemcpy(A.data(), d_A, sizeof(data_type)*A.size(), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
}

void To_Host(const data_type* d_A, data_type& A, const int offset) {
	cudaError(cudaMemcpy(&A, d_A + offset, sizeof(data_type), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
}


void cusolverError(cusolverStatus_t status, std::string file, int line) {
	if (status != CUSOLVER_STATUS_SUCCESS) {
		printf("cusolver error: %d, at file %s line %d\n", status, file.c_str(), line);
		exit(EXIT_FAILURE);
	}
}

void print_matrix(const std::vector<data_type> m, const int size) {
	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			printf("%.6f + %.6fi   ", m[j * size + i].x, m[j * size + i].y);
		}
		puts("");
	}
}

void print_matrix_real(const std::vector<data_type> m, const int size) {
	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			printf("%.3f\t", m[j * size + i].x);
		}
		puts("");
	}
}

void print_matrix(const data_type* m, const int size) {
	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			printf("%.6f + %.6fi   ", m[i + j * size].x, m[i + j * size].y);
		}
		puts("");
	}
}



void print_Hal(const vector<data_type> m, const int size) {
	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			printf("%.f\t", m[i + j * size].x * 5308);
		}
		puts("");
	}

}

float cuComplex_abs(const data_type c) {
	return sqrt(c.x * c.x + c.y * c.y);
}
