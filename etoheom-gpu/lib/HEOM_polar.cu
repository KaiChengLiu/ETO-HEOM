#include "HEOM_polar.h"

/*The array size is fixed to 3*/
void vector3_norm(float* vec) {
	float vabs;
	vabs = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	vec[0] = vec[0] / vabs;
	vec[1] = vec[1] / vabs;
	vec[2] = vec[2] / vabs;
}

/*Polar is a 3 * 4 matrix*/
void polar_mat_set(param& key) {
	vector<data_type> polar(12, make_cuComplex(0.0, 0.0));
	if (key.d_polar) cudaError(cudaFree(key.d_polar), __FILE__, __LINE__);

	for (int i = 0; i < 4; i++) {
		polar[i * 3 + 0] = make_cuComplex(0.0, 0.0);
		polar[i * 3 + 1] = make_cuComplex(sin(key.angle[i]), 0.0);
		polar[i * 3 + 2] = make_cuComplex(cos(key.angle[i]), 0.0);
	}
	key.d_polar = To_Device(polar);
}

void polar_mat_ranrot(param& key) {
	float a[3];
	float b[3];
	float c[3];

	vector<data_type> rotmat(9, make_cuComplex(0.0, 0.0));
	vector<data_type> tmp(12, make_cuComplex(0.0, 0.0));

	for (int i = 0; i < 3; i++) a[i] = (float)generateRandomGaussian(0, 10);
	for (int i = 0; i < 3; i++) c[i] = (float)generateRandomGaussian(0, 10);
	vector3_norm(c);
	float product = a[0] * c[0] + a[1] * c[1] + a[2] * c[2];
	a[0] = a[0] - product * c[0];
	a[1] = a[1] - product * c[1];
	a[2] = a[2] - product * c[2];
	vector3_norm(a);

	b[0] = c[1] * a[2] - c[2] * a[1];
	b[1] = c[2] * a[0] - c[0] * a[2];
	b[2] = c[0] * a[1] - c[1] * a[0];

	for (int i = 0; i < 3; i++) {
		rotmat[i * 3 + 0] = make_cuComplex(a[i], 0.0);
		rotmat[i * 3 + 1] = make_cuComplex(b[i], 0.0);
		rotmat[i * 3 + 2] = make_cuComplex(c[i], 0.0);
	}


	cuComplex* d_rotmat = To_Device(rotmat);
	cuComplex* d_tmp = To_Device(tmp);
	cublasHandle_t cublasH;
	cublasCreate(&cublasH);
	cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
		3, 4, 3,
		&ONE, d_rotmat, 3, key.d_polar, 3,
		&ZERO, d_tmp, 3), __FILE__, __LINE__);
	cudaError(cudaMemcpy(key.d_polar, d_tmp, 12 * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);

	cublasDestroy(cublasH);
	cudaError(cudaFree(d_rotmat), __FILE__, __LINE__);
	cudaError(cudaFree(d_tmp), __FILE__, __LINE__);
}


void compute_pulse_interaction(param& key, vector<data_type*>& M) {
	if (M.size() != 4) {
		cout << "The number of pulse polarization should be 4 instead of " << M.size() << '\n';
		exit(EXIT_FAILURE);
	}


	cublasHandle_t cublasH;
	cublasCreate(&cublasH);

	int size = key.sys_size;
	data_type E;
	vector<data_type> tmp(size * size, make_cuComplex(0.0, 0.0));
	data_type* d_tmp = To_Device(tmp);
	data_type* d_X[3];
	for (int i = 0; i < 3; i++) d_X[i] = To_Device(tmp);
	cudaError(cudaMemcpy(d_X[0], key.d_Xx, size * size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_X[1], key.d_Xy, size * size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_X[2], key.d_Xz, size * size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);

	vector<data_type> polar(12);
	To_Host(key.d_polar, polar);
	for (int i = 0; i < 4; i++) {
		cudaError(cudaMemset(M[i], 0.0, size * size * sizeof(data_type)), __FILE__, __LINE__);
		for (int j = 0; j < 3; j++) {
			E = polar[i * 3 + j];
			cudaError(cudaMemcpy(d_tmp, d_X[j], size * size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
			cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
			cublasError(cublasCaxpy_v2(cublasH, size * size, &E, d_tmp, 1, M[i], 1), __FILE__, __LINE__);
			cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

		}
	}

	cublasDestroy(cublasH);
	cudaError(cudaFree(d_tmp), __FILE__, __LINE__);
	for (int i = 0; i < 3; i++) cudaError(cudaFree(d_X[i]), __FILE__, __LINE__);
}
