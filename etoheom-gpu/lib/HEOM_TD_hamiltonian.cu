#include "HEOM_TD_hamiltonian.h"
#include "HEOM_kernel_func.h"


void build_V_matrix(cublasHandle_t cublasH, param& key, const data_type* d_X, data_type* d_V, const float t, const int idx) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;

	double E = key.pulses[idx].E0 * exp(-2.77 * (t - key.pulses[idx].tau0) * (t - key.pulses[idx].tau0) / key.pulses[idx].taup / key.pulses[idx].taup);

	// E* exp^ { i * omega * (t - tau0) }

	data_type e = make_cuComplex(cos(key.pulses[idx].w0 * (t - key.pulses[idx].tau0)) * E, sin(key.pulses[idx].w0 * (t - key.pulses[idx].tau0)) * E);
	cudaError(cudaMemcpy(d_V, d_X, total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
	// V = e * X
	cublasError(cublasCscal_v2(cublasH, total_size, &e, d_V, 1), __FILE__, __LINE__);
	//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
}

void build_V_dagger_matrix(cublasHandle_t cublasH, param& key, const data_type* d_X, data_type* d_V, const float t, const int idx) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;

	// E = E0 * exp^{-2 * ln4 * (t - tau0)^2 / taup^2}
	double E = key.pulses[idx].E0 * exp(-2.77 * (t - key.pulses[idx].tau0) * (t - key.pulses[idx].tau0) / key.pulses[idx].taup / key.pulses[idx].taup);

	// e = E * exp^{i * omega * (t - tau0)}
	data_type e = make_cuComplex(cos(key.pulses[idx].w0 * (t - key.pulses[idx].tau0)) * E, -sin(key.pulses[idx].w0 * (t - key.pulses[idx].tau0)) * E);
	cudaError(cudaMemcpy(d_V, d_X, total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
	// V = e * X
	cublasError(cublasCscal_v2(cublasH, total_size, &e, d_V, 1), __FILE__, __LINE__);
	//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
}