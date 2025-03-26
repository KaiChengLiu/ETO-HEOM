#include "HEOM_utilize.h"
#include "HEOM_kernel_func.h"
#include "HEOM_constant.h"
#include "HEOM_dynamics.h"
#include "HEOM_TD_hamiltonian.h"
#include "HEOM_param.h"
#include "HEOM_polar.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>


using namespace std;

int main(int argc, char** argv)
{
	string filename(argv[1]);
	param k;
	k.param_2des(filename);
	construct_ADO_set(k);

	cublasHandle_t cublasH;
	cublasCreate(&cublasH);

	vector<vector<data_type>> p;
	for (int i = 0; i < k.n_sample; i++) {
		//cout << "Now running sample " << i + 1 << '\n';
		//cout << "The program is running with " << k.ado.size() << " ADOs" << '\n';
		//cout << "The coherent time is " << k.pulses[1].tau0 - k.pulses[0].tau0 << " fs " << "and the population time is ";
		//if (k.pulses[1].tau0 >= k.pulses[0].tau0) cout << k.pulses[2].tau0 - k.pulses[1].tau0 << " fs" << '\n';
		//else cout << k.pulses[2].tau0 - k.pulses[0].tau0 << " fs" << '\n';

		vector<data_type> H(k.sys_size * k.sys_size);
		construct_Hal(k, H);
		//cout << "The disordered Hamiltonain is:" << '\n';
		//print_Hal(H, sys_size);
		//cout << '\n';
		data_type* d_H = To_Device(H);
		polar_mat_set(k);
		polar_mat_ranrot(k);


		vector<data_type> p1;
		twoD_spectrum_solver(k, d_H, 1, 1, 1, p1);

		vector<data_type> p2;
		twoD_spectrum_solver(k, d_H, 1, 1, 0, p2);

		vector<data_type> p3;
		twoD_spectrum_solver(k, d_H, 1, 0, 1, p3);

		vector<data_type> p_i;
		for (int j = 0; j < p1.size(); j++) {
			float real = p1[j].x - p2[j].x - p3[j].x;
			float imag = p1[j].y - p2[j].y - p3[j].y;
			p_i.push_back(make_cuComplex(real, imag));
		}
		p.push_back(p_i);

		cudaError(cudaFree(d_H), __FILE__, __LINE__);
	}

	vector<data_type> P(p[0].size(), make_cuComplex(0.0, 0.0));
	for (int i = 0; i < p.size(); i++) {
		for (int j = 0; j < p[i].size(); j++) {
			P[j].x += p[i][j].x / k.n_sample;
			P[j].y += p[i][j].y / k.n_sample;
		}
	}

	for (int i = 0; i < P.size(); i++) cout << P[i].x << " " << P[i].y << '\n';

	cublasError(cublasDestroy(cublasH), __FILE__, __LINE__);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

	k.param_free();
	return 0;
}



