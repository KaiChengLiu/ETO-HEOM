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
	k.param_dynamics(filename);
	construct_ADO_set(k);
	vector<int> sites = { 1 };

	vector<vector<float>> population(sites.size());
	
	dynamics_solver(k, sites, population);
	
	for (int i = 0; i < sites.size(); i++) {
		cout << "site" << " " << i + 1 << " " << "population: " << '\n' << '\n';
		for (int j = 0; j < population[0].size(); j++) {
			cout << population[i][j] << '\n';
		}
		cout << '\n';
	}
	
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

	k.param_free();
	return 0;
}



