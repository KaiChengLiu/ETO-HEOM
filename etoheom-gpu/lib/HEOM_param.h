#ifndef HEOM_PARAM_H
#define HEOM_PARAM_H

#include "cuda_runtime.h"
#include "cusolverDn.h"
#include "cublas_v2.h"
#include "device_launch_parameters.h"
#include "HEOM_kernel_func.h"
#include "HEOM_constant.h"
#include "HEOM_utilize.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>


using namespace std;

typedef struct pulse pulse;

struct pulse
{
	float E0; // pulse amplitude
	float tau0; // pulse central time
	float taup; // pulse width
	float w0; // pulse frequency
};

class param {

public:
	void param_2des(string filename);
	void param_dynamics(string filename);
	int sys_size = 0;
	int single_size = 0;
	vector<string> ado;
	unordered_map<string, int> ado_map;
	unordered_map<int, pair<int, int>> d_map;
	data_type* d_rho = nullptr;
	vector<data_type*> d_rho_ptr;
	data_type** d_rho_batch = nullptr;
	data_type* d_Hal = nullptr;
	data_type* d_Disorder = nullptr;
	data_type* d_S = nullptr;
	data_type* d_Xx = nullptr;
	data_type* d_Xy = nullptr;
	data_type* d_Xz = nullptr;
	data_type* d_X_site = nullptr;
	data_type* d_polar = nullptr;
	data_type* d_eigenstate = nullptr;
	int K = 0;
	int K_m = 0;
	int K_extra = 0;
	int L = 0;
	int n_sample = 0;
	data_type beta;
	vector<data_type> Hal;
	vector<data_type> Disorder;
	string bath_type;

	vector<float> angle;
	vector<vector<data_type>> alpha;
	vector<vector<data_type>> alpha_t;
	vector<vector<data_type>> gamma;
	vector<data_type> lambda;
	pulse pulses[3];
	float t_start = 0;
	float t_end = 0;
	float step_size = 0;
	float print_step = 0;

	void reset_rho();
	void reset_rho_batch();
	void param_free();
};

#endif
