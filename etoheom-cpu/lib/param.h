#ifndef PARAM_H
#define PARAM_H

#include "constant.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <utility>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



using namespace std;

typedef struct pulse pulse;

struct pulse
{
	double E0 = 0; // pulse amplitude
	double tau0 = 0; // pulse central time
	double taup = 0; // pulse width
	double w0 = 0; // pulse frequency
};

class param {

public:
	void param_2des(string filename);
	void param_dynamics(string filename);
	int single_size = 0;
	int sys_size = 0;
	unordered_map<int, pair<int, int>> d_map;
	vector<string> ado;
	unordered_map<string, int> ado_map;
	vector<gsl_matrix_complex*> rho;
	gsl_matrix_complex* Hal = nullptr;
	gsl_matrix_complex* Disorder = nullptr;
	vector<gsl_matrix_complex*> S;
	gsl_matrix_complex* X_site = nullptr;
	gsl_matrix_complex* X_x = nullptr;
	gsl_matrix_complex* X_y = nullptr;
	gsl_matrix_complex* X_z = nullptr;
	vector<gsl_matrix_complex*> X;
	gsl_matrix_complex* polar = nullptr;
	int K = 0;
	int K_m = 0;
	int L = 0;

	int n_sample = 0;
	unsigned long seed = 0;
	double beta = 0;
	string bath_type;
	vector<double> coupling_str;
	vector<double> wc;
	vector<double> angle;
	//vector<vector<double>> alpha;
	vector<vector<gsl_complex>> alpha;
	vector<vector<gsl_complex>> alpha_t;
	vector<vector<gsl_complex>> gamma;
	vector<double> lambda;
	pulse pulses[3];
	float t_start = 0;
	float t_end = 0;
	float step_size = 0;
	float print_step = 0;

	void reset_rho();
	void param_free();
};

#endif
