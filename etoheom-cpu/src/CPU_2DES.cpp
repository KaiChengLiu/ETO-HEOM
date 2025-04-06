#include <iostream>
#include <iomanip>
#include <time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "param.h"
#include "utilize.h"
#include "dynamics.h"
#include "polar.h"

using namespace std;

int main(int argc, char** argv)
{
	string filename(argv[1]);
	param k;
	k.param_2des(filename);


	construct_ADO_set(k);
  
	const gsl_rng_type* T;
	gsl_rng* r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, k.seed);

	//cout << "The program is running with importance value cutoff level " << k.cutoff_level << " and " << k.ado.size() << " ADOs" << '\n';
	//cout << "The coherent time is " << k.pulses[1].tau0 - k.pulses[0].tau0 << " fs " << "and the population time is ";
	//if (k.pulses[1].tau0 >= k.pulses[0].tau0) cout << k.pulses[2].tau0 - k.pulses[1].tau0 << " fs" << '\n';
	//else cout << k.pulses[2].tau0 - k.pulses[0].tau0 << " fs" << '\n';


	int sys_size = k.sys_size;
	int total_size = k.sys_size * k.sys_size;
	gsl_matrix_complex* H = gsl_matrix_complex_alloc(sys_size, sys_size);

	vector<vector<gsl_complex>> p;

	vector<gsl_matrix_complex*> OP(1);
	for (int i = 0; i < 1; i++) {
		OP[i] = gsl_matrix_complex_alloc(k.sys_size, k.sys_size);
		gsl_matrix_complex_set_all(OP[i], gsl_complex_rect(0.0, 0.0));
	}

	int data_size = 0;
	for (int i = 0; i < k.n_sample; i++) {
		//cout << "#Now running sample " << i + 1 << '\n';

		gsl_matrix_complex_set_all(H, gsl_complex_rect(0.0, 0.0));
		construct_Hal(k, H, r);
		//cout << "The disordered Hamiltonain is:" << '\n';
		//print_Hal(H);

		//cout << "The absolute dipole operator X_abs is:" << '\n';
		polar_mat_set(k);
		polar_mat_ranrot(k, r);
		compute_pulse_interaction(k);
		//print_matrix_real(k.X_site);

		vector<gsl_complex> p1;
		twoD_spectrum_solver(k, H, 1, 1, 1, p1);

		vector<gsl_complex> p2;
		twoD_spectrum_solver(k, H, 1, 1, 0, p2);

		vector<gsl_complex> p3;
		twoD_spectrum_solver(k, H, 1, 0, 1, p3);

		vector<gsl_complex> pi;
		for (int j = 0; j < p1.size(); j++) {
			double real = GSL_REAL(p1[j]) - GSL_REAL(p2[j]) - GSL_REAL(p3[j]);
			double imag = GSL_IMAG(p1[j]) - GSL_IMAG(p2[j]) - GSL_IMAG(p3[j]);
			pi.push_back(gsl_complex_rect(real, imag));
		}
		p.push_back(pi);
	}
	cblas_zscal(total_size, &ONE, gsl_matrix_complex_ptr(OP[0], 0, 0), 1);

	vector<gsl_complex> P(p[0].size(), gsl_complex_rect(0.0, 0.0));
	for (int i = 0; i < p.size(); i++) {
		for (int j = 0; j < p[i].size(); j++) {
		P[j].dat[0] += GSL_REAL(p[i][j]) / k.n_sample;
		P[j].dat[1] += GSL_IMAG(p[i][j]) / k.n_sample;
		}
	}

	for (int i = 0; i < P.size(); i++) cout << scientific << setprecision(6) << GSL_REAL(P[i]) << " " << GSL_IMAG(P[i]) << '\n';

	for(int i = 0; i < 1; i++) gsl_matrix_complex_free(OP[i]);
	k.param_free();
	return 0;
}


