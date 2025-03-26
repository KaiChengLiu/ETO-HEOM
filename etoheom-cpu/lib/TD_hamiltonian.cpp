#include "TD_hamiltonian.h"



void build_V_matrix(param& key, const gsl_matrix_complex* X, gsl_matrix_complex* V, const float t, const int idx) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;

	double E = key.pulses[idx].E0 * exp(-2.77 * (t - key.pulses[idx].tau0) * (t - key.pulses[idx].tau0) / key.pulses[idx].taup / key.pulses[idx].taup);

	// E* exp^ { i * omega * (t - tau0) }
	gsl_complex e = gsl_complex_polar(E, key.pulses[idx].w0 * (t - key.pulses[idx].tau0));
	// V = e * X
	gsl_matrix_complex_memcpy(V, X);
	gsl_matrix_complex_scale(V, e);
	
}

void build_V_dagger_matrix(param& key, const gsl_matrix_complex* X, gsl_matrix_complex* V, const float t, const int idx) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;

	// E = E0 * exp^{-2 * ln4 * (t - tau0)^2 / taup^2}
	double E = key.pulses[idx].E0 * exp(-2.77 * (t - key.pulses[idx].tau0) * (t - key.pulses[idx].tau0) / key.pulses[idx].taup / key.pulses[idx].taup);

	// E* exp^ { -i * omega * (t - tau0) }
	gsl_complex e = gsl_complex_polar(E, -key.pulses[idx].w0 * (t - key.pulses[idx].tau0));
	// V = e * X
	gsl_matrix_complex_memcpy(V, X);
	gsl_matrix_complex_scale(V, e);
}