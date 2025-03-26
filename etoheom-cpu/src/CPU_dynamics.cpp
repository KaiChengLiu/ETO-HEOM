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
	cout << "Now is CPU version" << '\n';
	cout << "Running " << filename << '\n';
	param k;
	k.param_dynamics(filename);
	cout << "The bath type is " << k.bath_type << '\n';
	construct_ADO_set(k);

	int sys_size = k.sys_size;
	int total_size = k.sys_size * k.sys_size;

	vector<int> sites = {1, 2};
	vector<vector<double>> population(sites.size(), vector<double>((k.t_end - k.t_start) / k.print_step + 1, 0.0));

	dynamics_solver(k, k.Hal, sites, population);

    vector<gsl_matrix_complex*> OP(1);
	for (int i = 0; i < OP.size(); i++) {
		OP[i] = gsl_matrix_complex_alloc(k.sys_size, k.sys_size);
		gsl_matrix_complex_set_all(OP[i], gsl_complex_rect(0.0, 0.0));
	}
    cblas_zscal(total_size, &ONE, gsl_matrix_complex_ptr(OP[0], 0, 0), 1);

	for (int i = 0; i < sites.size(); i++) {
        cout << "site" << " " << i + 1 << " " << "population: " << '\n' << '\n';
		for (int j = 0; j < population[i].size(); j++) {
			cout << population[i][j] << '\n';
		}
        cout << '\n';
	}

	k.param_free();
	return 0;
}


