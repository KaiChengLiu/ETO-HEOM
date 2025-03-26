#include "utilize.h"

void Build_ADO(param& key, std::string& current, const int cur_L) {
	if (current.length() == key.K * key.K_m) {
		if (cur_L <= key.L) {
			key.ado.push_back(current);
		}
		return;
	}

	for (char c = 0; c <= key.L; c++) {
		if (cur_L + c <= key.L) {
			current.push_back(c);
			Build_ADO(key, current, cur_L + c);
			current.pop_back();
		}
	}
}

void Build_ADO_map(param& key) {
	for (int i = 0; i < key.ado.size(); i++) {
		key.ado_map[key.ado[i]] = i;
	}
}

void construct_ADO_set(param& key) {
	string s1;
	Build_ADO(key, s1, 0);
	Build_ADO_map(key);
	for (int i = 0; i < key.ado.size(); i++) {
		gsl_matrix_complex* m = gsl_matrix_complex_alloc(key.sys_size, key.sys_size);
		gsl_matrix_complex_set_all(m, gsl_complex_rect(0.0, 0.0));
		key.rho.push_back(m);
	}
	gsl_matrix_complex_set(key.rho[0], 0, 0, gsl_complex_rect(1.0, 0.0));
}

void print_matrix_real(const gsl_matrix_complex* M) {
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			gsl_complex ele = gsl_matrix_complex_get(M, i, j);
			printf("%.4f\t", GSL_REAL(ele));
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void print_Hal(const gsl_matrix_complex* M) {
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			gsl_complex ele = gsl_matrix_complex_get(M, i, j);
			printf("%.f\t", GSL_REAL(ele) * 5308);
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void print_matrix(const gsl_matrix* M) {
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			double ele = gsl_matrix_get(M, i, j);
			std::cout << ele << " ";
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

double delta(int a, int b) {
	return a == b ? 1 : 0;
}
