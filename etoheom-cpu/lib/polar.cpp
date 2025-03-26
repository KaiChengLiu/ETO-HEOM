#include "polar.h"

void vector3_norm(double* vec) {
	double vabs;
	vabs = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	vec[0] = vec[0] / vabs;
	vec[1] = vec[1] / vabs;
	vec[2] = vec[2] / vabs;
}

void vector3_add(double* vec1, const double* vec2) {
	vec1[0] += vec2[0];
	vec1[1] += vec2[1];
	vec1[2] += vec2[2];
}

void polar_mat_set(param& key) {
	if (key.polar->size1 != 3 || key.polar->size2 != 4) {
		cout << "Polar matrix size should be 3 * 4" << '\n';
		exit(EXIT_FAILURE);
	}
	

	for (int i = 0; i < 4; i++) {
		gsl_matrix_complex_set(key.polar, 0, i, gsl_complex_rect(0.0, 0.0));
		gsl_matrix_complex_set(key.polar, 1, i, gsl_complex_rect(sin(key.angle[i]), 0.0));
		gsl_matrix_complex_set(key.polar, 2, i, gsl_complex_rect(cos(key.angle[i]), 0.0));
	}
}

void polar_mat_ranrot(param& key, gsl_rng* r) {
	if (key.polar->size1 != 3 || key.polar->size2 != 4) {
		cout << "Polar matrix size should be 3 * 4" << '\n';
		exit(EXIT_FAILURE);
	}

	double a[3] = {0};
	double b[3] = {0};
	double c[3] = {0};

	gsl_matrix_complex* rotmat = gsl_matrix_complex_alloc(3, 3);
	gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(3, 4);

	
	gsl_ran_dir_3d(r, c, c + 1, c + 2);
	vector3_norm(c);

	gsl_ran_dir_3d(r, a, a + 1, a + 2);


	double product = a[0] * c[0] + a[1] * c[1] + a[2] * c[2];
	a[0] = a[0] - product * c[0];
	a[1] = a[1] - product * c[1];
	a[2] = a[2] - product * c[2];
	vector3_norm(a);

	b[0] = c[1] * a[2] - c[2] * a[1];
	b[1] = c[2] * a[0] - c[0] * a[2];
	b[2] = c[0] * a[1] - c[1] * a[0];

	for (int i = 0; i < 3; i++) {
		gsl_matrix_complex_set(rotmat, 0, i, gsl_complex_rect(a[i], 0.0));
		gsl_matrix_complex_set(rotmat, 1, i, gsl_complex_rect(b[i], 0.0));
		gsl_matrix_complex_set(rotmat, 2, i, gsl_complex_rect(c[i], 0.0));
	}

	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, ONE, rotmat, key.polar, ZERO, tmp);
	gsl_matrix_complex_memcpy(key.polar, tmp);

	gsl_matrix_complex_free(rotmat);
	gsl_matrix_complex_free(tmp);
}

void compute_pulse_interaction(param& key) {
	if (key.X.size() != 4) {
		cout << "The number of pulse polarization should be 4 instead of " << key.X.size() << '\n';
		exit(EXIT_FAILURE);
	}
	int size = key.sys_size;
	gsl_complex E;
	gsl_matrix_complex* mtmp;
	gsl_matrix_complex* X[3];
	for (int i = 0; i < 3; i++) X[i] = gsl_matrix_complex_alloc(key.sys_size, key.sys_size);
	gsl_matrix_complex_memcpy(X[0], key.X_x);
	gsl_matrix_complex_memcpy(X[1], key.X_y);
	gsl_matrix_complex_memcpy(X[2], key.X_z);


	mtmp = gsl_matrix_complex_alloc(size, size);
	for (int i = 0; i < 4; i++) {
		gsl_matrix_complex_set_zero(key.X[i]);
		for (int j = 0; j < 3; j++) {
			E = gsl_matrix_complex_get(key.polar, j, i);
			gsl_matrix_complex_memcpy(mtmp, X[j]);
			gsl_matrix_complex_scale(mtmp, E);
			gsl_matrix_complex_add(key.X[i], mtmp);
		}
	}
	gsl_matrix_complex_free(mtmp);
	for (int i = 0; i < 3; i++) gsl_matrix_complex_free(X[i]);
}

