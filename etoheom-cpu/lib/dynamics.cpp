#include "dynamics.h"

void construct_Hal(param& key, gsl_matrix_complex* Hal, gsl_rng* r) {



	double sigma = 0;
	double tmp = 0;
	gsl_matrix_complex* deltaH = gsl_matrix_complex_alloc(key.sys_size, key.sys_size);
	gsl_matrix_complex_set_all(deltaH, gsl_complex_rect(0.0, 0.0));


	for (int i = 0; i < key.single_size; i++) {
		for (int j = 0; j < key.single_size; j++) {
			sigma = GSL_REAL(gsl_matrix_complex_get(key.Disorder, i, j));
			if (fabs(sigma) < 1e-10) tmp = 0;
			else tmp = gsl_ran_gaussian(r, sigma);
			gsl_matrix_complex_set(deltaH, i, j, gsl_complex_rect(tmp, 0.0));
		}
	}
	gsl_matrix_complex_memcpy(Hal, key.Hal);
	gsl_matrix_complex_add(Hal, deltaH);
  
  vector<double> Htmp(key.single_size, 0.0);
  for(int i = 1; i < key.single_size; i++) Htmp[i] = GSL_REAL(gsl_matrix_complex_get(Hal, i, i));
  sort(Htmp.begin(), Htmp.end());
  for(int i = 1; i < key.single_size; i++) gsl_matrix_complex_set(Hal, i, i, gsl_complex_rect(Htmp[i], 0.0));
  
	for (int i = key.single_size; i < key.sys_size; i++) {
		for (int j = key.single_size; j < key.sys_size; j++) {
			if (i == j) {
				int a = key.d_map[i].first;
				int b = key.d_map[i].second;
				gsl_matrix_complex_set(Hal, i, i, gsl_complex_add(gsl_matrix_complex_get(Hal, a, a), gsl_matrix_complex_get(Hal, b, b)));
			}
			else {
				int a = key.d_map[i].first;
				int b = key.d_map[i].second;
				int c = key.d_map[j].first;
				int d = key.d_map[j].second;
				double x1 = delta(a, c) * (1 - delta(b, d)) * GSL_REAL(gsl_matrix_complex_get(Hal, b, d));
				double x2 = delta(a, d) * (1 - delta(b, c)) * GSL_REAL(gsl_matrix_complex_get(Hal, b, c));
				double x3 = delta(b, c) * (1 - delta(a, d)) * GSL_REAL(gsl_matrix_complex_get(Hal, a, d));
				double x4 = delta(b, d) * (1 - delta(a, c)) * GSL_REAL(gsl_matrix_complex_get(Hal, a, c));
				gsl_matrix_complex_set(Hal, i, j, gsl_complex_rect(x1 + x2 + x3 + x4, 0.0));

			}
		}
	}
	//print_Hal(Hal);
}

void total_ADO_dynamics(const param& key, const vector<gsl_matrix_complex*>& rho_copy, vector<gsl_matrix_complex*>& drho, const gsl_matrix_complex* H) {
	
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;
	int K = key.K;
	int K_m = key.K_m;

	for (int i = 0; i < key.ado.size(); i++) {
		//L_s
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, ICNT, rho_copy[i], H, ZERO, drho[i]);
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, MINUS_ICNT, H, rho_copy[i], ONE, drho[i]);

		if (i != 0) {
			gsl_complex x = gsl_complex_rect(0.0, 0.0);
			for (int j = 0; j < K; j++) {
				for (int k = 0; k < K_m; k++) {
					int offset = j + k * K;
					if ((int)key.ado[i][offset] == 0) continue;
					gsl_complex y = gsl_complex_mul(MINUS_ONE, gsl_complex_mul_real(key.gamma[j][k], (int)key.ado[i][offset]));
					x = gsl_complex_add(x, y);
				}
			}
			cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(rho_copy[i], 0, 0), 1, gsl_matrix_complex_ptr(drho[i], 0, 0), 1);
		}

		for (int j = 0; j < K; j++) {
			for (int k = 0; k < K_m; k++) {
				int offset = j + k * K;
				if ((int)key.ado[i][offset] >= key.L) continue;
				string tmp(key.ado[i]);
				tmp[offset] += 1;
				if (key.ado_map.find(tmp) != key.ado_map.end()) {
					double x = sqrt(((double)key.ado[i][offset] + 1) * gsl_complex_abs(key.alpha[j][k]));
					//double x = 1.0;
					gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_mul_real(MINUS_ICNT, x), rho_copy[key.ado_map.at(tmp)], key.S[j], ONE, drho[i]);
					gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_mul_real(ICNT, x), key.S[j], rho_copy[key.ado_map.at(tmp)], ONE, drho[i]);
				}
			}
		}

		//Lower
		if (i != 0) {
			for (int j = 0; j < K; j++) {
				for (int k = 0; k < K_m; k++) {
					int offset = j + k * K;
					if ((int)key.ado[i][offset] <= 0) continue;
					string tmp(key.ado[i]);
					tmp[offset] -= 1;
					gsl_complex x = gsl_complex_rect(0.0, 0.0);
					if (key.ado_map.find(tmp) != key.ado_map.end()) {
						double y = sqrt((double)key.ado[i][offset] / gsl_complex_abs(key.alpha[j][k]));
						//double y = (double)key.ado[i][offset];
						x = gsl_complex_mul(ICNT, gsl_complex_mul_real(key.alpha[j][k], y));
						gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, x, key.S[j], rho_copy[key.ado_map.at(tmp)], ONE, drho[i]);
						x = gsl_complex_mul(MINUS_ICNT, gsl_complex_mul_real(key.alpha_t[j][k], y));
						gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, x, rho_copy[key.ado_map.at(tmp)], key.S[j], ONE, drho[i]);

					}
				}
			}
		}
	}
}

void total_ADO_dynamics_Ht(param& key, const vector<gsl_matrix_complex*>& rho_copy, vector<gsl_matrix_complex*>& drho, const gsl_matrix_complex* Ht) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;
	int K = key.K;
	int K_m = key.K_m;
	gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(sys_size, sys_size);

	for (int i = 0; i < key.ado.size(); i++) {
		//L_s
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, ICNT, rho_copy[i], Ht, ZERO, drho[i]);
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, MINUS_ICNT, Ht, rho_copy[i], ONE, drho[i]);

		if (i != 0) {
			gsl_complex x = gsl_complex_rect(0.0, 0.0);
			for (int j = 0; j < K; j++) {
				for (int k = 0; k < K_m; k++) {
					int offset = j + k * K;
					if ((int)key.ado[i][offset] == 0) continue;
					gsl_complex y = gsl_complex_mul(MINUS_ONE, gsl_complex_mul_real(key.gamma[j][k], (int)key.ado[i][offset]));
					x = gsl_complex_add(x, y);
				}
			}
			cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(rho_copy[i], 0, 0), 1, gsl_matrix_complex_ptr(drho[i], 0, 0), 1);
		}


		for (int j = 0; j < K; j++) {
			for (int k = 0; k < K_m; k++) {
				int offset = j + k * K;
				if ((int)key.ado[i][offset] >= key.L) continue;
				string tmp(key.ado[i]);
				tmp[offset] += 1;
				if (key.ado_map.find(tmp) != key.ado_map.end()) {
					double x = sqrt(((double)key.ado[i][offset] + 1) * gsl_complex_abs(key.alpha[j][k]));
					//double x = 1.0;
					gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_mul_real(MINUS_ICNT, x), rho_copy[key.ado_map.at(tmp)], key.S[j], ONE, drho[i]);
					gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_mul_real(ICNT, x), key.S[j], rho_copy[key.ado_map.at(tmp)], ONE, drho[i]);
				}
			}
		}


		//Lower
		if (i != 0) {
			for (int j = 0; j < K; j++) {
				for (int k = 0; k < K_m; k++) {
					int offset = j + k * K;
					if ((int)key.ado[i][offset] <= 0) continue;
					string tmp(key.ado[i]);
					tmp[offset] -= 1;
					gsl_complex x = gsl_complex_rect(0.0, 0.0);
					if (key.ado_map.find(tmp) != key.ado_map.end()) {
						double y = sqrt((double)key.ado[i][offset] / gsl_complex_abs(key.alpha[j][k]));
						//double y = (double)key.ado[i][offset];
						x = gsl_complex_mul(ICNT, gsl_complex_mul_real(key.alpha[j][k], y));
						gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, x, key.S[j], rho_copy[key.ado_map.at(tmp)], ONE, drho[i]);
						x = gsl_complex_mul(MINUS_ICNT, gsl_complex_mul_real(key.alpha_t[j][k], y));
						gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, x, rho_copy[key.ado_map.at(tmp)], key.S[j], ONE, drho[i]);

					}
				}
			}
		}
		/*
		if (key.K_m == 1) {
			//Temp correction
			for (int j = 0; j < K; j++) {
				gsl_complex x1 = gsl_complex_rect(0.0, 0.0);
				gsl_complex x2 = gsl_complex_div(gsl_complex_rect(2 * key.lambda[j] / key.beta / h_bar, 0.0), key.gamma[j][0]);
				for (int k = 0; k < K_m; k++) {
					x1 = gsl_complex_add(x1, gsl_complex_div(key.alpha[j][k], key.gamma[j][k]));
				}
				gsl_complex x = gsl_complex_sub(x1, x2);
				gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, ONE, key.S[j], rho_copy[i], ZERO, tmp);
				gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, MINUS_ONE, rho_copy[i], key.S[j], ONE, tmp);
				gsl_matrix_complex_scale(tmp, x);
				gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, ONE, key.S[j], tmp, ONE, drho[i]);
				gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, MINUS_ONE, tmp, key.S[j], ONE, drho[i]);
			}
		}
		*/
	}
	gsl_matrix_complex_free(tmp);
}

void dynamics_solver(param& key, gsl_matrix_complex* H, vector<int>& sites, vector<vector<double>>& population) {
	key.reset_rho();
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;
	if (sites.size() > sys_size) {
		cout << "number of output sites should not over the system size " << sys_size << "." << '\n';
		exit(EXIT_FAILURE);
	}
	for (int site : sites) {
		if (site > sys_size || site < 0) {
			cout << "site " << site << " is out of range" << '\n';
			exit(EXIT_FAILURE);
		}
	}

	vector<double> b = { 1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6 };

	vector<gsl_matrix_complex*> k1(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k1[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k1[i], gsl_complex_rect(0.0, 0.0));
	}
	vector<gsl_matrix_complex*> k2(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k2[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k2[i], gsl_complex_rect(0.0, 0.0));
	}
	vector<gsl_matrix_complex*> k3(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k3[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k3[i], gsl_complex_rect(0.0, 0.0));
	}
	vector<gsl_matrix_complex*> k4(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k4[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k4[i], gsl_complex_rect(0.0, 0.0));
	}

	float nstep = (key.t_end - key.t_start) / key.step_size;
	float t = key.t_start;
	int print_step = (int)(key.print_step / key.step_size);
	int ct = 0;

	for (int i = 0; i <= nstep; i++) {
		if (i % print_step == 0) {
			for (int j = 0; j < sites.size(); j++) {
				gsl_complex r = gsl_matrix_complex_get(key.rho[0], sites[j] - 1, sites[j] - 1);
				//files[j] << GSL_REAL(r) << " " << GSL_IMAG(r) << '\n';
				population[j][ct] += GSL_REAL(r);
			}
			ct++;
		}
		//f = drho/dt = f(rho)
		//k1 = f(rho)
		gsl_complex x;
		total_ADO_dynamics(key, key.rho, k1, H);

		x = gsl_complex_rect(key.step_size * 0.5, 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k1[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		total_ADO_dynamics(key, key.rho, k2, H);
		//for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_memcpy(rho_copy[j], key.rho[j]);
		x = gsl_complex_mul_real(x, -1);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k1[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);

		x = gsl_complex_rect(key.step_size * 0.5, 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k2[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		total_ADO_dynamics(key, key.rho, k3, H);
		//for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_memcpy(rho_copy[j], key.rho[j]);
		x = gsl_complex_mul_real(x, -1);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k2[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);

		x = gsl_complex_rect(key.step_size, 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k3[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		total_ADO_dynamics(key, key.rho, k4, H);
		//for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_memcpy(rho_copy[j], key.rho[j]);
		x = gsl_complex_mul_real(x, -1);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k3[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);

		//new rho
		x = gsl_complex_rect(key.step_size * b[0], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k1[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		x = gsl_complex_rect(key.step_size * b[1], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k2[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		x = gsl_complex_rect(key.step_size * b[2], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k3[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		x = gsl_complex_rect(key.step_size * b[3], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k4[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		t += key.step_size;
	}
	for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_free(k1[j]);
	for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_free(k2[j]);
	for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_free(k3[j]);
	for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_free(k4[j]);
}


void twoD_spectrum_solver(param& key, gsl_matrix_complex* H, const int nv1, const int nv2, const int nv3, vector<gsl_complex>& polarization) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;
	key.reset_rho();

	vector<double> b = { 1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6 };

	gsl_matrix_complex* X_site_trans_2 = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_transpose_memcpy(X_site_trans_2, key.X[1]);
	gsl_matrix_complex* X_site_trans_3 = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_transpose_memcpy(X_site_trans_3, key.X[2]);


	vector<gsl_matrix_complex*> k1(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k1[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k1[i], gsl_complex_rect(0.0, 0.0));
	}
	vector<gsl_matrix_complex*> k2(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k2[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k2[i], gsl_complex_rect(0.0, 0.0));
	}
	vector<gsl_matrix_complex*> k3(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k3[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k3[i], gsl_complex_rect(0.0, 0.0));
	}
	vector<gsl_matrix_complex*> k4(key.ado.size());
	for (int i = 0; i < key.ado.size(); i++) {
		k4[i] = gsl_matrix_complex_alloc(sys_size, sys_size);
		gsl_matrix_complex_set_all(k4[i], gsl_complex_rect(0.0, 0.0));
	}

	gsl_matrix_complex* V1 = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_set_all(V1, gsl_complex_rect(0.0, 0.0));
	gsl_matrix_complex* V2 = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_set_all(V2, gsl_complex_rect(0.0, 0.0));
	gsl_matrix_complex* V3 = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_set_all(V3, gsl_complex_rect(0.0, 0.0));
	gsl_matrix_complex* V = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_set_all(V, gsl_complex_rect(0.0, 0.0));
	gsl_matrix_complex* P = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_set_all(P, gsl_complex_rect(0.0, 0.0));
	//ofstream file1("site1.txt");
	//ofstream file2("site2.txt");
	float nstep = (key.t_end - key.t_start) / key.step_size;
	float t = key.t_start;
	int print_step = (int)(key.print_step / key.step_size);
	for (int i = 0; i <= nstep; i++) {
		gsl_matrix_complex_memcpy(V, H);
		// pulse 1
		if (nv1 == 1 && fabs(t - key.pulses[0].tau0) <= key.pulses[0].taup * 2) {
			build_V_matrix(key, key.X[0], V1, t, 0);
			gsl_matrix_complex_sub(V, V1);
		}
		// pulse 2
		if (nv2 == 1 && fabs(t - key.pulses[1].tau0) <= key.pulses[1].taup * 2) {
			build_V_dagger_matrix(key, X_site_trans_2, V2, t, 1);
			gsl_matrix_complex_sub(V, V2);
		}
		// pulse 3
		if (nv3 == 1 && fabs(t - key.pulses[2].tau0) <= key.pulses[2].taup * 2) {
			build_V_dagger_matrix(key, X_site_trans_3, V3, t, 2);
			gsl_matrix_complex_sub(V, V3);
		}
		//print_matrix_real(V);


		if (t >= key.pulses[2].tau0 && i % print_step == 0) {
			//print_matrix(key.rho[0]);
			gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, ONE, key.rho[0], key.X[3], ZERO, P);
			//cout << "t = " << t << '\n';
			//print_matrix(P);
			//puts("");
			///////////////////////////////////////////////
			// Tr(X*rho)
			gsl_complex tr = gsl_complex_rect(0.0, 0.0);
			for (int j = 0; j < sys_size; j++) {
				tr = gsl_complex_add(tr, gsl_matrix_complex_get(P, j, j));
			}
			//cout << GSL_REAL(tr) << '\n';
			polarization.push_back(tr);

			///////////////////////////////////////////////
		}

		gsl_complex x;
		total_ADO_dynamics_Ht(key, key.rho, k1, V);

		x = gsl_complex_rect(key.step_size * 0.5, 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k1[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		total_ADO_dynamics_Ht(key, key.rho, k2, V);
		//for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_memcpy(rho_copy[j], key.rho[j]);
		x = gsl_complex_mul_real(x, -1);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k1[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);

		x = gsl_complex_rect(key.step_size * 0.5, 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k2[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		total_ADO_dynamics_Ht(key, key.rho, k3, V);
		//for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_memcpy(rho_copy[j], key.rho[j]);
		x = gsl_complex_mul_real(x, -1);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k2[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);

		x = gsl_complex_rect(key.step_size, 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k3[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		total_ADO_dynamics_Ht(key, key.rho, k4, V);
		//for (int j = 0; j < key.ado.size(); j++) gsl_matrix_complex_memcpy(rho_copy[j], key.rho[j]);
		x = gsl_complex_mul_real(x, -1);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k3[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);

		//new rho
		x = gsl_complex_rect(key.step_size * b[0], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k1[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		x = gsl_complex_rect(key.step_size * b[1], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k2[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		x = gsl_complex_rect(key.step_size * b[2], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k3[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		x = gsl_complex_rect(key.step_size * b[3], 0.0);
		for (int j = 0; j < key.ado.size(); j++) cblas_zaxpy(total_size, &x, gsl_matrix_complex_ptr(k4[j], 0, 0), 1, gsl_matrix_complex_ptr(key.rho[j], 0, 0), 1);
		t += key.step_size;
	}

	for (int i = 0; i < key.ado.size(); i++) gsl_matrix_complex_free(k1[i]);
	for (int i = 0; i < key.ado.size(); i++) gsl_matrix_complex_free(k2[i]);
	for (int i = 0; i < key.ado.size(); i++) gsl_matrix_complex_free(k3[i]);
	for (int i = 0; i < key.ado.size(); i++) gsl_matrix_complex_free(k4[i]);
	gsl_matrix_complex_free(V1);
	gsl_matrix_complex_free(V2);
	gsl_matrix_complex_free(V3);
	gsl_matrix_complex_free(V);
	gsl_matrix_complex_free(P);
	gsl_matrix_complex_free(X_site_trans_2);
	gsl_matrix_complex_free(X_site_trans_3);
}





