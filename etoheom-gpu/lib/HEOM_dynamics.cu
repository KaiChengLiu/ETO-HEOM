#include "HEOM_dynamics.h"

int delta1(const int a, const int b) {
	return a == b ? 1 : 0;
}

void construct_Hal(param& key, vector<data_type>& Hal) {
	int sys_size = key.sys_size;
	float sigma = 0;
	float tmp = 0;
	vector<data_type> deltaH(key.single_size * key.single_size, make_cuComplex(0.0, 0.0));

	for (int i = 0; i < key.single_size; i++) {
		for (int j = 0; j < key.single_size; j++) {
			sigma = key.Disorder[i * key.single_size + j].x;
			if (fabs(sigma) < 1e-10) tmp = 0;
			else tmp = generateRandomGaussian(0, sigma);
			deltaH[i * key.single_size + j] = make_cuComplex(tmp, 0.0);
		}
	}

	Hal = key.Hal;
	for (int i = 0; i < key.single_size; i++) {
		for (int j = 0; j < key.single_size; j++) {
			Hal[i * sys_size + j] = cuCaddf(Hal[i * sys_size + j], deltaH[i * key.single_size + j]);
		}
	}

	for (int i = key.single_size; i < key.sys_size; i++) {
		for (int j = key.single_size; j < key.sys_size; j++) {
			if (i == j) {
				int a = key.d_map[i].first;
				int b = key.d_map[i].second;
				Hal[i * sys_size + j] = cuCaddf(Hal[a * sys_size + a], Hal[b * sys_size + b]);
			}
			else {
				int a = key.d_map[i].first;
				int b = key.d_map[i].second;
				int c = key.d_map[j].first;
				int d = key.d_map[j].second;
				double x1 = delta1(a, c) * (1 - delta1(b, d)) * Hal[b * sys_size + d].x;
				double x2 = delta1(a, d) * (1 - delta1(b, c)) * Hal[b * sys_size + c].x;
				double x3 = delta1(b, c) * (1 - delta1(a, d)) * Hal[a * sys_size + d].x;
				double x4 = delta1(b, d) * (1 - delta1(a, c)) * Hal[a * sys_size + c].x;
				Hal[i * sys_size + j] = make_cuComplex(x1 + x2 + x3 + x4, 0.0);
			}
		}
	}


}


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
	key.reset_rho();
}

std::unordered_map<std::string, int> Build_ADO_map(const std::vector<std::string> ado) {
	std::unordered_map<std::string, int> ado_map;
	for (int i = 0; i < ado.size(); i++) ado_map[ado[i]] = i;
	return ado_map;
}

void total_ADO_dynamics(cublasHandle_t& cublasH, const param& key, const data_type* d_rho_copy, data_type* d_drho) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;

	for (int i = 0; i < key.ado.size(); i++) {

		// Compute L_s operation
		cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
			sys_size, sys_size, sys_size,
			&ICNT, d_rho_copy + total_size * i, sys_size, key.d_Hal, sys_size,
			&ZERO, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
			sys_size, sys_size, sys_size,
			&MINUS_ICNT, key.d_Hal, sys_size, d_rho_copy + total_size * i, sys_size,
			&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		// Apply damping term
		if (i != 0) {
			for (int j = 0; j < key.K; j++) {
				for (int k = 0; k < key.K_m; k++) {
					int offset = j + k * key.K;
					if ((int)key.ado[i][offset] == 0) continue;

					// FIX: Store `cuCmulf(...)` result in a variable before passing to function
					cuComplex coeff = cuCmulf(MINUS_ONE, make_cuComplex(
						key.gamma[j][k].x * (int)key.ado[i][offset],
						key.gamma[j][k].y * (int)key.ado[i][offset]
					));

					cublasError(cublasCaxpy(cublasH, total_size, &coeff, d_rho_copy + total_size * i, 1, d_drho + total_size * i, 1), __FILE__, __LINE__);
					cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
				}
			}
		}

		// Handle upper-level terms
		for (int j = 0; j < key.K; j++) {
			for (int k = 0; k < key.K_m; k++) {
				int offset = j + k * key.K;
				if ((int)key.ado[i][offset] >= key.L) continue;
				std::string tmp(key.ado[i]);
				tmp[offset] += 1;
				if (key.ado_map.find(tmp) != key.ado_map.end()) {
					int target_idx = key.ado_map.at(tmp);
					cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
						sys_size, sys_size, sys_size,
						&MINUS_ICNT, d_rho_copy + total_size * target_idx, sys_size, key.d_S + total_size * j, sys_size,
						&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
					cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
					cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
						sys_size, sys_size, sys_size,
						&ICNT, key.d_S + total_size * j, sys_size, d_rho_copy + total_size * target_idx, sys_size,
						&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
					cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
				}
			}
		}

		// Handle lower-level terms
		if (i != 0) {
			for (int j = 0; j < key.K; j++) {
				for (int k = 0; k < key.K_m; k++) {
					int offset = j + k * key.K;
					if ((int)key.ado[i][offset] <= 0) continue;

					std::string tmp(key.ado[i]);
					tmp[offset] -= 1;
					if (key.ado_map.find(tmp) != key.ado_map.end()) {
						int target_idx = key.ado_map.at(tmp);

						// FIX: Store `cuCmulf(...)` result in a variable before passing to function
						cuComplex alpha_coeff = cuCmulf(ICNT, make_cuComplex(
							key.alpha[j][k].x * (int)key.ado[i][offset],
							key.alpha[j][k].y * (int)key.ado[i][offset]
						));

						cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
							sys_size, sys_size, sys_size,
							&alpha_coeff, key.d_S + total_size * j, sys_size, d_rho_copy + total_size * target_idx, sys_size,
							&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
						cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
						// FIX: Store `cuCmulf(...)` result in a variable before passing to function
						cuComplex alpha_t_coeff = cuCmulf(MINUS_ICNT, make_cuComplex(
							key.alpha_t[j][k].x * (int)key.ado[i][offset],
							key.alpha_t[j][k].y * (int)key.ado[i][offset]
						));

						cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
							sys_size, sys_size, sys_size,
							&alpha_t_coeff, d_rho_copy + total_size * target_idx, sys_size, key.d_S + total_size * j, sys_size,
							&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
						cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
					}
				}
			}
		}
	}
}

void total_ADO_dynamics_Ht(cublasHandle_t& cublasH, param& key, const data_type* d_rho_copy, data_type* d_drho, const data_type* d_Ht) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;

	// Allocate temporary memory on device
	data_type* d_tmp = nullptr;
	cudaError(cudaMalloc(&d_tmp, total_size * sizeof(data_type)), __FILE__, __LINE__);

	// Loop over all auxiliary density operators (ADO)
	for (int i = 0; i < key.ado.size(); i++) {

		// Compute L_s operation
		cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
			sys_size, sys_size, sys_size,
			&ICNT, d_rho_copy + total_size * i, sys_size, d_Ht, sys_size,
			&ZERO, d_drho + total_size * i, sys_size), __FILE__, __LINE__);

		cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
			sys_size, sys_size, sys_size,
			&MINUS_ICNT, d_Ht, sys_size, d_rho_copy + total_size * i, sys_size,
			&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);

		// Apply damping term
		if (i != 0) {
			data_type x = make_cuComplex(0.0, 0.0);
			for (int j = 0; j < key.K; j++) {
				for (int k = 0; k < key.K_m; k++) {
					int offset = k * key.K + j;
					if ((int)key.ado[i][offset] != 0) {
						// Store the multiplication result in a variable before adding
						data_type gamma_coeff = cuCmulf(MINUS_ONE, make_cuComplex(
							key.gamma[j][k].x * (int)key.ado[i][offset],
							key.gamma[j][k].y * (int)key.ado[i][offset]
						));
						x = cuCaddf(x, gamma_coeff);
					}
				}
			}
			// Ensure we pass a variable (not an rvalue) to cublasCaxpy
			cublasError(cublasCaxpy_v2(cublasH, total_size, &x, d_rho_copy + total_size * i, 1, d_drho + total_size * i, 1), __FILE__, __LINE__);
		}

		// Apply upper level terms
		for (int j = 0; j < key.K; j++) {
			for (int k = 0; k < key.K_m; k++) {
				std::string tmp(key.ado[i]);
				int offset = k * key.K + j;
				tmp[offset] += 1;
				if (key.ado_map.find(tmp) != key.ado_map.end()) {
					int target_idx = key.ado_map.at(tmp);

					// Store computed values before passing to functions
					data_type x = make_cuComplex(sqrt(((int)key.ado[i][offset] + 1) * cuComplex_abs(key.alpha[j][k])), 0.0);
					data_type minus_icnt_x = cuCmulf(MINUS_ICNT, x);
					data_type icnt_x = cuCmulf(ICNT, x);

					cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
						sys_size, sys_size, sys_size,
						&minus_icnt_x, d_rho_copy + total_size * target_idx, sys_size, key.d_S + total_size * j, sys_size,
						&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);

					cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
						sys_size, sys_size, sys_size,
						&icnt_x, key.d_S + total_size * j, sys_size, d_rho_copy + total_size * target_idx, sys_size,
						&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
				}
			}
		}

		// Apply lower level terms
		if (i != 0) {
			for (int j = 0; j < key.K; j++) {
				for (int k = 0; k < key.K_m; k++) {
					std::string tmp(key.ado[i]);
					int offset = k * key.K + j;
					tmp[offset] -= 1;
					if (key.ado_map.find(tmp) != key.ado_map.end()) {
						int target_idx = key.ado_map.at(tmp);

						// Compute scaling factor
						data_type y = make_cuComplex(sqrt((int)key.ado[i][offset] / cuComplex_abs(key.alpha[j][k])), 0.0);

						// Store results before using in function calls
						data_type alpha_x = cuCmulf(y, cuCmulf(ICNT, key.alpha[j][k]));
						data_type alpha_t_x = cuCmulf(y, cuCmulf(MINUS_ICNT, key.alpha_t[j][k]));

						cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
							sys_size, sys_size, sys_size,
							&alpha_x, key.d_S + total_size * j, sys_size, d_rho_copy + total_size * target_idx, sys_size,
							&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);

						cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
							sys_size, sys_size, sys_size,
							&alpha_t_x, d_rho_copy + total_size * target_idx, sys_size, key.d_S + total_size * j, sys_size,
							&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
					}
				}
			}
		}
	}

	// Free allocated memory
	cudaError(cudaFree(d_tmp), __FILE__, __LINE__);
}

void total_ADO_dynamics_Ht_batch(cublasHandle_t& cublasH, param& key, const data_type* d_rho_copy, data_type* d_drho, const data_type* d_Ht, const data_type* const* d_rho_batch, data_type** d_drho_batch, data_type* d_work_space, data_type** d_work_batch, data_type* d_rho_work_space, data_type** d_rho_work_batch, vector<cudaStream_t>& streams) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;

	data_type* d_tmp = nullptr;
	cudaError(cudaMalloc(&d_tmp, total_size * sizeof(data_type)), __FILE__, __LINE__);
	//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

	// L_s	
	cudaError(cudaMemcpy(d_work_space, d_Ht, total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	cublasError(cublasCgemmBatched(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
		sys_size, sys_size, sys_size,
		&ICNT, d_rho_batch, sys_size, d_work_batch, sys_size,
		&ZERO, d_drho_batch, sys_size, key.ado.size()), __FILE__, __LINE__);
	//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
	cublasError(cublasCgemmBatched(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
		sys_size, sys_size, sys_size,
		&MINUS_ICNT, d_work_batch, sys_size, d_rho_batch, sys_size,
		&ONE, d_drho_batch, sys_size, key.ado.size()), __FILE__, __LINE__);
	//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

	// Upper
	cudaError(cudaMemset(d_rho_work_space, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
	for (int j = 0; j < key.K; j++) {
		cublasError(cublasCcopy(cublasH, total_size, key.d_S + total_size * j, 1, d_work_space, 1), __FILE__, __LINE__);
		for (int k = 0; k < key.K_m; k++) {
			int offset = k * key.K + j;
			vector<data_type*> ptr_tmp;
			for (int i = 0; i < key.ado.size(); i++) {
				std::string tmp(key.ado[i]);
				tmp[offset] += 1;
				if (key.ado_map.find(tmp) != key.ado_map.end()) {
					int target_idx = key.ado_map.at(tmp);
					ptr_tmp.push_back(key.d_rho_ptr[target_idx]);
				}
				else ptr_tmp.push_back(d_rho_work_space);
			}
			cudaError(cudaMemcpy(d_rho_work_batch, ptr_tmp.data(), key.ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);

			cublasError(cublasCgemmBatched(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
				sys_size, sys_size, sys_size,
				&MINUS_ICNT, d_rho_work_batch, sys_size, d_work_batch, sys_size,
				&ONE, d_drho_batch, sys_size, key.ado.size()), __FILE__, __LINE__);
			//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
			cublasError(cublasCgemmBatched(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
				sys_size, sys_size, sys_size,
				&ICNT, d_work_batch, sys_size, d_rho_work_batch, sys_size,
				&ONE, d_drho_batch, sys_size, key.ado.size()), __FILE__, __LINE__);
			//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		}
	}

	// Damp
	for (int i = 0; i < key.ado.size(); i++) {
		if (i != 0) {
			data_type x = make_cuComplex(0.0, 0.0);
			for (int j = 0; j < key.K; j++) {
				for (int k = 0; k < key.K_m; k++) {
					int offset = k * key.K + j;
					if ((int)key.ado[i][offset] != 0) {
						x = cuCaddf(x, cuCmulf(MINUS_ONE, { key.gamma[j][k].x * (int)key.ado[i][offset], key.gamma[j][k].y * (int)key.ado[i][offset] }));
						//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
					}
				}
			}
			cublasError(cublasCaxpy_v2(cublasH, total_size, &x, d_rho_copy + total_size * i, 1, d_drho + total_size * i, 1), __FILE__, __LINE__);
			//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		}

		// Lower
		if (i != 0) {
			for (int j = 0; j < key.K; j++) {
				for (int k = 0; k < key.K_m; k++) {
					std::string tmp(key.ado[i]);
					int offset = k * key.K + j;
					tmp[offset] -= 1;
					data_type x = make_cuComplex(0.0, 0.0);
					if (key.ado_map.find(tmp) != key.ado_map.end()) {
						int target_idx = key.ado_map.at(tmp);
						x = cuCmulf(ICNT, key.alpha[j][k]);
						cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
							sys_size, sys_size, sys_size,
							&x, key.d_S + total_size * j, sys_size, d_rho_copy + total_size * target_idx, sys_size,
							&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
						//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
						x = cuCmulf(MINUS_ICNT, key.alpha_t[j][k]);
						cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
							sys_size, sys_size, sys_size,
							&x, d_rho_copy + total_size * target_idx, sys_size, key.d_S + total_size * j, sys_size,
							&ONE, d_drho + total_size * i, sys_size), __FILE__, __LINE__);
						//cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
					}
				}
			}
		}
	}
	// Free allocated memory
	cudaError(cudaFree(d_tmp), __FILE__, __LINE__);
}

void dynamics_solver(param& key, std::vector<int>& sites, vector<vector<float>>& population) {
	
	cublasHandle_t cublasH;
	cublasCreate(&cublasH);

	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;
	if (sites.size() > sys_size) {
		std::cout << "number of output sites should not over the system size " << sys_size << "." << '\n';
		exit(EXIT_FAILURE);
	}
	for (int site : sites) {
		if (site > sys_size || site < 0) {
			std::cout << "site " << site << " is out of range" << '\n';
			exit(EXIT_FAILURE);
		}
	}

	std::vector<float> b = { 1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6 };
	//cudaDeviceSynchronize();
	data_type* d_k1 = nullptr;
	cudaError(cudaMalloc(&d_k1, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k2 = nullptr;
	cudaError(cudaMalloc(&d_k2, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k3 = nullptr;
	cudaError(cudaMalloc(&d_k3, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k4 = nullptr;
	cudaError(cudaMalloc(&d_k4, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);

	float nstep = (key.t_end - key.t_start) / key.step_size;
	float t = key.t_start;
	int print_step = (int)(key.print_step / key.step_size);
	for (int i = 0; i <= nstep; i++) {

		if (i % print_step == 0) {
			for (int j = 0; j < sites.size(); j++) {
				data_type r;
				cudaError(cudaMemcpy(&r, key.d_rho + (sites[j] - 1) * (sys_size + 1), sizeof(data_type), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
				population[j].push_back(r.x);
			}
		}

		//f = drho/dt = f(rho)
		//k1 = f(rho)
		//cudaDeviceSynchronize();
		data_type coef;
		total_ADO_dynamics(cublasH, key, key.d_rho, d_k1);

		//k2 = f(rho + b[1][0]*h*k1)
		coef = make_cuComplex(key.step_size * 0.5, 0.0);
		cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k1, 1, key.d_rho, 1);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		total_ADO_dynamics(cublasH, key, key.d_rho, d_k2);
		coef = make_cuComplex(-key.step_size * 0.5, 0.0);
		cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k1, 1, key.d_rho, 1);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

		//k3 = f(rho + b[2][0]*h*k1 + b[2][1]*h*k2)
		coef = make_cuComplex(key.step_size * 0.5, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		total_ADO_dynamics(cublasH, key, key.d_rho, d_k3);
		coef = make_cuComplex(-key.step_size * 0.5, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

		//k4 = f(rho + b[3][0]*h*k1 + b[3][1]*h*k2+ b[3][2]*h*k3)
		coef = make_cuComplex(key.step_size, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		total_ADO_dynamics(cublasH, key, key.d_rho, d_k4);
		coef = make_cuComplex(-key.step_size, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

		//new rho
		coef = make_cuComplex(key.step_size * b[0], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k1, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[1], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[2], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[3], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k4, 1, key.d_rho, 1), __FILE__, __LINE__);
		cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

		t += key.step_size;
	}

	cudaFree(d_k1);
	cudaFree(d_k2);
	cudaFree(d_k3);
	cudaFree(d_k4);
	cublasDestroy(cublasH);
}


void twoD_spectrum_solver(param& key, const data_type* d_H, const int nv1, const int nv2, const int nv3, std::vector<data_type>& polarization) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;
	key.reset_rho();

	cublasHandle_t cublasH;
	cublasCreate(&cublasH);


	std::vector<float> b = { 1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6 };

	vector<data_type*> d_X(4);
	vector<data_type> X(total_size, make_cuComplex(0.0, 0.0));
	for (int i = 0; i < 4; i++) d_X[i] = To_Device(X);


	compute_pulse_interaction(key, d_X);

	dim3 tpb(ThreadPerBlock, ThreadPerBlock); // threads per block
	dim3 numblocks((sys_size + tpb.x - 1) / tpb.x, (sys_size + tpb.y - 1) / tpb.y);


	data_type* d_X_site_trans_2;
	cudaError(cudaMalloc(&d_X_site_trans_2, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_X_site_trans_2, d_X[1], total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	inplace_transpose << <numblocks, tpb >> > (d_X_site_trans_2, sys_size);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
	data_type* d_X_site_trans_3;
	cudaError(cudaMalloc(&d_X_site_trans_3, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_X_site_trans_3, d_X[2], total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	inplace_transpose << <numblocks, tpb >> > (d_X_site_trans_3, sys_size);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);



	//data_type* d_rho_copy = nullptr;
	//cudaError(cudaMalloc(&d_rho_copy, key.ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k1 = nullptr;
	cudaError(cudaMalloc(&d_k1, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k2 = nullptr;
	cudaError(cudaMalloc(&d_k2, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k3 = nullptr;
	cudaError(cudaMalloc(&d_k3, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k4 = nullptr;
	cudaError(cudaMalloc(&d_k4, total_size * key.ado.size() * sizeof(data_type)), __FILE__, __LINE__);

	data_type* d_V1 = nullptr;
	cudaError(cudaMalloc(&d_V1, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V1, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_V2 = nullptr;
	cudaError(cudaMalloc(&d_V2, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V2, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_V3 = nullptr;
	cudaError(cudaMalloc(&d_V3, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V3, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_V = nullptr;
	cudaError(cudaMalloc(&d_V, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);


	data_type* d_P = nullptr;
	cudaError(cudaMalloc(&d_P, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_ONE = nullptr;
	cudaError(cudaMalloc(&d_ONE, sys_size * sizeof(data_type)), __FILE__, __LINE__);
	init_cuComplex << <numblocks.x, tpb.x >> > (d_ONE, sys_size, ONE);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
	data_type tr;
	data_type* d_tr = nullptr;
	cudaError(cudaMalloc(&d_tr, sizeof(data_type)), __FILE__, __LINE__);

	//std::ofstream file1("site1.txt");
	//std::ofstream file2("site2.txt");

	int nstep = (key.t_end - key.t_start) / key.step_size;
	float t = key.t_start;
	int print_step = (int)(key.print_step / key.step_size);

	for (int i = 0; i <= nstep; i++) {
		cudaError(cudaMemcpy(d_V, d_H, total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
		// pulse 1
		if (nv1 == 1 && fabs(t - key.pulses[0].tau0) <= key.pulses[0].taup * 2) {
			build_V_matrix(cublasH, key, d_X[0], d_V1, t, 0);
			cublasError(cublasCaxpy_v2(cublasH, total_size, &MINUS_ONE, d_V1, 1, d_V, 1), __FILE__, __LINE__);
		}
		// pulse 2
		if (nv2 == 1 && fabs(t - key.pulses[1].tau0) <= key.pulses[1].taup * 2) {
			build_V_dagger_matrix(cublasH, key, d_X_site_trans_2, d_V2, t, 1);
			cublasError(cublasCaxpy_v2(cublasH, total_size, &MINUS_ONE, d_V2, 1, d_V, 1), __FILE__, __LINE__);
		}
		// pulse 3
		if (nv3 == 1 && fabs(t - key.pulses[2].tau0) <= key.pulses[2].taup * 2) {
			build_V_dagger_matrix(cublasH, key, d_X_site_trans_3, d_V3, t, 2);
			cublasError(cublasCaxpy_v2(cublasH, total_size, &MINUS_ONE, d_V3, 1, d_V, 1), __FILE__, __LINE__);
		}


		if (t >= key.pulses[2].tau0 && i % print_step == 0) {
			cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
				sys_size, sys_size, sys_size,
				&ONE, d_X[3], sys_size,
				key.d_rho, sys_size,
				&ZERO, d_P, sys_size), __FILE__, __LINE__);
			///////////////////////////////////////////////
			// Tr(X*rho)
			cublasError(cublasCdotu_v2(cublasH, sys_size, d_P, sys_size + 1, d_ONE, 1, d_tr), __FILE__, __LINE__);
			cudaError(cudaMemcpy(&tr, d_tr, sizeof(data_type), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
			polarization.push_back(tr);
			///////////////////////////////////////////////
		}

		data_type coef;

		//f = drho/dt
		//k1 = f(rho)
		total_ADO_dynamics_Ht(cublasH, key, key.d_rho, d_k1, d_V);

		//k2 = f(rho + b[1][0]*h*k1)
		coef = make_cuComplex(key.step_size * 0.5, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k1, 1, key.d_rho, 1), __FILE__, __LINE__);
		total_ADO_dynamics_Ht(cublasH, key, key.d_rho, d_k2, d_V);
		coef = cuCmulf(coef, MINUS_ONE);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k1, 1, key.d_rho, 1), __FILE__, __LINE__);

		//k3 = f(rho + b[2][0]*h*k1 + b[2][1]*h*k2)
		coef = make_cuComplex(key.step_size * 0.5, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);
		total_ADO_dynamics_Ht(cublasH, key, key.d_rho, d_k3, d_V);
		coef = cuCmulf(coef, MINUS_ONE);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);

		//k4 = f(rho + b[3][0]*h*k1 + b[3][1]*h*k2+ b[3][2]*h*k3)
		coef = make_cuComplex(key.step_size, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);
		total_ADO_dynamics_Ht(cublasH, key, key.d_rho, d_k4, d_V);
		coef = cuCmulf(coef, MINUS_ONE);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);

		//new rho
		coef = make_cuComplex(key.step_size * b[0], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k1, 1, key.d_rho, 1), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[1], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[2], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[3], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, total_size * key.ado.size(), &coef, d_k4, 1, key.d_rho, 1), __FILE__, __LINE__);
		t += key.step_size;
	}

	cublasError(cublasDestroy(cublasH), __FILE__, __LINE__);
	for (int i = 0; i < 4; i++) cudaError(cudaFree(d_X[i]), __FILE__, __LINE__);
	cudaError(cudaFree(d_X_site_trans_2), __FILE__, __LINE__);
	cudaError(cudaFree(d_X_site_trans_3), __FILE__, __LINE__);
	cudaError(cudaFree(d_k1), __FILE__, __LINE__);
	cudaError(cudaFree(d_k2), __FILE__, __LINE__);
	cudaError(cudaFree(d_k3), __FILE__, __LINE__);
	cudaError(cudaFree(d_k4), __FILE__, __LINE__);
	//cudaError(cudaFree(d_rho_copy), __FILE__, __LINE__);
	cudaError(cudaFree(d_V1), __FILE__, __LINE__);
	cudaError(cudaFree(d_V2), __FILE__, __LINE__);
	cudaError(cudaFree(d_V3), __FILE__, __LINE__);
	cudaError(cudaFree(d_V), __FILE__, __LINE__);
	cudaError(cudaFree(d_P), __FILE__, __LINE__);
	cudaError(cudaFree(d_ONE), __FILE__, __LINE__);
	cudaError(cudaFree(d_tr), __FILE__, __LINE__);

	//file1.close();
	//file2.close();

}

void propagation_Ht_batch(param& key, const data_type* d_H, const int nv1, const int nv2, const int nv3, std::vector<data_type>& polarization) {
	int sys_size = key.sys_size;
	int total_size = sys_size * sys_size;
	key.reset_rho();
	key.reset_rho_batch();

	cublasHandle_t cublasH = NULL;
	cublasCreate(&cublasH);

	int stream_num = 16;
	vector<cudaStream_t> streams(stream_num);
	for (int i = 0; i < stream_num; i++) cudaError(cudaStreamCreate(&streams[i]), __FILE__, __LINE__);

	vector<float> b = { 1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6 };

	vector<data_type*> d_X(4);
	vector<data_type> X(total_size, make_cuComplex(0.0, 0.0));
	for (int i = 0; i < 4; i++) d_X[i] = To_Device(X);


	compute_pulse_interaction(key, d_X);

	dim3 tpb(ThreadPerBlock, ThreadPerBlock); // threads per block
	dim3 numblocks((sys_size + tpb.x - 1) / tpb.x, (sys_size + tpb.y - 1) / tpb.y);


	data_type* d_X_site_trans_2;
	cudaError(cudaMalloc(&d_X_site_trans_2, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_X_site_trans_2, d_X[1], total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	inplace_transpose << <numblocks, tpb >> > (d_X_site_trans_2, sys_size);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
	data_type* d_X_site_trans_3;
	cudaError(cudaMalloc(&d_X_site_trans_3, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_X_site_trans_3, d_X[2], total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
	inplace_transpose << <numblocks, tpb >> > (d_X_site_trans_3, sys_size);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);

	/*
	vector<data_type> A(total_size);
	To_Host(d_X[0], A);
	print_matrix(A, sys_size);
	puts("");
	*/

	//data_type* d_rho_copy = nullptr;
	//cudaError(cudaMalloc(&d_rho_copy, key.ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k1 = nullptr;
	cudaError(cudaMalloc(&d_k1, key.ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k2 = nullptr;
	cudaError(cudaMalloc(&d_k2, key.ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k3 = nullptr;
	cudaError(cudaMalloc(&d_k3, key.ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_k4 = nullptr;
	cudaError(cudaMalloc(&d_k4, key.ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);

	vector<data_type*> d_k1_ptr;
	for (int i = 0; i < key.ado.size(); i++) d_k1_ptr.push_back(d_k1 + total_size * i);
	vector<data_type*> d_k2_ptr;
	for (int i = 0; i < key.ado.size(); i++) d_k2_ptr.push_back(d_k2 + total_size * i);
	vector<data_type*> d_k3_ptr;
	for (int i = 0; i < key.ado.size(); i++) d_k3_ptr.push_back(d_k3 + total_size * i);
	vector<data_type*> d_k4_ptr;
	for (int i = 0; i < key.ado.size(); i++) d_k4_ptr.push_back(d_k4 + total_size * i);

	data_type** d_k1_batch = nullptr;
	data_type** d_k2_batch = nullptr;
	data_type** d_k3_batch = nullptr;
	data_type** d_k4_batch = nullptr;
	cudaError(cudaMalloc(&d_k1_batch, key.ado.size() * sizeof(data_type*)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_k1_batch, d_k1_ptr.data(), key.ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);
	cudaError(cudaMalloc(&d_k2_batch, key.ado.size() * sizeof(data_type*)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_k2_batch, d_k2_ptr.data(), key.ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);
	cudaError(cudaMalloc(&d_k3_batch, key.ado.size() * sizeof(data_type*)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_k3_batch, d_k3_ptr.data(), key.ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);
	cudaError(cudaMalloc(&d_k4_batch, key.ado.size() * sizeof(data_type*)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_k4_batch, d_k4_ptr.data(), key.ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);


	data_type* d_work;
	cudaError(cudaMalloc(&d_work, total_size * sizeof(data_type)), __FILE__, __LINE__);
	vector<data_type*> d_work_ptr;
	for (int i = 0; i < key.ado.size(); i++) d_work_ptr.push_back(d_work);
	data_type** d_work_batch = nullptr;
	cudaError(cudaMalloc(&d_work_batch, key.ado.size() * sizeof(data_type*)), __FILE__, __LINE__);
	cudaError(cudaMemcpy(d_work_batch, d_work_ptr.data(), key.ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);

	data_type* d_rho_work;
	cudaError(cudaMalloc(&d_rho_work, key.ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
	vector<data_type*> d_rho_work_ptr;
	for (int i = 0; i < key.ado.size(); i++) d_rho_work_ptr.push_back(d_rho_work);
	data_type** d_rho_work_batch = nullptr;
	cudaError(cudaMalloc(&d_rho_work_batch, key.ado.size() * sizeof(data_type*)), __FILE__, __LINE__);
	for (int i = 0; i < key.ado.size(); i++) cudaError(cudaMemcpy(d_rho_work_batch, d_rho_work_ptr.data(), key.ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);


	data_type* d_V1 = nullptr;
	cudaError(cudaMalloc(&d_V1, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V1, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_V2 = nullptr;
	cudaError(cudaMalloc(&d_V2, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V2, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_V3 = nullptr;
	cudaError(cudaMalloc(&d_V3, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V3, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_V = nullptr;
	cudaError(cudaMalloc(&d_V, total_size * sizeof(data_type)), __FILE__, __LINE__);
	cudaError(cudaMemset(d_V, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);


	data_type* d_P = nullptr;
	cudaError(cudaMalloc(&d_P, total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type* d_ONE = nullptr;
	cudaError(cudaMalloc(&d_ONE, sys_size * sizeof(data_type)), __FILE__, __LINE__);
	init_cuComplex << <numblocks.x, tpb.x >> > (d_ONE, sys_size, ONE);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
	data_type tr;
	data_type* d_tr = nullptr;
	cudaError(cudaMalloc(&d_tr, sizeof(data_type)), __FILE__, __LINE__);

	//std::ofstream file1("site1.txt");
	//std::ofstream file2("site2.txt");

	int nstep = (key.t_end - key.t_start) / key.step_size;
	float t = key.t_start;
	int print_step = (int)(key.print_step / key.step_size);

	for (int i = 0; i <= nstep; i++) {
		cudaError(cudaMemcpy(d_V, d_H, total_size * sizeof(data_type), cudaMemcpyDeviceToDevice), __FILE__, __LINE__);
		// pulse 1
		if (nv1 == 1 && fabs(t - key.pulses[0].tau0) <= key.pulses[0].taup * 2) {
			build_V_matrix(cublasH, key, d_X[0], d_V1, t, 0);
			cublasError(cublasCaxpy_v2(cublasH, total_size, &MINUS_ONE, d_V1, 1, d_V, 1), __FILE__, __LINE__);
		}
		// pulse 2
		if (nv2 == 1 && fabs(t - key.pulses[1].tau0) <= key.pulses[1].taup * 2) {
			build_V_dagger_matrix(cublasH, key, d_X_site_trans_2, d_V2, t, 1);
			cublasError(cublasCaxpy_v2(cublasH, total_size, &MINUS_ONE, d_V2, 1, d_V, 1), __FILE__, __LINE__);
		}
		// pulse 3
		if (nv3 == 1 && fabs(t - key.pulses[2].tau0) <= key.pulses[2].taup * 2) {
			build_V_dagger_matrix(cublasH, key, d_X_site_trans_3, d_V3, t, 2);
			cublasError(cublasCaxpy_v2(cublasH, total_size, &MINUS_ONE, d_V3, 1, d_V, 1), __FILE__, __LINE__);
		}

		/*
		cout << "t = " << t << '\n';
		vector<data_type> X(total_size);
		To_Host(key.d_rho, X);
		print_matrix(X, 1);
		puts("");
		*/


		if (t >= key.pulses[2].tau0 && i % print_step == 0) {
			cublasError(cublasCgemm3m(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
				sys_size, sys_size, sys_size,
				&ONE, d_X[3], sys_size,
				key.d_rho, sys_size,
				&ZERO, d_P, sys_size), __FILE__, __LINE__);
			///////////////////////////////////////////////
			// Tr(X*rho)
			cublasError(cublasCdotu_v2(cublasH, sys_size, d_P, sys_size + 1, d_ONE, 1, d_tr), __FILE__, __LINE__);
			cudaError(cudaMemcpy(&tr, d_tr, sizeof(data_type), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
			polarization.push_back(tr);
			///////////////////////////////////////////////
		}

		data_type coef;

		//f = drho/dt
		//k1 = f(rho)
		total_ADO_dynamics_Ht_batch(cublasH, key, key.d_rho, d_k1, d_V, key.d_rho_batch, d_k1_batch, d_work, d_work_batch, d_rho_work, d_rho_work_batch, streams);

		//k2 = f(rho + b[1][0]*h*k1)
		coef = make_cuComplex(key.step_size * 0.5, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k1, 1, key.d_rho, 1), __FILE__, __LINE__);
		total_ADO_dynamics_Ht_batch(cublasH, key, key.d_rho, d_k2, d_V, key.d_rho_batch, d_k2_batch, d_work, d_work_batch, d_rho_work, d_rho_work_batch, streams);
		coef = cuCmulf(coef, MINUS_ONE);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k1, 1, key.d_rho, 1), __FILE__, __LINE__);

		//k3 = f(rho + b[2][0]*h*k1 + b[2][1]*h*k2)
		coef = make_cuComplex(key.step_size * 0.5, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);
		total_ADO_dynamics_Ht_batch(cublasH, key, key.d_rho, d_k3, d_V, key.d_rho_batch, d_k3_batch, d_work, d_work_batch, d_rho_work, d_rho_work_batch, streams);
		coef = cuCmulf(coef, MINUS_ONE);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);

		//k4 = f(rho + b[3][0]*h*k1 + b[3][1]*h*k2+ b[3][2]*h*k3)
		coef = make_cuComplex(key.step_size, 0.0);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);
		total_ADO_dynamics_Ht_batch(cublasH, key, key.d_rho, d_k4, d_V, key.d_rho_batch, d_k4_batch, d_work, d_work_batch, d_rho_work, d_rho_work_batch, streams);
		coef = cuCmulf(coef, MINUS_ONE);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);

		//new rho
		coef = make_cuComplex(key.step_size * b[0], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k1, 1, key.d_rho, 1), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[1], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k2, 1, key.d_rho, 1), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[2], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k3, 1, key.d_rho, 1), __FILE__, __LINE__);
		coef = make_cuComplex(key.step_size * b[3], 0.0);
		cublasError(cublasCaxpy_v2(cublasH, key.ado.size() * total_size, &coef, d_k4, 1, key.d_rho, 1), __FILE__, __LINE__);
		t += key.step_size;
	}

	cublasError(cublasDestroy(cublasH), __FILE__, __LINE__);
	for (int i = 0; i < stream_num; i++) cudaError(cudaStreamDestroy(streams[i]), __FILE__, __LINE__);
	for (int i = 0; i < 4; i++) cudaError(cudaFree(d_X[i]), __FILE__, __LINE__);
	cudaError(cudaFree(d_X_site_trans_2), __FILE__, __LINE__);
	cudaError(cudaFree(d_X_site_trans_3), __FILE__, __LINE__);
	cudaError(cudaFree(d_k1), __FILE__, __LINE__);
	cudaError(cudaFree(d_k2), __FILE__, __LINE__);
	cudaError(cudaFree(d_k3), __FILE__, __LINE__);
	cudaError(cudaFree(d_k4), __FILE__, __LINE__);
	cudaError(cudaFree(d_work), __FILE__, __LINE__);
	cudaError(cudaFree(d_k1_batch), __FILE__, __LINE__);
	cudaError(cudaFree(d_k2_batch), __FILE__, __LINE__);
	cudaError(cudaFree(d_k3_batch), __FILE__, __LINE__);
	cudaError(cudaFree(d_k4_batch), __FILE__, __LINE__);
	cudaError(cudaFree(d_work_batch), __FILE__, __LINE__);
	cudaError(cudaFree(d_rho_work), __FILE__, __LINE__);
	cudaError(cudaFree(d_rho_work_batch), __FILE__, __LINE__);
	cudaError(cudaFree(d_V1), __FILE__, __LINE__);
	cudaError(cudaFree(d_V2), __FILE__, __LINE__);
	cudaError(cudaFree(d_V3), __FILE__, __LINE__);
	cudaError(cudaFree(d_V), __FILE__, __LINE__);
	cudaError(cudaFree(d_P), __FILE__, __LINE__);
	cudaError(cudaFree(d_ONE), __FILE__, __LINE__);
	cudaError(cudaFree(d_tr), __FILE__, __LINE__);

	//file1.close();
	//file2.close();
}



