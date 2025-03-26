#include "HEOM_param.h"
double delta(int a, int b) {
	return a == b ? 1 : 0;
}

void param::param_2des(string filename) {
	fstream f(filename);
	string line, key_word;
	cublasHandle_t cublasH;
	cublasCreate(&cublasH);
	while (getline(f, line)) {
		if (line.empty() || line[0] == '#') continue;
		istringstream iss(line);
		iss >> key_word;

		if (key_word == "HEOM") {
			//cout << "Reading HEOM data" << '\n';
			getline(f, line);
			istringstream val(line);
			val >> this->K >> this->L;
			this->K_m = 0;
		}

		else if (key_word == "SIZE") {
			//cout << "Reading SIZE data" << '\n';
			getline(f, line);
			istringstream val(line);
			val >> this->single_size;
			this->sys_size = this->single_size + (this->single_size - 1) * (this->single_size - 2) / 2;

			this->sys_size = this->single_size + (this->single_size - 1) * (this->single_size - 2) / 2;
			int ct = this->single_size;
			for (int i = 1; i < this->single_size; i++) {
				for (int j = i; j < this->single_size; j++) {
					if (i != j) {
						this->d_map[ct] = make_pair(i, j);
						ct++;
					}
				}
			}
			
		}

		else if (key_word == "DISORDER") {
			getline(f, line);
			istringstream disorder_stream(line);
			disorder_stream >> this->n_sample;
			vector<data_type> Disorder(this->single_size * this->single_size);
			float val = 0;
			for (int i = 0; i < this->single_size; i++) {
				getline(f, line);
				istringstream row(line);
				int j = 0;
				while (row >> val) {
					Disorder[i + this->single_size * j] = make_cuComplex(val / h_bar.x, 0.0);
					j++;
				}
			}
			this->Disorder = Disorder;
			this->d_Disorder = To_Device(Disorder);
		}

		else if (key_word == "HAMILTONIAN") {
			//cout << "Reading HAMILTONIAN data" << '\n';
			vector<data_type> eigenstate(this->sys_size * this->sys_size, make_cuComplex(0.0, 0.0));
			this->d_eigenstate = To_Device(eigenstate);
			vector<data_type> Hal(this->sys_size * this->sys_size);
			float val;
			for (int i = 0; i < this->sys_size; i++) {
				getline(f, line);
				istringstream row(line);
				int j = 0;
				while (row >> val) {
					Hal[i + this->sys_size * j] = make_cuComplex(val / h_bar.x, 0.0);
					j++;
				}
			}
			this->Hal = Hal;
		}

		else if (key_word == "DIPOLE") {
			//cout << "Reading DIPOLE data" << '\n';
			int size = this->sys_size;
			int total_size = this->sys_size * this->sys_size;

			getline(f, line);
			istringstream ndipole_stream(line);
			int ndipole;
			ndipole_stream >> ndipole;
			vector<data_type> Xx(total_size, { 0.0, 0.0 });
			vector<data_type> Xy(total_size, { 0.0, 0.0 });
			vector<data_type> Xz(total_size, { 0.0, 0.0 });
			vector<data_type> X_site(total_size, { 0.0, 0.0 });
			for (int k = 0; k < ndipole; k++) {
				getline(f, line);
				istringstream dipole_stream(line);
				float x, y, z;
				int i, j;
				dipole_stream >> i >> j >> x >> y >> z;
				Xx[i + this->sys_size * j] = make_cuComplex(x, 0.0);
				Xy[i + this->sys_size * j] = make_cuComplex(y, 0.0);
				Xz[i + this->sys_size * j] = make_cuComplex(z, 0.0);
				X_site[i + this->sys_size * j] = make_cuComplex(sqrt(x * x + y * y + z * z), 0.0);
			}
			for (int i = 0; i < this->sys_size; i++) {
				for (int j = this->single_size; j < this->sys_size; j++) {
					int a = this->d_map[j].first;
					int b = this->d_map[j].second;
					float x1 = Xx[0 + size * a].x;
					float x2 = Xx[0 + size * b].x;
					float x = delta(b, i) * x1 + delta(a, i) * x2;
					Xx[i + size * j] = make_cuComplex(x, 0.0);
				}
			}
			for (int i = 0; i < this->sys_size; i++) {
				for (int j = this->single_size; j < this->sys_size; j++) {
					int a = this->d_map[j].first;
					int b = this->d_map[j].second;
					float x1 = Xy[0 + size * a].x;
					float x2 = Xy[0 + size * b].x;
					float x = delta(b, i) * x1 + delta(a, i) * x2;
					Xy[i + size * j] = make_cuComplex(x, 0.0);
				}
			}
			for (int i = 0; i < this->sys_size; i++) {
				for (int j = this->single_size; j < this->sys_size; j++) {
					int a = this->d_map[j].first;
					int b = this->d_map[j].second;
					float x1 = Xz[0 + size * a].x;
					float x2 = Xz[0 + size * b].x;
					float x = delta(b, i) * x1 + delta(a, i) * x2;
					Xz[i + size * j] = make_cuComplex(x, 0.0);
				}
			}
			for (int i = 0; i < this->sys_size; i++) {
				for (int j = 0; j < this->sys_size; j++) {
					float x, y, z;
					x = Xx[i + size * j].x;
					y = Xy[i + size * j].x;
					z = Xz[i + size * j].x;
					X_site[i + size * j] = make_cuComplex(sqrt(x * x + y * y + z * z), 0.0);
				}
			}
			//print_matrix_real(X_site, size);
			this->d_Xx = To_Device(Xx);
			this->d_Xy = To_Device(Xy);
			this->d_Xz = To_Device(Xz);
			this->d_X_site = To_Device(X_site);
		}

		else if (key_word == "POLARIZATION") {
			getline(f, line);
			istringstream polarization_stream(line);
			vector<float> angle(4, 0.0);
			polarization_stream >> angle[0] >> angle[1] >> angle[2] >> angle[3];
			for (int i = 0; i < 4; i++) this->angle.push_back(angle[i] * pi / 180);
		}

		else if (key_word == "TEMPERATURE") {
			//cout << "Reading TEMPERATURE data" << '\n';
			getline(f, line);
			istringstream temperature_stream(line);
			float temperature;
			temperature_stream >> temperature;
			this->beta = cuCdivf(ONE, cuCmulf({ temperature, 0.0 }, boltz_k));
		}

		else if (key_word == "BATHTYPE") {
			//cout << "Reading BATHTYPE data" << '\n';
			getline(f, line);
			istringstream bathtype_stream(line);
			bathtype_stream >> this->bath_type;
		}

		else if (key_word == "BATH") {
			//cout << "Reading BATH data" << '\n';
			if (this->bath_type == "artificial") {
				getline(f, line);
				istringstream bath_stream(line);
				float l;
				bath_stream >> l;
				for (int i = 0; i < this->K; i++) this->lambda.push_back(make_cuComplex(l, 0.0));
				int total_size = this->sys_size * this->sys_size;
				vector<data_type> S(total_size * this->K, { 0.0, 0.0 });
				for (int i = 0; i < this->K; i++) {
					S[i * total_size + (i + 1) + this->sys_size * (i + 1)] = { 1.0, 0.0 };
					this->Hal[(i + 1) + this->sys_size * (i + 1)].x += l;
				}


				this->d_S = To_Device(S);
				this->d_Hal = To_Device(this->Hal);
				int n_modes;
				getline(f, line);
				istringstream n_modes_stream(line);
				n_modes_stream >> n_modes;

				vector<data_type> alpha;
				vector<data_type> alpha_t;
				vector<data_type> gamma;
				for (int i = 0; i < n_modes; i++) {
					double a, b, g, w;
					getline(f, line);
					istringstream modes_stream(line);
					modes_stream >> a >> b >> g >> w;
					alpha.push_back(make_cuComplex(a, b));
					alpha_t.push_back(make_cuComplex(a, -b));
					gamma.push_back(make_cuComplex(g, w));
				}
				for (int i = 0; i < this->K; i++) {
					this->alpha.push_back(alpha);
					this->alpha_t.push_back(alpha_t);
					this->gamma.push_back(gamma);
				}
				this->K_m = this->gamma[0].size();
			}
		}

		else if (key_word == "PULSE") {
			//cout << "Reading PULSE data\n";
			int n;
			getline(f, line);
			istringstream npulse_stream(line);
			npulse_stream >> n;
			if (n != 3) {
				cout << "Only support 3 pulses simulation.\n";
				exit(EXIT_FAILURE);
			}
			for (int i = 0; i < n; i++) {
				getline(f, line);
				istringstream pulse_stream(line);
				pulse_stream >> this->pulses[i].E0 >> this->pulses[i].tau0 >> this->pulses[i].taup >> this->pulses[i].w0;
				this->pulses[i].w0 = this->pulses[i].w0 / h_bar.x; //cm^-1 to fs^-1
				this->pulses[i].E0 = this->pulses[i].E0 / h_bar.x; //cm^-1 to fs^-1
			}
		}

		else if (key_word == "TIME") {
			//cout << "Reading TIME data" << '\n';
			getline(f, line);
			istringstream time_stream(line);
			time_stream >> this->t_start >> this->t_end >> this->step_size >> this->print_step;
		}
	}
	cublasError(cublasDestroy(cublasH), __FILE__, __LINE__);
}

void param::param_dynamics(string filename) {
	{
		fstream f(filename);
		string line, key_word;
		cublasHandle_t cublasH;
		cublasCreate(&cublasH);
		while (getline(f, line)) {
			if (line.empty() || line[0] == '#') continue;
			istringstream iss(line);
			iss >> key_word;

			if (key_word == "HEOM") {
				//cout << "Reading HEOM data" << '\n';
				getline(f, line);
				istringstream val(line);
				val >> this->K >> this->L;
				this->K_m = 0;
			}

			else if (key_word == "SIZE") {
				//cout << "Reading SIZE data" << '\n';
				getline(f, line);
				istringstream val(line);
				val >> this->single_size;
				this->sys_size = this->single_size;
			}

			else if (key_word == "HAMILTONIAN") {
				//cout << "Reading HAMILTONIAN data" << '\n';
				vector<data_type> eigenstate(this->sys_size * this->sys_size, make_cuComplex(0.0, 0.0));
				this->d_eigenstate = To_Device(eigenstate);
				vector<data_type> Hal(this->sys_size * this->sys_size);
				float val;
				for (int i = 0; i < this->sys_size; i++) {
					getline(f, line);
					istringstream row(line);
					int j = 0;
					while (row >> val) {
						Hal[i + this->sys_size * j] = make_cuComplex(val / h_bar.x, 0.0);
						j++;
					}
				}
				this->Hal = Hal;
			}

			else if (key_word == "TEMPERATURE") {
				//cout << "Reading TEMPERATURE data" << '\n';
				getline(f, line);
				istringstream temperature_stream(line);
				float temperature;
				temperature_stream >> temperature;
				this->beta = cuCdivf(ONE, cuCmulf({ temperature, 0.0 }, boltz_k));
			}

			else if (key_word == "BATHTYPE") {
				getline(f, line);
				istringstream bathtype_stream(line);
				bathtype_stream >> this->bath_type;
			}

			else if (key_word == "BATH") {
				//cout << "Reading BATH data" << '\n';
				if (this->bath_type == "artificial") {
					getline(f, line);
					istringstream bath_stream(line);
					float l;
					bath_stream >> l;
					for (int i = 0; i < this->K; i++) this->lambda.push_back(make_cuComplex(l, 0.0));
					int total_size = this->sys_size * this->sys_size;
					vector<data_type> S(total_size * this->K, { 0.0, 0.0 });
					for (int i = 0; i < this->K; i++) {
						S[i * total_size + i + this->sys_size * i] = { 1.0, 0.0 };
						this->Hal[i + this->sys_size * i].x += l;
					}


					this->d_S = To_Device(S);
					this->d_Hal = To_Device(this->Hal);
					int n_modes;
					getline(f, line);
					istringstream n_modes_stream(line);
					n_modes_stream >> n_modes;

					vector<data_type> alpha;
					vector<data_type> alpha_t;
					vector<data_type> gamma;
					for (int i = 0; i < n_modes; i++) {
						double a, b, g, w;
						getline(f, line);
						istringstream modes_stream(line);
						modes_stream >> a >> b >> g >> w;
						alpha.push_back(make_cuComplex(a, b));
						alpha_t.push_back(make_cuComplex(a, -b));
						gamma.push_back(make_cuComplex(g, w));
					}
					for (int i = 0; i < this->K; i++) {
						this->alpha.push_back(alpha);
						this->alpha_t.push_back(alpha_t);
						this->gamma.push_back(gamma);
					}
					this->K_m = this->gamma[0].size();
				}
			}

			else if (key_word == "TIME") {
				//cout << "Reading TIME data" << '\n';
				getline(f, line);
				istringstream time_stream(line);
				time_stream >> this->t_start >> this->t_end >> this->step_size >> this->print_step;
			}
		}
		cublasError(cublasDestroy(cublasH), __FILE__, __LINE__);
	}
}

void param::reset_rho() {
	int total_size = this->sys_size * this->sys_size;
	if (!this->d_rho) {
		cudaError(cudaMalloc(&this->d_rho, this->ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
		for (int i = 0; i < this->ado.size(); i++) this->d_rho_ptr.push_back(this->d_rho + total_size * i);
	}
	cudaError(cudaMemset(this->d_rho, 0.0, this->ado.size() * total_size * sizeof(data_type)), __FILE__, __LINE__);
	data_type t = ONE;
	cudaMemcpy(this->d_rho, &t, sizeof(data_type), cudaMemcpyHostToDevice);
}

/*
void param::reset_rho() {
	int total_size = this->sys_size * this->sys_size;
	if (this->d_rho.size()) {
		for (int i = 0; i < this->ado.size(); i++) {
			data_type* tmp = nullptr;
			cudaError(cudaMalloc(&tmp, total_size * sizeof(data_type)), __FILE__, __LINE__);
			cudaError(cudaMemset(tmp, 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);
			this->d_rho.push_back(tmp);
		}
	}
	else {
		for (int i = 0; i < this->ado.size(); i++) {
			cudaError(cudaMemset(this->d_rho[i], 0.0, total_size * sizeof(data_type)), __FILE__, __LINE__);

		}
	}
	data_type* t = new data_type(ONE);
	cudaError(cudaMemcpy(this->d_rho[0], t, sizeof(data_type), cudaMemcpyHostToDevice), __FILE__, __LINE__);
	delete(t);
}
*/
void param::reset_rho_batch() {
	if (!this->d_rho_batch) cudaError(cudaMalloc(&this->d_rho_batch, this->ado.size() * sizeof(data_type*)), __FILE__, __LINE__);

	cudaError(cudaMemcpy(this->d_rho_batch, this->d_rho_ptr.data(), this->ado.size() * sizeof(data_type*), cudaMemcpyHostToDevice), __FILE__, __LINE__);
	cudaError(cudaDeviceSynchronize(), __FILE__, __LINE__);
}


void param::param_free() {
	if (this->d_Disorder) cudaError(cudaFree(this->d_Disorder), __FILE__, __LINE__);
	if (this->d_eigenstate) cudaError(cudaFree(this->d_eigenstate), __FILE__, __LINE__);
	if (this->d_Hal) cudaError(cudaFree(this->d_Hal), __FILE__, __LINE__);
	if (this->d_S) cudaError(cudaFree(this->d_S), __FILE__, __LINE__);
	//if (this->d_rho.size() > 0) for (int i = 0; i < this->ado.size(); i++) cudaError(cudaFree(this->d_rho[i]), __FILE__, __LINE__);
	if (this->d_rho) cudaError(cudaFree(this->d_rho), __FILE__, __LINE__);
	if (this->d_rho_batch) cudaError(cudaFree(this->d_rho_batch), __FILE__, __LINE__);
	if (this->d_Xx) cudaError(cudaFree(this->d_Xx), __FILE__, __LINE__);
	if (this->d_Xy) cudaError(cudaFree(this->d_Xy), __FILE__, __LINE__);
	if (this->d_Xz) cudaError(cudaFree(this->d_Xz), __FILE__, __LINE__);
	if (this->d_X_site) cudaError(cudaFree(this->d_X_site), __FILE__, __LINE__);
}
