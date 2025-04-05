#include "param.h"

double delta1(int a, int b) {
	return a == b ? 1 : 0;
}

void print_matrix_real1(const gsl_matrix_complex* M) {
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			gsl_complex ele = gsl_matrix_complex_get(M, i, j);
			printf("%.2f\t", GSL_REAL(ele));
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void param::param_2des(string filename) {

	fstream f(filename);
	string line, key_word;
	while (getline(f, line)) {
		if (line.empty() || line[0] == '#') continue;
		istringstream iss(line);
		iss >> key_word;

		if (key_word == "HEOM") {
			//cout << "Reading HEOM data" << '\n';

			int total_size = this->sys_size * this->sys_size;
			getline(f, line);
			istringstream val(line);
			val >> this->K >> this->L;
		}

		else if (key_word == "SIZE") {
			//cout << "Reading SIZE data" << '\n';
			getline(f, line);
			istringstream val(line);
			val >> this->single_size;
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
			//cout << "Reading DISORDER data" << '\n';
			getline(f, line);
			istringstream disorder_stream(line);
			disorder_stream >> this->n_sample >> this->seed;
			this->Disorder = gsl_matrix_complex_alloc(this->sys_size, this->sys_size);
			double val;
			for (int i = 0; i < this->sys_size; i++) {
				getline(f, line);
				istringstream row(line);
				int j = 0;
				while (row >> val) {
					gsl_matrix_complex_set(this->Disorder, i, j, gsl_complex_rect(val / h_bar, 0.0));
					j++;
				}
			}
		}

		else if (key_word == "HAMILTONIAN") {
			//cout << "Reading HAMILTONIAN data" << '\n';

			this->Hal = gsl_matrix_complex_alloc(this->sys_size, this->sys_size);
			gsl_matrix_complex_set_all(this->Hal, gsl_complex_rect(0.0, 0.0));
			double val;
			for (int i = 0; i < this->single_size; i++) {
				getline(f, line);
				istringstream row(line);
				int j = 0;
				while (row >> val) {
					gsl_matrix_complex_set(this->Hal, i, j, gsl_complex_rect(val / h_bar, 0.0));
					j++;
				}
			}
		}

		else if (key_word == "DIPOLE") {
			//cout << "Reading DIPOLE data" << '\n';
			getline(f, line);
			istringstream ndipole_stream(line);
			int ndipole;
			ndipole_stream >> ndipole;
			this->X_x = gsl_matrix_complex_alloc(this->sys_size, this->sys_size);
			gsl_matrix_complex_set_all(this->X_x, gsl_complex_rect(0.0, 0.0));
			this->X_y = gsl_matrix_complex_alloc(this->sys_size, this->sys_size);
			gsl_matrix_complex_set_all(this->X_y, gsl_complex_rect(0.0, 0.0));
			this->X_z = gsl_matrix_complex_alloc(this->sys_size, this->sys_size);
			gsl_matrix_complex_set_all(this->X_z, gsl_complex_rect(0.0, 0.0));


			for (int k = 0; k < ndipole; k++) {
				getline(f, line);
				istringstream dipole_stream(line);
				double x, y, z;
				int i, j;
				dipole_stream >> i >> j >> x >> y >> z;
				gsl_matrix_complex_set(this->X_x, i, j, gsl_complex_rect(x, 0.0));
				gsl_matrix_complex_set(this->X_y, i, j, gsl_complex_rect(y, 0.0));
				gsl_matrix_complex_set(this->X_z, i, j, gsl_complex_rect(z, 0.0));
			}

			this->X_site = gsl_matrix_complex_alloc(sys_size, sys_size);
			gsl_matrix_complex_set_all(this->X_site, gsl_complex_rect(0.0, 0.0));

			for (int i = 0; i < this->sys_size; i++) {
				for (int j = this->single_size; j < this->sys_size; j++) {
					int a = this->d_map[j].first;
					int b = this->d_map[j].second;
					double x1 = GSL_REAL(gsl_matrix_complex_get(this->X_x, 0, a));
					double x2 = GSL_REAL(gsl_matrix_complex_get(this->X_x, 0, b));
					double x = delta1(b, i) * x1 + delta1(a, i) * x2;
					gsl_matrix_complex_set(this->X_x, i, j, gsl_complex_rect(x, 0.0));
				}
			}
			for (int i = 0; i < this->sys_size; i++) {
				for (int j = this->single_size; j < this->sys_size; j++) {
					int a = this->d_map[j].first;
					int b = this->d_map[j].second;
					double x1 = GSL_REAL(gsl_matrix_complex_get(this->X_y, 0, a));
					double x2 = GSL_REAL(gsl_matrix_complex_get(this->X_y, 0, b));
					double x = delta1(b, i) * x1 + delta1(a, i) * x2;
					gsl_matrix_complex_set(this->X_y, i, j, gsl_complex_rect(x, 0.0));
				}
			}
			for (int i = 0; i < this->sys_size; i++) {
				for (int j = this->single_size; j < this->sys_size; j++) {
					int a = this->d_map[j].first;
					int b = this->d_map[j].second;
					double x1 = GSL_REAL(gsl_matrix_complex_get(this->X_z, 0, a));
					double x2 = GSL_REAL(gsl_matrix_complex_get(this->X_z, 0, b));
					double x = delta1(b, i) * x1 + delta1(a, i) * x2;
					gsl_matrix_complex_set(this->X_z, i, j, gsl_complex_rect(x, 0.0));
				}
			}

			for (int i = 0; i < this->sys_size; i++) {
				for (int j = 0; j < this->sys_size; j++) {
					gsl_complex ele_x, ele_y, ele_z;
					ele_x = gsl_matrix_complex_get(this->X_x, i, j);
					ele_y = gsl_matrix_complex_get(this->X_y, i, j);
					ele_z = gsl_matrix_complex_get(this->X_z, i, j);
					gsl_matrix_complex_set(this->X_site, i, j, gsl_complex_rect(sqrt(pow(ele_x.dat[0], 2) + pow(ele_y.dat[0], 2) + pow(ele_z.dat[0], 2)), 0.0));
				}
			}
		}

		else if (key_word == "POLARIZATION") {
			getline(f, line);
			istringstream polarization_stream(line);
			vector<double> angle(4, 0.0);
			polarization_stream >> angle[0] >> angle[1] >> angle[2] >> angle[3];
			for (int i = 0; i < 4; i++) this->angle.push_back(angle[i] * pi / 180);
			for (int i = 0; i < 4; i++) {
				gsl_matrix_complex* m = gsl_matrix_complex_alloc(this->sys_size, this->sys_size);
				gsl_matrix_complex_set_all(m, gsl_complex_rect(0.0, 0.0));
				this->X.push_back(m);
			}
			this->polar = gsl_matrix_complex_alloc(3, 4);
		}

		else if (key_word == "TEMPERATURE") {
			//cout << "Reading TEMPERATURE data" << '\n';
			getline(f, line);
			istringstream temperature_stream(line);
			double temperature;
			temperature_stream >> temperature;
			this->beta = 1 / (temperature * boltz_k);
		}

		else if (key_word == "BATH") {
			//cout << "Reading BATH data" << '\n';
			getline(f, line);
			istringstream bath_stream(line);
			double l;
			bath_stream >> l;
			for (int i = 0; i < this->K; i++) this->lambda.push_back(l);
			for (int i = 0; i < this->K; i++) gsl_matrix_complex_set(this->Hal, i + 1, i + 1, gsl_complex_rect(gsl_matrix_complex_get(this->Hal, i + 1, i + 1).dat[0] + this->lambda[i], 0.0));
			for (int i = 0; i < this->K; i++) {
				gsl_matrix_complex* m = gsl_matrix_complex_alloc(sys_size, sys_size);
				gsl_matrix_complex_set_all(m, gsl_complex_rect(0.0, 0.0));
				gsl_matrix_complex_set(m, i + 1, i + 1, gsl_complex_rect(1.0, 0.0));
				S.push_back(m);
			}
			int n_modes;
			getline(f, line);
			istringstream n_modes_stream(line);
			n_modes_stream >> n_modes;
			vector<gsl_complex> alpha;
			vector<gsl_complex> alpha_t;
			vector<gsl_complex> gamma;
			for (int i = 0; i < n_modes; i++) {
				double a, b, g, w;
				getline(f, line);
				istringstream modes_stream(line);
				modes_stream >> a >> b >> g >> w;
				alpha.push_back(gsl_complex_rect(a, b));
				alpha_t.push_back(gsl_complex_rect(a, -b));
				gamma.push_back(gsl_complex_rect(g, w));
			}
			for (int i = 0; i < this->K; i++) {
				this->alpha.push_back(alpha);
				this->alpha_t.push_back(alpha_t);
				this->gamma.push_back(gamma);
			}
			this->K_m = this->gamma[0].size();
			
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
				this->pulses[i].w0 = this->pulses[i].w0 / h_bar; //cm^-1 to fs^-1
				this->pulses[i].E0 = this->pulses[i].E0 / h_bar; //cm^-1 to fs^-1
			}
		}

		else if (key_word == "TIME") {
			//cout << "Reading TIME data" << '\n';
			getline(f, line);
			istringstream time_stream(line);
			time_stream >> this->t_start >> this->t_end >> this->step_size >> this->print_step;
		}

	}
}

void param::param_dynamics(string filename) {

	fstream f(filename);
	string line, key_word;
	while (getline(f, line)) {
		if (line.empty() || line[0] == '#') continue;
		istringstream iss(line);
		iss >> key_word;

		if (key_word == "HEOM") {
			//cout << "Reading HEOM data" << '\n';

			int total_size = this->sys_size * this->sys_size;
			getline(f, line);
			istringstream val(line);
			val >> this->K >> this->L;
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

			this->Hal = gsl_matrix_complex_alloc(this->sys_size, this->sys_size);
			gsl_matrix_complex_set_all(this->Hal, gsl_complex_rect(0.0, 0.0));
			double val;
			for (int i = 0; i < this->single_size; i++) {
				getline(f, line);
				istringstream row(line);
				int j = 0;
				while (row >> val) {
					gsl_matrix_complex_set(this->Hal, i, j, gsl_complex_rect(val / h_bar, 0.0));
					j++;
				}
			}
		}

		else if (key_word == "TEMPERATURE") {
			//cout << "Reading TEMPERATURE data" << '\n';
			getline(f, line);
			istringstream temperature_stream(line);
			double temperature;
			temperature_stream >> temperature;
			this->beta = 1 / (temperature * boltz_k);
		}

		else if (key_word == "BATH") {
			//cout << "Reading BATH data" << '\n';
			getline(f, line);
			istringstream bath_stream(line);
			double l;
			bath_stream >> l;
			for (int i = 0; i < this->K; i++) this->lambda.push_back(l);
			for (int i = 0; i < this->K; i++) gsl_matrix_complex_set(this->Hal, i, i, gsl_complex_rect(gsl_matrix_complex_get(this->Hal, i, i).dat[0] + this->lambda[i], 0.0));
			for (int i = 0; i < this->K; i++) {
				gsl_matrix_complex* m = gsl_matrix_complex_alloc(sys_size, sys_size);
				gsl_matrix_complex_set_all(m, gsl_complex_rect(0.0, 0.0));
				gsl_matrix_complex_set(m, i, i, gsl_complex_rect(1.0, 0.0));
				S.push_back(m);
			}
			int n_modes;
			getline(f, line);
			istringstream n_modes_stream(line);
			n_modes_stream >> n_modes;
			vector<gsl_complex> alpha;
			vector<gsl_complex> alpha_t;
			vector<gsl_complex> gamma;
			for (int i = 0; i < n_modes; i++) {
				double a, b, g, w;
				getline(f, line);
				istringstream modes_stream(line);
				modes_stream >> a >> b >> g >> w;
				alpha.push_back(gsl_complex_rect(a, b));
				alpha_t.push_back(gsl_complex_rect(a, -b));
				gamma.push_back(gsl_complex_rect(g, w));
			}
			for (int i = 0; i < this->K; i++) {
				this->alpha.push_back(alpha);
				this->alpha_t.push_back(alpha_t);
				this->gamma.push_back(gamma);
			}
			this->K_m = this->gamma[0].size();
		}

		else if (key_word == "TIME") {
			//cout << "Reading TIME data" << '\n';
			getline(f, line);
			istringstream time_stream(line);
			time_stream >> this->t_start >> this->t_end >> this->step_size >> this->print_step;
		}

	}
}

void param::reset_rho() {
	for (int i = 0; i < this->ado.size(); i++) gsl_matrix_complex_set_all(this->rho[i], gsl_complex_rect(0.0, 0.0));
	gsl_matrix_complex_set(this->rho[0], 0, 0, gsl_complex_rect(1.0, 0.0));
}



void param::param_free() {
	size_t ado_size = this->ado.size();
	int K = this->K;
	for (int i = 0; i < ado_size; i++) if(this->rho[i]) gsl_matrix_complex_free(this->rho[i]);
	for (int i = 0; i < K; i++) if(this->S[i]) gsl_matrix_complex_free(this->S[i]);
	if(!this->X.empty()) for (int i = 0; i < 4; i++) gsl_matrix_complex_free(this->X[i]);
	if(this->Hal) gsl_matrix_complex_free(this->Hal);
	if(this->X_site) gsl_matrix_complex_free(this->X_site);
	if(this->X_x) gsl_matrix_complex_free(this->X_x);
	if(this->X_y) gsl_matrix_complex_free(this->X_y);
	if(this->X_z) gsl_matrix_complex_free(this->X_z);
	if(this->polar) gsl_matrix_complex_free(this->polar);
	if(this->Disorder) gsl_matrix_complex_free(this->Disorder);
}
