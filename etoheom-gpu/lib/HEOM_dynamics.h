#include "HEOM_utilize.h"
#include "HEOM_TD_hamiltonian.h"
#include "HEOM_kernel_func.h"
#include "HEOM_param.h"
#include "HEOM_polar.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>

#ifndef HEOM_DYNAMICS_H
#define HEOM_DYNAMICS_H

void construct_Hal(param& key, vector<data_type>& Hal);

void Build_ADO(param& key, std::string& current, const int cur_L);

void Build_ADO_map(param& key);

void construct_ADO_set(param& key);

unordered_map<string, int> Build_ADO_map(const vector<string> ado);



void total_ADO_dynamics(cublasHandle_t& cublasH, const param& key, const data_type* d_rho_copy, data_type* d_drho);


void total_ADO_dynamics_Ht(cublasHandle_t& cublasH, param& key, const data_type* d_rho_copy, data_type* d_drho, const data_type* d_Ht);

void total_ADO_dynamics_Ht_batch(cublasHandle_t& cublasH, param& key, const data_type* d_rho_copy, data_type* d_drho, const data_type* d_Ht, const data_type* const* d_rho_batch, data_type** d_drho_batch, data_type* d_work_space, data_type** d_work_batch, data_type* d_rho_work_space, data_type** d_rho_work_batch, vector<cudaStream_t>& streams);


void dynamics_solver(param& key, std::vector<int>& sites, vector<vector<float>>& population);



void twoD_spectrum_solver(param& key, const data_type* d_H, const int nv1, const int nv2, const int nv3, std::vector<data_type>& polarization);

void propagation_Ht_batch(param& key, const data_type* d_H, const int nv1, const int nv2, const int nv3, std::vector<data_type>& polarization);


#endif
