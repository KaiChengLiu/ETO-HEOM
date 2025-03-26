#include "HEOM_utilize.h"
#include "HEOM_kernel_func.h"
#include "HEOM_param.h"
#include "cusolverDn.h"
#include <iostream>
#include <string>
#include <cmath>
#include <complex>

#ifndef HEOM_TD_HAMILTONIAN_H
#define HEOM_TD_HAMILTONIAN_H

void build_V_matrix(cublasHandle_t cublasH, param& key, const data_type* d_X, data_type* d_V, const float t, const int idx);

void build_V_dagger_matrix(cublasHandle_t cublasH, param& key, const data_type* d_X, data_type* d_V, const float t, const int idx);

#endif