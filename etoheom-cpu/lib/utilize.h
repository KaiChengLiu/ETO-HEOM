#ifndef UTILIZE_H
#define UTILIZE_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>

#include "param.h"
#include "constant.h"

void Build_ADO(param& key, std::string& current, const int cur_L);

void Build_ADO_map(param& key);

void construct_ADO_set(param& key);

void print_matrix_real(const gsl_matrix_complex* M);

void print_Hal(const gsl_matrix_complex* M);

void print_matrix(const gsl_matrix* M);

double delta(int a, int b);

#endif 
