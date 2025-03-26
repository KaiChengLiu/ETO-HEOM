#ifndef TD_HAMILTONIAN_H
#define TD_HAMILTONIAN_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "param.h"


void build_V_matrix(param& key, const gsl_matrix_complex* X, gsl_matrix_complex* V, const float t, const int idx);

void build_V_dagger_matrix(param& key, const gsl_matrix_complex* X, gsl_matrix_complex* V, const float t, const int idx);

#endif