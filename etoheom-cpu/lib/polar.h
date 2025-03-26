#ifndef POLAR_H
#define POLAR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include "constant.h"
#include "param.h"

using namespace std;

/*The array size is fixed to 3*/
void vector3_norm(double* vec);

/*Polar is a 3 * 4 matrix*/
void polar_mat_set(param& key);

void polar_mat_ranrot(param& key, gsl_rng* r);

void compute_pulse_interaction(param& key);

#endif POLAR_H
