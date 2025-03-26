#ifndef POLAR_H
#define POLAR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

#include "HEOM_constant.h"
#include "HEOM_param.h"
#include "HEOM_utilize.h"

using namespace std;

/*The array size is fixed to 3*/
void vector3_norm(float* vec);

/*Polar is a 3 * 4 matrix*/
void polar_mat_set(param& key);

void polar_mat_ranrot(param& key);

void compute_pulse_interaction(param& key, vector<data_type*>& M);

#endif
