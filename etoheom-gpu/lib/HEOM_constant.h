#ifndef HEOM_CONSTANT_H
#define HEOM_CONSTANT_H
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

using data_type = cuComplex;
extern const int ThreadPerBlock;
extern const float pi;

extern const data_type ONE;
extern const data_type ZERO;
extern const data_type MINUS_ONE;
extern const data_type MINUS_ICNT;
extern const data_type ICNT;
extern const data_type boltz_k; // cm^-1 / K
extern const data_type h_bar; // cm^-1 * fs
#endif
