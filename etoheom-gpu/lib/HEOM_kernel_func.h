#ifndef HEOM_KERNEL_FUNC_H
#define HEOM_KERNEL_FUNC_H
#include "HEOM_constant.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

__global__ void mul_int_scalarKernel(const int a, const data_type* b, data_type* c, const int total_size);

__global__ void mul_double_scalarKernel(const double a, const data_type* b, data_type* c, const int total_size);

__global__ void init_cuComplex(data_type* data, const int total_size, data_type val);

__global__ void concatenateArrays(const data_type* arrays, data_type* result, const int i, const int sys_size);

__global__ void copyUpperToLower(data_type* matrix, const int sys_size);

__global__ void elementwise_pow_kernel(data_type* arrays, const double p, const int size);

__global__ void transpose(const data_type* in, data_type* out, const int sys_size);

__global__ void inplace_transpose(data_type* matrix, int sys_size);

__global__ void build_I(data_type* d_I, int sys_size);

__global__ void trace(data_type* d_A, int sys_size, data_type* res);

__global__ void copyKernel(data_type* d_work, const data_type* d_V, int total_size, int batch_size);

#endif