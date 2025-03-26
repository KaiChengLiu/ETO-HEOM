#include "HEOM_kernel_func.h"

__global__ void mul_int_scalarKernel(const int a, const data_type* b, data_type* c, const int total_size) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < total_size) {
		c[idx].x = b[idx].x * a;
		c[idx].y = b[idx].y * a;
	} 
}

__global__ void mul_double_scalarKernel(const double a, const data_type* b, data_type* c, const int total_size) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < total_size) {
		c[idx].x = b[idx].x * a;
		c[idx].y = b[idx].y * a;
	}
}

__global__ void init_cuComplex(data_type* data, const int size, data_type val) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size) {
		data[idx] = val;
	}
}

__global__ void concatenateArrays(const data_type* arrays, data_type* result, const int i, const int sys_size) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < sys_size * sys_size) {
		result[i * sys_size * sys_size + idx] = arrays[idx];
	}
}

__global__ void copyUpperToLower(data_type* matrix, const int sys_size) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row < sys_size && col < sys_size && col > row) {
		matrix[row * sys_size + col].x = matrix[col * sys_size + row].x;
		matrix[row * sys_size + col].y = matrix[col * sys_size + row].y;
	}
}
__global__ void elementwise_pow_kernel(data_type* arrays, const double p, const int size) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < size) {
		arrays[idx].x = pow(arrays[idx].x, p);
		arrays[idx].y = pow(arrays[idx].y, p);
	}
}

__global__ void transpose(const data_type* in, data_type* out, const int sys_size) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row < sys_size && col < sys_size) {
		out[col * sys_size + row].x = in[row * sys_size + col].x;
		out[col * sys_size + row].y = in[row * sys_size + col].y;
	}
}

__global__ void inplace_transpose(data_type* matrix, int sys_size) {
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int row = blockIdx.y * blockDim.y + threadIdx.y;

	if (col < sys_size && row < sys_size && col > row) {
		data_type temp = matrix[row * sys_size + col];
		matrix[row * sys_size + col] = matrix[col * sys_size + row];
		matrix[col * sys_size + row] = temp;
		
	}
}

__global__ void build_I(data_type* d_I, int sys_size) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	if (row < sys_size && col < sys_size) {
		if (row == col) {
			d_I[row * sys_size + col] = make_cuComplex(1.0f, 0.0f);
		}
		else {
			d_I[row * sys_size + col] = make_cuComplex(0.0f, 0.0f);
		}
	}
}

__global__ void trace(data_type* d_A, int sys_size, data_type* d_res) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	if (row < sys_size && col < sys_size && row == col) {
		d_res[row].x = d_A[row * sys_size + col].x;
		d_res[row].y = d_A[row * sys_size + col].y;
	}
}

__global__ void copyKernel(data_type* d_work, const data_type* d_V, int total_size, int batch_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < batch_size) {
        // Calculate the start position in d_work for each batch
        data_type* dst = d_work + idx * total_size;
        // Copy data from d_V to the corresponding position in d_work
        for (int j = 0; j < total_size; ++j) {
            dst[j] = d_V[j];
        }
    }
}
