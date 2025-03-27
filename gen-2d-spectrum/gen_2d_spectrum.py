import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import griddata, RegularGridInterpolator
from scipy.fft import fft2, fftshift, fftfreq
from scipy.ndimage import gaussian_filter
        
size = 2048
size1 = 1000
size2 = 61

time = 1000
t_size = 0
tau_size = 0
t_step = 10
tau_step = 10
T = 0
tau_start = -600
tau_end = 600

non_rephaasing_tau_start = tau_start
non_rephaasing_tau_end = 0
non_rephaasing_data = np.zeros((size1, size2), dtype=np.complex128)  

#loading non-rephasing data
for i, tau in enumerate(np.arange(non_rephaasing_tau_start, non_rephaasing_tau_end + 1, tau_step)):
    filename = f"./2d-output/out_{tau}_{T}.out"
    d = []
    with open(filename, 'r') as file:
        # Convert the file content into an array of complex numbers
        complex_numbers = np.loadtxt(filename, dtype=complex, converters={0: lambda x: complex(float(x), 0), 1: lambda x: complex(0, float(x))})

        # Combine the real and imaginary parts into complex numbers 
        d = np.array([complex(r.real, i.imag) for r, i in complex_numbers])
        d = d * 1j
        for j in range(size1):
            if j < len(d): non_rephaasing_data[j][i] = d[j]
            else: non_rephaasing_data[j][i] = 0

            
    if(i == 0): t_size = len(d)
    tau_size += 1
    
rephaasing_data = np.zeros((size1, size2), dtype=np.complex128) 
rephasing_tau_start = 0
rephasing_tau_end = tau_end

#loading rephasing data
for i, tau in enumerate(np.arange(rephasing_tau_start, rephasing_tau_end + 1, tau_step)):
    filename = f"./2d-output/out_{tau}_{T}.out"
    d = []
    with open(filename, 'r') as file:
        # Convert the file content into an array of complex numbers
        complex_numbers = np.loadtxt(filename, dtype=complex, converters={0: lambda x: complex(float(x), 0), 1: lambda x: complex(0, float(x))})

        # Combine the real and imaginary parts into complex numbers
        d = np.array([complex(r.real, i.imag) for r, i in complex_numbers])
        d = d * 1j
        for j in range(size1):
            if j < len(d): rephaasing_data[j][i] = d[j]
            else: rephaasing_data[j][i] = 0

#2DFFT
fft_non_rephaasing_data = fftshift(fft2(non_rephaasing_data)).real
fft_rephaasing_data = fftshift(fft2(rephaasing_data)).real
fft_data = fft_rephaasing_data + fft_non_rephaasing_data


# Apply the Gaussian filter
fft_data = gaussian_filter(fft_data, sigma=(8, 0.8))  # Adjust sigma as needed


'''
#arcsinh signal amplifier, use as need
# arcsinh scaling
fft_data = fft_data / np.max(np.abs(fft_data)) * 3
fft_data = np.log(fft_data + np.sqrt(1 + fft_data ** 2))
'''

# set inpolation gactor
interp_factor = 2
new_size1 = fft_data.shape[0] * interp_factor
new_size2 = fft_data.shape[1] * interp_factor

# produce 2D grid of origninal data
x = np.arange(fft_data.shape[0])
y = np.arange(fft_data.shape[1])

# produce new grid of interpolated data
new_x = np.linspace(0, fft_data.shape[0] - 1, new_size1)
new_y = np.linspace(0, fft_data.shape[1] - 1, new_size2)

# use RegularGridInterpolator to do bilinear interpolation
interpolator = RegularGridInterpolator((x, y), fft_data)

# produce new points on new grid
new_grid = np.meshgrid(new_x, new_y, indexing='ij')
new_points = np.array([new_grid[0].ravel(), new_grid[1].ravel()]).T

# do interpolation on new points
interpolated_data = interpolator(new_points).reshape(new_size1, new_size2)
fft_data = interpolated_data


t = np.linspace(0, time + (new_size1 - t_size) * t_step, new_size1)
tau = np.linspace(0, tau_end + (new_size2 - tau_size) * tau_step, new_size2)

t = 5308 * 2 * np.pi * fftshift(fftfreq(len(t), t[1] - t[0]))
tau = 5308 * 2 * np.pi * fftshift(fftfreq(len(tau), tau[1] - tau[0]))

# draw high resolution contour map
# interpolating original axises onto a high-resolution grid
high_res_factor = 12
tau_high_res = np.linspace(tau.min(), tau.max(), size2 * high_res_factor)
t_high_res = np.linspace(t.min(), t.max(), size1 * high_res_factor)
X_high_res, Y_high_res = np.meshgrid(tau_high_res, t_high_res)

# interpolating original data onto a high-resolution grid
X, Y = np.meshgrid(tau, t)
points = np.array([X.flatten(), Y.flatten()]).T
values = fft_data.flatten()
fft_data_high_res = griddata(points, values, (X_high_res, Y_high_res), method='cubic')
fft_data = fft_data_high_res

# signal normalization
fft_data = fft_data / np.max(np.abs(fft_data))

# plotting a smooth contour map
lv = np.linspace(np.min(fft_data_high_res), np.max(fft_data_high_res), 10)
norm = Normalize(vmin=-1.0, vmax=1.0)
plt.contour(X_high_res, Y_high_res[::-1], fft_data_high_res, colors=['#000', '#000'], levels=lv, norm=norm, linewidths=1)
img = plt.imshow(fft_data, extent=(tau.min(), tau.max(), t.min(), t.max()), cmap='jet', norm = norm)
plt.colorbar(img, label='Magnitude', ticks=[-1.0, -0.5, 0, 0.5, 1.0])
plt.xlabel("$\omega_{\\tau} \ cm^{-1}$")
plt.ylabel("$\omega_t \ cm^{-1}$")
plt.xlim(0, 900)
plt.ylim(0, 900)
plt.title(f'pulse width = 100 fs')
plt.savefig('my_spectrum.png')
#show the figure as need
#plt.show()
