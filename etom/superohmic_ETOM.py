import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit

# Define constants
hbar = 5308              # Planck constant (cm^-1 * fs)
w_c = 100 / hbar         # Cutoff frequency (fs^-1)
coupling_str = 1.0       # Coupling strength (dimensionless)
kt = 0.5 * w_c           # Thermal energy (fs^-1)
beta = 1 / kt            # Inverse temperature (fs)
N = 3                    # Number of effective oscillators
init_guess = 0.01        # Initial guess for fit parameters
save_fig = False      

# Define spectral density
# You can change different J(w) to meet your work
def J(w):
    """
    Defines the spectral density J(w) = w^3 / w_c^2 * exp(-w / w_c).
    """
    return w**3 / w_c**2 * np.exp(-w / w_c)

def integrand_real(w, t):
    return J(w) * np.cos(w * t) / np.tanh(0.5 * beta * w)

def integrand_imag(w, t):
    return -J(w) * np.sin(w * t)

def original_function_real(t):
    integral, _ = quad(integrand_real, 0, np.inf, args=(t,), epsabs=1e-10, epsrel=1e-10, limit=10_000_000)
    return integral / np.pi

def original_function_imag(t):
    integral, _ = quad(integrand_imag, 0, np.inf, args=(t,), epsabs=1e-10, epsrel=1e-10, limit=10_000_000)
    return integral / np.pi

# Define fit function using multiple bases
def fit_function_real(t, *params):
    """
    Fitting function for the real part of the correlation function.
    We use a sum of N damped cosines:
        sum_{i=1 to N} [a_i cos(w_i t) exp(-g_i t)]
    
    params is expected to contain 3*N parameters in the order:
    (a_0, g_0, w_0, a_1, g_1, w_1, ..., a_(N-1), g_(N-1), w_(N-1)).
    """
    result = 0.0
    for j in range(N):
        a = params[3 * j]
        g = params[3 * j + 1]
        w = params[3 * j + 2]
        result += a * np.cos(w * t) * np.exp(-g * t)
    return result

def fit_function_imag(t, *params):
    """
    Fitting function for the imaginary part of the correlation function.
    We reuse (g, w) from the real part fit and only fit b_i for amplitude:
        sum_{i=1 to N} [b_i cos(w_i t) exp(-g_i t)]
    
    params is expected to contain N parameters (b_0, b_1, ..., b_(N-1)).
    g_values and w_values are taken from the real part fit.
    """
    result = 0.0
    for j in range(N):
        b = params[j]
        g = g_values[j]
        w = w_values[j]
        result += b * np.cos(w * t) * np.exp(-g * t)
    return result

# Generate sample data
t_values = np.linspace(0, 500, 1000)
y_values_real = np.array([original_function_real(t) for t in t_values])
y_values_imag = np.array([original_function_imag(t) for t in t_values])

# Perform the fit for the real part
initial_guess_real = [init_guess, init_guess, init_guess] * N  # (a_0, g_0, w_0, ..., a_(N-1), g_(N-1), w_(N-1))
popt_real, pcov_real = curve_fit(
    fit_function_real, t_values, y_values_real,
    p0=initial_guess_real, maxfev=1_000_000_000
)

# Extract (g_i, w_i) from the real-part fit
g_values = popt_real[1::3]  # every 3rd element, starting from index 1
w_values = popt_real[2::3]  # every 3rd element, starting from index 2

# Perform the fit for the imaginary part
initial_guess_imag = [init_guess] * N
popt_imag, pcov_imag = curve_fit(
    fit_function_imag, t_values, y_values_imag,
    p0=initial_guess_imag, maxfev=1_000_000_000,
    bounds=(-np.inf, np.inf)
)


# Print the ETOM parameters
for j in range(N):
    a = popt_real[3 * j]
    b = popt_imag[j]
    g = popt_real[3 * j + 1]
    w = popt_real[3 * j + 2]
    
    # Print pairs for positive and negative frequencies
    print(f"{coupling_str * a / 2:.8f} {coupling_str * b / 2:.8f} {g:.8f} {w:.8f}")
    print(f"{coupling_str * a / 2:.8f} {coupling_str * b / 2:.8f} {g:.8f} {-w:.8f}")

# Plot the Real Part

fontsize = 30
labelpad = 12
labelsize = 16

def format_func(value, tick_number):
    return f"{value:.2e}"

# Plot the Real Part
plt.figure(figsize=(10, 8))
plt.scatter(
    t_values[::20] * w_c,
    y_values_real[::20],
    label=fr'Analytic $(k_bT = {kt / w_c}\omega_c)$',
    s=10,
    color='black'
)
plt.plot(
    t_values * w_c,
    fit_function_real(t_values, *popt_real),
    label=fr'ETOM $(k_bT = {kt / w_c}\omega_c)$',
    color='blue'
)
plt.xlabel(r"$t\ (1 / \omega_c)$", fontsize=fontsize, labelpad=labelpad)
plt.ylabel(r"$\mathrm{Re}\{C(t)\} \,/\, \eta$", fontsize=fontsize, labelpad=labelpad)
plt.legend(fontsize=fontsize - 5)
plt.tick_params(axis='both', which='major', labelsize=labelsize)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
if save_fig: plt.savefig('superohmic_TCF_real.png')
plt.show()

# Plot the Imaginary Part
plt.figure(figsize=(10, 8))
plt.scatter(
    t_values[::20] * w_c,
    y_values_imag[::20],
    label=fr'Analytic $(k_bT = {kt / w_c}\omega_c)$',
    s=10,
    color='black'
)
plt.plot(
    t_values * w_c,
    fit_function_imag(t_values, *popt_imag),
    label=fr'ETOM $(k_bT = {kt / w_c}\omega_c)$',
    color='blue'
)
plt.xlabel(r"$t\ (1 / \omega_c)$", fontsize=fontsize, labelpad=labelpad)
plt.ylabel(r"$\mathrm{Im}\{C(t)\} \,/\, \eta$", fontsize=fontsize, labelpad=labelpad)
plt.legend(fontsize=fontsize - 5)
plt.tick_params(axis='both', which='major', labelsize=labelsize)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
if save_fig: plt.savefig('superohmic_TCF_imag.png')
plt.show()


import os

# Assume we already have:
#   l: float, the reorganization energy 
#   N: int, the number of effective oscillators
#   popt_real: array from curve_fit (includes a_j, g_j, w_j for each oscillator)
#   popt_imag: array from curve_fit (includes b_j for each oscillator)

# Prepare the lines to be inserted
new_data_lines = []

# First line: reorganization energy
new_data_lines.append(f"{6 * coupling_str * w_c / np.pi}\n")

# Second line: number of ETOM modes
new_data_lines.append(f"{2 * N}\n")

# Next pairs of lines:
#   (a/2) (b/2) g w
#   (a/2) (b/2) g -w
for j in range(N):
    a = popt_real[3 * j]
    b = popt_imag[j]
    g = popt_real[3 * j + 1]
    w = popt_real[3 * j + 2]

    new_data_lines.append(f"{coupling_str * a / 2:.8f} {coupling_str * b / 2:.8f} {g:.8f} {w:.8f}\n")
    new_data_lines.append(f"{coupling_str * a / 2:.8f} {coupling_str * b / 2:.8f} {g:.8f} {-w:.8f}\n")

# Read the existing file
key_file_path = "./key.key-tmpl"
with open(key_file_path, "r") as f:
    lines = f.readlines()

# Locate the lines containing BATH and DIPOLE
bath_index = None
dipole_index = None

for i, line in enumerate(lines):
    content = line.strip()
    if content == "BATH" and bath_index is None:
        bath_index = i
    elif content == "DIPOLE" and dipole_index is None:
        dipole_index = i

# Check if we found both keywords and if BATH appears before DIPOLE
if bath_index is None or dipole_index is None or bath_index >= dipole_index:
    print("Error: Could not find correct BATH and DIPOLE order in the file.")
else:
    # 4) Remove everything between BATH and DIPOLE and insert new data
    
    # Keep lines from the start of the file up to BATH (including BATH)
    updated_lines = lines[:bath_index+1]
    
    # Insert our new data lines
    updated_lines.extend(new_data_lines)
    updated_lines.append("\n")
    
    # Keep DIPOLE and everything after it
    updated_lines.append(lines[dipole_index])
    updated_lines.extend(lines[dipole_index+1:])
    
    # 5) Write back to the file
    with open(key_file_path, "w") as f:
        f.writelines(updated_lines)
    
    print("Successfully updated key.key-tmpl! Content between BATH and DIPOLE removed, and new data inserted.")