import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit

# Define constants
hbar = 5308
w_c = 120 / hbar  # Adjust based on your requirements
coupling_str = 1
temp = 77
beta = 1 / 0.695 / temp  # Adjust based on your requirements
N = 1  # Number of bases
init_guess = 1


# Define spectral density
# You can change different J(w) to meet your work
def J(w):
  return coupling_str * w * np.exp(-w / w_c)

# Define integrand functions
def integrand_real(w, t):
    return J(w) * np.cos(w * t) / np.tanh(0.5 * beta * hbar * w)

def integrand_imag(w, t):
    return -J(w) * np.sin(w * t)

# Define original function (result of integration)
def original_function_real(t):
    integral, _ = quad(integrand_real, 0, np.inf, args=(t,))
    return integral

def original_function_imag(t):
    integral, _ = quad(integrand_imag, 0, np.inf, args=(t,))
    return integral

# Define fit function using multiple bases
def fit_function_real(t, *params):
    result = 0.0
    for i in range(N):
        a = params[3 * i]
        g = params[3 * i + 1]
        w = params[3 * i + 2]
        result += a * np.cos(w * t) * np.exp(-g * t)
    return result

t_values = np.linspace(0, 1000, 5000)
y_values_real = np.array([original_function_real(t) for t in t_values])

# Perform curve fitting
initial_guess_real = [init_guess, init_guess, init_guess] * N  # Initial parameter guesses
initial_guess_imag = [init_guess] * N 

popt_real, pcov_real = curve_fit(fit_function_real, t_values, y_values_real, p0=initial_guess_real, maxfev=10000000)

g_values = popt_real[1::3]
w_values = popt_real[2::3]

def fit_function_imag(t, *params):
    result = 0.0
    for i in range(N):
        a = params[i]
        w = w_values[i]
        g = g_values[i]  # Use the fitted g values
        result += a * np.cos(w * t) * np.exp(-g * t)
    return result

# Generate sample data
y_values_imag = np.array([original_function_imag(t) for t in t_values])

popt_imag, pcov_imag = curve_fit(fit_function_imag, t_values, y_values_imag, p0=initial_guess_imag, maxfev=10000000)

#for i in range(N):
#    print(f'a = {popt_real[3 * i] / 2}')
#    print(f'b = {popt_imag[i] / 2}')
#    print(f'g = {popt_real[3 * i + 1]}')
#    print(f'w = {popt_real[3 * i + 2]}')

import re
# Read the content of the key.txt file
key_file_path = "../2d_input/key.key-tmpl"

with open(key_file_path, 'r') as file:
    lines = file.readlines()

def find_keyword_index(keyword, lines):
    pattern = rf'^{keyword}\b'
    for i, line in enumerate(lines):
        if re.search(pattern, line):
            return i
    return None

# Find the index of the line containing the keyword 'BATHTYPE'
bathtype_index = find_keyword_index('BATHTYPE', lines)
if bathtype_index is not None:
    # Overwrite 'etom' after the 'BATHTYPE' line if something exists
    lines[bathtype_index + 1] = "etom\n"
else:
    print("Keyword 'BATHTYPE' not found in the file.")

# Find the index of the line containing the keyword 'BATH'
bath_index = find_keyword_index('BATH', lines)
if bath_index is not None:
    # Prepare the popt_real and popt_imag data
    bath_data = []
    bath_data.append(f"{w_c * coupling_str / np.pi:.6f}\n")
    bath_data.append(f"{N * 2}\n")
    for i in range(N):
        a = popt_real[3 * i]
        b = popt_imag[i]
        g = popt_real[3 * i + 1]
        w = popt_real[3 * i + 2]
        bath_data.append(f"{a / 2:.8f} {-b / 2:.8f} {g:.8f} {-w:.8f}\n")
        bath_data.append(f"{a / 2:.8f} {b / 2:.8f} {g:.8f} {w:.8f}\n")
    
    # Overwrite all bath_data after the 'BATH' line
    for i, data in enumerate(bath_data):
        lines[bath_index + 1 + i] = data
else:
    print("Keyword 'BATH' not found in the file.")

# Write the modified content back to the key.txt file
with open(key_file_path, 'w') as file:
    file.writelines(lines)
