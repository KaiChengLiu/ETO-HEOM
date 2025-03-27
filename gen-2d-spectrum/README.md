# Generate 2D Spectrum & Visualization Tool

This Python script processes simulated 2D electronic spectra from time-domain outputs, performs Fourier transforms, interpolates, and visualizes the frequency-domain 2D spectrum.

---

## Overview

The script:

1. Loads **non-rephasing** and **rephasing** time-domain data from `./2d_output/out_{tau}_{T}.out` files.
2. Constructs a 2D time-domain signal matrix for each.
3. Performs 2D FFT (Fast Fourier Transform) and combines both signals.
4. Applies optional smoothing (Gaussian filter) or signal amplification (arcsinh).
5. Interpolates to a higher resolution using bilinear interpolation.
6. Converts time axes to frequency axes using `fftshift` and `fftfreq`.
7. Re-interpolates the frequency-domain data to a higher resolution grid.
8. Normalizes and plots a smooth contour map of the 2D spectrum.
9. Saves the spectrum as `my_spectrum.{format}`.

---

## Adjustable Parameters

You can modify the following parameters to control the behavior and output:

```python
size = 2048        # Unused placeholder (can be removed or re-purposed)
size1 = 1000       # Number of time points (t-axis)
size2 = 61         # Number of delay points (tau-axis)
time = 1000        # Total simulation time
t0 = 0             # Initial pulse time offset (not used directly)
t_step = 10        # Time step size (fs)
tau_step = 10      # Tau step size (fs)
T = 0              # Population time (indexed in filename)

tau_start = -600   # Start of tau window (fs)
tau_end = 600      # End of tau window (fs)
```

### Interpolation

- `interp_factor`: Interpolation factor for increasing resolution before frequency transform (default = 2)
- `high_res_factor`: Interpolation factor after FFT for final plotting (default = 12)

### Optional Filters (can be uncommented as needed)

- **Gaussian filter** to smooth data:

```python
fft_data = gaussian_filter(fft_data, sigma=(8, 0.8))
```

- **Arcsinh signal amplifier** for enhancing dynamic range:

```python
fft_data = fft_data / np.max(np.abs(fft_data)) * 3
fft_data = np.log(fft_data + np.sqrt(1 + fft_data ** 2))
```

---

## Output

- The final plot is saved as `my_spectrum.png`.
- X-axis: $\omega_{\tau}$ (cm⁻¹)
- Y-axis: $\omega_t$ (cm⁻¹)
- Colormap: `'jet'`, normalized between -1 and 1

---

## Example Output

The plot shows a smooth, normalized, high-resolution 2D spectrum constructed from both rephasing and non-rephasing signals.

---

## For More Information

A Jupyter Notebook version is available in the same directory for interactive exploration and parameter tuning. Please refer to it for:

- Step-by-step explanation of each stage
- Visualization of intermediate results
- Parameter customization and testing
```
