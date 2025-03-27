# ETOM (Effective Thermal Oscillator Model Generator)

This directory contains Python modules for simulating the time correlation function (TCF) using the effective thermal oscillators (ETOs) under different spectral densities. These modules extract ETOM parameters that can be directly used for 2DES simulations.

---

## Included Python Scripts

The folder includes the following three modules:

- **`debye_lorentz_ETOM.py`**  
  Simulates TCF under the Debye–Lorentz spectral density using the ETOM method.

- **`ohmic_ETOM.py`**  
  Simulates TCF under the Ohmic spectral density using the ETOM method.

- **`superohmic_ETOM.py`**  
  Simulates TCF under the Super-Ohmic spectral density using the ETOM method.

Each script computes ETOM parameters and outputs them in a format that can be automatically inserted into the `BATH` section of the `key.key-tmpl` input file for HEOM-based 2DES simulations.

---

## Output

The generated ETOM parameters are automatically written to the `BATH` section of `key.key-tmpl`, so users do not need to manually edit them.

---

## How to Use

For detailed instructions and example usage, refer to the corresponding `.ipynb` (Jupyter Notebook) files provided in this folder. These notebooks demonstrate how to:

- Define spectral density parameters
- Run the ETOM modules
- Visualize the generated thermal correlation functions (TCFs)
- Automatically update `key.key-tmpl` with output parameters

---
