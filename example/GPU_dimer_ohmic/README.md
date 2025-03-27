# Workflow Overview for 2DES Simulation Pipeline

This README provides the step-by-step workflow for generating, running, and visualizing 2D electronic spectroscopy (2DES) simulations using HEOM and GPU/CPU implementations.

---

## Step-by-Step Instructions

### 1. Generate ETOM Data

Run the appropriate ETOM module based on your spectral density type:

```bash
python {BATHTYPE}_ETOM.py
```

Replace `{BATHTYPE}` with one of the supported types:
- `debye_lorentz`
- `ohmic`
- `superohmic`

The script will generate ETOM bath parameters and automatically insert them into your input template.

📄 For more information, refer to:  
`../etom/README.md`

---

### 2. Edit the Input Template

Update `key.key-tmpl` to match your system configuration, such as system size, Hamiltonian, dipole directions, etc.

📄 For detailed guidance, see:  
`../script/IO-scripts/README.md`

---

### 3. Generate Input Files

Execute the input generation script:

```bash
./gen_2d_data.sh
```

This will create a set of simulation input files:

```
./2d-input/key_{TAU}_{T}.key
```

These files include time-domain configurations for each TAU and T combination.

---

### 4. Submit Simulation Jobs

- **For CPU version**, submit the batch job using:

  ```bash
  ./submit_jobs.sh
  ```

  This will generate and submit PBS job scripts, producing output files like:

  ```
  ./2d-output/out_{TAU}_{T}.out
  ```

- **For GPU version**, directly run:

  ```bash
  ./Run_GPU_2DES.sh
  ```

  This script executes the GPU-accelerated 2DES simulations and generates the same output format.

---

### 5. Generate 2D Spectrum

Visualize and process the simulation results using:

```bash
python gen_2d_data.py
```

This script performs 2D FFT and interpolation to produce the final spectrum plot.

📄 For configuration and visualization details, refer to:  
`../gen-2d-data/README.md`

---

### 6. Clean Up Generated Files (Optional)

If needed, you can delete all generated `.key` and `.out` files using:

```bash
./clean_2d_data.sh
```

This will clear:

- `./2d-input/key_*.key`
- `./2d-output/out_*.out`

---

This modular workflow allows for high flexibility in simulation setup, automated batch processing, and high-resolution spectral visualization.
