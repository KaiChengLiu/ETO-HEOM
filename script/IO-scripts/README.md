# IO-script

This repository generates input files for 2D electronic spectroscopy (2DES) simulations using pre-defined templates and configurations.

---

## Workflow

Follow the instructions below to generate and manage input files:

---

### 1. Modify `key.key-tmpl`

This is the template input file containing all system-specific parameters. Update the following sections to reflect your system configuration:

- **SIZE**: Define the system size including the ground and first excited states.

- **HEOM**: Define the HEOM configuration:  
  Format — `(SITE_NUMBER) (TRUNCATION_LEVEL)`

- **HAMILTONIAN**: Define the system Hamiltonian (ground + excited states).

- **DISORDER**: Static disorder matrix of the Hamiltonian.  
  First row: number of samples and random seed.

- **BATH**: Bath information, automatically generated by `{spectral_density_type}_ETOM.py` in `../bath_model`.

- **DIPOLE**:  
  - First row: number of transition dipole vectors  
  - Following rows: dipole direction and amplitude in XYZ components

- **POLARIZATION**: Define polarization angles for the four pulses.  
  *(You typically don't need to modify this.)*

- **PULSE**:  
  - First row: number of pulses (currently only supports 3)  
  - Next 3 rows: amplitude, central time, width, and frequency  
  - Use placeholders `TAU1`, `TAU2`, `TAU3` that will be replaced automatically

- **TIME**:  
  - Format — start time, end time (`T_END` placeholder), time step, and sampling interval  
  - `T_END` will be replaced during script execution

---

### 2. Modify `gen_2d_data.sh`

This shell script controls how the input files are generated. Key parameters you can adjust:

```bash
t0=210                # Initial central time for the first pulse
propagate_time=600    # Duration to propagate after the last pulse

tau_step=10           # Step size for tau2 (in fs)
tau_bound=600         # Upper and lower bounds for tau2 offset (in fs)

T=(0)                 # Array of population times (in fs). You can add more values.
input_file="key.key-tmpl"  # Path to your template file
```

The script will:

- Loop over values of `tau2` from `-tau_bound` to `+tau_bound`
- Automatically compute `tau1`, `tau2`, `tau3`, and `T_END` based on `t0`, `tau_step`, and `T`
- Replace placeholders (`TAU1`, `TAU2`, `TAU3`, `T_END`) in the template
- Output files as `./2d-input/key_{tau2}_{T}.key`

You can modify this script to include more values of population time `T` or adjust the pulse timing logic.

---

### 3. Run `gen_2d_data.sh`

Execute the script to generate the input files:

```bash
./gen_2d_data.sh
```

Generated files will be saved in the `./2d-input` directory.

---

### 4. Clean Input Files: `clean_2d_data.sh`

Use this script to remove all generated input files and clean the `./2d-input` directory.

```bash
./clean_2d_data.sh
```

The script deletes all files in the `./2d-input` folder:
