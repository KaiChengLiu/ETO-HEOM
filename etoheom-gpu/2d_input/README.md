Files in this directory is the input files of 2DES simulation, each
file represent a conherent time and population time.

follow the instructions below:
1. modify .key-tmpl file according to your system.
2. modify gen_2d.sh, especially the "of_name"
3. run gen_2d.sh input files.


SIZE: Define the system size including ground state and first excited states.

HEOM: Define the HEOM configuration by (SITE_NUMBER) (TRUNCATION_LEVEL)

HAMILTONIAN: Define the system Hamiltonian including ground state and first excited states.

DISORDER: Define the static disorder matrix of Hamiltonian. The first two number in first row are sampling times and random seed

BATHTYPE: Define the bath type, currently only support "artificial".

BATH: Define bath information, it will be automatically written by {spectral_density_type}_ETOM.py in ../bath_model.

DIPOLE: Define dipole information, number in first row is the number of transition diploe vector.
        The following rows first represent the dipole direction and its amplitude in XYZ component.

POLARIZATION: Define the angles account for the four pulses polarization. You don't have to change it in normal case.

PULSE: Define the pulse profile, number in first row is the number of pulse, now only support 3 pulses.
       The following 3 rows represent the amplitude, central time, width and frequency of pulse. The key
       words TAU1, TAU2 and TAU3 will be replace to correspond values by gen_2d.sh.

TIME: Define the simulation and sample time. first two number are start and end time of simulation where T_END  
      be replace to correspond values by gwn_2d.sh. Third number is time step and fourth number in sampling interval.
