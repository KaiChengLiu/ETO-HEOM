# ETO-HEOM: Efficient Simulation of 2D Electronic Spectroscopy by Hierarchical Equation of Motion

**ETO**: Effective Thermal Oscillator  
**HEOM**: Hierarchical Equation of Motion

---

## Description

This code simulates two-dimensional electronic spectroscopy (2DES) of multilevel excitonic systems using the Hierarchical Equation of Motion (HEOM) with the ETO (Effective Thermal Oscillator) model. Compared to traditional HEOM approaches, the ETO model provides significantly improved computational efficiency, especially for arbitrary spectral densities and low temperature scene.

The simulation framework uses HEOM to compute numerically exact quantum dynamics of open quantum systems. The density matrix is propagated using a fourth-order Runge-Kutta (RK4) integrator.

The code also supports:
- Multiple laser pulses
- Finite pulse width effects

These features enable the simulation of nonlinear optical signals such as three-pulse photon-echo peak shift measurements and other four-wave mixing spectroscopies.

---

## CPU and GPU Versions

- **CPU Version**: Suitable for small systems or testing scenarios.  
- **GPU Version**: Recommended for large-scale systems, offering significant acceleration on supported NVIDIA GPUs.

Both versions share the same input/output format and simulation framework.

---

## Installation and Job Setup

For instructions on installation, configuration, and how to set up and process jobs, please refer to the `INSTALL` file.

---

## Directory Overview

