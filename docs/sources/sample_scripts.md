# Sample Scripts

This page describes the collection of sample scripts contained
in the [sample scripts folder][scripts].

[scripts]: https://github.com/ComCatLab/welcome-guide/tree/main/sample_scripts

## Python Scripts

`sample_scripts/python_scripts/run.py`

Run a VASP relaxation using the internal quasi-Newton optimization algorithm (RMM-DIIS) in VASP.

`sample_scripts/python_scripts/run_ase.py`

Run a VASP relaxation using the BFGSLineSearch optimization routine in ASE.

## bash profile files (`.bashrc`)

`sample_scripts/bash_profiles/bashrc`

A standard bash profile for use on remote clusters defining
utility functions and environment variables for software.
**Don't forget to replace 'username' with your username.**

## Slurm submission files

`sample_scripts/slurm/run_vasp.sh`

Submit a VASP calculation. The following variables should
be updated:

- Slurm options (e.g., time, memory, notification email)
- The name of the Python script used to execute the calculation (if any)

`sample_scripts/slurm/run_espresso.sh`

Submit a Quantum Espresso calculation. The following variables should
be updated:

- Slurm options (e.g., time, memory, notification email)
- The name of the Python script used to execute the calculation (if any)

`sample_scripts/slurm/run_gaussian.sh`

Submit a Gaussian calculation. The following variables should
be updated:

- Slurm options (e.g., time, memory, notification email)
- The name of the Python script used to execute the calculation (if any)

`sample_scripts/slurm/run_orca.sh`

Submit an ORCA calculation. The following variables should
be updated:

- Slurm options (e.g., time, memory, notification email)
- The name of the Python script used to execute the calculation (if any)

## Software Input Files

- VASP
- Quantum Espresso
- Gaussian
- ORCA
