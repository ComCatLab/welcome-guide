# Workflows

It takes several steps to go from a structure to a publishable result.
This page contextualizes common steps and highlights the infrastructure required to execute those steps.

1. Obtain a structure file.

    Structure files can be:

    - generated within various Python packages (e.g., ASE, Pymatgen)
    - obtained from the supplementary information of publications
    - retrieved from databases (e.g., Materials Project)

2. Perform structure manipulations.

    Depending on the project, this can any combination of the following:

    - modifying the structure (e.g., creating defects, substituting atoms, creating surfaces)
    - adding adsorbates
    - and more...

3. Creating input files.

    There are several ways to do this. Computational codes like VASP and Gaussian fully specify
    their input file formats and how to configure input files for various calculations. However,
    there are also several utilities that can help with this process.

    - ASE: provides Python interface to input file generation for a number of computational chemistry codes
    - MaterialsCloud: provides a graphical interface for generating input files for a number of computational chemistry codes

Typically, your workflow will consist of generating structure files (often a `.traj` file),
configuring calculation parameters (using ASE calculators in a python script), transferring those
files to a `/project` subdirectory on a computing cluster, configuring scheduler parameters (using
`#SBATCH` directives in a shell script), submitting your calculation using a scheduler, monitoring the
calculation using various `slurm` commands (e.g., `seff`), compiling the results, and analyzing the
results.
