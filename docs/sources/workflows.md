# Workflows

Outline the general workflow of going from idea to result with links to relevant software pages.

Typically, your workflow will consist of generating structure files (often a .traj file), configuring calculation parameters (using ASE calculators in a python script), transferring those files to a /project subdirectory on a computing cluster, configuring scheduler parameters (using #SBATCH directives in a shell script), submitting your calculation using a scheduler, monitoring the calculation using various slurm commands (e.g., seff), compiling the results, and analyzing the results.
