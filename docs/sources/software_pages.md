# Software Pages

## Digital Research Alliance (formerly Compute Canada)

[Documentation][dra]

This is the national organization that manages our main computational resources.
Their documentation is comprehensive, including details about specific
[clusters][clusters], [available software][software], and [how to get started][get-started].
To use these resources, [you'll need to create a CCCDB account](tutorials/cccdb.md).

Related:

- [Creating a CCCDB account](tutorials/cccdb.md)
- [Transferring data to remote clusters](tutorials/data_transfer.md)
- [Setting up your cluster environment](tutorials/cluster_setup.md)

## Slurm Workload Manager

[Documentation][slurm]

Slurm is the job scheduler that runs on our remote clusters. You can find a useful
cheatsheet [here][slurm-cheatsheet]. `sbatch`, `squeue`, and `sacct` are commonly
used for submitting and monitoring jobs.

## bash

[Documentation][bash]

bash is a shell program and shell-scripting language. Notably, it is the default 
on most compute clusters. Navigating and interacting with the file system on remote
clusters is done using bash. For example, the command-line functions `ssh` and `scp`
are used to connect to and transfer files between the Alliance clusters.

## Git

[Documentation][git]

Git is a version control system (VCS) that allows you to keep track of multiple versions
of files at once. With git it's easy to revert back to an old version of a file or merge
different versions of files. This is very helpful when you're trying to determine what
change resulted in a particular error. Chapters 1 to 3 provide a good
enough background to get started. See Section 5.1 for a description of
our git workflow, the ["Integration-Manager Workflow"][git-workflow].

## Python

[Documentation][python]

Python is a relevant programming language in which many useful utilities for our
computational workflows are written. Python files can be identified by the `.py` extension.

## Atomic Simulation Environment (ASE)

[Documentation][ase]

The Atomic Simulation Environment (ASE) is a set of tools and Python modules for
setting up, manipulating, running, visualizing and analyzing atomistic simulations.

## Pymatgen

[Documentation][pymatgen]

Pymatgen (Python Materials Genomics) is an open-source Python library for materials analysis.

## Vienna Ab Initio Simulation Package (VASP)

[Documentation][vasp]

VASP is a proprietary, plane-wave-based code that we use for studying periodic systems.

## Quantum Espresso

[Documentation][espresso]

Quantum Espresso is an open-source plane-wave-based code that we use for studying
periodic systems.

## Gaussian

[Documentation][gaussian]

Gaussian is a proprietary, gaussian-orbital based code that we mainly use for studying
molecular systems.

## ORCA

[Documentation][orca]

ORCA is a free-for-academic gaussian-orbital based code that we mainly use for studying
molecular systems.

## Markdown

[Documentation][markdown]

Markdown is a markup language commonly used for writing documentation. A common extension
for Markdown files is `.md`. The documentation for this Welcome Guide is written in
Markdown. Markdown files can easily be converted into the pages of a static website using
a tool like [mkdocs][mkdocs].

## Lmod

[Documentation][lmod]

Lmod is a software for managing software environments. The necessary commands
and variables required for a software to be used are specified in 'module files'
that are written in the [Lua programming language][lua]. Modules are managed using
commands like `module load`, `module purge`, and `module unload`.

## WIP

- Materials Cloud
- Materials Project
- Globus
- FileZilla
- Cyberduck
- bader
- Multiwfn
- Lobster
- chargemol
- Fireworks
- jobflow
- Custodian
- atomate2
- catmap
- matplotlib
- ccu (e.g., FancyPlots)
- autojob

[vasp]: https://www.vasp.at/wiki/index.php/Main_page
[ase]: https://wiki.fysik.dtu.dk/ase/index.html
[slurm]: https://slurm.schedmd.com/documentation.html
[dra]: https://docs.alliancecan.ca/wiki/Technical_documentation
[clusters]: https://docs.alliancecan.ca/wiki/National_systems#Compute_clusters
[software]: https://docs.alliancecan.ca/wiki/Available_software
[get-started]: https://docs.alliancecan.ca/wiki/Getting_started
[pymatgen]: https://pymatgen.org
[python]: https://www.python.org
[bash]: https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html
[slurm-cheatsheet]: https://slurm.schedmd.com/pdfs/summary.pdf
[gaussian]: https://gaussian.com/man/
[espresso]: https://www.quantum-espresso.org
[orca]: https://www.orcasoftware.de/tutorials_orca/index.html
[markdown]: https://www.markdownguide.org
[mkdocs]: https://www.mkdocs.org
[lmod]: https://lmod.readthedocs.io/en/latest/
[lua]: https://www.lua.org/docs.html
[git]: https://git-scm.com/book/en/v2
[git-workflow]: https://www.git-scm.com/book/en/v2/ch00/wfdiag_b
