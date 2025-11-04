# Software Pages

## Digital Research Alliance (formerly Compute Canada)

[Documentation][dra]

This is the national organization that manages our main computational resources.
Their documentation is comprehensive, including details about specific
[clusters][clusters], [available software][software], and [how to get started][get-started].
To use these resources, [you'll need to create a CCDB account](tutorials/ccdb.md).

Related:

- [Creating a CCDB account](tutorials/ccdb.md)
- [Transferring data to remote clusters](tutorials/data_transfer.md)
- [Setting up your cluster environment][cluster-setup]

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

## ssh

[Documentation][ssh]

SSH, or Secure Shell, is a network protocol that provides a secure, encrypted connection
between two computers. In the context of ComCat Lab, we use it to log in and execute
commands on DRA clusters, transfer files to and from clusters (with protocols like SFTP
and SCP), and "tunnel" other applications (e.g., VSCode). It is highly recommended that
you [create a password-protected SSH key](#create-ssh-key).

## Git

[Documentation][git]

Git is a version control system (VCS) that allows you to keep track of multiple versions
of files at once. With git it's easy to revert back to an old version of a file or merge
different versions of files. This is very helpful when you're trying to determine what
change resulted in a particular error. Chapters 1 to 3 provide a good
enough background to get started. See Section 5.1 for a description of
our git workflow, the ["Integration-Manager Workflow"][git-workflow].
You can find a guide to git best practices [here][git-best-practices].

## Miniconda

[Documentation][miniconda]

Miniconda is a software distribution that contains `conda`, a package and environment
manager for your command line interface, Python, and their dependencies. Miniconda
is a subset of Anaconda with the goal of reducing the memory footprint of the latter.
While we advocate the use of Miniconda/Anaconda on your local machine, the DRA
discourages their use in favour of using Lmod modules and Python virtual environments.

## Python

[Documentation][python]

Python is a relevant programming language in which many useful utilities for our
computational workflows are written. Python files can generally be identified by the
`.py` extension. When specific Python packages (software) is required, the use of
[virtual environments][venvs] is handy.

!!! Tip

    You can create a virtual environment by running the following command in your
    terminal:

    ```shell
    python -m venv name-of-virtual-environment
    ```

## Atomic Simulation Environment (ASE)

[Documentation][ase]

The Atomic Simulation Environment (ASE) is a set of tools and Python modules for
setting up, manipulating, running, visualizing and analyzing atomistic simulations.

## Pymatgen

[Documentation][pymatgen]

Pymatgen (Python Materials Genomics) is an open-source Python library for materials analysis.

## cclib

[Documentation][cclib]

cclib is an open source library, written in Python, for parsing and interpreting the results
of computational chemistry packages.

## Vienna Ab Initio Simulation Package (VASP)

[Documentation][vasp]

VASP is a proprietary, plane-wave-based code that we use for studying periodic systems.

## Quantum Espresso

[Documentation][espresso]

From their webpage:

> Quantum Espresso is an integrated suite of Open-Source computer codes for electronic-structure
> calculations and materials modeling at the nanoscale. It is based on density-functional theory,
> plane waves, and pseudopotentials.

## Gaussian

[Documentation][gaussian]

Gaussian is a proprietary, gaussian-orbital based code that we mainly use for studying
molecular systems. The Alliance also maintains a [documentation page][alliance-gaussian]
for running Gaussian on Alliance clusters.

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

## **Comp**utational **Chem**istry **Util**itie**s**

[Documentation][ccu]

CompChemUtils (Computational Chemistry Utilities; AKA `ccu`) is a Python package
containing several useful Python classes and routines for computational
chemistry such as adsorbate complex creation, defect creation, free energy
diagram generation, and more.

## Autojob

[Documentation][autojob]

`autojob` is a framework for executing high throughput calculations with the
flexibility to granularly resubmit and modifying calculations between
execution. Autojob also provides a codeless interface for generating input
directories for calculations.

## WIP

- Materials Cloud
- Materials Project & API
- VIM
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
[alliance-gaussian]: https://docs.alliancecan.ca/wiki/Gaussian
[espresso]: https://www.quantum-espresso.org
[orca]: https://www.faccts.de/docs/orca/6.0/manual/index.html
[markdown]: https://www.markdownguide.org
[mkdocs]: https://www.mkdocs.org
[lmod]: https://lmod.readthedocs.io/en/latest/
[lua]: https://www.lua.org/docs.html
[git]: https://git-scm.com/book/en/v2
[git-workflow]: https://www.git-scm.com/book/en/v2/ch00/wfdiag_b
[git-best-practices]: https://about.gitlab.com/topics/version-control/version-control-best-practices/
[ccu]: http://python-comp-chem-utils.rtfd.io/
[autojob]: http://python-autojob.readthedocs.io/
[cclib]: https://cclib.github.io
[cluster-setup]: https://github.com/ComCatLab/cluster-setup
[miniconda]: https://www.anaconda.com/docs/getting-started/miniconda/main
[venvs]: https://docs.python.org/3/tutorial/venv.html
[ssh]: https://docs.alliancecan.ca/wiki/SSH
