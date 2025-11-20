# Setting Up Your First Calculation

This tutorial explains how to set up your first calculations to run
on a Digital Research of Alliance cluster. In particular, this tutorial
will walk through how to run your first calculations on our home cluster,
Fir!

## System Requirements

- MacOS 26.0.1 or later
- A [valid CCDB account](../tutorials/ccdb.md)
- A VASP license

## Prerequisites

The prerequisites for this tutorial can be satisfied by completing the
[Basic Setup](../onboarding/basic_setup.md),
[Software Development Setup](../onboarding/development.md) and
[Cluster Setup](../onboarding/cluster_setup.md) tutorials. In particular,
it will be assumed that you have already:

- set up your local machine with Python
- set up VSCode
- created an SSH key
- created a repository in which to house your calculation files

Further, these instructions do not assume any familiarity with shell scripting
or Python; however, the following reference pages may be of use:

- [The Python Tutorial](https://docs.python.org/3/tutorial/index.html)
- [`scp`](https://docs.alliancecan.ca/wiki/Transferring_data#SCP)
- [`ssh`](https://docs.alliancecan.ca/wiki/SSH)
- [Linux introduction](https://docs.alliancecan.ca/wiki/Linux_introduction)

## Objectives

- to organize a calculation directory
- to obtain structure files from the [Materials Project][mat-pro]
- to submit your first calculation

## Step-by-Step Instructions

The following steps will make heavy use of VSCode, but the steps can
analogously be performed from the command line without VSCode. Many steps
will require you to type commands from the command line. This can be
done within the Terminal app, but they can also be executed from the terminal
subwindow of VSCode.

(screenshot of terminal subwindow)

If the VSCode terminal subwindow is ever not visible, you can always make it
visible by typing ``ctrl-` `` from within VSCode or selecting "Terminal" from
the "View" submenu.

(screenshot of view subwindow)

1. **Obtain a structure file for your calculation.**

    Depending on your desired application, this can be done any number of ways.
    As noted in [Computational Catalysis in a Nutshell](../nutshell.md), you
    can create the structure within ASE, obtain a `.cif` file from a paper, or
    download a structure from [Materials Project][mat-pro].

2. **Create a virtual environment in your project folder.**

    From the root of your project directory, type the following command:

    ```shell
    python -m venv .venv
    ```

    or by using the "Python: Create Environment" command from the VSCode
    command palette (accessible via the keyboard shortcut `option-x` or from
    the "View" submenu).

    ??? info "Explanation"

        This will create a virtual environment in your repository directory
        named `.venv`.

    ??? note

        Although it is possible to use a global environment (like that created
        in [Basic Setup](../onboarding/basic_setup.md)), different projects will likely
        require different packages, and creating separate environments ensures
        that modifications in the environment for one project do not affect the
        environment of another project.

    !!! tip

        You may receive a notification remarking on the creation of a virtual
        environment and prompting you to set the environment as the default
        for the workspace. Do so! This will prevent you from having to activate
        your virtual environment every time you open the terminal subwindow of
        VSCode, and it will ensure that when you run Python files from within
        your project workspace, you use the correct Python interpreter with
        all of your desired packages installed.

    Now, add the following line to your `.gitignore` file:

    ```text
    .venv
    ```

    ??? info "Explanation"

        This will prevent Git from tracking any files in this virtual
        environment.

3. **Activate the virtual environment and install some necessary packages.**

   ```shell
    source .venv/bin/activate
    pip install ase matplotlib ipython ruff pymatgen matplotlib mp-api python-autojob comp-chem-utils 
   ```

   For posterity, it can be useful to record the state of your virtual
   environment, so that you can replicate it on other machines. To do so,
   create a `requirements.txt` file like so:

   ```shell
   pip freeze > requirements.txt
   ```

   ??? info "Explanation"

        This command will record the names and versions of every package
        installed in the current Python environment into a `requirements.txt`
        named.

4. **Create folders to organize your files.**

   For starters, it is recommended to create folders for structures and
   calculations:

   ```shell
   mkdir structure calculations
   ```

[mat-pro]: https://next-gen.materialsproject.org
