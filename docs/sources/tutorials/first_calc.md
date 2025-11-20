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

    and commit your changes.

    ```shell
    git add .gitignore
    git commit -m 'Ignore virtual environment folder'
    ```

    !!! tip

        You can perform this step in the VSCode GUI as well!

        (TODO: add screenshots)

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
            named.commit your changes.

    You should also commit these changes:

    ```shell
    git add requirements.txt
    git commit -m 'Record Python dependencies'
    ```

4. **Create folders to organize your files.**

    For starters, it is recommended to create folders for structures and
    calculations. This can be done from the command line like so:

    ```shell
    mkdir structures calculations
    ```

    or using the VSCode GUI:

    (TODO: add screenshot of VSCode create folder button)

    Copy the structure that you obtained from step 1 into the `structures/`
    directory and commit the change:

    ```shell
    git add structures
    git commit -m 'Add first structure file'
    ```

5. **Create your first calculation directory.**

    Our first calculation will be a relaxation. Assuming that the structure
    that you obtained in step 1 is named `structure`, create a subfolder in
    the `calcuations/` directory called `structure_relaxation`.

6. **Populate your first calculation directory.**

    To run calculations on DRA clusters, you will generally need three files:

    - *a structure file*: This describes the structure on which the calculation
      will be performed.
    - *a SLURM submission script*: This script specifies both the resource
      requirements for the job as well as the commands to execute during the
      job. Generally, the commands to execute include loading necessary
      [modules][dra-modules], calling a Python script to execute the
      calculation, and performing any necessary clean up once the job
      completes.
    - *a Python script*: This script actually runs the calculation. This
      generally leverages an interface (e.g., [ASE][ase]) with an underlying
      computational code (e.g., [VASP][vasp]).

    Copy the structure file from your `structures` directory into the
    `structure_relaxation` folder.

    (TODO: screenshot of directory structure with new file)

    In addition to a structure file, you will also need a SLURM submission
    script and a Python script to perform the calculation. Samples for each
    of these scripts can be obtained from [here](../samples/slurm.md) and
    [here](../samples/python.md), respectively.

    In the SLURM submission script, you will need to modify the lines which
    specify your resource requirements to match your job:

    ```shell
    ...
    #SBATCH --account=def-samiras
    #SBATCH --job-name=JOB_NAME
    #SBATCH --mem-per-cpu=1000MB
    #SBATCH --nodes=2
    #SBATCH --ntasks-per-node=24
    #SBATCH --time=23:00:00
    #SBATCH --mail-user=SFU_ID@sfu.ca
    #SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT,TIME_LIMIT_90
    ...
    ```

    In the Python script, you will need to change both the name of the file
    read by `ase.io.read`:

    ```python
    ...
    atoms = ase.io.read("in.traj")
    ...
    ```

    as well as the parameters used to configure the ASE calculator:

    ```python
    ...
        atoms.calc = Vasp(
        gga="PE",
        gamma=False,
        ibrion=1,
        isif=2,
        ismear=0,
        kpts=(1, 1, 1),
        nelm=60,
        nsw=100,
        prec="Accurate",
    )
    ...
    ```

    !!! tip

        You may find the following links helpful:

        - [SLURM submission parameters](https://slurm.schedmd.com/sbatch.html)
        - [ASE calculator parameters](https://ase-lib.org/ase/calculators/calculators.html)

    !!! tip

        Confirm that your modifications make sense with a group member.

    Once you have edited the SLURM and Python scripts, save and commit your
    changes:

    ```shell
    git add calculations
    git commit -m 'Create calc folder'
    ```

7. **Push your changes.**

    (TODO: add screenshot)

    ```shell
    git push
    ```

    Now, your calculation files are accessbile remotely.

8. **Connect to your favorite cluster.**

    Now, `ssh` into the cluster on which you would like to perform the
    calculations:

    ```shell
    ssh -Y <dra-username>@fir.alliancecan.ca
    ```

    or, if you created aliases:

    ```shell
    fir
    ```

    and authenticate with 2FA.

9. **Start the SSH agent.**

    ```shell
    eval $(ssh-agent) && ssh-add path/to/SSH/private/key
    ```

10. **Clone your calculation repository.**

    From the command-line, type:

    ```shell
    cd ~/ && git clone git@github.com:GH_USERNAME/GH_REPO.git
    ```

    !!! tip

        Dont' forget to replace `GH_USERNAME` and `GH_REPO` with your GitHub
        username and the name of your project repository, respectively.

    ??? into "Explanation"

        This will download the current state of your Git repository from GitHub
        into a folder named after the repository.

    You should now be able to see all of the files from your repository:

    ```shell
    cd GH_REPO
    ls
    ```

11. **Navigate to your calculation directory and submit the calculation.**

    `cd` into the calculation directory that you created in step 5. For
    example, assuming that you named this directory `structure_relaxation`:

    ```shell
    cd calculations/structure_relaxation
    ```

    Using the `ls` command, you should see the three files that you created
    in step 6.

    You can now submit the calculation using `sbatch`. Assuming that you named
    your SLURM submission script `run.sh`:

    ```shell
    sbatch vasp.sh
    ```

    The SLURM job ID should be printed to the terminal with a message
    indicating that your job was submitted. If you receive an error, there
    is likely an issue with your submission script.

12. **Monitor the calculation.**

    The SLURM scheduler provides several utilities for monitoring submitted
    jobs and viewing the statistics of completed jobs. In particular, the
    [`sacct`][sacct] and `seff` commands can be useful. The simplest
    way to invoke these commands is like so:

    ```shell
    sacct -j JOB_ID
    seff JOB_ID
    ```

    where `JOB_ID` is the SLURM job ID of the job.

[mat-pro]: https://next-gen.materialsproject.org
[dra-modules]: https://docs.alliancecan.ca/wiki/Utiliser_des_modules/en
[ase]: https://ase-lib.org
[vasp]: https://ase-lib.org/ase/calculators/vasp.html
[sacct]: https://slurm.schedmd.com/sacct.html
