# :sparkles: Make your own Quantum ESPRESSO compilation :sparkles:

## Overview

This guide will walk through compiling Quantum Espresso with the following
options

- compiled with Intel OneAPI libraries
- signal trapping: This enables QE to recognize common signals sent by job schedulers
  (e.g., SIGTERM) and exit gracefully
- exit statuses
- [libxc][libxc] support: This expands the library of functionals that can be used
- [Environ][Environ]: This compiles Quantum Espresso with the Environ library,
  enabling the use of advanced solvent models (e.g., SCCS, SSCS)

## Tips

You may need to replace `def-samiras` with a valid project directory. To see
which project directories are available, run:

```shell
ls ~/projects
```

!!! note

    As of 2024/10/20, the installation instructions for Environ on the
    documentation are out-of-date for Environ 3+. The GitHub repository contains
    sparse instructions for Environ 3+ [here][Environ3-installation].

A good practice for managing the installed software is to store the raw source
code in a separate folder from where it is installed. For example, one can
keep the pre-compiled source code in a `software_support` folder and install
the software in a `software` folder. In this way, if one ever has to
compile the software, they can simply delete the associated subdirectory
in `software` without having to worry about the source code. Additionally,
for large source codebases, one can archive the directory in `software_support`
to conserve the file quota. This can be done for the Intel OneAPI,
QuantumEspresso, and libxc software in this tutorial.

## Step-by-Step

👩‍💻  Here's a (not really) quick step-by-step to compile
Quantum Espresso (QE) 7.3.1  👩‍💻

1. **Download the Intel OneAPI libraries** from [here][oneapi-libs].

2. **Copy the folder/archive of Intel OneAPI files** to where you
    will install QE. This may be in your home directory (e.g., `/home/USER`) or
    a subdirectory of your project folder (e.g.,
    `/home/USER/projects/def-samiras/USER/software_support`).

    ```shell
    scp -r oneapi_archive USER@cedar.computecanada.ca:/home/USER/
    ```

    where `USER` is your Digital Research Alliance (DRA) username.

    !!! note

        These libraries are from the [OneAPI suite][oneapi-suite], however,
        only the current year's libraries are available for download for free. The
        2024 libraries were tested and they did not work for QE compilation on
        Cedar. The 2022 libraries have been confirmed to have worked. Both the
        2022 and 2023 libraries are included in the shared folder from step 1.

3. **Install the Base Toolkit.** (This could take up to 30 minutes.)
    Navigate to where you copied the folder on Cedar

    ```shell
    $SHELL l_BaseKit_p_2022.1.1.119.sh
    ```

    This will install the Base Toolkit under the directory
    `/home/USER/intel/oneapi` (where `USER` is your DRA username).

    Optionally, you can change the installation location by opting to
    customize the installation (select "Accept & customize installation" prior
    to starting the installation process).

    ![Base Toolkit Installer Homepage](images/basetoolkit_welcome.png)

    You may want to do this if you would like the libraries to be accessible
    for others, in which case, change the default installation directory to a
    subdirectory of the project folder. For example,

    ![Base Toolkit Customization](images/basetoolkit_customization.png)

    !!! warning

        You may receive warnings about the operating system being
        "Unknown" or missing packages required for the Intel® VTune(TM) Profiler,
        but you can ignore these.

4. **Install the HPC Toolkit.** (This shouldn't take more than a few minutes.)

    ```shell
    $SHELL l_HPCKit_p_2022.1.1.97.sh
    ```

    Again, you may want to specify a custom installation location. In that case,
    it is reasonable to create a `oneapi` directory as a subdirectory of the
    custom location specified for the Base Toolkit and use this as the
    installation location.

    At this point, the directory that was chosen as the installation location
    should be populated with the Intel libraries, and there should be a file
    called `setvars.sh`, which when sourced (e.g., `source setvars.sh`), should
    give you a list of modules loaded in the environment. Check this step prior
    to QE compilation.

5. **Download Quantum Espresso 7.3.1**. You can get QE [here][qe-7.3.1].
    Alternatively, you can just go to their homepage, register your email,
    and go to Downloads > to find the latest version.

6. **Copy the downloaded `QE.tar` file to Cedar and extract it**.

7. **Forcibly purge your loaded modules** and reload the Gentoo Linux module.

    ```shell
    module --force purge
    module load gentoo
    ```

    !!! important

        The main problem with Cedar seems to be some
        library/ies dependency, so make sure you run module --force purge before
        moving on.

    !!! note

        The Gentoo Linux module provides access to the `git` CLI utility
        which is needed to configure the Environ module.

8. **Setup the Intel environment.**

    ```shell
    source /path/to/setvars.sh
    ```

9. **Locate the root directory for libxc.**
    On Cedar, this directory is found at
    `/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcc12/libxc/6.2.2`

    !!! note

        The latest version of libxc (7.0.0) is not installed on
        Cedar. However, installing libxc is quite straightforward. Instructions can
        be found [here][libxc-installation].

10. Ensure that the local language is set to the standard, i.e. ”C”.

    ```shell
    export LC_ALL=C
    ```

11. **Obtain a copy of Environ 3+** and copy it to the QE folder. You can
    either clone it from the git repository by running this command from
    inside of the QE folder.

    ```shell
    git clone https://github.com/environ-developers/Environ.git
    ```

    or you can [download an archive][Environ-releases], copy it to Cedar, and
    then extract it into the QE folder.

    !!! warning

        The `Environ` folder should be inside of the Quantum Espresso
        folder (e.g., `qe-X.Y.Z/Environ`).

12. **Configure the QE compilation.**

    For clarity, we define variables for the location of the Intel Base
    Toolkit and where you would like to install the QE binaries. If you
    didn't modify the default directory where the Intel libraries are
    installed, then you should define `intel_dir` as:

    ```shell
    intel_dir=/home/$USER/intel/oneapi
    ```

    If you used a custom directory, your variable definition may look like
    this:

    ```shell
    intel_dir=/home/$USER/projects/def-samiras/$USER/software/intel-2022.1.1
    ```

    Next, specify a path (outside of the current directory) where you would
    like to install the QE executables.

    ```shell
    espresso_dir=/home/$USER/projects/def-samiras/$USER/software/espresso-X.Y.Z
    ```

    !!! note

        These directories must be **absolute** paths (i.e., starting with `/`).

    Run the `configure` script from inside the `qe-X.Y.Z` folder (where `X.Y.Z`
    is the QE version number).

    ```shell
    ./configure LIBDIRS="$intel_dir/mkl $intel_dir/mpi $intel_dir/compiler" --enable-parallel --with-scalapack=intel FC=ifort F90=ifort mpif90=mpiifort CC=icc mpicc=mpiicc --enable-signals --enable-exit-status --prefix=$espresso_dir
    ```

    !!! note

        This will take a while, (between minutes and hour-ish), so be sure
        to have some time at this step. Keep an eye on what will come up at the
        screen, because this will inform whether or not QE found the libraries
        you told it to, or if it skipped any of them. It will also let you know
        if the parallel compilation was successfully identified (QE has different
        compilations for serial and parallel execution). Make sure it found all
        dependencies you need before moving on. If you are unable to see all the
        output from the command in your terminal, it may be useful to redirect
        the standard output and standard error from the `./configure` comand to a
        file:

        ```shell
        ./configure LIBDIRS="$intel_dir/mkl $intel_dir/mpi $intel_dir/compiler" --enable-parallel --with-scalapack=intel FC=ifort F90=ifort mpif90=mpiifort CC=icc mpicc=mpiicc --enable-signals --enable-exit-status --prefix=$espresso_dir >&espresso.log
        ```

    This command will redirect all the configuration information to
    `espresso.log`.

13. **Modify the `make.inc` file to configure libxc.**

    This includes:

    - adding `-D__LIBXC` to `DFLAGS`
    - adding `-I/path/to/libxc/include` to `IFLAGS`
    - setting `LD_LIBS=-L/path/to/libxc/lib -lxcf03 -lxc`

    Note that `/path/to/libxc` should be the path determined in Step 9.

14. **Configure Environ.**

    Change the current directory to inside the `Environ` folder and run:

    ```shell
    ./configure LIBDIRS="$intel_dir/mkl $intel_dir/mpi $intel_dir/compiler" --enable-parallel FC=ifort F90=ifort mpif90=mpiifort CC=icc mpicc=mpiicc --prefix=$espresso_dir
    ```

    !!! note

        It is very important to use the same libraries and parallelization
        flags for Environ as for QE.

15. **Compile Environ.**

    ```shell
    make -jN compile
    ```

    where `N` is the number of processors to use for parallel compilation.

16. **Build the QE executables.**

    ```shell
    make all
    ```

    :star2: The [QE manual][qe-manual] is really good and useful and easy to
    follow, in case you want to test different things on your own.

    :red_circle: You can also compile their default version, which is mostly
    foolproof and easy, and test it. It usually finds its own way and works
    well. Beware that there is a problem with library dependencies in the
    cluster, so if you compile the default version, it will have the same
    libraries as the cluster version and likely have the same issues!

    !!! note

        To speed up the installation, you can spawn an interactive
        job with multiple cores and run `make` in parallel. For example, if you
        want to spawn a job with 16 cores, do the following:

        ```shell
        salloc --mem=32GB --ntasks-per-node=16 --nodes=1 --account=def-samiras --time=00:10:00
        make -j16 all
        ```

    I've found that this scales linearly, with 32 cores taking about 2 minutes
    to complete.

17. **Install the executables.**

    ```shell
    make install
    ```

    This step copies the built executables to the directory specified by the
    `--prefix` option (the `espresso_dir` variable defined in Step 11).

## Modulefiles

Recall that the `setvars.sh` script **must** be sourced prior to executing QE
in order to setup the Intel environment for Quantum Espresso. This has a few
disadvantages, the most significant of which being that there is no
straightforward way to undo the changes. Additionally, there is the issue of
convenience and portability due to having to specify the *exact* location of
the `setvars.sh` script in order to source it. This section covers how to
create modulefiles for Quantum Espresso and the Intel OneAPI libraries.

### Intel OneAPI

Modulefiles solve the above problems and offer some additional benefits. In
short, modulefiles provide a convenient way to dynamically change the user's
environment. This may involve modifying environment variables, defining shell
functions, or loading other modulefiles. The Digital Research Alliance clusters
use [Lmod][lmod] to manage modulefiles. Lmod provides the `sh_to_modulefile`
utility. See [this tip](../resources/snippets.md#create-a-modulefile-from-a-script)
for how to convert `setvars.sh` into a modulefile.

Once you have created the modulefile, ensure that the file
is in your `MODULEPATH` (this is controlled by the command `module use`). A
reasonable strategy is to name the modulefile after the version of the
libraries (e.g., 2022.1.1) and to place this file in a subdirectory of a path
in `MODULEPATH`. The name of the subdirectory will be the name used to load
the module (i.e., set up the environment), so a reasonable name would be
something like `intel`. Additionally, it's useful to add the following
metadata to the modulefile:

```text title="intel/2022.1.1.lua"
help([[
Description
===========
Intel C, C++ & Fortran compilers (classic and oneAPI)


More information
================
 - Homepage: https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html
]])
whatis("Description: Intel C, C++ & Fortran compilers (classic and oneAPI)")
whatis("Homepage: https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html")
whatis("URL: https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html")
whatis("Version: 2022.1.1")
whatis("Keywords: HPC")
conflict("intel")
```

This mostly just descriptive information about the module. The last line
ensures that you don't load conflicting modules and is especially helpful
since we've had troubles with the default DRA libraries. An error will
be thrown if you try to load this module when another version is already
loaded.

If you placed the created file in a directory called `intel` in your `MODULEPATH`
and called the file `2022.1.1` (after the version of the Intel OneAPI
libraries), then you perform the required set up with the command:

```shell
module load intel/2022.1.1
```

### Quantum Espresso

You can also define a modulefile for your newly installed version of
Quantum Espresso. Essentially, the only necessary change that must be made
is to prepend the directory containing QE executables to the `PATH` environment
variable. A sample modulefile is shown below:

```text title="espresso/7.3.1.lua"
help([[
For detailed instructions, go to:
    https://www.quantum-espresso.org/documentation/package-specific-documentation/

]])

whatis("Version: 7.3.1")
whatis("Keywords: QuantumEspresso")

conflict("quantumespresso")
conflict("espresso")

depends_on("intel/2022.1.1")
prepend_path("PATH", "espresso_dir/bin")
```

Here, `espresso_dir` should be replaced with the path used to define
`espresso_dir` in Step 11. Also, note that the DRA module is named
`quantumespresso` whereas this module is named `espresso`. The `conflict`
commands reflect this fact. This is beneficial (but not necessary) since ASE's
Quantum Espresso calculator class is called `Espresso`, so both Python and
SLURM submission scripts can be templated with the same variable. Finally,
note that this modulefile "depends_on" the `intel/2022.1.1` modulefile. This
ensures that the `intel/2022.1.1` module is loaded prior to loading the
`espresso` module Optionally, one can explicitly load the intel module by
replacing that `depends_on` with `load`.

## Running a sample calculation

The following commands should be added to your SLURM submission script to run
calculations with QE.

```shell
module --force purge
source /home/USER/intel/oneapi/setvars.sh
export PATH=$PATH:/home/USER/qe-7.3.1/bin/
export PATH=$PATH:/home/USER/intel/oneapi/mpi/2021.5.0/bin
/home/YOUR$USERNAME/intel/oneapi/mpi/2021.5.0/bin/mpirun pw.x < espresso.in > espresso.out
```

where `USER` is your DRA username.

If you've defined modules as indicated [above](#modulefiles), add the
following commands instead.

```shell
module --force purge
module load gentoo intel/2022.1.1 espresso/7.3.1
mpirun pw.x
```

If running Quantum Espresso with ASE in a script called `run.py`, add the
following commands instead.

```shell
module --force purge
module load gentoo intel/2022.1.1 espresso/7.3.1

# Load your Python environment here, if necessary

python3 run.py
```

If you changed the `depends_on` command to `load`, then you don't need to load
the `intel` module.

## Benchmarking

Finally, benchmark and test its performance using their provided examples or
your job. What to look for:

- The performance should be equivalent to their examples for simple structures
  with similar parameters. For unit cells, it is useless to scale up too many
  cores, so test it with few at first (1 core vs 2 vs 4, for example). Check
  the last lines of the out file the total CPU time and resources. If the
  performance is worse with more cores, the compilation is not good.
- Check the performance of a decent sized job with your compilation against
  the one available in Cedar, the time should be around half but it will also
  depend on your own system (it was half for mine with a ~100 atoms
  supercell). If it is similar to Cedar, or worse, again, the compilation it
  not good.

[oneapi-libs]: https://1sfu-my.sharepoint.com/:u:/g/personal/lebarbos_sfu_ca/EYL5CdD_j99OmZirKP2MV4cBTJ-FPv0yE_fB28Jmnr_RVQ?e=xnV2Hd
[oneapi-suite]: https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.dhpux3
[qe-7.3.1]: https://www.quantum-espresso.org/rdm-download/488/v7-3-1/b44218b83782ba64fb5b47781e0fbcfb/qe-7.3.1-ReleasePack.tar.gz
[qe-manual]: https://www.quantum-espresso.org/Doc/user_guide_PDF/user_guide.pdf
[libxc]: https://libxc.gitlab.io/
[libxc-installation]: https://libxc.gitlab.io/installation/
[lmod]: https://lmod.readthedocs.io/en/latest/
[Environ]: https://environ.readthedocs.io/en/latest/index.html
[Environ3-installation]: https://github.com/environ-developers/Environ/blob/documentation/Doc/sphinx/source/install/install_3.rst
[Environ-releases]: https://github.com/environ-developers/Environ/releases
