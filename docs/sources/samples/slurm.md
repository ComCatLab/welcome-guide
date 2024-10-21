# Slurm Submission Files (WIP)

## Submit a DFT calculation

=== "VASP"

    ``` py title="samples/slurm/vasp.sh"
    --8<-- "./samples/slurm/vasp.sh"
    ```

=== "Espresso"

    ``` py title="samples/slurm/espresso.sh"
    --8<-- "./samples/slurm/espresso.sh"
    ```

    !!! important "Reminder"
        This script assumes that you are using a self-compiled version of
        Quantum Espresso and have created a corresponding module named
        `espresso`. See [this tutorial](../tutorials/espresso_compilation.md)
        for how to compile Quantum Espresso and create the necessary
        modulefile.

=== "Gaussian"

    ``` py title="samples/slurm/gaussian.sh"
    --8<-- "./samples/slurm/gaussian.sh"
    ```

=== "ORCA"

    ``` py title="samples/slurm/orca.sh"
    --8<-- "./samples/slurm/orca.sh"
    ```

!!! Reminder

    Don't forget to replace `JOB_NAME`, `SFU_ID`, and `PYTHON_SCRIPT` with
    appropriate values in addition to setting your desired SLURM parameters.
    Also, if you don't define the path to a Python virtual environment in your
    `.bashrc` file, then you should replace `$COMP_CHEM_ENV` with the path to
    the `activate` script (usually, `path-to-environment/bin/activate`).
