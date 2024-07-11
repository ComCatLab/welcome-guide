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
