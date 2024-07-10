# Troubleshooting

This page contains potential solutions to whatever may currently be troubling you. Your
mileage may vary.

## VASP

### VASP internal routines  have requested a change of the k-point set

#### Context

Running phonon calculations for nitrate reduction intermediates on M-BHT slabs.

#### Environment

- **Cluster**: Cedar
- **Date**: July 9, 2024
- **Software Versions**: Python 3.11.9, ASE 3.23.0, VASP 5.4.4

#### Relevant Files

Input files for the calculation.

=== "INCAR"

    ```text
    --8<-- "./docs/sources/resources/files/kpt_set_change_error/INCAR"
    ```

=== "run.py"

    ```py
    --8<-- "./docs/sources/resources/files/kpt_set_change_error/run.py"
    ```

=== "vasp.sh"

    ```py
    --8<-- "./docs/sources/resources/files/kpt_set_change_error/vasp.sh"
    ```

#### The Error

<!-- markdownlint-disable-next-line MD046 -->
```text title="vasp.out"
 -----------------------------------------------------------------------------
|                                                                             |
|     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###     |
|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
|     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       #      |
|     E        R   R    R   R    O     O  R   R                               |
|     E        R    R   R    R   O     O  R    R      ###     ###     ###     |
|     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |
|                                                                             |
|      VASP internal routines  have requested a change of the k-point set.    |
|      Unfortunately this is only possible if NPAR=number of nodes.           |
|      Please remove the tag NPAR from the INCAR file and restart the         |
|      calculations.                                                          |
|                                                                             |
|      ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ..., BYE!!! <----       |
|                                                                             |
 -----------------------------------------------------------------------------
```

This is an error due to symmetry breaking (see [here][ibrion-error]).

[ibrion-error]: https://mattermodeling.stackexchange.com/a/9013

#### The Solution

Set `ISYM=0`. in your `INCAR` file.
