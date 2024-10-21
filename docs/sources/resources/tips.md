<!-- markdownlint-disable MD024 -->
# Tips

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

### ASE is not copying the `vdw-kernel.bindat` to the job directory (WIP)

#### Context

Running any calculation with VDW corrections.

#### Environment

- **Cluster**: Cedar
- **Date**: July 9, 2024
- **Software Versions**: Python 3.11.9, ASE 3.23.0, VASP 5.4.4

#### The Error

You set the `ivdw` keyword argument when configuring a VASP calculator
in ASE, but you notice that the `vdw-kernel.bindat` file is not copied to the
calculation directory.

#### The Solution

As of ASE 3.23.0, the `vdw-kernel.bindat` is only copied if you also set the
keyword argument `luse_vdw=True` for the VASP calculator.

## Gaussian

### Frequency Calculations

Watch out for internal rotations! If performing frequency calculations in
Gaussian, you may get the following warning:

    Warning -- explicit consideration of  23 degrees of freedom as
            vibrations may cause significant error

See [here][gaussian-freq-error] for a resolution.

[gaussian-freq-error]: http://www.ccl.net/chemistry/resources/messages/2005/04/01.002-dir/

## Modulefiles

### Create a modulefile from a script

[Lmod](../software_pages.md#lmod) provides a handy utility
[`sh_to_modulefile`][sh_to_modulefile] for creating modulefiles from scripts.
`sh_to_modulefile` records the initial environment, runs the script, and
compares the final environment to determine the changes. It then converts these
changes to Lua commands and prints them to the terminal. Given a script,
`my_script.sh`, one can run the following command to generate an Lmod
modulefile written in Lua.

```shell
$LMOD_DIR/sh_to_modulefile my_script.sh > my_script.lua
```

The command `$LMOD_DIR/sh_to_modulefile` calls the utility by using the
Lmod-defined environment variable `LMOD_DIR`, which points to the directory
in which Lmod is installed. `my_script.sh` is the name of the script to be
converted into a modulefile. The output of the command is redirected (`>`)
to the file `my_script.lua`.

[sh_to_modulefile]: https://lmod.readthedocs.io/en/latest/260_sh_to_modulefile.html

## SLURM

### Attaching to Existing Jobs

Say, you have a window with an interactive job running (e.g., `salloc`) and
you would like to open another window with that same allocation. If the job
ID for the interactive job is JOBID, you can connect to the interactive job
with the command

```shell
srun --pty --overlap --jobid JOBID bash
```

## General

### Property Calculation

When calculating zero-point energies using computational codes, it is worth it
to examine whether there is an internal routine. In the case of VASP, one can
use the INCAR tag `IBRION=5` (or 6) to perform such a calculation using finite
differences. The benefit in this case over using the ASE reliant method
`ccu.thermo.vibration.run_vibration` is that all code-specific data is
retained. For example, one can also go back and collect IR frequency data
from the calculated dipoles of each image.

Additionally, one should also check whether IR frequency data is any more
expensive than ZPE data. For example, in VASP, dipole moments are calculated
at every optimization point during a phonon frequency calculation. This means
that if the phonon calculation is executed in VASP, IR frequencies are
obtained for free. This benefit highlights a major advantage to executing
frequency calculations in VASP as opposed to with an external routine like
ASE’s `ase.vibration.vibration.Vibration`. ASE does not archive all calculated
results for each image, but instead, only records calculated forces.
