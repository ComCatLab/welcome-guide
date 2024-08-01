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
