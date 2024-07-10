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

``` title="INCAR"
INCAR created by Atomic Simulation Environment
 ENCUT = 550.000000
 POTIM = 0.015000
 SIGMA = 0.040000
 SMASS = -3.000000
 EDIFF = 1.00e-10
 EDIFFG = -1.00e-02
 ALGO = Normal
 GGA = PE
 PREC = Accurate
 IBRION = 5
 ISIF = 2
 ISMEAR = 0
 ISPIN = 2
 KPAR = 4
 LDAUTYPE = 3
 LMAXMIX = 4
 LORBIT = 11
 NELM = 250
 NFREE = 2
 NSW = 1
 IVDW = 11
 NCORE = 4
 LDAU = .TRUE.
 LDIPOL = .FALSE.
 LPLANE = .TRUE.
 LREAL = Auto
 LDAUL = 2 -1 -1 -1 -1
 LDAUU = 5.200 0.000 0.000 0.000 0.000
 LDAUJ = 0.000 0.000 0.000 0.000 0.000
```

#### The Error

```title="vasp.out"
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
