# Python Scripts

This page describes the collection of Python sample scripts contained
in the [Python sample scripts folder][scripts].

## Run a VASP relaxation using the internal quasi-Newton optimization algorithm (RMM-DIIS) in VASP

``` py title="samples/python/run.py"
--8<-- "./samples/python/run.py"
```

!!! important "Reminder"
    Replace `"in.traj"` with the name of your structure file.

## Run a VASP relaxation using the BFGSLineSearch optimization routine in ASE

``` py title="samples/python/run_ase.py"
--8<-- "./samples/python/run_ase.py"
```

!!! important "Reminder"
    Replace `"in.traj"` with the name of your structure file.

[scripts]: https://github.com/ComCatLab/welcome-guide/tree/main/samples/python
