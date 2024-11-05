# Python Scripts

This page describes the collection of Python sample scripts contained
in the [Python sample scripts folder][scripts].

## Run a VASP relaxation using the internal quasi-Newton optimization algorithm (RMM-DIIS) in VASP

``` py title="samples/python/run_vasp.py"
--8<-- "./samples/python/run_vasp.py"
```

!!! important "Reminder"
    Replace `"in.traj"` with the name of your structure file.

## Run a mixed-basis Gaussian relaxation using the internal Berny optimization algorithm in Gaussian

``` py title="samples/python/run_gaussian.py"
--8<-- "./samples/python/run_gaussian.py"
```

!!! important "Reminder"
    Replace `"in.traj"` with the name of your structure file, and read
    [this][gaussian-alliance] article about selecting which executable to pass
    to the `command` keyword argument to the `Gaussian` constructor. Note
    that the syntax is *slightly* different
    (`g16 < Gaussian.com vs. G16 Gaussian.com`).

## Run a VASP relaxation using the BFGSLineSearch optimization routine in ASE

``` py title="samples/python/run_ase.py"
--8<-- "./samples/python/run_ase.py"
```

!!! important "Reminder"
    Replace `"in.traj"` with the name of your structure file.

## Run a VASP relaxation using [`ccu`][ccu]

``` py title="samples/python/run_ccu.py"
--8<-- "./samples/python/run_ccu.py"
```

`ccu` is a set of tools for computational chemistry workflows. In particular,
[`run_relaxation`][ccu-run-relaxation] is a wrapper function around
`ase.atoms.Atoms.get_potential_energy()` that handles the logging and archiving
of a calculation's final results.

!!! important "Reminder"
    Replace `"in.traj"` with the name of your structure file.

## Run an AIMD calculation with Quantum Espresso

``` py title="samples/python/run_espresso.py"
--8<-- "./samples/python/run_espresso.py"
```

!!! important "Reminder"
    Replace `"in.traj"` with the name of your structure file, and don't forget to
    change the string passed to the `pseudo_dir` keyword argument of the
    `EspressoProfile` constructor to the absolute path to where your Quantum
    Espresso pseudopotentials are located.

[scripts]: https://github.com/ComCatLab/welcome-guide/tree/main/samples/python
[ccu]: https://python-comp-chem-utils.readthedocs.io/en/latest/
[ccu-run-relaxation]: https://python-comp-chem-utils.readthedocs.io/en/latest/reference/ccu.html#ccu.relaxation.run_relaxation
[gaussian-alliance]: https://docs.alliancecan.ca/wiki/Gaussian#G16_(G09,_G03)
