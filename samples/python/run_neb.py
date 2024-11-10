"""This script needs testing and refactoring!"""

from pathlib import Path
from typing import TypeVar

import matplotlib as mpl

mpl.use("Agg")
import subprocess
import sys

import ase
from ase import Atoms
import ase.calculators.vasp as vasp_calculator
from ase.io import read
import ase.io.vasp
from ase.neb import NEB

_S = TypeVar("_S")
_T = TypeVar("_T")


def sort_species(to_sort: Atoms, symbols: list[str]) -> None:
    """Sort the atoms in an Atoms object.

    This method modifies `to_sort` in place.

    Args:
        to_sort: An Atoms object whose atoms are to be sorted.
        symbols: A list of strings indicating the order in which the atoms
            in `to_sort` will be sorted.
    """
    old_symbols = to_sort.get_chemical_symbols()
    old_positions = to_sort.get_positions()
    symbols_and_positions = list(zip(old_symbols, old_positions, strict=False))
    new_symbols = []
    new_positions = []

    for symbol_to_match in symbols:
        for symbol, position in symbols_and_positions:
            if symbol == symbol_to_match:
                new_symbols.append(symbol)
                new_positions.append(position)

    to_sort.set_chemical_symbols(new_symbols)
    to_sort.set_positions(new_positions)


def dict_to_list(d: dict[_S, _T]) -> list[tuple[_S, _T]]:
    """Convert a dictionary to a list of 2-tuples.

    Args:
        d: A dictionary to convert.

    Returns:
        A list of 2-tuples (`key`, `value`) where each `value` corresponds
        to `key` in `d`.
    """
    _list: list[tuple[_S, _T]] = []
    for name, value in d.items():
        _list.append((name, value))
    return _list


def species(to_count: Atoms) -> list[tuple[str, int]]:
    """Count how many of each symbol exists in an Atoms object.

    Args:
        to_count: An Atoms object.

    Returns:
        A list of 2-tuples (`symbol`, `count`) where `count` is the number of
        symbol in `to_count`.
    """
    _dict = {}
    for i in to_count.get_chemical_symbols():
        if i not in _dict:
            _dict[i] = 1
        else:
            _dict[i] += 1
    _list = dict_to_list(_dict)
    return _list


try:
    submitdir = sys.argv[1]
except IndexError:
    submitdir = ""
if submitdir != "":
    submitdir += "/"

# SCRIPT STARTS HERE

# total number of images (fixed initial and final ones also count)
nimages = 12
# options are 'idpp' or ''
idpp = "idpp"
# if False, it will only generate the interpolated POSCAR files
run = True


POSCAR = True
for i in range(nimages):
    if Path(f"{i:0>2}/POSCAR").exists():
        pass
    else:
        POSCAR = False

if not POSCAR:
    initial = read("initial.traj")
    final = read("final.traj")
    images = [initial]
    for _ in range(nimages - 2):
        image = initial.copy()
        images.append(image)

    images.append(final)

    neb = NEB(images, climb=False, k=0.1)
    if idpp:
        neb.interpolate(idpp)
    else:
        neb.interpolate()

    for index, image in enumerate(neb.images):
        Path(f"{index:0>2}").mkdir()
        _species = species(image)
        # sort_species is used to order positions according to _species list,
        # which will be defined as symbol_count later
        sort_species(to_sort=image, symbols=[x for x, y in _species])
        # symbol count is used to write POSCARS in compact notation e.g. H C O
        # rather than H C O C H (3 species only)
        # this is relevant because otherwise vasp will take the former example
        # as 5 species instead of 3
        # and will crash because potcar only has 3 available species.
        ase.io.vasp.write_vasp(
            f"{index:0>2}/POSCAR", image, symbol_count=species(image)
        )

extra_string = ""
if POSCAR:
    extra_string += "POSCARs read from ##/POSCAR files."
else:
    extra_string += (
        f"POSCARs not read. Generating POSCARs from ASE NEB({idpp}) "
        "interpolation method."
    )

subprocess.call(  # noqa: S602
    f"echo 'computing {nimages!s} images. {extra_string}'",
    shell=True,
    stdout=None,
    cwd=".",
)


if run:
    atoms = read("00/POSCAR")
    calc = vasp_calculator.Vasp(
        encut=400,
        xc="PBE",
        gga="PE",
        ncore=8,
        isif=2,
        images=nimages - 2,  # start NEB
        spring=-5.0,
        ichain=0,
        lclimb=True,  # end NEB
        ivdw=11,
        kpts=(1, 1, 1),
        gamma=True,  # Gamma-centered (defaults to Monkhorst-Pack)
        ismear=0,
        sigma=0.1,
        nelm=250,
        algo="fast",
        ibrion=1,  # -1 for no relaxation with vasp, 1 otherwise
        ediffg=-0.01,  # forces
        ediff=1e-8,  # energy conv.
        prec="Accurate",
        nsw=500,  # don't use the VASP internal relaxation, only use ASE
        lreal="Auto",
        ispin=1,
    )
    atoms.set_calculator(calc)
    e = atoms.get_potential_energy()
    print("final energy", e)
    with Path("final.e").open(mode="w", encoding="utf-8") as f:
        f.write(str(e) + "\n")
