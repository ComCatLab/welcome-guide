import logging

from ase.calculators.vasp import Vasp
import ase.io
from ccu.relaxation import run_relaxation

logging.basicConfig(level=logging.DEBUG)

atoms = ase.io.read("in.traj")

calc = Vasp(
    gga="PE",
    gamma=False,
    ibrion=1,
    isif=2,
    ismear=0,
    kpts=(1, 1, 1),
    nelm=60,
    nsw=100,
    prec="Accurate",
)

atoms.calc = calc
run_relaxation(atoms, run_bader=False, run_chargemol=False)
