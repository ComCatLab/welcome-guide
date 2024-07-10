import logging

import ase.io
from ase.calculators.vasp import Vasp
from ccu.relaxation import run_relaxation

logging.basicConfig(level=logging.DEBUG)

ldau_luj = {
    "H": {"L": -1, "U": 0, "J": 0},
    "N": {"L": -1, "U": 0, "J": 0},
    "O": {"L": -1, "U": 0, "J": 0},
    "S": {"L": -1, "U": 0, "J": 0},
    "Ag": {"L": 2, "U": 3.87, "J": 0},
    "Co": {"L": 2, "U": 5.2, "J": 0},
    "Mn": {"L": 2, "U": 5.3, "J": 0},
    "Ni": {"L": 2, "U": 6.1, "J": 0},
}

atoms = ase.io.read("Co-BHT_on_Co_vertical_1.traj")

calc = Vasp(
    algo='Normal',
    ediff=1e-8,
    ediffg=-0.01,
    encut=550,
    gga='PE',
    gamma=False,
    ibrion=5,
    isif=2,
    ismear=0,
    ispin=2,
    isym=0,
    ivdw=11,
    kpts=(4, 4, 1),
    kpar=2,
    ldautype=3,
    ldau=True,
    ldipol=False,
    lmaxmix=4,
    lorbit=11,
    lplane=True,
    lreal='Auto',
    npar=2,
    nelm=250,
    nfree=2,
    nsw=1,
    potim=0.015,
    prec='Accurate',
    sigma=0.04,
    smass=-3,
    ldau_luj=ldau_luj,
)

atoms.calc = calc
run_relaxation(atoms, run_bader=False, run_chargemol=False)
