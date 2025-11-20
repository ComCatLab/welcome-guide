import logging

from ase.calculators.vasp.vasp import Vasp
import ase.io
from ccu.workflows.calculation import run_calculation

logging.basicConfig(level=logging.DEBUG)

# Replace in.traj with the name of your structure file
atoms = ase.io.read("in.traj")

# see https://www.vasp.at/wiki/index.php/Category:INCAR_tag
# for details on what each of these keywords mean
# if not found, check https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#module-ase.calculators.vasp
atoms.calc = Vasp(
    algo="Normal",
    ediff=1e-8,
    ediffg=-1e-2,
    encut=450,
    gga="PE",
    ivdw=12,
    kpts=(1, 1, 1),
    nelm=50,
    nsw=50,
    prec="Accurate",
)
run_calculation(atoms)
