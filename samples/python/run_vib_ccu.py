import logging

from ase.calculators.vasp import Vasp
import ase.io
from ccu.workflows.vibration import run_vibration

logging.basicConfig(level=logging.DEBUG)

# Replace in.traj with the name of your structure file
atoms = ase.io.read("in.traj")

# see https://www.vasp.at/wiki/index.php/Category:INCAR_tag
# for details on what each of these keywords mean
# if not found, check https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#module-ase.calculators.vasp
atoms.calc = Vasp(
    algo="Fast",
    encut=450,
    gga="PE",
    gamma=False,
    ibrion=1,
    isif=2,
    ismear=0,
    ncore=4,
    nelm=60,
    nsw=100,
    prec="Accurate",
    sigma=0.1,
    kpts=(4, 4, 1),
)

# This will only vibrate O atoms; redefine indices to control which
# atoms will be vibrated. If indices is not specified, only those non-
# fixed atoms will be vibrated
indices = [a.index for a in atoms if str(a.symbol) == "O"]
run_vibration(atoms, nfree=4, indices=indices)
