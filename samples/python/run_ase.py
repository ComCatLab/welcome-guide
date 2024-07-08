from pathlib import Path

from ase import io
from ase.calculators.vasp import Vasp
from ase.optimize.bfgslinesearch import BFGSLineSearch
from numpy.linalg import norm

# Replace in.traj with the name of your structure file
atoms = io.read("in.traj")

# see https://www.vasp.at/wiki/index.php/Category:INCAR_tag
# for details on what each of these keywords mean
# if not found, check https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#module-ase.calculators.vasp
calc = Vasp(
    algo="Fast",
    ediff=1e-8,
    encut=450,
    gga="PE",
    gamma=False,
    ibrion=1,
    isif=2,
    ismear=0,
    ispin=1,
    ivdw=11,
    ldipol=False,
    lplane=True,
    lreal="Auto",
    ncore=4,
    nelm=60,
    # Note that nsw must equal zero in order to properly  use the ASE optimizer
    nsw=0,
    prec="Accurate",
    sigma=0.04,
    kpts=(4, 4, 1),
    dipol=None,
)

atoms.calc = calc

# Create the ASE optimizer and run the optimization
dyn = BFGSLineSearch(atoms, logfile="BFGS_output.log", trajectory="relax.traj", restart='relax.pckl')
dyn.run(fmax=0.01, steps=150)

e = atoms.get_potential_energy()
f = norm(max(atoms.get_forces(), key=norm))

# Write final structure to file
atoms.write('final.traj')

# Print final energy and max force to standard output
print(f'final energy {e}')
print(f'max force: {f}')

# Write final energy to file
with Path('final.e').open(mode='x', encoding='utf-8') as file:
    file.write(f'{e}\n')

