from pathlib import Path

from ase import io
from ase.calculators.gaussian import Gaussian
from numpy.linalg import norm

# Replace in.traj with the name of your structure file
atoms = io.read("in.traj")

# see https://gaussian.com/keywords/
# for details on what each of these keywords mean
# if not found, check https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html
calc = Gaussian(
    command="g16 < Gaussian.com > Gaussian.log",
    label="Gaussian",
    chk="Gaussian",
    save=None,
    xc="pbe0",
    basis="gen",
    basis_set="-H -N -C -O\nDef2SVP\n****\n-Cu -Fe\nDef2TZVP\n****\n",
    charge=0,
    mult=3,
    scf="QC,MaxCycle=100,IntRep",
    opt="Tight,MaxCycles=250,CalcFC",
    mem="40GB",
)

atoms.calc = calc

e = atoms.get_potential_energy()
f = norm(max(atoms.get_forces(), key=norm))

# Write final structure to file
atoms.write("final.traj")

# Print final energy and max force to standard output
print(f"final energy {e}")
print(f"max force: {f}")

# Write final energy to file
with Path("final.e").open(mode="x", encoding="utf-8") as file:
    file.write(f"{e}\n")
