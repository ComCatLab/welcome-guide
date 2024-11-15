from pathlib import Path

from ase import io
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.vibrations.vibrations import Vibrations
from numpy.linalg import norm

# Replace in.traj with the name of your structure file
atoms = io.read("in.traj")

# see https://www.vasp.at/wiki/index.php/Category:INCAR_tag
# for details on what each of these keywords mean
# if not found, check https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#module-ase.calculators.vasp
calc = Vasp(
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

print("Running vibration calculation")
forces = atoms.get_forces(apply_constraint=True)
atoms.write("final.traj")

if forces is None:
    msg = "Unable to perform vibrational calculation due to missing forces"
    raise RuntimeError(msg)

print(">>> BEGIN print full force")

for force, atom in zip(forces, atoms, strict=True):
    print("%s %s: %s, %s, %s", atom.index, atom.symbol, *force)

print("norm of force: %s", norm(forces))
print("max force: %s", norm(max(forces, key=norm)))
print("<<< END print full force")

fixed = []
if atoms.constraints:
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            fixed = constraint.index
            break

print("Fixed indices: %s", fixed)

to_vibrate = [atom.index for atom in atoms if atom.index not in fixed]
print("Indices to vibrate: %s", to_vibrate)

vib = Vibrations(
    atoms=atoms, nfree=4, name="vib", delta=0.015, indices=to_vibrate
)
vib.run()

with Path("vib.txt").open(mode="w", encoding="utf-8") as file:
    vib.summary(log=file)

for mode in range(len(to_vibrate) * 3):
    vib.write_mode(mode)

zpe = vib.get_zero_point_energy()

print("Successfully ran vibration calculation. ZPE: %s", zpe)
