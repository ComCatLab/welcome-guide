from pathlib import Path

from ase import io
from ase.calculators.espresso import Espresso
from ase.calculators.espresso import EspressoProfile
from numpy.linalg import norm

# Replace in.traj with the name of your structure file
atoms = io.read("in.traj")

# see https://www.quantum-espresso.org/Doc/INPUT_PW.html
# for details on what each of these keywords mean
# if not found, check https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html
input_data = {
    "calculation": "md",
    "verbosity": "low",
    "restart_mode": "from_scratch",
    "tstress": False,
    "tprnfor": False,
    # This should be set to less than the job time limit
    "max_seconds": 75000,
    "occupations": "smearing",
    "smearing": "gaussian",
    "degauss": 0.1,
    "vdw_corr": "grimme-d3",
    "dftd3_version": 4,
}

# Ensure that you define an entry for each atom type in your structure
pseudopotentials = {
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Co": "Co_pbe_v1.2.uspp.F.UPF",
    "S": "s_pbe_v1.4.uspp.F.UPF",
}

profile = EspressoProfile(
    # Replace the path in pseudo_dir with the path to your QE pseudopotential
    # folder
    pseudo_dir="/home/ugognw/projects/def-samiras/ugognw/software_support/espresso/SSSP_1.3.0_PBE_precision",
    command="mpiexec pw.x -nband 2 -ntg 2",
)

calc = Espresso(
    input_data=input_data,
    pseudopotentials=pseudopotentials,
    profile=profile,
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
