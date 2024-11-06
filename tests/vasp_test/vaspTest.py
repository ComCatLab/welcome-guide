from ase.build import molecule
from ase.calculators.vasp import Vasp
import pytest


@pytest.mark.calculator
def test_vasp_calculator() -> None:
    atoms = molecule("N2")
    atoms.center(vacuum=5)

    atoms.pbc = True  # VASP cannot handle non-periodic boundary conditions
    calc = Vasp(
        xc="pbe",  # Select exchange-correlation functional
        encut=400,  # Plane-wave cutoff
        kpts=(1, 1, 1),  # k-points
    )

    atoms.calc = calc
    e = atoms.get_potential_energy()  # This call will start the calculation
    print(f"Potential energy: {e:.2f} eV")
