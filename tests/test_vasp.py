import os
from pathlib import Path

from ase.build import molecule
from ase.calculators.vasp import Vasp
import pytest

VASP_INSTALLED = (
    os.environ.get("ASE_VASP_COMMAND")
    or os.environ.get("VASP_COMMAND")
    or os.environ.get("VASP_SCRIPT")
)


@pytest.mark.calculator
@pytest.mark.skipif(
    not VASP_INSTALLED,
    reason="The environment is not configured for VASP. "
    "See the ASE documentation for details.",
)
def test_should_perform_vasp_calculation(tmp_path: Path) -> None:
    atoms = molecule("N2")
    atoms.center(vacuum=5)

    atoms.pbc = True  # VASP cannot handle non-periodic boundary conditions
    calc = Vasp(
        xc="pbe",  # Select exchange-correlation functional
        encut=400,  # Plane-wave cutoff
        kpts=(1, 1, 1),  # k-points
        directory=tmp_path,
    )

    atoms.calc = calc
    e = atoms.get_potential_energy()  # This call will start the calculation
    print(f"Potential energy: {e:.2f} eV")
    assert e is not None
