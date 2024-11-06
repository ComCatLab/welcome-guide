from ase import Atoms
import pytest

try:
    from gpaw import GPAW
except ImportError:
    GPAW = None


@pytest.mark.calculator
@pytest.mark.skipif(GPAW is None, reason="GPAW is not installed")
def test_gpaw_calculator() -> None:
    atoms = Atoms("N2", positions=[[0, 0, -1], [0, 0, 1]])
    atoms.write("myatoms.traj")
    calc = GPAW(mode="lcao", basis="dzp", txt="gpaw.txt", xc="LDA")

    atoms.calc = calc

    atoms.center(vacuum=3.0)
    print(atoms)

    e = atoms.get_potential_energy()
    print("The Energy is: ", e)

    f = atoms.get_forces()
    print("Forces \n", f)
