from collections.abc import Generator
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator

try:
    from gpaw import GPAW
except ImportError:
    GPAW = None


@pytest.fixture(name="n2")
def fixture_n2(gpaw_calc: Calculator) -> Atoms:
    n2 = Atoms("N2", positions=[[0, 0, -1], [0, 0, 1]])
    n2.center(vacuum=3.0)
    n2.calc = gpaw_calc
    return n2


@pytest.fixture(name="gpaw_calc")
def fixture_gpaw_calc(tmp_path: Path) -> Generator[Calculator, None, None]:
    output_file = tmp_path.joinpath("gpaw.txt")
    with output_file.open(mode="w", encoding="utf-8") as file:
        yield GPAW(mode="lcao", basis="dzp", txt=file, xc="LDA")


@pytest.mark.calculator
@pytest.mark.skipif(GPAW is None, reason="GPAW is not installed.")
def test_should_calculate_n2_energy_and_forces(n2: Atoms) -> None:
    e = n2.get_potential_energy()
    f = n2.get_forces()

    print("The Energy is: ", e)
    print("Forces \n", f)

    assert None not in (e, f)


@pytest.mark.calculator
@pytest.mark.skipif(GPAW is None, reason="GPAW is not installed.")
def test_should_calculate_energy_vs_n2_bond_length(n2: Atoms, tmp_path: Path) -> None:
    # Distances are in picometers
    start = 50
    stop = 300
    dstep = 5
    energies = []
    distances = []

    # Move the second N in lengths of dstep from 0.5 Å separation to 3 Å
    # separation and updates energies and forces
    for d in range(start, stop, dstep):
        n2.positions[1, 2] = n2.positions[0, 2] + d
        e = n2.get_potential_energy()
        f = n2.get_forces()

        energies.append(e)
        distances.append(d)

        print(n2.positions)
        print("The Energy is: ", e, " eV")
        print("The Forces on each atom are: ", f, " ev/A")

    ax = plt.gca()
    ax.plot(distances, energies)
    ax.set_xlabel("Distance [A]")
    ax.set_ylabel("Total Energy [eV]")
    figure = tmp_path.joinpath("N2_distances_vs_energy.png")
    plt.savefig(figure, bbox_inches="tight")
    assert figure.exists()
