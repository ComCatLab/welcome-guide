# Standard imports for ase
from pathlib import Path

from ase import Atoms
from ase.build import fcc111
import numpy as np
import pytest


@pytest.fixture(name="water_layer")
def fixture_water_layer() -> Atoms:
    p = np.array(
        [
            [0.27802511, -0.07732213, 13.46649107],
            [0.91833251, -1.02565868, 13.41456626],
            [0.91865997, 0.87076761, 13.41228287],
            [1.85572027, 2.37336781, 13.56440907],
            [3.13987926, 2.3633134, 13.4327577],
            [1.77566079, 2.37150862, 14.66528237],
            [4.52240322, 2.35264513, 13.37435864],
            [5.16892729, 1.40357034, 13.42661052],
            [5.15567324, 3.30068395, 13.4305779],
            [6.10183518, -0.0738656, 13.27945071],
            [7.3856151, -0.07438536, 13.40814585],
            [6.01881192, -0.08627583, 12.1789428],
        ]
    )
    c = np.array(
        [[8.490373, 0.0, 0.0], [0.0, 4.901919, 0.0], [0.0, 0.0, 26.93236]]
    )
    water_layer = Atoms("4(OH2)", positions=p, cell=c, pbc=[1, 1, 0])
    water_layer.cell = [water_layer.cell[1, 1], water_layer.cell[0, 0], 0.0]
    water_layer.rotate(90, "z", center=(0, 0, 0))
    water_layer.wrap()
    return water_layer


@pytest.fixture(name="lattice_constant")
def fixture_lattice_constant() -> float:
    lattice_constant = 3.55
    return lattice_constant


@pytest.fixture(name="nickel_slab")
def fixture_nickel_slab(lattice_constant: float) -> Atoms:
    nickel_slab = fcc111(
        "Ni", size=[2, 4, 3], a=lattice_constant, orthogonal=True
    )
    return nickel_slab


def test_should_create_nickel_slab(nickel_slab: Atoms, tmp_path: Path) -> None:
    nickel_slab.center(vacuum=5.0, axis=2)

    # Add silver atom
    h = 1.9
    relative = (1 / 6, 1 / 6, 0.5)
    absolute = (*np.dot(relative, nickel_slab.cell), 0, 0, h)
    nickel_slab.append("Ag")
    nickel_slab.positions[-1] = absolute

    # Write file
    xyz_file = tmp_path.joinpath("slab.xyz")
    nickel_slab.write(xyz_file)

    assert xyz_file.exists()


def test_should_create_slab_with_water_layer(
    nickel_slab: Atoms,
    water_layer: Atoms,
    tmp_path: Path,
) -> None:
    # Add water layer to slab
    water_layer.set_cell(nickel_slab.cell, scale_atoms=True)

    # Shift water layer
    zmin = water_layer.positions[:, 2].min()
    zmax = nickel_slab.positions[:, 2].max()
    water_layer.positions += (0, 0, zmax - zmin + 1.5)

    # Add water layer to slab
    interface = nickel_slab + water_layer
    interface.center(vacuum=6, axis=2)

    # Write file
    traj_file = tmp_path.joinpath("NiH2O.traj")
    interface.write(traj_file)

    assert traj_file.exists()
