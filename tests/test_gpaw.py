from ase import Atoms
from ase.io import iread
from ase.io.trajectory import Trajectory
import matplotlib.pyplot as plt
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


@pytest.mark.calculator
@pytest.mark.skipif(GPAW is None, reason="GPAW is not installed")
def test_should_perform_gpaw_calculation() -> None:
    # Creates an N2 molecule with z positions of the N atoms
    atoms = Atoms("N2", positions=[[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]])

    # Loads the N2 into the ase gui and displays them
    # view(atoms)

    atoms.write("myatoms.traj")

    calc = GPAW(mode="lcao", basis="dzp", txt="gpaw.txt", xc="LDA")

    atoms.calc = calc

    atoms.center(vacuum=3.0)

    print(atoms)

    dstep = 0.05
    energies = []
    distances = []

    ## Prints initial position and sets the position of the second N to be 0.05 A from the first N and relcaulate energy and forces
    # print(atoms.positions)
    # atoms.positions[1, 2] = atoms.positions[0, 2] + 0.05
    # print(atoms.positions)

    traj = Trajectory("mytrajectory.traj", "w")

    # Move the second N in lengths of dstep from 0.5 A separation to 5 A separation and updates energies and forces
    for i in range(int((3.0 - 0.5) / dstep)):
        d = 0.5 + i * dstep
        atoms.positions[1, 2] = atoms.positions[0, 2] + d
        print(atoms.positions)
        e = atoms.get_potential_energy()
        print("The Energy is: ", e, " eV")

        f = atoms.get_forces()
        print("The Forces on each atom are: ", f, " ev/A")
        traj.write(atoms)

    # This extracts the energies and distances from the trajectory file and appends them to the appropriate arrays
    for atoms in iread("mytrajectory.traj"):
        energies.append(atoms.get_potential_energy())
        distances.append(atoms.positions[1, 2] - atoms.positions[0, 2])

    ax = plt.gca()
    ax.plot(distances, energies)
    ax.set_xlabel("Distance [A]")
    ax.set_ylabel("Total Energy [eV]")
    plt.savefig("N2_distances_vs_energy.png", bbox_inches="tight")
