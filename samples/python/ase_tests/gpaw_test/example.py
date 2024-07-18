#!/usr/bin/python3

from ase import Atoms
atoms = Atoms('N2', positions=[[0, 0, -1], [0, 0, 1]])

# from ase.visualize import view
# view(atoms)

from ase.io import write
write('myatoms.traj', atoms)

from gpaw import GPAW

calc = GPAW(mode='lcao', basis='dzp', txt='gpaw.txt', xc='LDA')

atoms.calc = calc

atoms.center(vacuum=3.0)
print(atoms)

e = atoms.get_potential_energy()
print('The Energy is: ', e)

f = atoms.get_forces()
print('Forces \n', f)
