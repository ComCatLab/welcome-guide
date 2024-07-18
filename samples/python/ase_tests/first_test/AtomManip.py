# Standard imports for ase
from math import sqrt
import numpy as np
from ase import Atoms
from ase import visualize
from ase.visualize import view
from ase.io import read
from ase.build import fcc111

a = 3.55
atoms = Atoms('Ni4',
              cell = [sqrt(2) * a, sqrt(2) * a, 1.0, 90, 90, 120],
              pbc = (1, 1, 0),
              scaled_positions = [(0, 0, 0),
                                (0.5, 0, 0),
                                (0, 0.5, 0),
                                (0.5, 0.5, 0)])
atoms.center(vacuum = 5.0, axis = 2)

atoms.write('slab.xyz')

h = 1.9
relative = (1/6, 1/6, 0.5)
absolute = np.dot(relative, atoms.cell) + (0, 0, h)
atoms.append('Ag')
atoms.positions[-1] = absolute

#view(atoms)
W = read('WL.traj')

slab = fcc111('Ni', size=[2, 4, 3], a=3.55, orthogonal=True)

W.cell = [W.cell[1, 1], W.cell[0, 0], 0.0]
W.rotate(90, 'z', center = (0, 0, 0))
W.wrap()
W.set_cell(slab.cell, scale_atoms = True)

zmin = W.positions[:, 2].min()
zmax = slab.positions[:, 2].max()

W.positions += (0, 0, zmax - zmin + 1.5)

interface = slab + W
interface.center(vacuum = 6, axis = 2)
interface.write('NiH2O.traj')
