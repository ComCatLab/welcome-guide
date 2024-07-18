from ase.build import molecule

atoms = molecule('N2')
atoms.center(vacuum=5)

atoms.pbc=True # VASP cannot handle non-periodic boundary conditions
from ase.calculators.vasp import Vasp

calc = Vasp(xc = 'pbe', # Select exchange-correlation functional
            encut = 400, # Plane-wave cutoff
            kpts = (1, 1, 1) #k-points
            )

atoms.calc = calc
e = atoms.get_potential_energy() # This call will start the calculation
print('Potential energy: {:.2f} eV'.format(e))
