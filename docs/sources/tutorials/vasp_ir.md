# Calculate IR Intensities with VASP (WIP)

This tutorial covers how to calculate IR spectra using VASP.

1. Relax the structure.

2. Conduct a phonon calculation.

   It is important to set the following `INCAR` tags:

   - `EDIFF`: should be no larger than 10<sup>-6</sup>
   - `IBRION`: set to 5 (do not use symmetry to reduce the number of
      directions) or 6  (use symmetry to reduce the number of directions)
   - `LEPSILON`: must be set to `.TRUE.` in order to calculate the Born
      effective field tensor
   - `NFREE`: set to either 2 or 4; determines how many displacements in
     each direction are used for the finite differences
   - `NSW`: must be set to exactly 1 in order for the finite differences
     step to be calculated
   - `POTIM`: should be set to a very small number; the default is 0.015
     but a valid number can be obtained by progressively increasing from 0.015
     and comparing with the results from a DFPT phonon calculation
