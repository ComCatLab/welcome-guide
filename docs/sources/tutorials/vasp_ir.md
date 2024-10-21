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
  
   Additionally, the calculation must be a gamma-point only calculation since
   [as of 4/17/2024, VASP DFPT only works for $$q = \Gamma$$][MM-answer].
   Also note that, depending on how you run VASP, a separate executable may be
   selected based on the k-points chosen for the calculation. DFPT calculations
   are not supported with the `vasp_gam` executable.

## Resources

[VASP Wiki for IBRION INCAR Tag][VASP-IBRION]
[Stack Exchange Question with Instructions][MM-question]

[MM-question]: https://mattermodeling.stackexchange.com/questions/12494/phonon-calculations-using-dfpt
[VASP-IBRION]: https://www.vasp.at/wiki/index.php/IBRION#Computing_the_phonon_modes
[MM-answer]: https://mattermodeling.stackexchange.com/a/12496/5260
