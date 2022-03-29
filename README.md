# MolecularMechanics
Final project: Molecular Mechanics

This program was written as part of the Introduction to Scientific Computing and Programming course at the University of Amsterdam.
In this program, the minimization of energy from a given molecular sturcture is calculated. The molecular structure consists only of hydrogen and carbon atoms. The main program is combined with four other modules. In the AtomMod.f90, the input data is read from the text file and stored, and the variables are defined. In Calculationmod.f90, the bonds are assigned, and the angles are calculated. In Energymod.f90, the total energy and its contributions are calculated. In Minimizationmod.f90, the energy is minimized by making small random steps in the atom coordinates using the Metropolis algorithm. 
