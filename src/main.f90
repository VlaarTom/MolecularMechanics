!! Author:              Tom Vlaar
!! Studentnumber:       14101971
!! University:          University of Amsterdam
!! Last change date:    2022-03-22
!!
!! Molecular Mechanics project
!! Scientific Computing and Programming Course 
!!
!! In this program, the minimization of energy from a given molecular structure is calculated. The molecular 
!! structure consists only of hydrogen and carbon atoms. This main program is combined with four other modules: 
!! Atommod.90, Calculationmod.f90, Energymod.f90 and Minimizationmod.90.

program MolecularMechanics
    use AtomModule
    use CalculationModule
    use EnergyModule    
    use MinimizationModule
implicit none

type(Atom), ALLOCATABLE     :: Molecule(:)
type(Bonds), ALLOCATABLE    :: BondsArray(:)
type(Bonds)                 :: Bond
type(Parameter)             :: Variable
real*8, ALLOCATABLE         :: AnglesArray(:)
real*8, ALLOCATABLE         :: TorsionAnglesArray(:)
real*8                      :: Energy
integer                     :: TotalSamples, AcceptedSamples

call InputReader ('c4h10.xyz', Molecule)
call ParamaterInput (Variable)
call AssignBonds (Molecule, BondsArray, Bond, Variable)
call BendAngle (Molecule, BondsArray, AnglesArray)
call TorsionAngle (Molecule, BondsArray, TorsionAnglesArray)
Energy = TotalEnergy (Molecule, Variable, BondsArray, AnglesArray, TorsionAnglesArray)

print *, "The initial energy of the molecule:", Energy, "J/mol"
print *, "---------------------------------------------------------------------------------"

Call Minimization(Energy, Molecule, Variable, BondsArray, AnglesArray, TorsionAnglesArray, &
TotalSamples, AcceptedSamples)

print *, "The minimized energy of the molecule:", Energy,"J/mol"
print *, "Number of new positions tested: ", TotalSamples 
print *, "Number of new positions accepted: ", AcceptedSamples

end program MolecularMechanics