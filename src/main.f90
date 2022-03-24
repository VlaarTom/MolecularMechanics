program MolecularMechanics

use AtomModule
use CalculationModule
use EnergyModule

implicit none

type(Atom), ALLOCATABLE :: Molecule(:)
type(Bonds), ALLOCATABLE :: BondsArray(:)
type(Bonds) :: Bond
type(Parameter) :: Variable


call InputReader ('c4h10.xyz', Molecule)
call ParamaterInput(Variable)
call AssignBonds (Molecule, BondsArray, Bond, Variable)
call BendAngle (molecule, BondsArray)!, Angles)

print *, StretchEnergy (BondsArray, Variable)
print *, Molecule
print *, BondsArray
end program MolecularMechanics