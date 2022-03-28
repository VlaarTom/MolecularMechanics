program MolecularMechanics

use AtomModule
use CalculationModule
use EnergyModule

implicit none

type(Atom), ALLOCATABLE :: Molecule(:)
type(Bonds), ALLOCATABLE :: BondsArray(:)
type(Bonds) :: Bond
type(Parameter) :: Variable
real*8, ALLOCATABLE :: AnglesArray(:)
real*8, ALLOCATABLE :: TorsionAnglesArray(:)



call InputReader ('c4h10.xyz', Molecule)
call ParamaterInput (Variable)
call AssignBonds (Molecule, BondsArray, Bond, Variable)
call BendAngle (molecule, BondsArray, AnglesArray)
call TorsionAngle (Molecule, BondsArray, TorsionAnglesArray)

print *, TotalEnergy (Molecule, Variable, BondsArray, AnglesArray, TorsionAnglesArray)
end program MolecularMechanics