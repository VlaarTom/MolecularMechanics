module EnergyModule

    use AtomModule
    use CalculationModule

implicit none
PRIVATE
PUBLIC StretchEnergy!, BendingEnergy

    
contains
real*8 function StretchEnergy (BondsArray, Variable)
!The stretch energy function is a function of all pairs of atoms in the molecule.

    type(Bonds), INTENT(IN), allocatable    :: BondsArray(:)
    type(Parameter), INTENT(IN)             :: Variable
    integer                                 :: i

    StretchEnergy = 0.0
    do i = 1, size(BondsArray)
        if (BondsArray(i)%Elements == "CH") then
            StretchEnergy = StretchEnergy + Variable%ForceConstantCH*((BondsArray(i)%Length - Variable%CHBond)**2)
        else
            StretchEnergy = StretchEnergy + Variable%ForceConstantCC*((BondsArray(i)%Length - Variable%CCBond)**2) 
        endif
    enddo
end function

!real*8 function BendingEnergy (BondsArray)
end module EnergyModule