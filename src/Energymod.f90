module EnergyModule
    use AtomModule
    use CalculationModule
implicit none
PRIVATE
PUBLIC StretchEnergy, BendingEnergy, TorsionEnergy, NonBondedEnergy, TotalEnergy

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

real*8 function BendingEnergy (AnglesArray, Variable)
!The BendingEnergy function calculates the bending energy of all triples of atoms. 
    real*8, INTENT(IN), ALLOCATABLE :: AnglesArray(:)
    type(Parameter), INTENT(IN)     :: Variable
    integer                         :: i
    
    BendingEnergy = 0.0
    do i = 1, size(AnglesArray)
        if (abs(acos(AnglesArray(i)) - Variable%EquiAngle) < 1) then
            BendingEnergy = BendingEnergy + Variable%ForceConstantAngle*((acos(AnglesArray(i)) - Variable%EquiAngle)**2)
        endif
    enddo
end function

real*8 function TorsionEnergy (TorsionAnglesArray, Variable)
!The TorsionEnergy function calculates the energy of all torsional angles.
    real*8, INTENT(IN), ALLOCATABLE :: TorsionAnglesArray(:)
    type(Parameter), INTENT(IN)     :: Variable
    integer                         :: i

    TorsionEnergy = 0.0
    do i = 1, size(TorsionAnglesArray)
        TorsionEnergy = TorsionEnergy + 0.5 * Variable%V1 * (1 + cos(Variable%n * TorsionAnglesArray(i) - Variable%gamma))
    enddo
end function

real*8 function NonBondedEnergy (Molecule, Variable)
!The NonBondedEnergy function calculates the non bonded energy of the atoms which can be divided into the Van-der Waals 
!interactions and the electrostatic interactions. The if statements determine the use of the correct variables with the 
!corresponding atom interactions. The j > i statements make sure that the Van der Waals interactions are counted double 
!and the electrostatic interactions once. 
    type(atom), INTENT(IN), ALLOCATABLE :: Molecule(:)
    type(Parameter), INTENT(IN)         :: Variable
    integer                             :: i, j

    NonBondedEnergy = 0.0

    do i = 1, size(Molecule)
        do j = i, size(Molecule)
            if (i /= j) then
                if (Molecule(i)%Element == "H" .and. Molecule(j)%Element == "H") then
                    if (j > i) then
                        NonBondedEnergy = NonBondedEnergy + &
                        Variable%Avogadro*(Variable%HCharge**2)*Variable%CoulombConstant / &
                        (BondLength(Molecule,i,j)/(10.0**10))
                    endif
                    NonBondedEnergy = NonBondedEnergy + &
                    4*Variable%HWellDepth*(((Variable%HSigma/BondLength(Molecule,i,j))**12) &
                    -((Variable%HSigma/BondLength(Molecule,i,j)))) 
                elseif (Molecule(i)%Element == "C" .and. Molecule(j)%Element == "C") then
                    if (j > i) then
                        NonBondedEnergy = NonBondedEnergy + &
                        Variable%Avogadro*(Variable%HCharge**2)*Variable%CoulombConstant / &
                        (BondLength(Molecule,i,j)/(10.0**10))
                    endif
                    NonBondedEnergy = NonBondedEnergy + &
                    4*Variable%HWellDepth*(((Variable%HSigma/BondLength(Molecule,i,j))**12) &
                    -((Variable%HSigma/BondLength(Molecule,i,j)))) 
                else
                    if (j > i) then
                        NonBondedEnergy = NonBondedEnergy + &
                        Variable%Avogadro*(Variable%HCharge**2)*Variable%CoulombConstant / &
                        (BondLength(Molecule,i,j)/(10.0**10))
                    endif
                    NonBondedEnergy = NonBondedEnergy + ((4*((Variable%CWellDepth/2) + (Variable%HWellDepth/2))) * &
                    ((((Variable%CSigma/2)+(Variable%HSigma/2))/BondLength(Molecule,i,j)**12) - &
                    ((((Variable%CSigma/2)+(Variable%HSigma/2))/BondLength(Molecule,i,j)**6))))
                endif
            endif
        enddo
    enddo         
end function

real*8 function TotalEnergy (Molecule, Variable, BondsArray, AnglesArray, TorsionAnglesArray)
!The TotalEnergy function calculates the total energy of the molecule by the sum of its contributions.
    type(Atom), INTENT(IN), ALLOCATABLE     :: Molecule(:)
    type(Parameter), INTENT(IN)             :: Variable
    type(Bonds), INTENT(IN), ALLOCATABLE    :: BondsArray(:)
    real*8, INTENT(IN), ALLOCATABLE         :: AnglesArray(:)
    real*8, INTENT(IN), ALLOCATABLE         :: TorsionAnglesArray(:)

    TotalEnergy = StretchEnergy(BondsArray, Variable) + BendingEnergy(AnglesArray, Variable) &
    + TorsionEnergy(TorsionAnglesArray, Variable) + NonBondedEnergy(Molecule, Variable)
end function
end module EnergyModule