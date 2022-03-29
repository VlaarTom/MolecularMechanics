module MinimizationModule
    use AtomModule
    use CalculationModule
    use EnergyModule
implicit none
PRIVATE
public Minimization, NewCoordinates

contains
subroutine Minimization(Energy, Molecule, Variable, BondsArray, AnglesArray, TorsionAnglesArray, &
    TotalSamples, AcceptedSamples)
!Subroutine Minimization determines wheter the new energy value is lower than the previous value.
!If the value is lower, the new energy value and coordinates are accepted, if the value is higher than
!the previous value, the value is rejected and new coordinates are calculated. The algorithm is terminated
!by calculating a rejection rate.
    real*8, INTENT(INOUT)                   :: Energy
    type(Atom), INTENT(INOUT), ALLOCATABLE  :: Molecule(:)
    type(Atom), ALLOCATABLE                 :: MoleculeNew(:)
    type(Parameter), INTENT(IN)             :: Variable
    type(Bonds), INTENT(IN), ALLOCATABLE    :: BondsArray(:)
    real*8, INTENT(IN), ALLOCATABLE         :: AnglesArray(:), TorsionAnglesArray(:)
    real*8                                  :: q
    real*8                                  :: DeltaEnergy, NewEnergy, Probability
    integer, INTENT(OUT)                    :: TotalSamples, AcceptedSamples 
    integer                                 :: RejectionRate

    TotalSamples = 0
    AcceptedSamples = 0
    RejectionRate = 0

    allocate(MoleculeNew(size(Molecule)))
    MoleculeNew = Molecule
    do
        Call NewCoordinates(MoleculeNew, Variable)
        NewEnergy = TotalEnergy(MoleculeNew, Variable, BondsArray, AnglesArray, TorsionAnglesArray)
        DeltaEnergy = NewEnergy - Energy
        if (DeltaEnergy < 0) then
            Molecule = MoleculeNew
            Energy = NewEnergy 
            AcceptedSamples = AcceptedSamples + 1       
        else
            Probability = exp((-DeltaEnergy)/(Variable%Boltzmann*Variable%Temp))
            call random_number(q)
            if (q < Probability) then
                Molecule = MoleculeNew
                Energy = NewEnergy
                AcceptedSamples = AcceptedSamples + 1             
            endif
        endif
        TotalSamples = TotalSamples + 1
        if (AcceptedSamples > 1) then
            RejectionRate = (AcceptedSamples+TotalSamples)/AcceptedSamples
            if (RejectionRate >1000) EXIT
        endif
    enddo
end subroutine

subroutine NewCoordinates(MoleculeNew, Variable)
!Subroutine NewCoordinates calculates the new random coordinates for one atom at the time.
    type(Atom), INTENT(INOUT), ALLOCATABLE  :: MoleculeNew(:)
    type(Parameter), INTENT(IN)             :: Variable
    real*8, DIMENSION(3)                    :: NewCoord
    real*8, DIMENSION(3)                    :: Random
    real*8                                  :: RandomI
    integer                                 :: i
    integer                                 :: j 

    do 
        call random_number(RandomI)
        i = int(RandomI * size(MoleculeNew))
        if (i > 0) exit
    enddo

    call random_number(NewCoord(1))
    call random_number(NewCoord(2))
    call random_number(NewCoord(3))
    do j = 1,3
        call random_number(Random(j))
        if(Random(j) <= 0.5)then
            NewCoord(j) = -NewCoord(j)
        endif
    enddo

    MoleculeNew(i)%x = MoleculeNew(i)%x + Variable%r*NewCoord(1)
    MoleculeNew(i)%y = MoleculeNew(i)%y + Variable%r*NewCoord(2)
    MoleculeNew(i)%z = MoleculeNew(i)%z + Variable%r*NewCoord(3)
end subroutine
end module MinimizationModule