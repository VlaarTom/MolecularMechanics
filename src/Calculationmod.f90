module CalculationModule

    use AtomModule

implicit none
PRIVATE
PUBLIC Bonds, BondLength, BendAngle, AssignBonds, AddBond

type Bonds
    CHARACTER(2)    :: Elements
    real*8          :: Length
    integer         :: Atom1, Atom2
end type    

contains
real*8 function BondLength (Molecule, Atom1, Atom2)
!Bondlength Function calculates the distance between two atoms 

    type(Atom), INTENT(IN), ALLOCATABLE  :: Molecule(:)
    real*8                               :: DeltaX, DeltaY, DeltaZ
    integer, INTENT(IN)                  :: Atom1, Atom2

    DeltaX = Molecule(Atom1)%x - Molecule(Atom2)%x
    DeltaY = Molecule(Atom1)%Y - Molecule(Atom2)%Y
    DeltaZ = Molecule(Atom1)%Z - Molecule(Atom2)%z

    Bondlength = sqrt(DeltaX**2 + DeltaY**2 + DeltaZ**2)
end function


subroutine AssignBonds (Molecule, BondsArray, Bond, Variable)
    !AssignBonds subroutine determines which atoms are bonded and adds the bonds to the BondsArray.
    !The CheckBonds array stores which bonds have been made so that no duplicates occur.
            
    type(Atom), INTENT(IN), ALLOCATABLE      :: Molecule(:)
    type(Parameter), INTENT(IN)              :: Variable
    type(Bonds), INTENT(INOUT), ALLOCATABLE  :: BondsArray(:)
    type(Bonds)                              :: Bond
    integer                                  :: Atom1, Atom2
    integer, ALLOCATABLE                     :: CheckBonds(:,:)
        
    allocate(CheckBonds(size(Molecule), size(Molecule)))
    CheckBonds = 0
        
    do Atom1 = 1, size(Molecule)
        do Atom2 = 1, size(Molecule)
            if (Atom1 /= Atom2) then
                if (Molecule(Atom1)%Element == "C" .and. Molecule(Atom2)%Element == "C" &
                    .and. (abs(Bondlength(Molecule, Atom1, Atom2)) - Variable%CCBond  < 0.2) &
                    .and. CheckBonds(Atom1, Atom2) == 0) then
                            
                    Bond%Elements = "CC"
                    Bond%Length = BondLength(Molecule, Atom1, Atom2)
                    Bond%Atom1 = Atom1
                    Bond%Atom2 = Atom2
                    CheckBonds(Atom1, Atom2) = 1
                    CheckBonds(Atom2, Atom1) = 1
                    call AddBond (Bond, BondsArray)
                            
                elseif (Molecule(Atom1)%Element == "C" .and. Molecule(Atom2)%Element == "H" &
                .and. abs(Bondlength(Molecule, Atom1, Atom2)) - Variable%CHBond < 0.2) then
        
                    Bond%Elements = "CH"
                    Bond%Length = BondLength(Molecule, Atom1, Atom2)
                    Bond%Atom1 = Atom1
                    Bond%Atom2 = Atom2
                    Call AddBond (Bond, BondsArray)
                endif
            endif
        enddo
    enddo
end subroutine
        

subroutine AddBond (Bond, BondsArray)
!Subroutine AddBond stores the bonds made in the AssignBonds subroutine in an array
!A DummyArray is created to temporarily store the added bond
    type(Bonds), INTENT(IN)                 :: Bond
    type(Bonds), INTENT(INOUT), ALLOCATABLE :: BondsArray(:)
    type(bonds), ALLOCATABLE                :: DummyArray(:)
    integer                                 :: i, j
        
    if (allocated(BondsArray)) then
        i = size(BondsArray)
        allocate(DummyArray(i+1))
        do j = 1, i
            DummyArray(j) = BondsArray(j)
        enddo
        DummyArray(i+1) = Bond
        deallocate(BondsArray)
        call move_alloc (DummyArray, BondsArray)
    else
        allocate(BondsArray(1))
        BondsArray(1) = Bond
    endif
end subroutine


subroutine BendAngle (molecule, BondsArray)!, Angles)
    type(Atom), INTENT(IN), ALLOCATABLE :: Molecule(:)
    type(Bonds), INTENT(IN)             :: BondsArray(:)
    real*8, ALLOCATABLE                 :: Angles(:) 
    real*8                              :: Angle
    integer                             :: i, j, k, l

    do i = 1, size(Molecule)
        if (Molecule(i)%Element == "C") then
            do j = 1, size(BondsArray)
                do k = j, size(BondsArray)
                    if (j /= k .and. BondsArray(i)%Atom1 == i .and. BondsArray(k)%Atom1 == i) then
                        print *, Molecule(i)%Element, Molecule(j)%Element, Molecule(k)%Element
                        
                    endif
                    if (j /= k .and. BondsArray(i)%Atom2 == i .and. BondsArray(k)%Atom1 == i) then
                        print *, Molecule(i)%Element, Molecule(j)%Element, Molecule(k)%Element
                    endif
                enddo
            enddo
        endif
    enddo
end subroutine
end module CalculationModule