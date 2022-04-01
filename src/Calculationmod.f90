module CalculationModule
    use AtomModule
implicit none
PRIVATE
PUBLIC Bonds, BondLength, BendAngle, AssignBonds, AddBond, AddAngle, TorsionAngle

type Bonds
    CHARACTER(2)    :: Elements
    real*8          :: Length
    integer         :: Atom1, Atom2
end type

contains
real*8 function BondLength (Molecule, Atom1, Atom2)
!Bondlength Function calculates the distance between two atoms using their cartesian coordinates 
    type(Atom), INTENT(IN), ALLOCATABLE  :: Molecule(:)
    real*8                               :: DeltaX, DeltaY, DeltaZ
    integer, INTENT(IN)                  :: Atom1, Atom2

    DeltaX = Molecule(Atom1)%x - Molecule(Atom2)%x
    DeltaY = Molecule(Atom1)%y - Molecule(Atom2)%y
    DeltaZ = Molecule(Atom1)%z - Molecule(Atom2)%z

    Bondlength = sqrt(DeltaX**2 + DeltaY**2 + DeltaZ**2)
end function

subroutine AssignBonds (Molecule, BondsArray, Bond, Variable)
!AssignBonds subroutine determines which atoms are bonded with a certain threshold and adds the bonds 
!to the BondsArray. The CheckBonds array stores which bonds have been made so that no duplicates occur.        
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
!Subroutine AddBond stores the bonds made in the AssignBonds subroutine in an array.
!A DummyArray is created to temporarily store the added bond.
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

subroutine BendAngle (Molecule, BondsArray, AnglesArray)
!Subroutine BendAngle calculates the angle between three atoms. The first if statement makes sure the angle is being calculated 
!with a carbon atom as centre atom. The other statements make sure that the other two atoms are connected to the same centre atom.
    type(Atom), INTENT(IN), ALLOCATABLE :: Molecule(:)
    type(Bonds), INTENT(IN)             :: BondsArray(:)
    real*8, ALLOCATABLE                 :: AnglesArray(:) 
    real*8                              :: Angle
    integer                             :: i, j, k

    do i = 1,size(Molecule)
        if (Molecule(i)%element == 'C') then
           do j = 1,size(BondsArray)
              do k = j,size(BondsArray)
                 if (j /= k .and. BondsArray(j)%Atom1 == i .and. BondsArray(k)%Atom1 == i) then   
                    Angle = -((BondLength(Molecule, BondsArray(j)%Atom2, BondsArray(k)%Atom2))**2 - &
                            (BondsArray(j)%length)**2 - (BondsArray(k)%length)**2) / (2*(BondsArray(j)%length)*   &
                            (BondsArray(k)%length))
                    call AddAngle(Angle, AnglesArray)
                 endif
                 if (j /= k .and. BondsArray(j)%Atom2 == i .and. BondsArray(k)%Atom1 == i) then
                    Angle = -((BondLength(Molecule, BondsArray(j)%Atom2, BondsArray(k)%Atom2))**2 - &
                            ((BondsArray(j)%length)**2 - (BondsArray(k)%length)**2) / (2*(BondsArray(j)%length)*  &
                            (BondsArray(k)%length)))
                    call AddAngle(Angle, AnglesArray)
                 endif
              enddo
           enddo
        endif
     enddo
end subroutine

subroutine AddAngle (Angle, AnglesArray)
!Subroutine AddAngle stores the angle made in the BendAngle subroutine in an array.
!A DummyArray is created to temporarily store the added angle.
    real*8, INTENT(IN)                 :: Angle
    real*8, INTENT(INOUT), ALLOCATABLE :: AnglesArray(:)
    real*8, ALLOCATABLE                :: DummyArray(:)
    integer                                 :: i, j
            
    if (allocated(AnglesArray)) then
        i = size(AnglesArray)
        allocate(DummyArray(i+1))
        do j = 1, i
            DummyArray(j) = AnglesArray(j)
        enddo
        DummyArray(i+1) = Angle
        deallocate(AnglesArray)
        call move_alloc (DummyArray, AnglesArray)
    else
        allocate(AnglesArray(1))
        AnglesArray(1) = Angle
    endif
end subroutine
   
subroutine TorsionAngle (Molecule, BondsArray, TorsionAnglesArray)
!Subroutine TorsionAngle calculates the torsionangle over four atoms, with two carbon atoms 
!as the centre atoms. First, via the if statements is determined which connections can be made and makes sure no duplicate 
!angles are calculated. Of these connections, vectors are made which are used to calculate the torsion angle which are stored 
!in an array.
    type(Atom), INTENT(IN), ALLOCATABLE :: Molecule(:)
    type(Bonds), INTENT(IN)             :: BondsArray(:)
    real*8, ALLOCATABLE                 :: TorsionAnglesArray(:)
    real*8                              :: Phi, x, y
    real*8, dimension(3)                :: b1, b2, b3, n1, n2, m1
    integer                             :: i, j, k

    do i = 1,size(BondsArray)
        if (BondsArray(i)%Elements == 'CC') then
            do j = 1,size(BondsArray)
                do k = j,size(BondsArray)
                    if (j /= k .and. BondsArray(j)%Atom1 == bondsArray(i)%Atom1 &
                    .and. BondsArray(j)%Atom2 /= bondsArray(i)%Atom2 &
                    .and. BondsArray(k)%Atom1 == BondsArray(i)%Atom2 & 
                    .and. BondsArray(k)%Atom2 /= BondsArray(i)%Atom1) then 

                        b1(1) = Molecule(BondsArray(j)%Atom2)%x - Molecule(BondsArray(i)%Atom1)%x
                        b1(2) = Molecule(BondsArray(j)%Atom2)%y - Molecule(BondsArray(i)%Atom1)%y
                        b1(3) = Molecule(BondsArray(j)%Atom2)%z - Molecule(BondsArray(i)%Atom1)%z

                        b2(1) = Molecule(BondsArray(i)%Atom1)%x - Molecule(BondsArray(i)%Atom2)%x
                        b2(2) = Molecule(BondsArray(i)%Atom1)%y - Molecule(BondsArray(i)%Atom2)%y
                        b2(3) = Molecule(BondsArray(i)%Atom1)%z - Molecule(BondsArray(i)%Atom2)%z

                        b3(1) = Molecule(BondsArray(k)%Atom2)%x - Molecule(BondsArray(i)%Atom2)%x
                        b3(2) = Molecule(BondsArray(k)%Atom2)%y - Molecule(BondsArray(i)%Atom2)%y
                        b3(3) = Molecule(BondsArray(k)%Atom2)%z - Molecule(BondsArray(i)%Atom2)%z
                        
                        n1 = CrossProduct(b1, b2)
                        n2 = CrossProduct(b2, b3)
                        m1 = CrossProduct(n1, b2)
                        x = DOT_PRODUCT(n1, n2)
                        y = DOT_PRODUCT(m1, n2)
                        Phi = atan2(y,x)
                        call AddTorsionAngle (Phi, TorsionAnglesArray)
                         
                    elseif (j /= k .and. BondsArray(j)%Atom2 == bondsArray(i)%Atom1 &
                    .and. BondsArray(j)%Atom2 /= bondsArray(i)%Atom2 &
                    .and. BondsArray(k)%Atom1 == BondsArray(i)%Atom2 & 
                    .and. BondsArray(k)%Atom2 /= BondsArray(i)%Atom1) then

                        b1(1) = Molecule(BondsArray(j)%Atom2)%x - Molecule(BondsArray(i)%Atom1)%x
                        b1(2) = Molecule(BondsArray(j)%Atom2)%y - Molecule(BondsArray(i)%Atom1)%y
                        b1(3) = Molecule(BondsArray(j)%Atom2)%z - Molecule(BondsArray(i)%Atom1)%z

                        b2(1) = Molecule(BondsArray(i)%Atom1)%x - Molecule(BondsArray(i)%Atom2)%x
                        b2(2) = Molecule(BondsArray(i)%Atom1)%y - Molecule(BondsArray(i)%Atom2)%y
                        b3(3) = Molecule(BondsArray(i)%Atom1)%z - Molecule(BondsArray(i)%Atom2)%z

                        b3(1) = Molecule(BondsArray(k)%Atom2)%x - Molecule(BondsArray(i)%Atom2)%x
                        b3(2) = Molecule(BondsArray(k)%Atom2)%y - Molecule(BondsArray(i)%Atom2)%y
                        b3(3) = Molecule(BondsArray(k)%Atom2)%z - Molecule(BondsArray(i)%Atom2)%z

                        n1 = CrossProduct(b1, b2)
                        n2 = CrossProduct(b2, b3)
                        m1 = CrossProduct(n1, b2)
                        x = DOT_PRODUCT(n1, n2)
                        y = DOT_PRODUCT(m1, n2)
                        Phi = atan2(y,x)
                        call AddTorsionAngle (Phi, TorsionAnglesArray)
    
                    endif
                enddo
            enddo
        endif
    enddo
end subroutine

function CrossProduct (x, y)
!The CrossProduct function calculates the cross product of two vectors  
    real*8, DIMENSION(3)    :: CrossProduct, x, y

    CrossProduct(1) = (x(2)*y(3)) - (x(3)*y(2))
    CrossProduct(2) = (x(1)*y(3)) - (x(3)*y(1))
    CrossProduct(3) = (x(1)*y(2)) - (x(2)*y(1))
end function

subroutine AddTorsionAngle (Phi, TorsionAnglesArray)
!Subroutine AddTorsionAngle stores the angle made in the TorsionAngle subroutine in an array.
!A DummyArray is created to temporarily store the added torsion angle.
    real*8, INTENT(IN)                 :: Phi
    real*8, INTENT(INOUT), ALLOCATABLE :: TorsionAnglesArray(:)
    real*8, ALLOCATABLE                :: DummyArray(:)
    integer                            :: i, j
                
    if (allocated(TorsionAnglesArray)) then
        i = size(TorsionAnglesArray)
        allocate(DummyArray(i+1))
        do j = 1, i
            DummyArray(j) = TorsionAnglesArray(j)
        enddo
        DummyArray(i+1) = Phi
        deallocate(TorsionAnglesArray)
        call move_alloc (DummyArray, TorsionAnglesArray)
    else
        allocate(TorsionAnglesArray(1))
        TorsionAnglesArray(1) = Phi
    endif
end subroutine
end module CalculationModule