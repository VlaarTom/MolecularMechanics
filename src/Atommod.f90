module AtomModule
implicit none
PRIVATE
PUBLIC Atom, Parameter, InputReader, ParamaterInput

type Atom
    Character(1)    :: Element
    real*8          :: x, y, z
end type

type Parameter
    real*8 :: ForceConstantCC, ForceConstantCH, CCBond, CHBond
end type

contains
subroutine InputReader (filename, Molecule)
    CHARACTER(len=*), INTENT(IN) :: filename
    type(Atom), INTENT(INOUT), ALLOCATABLE :: Molecule (:)
    integer :: fu, AtomNumber, i

    open(newunit = fu, file = filename, action = 'read')
    read(fu,*) AtomNumber
    allocate(Molecule(AtomNumber))

    do i = 1, AtomNumber
    read(fu,*) Molecule(i)%Element, Molecule(i)%x, Molecule(i)%y, Molecule(i)%z
    enddo

    close(fu)
end subroutine

subroutine ParamaterInput (Variable)
    type(Parameter), INTENT(INOUT) :: Variable

    Variable%ForceConstantCC    = 317   !kcal/(mol Angstrom**2)
    Variable%ForceConstantCH    = 300   !?????
    Variable%CCBond             = 1.507 !Angstrom
    Variable%CHBond             = 1.094 !Angstrom
end subroutine

end module AtomModule