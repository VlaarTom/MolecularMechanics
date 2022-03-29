module AtomModule
implicit none
PRIVATE
PUBLIC Atom, Parameter, InputReader, ParamaterInput

type Atom
    Character(1)    :: Element
    real*8          :: x, y, z
end type

type Parameter
    real*8 :: ForceConstantCC, ForceConstantCH, CCBond, CHBond, ForceConstantAngle, EquiAngle, &
    V1, n, gamma, CCharge, HCharge, CoulombConstant, Avogadro, CWellDepth, HWellDepth, CSigma, HSigma, &
    r, Boltzmann, Temp
end type

contains
subroutine InputReader (filename, Molecule)
!Subroutine InputReader reads the text file and stores the coordinates in the Molecule array.
    CHARACTER(len=*), INTENT(IN) :: filename
    type(Atom), INTENT(OUT), ALLOCATABLE :: Molecule (:)
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

    Variable%ForceConstantCC    = 317                       !kcal/(mol Angstrom**2)
    Variable%ForceConstantCH    = 300                       !?????
    Variable%CCBond             = 1.507                     !Angstrom
    Variable%CHBond             = 1.094                     !Angstrom
    Variable%ForceConstantAngle = 150                       !Unknown
    Variable%EquiAngle          = 2                         !Unknown
    Variable%V1                 = 50                        !Size of the barrier for rotation
    Variable%n                  = 1                         !Potential shift 
    Variable%gamma              = 0                         !Unknown
    Variable%CCharge            = 1.602176634/(10.0**19)    !Unknown
    Variable%HCharge            = 1.602176634/(10.0**19)    !Unknown
    Variable%CoulombConstant    = 8.9875517923*(10**9)      !kg*m**3*S**-2*C**-2
    Variable%Avogadro           = 6.02214076*(10.0**23)     !mol**-1
    Variable%CWellDepth         = 0.066                     !kcal/mol
    Variable%HWellDepth         = 0.030                     !kcal/mol
    Variable%CSigma             = 1                         !Unknown
    Variable%HSigma             = 1                         !Unknown
    Variable%r                  = 0.00000001                !Search radius
    Variable%Boltzmann          = 1.3806485/(10.0**23)      !m**2*kg/s**-2*K**-1
    Variable%Temp               = 300                       !K
end subroutine
end module AtomModule