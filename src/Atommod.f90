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

    Variable%ForceConstantCC    = 1.326e6                   !J/mol
    Variable%ForceConstantCH    = 1.255e6                   !J/mol
    Variable%CCBond             = 1.507                     !Angstrom
    Variable%CHBond             = 1.094                     !Angstrom
    Variable%ForceConstantAngle = 475720.8                  !J/mol
    Variable%EquiAngle          = 2                         !Unknown
    Variable%V1                 = 5857.6                    !J/mol 
    Variable%n                  = 2                         !Potential shift 
    Variable%gamma              = 0                         !Unknown
    Variable%CCharge            = -0.344*(1.602176634e-19)  !Charge
    Variable%HCharge            = 0.078*(1.602176634e-19)   !Charge
    Variable%CoulombConstant    = 8.9875517923e9            !kg*m**3*S**-2*C**-2
    Variable%Avogadro           = 6.02214076e23             !mol e-1
    Variable%CWellDepth         = 457.7296                  !J/mol
    Variable%HWellDepth         = 65.6888                   !J/mol
    Variable%CSigma             = 1.9080                    !Unknown
    Variable%HSigma             = 0.6000                    !Unknown
    Variable%r                  = 1e-9                      !Search radius
    Variable%Boltzmann          = 1.3806485e-23             !m**2*kg/s**-2*K**-1
    Variable%Temp               = 550                       !K
end subroutine
end module AtomModule