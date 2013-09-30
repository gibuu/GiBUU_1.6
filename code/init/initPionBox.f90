!*****************************************************************************
!****m* /initPionBox
! NAME
! module initPionBox
! PURPOSE
! Initializes pions for a box of pions
!*****************************************************************************
module initPionBox

  implicit none

  Private

  !*************************************************************************
  !****g* initPionBox/nDens
  ! SOURCE
  !
  real, save :: nDens = 1.0
  ! PURPOSE
  ! particle density [fm^-3]
  !*************************************************************************

  !*************************************************************************
  !****g* initPionBox/ChargeSelection
  ! SOURCE
  !
  integer, save :: ChargeSelection = 0
  ! PURPOSE
  ! define the type of the charge selection:
  ! * 0: only pi0
  ! * 1: 50% pi+, 50% pi-
  ! * 2: 33% for +,0,-
  !*************************************************************************

  !*************************************************************************
  !****g* initPionBox/pInit
  ! SOURCE
  !
  real, save :: pInit = 0.5
  ! PURPOSE
  ! initial momentum of particles [GeV/c]
  !*************************************************************************

  Public :: initializePionBox

contains

  !*************************************************************************
  !****s* initPionBox/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'initBox'.
  !*************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput
    use callstack, only: traceBack

    !***********************************************************************
    !****n* initPionBox/PionBox
    ! NAME 
    ! NAMELIST PionBox
    ! PURPOSE
    ! Includes the input parameters:
    ! * nDens
    !***********************************************************************
    NAMELIST /PionBox/ nDens,ChargeSelection,pInit

    integer :: ios

    call Write_ReadingInput('PionBox',0)
    rewind(5)
    read(5,nml=PionBox,iostat=ios)
    call Write_ReadingInput('HiLeptonNucleus',0,ios)

    write(*,'(A,F8.3,A)') ' particle density=',nDens,' fm^(-3)'
    write(*,*) 'charge selection: ',ChargeSelection
    write(*,'(A,F8.3,A)') ' initial momentum=',pInit,' GeV/c'

    call Write_ReadingInput('PionBox',1)

  end subroutine initInput

  !*************************************************************************
  !****s* initPionBox/initializePionBox
  ! NAME
  ! subroutine initializePionBox(part)
  ! PURPOSE
  ! Initialize nucleons in a box
  !*************************************************************************
  subroutine initializePionBox(part)

    use particleDefinition
    use IdTable, only: pion
    use constants, only: mPi
    use densityModule, only: gridsize,get_densitySwitch
    use output, only: Write_InitStatus
    use callstack, only: traceBack
    use insertion, only: GarbageCollection
    use collisionNumbering, only: real_firstnumbering

    type(particle), dimension(:,:),intent(inOut) :: part

    integer :: numberPions, dummy, iEns,iPart, nEns
    real, dimension(1:3) :: pSum

    call Write_InitStatus('box of pions',0)
    call initInput

    dummy = get_densitySwitch() ! force density module to be initialized

    numberPions = NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*nDens)
    nEns = size(part(:,1))

    write(*,*) ' Number of pions per ensemble: ', numberPions
    write(*,'(A,3F9.3)') '  Gridsize   =',gridsize(:)
    write(*,*) ' Size of box= (8*Gridsize) = ',8.*gridsize(1)*gridsize(2)*gridsize(3)
    write(*,*) ' Number Ensembles             =',nEns
    write(*,*) ' Number Particles per Ensemble=',size(part(1,:))

    if (numberPions > size(part(1,:))) then
       call traceback('particle vector too small!')
    end if

    do iEns=1,nEns
       pSum = 0
       do iPart=1,numberPions
          call setToDefault(part(iEns,iPart))
          call setNumber(part(iEns,iPart)) ! give it a unique number

          part(iEns,iPart)%event = real_firstnumbering()
          part(iEns,iPart)%ID = pion
          part(iEns,iPart)%mass = mPi
          call ChooseCharge
          call ChoosePosition
          call ChooseMomentum

          pSum = pSum + part(iEns,iPart)%momentum(1:3)

       end do

       pSum = pSum/numberpions

       ! correct for 'resting box'

!!$       do iPart=1,numberPions
!!$          part(iEns,iPart)%momentum(1:3) = part(iEns,iPart)%momentum(1:3) - pSum
!!$          part(iEns,iPart)%momentum(0) = sqrt(mPi**2+sum(part(iEns,iPart)%momentum(1:3)**2))
!!$       end do

    end do
    
    call GarbageCollection(part)
    call Write_InitStatus('box of pions',1)

  contains

    subroutine ChooseCharge
      use random, only: rn
      real :: r

      select case (ChargeSelection)
      case (0)

         part(iEns,iPart)%charge = 0

      case (1)

         r = rn() * 2.0
         if (r < 1.0) then
            part(iEns,iPart)%charge = -1
         else
            part(iEns,iPart)%charge =  1
         end if

      case (2)
         
         r = rn() * 3.0
         if (r < 1.0) then
            part(iEns,iPart)%charge = -1
         else if (r < 2.0) then
            part(iEns,iPart)%charge =  0
         else
            part(iEns,iPart)%charge =  1
         end if

      case DEFAULT
         write(*,*) 'wrong ChargeSelection = ',ChargeSelection
         call traceback('correct input')
      end select

    end subroutine ChooseCharge

    subroutine ChoosePosition
      use random, only: rn

      part(iEns,iPart)%position(1)=(1.-2.*rn())*gridSize(1)
      part(iEns,iPart)%position(2)=(1.-2.*rn())*gridSize(2)
      part(iEns,iPart)%position(3)=(1.-2.*rn())*gridSize(3)

    end subroutine ChoosePosition

    subroutine ChooseMomentum
      use random, only: rnOmega
      
      part(iEns,iPart)%momentum(1:3) = rnOmega() * pInit
      part(iEns,iPart)%momentum(0) = sqrt(mPi**2+pInit**2)

    end subroutine ChooseMomentum

  end subroutine initializePionBox

end module initPionBox
