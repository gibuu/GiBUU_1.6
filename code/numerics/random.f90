!***************************************************************************
!****m* /random
! NAME
! module random
! PURPOSE
! Contains all information and routines, which are necessary for the
! random number generation.
!
! NOTES
! Formerly we used the random number generator RAN3 from Numerical Recipes.
! We switched to the generator "PYR" provided by Pythia, which is a rewrite
! of the RANMAR routine from the Cernlib.
!
! Contrary to the old RAN3 behaviour, the new rn() excludes the endpoints
! 0 and 1 !
!***************************************************************************
module random
  implicit none
  private


  !***********************************************************************
  !****g* random/Seed
  ! PURPOSE
  ! Random Seed (used to initialize the random number generator),
  ! accessible through the namelist 'initRandom'.
  ! If Seed is zero (default), then it is set via "SYSTEM_CLOCK()".
  !
  ! SOURCE
  !
  integer,save ::  Seed = 0
  !***********************************************************************


  !***********************************************************************
  !****g* random/resetRandom
  ! PURPOSE
  ! Reread random generator, used by setRandom, useful for debugging.
  !
  ! SOURCE
  !
  logical,save :: resetRandom = .false.
  !***********************************************************************


  !*************************************************************************
  !****f* random/rnExp
  ! NAME
  ! real function rnExp(A)
  ! real function rnExp(A,x1,x2)
  ! PURPOSE
  ! Evaluates a random number x according to exp(-A*x) distribution.
  ! If given, x is restricted to lie between x1 and x2
  ! INPUTS
  ! * real :: A -- slope parameter
  ! * real :: x1 -- minimal x-value
  ! * real :: x2 -- maximal x-value
  ! RESULT
  ! random number
  !*************************************************************************
  interface rnExp
     module procedure rnExp1, rnExp2
  end interface

  logical,save :: first=.true.

  public :: rn,rnFlat,rnExp,rnOmega,rnOmega_anis,rnGauss,rnCos,rn_openInterval,rnOmega_angles,setRandom, rn_truefalse

contains


  !*************************************************************************
  !****f* random/InitRan
  ! NAME
  ! subroutine InitRan
  ! PURPOSE
  ! Reads random seed out of namelist 'initRandom' in jobcard and initializes
  ! the random number generator
  !*************************************************************************
  subroutine InitRan
    use output, only: Write_ReadingInput

    integer :: ios
    integer(8) :: seed8, rate

    !***********************************************************************
    !****n* random/initRandom
    ! NAME
    ! NAMELIST /initRandom/
    ! PURPOSE
    ! Includes the input variables:
    ! * Seed
    ! * resetRandom
    !***********************************************************************
    NAMELIST /initRandom/ Seed,resetRandom

    call Write_ReadingInput('initRandom',0)
    rewind(5)
    read(5,nml=initRandom,IOSTAT=ios)
    call Write_ReadingInput('initRandom',0,ios)

    if (Seed == 0) then
      call system_clock (seed8, rate)              ! get 8-byte clock value
      write(*,*) 'Resetting Seed via system clock: ', seed8
      write(*,*) '                      precision: ',rate,'counts/sec'
      seed = int(mod(seed8,int(huge(seed),8)),4)   ! convert to 4-byte seed
    end if

    write(*,*) 'Seed: ', Seed
    write(*,*) 'resetRandom: ', resetRandom

    call InitPYR (Seed)

    call Write_ReadingInput('initRandom',1)

  end subroutine InitRan


  !*************************************************************************
  !****s* random/SetRandom
  ! NAME
  ! subroutine SetRandom
  ! PURPOSE
  ! write out/read random number generators
  ! NOTES
  ! This routine provides a shortcut for MC debugging:
  ! if an error occurs in run 12345 of 99999 runs/energy after 3 days
  ! of CPU time, the idea is to reset the random generators by reading
  ! in some files, that already run 1 reproduces this error.
  !
  ! In order to avoid some re-ordering of the random number lists, only
  ! in run number 3 the previously written files are read.
  ! The code stops in run 4.
  !*************************************************************************
  subroutine setRandom

    logical,external :: ReadPYR
    integer,save :: nRead = 0

    common/DebugWRT/ DoWRT
    logical DoWRT
    save/DebugWRT/

    if (first) then
       ! Initializes the random number generator
       ! when it's called for the first time
       call InitRan
       first=.false.
    end if

    if (resetRandom) then
       nread = nread+1
       DoWRT = .false.
       if (nRead.lt.3) return
       DoWRT = .true.

       write(*,*) '=== ATTENTION: reset of random generators!!! ==='

       if (nRead.gt.3) then
          write(*,*) 'nRead>3. stop.'
          stop
       endif

       if (.not.ReadPYR()) then
          write(*,*) 'ReadPYR failed. stop.'
          stop
       endif
    else
       call DumpPYR()  ! write the PYTHIA ran.gen. to file
    endif

    !    write(*,*) rn(0), pyr(0)

  end subroutine setRandom


  !*************************************************************************
  !****f* random/rn_openInterval
  ! NAME
  ! real function rn_openInterval()
  ! PURPOSE
  ! evaluates a random number in (0,1)
  ! USAGE
  ! (real)=rn_openInterval()
  ! Notes
  ! Finds random number which is in [1E-8,1-1E-8], which is approximately (0,1).
  !*************************************************************************
  real function rn_openInterval()
    real,parameter :: eps=1E-8
    do
       ! Find random number which is not 1 or 0:
       rn_openInterval=rn()
       if(abs(rn_openInterval-1.).gt.eps.and.abs(rn_openInterval).gt.eps) return
    end do
  end function rn_openInterval


  !*************************************************************************
  !****f* random/rn
  ! NAME
  ! real function rn()
  ! PURPOSE
  ! evaluates a random number in (0,1)
  ! USAGE
  ! (real)=rn()
  ! NOTES
  ! ckecks whether the random generator should be initialised, otherwise
  ! it just calls PYR. Endpoints are excluded.
  !*************************************************************************
  real function rn()
    real,external :: PYR ! prototype

    if (first) then
       ! Initializes the random number generator
       ! when it's called for the first time
       call InitRan
       first=.false.
    end if

    rn=PYR(0)
    return
  end function rn


  !*************************************************************************
  !****f* random/rnFlat
  ! NAME
  ! real function rnFlat(xMin,xMax)
  ! PURPOSE
  ! draws a random number from a flat (uniform) distribution on (xMin,xMax)
  ! INPUTS
  ! * real :: xMin,xMax -- minimum and maximum value
  !*************************************************************************
  real function rnFlat (xMin,xMax)
    real, intent(in) :: xMin, xMax
    rnFlat = XMin + rn()*(xMax-xMin)
  end function rnFlat


  !*************************************************************************
  !****f* random/rnGauss
  ! NAME
  ! real function rnGauss(StdDev, Mean)
  ! PURPOSE
  ! evaluates a random number in according a Gauss distribution
  ! INPUTS
  ! * real :: StdDev -- standard deviation
  ! * real :: Mean   -- mean value of distribution
  ! RESULT
  ! random number
  !*************************************************************************
  real FUNCTION rnGauss(StdDev, Mean)
    real,intent(in) :: StdDev, Mean
    real,save :: x1,x2
    integer,save :: isCached = 0
    real v1,v2, w, y

    if (isCached.ne.0) then
       isCached = 0
       rnGauss = x2*StdDev+Mean
       return
    endif

    do
      v1 = 2.0*rn() - 1.0
      v2 = 2.0*rn() - 1.0
      w = v1**2 + v2**2
      if (w.le.1.0) exit
    end do

    y = sqrt( (-2.0*log(w))/w )
    x1 = v1 * y
    x2 = v2 * y

    isCached = 1
    rnGauss = x1*StdDev+Mean
  end FUNCTION rnGauss


  !*************************************************************************
  !****f* random/rnCos
  ! NAME
  ! real function rnCos()
  ! PURPOSE
  ! Evaluates a random number x according to cos(x) distribution.
  ! So cos(x) is assumed to be isotropic in [-1,1]
  ! INPUTS
  ! * NONE
  ! RESULT
  ! random number
  !*************************************************************************
  FUNCTION rnCos() Result(x)
    real :: x
    real :: cosX
    cosX=1.-2.*rn()
    x=acos(cosX)
  end FUNCTION rnCos


  !*************************************************************************
  ! cf. interface random/rnExp
  !*************************************************************************
  FUNCTION rnExp1(A) Result(x)
    real, intent(in) :: A
    real :: x
    x= -Log(rn())/A
  end FUNCTION rnExp1
  !-------------------------------------------------------------------------
  FUNCTION rnExp2(A,x1,x2) Result(x)
    real, intent(in) :: A,x1,x2
    real :: x,r1,r2
    r1 = exp(-A*x1)
    r2 = exp(-A*x2)
    x = -log(r2+rn()*(r1-r2))/A
  end FUNCTION rnExp2
  !*************************************************************************


  !*************************************************************************
  !****f* random/rnOmega
  ! NAME
  ! function rnOmega()
  ! PURPOSE
  ! Generates a unit vector with random angle (isotropical distribution).
  ! INPUTS
  ! * NONE
  ! RESULT
  ! real :: rnOmega(3) - random 3-dim. unit-vector
  !*************************************************************************
  function rnOmega()
    use constants, only : twopi
    real:: rnOmega(3)
    real:: phi,cost,sint
    phi=twopi*rn()                      ! phi
    cost=2.*rn()-1.                     ! cos(theta)
    sint=sqrt(1.-cost**2)               ! sin(theta)
    rnOmega(1:3) = (/ sint*cos(phi), sint*sin(phi), cost /)
  end function rnOmega


  !*************************************************************************
  !****f* random/rnOmega_anis
  ! NAME
  ! function rnOmega_anis (B)
  ! PURPOSE
  ! Generates a unit vector with random angle (anisotropical distribution).
  ! cos(theta) is chosen according to a probability distribution ~ 1 + B*cos**2(theta).
  ! INPUTS
  ! * real, intent(in) :: B  ---  parameter of prob. distr.
  ! RESULT
  ! real :: rnOmega_anis(3) - random 3-dim. unit-vector
  !*************************************************************************
  function rnOmega_anis (B)
    use constants, only : twopi
    real :: rnOmega_anis(3)                   ! return value: 3-vector
    real, intent(in) :: B                     ! parameter of prob. distr.
    real :: phi,cost,sint,x

    phi=twopi*rn()

    ! Choose cos(theta) according to f(cost)~(1+B*cost**2)
    do
       cost = 1 - 2.*rn()
       x = rn()
       if (x < (1+B*cost**2)/(2*(1+B/3.))) exit
    end do

    sint = sqrt(1.-cost**2)   ! sin(theta)
    rnOmega_anis(1:3) = (/ sint*cos(phi), sint*sin(phi), cost /)

  end function rnOmega_anis


  !*************************************************************************
  !****s* random/rnOmega_angles
  ! NAME
  ! subroutine rnOmega(theta,phi)
  ! PURPOSE
  ! Generates random angles theta and phi (istropical distribution) in degree.
  ! INPUTS
  ! * NONE
  ! OUTPUT
  ! * real :: theta, phi
  !*************************************************************************
  subroutine rnOmega_angles(theta,phi)
    use constants, only: pi
    real,intent(out) :: phi,theta
    real :: cost
    phi=rn()   *360.                   ! phi
    cost=2.*rn()-1.                     ! cos(theta)
    theta=acos(cost)*180./pi
  end subroutine rnOmega_angles
 !*************************************************************************


  !*************************************************************************
  !****f* random/rn_trueFalse
  ! NAME
  ! function rn_trueFalse
  ! PURPOSE
  ! Gives true or false with same propability
  ! INPUTS
  ! * NONE
  ! RESULT
  ! logical
  !*************************************************************************
  function rn_truefalse()
    logical:: rn_trueFalse
    rn_trueFalse=(rn().gt.0.5)
  end function rn_truefalse
  !*************************************************************************




end module random
