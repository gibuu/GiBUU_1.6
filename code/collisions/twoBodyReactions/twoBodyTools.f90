!*******************************************************************************
!****m* /twoBodyTools
! NAME
! module twoBodyTools
! PURPOSE
! This module contains some auxiliary routines for two-body collisions.
!*******************************************************************************
module twoBodyTools

  implicit none
  private

  public :: sqrtS_free, pCM, pCM_sqr, p_lab
  public :: convertToAntiParticles, ranCharge, get_pInitial, velocity_correction
  public :: searchInInput, MomentumTransfer2, IsChargeExchange


  !*************************************************************************
  !****f* twoBodyTools/pCM
  ! NAME
  ! real function pCM(srts, mass1, mass2)
  !
  ! real function pCM(srts, mass1, mass2, flagOK)
  ! PURPOSE
  ! This routine evaluates the CM-momentum of two particles, asssuming 
  ! vacuum kinematics.
  !
  ! The formula used here looks different than the standard formula
  ! given in PDG, but the are actually the same. In particular,
  ! pCM is of course invariant with respect to an exchange of
  ! 'mass1' and 'mass2'.
  !
  ! If the no logical 'falgOK' is given, the routine aims for performance, 
  ! and does not check whether the argument of the sqrt is positive. 
  ! If you need to check this, give this logical. 
  !
  ! INPUTS
  ! * real :: mass1, mass2  ! masses of both particles
  ! * real :: srts          ! SQRT(s)
  !
  ! OUTPUT
  ! * real :: pCM ! center of mass momentum
  !*************************************************************************
  Interface pCM
     Module Procedure pCM_1, pCM_2
  end Interface


contains

  !*************************************************************************
  !****f* twoBodyTools/velocity_correction
  ! NAME
  ! real function velocity_correction(pair)
  ! PURPOSE
  ! Evaluates relative velocity/relative velocity in vacuum
  ! INPUTS
  ! *  type(particle), dimension(1:2) :: pair
  ! OUTPUT
  ! function value
  !*************************************************************************
  real function velocity_correction(pair)
    use particleDefinition
    use lorentzTrafo
    type(particle), dimension(1:2) :: pair
    real, dimension(1:2,0:3) :: velocity
    real, dimension(0:3)     :: dummy, momCM_1,momCM_2
    real, dimension(1:3)     ::beta_12,beta_12_vacuum,beta
    real :: beta_12_squared
    logical, parameter :: debug=.false.

    !*********************************************+
    ! Evaluate real relative velocity in CM-Frame
    !*********************************************+

    ! Define 4-velocities of both particles
    ! First particle:
    velocity(1,0)=1./Sqrt(1.-Dot_Product(pair(1)%velocity(1:3),pair(1)%velocity(1:3)))
    velocity(1,1:3)=pair(1)%velocity(1:3)*velocity(1,0)
    ! Second particle:
    velocity(2,0)=1./Sqrt(1.-Dot_Product(pair(2)%velocity(1:3),pair(2)%velocity(1:3)))
    velocity(2,1:3)=pair(2)%velocity(1:3)*velocity(2,0)

    !Evaluate the velocity of the CM-frame of both particles:
    beta(1:3)=(pair(1)%momentum(1:3)+pair(2)%momentum(1:3))/ &
         &    (pair(1)%momentum(0)+pair(2)%momentum(0))

    !Boost everything to the CM-frame
    ! Need this dummy because otherwise my compiler throws warnings : "In call to LORENTZ, an array temporary was created for argument #2"
    
    dummy=velocity(1,0:3)
    call lorentz(beta,dummy, 'velocity_correction')
    velocity(1,0:3)=dummy

    dummy=velocity(2,0:3)
    call lorentz(beta,dummy, 'velocity_correction')
    velocity(2,0:3)=dummy

    beta_12(1:3)=velocity(1,1:3)/velocity(1,0)-velocity(2,1:3)/velocity(2,0)

    beta_12_squared=Dot_Product(beta_12,beta_12)
    if(beta_12_squared .lt. 1.e-06) then
       velocity_Correction=1.
       return
    end if
   
    !*********************************************+
    ! Evaluate vacuum relative velocity in CM-Frame
    !*********************************************+

    !Boost momenta to the CM-frame
    
    momCM_1=pair(1)%momentum
    call lorentz(beta,momCM_1, 'velocity_correction')

    momCM_2=pair(2)%momentum
    call lorentz(beta,momCM_2, 'velocity_correction')

    beta_12_vacuum(1:3)=momCM_1(1:3)/SQRT(pair(1)%mass**2+Dot_product(momCM_1(1:3),momCM_1(1:3))) &
         &             -momCM_2(1:3)/SQRT(pair(2)%mass**2+Dot_product(momCM_2(1:3),momCM_2(1:3)))

    velocity_Correction=SQRT(Dot_Product(beta_12_vacuum,beta_12_vacuum)/beta_12_squared)

    If(debug) write(*,*) 'Velo-Correction=',velocity_Correction

  end function velocity_correction



  !*************************************************************************
  !****f* twoBodyTools/get_PInitial
  ! NAME
  ! real function get_PInitial(teilchenIn)
  ! PURPOSE
  ! return CM momentum of 2 particles
  ! INPUTS
  ! * type(particle),intent(in),dimension(1:2) :: teilchenIn
  ! * integer, intent(in) :: mode
  ! NOTES
  ! Return value depends on mode:
  ! * mode=0 : getP_Pinitial returns CM-momentum of two particles
  ! * mode=1 : getP_Pinitial returns old ikodama=1 prescription [not yet implemented]
  ! * mode=2 : getP_Pinitial returns old ikodama=2 prescription [not yet implemented]
  !*************************************************************************
  real function get_PInitial(teilchenIn,mode)
    use particleDefinition
    use lorentzTrafo
    type(particle), intent(in), dimension(1:2):: teilchenIn
    integer, intent(in) :: mode
    ! mode=0 ! getP_Pinitial returns CM-momentum of two particles
    ! mode=1 ! getP_Pinitial returns old ikodama=1 prescription [not yet implemented]
    ! mode=2 ! getP_Pinitial returns old ikodama=2 prescription [not yet implemented]
    real, dimension(1:3) :: betaCM
    real, dimension(0:3) :: momentum_Total,mom
    integer :: k
    
    If(mode.eq.0) then
       Do k=0,3
          momentum_Total(k)=Sum(teilchenIn(:)%momentum(k))
       end do
       
       betacm = lorentzCalcBeta (momentum_Total, 'get_PInitial')
       
       mom(0:3)=teilchenIN(1)%momentum(0:3)
       call lorentz(betacm,mom,'twoBodyTools')

       get_PInitial=SQRT(Dot_Product(mom(1:3),mom(1:3)))
    else
       Write(*,*) 'Mode not yet Implemented in getPinitial:', mode
       stop
    end if

  end function get_PInitial



  !*************************************************************************
  !****f* twoBodyTools/sqrtS_Free
  ! NAME
  ! function sqrtS_Free(teilchen)
  ! PURPOSE
  ! return the "free" sqrt(s) of 2 particles
  ! INPUTS
  ! * type(particle),intent(in),dimension(1:2) :: teilchen
  ! NOTES
  ! The name "free Sqrt(s)" is according to standard BUU prescription the following:
  ! Transform everything to CM system and define there the sqrt(s) by neglecting the 
  ! potentials. Therefore one wants to take care of the fact, that also the Xsections 
  ! are measured in vacuum, therefore without potentials.
  ! RESULT
  ! * Real : The free Sqrt(s)
  !*************************************************************************
  function sqrtS_Free(teilchen)
    use particleDefinition
    use lorentzTrafo
    real :: sqrtS_Free
    type(particle),intent(in),dimension(1:2) :: teilchen
    real, dimension(1:3) :: betaCM
    real, dimension(0:3) :: momentum_Total,mom
    integer :: k
    real :: mom2

    Do k=0,3
       momentum_Total(k)=Sum(teilchen(:)%momentum(k))
    end do

!    write(*,*)' part1: ', teilchen(1)%id,  teilchen(1)%charge
!    write(*,*) teilchen(1)%momentum
!    write(*,*)' part2: ', teilchen(2)%id,  teilchen(2)%charge
!    write(*,*) teilchen(2)%momentum

    betacm = lorentzCalcBeta (momentum_Total, 'sqrtS_Free')

    mom(0:3)=teilchen(1)%momentum(0:3)
    call lorentz(betacm,mom,'twoBodyTools')

    ! momentum(1:3) is CM-Momentum
    ! In CM-frame without potentials: s=(E_1+E_2)**2 with E=SQRT[m**2+pCM**2]. 
    ! Therefore:
    mom2 = Dot_Product(mom(1:3),mom(1:3))
    sqrtS_free=sqrt(teilchen(1)%mass**2+mom2) &
         &    +sqrt(teilchen(2)%mass**2+mom2)
  end function sqrtS_Free


  !*************************************************************************
  !****s* twoBodyTools/ranCharge
  ! NAME
  ! subroutine ranCharge(izmin,izmax,iztot,izout,flag)
  ! PURPOSE 
  ! This subroutine distributes the charges randomly to a given total charge iztot
  ! INPUTS
  ! * integer, intent(in)   :: izmin(:) ! Vector of minimal charges
  ! * integer, intent(in)   :: izmax(:) ! Vector of maximal charges
  ! * integer, intent(in)   :: iztot    ! total charge
  ! OUTPUT
  ! * real, dimension(:)    :: izout    ! Random charge configuration
  ! * logical, intent(out)  :: flag     ! =.false. if procedure failed
  !*************************************************************************
  subroutine rancharge(izmin,izmax,iztot,izout,flag)
    use random
    integer, intent(in)  :: izmin(:),izmax(:), iztot
    integer, intent(out) :: izout(:)
    logical, intent(out) :: flag

    integer :: n,totmax,totmin,numc,numf,i,i2,j,d,iztv
    integer, parameter :: numfm = 100
    integer :: ifo(numfm)

    ! Check Input
    If((size(izout,dim=1).ne.size(izmax,dim=1)).or.(size(izout,dim=1).ne.size(izmin,dim=1))) then
       write(*,*) 'Critical error in ranCharge, dimensions do not fit'
    else
       n=(size(izout,dim=1))
    end if

    totmax=0
    totmin=0
    numc=1

    do i=1,n
       totmax=totmax+izmax(i)
       totmin=totmin+izmin(i)
       numc=numc*(izmax(i)-izmin(i)+1)
    end do

    if(totmin.gt.iztot.or.totmax.lt.iztot) then
       flag=.false.
       return
    else
       flag=.true.
    end if

    numf=0
    do i=0,numc-1
       i2=i
       iztv=0
       d=numc
       do j=1,n
          d=d/(izmax(j)-izmin(j)+1)
          izout(j)=izmin(j)+i2/d
          i2=i2-(i2/d)*d
          iztv=iztv+izout(j)
       end do
       !   *         write(*,*)i,'charges:',izout,numf
       if(iztv.eq.iztot) then
          numf=numf+1
          if(numf.gt.numfm) then
             write(*,*)'problems in rancharge numf.gt.numfm'
             numf=numfm
          end if
          ifo(numf)=i
       end if
    end do
    ! *monte-carlo decision
    if(numf.eq.0) then
       write(*,*)'problems in rancharge numf=0'
       flag=.false.
       return
    end if
    i=ifo(int(rn()*numf)+1)
    i2=i
    d=numc
    do j=1,n
       d=d/(izmax(j)-izmin(j)+1)
       izout(j)=izmin(j)+i2/d
       i2=i2-(i2/d)*d
    end do
    return
  end subroutine rancharge



  !*************************************************************************
  !****s* twoBodyTools/searchInInput
  ! NAME
  ! subroutine searchInInput(teilchenIn,id1,id2,particle_1,particle_2,failFlag)
  ! PURPOSE
  ! Given two IDs and a particleVector of dimension 2, this routine tries to find
  ! the given IDs in the particleVector. 
  ! INPUTS
  ! * type(particle),dimension(1,2), intent(in) :: teilchenIn
  ! * integer, intent(in) :: id1,id2
  ! OUTPUT
  ! * type(particle), intent(out) :: particle_1 ! copy of particle which has id1
  ! * type(particle), intent(out) :: particle_2 ! copy of particle which has id2
  ! * logical, intent(out) :: failFlag ! .true. if id1 or id2 could not be found in the particleVector
  !*************************************************************************
  subroutine searchInInput(teilchenIn,id1,id2,particle_1,particle_2,failFlag)
    use particleDefinition, only : particle
    type(particle),dimension(1:2), intent(in) :: teilchenIn
    integer, intent(in) :: id1,id2
    type(particle), intent(out) :: particle_1, particle_2
    logical, intent(out) :: failFlag

    If((teilchenIn(1)%ID.eq.id1).and.(teilchenIn(2)%ID.eq.id2)) then
       particle_1=teilchenIN(1)
       particle_2=teilchenIN(2)
       failflag=.false.
    else If((teilchenIn(1)%ID.eq.id2).and.(teilchenIn(2)%ID.eq.id1)) then
       particle_1=teilchenIN(2)
       particle_2=teilchenIN(1)
       failflag=.false.
    else 
       failflag=.true.
    end if

  end subroutine searchInInput



  !*************************************************************************
  ! cf. interface pCM
  !*************************************************************************
  real function pCM_1(srts, mass1, mass2)
    real, intent(in) :: srts, mass1, mass2
    real :: s,mass12
    s=srts**2
    mass12=mass1**2
    pCM_1 = sqrt(max(0.,(s+mass12-mass2**2)**2/(4.*s)-mass12))
  end function pCM_1
  !-------------------------------------------------------------------------
  real function pCM_2(srts, mass1, mass2, flagOK)
    real, intent(in) :: srts, mass1, mass2
    real :: s,mass12
    logical, intent(out) :: flagOK
    s=srts**2
    mass12=mass1**2
    pCM_2 = (s+mass12-mass2**2)**2/(4.*s)-mass12
    if (pCM_2.ge.0.) then
       flagOK = .true.
       pCM_2 = sqrt(pCM_2)
    else
       flagOK = .false.
       pCM_2 = 0.
    end if
  end function pCM_2
  !*************************************************************************


  !*************************************************************************
  !****f* twoBodyTools/pCM_sqr
  ! NAME
  ! real function pCM_sqr(s, msqr1, msqr2)
  ! PURPOSE
  ! This routine evaluates the CM-momentum of two particles, asssuming 
  ! vacuum kinematics.
  !
  ! Attention: Takes all arguments as squared quantities!!!
  !
  ! This routine does almost the same as 'pCM', but without taking the sqrt.
  ! Useful if a validity check is needed.
  !
  ! INPUTS
  ! * real :: msqr1, msqr2  ! squared masses of both particles
  ! * real :: s             ! Mandelstam s
  !
  ! OUTPUT
  ! * real :: pCM_sqr ! center of mass momentum squared
  !*************************************************************************
  real function pCM_sqr(s, msqr1, msqr2)
    real :: s,msqr1,msqr2
    pCM_sqr = (s+msqr1-msqr2)**2/(4.*s)-msqr1
  end function pCM_sqr



  !*************************************************************************
  !****f* twoBodyTools/p_lab
  ! NAME
  ! real function p_lab(srts, mass1, mass2)
  ! PURPOSE
  ! This routine evaluates the momentum of the 1-st particle in the rest 
  ! frame of the 2-nd particle, asssuming vacuum kinematics.
  !
  ! INPUTS
  ! * real :: mass1, mass2 ! masses of both particles
  ! * real :: srts ! SQRT(s)
  !
  ! OUTPUT
  ! * real :: p_lab ! laboratory momentum
  !*************************************************************************
  real function p_lab(srts, mass1, mass2)
    real :: srts,s,mass1,mass2,mass12
    s=srts**2
    mass12=mass1**2
    p_lab=sqrt(max(1.e-10,((s-mass12-mass2**2)/(2.*mass2))**2-mass12))
  end function p_lab



  !*************************************************************************
  !****f* twoBodyTools/MomentumTransfer2
  ! NAME
  ! real function MomentumTransfer2(part1,part2)
  ! PURPOSE
  ! This routine evaluates the momentum transfer squared, i.e. -t=-(p1-p2)^2
  !
  ! INPUTS
  ! * type(particle), intent(in) :: part1,part2   ! Incoming and outgoing particles, respectively
  ! OUTPUT
  ! * function value
  !*************************************************************************
  real function MomentumTransfer2(part1,part2)
    use particleDefinition
    type(particle), intent(in) :: part1,part2
    real :: momentum_transfer(0:3)
    momentum_transfer=part1%momentum-part2%momentum
    MomentumTransfer2=dot_product(momentum_transfer(1:3),momentum_transfer(1:3)) &
                    &-momentum_transfer(0)**2
  end function MomentumTransfer2


  !*************************************************************************
  !****s* twoBodyTools/convertToAntiParticles
  ! NAME
  ! subroutine convertToAntiParticles(a)
  !
  ! PURPOSE
  ! Converts a given preEvent to AntiParticles
  !
  ! INPUTS
  ! * type(preEvent), intent(inOut) :: a
  ! 
  ! OUTPUT
  ! * type(preEvent), intent(inOut) :: a
  !*************************************************************************
  subroutine convertToAntiParticles(a)
    use IdTable, only: isMeson, isBaryon, getAntiMeson
    use preEventDefinition, only : preEvent
    type(preEvent), intent(inOut), dimension(:) :: a
    integer :: newID, newCharge, i

    do i=lBound(a,dim=1),uBound(a,dim=1)
       
       If(a(i)%ID.ne.0) then
          if(isMeson(a(i)%ID)) then
             call getAntiMeson(a(i)%ID,a(i)%charge,newID,newCharge)
             a(i)%ID=newID
             a(i)%charge=newCharge
          else if(isBaryon(a(i)%ID)) then
             a(i)%charge=-a(i)%charge
             a(i)%antiparticle=.true.
          end if
       end if

    end do
  end subroutine convertToAntiParticles


  !*************************************************************************
  !****f* twoBodyTools/IsChargeExchange(part1,part2,part3,part4)
  ! NAME
  ! logical function IsChargeExchange(part1,part2,part3,part4)
  ! PURPOSE
  ! Compare incoming particles part1, part2 with outgoing particles
  ! part3, part4 and return .TRUE., if this is a charge exchange reaction.
  ! INPUTS
  ! * type(particle), intent(in) :: part1,part2,part3,part4
  ! OUTPUT
  ! * function value
  !*************************************************************************
  logical function IsChargeExchange(part1,part2,part3,part4)
    use particleDefinition
    type(particle), intent(in) :: part1,part2,part3,part4

    IsChargeExchange=.false.

    if( part1%Id.eq.part2%Id .and. (part1%antiParticle.eqv.part2%antiParticle) ) then

        if(  part3%Id.eq.part1%Id .and. part4%Id.eq.part1%Id .and. &
          &  (part3%antiParticle.eqv.part1%antiParticle) .and. &
          &  (part4%antiParticle.eqv.part1%antiParticle) .and. &
          & .not.(part3%charge.eq.part1%charge.or.part3%charge.eq.part2%charge) ) IsChargeExchange=.true.

    else

        if( part3%Id.eq.part1%Id .and. part4%Id.eq.part2%Id .and. &
          & (part3%antiParticle.eqv.part1%antiParticle) .and. &
          & (part4%antiParticle.eqv.part2%antiParticle) .and. &
          & part3%charge.ne.part1%charge )   then 
             IsChargeExchange=.true.
        else if( part4%Id.eq.part1%Id .and. part3%Id.eq.part2%Id .and. &
          &      (part4%antiParticle.eqv.part1%antiParticle) .and. &
          &      (part3%antiParticle.eqv.part2%antiParticle) .and. &
          &      part4%charge.ne.part1%charge ) then
             IsChargeExchange=.true.
        end if

    end if

  end function IsChargeExchange



end module twoBodyTools
