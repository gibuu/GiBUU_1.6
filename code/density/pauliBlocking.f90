!*****************************************************************************
!****m* /pauliBlockingModule
! NAME
! module PauliBlockingModule
! PURPOSE
! Contains all information and routines, which are necessary for the
! Pauli-Blocking of the test-particles
!*****************************************************************************
module PauliBlockingModule
  use histf90
  use CallStack, only: Traceback

  implicit none

  Private

  !***************************************************************************
  !****g* pauliBlockingModule/pauliSwitch
  ! SOURCE
  !
  integer, save :: pauliSwitch=1
  ! 
  ! PURPOSE
  ! * 0 : No Pauli blocking
  ! * 1 : dynamic Pauli blocking (use actual phase space densities)
  ! * 2 : analytic Pauli blocking (use ground state assumption) 
  !   (not possible for Heavy Ions!)
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/densDepMomCutFlag
  ! SOURCE
  !
  logical, save :: densDepMomCutFlag=.false.
  ! PURPOSE
  ! if .true. - the radius in momentum space for selecting nucleons around
  ! given nucleon will depend on local Fermi momentum
  ! NOTES
  ! Used only for dynamic pauli blocking. 
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/Gauss
  ! SOURCE
  !
  real, save :: Gauss=1.0
  ! PURPOSE
  ! Smearing for dynamic pauli blocking
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/cutGauss
  ! SOURCE
  !
  real, save :: cutGauss=2.2
  ! PURPOSE
  ! Cutoff for gauss Smearing
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/cutMom
  ! SOURCE
  !
  real, save :: cutMom=0.08
  ! PURPOSE
  ! * for densDepMomCutFlag=.false. --- 
  !   radius of phase space box in momentum space
  ! * for densDepMomCutFlag=.true. --- 
  !   minimum radius of phase space box in momentum space
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/cutPos
  ! SOURCE
  !
  real, save :: cutPos=1.86
  ! PURPOSE
  ! Radius of phase space box in position space
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/nGridPos
  ! SOURCE
  !
  integer, save :: nGridPos=30
  ! PURPOSE
  ! number of points in position space to save weights on
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/ensembleJump
  ! SOURCE
  !
  integer, save :: ensembleJump=5
  ! PURPOSE
  ! Parameter for speedup. Only every "ensemblejump"th ensemble is considered 
  ! to evaluate the probability for pauli blocking.
  !***************************************************************************

  !***************************************************************************
  !****g* pauliBlockingModule/DoHistogram
  ! SOURCE
  !
  logical, save :: DoHistogram=.false.
  ! PURPOSE
  ! if .true., a historgram is filled representing the blocking probability
  ! as function of the fermi momentum.
  ! You have to call 'WriteBlockMom' explicitely for writing the histogram
  !***************************************************************************

  real, save :: deltaPos   !difference in position space between different array elements in weights
  real, save, allocatable, dimension(:) :: weights  !weights for dynamic pauli blocking

  type(histogram),   save :: hBlockMom

  Public :: pauliBlocking, checkPauli, WriteBlockMom

contains

  !*****************************************************************************
  !****f* pauliBlockingModule/pauliBlocking
  ! NAME
  ! logical function pauliBlocking(momentum,position,Nukcharge,Teilchen,probabilityOut,weight)
  !
  ! PURPOSE
  ! * evaluates Pauli blocking for neutrons or protons
  ! * returns .true. if position and momentum are blocked by Pauli principle 
  ! * returns .false. if no Pauli blocking
  !
  ! INPUTS
  ! * real, dimension(0:3)  :: momentum
  ! * real, dimension(1:3)  :: position  
  ! * integer               :: Nukcharge -- Charge of nucleus
  ! * type(particle), dimension(:,:), OPTIONAL :: Teilchen -- 
  !   Full Particle Vector
  ! * real, OPTIONAL :: weight -- particle weight for documentation purposes
  !
  ! OUTPUT
  ! * real, OPTIONAL :: probabilityOut -- Probability of blocking, 
  !   i.e. occupation number
  !*****************************************************************************
  logical function pauliBlocking (momentum, position, Nukcharge, Teilchen, probabilityOut, weight)

    use particledefinition
    use constants, only: pi, hbarc

    real, dimension(0:3), intent(in)           :: momentum
    real, dimension(1:3), intent(in)           :: position  
    integer, intent(in)                        :: Nukcharge
    type(particle), dimension(:,:), intent(in), optional :: Teilchen
    real, intent(out), optional :: probabilityOut
    real, intent(in), optional :: weight

    real :: probability
    logical, save :: initSwitch=.true.

    !Check Input
    If ((Nukcharge/=1).and.(nukcharge/=0)) then
       Write(*,*) 'Wrong input in PauliBlocking. Nukcharge=', nukcharge
       call Traceback()
    end if

    !Read out jobcard
    if (initSwitch) then
       call init
       initSwitch=.false.
       if (DoHistogram) call CreateHist(hBlockMom, "pFermi", 0.000,0.300,0.005)
    end if

    !Evaluate Pauli Blocking
    select case(PauliSwitch)
    case(0) !Switched off
       pauliBlocking=.false.
    case(1) !Dynamic Pauli Blocking
       if(.not.Present(Teilchen)) then
          call TraceBack('Error: For dynamic pauliblocking the real particle vector must be input. stop!')
       end if
       call dynamicPauli(pauliBlocking,probability)
       if(present(probabilityOut)) probabilityOut=probability
    case(2) !Analytic Pauli Blocking
       call analyticPauli(pauliBlocking)
    case default
       write(*,*) 'This PauliSwitch is not valid:',PauliSwitch,'Choose 1 or 2!'
       call Traceback('Stop program')
    end select

  contains

    !***************************************************************************
    ! Analytic Pauli Blocking according to Local Density Approximation
    ! Makes only sense if the nucleus stays in groundstate
    ! Not for  heavy ions!!
    !***************************************************************************
    subroutine analyticPauli (pauliblocking)

      use dichteDefinition
      use densityModule, only: densityAt, FermiMomAt
      use lorentzTrafo, only: lorentz
      use minkowski, only: abs3

      real :: momAbs ! absolute momentum
      real :: fermiMomentum ! absolute momentum
      real :: beta(1:3) !beta value to boost to lrf
      logical,intent(out) :: pauliblocking
      real, dimension(0:3) :: momSave
      type(dichte) :: density

      density=densityAt(position)
      select case (nukcharge)
      case (0)
         beta=density%neutron(1:3)/density%neutron(0)
      case (1)
         beta=density%proton(1:3)/density%proton(0)
      end select
      fermiMomentum=FermiMomAt(position,nukcharge)
      momSave=momentum
      call lorentz(beta,momSave,'pauliBlocking') !boost momentum to lrf
      momAbs=abs3(momSave)

      pauliBlocking=(momAbs<fermiMomentum)

      if (DoHistogram) then
         if (present(weight)) then
            if (pauliBlocking) then
               call AddHist(hBlockMom,fermiMomentum, weight,weight)
            else
               call AddHist(hBlockMom,fermiMomentum, weight)
            end if
         end if
      end if

      !      write(*,'(A,4f11.3,L3)') 'Pauli :  ',position,fermiMomentum,pauliBlocking
      !      write(*,*) beta

    end subroutine analyticPauli

    !***************************************************************************
    subroutine dynamicPauli (pauliBlocking, probability)

      use dichteDefinition
      use densityModule, only: densityAt
      use lorentzTrafo, only: lorentz
      use random, only: rn
      use output, only: Write_InitStatus
      use IdTable, only: nucleon
      use minkowski, only: abs4

      logical, intent(out) :: pauliBlocking
      real, intent(out)    :: probability

      real :: dens   !density of considered nucleon species in local rest frame (=lrf)
      real :: momAbs ! absolute momentum
      real :: fermiMomentum ! Fermi momentum
      real :: distPosSQR, distMomSQR !Squares of momentum and position distances
      real :: MomentumCut
      real :: cutMomSqr,cutPosSqr !Squares of momentum and position cut off
      real :: beta(1:3) !beta value to boost to lrf
      real, dimension(0:3) :: momSave
      integer :: k,index,j
      logical, save :: initWeightsFlag=.true.
      type(dichte) :: density      

      if (initWeightsFlag) then
         call Write_InitStatus('Dynamic Pauli Blocking',0)
         Write(*,*) 'distance grid has grid spacing:', deltaPos
         write(*,*) 'Start numerical integrations for the weight factors'
         call initWeights
         call Write_InitStatus('Dynamic Pauli Blocking',1)
         initWeightsFlag=.false.
      end if

      if (Size(Teilchen,dim=1)<ensembleJump) ensembleJump=Size(Teilchen,dim=1)

      if (densDepMomCutFlag) then
         density=densityAt(position)
         if (nukcharge==0) then    !Evaluate invariant Neutron Density
            If (density%neutron(0)<1E-08) then
               pauliBlocking=.false.
               probability=0.
               return
            end if
            dens = abs4(density%neutron)
            beta = density%neutron(1:3)/density%neutron(0)
         else if (nukcharge==1) then !Evaluate invariant Proton Density
            If (density%proton(0)<1E-08) then
               pauliBlocking=.false.
               probability=0.
               return
            end if
            dens = abs4(density%proton)
            beta = density%proton(1:3)/density%proton(0)
         end if
         momSave=momentum
         call lorentz(beta,momSave,'pauliBlocking') !boost momentum to lrf
         fermiMomentum=(3.*pi**2*dens)**(1./3.)*hbarc
         momAbs=sqrt(momSave(1)**2+momSave(2)**2+momSave(3)**2)
         MomentumCut=max(cutMom,fermiMomentum-momAbs)
      else
         MomentumCut=cutMom
      end if

      cutMomSqr=MomentumCut**2
      cutPosSqr=(cutPos+cutGauss)**2
      probability=0.
      Do k=1,Size(Teilchen,dim=2) !loop over all particles in one ensemble
         Do j=1,Size(Teilchen,dim=1) !loop over all ensembles
            If (Mod(j,ensembleJump)/=0) cycle !only consider every 5th ensemble
            !Check wether particle is a real nucleon of the considered charge
            if (Teilchen(j,k)%perturbative.or.(Teilchen(j,k)%id/=nucleon) &
                .or.(Teilchen(j,k)%charge/=Nukcharge).or.(Teilchen(j,k)%antiParticle)) cycle
            !distance to particles in momentum space
            distMomSQR=(Teilchen(j,k)%momentum(1)-momentum(1))**2 &
                      +(Teilchen(j,k)%momentum(2)-momentum(2))**2 &
                      +(Teilchen(j,k)%momentum(3)-momentum(3))**2
            If (distMomSQR>cutMomSQR) cycle
            !distance to particles in position space
            distPosSQR=(Teilchen(j,k)%position(1)-position(1))**2 &
                      +(Teilchen(j,k)%position(2)-position(2))**2 &
                      +(Teilchen(j,k)%position(3)-position(3))**2
            If (distPosSQR>=cutPosSqr) cycle
            index=NINT(Sqrt(distPosSQR)/deltaPos)
            !linear interpolation of weights
            !probability=probability&
            !&+weights(index)+(weights(index+1)-weights(index))/deltaPos*(distPos-index*deltaPos)   
            probability=probability+weights(index)   
         End do
      End do
      probability=probability/float(size(Teilchen,dim=1)/ensembleJump) !divide by number ensembles
      if (densDepMomCutFlag) probability=probability*(cutMom/MomentumCut)**3
      !      write(*,*) ' Pos, MomCut, Prob:', &
      !                 & sqrt(dot_product(position(:),position(:))),MomentumCut,probability

      pauliblocking = (probability > rn())

    end subroutine dynamicPauli

    !***************************************************************************
    ! Evaluates Integral over box in position space to get the contribution
    ! of a particle which is rdist away.
    ! wheights(i) is array of those integrals with rdist=i*deltaPos
    !***************************************************************************
    subroutine initWeights
      use errorFunction, only: errorFunc

      integer :: k,m,n
      real :: rdist, integral 
      integer, parameter :: stepsR=2000        !integration steps for radial integration
      integer, parameter :: stepsTheta=2000    !integration steps for theta integration
      real :: distance !distance between vec(rdist) and vec(r)
      real :: Costheta,r,volumeBox,norm,deltaCosTheta,deltaR,deltaIntegral

      ! Volume of phase space box; factor two for spin:
      volumeBox=2.*(4./3.*pi*cutPos**3)*(4./3.*pi*cutMom**3)/((2.*pi)**3)/(hbarc**3)
      norm=errorFunc(cutGauss/sqrt(2.)/gauss)-cutGauss*sqrt(2./pi)/gauss*Exp(-cutGauss**2/2/gauss**2)
      deltaCosTheta=1./float(stepsTheta)
      deltaR=cutPos/float(StepsR)
      deltaPos=(cutPos+cutGauss)/float(nGridPos)
      allocate(weights(0:nGridPos))
      Do k=0,Size(weights)-1
         rdist=float(k)*deltaPos                  !rdist=0...gaussCut+positionCut
         Integral=0.
         Do n=-stepstheta,stepstheta-1            !Integral over Cos(theta)=-1,...1
            Costheta=(float(n)+0.5)*deltaCosTheta
            Do m=0,stepsR-1                       !integral over r=0...posCut
               r=(float(m)+0.5)*deltaR
               distance=SQRT(r**2+rdist**2-2.*r*rdist*Costheta)     !ABS(vec(r)-vec(rdist))
               If (distance.le.cutGauss) then
                  DeltaIntegral=r**2/exp((distance**2)/2./gauss**2)
                  Integral=Integral+DeltaIntegral
               end if
            end do
         End do
         Integral=Integral*2.*pi*deltaCosTheta*deltaR/(2.*pi*gauss**2)**(3./2.)
         weights(k)=integral/norm/volumeBox
      end do

      !Output
      Write(*,*) ' Print weights to pauliWeights.dat'
      open(15,file='pauliWeights.dat')
      write(15,*) '# Volume of phase space box=',volumeBox
      write(15,*) '# distance , weigth factor of nucleon at this distance'
      Do k=0,Size(weights)-1
         Write(15,*) deltaPos*float(k),weights(k)
      End do
      close(15)

    end subroutine initWeights

  end function pauliBlocking


  !***************************************************************************
  !****s* pauliBlockingModule/WriteBlockMom
  ! NAME
  ! subroutine WriteBlockMom(mul)
  ! PURPOSE
  ! Write Histogram 'PauliBlocking.BlockMom.dat'
  !***************************************************************************
  subroutine WriteBlockMom (mul)
    real, intent(in) :: mul
    if (DoHistogram) &
      call WriteHist (hBlockMom, mul=mul, add=1e-20, file='PauliBlocking.BlockMom.dat', dump=.true.)
  end subroutine WriteBlockMom


  !***************************************************************************
  !****s* pauliBlockingModule/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads data out of jobcard 'initPauli'
  !***************************************************************************
  subroutine init
    use output, only: Write_ReadingInput
    use nucleusDefinition
    use nucleus, only: getTarget
    use inputGeneral, only: eventType
    use eventtypes, only: elementary

    !***************************************************************************
    !****n* pauliBlockingModule/initPauli
    ! NAME
    ! NAMELIST initPauli
    ! PURPOSE
    ! Includes the input switches and variables:
    ! * pauliSwitch
    ! * densDepMomCutFlag
    ! * Gauss
    ! * cutGauss
    ! * cutMom
    ! * cutPos
    ! * nGridPos 
    ! * ensembleJump
    ! * DoHistogram
    !***************************************************************************
    
    NAMELIST /initPauli/ pauliSwitch, densDepMomCutFlag, Gauss, cutGauss, &
                         cutMom, cutPos, nGridPos, ensembleJump, DoHistogram

    character(len=23), dimension(0:2), parameter :: NN = (/ 'No Pauli blocking      ', &
                                                            'Dynamic Pauli blocking ', &
                                                            'Analytic Pauli blocking' /)
    integer :: ios
    type(tNucleus), pointer :: targetNuc

    call Write_ReadingInput('initPauli',0)
    rewind(5)
    read(5,nml=initPauli,iostat=ios)
    call Write_ReadingInput('initPauli',0,ios)

    if (pauliSwitch<0 .or. pauliSwitch>2) then
       write(*,'(A,i3," = ",A)') '  pauliSwitch        = ',pauliSwitch,' STOP !!!'
       call Traceback()
    end if

    write(*,'(A,i3," = ",A)') '  pauliSwitch        = ',pauliSwitch,NN(pauliSwitch)
    
    targetNuc => getTarget()
    if (eventType==elementary .or. targetNuc%mass==1) then
      pauliSwitch = 0
      write (*,*) 'pauliSwitch is set to 0 for elementary target'
    end if

    write(*,*) ' densDepMomCutFlag  = ',densDepMomCutFlag
    write(*,*) ' Gauss         [fm] = ',Gauss
    write(*,*) ' cutGauss      [fm] = ',cutGauss
    write(*,*) ' cutMom     [GeV/c] = ',cutMom
    write(*,*) ' cutPos        [fm] = ',cutPos
    write(*,*) ' nGridPos           = ',nGridPos
    write(*,*) ' ensembleJump       = ',ensembleJump
    write(*,*) ' DoHistogram        = ',DoHistogram

    call Write_ReadingInput('initPauli',1)

  end subroutine init


  !***************************************************************************
  !****s*  pauliBlockingModule/checkPauli
  ! NAME
  ! subroutine checkPauli(teilchen,realParticles,collisionFlag)
  !
  ! PURPOSE
  ! Administrates the Pauli-Blocking decision
  !
  ! INPUTS
  ! * type(particle),dimension(:)   :: teilchen      -- outgoing particles
  ! * type(particle),dimension(:,:) :: realParticles -- real particle vector
  ! 
  ! OUTPUT
  ! * logical ::  collisionFlag -- 
  !  .true. = momentum of particle "teilchen" IS NOT pauli blocked, 
  !  .false.= momentum of particle "teilchen" IS pauli blocked
  !
  ! NOTE
  ! The position and momentum of the outgoing particles must be set!
  !***************************************************************************
  subroutine checkPauli (teilchen, realParticles, collisionFlag)

    use IdTable, only : nucleon
    use particleDefinition

    type(particle), intent(in), dimension(:)   :: teilchen
    type(particle), intent(in), dimension(:,:) :: realParticles
    logical, intent(out)                       :: collisionFlag

    integer :: i

    collisionFlag=.true.

    ! Check for pauli-blocking of particles in final state
    do i=lbound(teilchen,dim=1),ubound(teilchen,dim=1)
       if (teilchen(i)%ID==nucleon .and. .not.(teilchen(i)%Antiparticle)) then
          if (pauliBlocking(teilchen(i)%momentum,teilchen(i)%position,teilchen(i)%charge, &
                            realParticles,weight=teilchen(i)%perweight)) then
             collisionFlag=.false.
             return
          end if
       end if
    End do

  end subroutine checkPauli


end module PauliBlockingModule
