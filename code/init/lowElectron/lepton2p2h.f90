!***************************************************************************
!****m* /lepton2p2h
! NAME
! module lepton2p2h
!
! PURPOSE
! Do all the internals for 2p2h scattering:
!
! EM:
! * ell N1 N2 --> ell' N1' N2'  == gamma* N1 N2 --> N1' N2'
! * ell N1 N2 --> ell' N Delta  == gamma* N1 N2 --> N Delta
! NC:
! * nu  N1 N2 --> nu'  N1' N2'
! * nu  N1 N2 --> nu'  N Delta
! CC:
! * nu  N1 N2 --> ell- N1' N2' (sum of hadronic charges increases by +1)
! * nu  N1 N2 --> ell- N Delta (  -- " --                              )
!
! antiEM, antiNC and antiCC are the same as EM, NC, CC.
!***************************************************************************
module lepton2p2h
  use particleDefinition
  use eN_eventDefinition
  use leptonicID
  use CALLSTACK, only: TRACEBACK

  implicit none

  PUBLIC :: lepton2p2h_DoQE
  PUBLIC :: lepton2p2h_DoDelta

  !*************************************************************************
  !****g* lepton2p2h/ME_Version
  ! PURPOSE
  ! indicate the type of matrix element parametrisation
  !
  ! SOURCE
  integer, save :: ME_Version = 8
  !
  ! possible values:
  ! * 1: Monopole parametrisation (cf. const*ME_Norm_XX/(1+ Q^2/ME_Mass_XX^2)^2 ) ! const for CC  fitted to MiniBooNE is 3.5e-6
  ! * 2: Monopole parameterization in transverse part only ! const for CC  fitted to MiniBooNE is 1e-4
  ! * 3: Monopole with suppressed low Q2  ( ME_Norm_XX*Q2/(1+ Q^2/ME_Mass_XX^2)^2 )
  ! * 4: const ME_Norm_XX  ! const for CC  fitted to MiniBooNE is 1.8e-6
  ! * 5: const with transverse part only
  ! * 6: fall with W
  ! * 7: as 4  and decreasing with Enu
  ! * 8: as 5  and decreasing with Enu
  ! * 9: exponential fall with Q2
  ! * 10: transverse,  exponential fall with Q2
  ! * 11:
  ! * 12: "Dipole transverse" transverse,  fall with Q2 as 4-th power
  ! * 13: transverse,  fall with Enu
  ! * 14: transverse,  fall with qz^2
  ! * 15: like Bosted arXiV:1203.2262  and transverse
  !*************************************************************************

  !*************************************************************************
  !****g* lepton2p2h/ME_Norm_QE
  ! PURPOSE
  ! Parametrisation of matrix element, parameter no. 1,
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Norm_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !*************************************************************************

  !*************************************************************************
  !****g* lepton2p2h/ME_Norm_Delta
  ! PURPOSE
  ! Parametrisation of matrix element, parameter no. 1,
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Norm_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !*************************************************************************

  !*************************************************************************
  !****g* lepton2p2h/ME_Mass_QE
  ! PURPOSE
  ! Parametrisation of matrix element, parameter no. 2
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Mass_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !*************************************************************************

  !*************************************************************************
  !****g* lepton2p2h/ME_Mass_Delta
  ! PURPOSE
  ! Parametrisation of matrix element, parameter no. 2
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Mass_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !*************************************************************************

  logical, save:: initflag = .true.

contains
  !*************************************************************************
  !****is* lepton2p2h/readInput
  ! NAME
  ! subroutine readInput
  !*************************************************************************
  subroutine readInput
    use output

    integer :: ios, i

    !***********************************************************************
    !****n* lepton2p2h/Lepton2p2h
    ! NAME
    ! NAMELIST /Lepton2p2h/
    ! PURPOSE
    ! Includes parameters for QE events:
    ! * ME_Version
    ! * ME_Norm_QE
    ! * ME_Norm_Delta
    ! * ME_Mass_QE
    ! * ME_Mass_Delta
    !***********************************************************************
    NAMELIST /lepton2p2h/ ME_Version, &
                          ME_Norm_QE, ME_Norm_Delta, &
                          ME_Mass_QE, ME_Mass_Delta

    if(.not.initFlag) return

    call Write_ReadingInput('lepton2p2h',0)
    rewind(5)
    read(5,nml=lepton2p2h,IOSTAT=ios)
    call Write_ReadingInput("lepton2p2h",0,ios)

    select case (ME_Version)
    case (1)
       write(*,'(A)') 'ME1 monopole  parametrization of the "form factor":'
       write(*,'(A)') '  [parametrisation of the ME squared: 7.0e-6 * ME_Norm_XX * (1+Q^2/M^2)^(-2) ]'
    case (2)
       write(*,'(A)') 'ME2  monopole parameterization in transverse part only'
    case (3)
       write(*,'(A)') 'ME3  monopole parameterization with suppressed low Q2 ( ME_Norm_XX*Q2/(1+ Q^2/ME_Mass_XX^2)^2 )'
    case (4)
       write(*,'(A)') 'ME4  const =  4.0e-6 * ME_Norm_XX'
    case (5)
       write(*,'(A)') 'ME5  const =  2e4 * ME_Norm_XX in transverse part only'
    case (6)
       write(*,'(A)') 'ME6   decreasing with W '
    case (7)
       write(*,'(A)') 'ME7   2.0e-6 * ME_Norm_XX  decreasing with Enu'
    case (8)
       write(*,'(A)') 'ME8  4.8e4 * 0.635 / Enu^2 * ME_Norm_XX in transverse part only   decreasing with Enu'
    case (9)
       write(*,'(A)') 'ME9  exponential fall with Q2'
    case (10)
       write(*,'(A)') 'ME10  transverse part, exponential fall with Q2'
    case (11)
       write(*,'(A)') 'ME11  dipole parameterization of the "form factor" '
    case (12)
       write(*,'(A)') 'ME12  dipole parameterization of the "form factor" ME=8e4*(1+Q2/MA2)^{-4},transverse part only'
    case (13)
       write(*,'(A)') 'ME13  4.8e4 / Enu  * ME_Norm_XX in transverse part only   decreasing with Enu'
       write(*,'(A)') '               (similar to ME8 but to fir versus reconstructed)'
    case (14)
       write(*,'(A)') 'ME14  fall down with qz'
    case (15)
       write(*,'(A)') 'parameterization like Bosted arXiV:1203.2262  and transverse'
    case DEFAULT
       write(*,*) 'ME_Version = ',ME_Version
       call TRACEBACK('wrong value for ME_Version')
    end select

   write(*,'(A)') 'parameters ( N N final state ):  [i=EM,CC,NC]'
   do i=1,3
      write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
           & ME_Norm_QE(i),ME_Mass_QE(i)
   end do
   write(*,'(A)') 'parameters ( N Delta final state ):'
   do i=1,3
      write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
           & ME_Norm_Delta(i),ME_Mass_Delta(i)
   end do

    call Write_ReadingInput('lepton2p2h',1)

    initFlag=.false.

  end subroutine readInput

  !*************************************************************************
  !****s* lepton2p2h/lepton2p2h_SelectN2
  ! NAME
  ! subroutine lepton2p2h_SelectN2(eN)
  !
  ! PURPOSE
  ! Finds the second nucleon for the 2p2h collision
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(electronNucleon_event) :: eN -- a second nucleon is added
  !
  ! NOTES
  ! * The seond particle is generated analytically, not by selecting
  !   a testparticle from the real particle vector.
  ! * This is at a very basic level. You may add more sophisticated features
  !   as eq. two-particle correlatione etc.
  ! * A threshold check  Wfree>(2*mN+1MeV) is performed
  !*************************************************************************
  subroutine lepton2p2h_SelectN2(eN,flagOK)
    use mediumDefinition
    use mediumModule, only: mediumAt
    use densitymodule, only : FermiMomAt
    use random, only : rn, rnOmega
    use constants, only : mN
    use energyCalc, only : energyDetermination
    use minkowski, only : abs4
    use lorentzTrafo, only : lorentz

    type(electronNucleon_event), intent(inout) :: eN
    logical, intent(out) :: flagOK

    type(medium)   :: media
    type(particle) :: partN2
    real :: p,pF
    type(particle), dimension(2) :: nucleon
    real, dimension(0:3) :: momentum
    integer :: i

    ! 0) Set some defaults:
    call setToDefault(partN2)
    flagOK = .false.

    partN2%ID = 1
    partN2%mass=mN

    ! 1) select charge:
    media=mediumAt(eN%nucleon%position)
    if (rn()*media%density.gt.media%densityProton) then
       partN2%charge = 0
    else
       partN2%charge = 1
    end if

    ! 2) select position:
    partN2%position = eN%nucleon%position

    ! 3) select 3-momentum:
    pF = FermiMomAt(partN2%position)
    p = pF * rn()**(1./3.)
    partN2%momentum(1:3) = p * rnOmega()
    partN2%momentum(0) = sqrt(mN**2+p**2)

    call energyDetermination(partN2)

    ! 4) change the eN information:

    eN%nucleon2 = partN2

    nucleon(1) = eN%nucleon
    nucleon(2) = eN%nucleon2
    momentum = eN%boson%momentum+nucleon(1)%momentum+nucleon(2)%momentum
    eN%betacm = momentum(1:3)/momentum(0)
    eN%W = abs4(momentum)

    ! we calculate Wfree in the CM system:
    do i=1,2
       nucleon(i)%position = 999999999
       call lorentz(eN%betacm,nucleon(i)%momentum)
       nucleon(i)%momentum(0) = FreeEnergy(nucleon(i))
       call lorentz(-eN%betacm,nucleon(i)%momentum)
    end do

    momentum = eN%boson%momentum+nucleon(1)%momentum+nucleon(2)%momentum
    eN%W_free = abs4(momentum)

    if (eN%W_free.le.2*mN+0.001) return ! ===> failure

    ! 5) abuse of 'offshellParameter' for storage of density,
    !    needed for the 'cross section' calculation

    eN%nucleon2%offshellParameter = media%density
    flagOK = .true.

  end subroutine lepton2p2h_SelectN2

  !*************************************************************************
  !****f* lepton2p2h/lepton2p2h_XS
  ! NAME
  ! real function lepton2p2h_XS(eN,outPart,DoQE)
  ! PURPOSE
  ! calculate the electron induced 2p2h-QE cross section
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  ! * type(particle),dimension(:) :: OutPart -- the outgoing particles
  ! * logical :: DoQE -- .true. for NN final state, .false. for N Delta
  ! OUTPUT
  ! * the function value
  ! NOTES
  ! * One has to give a realistic parametrization of the matrix element
  ! * If one randomly selects the position of the second particle, one
  !   has to respect this in the XS calculation (and maybe not to
  !   multiply it with the density at the position)
  !*************************************************************************
  real function lepton2p2h_XS(eN,outPart,DoQE)
    use minkowski, only : SP,abs4,abs4Sq
    use constants, only : twopi,mN,hbarc
    use twoBodyTools, only: pCM_sqr

    type(electronNucleon_event), intent(in) :: eN
    type(particle),dimension(:), intent(in) :: OutPart
    logical,                     intent(in) :: DoQE

    real :: mf1_2,mf2_2, k1 ! absolute value of the 3-momentum of the outgoing lepton
    real :: A2,pcm2
    real :: ME

    mf1_2 = abs4Sq(outPart(1)%momentum)
    mf2_2 = abs4Sq(outPart(2)%momentum)

    ! in analogy to dSigmadcosTheta_l_dE_l_BW_eN or dSdO_fdE_fdO_k_med_eN:

    k1=max(  0.0, sqrt( Dot_Product(en%lepton_out%momentum(1:3),en%lepton_out%momentum(1:3)) )  )

    lepton2p2h_XS=k1/(&
         twopi**5*8.*SP(eN%lepton_in%momentum,eN%nucleon%momentum))

    A2 = abs4Sq(eN%boson%momentum+eN%nucleon%momentum+eN%nucleon2%momentum)
    pcm2=pCM_sqr(A2, mf1_2, mf2_2)

    lepton2p2h_XS=lepton2p2h_XS &
         * sqrt(pcm2/(16.0*A2)) & ! <-- the deltas
         * 2*twopi              & ! <-- the angular integration
         * eN%nucleon2%offshellParameter ! <-- abuse !! = media%density

    ! Now we have to calculate the Matrixelement:
    select case (ME_Version)
    case (1)
       ME = ME_Monopole(eN)
    case (2)
       ME = ME_Monopole_transverse(eN)
    case (3)
       ME = ME_Monopole(eN)*eN%Qsquared
    case (4)
       ME = ME_const(eN)
    case (5)
       ME = ME_transverse(eN)
    case (6)
       ME = ME_fallWithW(eN)
    case (7)
       ME = ME_const(eN)*min(1.0, 5.0/( eN%lepton_in%momentum(0) + 1.43 )**2 )
    case (8)
       ME = ME_transverse(eN)*0.635/( eN%lepton_in%momentum(0) )**2 !!!*1.17/eN%lepton_in%momentum(0)/(0.42+eN%lepton_in%momentum(0))**2
    case (9)
       ME = ME_const(eN)*1.4*exp(-ME_Mass_QE(abs(eN%idProcess))*eN%QSquared)
    case (10)
       ME = ME_transverse(eN)*exp(-ME_Mass_QE(abs(eN%idProcess))*eN%QSquared)
    case (11)
       ME = ME_Dipole(eN)
    case (12)
       ME = ME_Dipole_transverse(eN)
    case (13)
       ME = ME_transverse(eN)*0.25 / eN%lepton_in%momentum(0)
    case (14)
       ME = ME_von_qz2_transverse(eN)
    case (15)
       ME = ME_Bosted_transverse(eN)
    end select

    lepton2p2h_XS=lepton2p2h_XS* ME / 2./eN%nucleon2%momentum(0)

    ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2= (197/1000)**2 * 10 mb
    ! Now the cross section is given in units of mb/GeV:
    lepton2p2h_XS = lepton2p2h_XS*hbarc**2*10.

    ! Symmetry-Factor:
    if (IsSamePart(OutPart(1),OutPart(2))) lepton2p2h_XS = lepton2p2h_XS *0.5
    if (IsSamePart(eN%nucleon,eN%nucleon2)) lepton2p2h_XS = lepton2p2h_XS *0.5

  contains

    !***********************************************************************
    !****f* lepton2p2h_XS/ME_Monopole
    ! NAME
    ! real function ME_Monopole(eN)
    ! PURPOSE
    ! calculate the matrix element according the monopole parametrisation of the "form factor"
    ! so that the ME falls down as dipole
    !***********************************************************************
    real function ME_Monopole(eN)
      type(electronNucleon_event), intent(in) :: eN

      integer :: iP
      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===

         ME_Monopole=7.0e-6*ME_Norm_QE(iP)*(1.+eN%QSquared/ME_Mass_QE(iP)**2)**(-2)

      else           !=== N Delta final state ===

         ME_Monopole=ME_Norm_Delta(iP)*(1.+eN%QSquared/ME_Mass_Delta(iP)**2)**(-2)

      end if

    end function ME_Monopole

    !***********************************************************************
    !****f* lepton2p2h_XS/ME_Monopole_transverse
    ! NAME
    ! real function ME_dipole_transverse(eN)
    ! PURPOSE
    ! calculate the 2p2h matrix element according to  W_1(g_munu -q_um q_nu /Q2) * L^munu
    ! so that the contribution is only to the transverse part
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !***********************************************************************
    real function ME_Monopole_transverse(eN)
      use minkowski, only : metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME=3.4e4*ME_Norm_QE(iP)*(1+eN%QSquared/ME_Mass_QE(iP)**2)**(-2)/eN%lepton_in%momentum(0)
      else           !=== N Delta final state ===
         ME=3.4e4*ME_Norm_Delta(iP)*(1+eN%QSquared/ME_Mass_Delta(iP)**2)**(-2)/eN%lepton_in%momentum(0)
      end if


      do mu=0,3
        do nu=0,3
           dummy(mu,nu)=eN%boson%momentum(mu)*eN%boson%momentum(nu)/eN%QSquared
        end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%momentum,eN%lepton_out%momentum)

      ME_Monopole_transverse = Contract(hadronTens,leptonTens)

    end function ME_Monopole_transverse




    real function ME_const(eN)

      type(electronNucleon_event), intent(in) :: eN
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME_const=1.0e-5*ME_Norm_QE(iP)
      else           !=== N Delta final state ===
         ME_const=1.0e-5*ME_Norm_Delta(iP)
      end if

    end function ME_const






    real function ME_transverse(eN)
      use minkowski, only : metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME=4.8e4*ME_Norm_QE(iP)
      else           !=== N Delta final state ===
         ME=4.8e4*ME_Norm_Delta(iP)
      end if

      do mu=0,3
        do nu=0,3
           dummy(mu,nu)=eN%boson%momentum(mu)*eN%boson%momentum(nu)/eN%QSquared
        end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%momentum,eN%lepton_out%momentum)

      ME_transverse = Contract(hadronTens,leptonTens)

    end function ME_transverse



    real function ME_fallWithW(eN)
      use minkowski, only : metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      real :: ME, W2 !, sigma

      W2=abs4Sq(eN%boson%momentum+(eN%nucleon%momentum+eN%nucleon2%momentum)/2.)
!       sigma=ME_Mass_QE(abs(eN%idProcess)) ! width of the gaussian distribution
      ME=0.

      if (DoQE) then !=== N N final state ===
        if (W2>0)   ME=ME_Norm_QE(abs(eN%idProcess)) * 1.0e-5 * exp(-ME_Mass_QE(abs(eN%idProcess))*W2)
        !! gauss distribution in W with center mN
        !!if (W2>0) ME=ME_Norm_QE(abs(eN%idProcess)) * 1.0e-5 * exp( - (sqrt(W2)-mN)**2 / 2./ sigma**2 ) / sigma/sqrt(twopi)
      else           !=== N Delta final state ===
        if (W2>0) ME=ME*exp(-ME_Mass_Delta(abs(eN%idProcess))*W2**2)
      end if

      ME_fallWithW = ME

    end function ME_fallWithW


    !***********************************************************************
    !****f* lepton2p2h_XS/ME_Dipole
    ! NAME
    ! real function ME_Dipole(eN)
    ! PURPOSE
    ! calculate the 2p2h matrix element according the dipole parametrisation
    !***********************************************************************
    real function ME_Dipole(eN)
      type(electronNucleon_event), intent(in) :: eN

      integer :: iP
      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===

         ME_Dipole=1.2e-5*ME_Norm_QE(iP)*(1.+eN%QSquared/ME_Mass_QE(iP)**2)**(-4)

      else           !=== N Delta final state ===

         ME_Dipole=ME_Norm_Delta(iP)*(1.+eN%QSquared/ME_Mass_Delta(iP)**2)**(-4)

      end if

    end function ME_Dipole

    !***********************************************************************
    !****f* lepton2p2h_XS/ME_Dipole_transverse
    ! NAME
    ! real function ME_Dipole_transverse(eN)
    ! PURPOSE
    ! calculate the 2p2h matrix element according to  W_1(g_munu -q_um q_nu /Q2) * L^munu
    ! so that the contribution is only to the transverse part
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !***********************************************************************
    real function ME_Dipole_transverse(eN)
      use minkowski, only : metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME=8.0e4*ME_Norm_QE(iP)*(1. + eN%QSquared/ME_Mass_QE(iP)**2)**(-4)
      else           !=== N Delta final state ===
         ME=8.0e4*ME_Norm_Delta(iP)*(1+eN%QSquared/ME_Mass_Delta(iP)**2)**(-4)
      end if


      do mu=0,3
        do nu=0,3
           dummy(mu,nu)=eN%boson%momentum(mu)*eN%boson%momentum(nu)/eN%QSquared
        end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%momentum,eN%lepton_out%momentum)

      ME_Dipole_transverse = Contract(hadronTens,leptonTens)

    end function ME_Dipole_transverse


    !***********************************************************************
    !****f* lepton2p2h_XS/ME_von_qz2_transverse
    ! NAME
    ! real function ME_von_qz2_transverse(eN)
    ! PURPOSE
    ! calculate the 2p2h matrix element according to  W_1(g_munu -q_um q_nu /Q2) * L^munu
    ! so that the contribution is only to the transverse part
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !***********************************************************************
    real function ME_von_qz2_transverse(eN)
      use minkowski, only : metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME, qz2
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      qz2=eN%QSquared+eN%boson%momentum(0)**2

      if (DoQE) then !=== N N final state ===
         ME=1.0e7*ME_Norm_QE(iP)/(1.+(qz2/ME_Mass_QE(iP)**2))
      else           !=== N Delta final state ===
         ME=1.0e7*ME_Norm_Delta(iP)/(1.+(qz2/ME_Mass_Delta(iP)**2))
      end if


      do mu=0,3
        do nu=0,3
           dummy(mu,nu)=eN%boson%momentum(mu)*eN%boson%momentum(nu)/eN%QSquared
        end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%momentum,eN%lepton_out%momentum)

      ME_von_qz2_transverse = Contract(hadronTens,leptonTens)

    end function ME_von_qz2_transverse




    !***********************************************************************
    !****f* lepton2p2h_XS/ME_Bosted_transverse
    ! NAME
    ! real function ME_Bosted_transverse(eN)
    ! PURPOSE
    ! calculate the 2p2h matrix element according to  W_1(g_munu -q_um q_nu /Q2) * L^munu
    ! so that the contribution is only to the transverse part
    ! W1 is taken from Bosted, Mamyan arXiv:1303.2262
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !***********************************************************************
    real function ME_Bosted_transverse(eN)
      use minkowski, only : metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      real :: p0,p1,p2,p3,p4,p5,p18,p19,f, Q2,omega
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC
      p0=0.005138
      p1=0.980710
      p2=0.046379
      p3=1.643300
      p4=6.982600
      p5=-0.226550
      p18=215.0
      p19=-0.045536

      Q2=eN%QSquared
      omega=eN%boson%momentum(0)

      f=( 1. + max(0.3,Q2)/p3 )**p4 * omega**p5 * (1. + p18*12.0**(1. + p19*Q2/2./mN/omega))
      if (DoQE) then !=== N N final state ===
         ME=ME_Norm_QE(iP)*1.0e20*p0*exp( -(eN%W-p1)**2/p2 )/f
         !write (*,*) '2p2h: ME15=',ME, '     f=',f
      else           !=== N Delta final state ===
         ME=ME_Norm_Delta(iP)*1.0e20*p0*exp( -(eN%W-p1)**2/p2 )/f
      end if


      do mu=0,3
        do nu=0,3
           dummy(mu,nu)=eN%boson%momentum(mu)*eN%boson%momentum(nu)/eN%QSquared
        end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%momentum,eN%lepton_out%momentum)

      ME_Bosted_transverse = Contract(hadronTens,leptonTens)

    end function ME_Bosted_transverse


  end function lepton2p2h_XS



  !*************************************************************************
  !****s* lepton2p2h/lepton2p2h_FinalState
  ! NAME
  ! subroutine lepton2p2h_FinalState(eN,outPart,DoQE,flagOK)
  ! PURPOSE
  ! Generate the final state of the electron 2p2h event
  !*************************************************************************
  subroutine lepton2p2h_FinalState(eN,outPart,DoQE,flagOK)
    use minkowski, only : abs4
    use mediumDefinition
    use mediumModule, only: mediumAt
    use collisionNumbering, only: pert_numbering
    use master_2Body, only: setKinematics
    use propagation, only: updateVelocity
    use IDtable, only: nucleon, delta
    use random, only: rn
    use particleProperties, only: hadron
    use baryonWidthMedium, only: get_MediumSwitch_coll

    type(electronNucleon_event), intent(in)    :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    logical,                     intent(in)    :: DoQE
    logical,                     intent(out)   :: flagOK

    type(particle), dimension(2) :: pairIN
    type(medium)         :: media
    real, dimension(1:3) :: betaToLRF
!     real, dimension(0:3) :: momentum
    integer :: i, ChargeIn

    flagOK = .false.

    if (size(OutPart).lt.2) call TRACEBACK('OutPart array too small.')

    pairIN = (/ eN%boson, eN%nucleon /)
    media=mediumAt(eN%nucleon%position)

    ChargeIn = eN%nucleon%Charge+eN%nucleon2%Charge

    call setToDefault(OutPart)

    if (DoQE) then !=== N N final state ===
       OutPart%ID =     (/ nucleon          , nucleon /)

       select case (eN%idProcess)
       case DEFAULT ! == EM, NC, antiEM, antiNC
          OutPart%Charge = (/ eN%nucleon%Charge, eN%nucleon2%Charge /)

       case(2)      ! == CC
          select case(ChargeIn)
          case (0)
             OutPart%Charge = (/ 0, 1 /)
          case (1)
             OutPart%Charge = (/ 1, 1 /)
          case (2)
             return ! ==> failure
          case DEFAULT
             call TRACEBACK('ChargeIn not allowed')
          end select

       case(-2)     ! == antiCC
          select case(ChargeIn)
          case (0)
             return ! ==> failure
          case (1)
             OutPart%Charge = (/ 0, 0 /)
          case (2)
             OutPart%Charge = (/ 0, 1 /)
          case DEFAULT
             call TRACEBACK('ChargeIn not allowed')
          end select

       end select

    else           !=== N Delta final state ===
       OutPart%ID =     (/ nucleon          , delta /)

       select case (eN%idProcess)
       case DEFAULT ! == EM, NC, antiEM, antiNC
          OutPart(1)%Charge = nint(rn())
          OutPart(2)%Charge = eN%nucleon%Charge+eN%nucleon2%Charge - OutPart(1)%Charge
       case(2)      ! == CC
          select case(ChargeIn)
          case (0)
             OutPart(1)%Charge = nint(rn())
             OutPart(2)%Charge = 1 - OutPart(1)%Charge
          case (1)
             OutPart(1)%Charge = 1
             if (rn()<1./3.) OutPart(1)%Charge = 0  ! Delta++ n : Delta+ p = 1:2
             !based on counting diagrams, but better ideas are needed
             OutPart(2)%Charge = 2 - OutPart(1)%Charge
          case (2)
             OutPart%Charge = (/ 1, 2 /)
          case DEFAULT
             call TRACEBACK('ChargeIn not allowed')
          end select
!          call TRACEBACK('CC not yet implemented')
       case(-2)     ! == antiCC
          select case(ChargeIn)
          case (2)
             OutPart(1)%Charge = nint(rn())
             OutPart(2)%Charge = 1 - OutPart(1)%Charge
          case (1)
             OutPart(1)%Charge = 0
             if (rn()<3./4.) OutPart(1)%Charge = 1  ! Delta0 n : Delta- p = 4:3
             !based on counting diagrams, but better ideas are needed
             OutPart(2)%Charge =  - OutPart(1)%Charge
          case (0)
             OutPart%Charge = (/ 0, -1 /)
          case DEFAULT
             call TRACEBACK('ChargeIn not allowed')
          end select
!          call TRACEBACK('antiCC not yet implemented')
       end select

       ! The following is in order to avoid problems in massass:
       if (.not.get_MediumSwitch_coll()) then
          ! minimal value: 0.938 + Delta-MinMass + epsilon
          if (eN%W_free .lt. hadron(1)%mass+hadron(2)%minmass+0.005) then
             return ! ==> failure
          end if
       end if


    end if

    OutPart%antiparticle=.false.
    OutPart%perturbative=.true.

    OutPart%perWeight=0. ! perturbative weight = XS (here only dummy)

    do i=1,2
       OutPart(i)%position=eN%nucleon%position
       OutPart(i)%event=pert_numbering(eN%nucleon)
    end do

    betaToLRF=0.
!     momentum = eN%boson%momentum+eN%nucleon%momentum+eN%nucleon2%momentum

    call setKinematics (eN%W, eN%W_free, betaToLRF, eN%betacm, media, pairIn, OutPart(1:2), flagOK, .false.)

    if (.not.flagOK) return ! ==> failure

    call updateVelocity(OutPart)

  end subroutine lepton2p2h_FinalState

  !*************************************************************************
  !****s* lepton2p2h/lepton2p2h_DoQE
  ! NAME
  ! subroutine lepton2p2h_DoQE(eN,outPart,XS)
  !
  ! PURPOSE
  ! Do all the electron induced 2p2h-QE scattering gamma* N1 N2 -> N1' N2'
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: OutPart -- the two produced nucleons
  ! * real :: XS -- the cross section
  !*************************************************************************
  subroutine lepton2p2h_DoQE(eN,outPart,XS)
    use output, only: WriteParticle

    type(electronNucleon_event), intent(inout) :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    real, intent(out) :: XS

    logical :: flagOK

    if(initFlag) call readInput

    XS = 0.0

    call lepton2p2h_SelectN2(eN,flagOK)
    if (.not.flagOK)  return ! ==> failure

!    call write_electronNucleon_event(eN)

    call lepton2p2h_FinalState(eN,outPart,.true.,flagOK)

!!$    if (.not.flagOK) then
!!$       call write_electronNucleon_event(eN)
!!$       write(*,*) 'Failure'
!!$       stop
!!$    end if

    if (.not.flagOK) return ! ==> failure

    XS = lepton2p2h_XS(eN,outPart,.true.)
    outPart%perWeight=XS

!    call WriteParticle(6,1,outPart)
!    stop

  end subroutine lepton2p2h_DoQE

  !*************************************************************************
  !****s* lepton2p2h/lepton2p2h_DoDelta
  ! NAME
  ! subroutine lepton2p2h_DoDelta(eN,outPart,XS)
  !
  ! PURPOSE
  ! Do all the electron induced 2p2h-QE scattering gamma* N1 N2 -> N Delta
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: OutPart -- the two produced hadrons
  ! * real :: XS -- the cross section
  !*************************************************************************
  subroutine lepton2p2h_DoDelta(eN,outPart,XS)
    use output, only: WriteParticle

    type(electronNucleon_event), intent(inout) :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    real, intent(out) :: XS

    logical :: flagOK

    if(initFlag) call readInput

    XS = 0.0

    call lepton2p2h_SelectN2(eN,flagOK)
    if (.not.flagOK) return ! ==> failure

    call lepton2p2h_FinalState(eN,outPart,.false.,flagOK)
    if (.not.flagOK) return ! ==> failure

    XS = lepton2p2h_XS(eN,outPart,.false.)
    outPart%perWeight=XS

  end subroutine lepton2p2h_DoDelta


end module lepton2p2h
