!***************************************************************************
!****m* /resonanceCrossSections
! NAME
! module resonanceCrossSections
! PURPOSE
! Includes all rottines to evaluate "baryon meson -> R -> X" cross sections
! with resonances R in the intermediate state.
!***************************************************************************


module resonanceCrossSections
  implicit none
  PRIVATE

  public :: barMes2resonance,barMes_R_barMes,sigma_npi_n2pi_resonances

  logical, save :: initInput_Flag=.true.

  !*************************************************************************
  !****g* resonanceCrossSections/fullPropagator
  ! PURPOSE
  ! Includes also the real parts in the resonance propagator. In former works
  ! (i.e. in the old Efffenberger code this has
  ! been neglected
  ! It should be set to .true. only if mediumSwitch_coll=.true. in the
  ! namelist width_Baryon
  !
  ! SOURCE
  !
  logical, save :: fullPropagator=.false.
  !*************************************************************************

contains


  subroutine readInput
    use output, only :Write_ReadingInput

    integer :: ios

    !***********************************************************************
    !****n* resonanceCrossSections/ResonanceCrossSections
    ! NAME
    ! NAMELIST ResonanceCrossSections
    ! PURPOSE
    ! Includes parameters:
    ! * fullPropagator
    !***********************************************************************

    NAMELIST /ResonanceCrossSections/ fullPropagator

    call Write_ReadingInput("ResonanceCrossSections",0)
    rewind(5)
    read(5,nml=ResonanceCrossSections,IOSTAT=ios)
    call Write_ReadingInput("ResonanceCrossSections",0,ios)
    write(*,*) '  Use the full self energy in the propagator of the resonances?',fullPropagator
    write(*,'(A)') '   (It should be TRUE only if mediumSwitch_coll=.true. in the namelist width_Baryon)'

    call Write_ReadingInput("ResonanceCrossSections",1)
  end subroutine readInput


  !*************************************************************************
  !****f* resonanceCrossSections/barMes2resonance
  ! NAME
  ! function barMes2resonance (idMeson_ini, idBaryon_ini, chargeMeson_ini, chargeBaryon_ini, propagated,
  ! mediumAtCollision, momentumLRF, masses, mesonMass, baryonMass, position, perturbative, srts) result (sigma)
  ! PURPOSE
  ! Evaluates baryon meson -> Resonance cross sections by summing over the resonances final state.
  ! One gets baryon meson -> X by summing over the (baryon meson -> R-> X) channels for all resonances.
  ! INPUTS
  ! * integer, intent(in) :: idMeson_ini,idBaryon_ini
  ! * integer, intent(in) :: chargeMeson_ini,chargeBaryon_ini
  ! * logical, intent(in) :: propagated
  !   if .true., then onsider only propagated particles in final state, otherwise
  !   consider all particles in final state
  ! * type(medium), intent(in) :: mediumAtCollision      ! Medium information
  ! * real, intent(in),dimension(0:3) :: momentumLRF     ! Total  Momentum in LRF=Momentum of resonance in LRF
  ! * real, intent(in)     ::  mesonMass, baryonMass     ! masses of colliding particles.
  ! * real, intent(in), dimension(1:3) :: position       ! Position where resonance is produced
  ! * logical, intent(in)              :: perturbative   ! flag for the perturbative nature of the resonance
  ! * real, intent(in), optional          :: srts        ! sqrt(s) of the collision -- needed in RMF mode
  ! RESULT
  ! * real, intent(out),  dimension(Delta:nbar) :: sigma   --- Cross section for each intermediate baryon resonance in units of mB. The index denotes the resonance.
  ! * real, intent(out),  dimension(Delta:nbar) :: masses  --- Masses of the produced resonances; index denotes the resonance.
  !*************************************************************************
  function barMes2resonance (idMeson_ini, idBaryon_ini, chargeMeson_ini, chargeBaryon_ini, propagated, mediumAtCollision, &
                             momentumLRF, masses, mesonMass, baryonMass, position, perturbative, srts) result (sigma)

    use mediumDefinition, only: medium, vacuum
    use baryonWidthMedium, only: WidthBaryonMedium, partialWidthBaryonMedium
    use particleProperties, only: hadron, isNonExotic
    use ClebschGordan, only: clebschSquared
    use twoBodyTools, only : resonanceMass
    use constants, only: pi, GeVSquared_times_mb
    use IdTable, only: Delta, nbar, isMeson, isBaryon
    use RMF, only : getRMF_flag
    use selfenergy_baryons, only: get_realPart,selfEnergy_imag
    use minkowski, only: abs4

    integer, intent(in) :: idMeson_ini, idBaryon_ini
    integer, intent(in) :: chargeMeson_ini, chargeBaryon_ini
    logical, intent(in) :: propagated
    type(medium), intent(in) :: mediumAtCollision
    real, intent(in), dimension(0:3) :: momentumLRF
    real, intent(in)                 :: mesonMass, baryonMass
    real, intent(in), dimension(1:3) :: position
    logical, intent(in)              :: perturbative
    real, intent(in), optional       :: srts

    real, dimension(Delta:nbar) :: sigma
    real, dimension(Delta:nbar), intent(out) :: masses ! masses of resonances

    integer :: resID, chargeRes
    real :: momCm, gamma_In, gammaTot, spinFactor, isoFactor
    real :: i1,i2,i3,iz1,iz2 ! Isospins & z-components
    real :: absMom_LRF, fullMass, PI_REAL, PI_IMAG
    logical, parameter :: debug=.false.

    if (initInput_Flag) then
       call readInput()
       initInput_Flag=.false.
    end if

    ! Check Input
    If (.not.(isMeson(idMeson_ini).and.isBaryon(idBaryon_ini))) then
       Write(*,*) 'Error in in input of "resonanceXsection" '
       write(*,*) idMeson_ini,idBaryon_ini
       stop
    End if

    sigma=0.
    masses=0.

    absMom_LRF=SQRT(Dot_product(momentumLRF(1:3),momentumLRF(1:3)))

    if (fullPropagator) fullmass = abs4(momentumLRF)

    ! Loop over intermediate resonances R : m B->  R -> X
    Do resId=Delta,nbar
       If (.not.hadron(resId)%usedForXsections) cycle ! Exclude resonances

       If (propagated.and.(.not.(hadron(resId)%propagated))) cycle

       ! Formula(2.52) Effenberger Phd.
       If (abs(hadron(idMeson_ini)%strangeness+hadron(idBaryon_ini)%strangeness-hadron(resID)%strangeness)>0.01) cycle ! No strangeness conservation
       If (abs(hadron(idMeson_ini)%charm+hadron(idBaryon_ini)%charm-hadron(resID)%charm)>0.01) cycle ! No charm conservation

       ! Evaluate isoSpin factors
       ! Convert Charges to isoSpinFactors in SU(4)
       ! Q=Y/2+I_3 with Y=S+C+B
       ! =>I_3=Q-(S+B+C)/2.
       ! for Initial State :
       i1 = real(hadron(idMeson_ini)%isoSpinTimes2)/2.
       iz1= real(chargeMeson_ini) - 0.5*real(hadron(idmeson_ini)%strangeness) - 0.5*real(hadron(idmeson_ini)%charm)
       i2 = real(hadron(idBaryon_ini)%isoSpinTimes2)/2.
       iz2= real(chargeBaryon_ini) - 0.5 - 0.5*real(hadron(idbaryon_ini)%strangeness) - 0.5*real(hadron(idbaryon_ini)%charm)
       i3 = real(hadron(resID)%isoSpinTimes2)/2.

       If (debug) then
         write(*,*)
         write (*,*) i1,i2,i3,iz1,iz2
       end if

       isoFactor=ClebschSquared(i1,i2,i3,iz1,iz2)

       ! Evaluate the spin factors
       spinFactor=(2.*hadron(resId)%spin+1.)/(2.*hadron(idBaryon_ini)%spin+1.)/(2.*hadron(idMeson_ini)%spin+1.)
       ! Evaluate the partial decay width for incoming and outgoing channel
       ! use no in-Width since in experiment I will also have no sharp particle

       ! Evaluate mass of the resonance and store it
       chargeRes=chargeMeson_ini+chargeBaryon_ini
       if (.not. getRMF_flag()) then
         masses(resID)=resonanceMass(resID,chargeRes,mediumAtcollision,momentumLRF,position,perturbative)
       else if( present(srts) ) then
         masses(resID)=srts
       else
         write(*,*)' In barMes2resonance: srts must present in RMF mode !!!'
         stop
       end if

       If (debug) write(*,*) masses(resID), resID
       If (masses(resID)<=0) then
          write(*,*) 'Mass of resonance less or equal to zero in resonanceCrossSections/barMes2resonance'
          write(*,*) 'Severe Error! stopping'
          write(*,*) idMeson_ini,idBaryon_ini,chargeMeson_ini, &
               chargeBaryon_ini,propagated,momentumLRF,sigma,masses
       end if

       !Evaluate CM-Momentum
       momCM=Sqrt(max((masses(resID)**2-(mesonMass+baryonMass)**2)*(masses(resID)**2-(mesonMass-baryonMass)**2) &
                      /4./masses(resId)**2,1E-8))
       If (debug) write(*,*) 'momCM=',momCM,'mass=',masses(ResID),'mesmass=', mesonMass,'barmass=', baryonMass

       ! use inWidth since there are both masses of the incoming particles given
       ! And don't use medium since only the out-width shall be dressed due to the final states
       gamma_In = partialWidthBaryonMedium (resID,masses(resID),.true.,IDmeson_ini, &
                                            IDbaryon_ini,momentumLRF,vacuum,baryonMass,mesonMass)
       gammaTot = WidthBaryonMedium (resID,masses(resID),momentumLRF,mediumATcollision)

       If ((Abs(gamma_In)<1E-6).or.(Abs(gammaTot)<1E-6)) then
          sigma(resID)=0.
          cycle
       end if

       If (resID==Delta.and.debug) write(300,'(6F8.4)')  momCM, Gamma_In,GammaTot, masses(resID)

       if (fullPropagator.and.isNonExotic(resID)) then
          PI_imag=selfenergy_imag(resID,absMom_LRF,momentumLRF(0),mediumAtCollision)
          PI_real=   get_RealPart(resID,absMom_LRF,fullmass      ,mediumAtCollision)
          sigma(resID) = isoFactor*spinFactor*4.*pi/momCm**2*masses(ResID)*Gamma_In*(-PI_imag) &
                         / ((fullMass**2-hadron(resID)%mass**2- PI_real)**2+PI_imag**2) / GeVSquared_times_mb
       else
          ! neglect real part: old Effenberger ansatz
          sigma(resID) = isoFactor*spinFactor*4.*pi/momCm**2*masses(ResID)**2*Gamma_In*GammaTot &
                         / ((masses(resID)**2-hadron(resID)%mass**2)**2+gammaTot**2*masses(resID)**2) / GeVSquared_times_mb
       end if

    End do
  end function barMes2resonance



  !*************************************************************************
  !****f* resonanceCrossSections/barMes_R_barMes
  ! NAME
  ! real function barMes_R_barMes (idMeson_ini, idBaryon_ini, idMeson_Final,
  ! idBaryon_Final, chargeMeson_ini, chargeBaryon_ini, chargeMeson_Final,
  ! chargeBaryon_Final, background, propagated, MediumAtCollision,
  ! momentumLRF, mesonMass_ini, baryonMass_ini, position, perturbative, srts)
  !
  ! PURPOSE
  ! Evaluates contribution of resonances to baryon meson -> baryon Meson
  ! One gets B m -> B m by summing over the (b M -> R->b M) channels for all
  ! resonances.
  !
  ! INPUTS
  ! * integer, intent(in) :: idMeson_ini,idBaryon_ini,idMeson_Final,idBaryon_Final --- ID's of particles in initial and final state.
  ! * integer, intent(in) :: chargeMeson_ini,chargeBaryon_ini,chargeMeson_Final,chargeBaryon_Final --- Charges of particles in initial and final state.
  ! * logical, intent(in) :: background ---
  !   .true. = Compute background cross section : Therefore consider only non-propagated particles in intermediate state;
  !   .false. = Full Xsection : Consider all particles in intermediate state
  ! * logical, intent(in) :: propagated ---
  !   .true. = Consider only propagated particles in intermediate state;
  !   .false. = Consider all particles in intermediate state
  !
  ! Information for In-Medium modifications :
  ! * type(medium), intent(in) :: mediumAtCollision      ! Medium information
  ! * real, intent(in),dimension(0:3) :: momentumLRF     ! Total  Momentum in LRF=Momentum of resonance in LRF
  ! * real, intent(in) :: mesonMass_ini, baryonMass_ini ! masses of colliding particles
  ! * real, intent(in),dimension(1:3) :: position        ! Position where resonance is produced
  ! * logical, intent(in)              :: perturbative   ! flag for the perturbative nature of the resonance
  ! * real, intent(in), optional                :: srts  ! sqrt(s) of the collision -- needed in RMF mode
  !

  ! NOTES
  ! Note, that it is not useful to set both background and propagated to .true. : Result=0.
  ! RESULT
  ! * resonance contribution to baryon meson -> baryon meson (in mb)
  !*************************************************************************

  real function barMes_R_barMes (idMeson_ini, idBaryon_ini, idMeson_Final, idBaryon_Final, &
                                 chargeMeson_ini, chargeBaryon_ini, chargeMeson_Final, chargeBaryon_Final, &
                                 background, propagated, MediumAtCollision, momentumLRF, &
                                 mesonMass_ini, baryonMass_ini, position,perturbative, srts)

    use baryonWidthMedium, only: WidthBaryonMedium, partialWidthBaryonMedium
    use mediumDefinition, only: medium, vacuum
    use twoBodyTools, only : resonanceMass
    use particleProperties, only: hadron, isNonExotic
    use ClebschGordan, only: clebschSquared
    use constants, only: pi, GeVSquared_times_mb
    use RMF, only : getRMF_flag
    use selfenergy_baryons, only: get_realPart, selfEnergy_imag
    use minkowski, only: abs4
    use IdTable, only: nbar, isBaryon, isMeson

    integer, intent(in)            :: idMeson_ini,idBaryon_ini,idMeson_Final,idBaryon_Final                       ! ID's of particles in initial and final state
    integer, intent(in)            :: chargeMeson_ini,chargeBaryon_ini,chargeMeson_Final,chargeBaryon_Final       ! Charges of particles in initial and final state

    logical, intent(in) :: background
    ! .true. = Compute background cross section : Therefore consider only non-propagated particles in intermediate state
    ! .false. = Full Xsection : Consider all particles in intermediate state

    logical, intent(in) :: propagated
    ! .true. = Consider only propagated particles in intermediate state
    ! .false. = Consider all particles in intermediate state

    type(medium), intent(in)          :: mediumATcollision      ! medium informations at the point of collision
    real, intent(in), dimension(0:3)  :: momentumLRF            ! total momentum in LRF
    real, intent(in)   :: mesonMass_ini, baryonMass_ini ! masses of colliding particles
    real, intent(in),dimension(1:3)  :: position                ! position of resonance
    logical, intent(in)              :: perturbative         ! flag for the perturbative nature of the resonance
    real, intent(in), optional       :: srts  ! sqrt(s) of the collision -- needed in RMF mode

    real :: sigma
    real :: momCM ! CM momentum
    integer :: resID ! Id of resonance
    real :: Gamma_In,Gamma_Out,Gamma_Tot
    real :: spinFactor,isoFactor
    real, dimension(1:3) :: i ! Isospins
    real, dimension(1:2) :: iz ! z-components of Isospins
    real :: massRes
    integer :: chargeRes
    logical,parameter ::  debug =.false.
    real :: absMom_LRF, fullMass,PI_REAL,PI_IMAG

    if(initInput_Flag) then
       call readInput()
       initInput_Flag=.false.
    end if

    ! Check Input
    If (.not.(isMeson(idMeson_ini).and.isMeson(idMeson_Final) &
         & .and.(isBaryon(idBaryon_ini)).and.(isBaryon(idBaryon_Final)))) then
       Write(*,*) 'Error in in input of "resonanceXsection" '
       write(*,*) idMeson_ini,idBaryon_ini,idMeson_Final, idBaryon_Final
       stop
    End if

    If(abs(chargeMeson_ini+chargeBaryon_ini-(chargeMeson_Final+chargeBaryon_Final)).gt.0.01) then
       ! No charge conservation
       barMes_R_barMes=0.
       return
    end if

    If(abs(hadron(idMeson_ini)%strangeness+hadron(idBaryon_ini)%strangeness &
         -(hadron(idMeson_Final)%strangeness+hadron(idBaryon_Final)%strangeness)).gt.0.01) then
       ! No strangeness conservation
       barMes_R_barMes=0.
       return
    end if
    If(abs(hadron(idMeson_ini)%charm+hadron(idBaryon_ini)%charm &
         -(hadron(idMeson_Final)%charm+hadron(idBaryon_Final)%charm)).gt.0.01) then
       ! No charm conservation
       barMes_R_barMes=0.
       return
    end if

    chargeRes=chargeMeson_ini+chargeBaryon_ini

    sigma=0.

    absMom_LRF=SQRT(Dot_product(momentumLRF(1:3),momentumLRF(1:3)))
    if(fullPropagator) then
       fullmass  =abs4(momentumLRF)
    end if

    ! Loop over intermediate resonances R : m B->  R -> m' B'
    Do resId=1,nbar
       if( .not.getRMF_flag() ) then
         massRes=resonanceMass(resID,chargeRes,mediumAtcollision,momentumLRF,position,perturbative)
       else if( present(srts) ) then
         massres=srts
       else
         write(*,*)' In barMes_R_barMes: srts must present in RMF mode !!!'
         stop
       end if
       If(debug) write(*,*) massres

       !Evaluate CM-Momentum

       momCM=Sqrt(max((massRes**2-(mesonMass_ini+baryonMass_ini)**2) &
            & *(massRes**2-(mesonMass_ini-baryonMass_ini)**2)/4./massRes**2,1E-8))

       If(.not.hadron(resId)%usedForXsections) cycle   ! Exclude resonances

       If(background.and.hadron(resId)%propagated) cycle ! Exclude propagated resonances for background

       If(propagated.and.(.not.(hadron(resId)%propagated))) cycle ! Exclude non-propagated resonances, use only the propagated ones

       ! Formula(2.52) Effenberger Phd.
       If(ABS(hadron(idMeson_ini)%strangeness+hadron(idBaryon_ini)%strangeness-hadron(resID)%strangeness).gt.0.01) then
          ! No strangeness conservation
          cycle
       end if
       If(ABS(hadron(idMeson_ini)%charm+hadron(idBaryon_ini)%charm-hadron(resID)%charm).gt.0.01) then
          ! No charm conservation
          cycle
       end if
       ! Evaluate isoSpin factors
       ! Convert Charges to isoSpinFactors in SU(4)
       ! Q=Y/2+I_3 with Y=S+C+B
       ! =>I_3=Q-(S+B+C)/2.
       ! (a) Initial State
       i(1) = float(hadron(idMeson_ini)%isoSpinTimes2)/2.
       iz(1)= float(chargeMeson_ini)-0.5*float(hadron(idmeson_ini)%strangeness)-0.5*float(hadron(idmeson_ini)%charm)
       i(2) = float(hadron(IDBaryon_ini)%isoSpinTimes2)/2.
       iz(2)= float(chargeBaryon_ini)-0.5 &
            & -0.5*float(hadron(idbaryon_ini)%strangeness)-0.5*float(hadron(idbaryon_ini)%charm)
       i(3)=  float(hadron(resID)%isoSpinTimes2)/2.
       If(debug) Write(*,*) 'In resonanceCrossSections'
       If(debug) Write(*,*) i(1),i(2),i(3),iz(1),iz(2)

       isoFactor=ClebschSquared(i(1),i(2),i(3),iz(1),iz(2))
       ! (b) Final State
       i(1) = float(hadron(idMeson_Final)%isoSpinTimes2)/2.
       iz(1)= float(chargeMeson_Final)-0.5*float(hadron(idmeson_Final)%strangeness)-0.5*float(hadron(idmeson_Final)%charm)
       i(2) = float(hadron(IDBaryon_Final)%isoSpinTimes2)/2.
       iz(2)= float(chargeBaryon_Final)-0.5-0.5*float(hadron(idbaryon_Final)%strangeness)-0.5*float(hadron(idbaryon_Final)%charm)
       i(3)=  float(hadron(resID)%isoSpinTimes2)/2.
       isoFactor=isoFactor*ClebschSquared(i(1),i(2),i(3),iz(1),iz(2))

       ! Evaluate the spin factors
       spinFactor=(2.*hadron(resId)%spin+1.)/(2.*hadron(idBaryon_ini)%spin+1.)/(2.*hadron(idMeson_ini)%spin+1.)

       ! Evaluate the partial decay width for incoming and outgoing channel
       ! use inWidth since there are both masses of the incoming particles given
       ! And don't use medium since only the out-width shall be dressed due to the final states

       Gamma_In = partialWidthBaryonMedium(resID,massRes,.true.,IDmeson_ini,IDbaryon_ini, &
                                           momentumLRF,vacuum,baryonMass_ini,mesonMass_ini)

       Gamma_Out = partialwidthBaryonMedium(resID,massRes,.false.,IDmeson_Final,IDbaryon_Final, &
                                            momentumLRF,mediumAtcollision)

       Gamma_Tot = WidthBaryonMedium(resID,massRes,momentumLRF,mediumAtcollision)

       ! Evaluate cross section
       if(fullPropagator.and.isNonExotic(resID)) then
          PI_imag=selfenergy_imag(resID,absMom_LRF,momentumLRF(0),mediumAtCollision)
          PI_real=   get_RealPart(resID,absMom_LRF,fullmass      ,mediumAtCollision)
          sigma = sigma + isoFactor*spinFactor*4.*pi/momCm**2*massRes**2*Gamma_In*Gamma_Out &
                          / ((fullMass**2-hadron(resID)%mass**2-PI_real)**2+PI_imag**2) / GeVSquared_times_mb
       else
          ! neglect real part: old Effenberger ansatz
          sigma = sigma + isoFactor*spinFactor*4.*pi/momCm**2*massRes**2*Gamma_In*Gamma_Out &
                          / ((massRes**2-hadron(resID)%mass**2)**2+Gamma_Tot**2*massRes**2) / GeVSquared_times_mb
       end if

    End do

    barMes_R_barMes=sigma

  end function barMes_R_barMes


  !*************************************************************************
  !****s* resonanceCrossSections/sigma_npi_n2pi_resonances
  ! NAME
  ! subroutine sigma_npi_n2pi_resonances(srts,charge_iniPion,background,sigmaTotal)
  ! PURPOSE
  ! Evaluates contribution of resonances to proton Pion -> Nucleon Pion Pion in the VACUUUM assumption [Gamma's are not the medium modified ones]
  ! One gets pion p -> pi pi N by summing over the (pion p -> R->pion Delta), (pion p -> R-> N rho),
  ! (pion p -> R-> N sigma) and (pion p -> R-> pion P11_1440) channels.
  ! All resonances but the P_11(1440) decay fully into N pi. Therefore the final states are in the end
  ! N Pi Pi states. For the P_11(1440) the ratio of N pi pi final states is given by its decay ratio into N pi.
  !
  ! INPUTS
  ! * logical, intent(in) :: background ---
  !   .true. = Compute background cross section : Therefore consider only non-propagated particles in intermediate state;
  !   .false. = Full Xsection : Consider all particles in intermediate state
  ! * real, intent(in) :: srts ---
  !   sqrt(s) in the process
  ! * integer, intent(in) :: charge_iniPion ---
  !   charge of incoming pion
  !
  ! RESULT
  ! * real, intent(out),dimension(-2:2) :: sigmaTotal ---
  !   cross sction in mB
  !
  ! Meaning of index in sigmaTotal :
  ! * -2 : pi+ pi- in final state
  ! * -1 : pi0 pi- in final state
  ! *  0 : pi0 pi0 in final state
  ! *  1 : pi+ pi0 in final state
  ! *  2 : pi+ pi+ in final state
  !*************************************************************************


  subroutine sigma_npi_n2pi_resonances(srts,charge_iniPion,background,sigmaTotal)
      use baryonWidth, only: partialwidthBaryon, FullWidthBaryon
      use particleProperties, only: hadron, isNonExotic
      use IdTable, only: nucleon, Delta, P11_1440, nbar, pion, rho, sigmaMeson
      use constants, only: pi, mN, mPi, GeVSquared_times_mb
      use selfenergy_baryons, only: get_realPart, selfEnergy_imag
      use mediumDefinition

      logical, intent(in) :: background
     ! .true. = Compute background cross section : Therefore consider only non-propagated particles in intermediate state
     ! .false. = Full Xsection : Consider all particles in intermediate state

      real, intent(in) :: srts
      integer, intent(in)  :: charge_iniPion
      real, intent(out),dimension(-2:2) :: sigmaTotal
      ! Meaning of index in sigmaTotal :
      ! -2 : pi+ pi- in final state
      ! -1 : pi0 pi- in final state
      !  0 : pi0 pi0 in final state
      !  1 : pi+ pi0 in final state
      !  2 : pi+ pi+ in final state

      real, dimension(1:4) :: gamma_Out,sigma
      real :: gammaTot
      real :: momCM ! CM momentum
      integer :: resID,decID
      real :: gamma_In
      real :: spinFactor
      real :: absMom_LRF
      real :: fullMass,PI_REAL,PI_IMAG,energy
      type(medium) :: med


      if(initInput_Flag) then
         call readInput()
         initInput_Flag=.false.
      end if


      absMom_LRF=sqrt(((srts**2-mN**2-mPi**2)/(2*mN))**2-mPi**2)
      if(fullPropagator) then
         fullmass  =srts
         energy=sqrt(srts**2+absMom_LRF**2)
         med%temperature    =0.
         med%useMedium      =.true.
         med%density        = 0
         med%densityProton  = 0
         med%densityNeutron = 0
      end if

      !check input
      If(.not.(abs(charge_iniPion).le.1)) then
         Write(*,*) 'Problems with input in npi_n2pi_resonances', charge_iniPion
      end if

      !Evaluate CM-Momentum
      momCM=Sqrt(max((srts**2-(mPi+mN)**2) &
           &                      *(srts**2-(mPi-mN)**2)  /4./srts**2,1E-8))


      ! Evaluate cross section for all channels
      sigmaTotal(-2:2)=0.
      ! Loop over intermediate resonances R : m B->  R -> m' B'
      Do resId=1,nbar
         If(.not.hadron(resId)%usedForXsections) cycle   ! Exclude resonances
         If(background.and.hadron(resId)%propagated) cycle ! Exclude propagated resonances for background contribution

         ! Evaluate the spin factors
         spinFactor=(2.*hadron(resId)%spin+1.)/(2.*hadron(nucleon)%spin+1.)/(2.*hadron(pion)%spin+1.)
         ! Evaluate the partial decay width for incoming and outgoing channel
         gamma_In=partialwidthBaryon(resID,srts,.true.,pion,nucleon)
         ! Get pion N -> pi pi N by summing just over the ...
         ! pion N -> pion Delta, pion N -> N rho, pion N -> N sigma, pion N -> pion P11_1440
         ! ... channels. All resonances but the P_11(1440) decay fully into N pi.
         ! We correct for the P_11(1440) later.
         gamma_Out(1)=partialwidthBaryon(resID,srts,.false.,pion,Delta)
         gamma_Out(2)=partialwidthBaryon(resID,srts,.false.,rho,nucleon)
         gamma_Out(3)=partialwidthBaryon(resID,srts,.false.,sigmaMeson,nucleon)
         gamma_Out(4)=partialwidthBaryon(resID,srts,.false.,pion,P11_1440)
         ! Formula(2.52) Effenberger Phd.for each channel

         gammaTot=FullWidthBaryon(resID,srts)

         do decId=1,4
            if(fullPropagator.and.isNonExotic(resID)) then
               PI_imag=selfenergy_imag(resID,absMom_LRF,energy        ,med)
               PI_real=   get_RealPart(resID,absMom_LRF,fullmass      ,med)
               sigma(decID) = spinFactor*4.*pi/momCm**2*srts**2*Gamma_In*Gamma_Out(decID) &
                              / ((fullMass**2-hadron(resID)%mass**2-PI_real)**2+PI_imag**2) / GeVSquared_times_mb
            else
               ! neglect real part: old Effenberger ansatz
               sigma(decID) = spinFactor*4.*pi/momCm**2*srts**2*Gamma_In*Gamma_Out(decID) &
                              / ((srts**2-hadron(resID)%mass**2)**2+gammaTot**2*srts**2) / GeVSquared_times_mb
            end if
         end do
         ! correction for N(1440)
         ! multiply the cross section by the N_1440 decay ratio into Nucleon Pion
         ! P11(1440) is the only resonance which has different decay channels than n pi
         ! other channels of P11(1440) might lead to 3 pion nucleon final states.
         sigma(4)=sigma(4)*hadron(P11_1440)%decays(1)

         ! differentiate different isospinChannels and make isospin calculus for the different channnels
         Select Case(hadron(resID)%isoSpinTimes2)
         Case(1) ! Isospin=1/2
            ! not possible to scatter pi^{+} proton via a I=1/2 channel, only contributions of pi^{-} and pi^{0}
            If(charge_iniPion.eq.0) then
               sigmaTotal(0)=sigmaTotal(0)  + 1./3.* (2./9.*sigma(1)+1./3.*sigma(3)+ 1./9.*sigma(4))
               sigmaTotal(-2)=sigmaTotal(-2)+ 1./3.* (5./9.*sigma(1)+1./3.*sigma(2)+ 2./3.*sigma(3)+4./9.*sigma(4))
               sigmaTotal(1)=sigmaTotal(1)  + 1./3.* (2./9.*sigma(1)+2./3.*sigma(2)+ 4./9.*sigma(4))
            else if(charge_iniPion.eq.-1) then
               sigmaTotal(0)=sigmaTotal(0)  + 2./3.* (2./9.*sigma(1)+1./3.*sigma(3)+ 1./9.*sigma(4))
               sigmaTotal(-2)=sigmaTotal(-2)+ 2./3.* (5./9.*sigma(1)+1./3.*sigma(2)+ 2./3.*sigma(3)+4./9.*sigma(4))
               sigmaTotal(-1)=sigmaTotal(-1) + 2./3.* (2./9.*sigma(1)+2./3.*sigma(2)+ 4./9.*sigma(4))
            end if
         case(3) ! Isospin=3/2
            If(charge_iniPion.eq.0) then
                sigmaTotal(0)=sigmaTotal(0)  +  2./3.* (2./45.*sigma(1)+2./9.*sigma(4))
                sigmaTotal(-2)=sigmaTotal(-2)+  2./3.* (26./45.*sigma(1)+2./3.*sigma(2)+2./9.*sigma(4))
                sigmaTotal(1)=sigmaTotal(1)  +  2./3.* (17./45.*sigma(1)+1./3.*sigma(2)+5./9.*sigma(4))
            else if(charge_iniPion.eq.-1) then
                sigmaTotal(0)=sigmaTotal(0)+   1./3.*(2./45.*sigma(1)+2./9.*sigma(4))
                sigmaTotal(-2)=sigmaTotal(-2)+ 1./3.*(26./45.*sigma(1)+2./3.*sigma(2)+2./9.*sigma(4))
                sigmaTotal(-1)=sigmaTotal(-1)+  1./3.*(17./45.*sigma(1)+1./3.*sigma(2)+5./9.*sigma(4))
            else if (charge_iniPion.eq.1) then
                sigmaTotal(1)=sigmaTotal(1)+    13./15.*sigma(1)+sigma(2)+1./3.*sigma(4)
                sigmaTotal(2)=sigmaTotal(2)+     2./15.*sigma(1)+2./3.*sigma(4)
            end if
         end select
      End do
   end subroutine sigma_npi_n2pi_resonances

  end module resonanceCrossSections
