!***************************************************************************
!****m* /mesonWidth
! NAME
! module mesonWidth
! PURPOSE
! When this module is initialized then all information for the VACUUM width 
! is once calculated for all meson resonances and then stored into the
! field gammaField, which is of type gammaFieldType. This is done by 
! initWidth. Afterwards this field is used to return full and 
! partial width of the meson resonances in the vacuum by the subroutine 
! "partialWidthMeson, fullWidthMeson"
! USES
! MesonWidthVacuum
!***************************************************************************
module mesonWidth
  use particleProperties, only: nDecays

  implicit none

  Private

  ! Parameters to store the partial widths depending on mass  in the range 
  ! (minimalmass, minimalMass+maxIndex*deltaMass)
  integer, parameter :: maxIndexMass=2000
  real, parameter :: deltaMass=0.004

  ! The type used for storage:
  type tGammaField 
     real :: gammatotal                ! total decay rate
     real, dimension(nDecays) :: ratio ! ratio of different  decay channels
  end type tGammaField

  ! Field which holds all the information for concerning  the vacuum width, 
  ! initialized in "initWidth"
  ! First Index: ID of Resonance
  ! Second Index: Mass Index in the range (0,maxIndexMass)
  type(tGammaField), dimension(:,:),Allocatable,save  :: gammaField 

  ! Flag to check wether this module is initialized by initWidth
  logical, save :: initFlag = .true.  


  ! Switching debugging infos off and on
  logical, parameter :: debug=.false.

  ! public subroutines :
  PUBLIC :: partialWidthMeson
  PUBLIC :: fullWidthMeson
  PUBLIC :: decayWidthMeson
  PUBLIC :: cleanUp
  PUBLIC :: GetMinMaxMass
  PUBLIC :: GetMaxQ

contains


  !*************************************************************************
  !****f* mesonWidth/partialWidthMeson
  ! NAME
  ! function partialWidthMeson (ID, mass, IDout1, IDout2, IDout3, iQ1, iQ2, iQ3)
  ! PURPOSE
  ! This function calculates the partial width (energy dependent) of all 
  ! meson resonances.
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real :: mass -- sqrt(p_mu p^mu) = mass of the resonance (offshell)
  ! * integer :: IDout1, IDout2 -- IDs of decay products 
  !   (selecting channel of interest)
  ! * integer, OPTIONAL :: IDout3 -- ID of third decay product
  ! * integer, OPTIONAL :: iQ1, iQ2, iQ3 -- Charges of decay products 
  !   (only relevant for 3-body decays)
  !*************************************************************************
  real function partialWidthMeson (ID, mass, IDout1, IDout2, IDout3, iQ1, iQ2, iQ3)
    use DecayChannels, only : Decay2bodyMeson
    use particleProperties, only: hadron, nDecays 
    use idTable, only: pion, eta, photon, isMeson
    use CallStack, only: TRACEBACK

    real,   intent(in) :: mass
    integer,intent(in) :: ID
    integer,intent(in) :: IDout1, IDout2

    ! Only for three-body-channels necessary
    integer, intent(in),optional :: IDout3
    integer, intent(in),optional :: iQ1,iQ2,iQ3

    real  :: Width       
    integer :: massIndex
    integer :: i, dId, dId_wished

    If (initFlag) call initWidth

    !*****************************************
    ! (1) Check Input : 
    if (.not.IsMeson(ID)) call TRACEBACK()
    if (.not.((IsMeson(IDout1).or.(IDout1.eq.photon)))) call TRACEBACK()
    if (.not.((IsMeson(IDout2).or.(IDout2.eq.photon)))) call TRACEBACK()
    If(Present(IDout3)) then
       if (.not.((IsMeson(IDout3).or.(IDout3.eq.photon)))) call TRACEBACK()
       if (.not.(Present(iQ1).and.Present(iQ2).and.Present(iQ3))) call TRACEBACK()
    end If

    !*****************************************
    ! (2) Evaluate the index of the given mass:
    massIndex=NINT((mass-hadron(ID)%minmass)/deltaMass)

    If (mass.le.hadron(ID)%minmass) then
       ! If mass is lower than minimalMas than gamma should be zero
       partialwidthMeson=0.

    else If (massIndex.ge.0) then
       ! Assume constant partial width at very high mass:
       If (massIndex.gt.maxIndexMass) then
!          write(*,*) 'Warning in mesonResonanceWidth/partialWidth'
!          Write(*,*) 'Mass of resonance is out of bounds. Mass=', mass
!          write(*,*) 'ID=',ID
          massIndex=maxIndexMass
       end if

       !*****************************************
       width=0.

       ! (3) decide on three or two body decay
       if (.not. Present(IDout3)) then  !two body
          ! Loop over all decay channels and search for demanded final state
          do i=1,nDecays
             dID = hadron(ID)%decaysID(i)
             if (dId<=0) cycle !  not 2Body
             If ((Decay2BodyMeson(dId)%ID(1)==IDout1 .and. Decay2BodyMeson(dId)%ID(2)==IDout2) .or. &
                 (Decay2BodyMeson(dId)%ID(2)==IDout1 .and. Decay2BodyMeson(dId)%ID(1)==IDout2)) &
                Width=width+gammaField(ID,massIndex)%ratio(i)*gammaField(ID,massIndex)%gammaTotal  
          end do

       else   !three body decay

          dId_wished = 999

          If  ((IDout1.eq.pion).and.(IDout2.eq.pion).and.(IDout3.eq.pion)) then 
             ! channel 2 or 3 :  pion pion pion channels
             If ((iQ1+iQ2+iQ3).eq.0) then
                If((iQ1.ne.0).or.(iQ2.ne.0).or.(iQ3.ne.0)) then
                   dId_wished = -2  ! piPlus, piMinus, piNull channel
                else
                   dId_wished = -3  ! piNull, piNull, piNull channel
                end if
             end if
          else If  ((IDout1.eq.eta).or.(IDout2.eq.eta).or.(IDout3.eq.eta)) then 
             if ((IDout1+IDout2+IDout3).eq.(2*pion+eta)) then
                ! eta pion pion channel
                If((iQ1.eq.0).and.(iQ2.eq.0).and.(iQ3.eq.0)) then
                   dId_wished = -1  ! pi0 pi0 eta
                else If(iQ1+iQ2+iQ3.eq.0) then
                   dId_wished = -4  ! pi+ pi- eta 
                end if
             end if
          end if

          do i=1,nDecays
             dID = hadron(ID)%decaysID(i) ! 3Body: <0
             if (dId.ne.dId_wished) cycle
             Width=width+gammaField(ID,massIndex)%ratio(i)*gammaField(ID,massIndex)%gammaTotal 
          end do

       end if

       partialwidthMeson=width

    else
       write(*,*) 'strange error in partialWidthMeson. STOP!', massIndex, mass,hadron(ID)%minmass,ID
       stop
       
    end if
    
  end function partialWidthMeson
  


  !*************************************************************************
  !****f* mesonWidth/FullWidthMeson
  ! NAME
  ! real function FullWidthMeson(ID,mass)
  ! PURPOSE
  ! This function calculates the full width (energy dependent) of all meson 
  ! resonances.
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real :: mass -- sqrt(p_mu p^mu) = mass of the resonance (offshell)
  !*************************************************************************
  real function FullWidthMeson(ID,mass)
    use IdTable, only: isMeson
    use particleProperties, only: hadron
    use CallStack, only: TRACEBACK

    integer, intent(in) :: ID
    real,    intent(in) :: mass

    integer :: down
    real    :: mass_down,weight

    If (initFlag) call initWidth

    ! (1) Check Input
    If (.not. isMeson(ID)) call TRACEBACK()

    ! (2) Return the full width
    If (mass <= hadron(ID)%minmass) then
       ! If mass is lower than minimalMass, then gamma should be zero
       FullWidthMeson=0.
    else If (mass > hadron(ID)%minmass+deltaMass*(maxIndexMass-1)) then
!       write(*,*) 'Warning in mesonResonanceWidth/widht'
!       Write(*,*) 'Mass of resonance is out of bounds. Mass=', mass
!       write(*,*) 'ID=',ID
       FullWidthMeson=gammaField(ID,maxIndexMass)%gammaTotal
    else if (mass > hadron(ID)%minmass) then
       ! Do linear interpolation between next two grid points "down" and "down+1"
       down = floor((mass-hadron(ID)%minmass)/deltaMass)
       mass_down = hadron(ID)%minmass+float(down)*deltaMass
       weight = (mass-mass_down)/deltaMass ! weight for the interpolation
       FullWidthMeson = gammaField(ID,down  )%gammaTotal * (1.-weight) &
                      + gammaField(ID,down+1)%gammaTotal * weight
    else
      print *,"problem in fullWidthMeson:", ID, mass
      call TRACEBACK()
    end if

    !If(debug) Print *, "In mesonWidth", mass,FullWidthMeson

  end function fullWidthMeson


  !*************************************************************************
  !****s* mesonWidth/initWidth
  ! NAME
  ! subroutine initWidth
  ! PURPOSE
  ! Stores the vacuum width of each meson to the field "gammaField".
  ! Should be called only once.
  !*************************************************************************
  subroutine initWidth
    use mesonWidthVacuum, only: vacuumWidth
    use particleProperties, only: nDecays, hadron, get_rho_dilep
    use IdTable, only: pion, rho, nMes

    integer :: massIndex, ID
    real :: mass, gammaTotal
    real, dimension(nDecays) :: ratio 

    ! allocate the field which holds the decay ratio information for each 
    ! meson, depending on mass, 
    ! first index: meson ID
    ! second index : mass 
    Allocate(gammaField(pion:pion+nMes-1,0:maxIndexMass))

    ! Initialize the gamma fields for the mesons by calling vacuumWidth
    do ID=pion,pion+nMes-1
       If(debug) print *, "Resonance=", ID
       do massIndex=0,MaxIndexMass
          mass=real(massIndex)*deltaMass+hadron(ID)%minmass
          gammaTotal = vacuumWidth (mass, ID, ratio)
          gammaField(ID,MassIndex)%gammaTotal=gammaTotal
          gammaField(ID,MassIndex)%ratio=ratio
          If (debug) Write(100+Id,'(6F14.9)') mass, gammaTotal 
          if (.not. (ID==rho .and. get_rho_dilep()) .and. (abs(sum(ratio)-1)>1E-6) .and. (abs(sum(ratio))>1E-6)) then
             Write(*,*) 'Problem in mesonWidth/initWidth'
             write(*,*) "Ratios don't add up to 1! ",gammaTotal
             write(*,*) "Ratio :  ", ratio
             write(*,*) 'ResonanceID:', ID
          end if
       end do
    end do

    initFlag=.false.

  end subroutine initWidth


  !*************************************************************************
  !****s* mesonWidth/cleanUp
  ! subroutine cleanUp
  ! PURPOSE
  ! Deallocate all fields
  !*************************************************************************
  subroutine cleanUp()
    if (initFlag .or. .not. Allocated(gammaField)) return
    DeAllocate(gammaField)
  end subroutine


  !*************************************************************************
  !****f* mesonWidth/decayWidthMeson
  ! NAME
  ! function decayWidthMeson (ID, mass, ch) result(decayWidth)
  ! PURPOSE
  ! This function returns the partial out width (energy dependent) for all decay channels.
  !
  ! The dStar resonances are treated explicitly since they are the only resonances 
  ! which have decay ratios which depend on the charge of the resonance.
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real    :: mass       -- baremass of the resonance (offshell)
  ! * integer :: ch         -- charge of the resonance
  ! OUTPUT
  ! * real, dimension(nDecays)  :: decayWidth -- widths
  !*************************************************************************
  function decayWidthMeson (ID, mass, ch) result(decayWidth)
    use IDTable, only: dStar, dStarBar, IsMeson
    use particleProperties, only: hadron, nDecays
    use decayChannels, only: decay2BodyMeson, decay3BodyMeson
    use CallStack, only: TRACEBACK

    integer,intent(in) :: ID
    real,   intent(in) :: mass
    integer,intent(in) :: ch
    real, dimension(nDecays) :: decayWidth

    integer :: i, dID, down
    real :: thr, mass_down, weight

    If (initFlag) call initWidth

    ! Check Input
    If (.not. isMeson(ID)) call TRACEBACK()

    ! Initialize to zero
    decayWidth = 0.

    ! Decays of the D* Mesons:
    ! Treated seperately since the D* is very different in the decays in the different charge modes.
    select case (ID)
    case (dStar)
       If (ch == 0) then
          decayWidth(1) = 0.002*0.381 ! photon+dMeson
          decayWidth(2) = 0.002*0.619 ! pion  +dMeson
       else
          decayWidth(1) = 0.096/1000*0.016 ! photon+dMeson
          decayWidth(2) = 0.096/1000*0.984 ! pion  +dMeson
       end if
       return
    case (dStarBar)
       If(ch == 0) then
          decayWidth(1) = 0.002*0.381 ! photon+dBar
          decayWidth(2) = 0.002*0.619 ! pion  +dBar
       else
          decayWidth(1) = 0.096/1000.*0.016 ! photon+dBar
          decayWidth(2) = 0.096/1000.*0.984 ! pion  +dBar
       end if
       return
    end select

    do i=1,nDecays
      dID = hadron(ID)%decaysID(i)
      if (dID==0) then
        cycle
      else if (dID>0) then
        thr = decay2BodyMeson(dId)%threshold
      else
        thr = decay3BodyMeson(-dId)%threshold
      end if
      If (mass < thr) then
        decayWidth(i) = 0.
      else If (mass > hadron(ID)%minmass+deltaMass*(maxIndexMass-1)) then
        decayWidth(i) = gammaField(ID,maxIndexMass)%gammaTotal * gammaField(ID,maxIndexMass)%ratio(i)
      else
        ! Do linear interpolation between next two grid points "down" and "down+1"
        down = floor((mass-hadron(ID)%minmass)/deltaMass)
        mass_down = hadron(ID)%minmass+float(down)*deltaMass
        weight = (mass-mass_down)/deltaMass ! weight for the interpolation
        decayWidth(i) = gammaField(ID,down  )%gammaTotal * gammaField(ID,down  )%ratio(i) * (1.-weight) &
                      + gammaField(ID,down+1)%gammaTotal * gammaField(ID,down+1)%ratio(i) * weight
      end if
    end do

  end function decayWidthMeson


  !*************************************************************************
  !****s* mesonWidth/GetMinMaxMass
  ! NAME 
  ! subroutine GetMinMaxMass(ID,MinMass,MaxMass,InMedium)
  ! PURPOSE
  ! return values of minimal and maximal mass according the mass tabulation
  ! INPUTS
  ! * integer :: ID -- ID of particle
  ! * logical :: InMedium -- Flag to override minimal mass of vector mesons
  ! OUTPUT
  ! * real :: MinMass -- minimal mass value
  ! * real :: MaxMass -- maximal mass value
  ! NOTES
  ! This returns the minimal mass as stored as default value. 
  ! For in-medium vector mesons, the minimal mass is reduced to zero.
  !
  ! The maximal mass is given by the size of the array times the bin width.
  ! All masses are restricted by an upper bound (=3 GeV).
  !*************************************************************************
  subroutine GetMinMaxMass(ID,MinMass,MaxMass,InMedium)
    use particleProperties, only: hadron
    use IdTable, only: rho,omegaMeson,phi

    integer, intent(in) :: ID
    logical, intent(in) :: InMedium
    real, intent(out) :: MinMass,MaxMass
    
    MinMass = hadron(ID)%minmass
    if (inMedium) then
       select case (ID)
       case (rho,omegaMeson,phi)
          MinMass = 0.01
       end select
    end if

    MaxMass = MinMass + deltaMass*(maxIndexMass-1)

    MaxMass = 3.0
!    if (ID==107) maxMass = 2.0

  end subroutine GetMinMaxMass

  !*************************************************************************
  !****s* mesonWidth/GetMaxQ
  ! NAME
  ! subroutine GetMaxQ(ID,mass0,gamma0,gammaDelta,BinM,BinMaxQ)
  ! PURPOSE
  ! Calculate the maximal values of the Q weight for bins according BinM
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real :: mass0 -- pole mass
  ! * real :: gamma0 -- width at pole mass
  ! * real :: gammaDelta -- additional width to be added during calculations
  ! * real, dimension(:) :: BinM -- array with boundaries for M binning
  ! OUTPUT
  ! * real, dimension(:) :: BinMaxQ -- the maximal Q values for each bin.
  !
  ! NOTES
  ! * The size of BinMaxQ has to be at least the size of BinM minus 1.
  ! * It first calculates Q at the boundaries, then it iterates over the 
  !   tabulated width values in order to take into account, that the Q
  !   value may be larger inbetween the boundaries.
  ! * if the Q value is maximal at the upper bound, we store its value as
  !   -Q.
  !
  !************************************************************************
  subroutine GetMaxQ(ID,mass0,gamma0,gammaDelta,BinM,BinMaxQ)
    use MassAssInfoDefinition, only: MassAssInfoQ

    integer, intent(in) :: ID
    real, intent(in) :: mass0,gamma0,gammaDelta
    real, dimension(:),intent(in)  :: BinM
    real, dimension(:),intent(out) :: BinMaxQ

    integer :: iB, iM,iM1,iM2
    real :: mass,gamma, Q

    ! 1) calculate Q at the boundaries:

    do iB=1,size(BinM)
       gamma = FullWidthMeson(ID,BinM(iB))
       BinMaxQ(iB) = MassAssInfoQ(mass0,gamma0+gammaDelta,BinM(iB),gamma+gammaDelta)
    end do

    do iB=1,size(BinM)-1
       if (BinMaxQ(iB+1).gt.BinMaxQ(iB)) BinMaxQ(iB) = -BinMaxQ(iB+1) ! sign!!
    end do

    ! 2) calculate Q between the boundaries:

    do iB=1,size(BinM)-1
       iM1 = INT((BinM(iB)  -BinM(1))/deltaMass)
       iM2 = INT((BinM(iB+1)-BinM(1))/deltaMass)

       do iM=iM1+1,iM2
          mass  = iM*deltaMass+BinM(1)
          gamma = gammaField(ID,iM)%gammaTotal
          Q = MassAssInfoQ(mass0,gamma0+gammaDelta,mass,gamma+gammaDelta)
          if (Q.gt.abs(BinMaxQ(iB))) BinMaxQ(iB) = Q
       end do
    end do

  end subroutine GetMaxQ


end module mesonWidth
