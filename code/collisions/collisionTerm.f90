!*******************************************************************************
!****m* /collisionTerm
! NAME
! module collisionTerm
! PURPOSE
! Includes the collision term.
!*******************************************************************************
module collisionTerm

  implicit none
  Private
  
  !*****************************************************************************
  !****g* collisionTerm/useStatistics
  ! SOURCE
  ! 
  logical, save :: useStatistics=.false.
  ! PURPOSE
  ! Generate statistical information using the module statistics.
  !*****************************************************************************


  !*****************************************************************************
  !****g* collisionTerm/noNucNuc
  ! SOURCE
  !
  logical, save :: noNucNuc=.false.
  ! PURPOSE
  ! Switch off perturbative  NN reactions.
  !*****************************************************************************


  !*****************************************************************************
  !****g* collisionTerm/energyCheck
  ! SOURCE
  ! 
  real, save :: energyCheck=0.01   
  ! PURPOSE
  ! Precision of energy check for each collision in GeV.
  !*****************************************************************************


  !*****************************************************************************
  !****g* collisionTerm/oneBodyProcesses
  ! SOURCE
  ! 
  logical,save :: oneBodyProcesses=.true. 
  ! PURPOSE
  ! Switch to turn off one-body-induced processes.
  !*****************************************************************************


  !*****************************************************************************
  !****g* collisionTerm/oneBodyAdditional
  ! SOURCE
  ! 
  logical,save :: oneBodyAdditional=.true. 
  ! PURPOSE
  ! Switch to turn off additional Pythia one-body-induced processes.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/twoBodyProcesses
  ! SOURCE
  ! 
  logical,save :: twoBodyProcesses=.true. 
  ! PURPOSE
  ! Switch to turn off two-body-induced processes.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/twoBodyProcessesRealReal
  ! SOURCE
  ! 
  logical,save :: twoBodyProcessesRealReal=.true. 
  ! PURPOSE
  ! Switch to turn of two-body-induced processes between two real particles.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/twoBodyProcessesRealPert
  ! SOURCE
  ! 
  logical,save :: twoBodyProcessesRealPert=.true. 
  ! PURPOSE
  ! Switch to turn off two-body-induced processes between a 
  ! real and a perturbative particle.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/threeBodyProcesses
  ! SOURCE
  ! 
  logical,save :: threeBodyProcesses=.true. 
  ! PURPOSE
  ! Switch to turn off three-body-induced processes.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/twoPlusOneBodyProcesses
  ! SOURCE
  ! 
  logical,save :: twoPlusOneBodyProcesses=.false. 
  ! PURPOSE
  ! Switch to turn off 2+1 body processes 
  ! (two really colliding particles plus one nearby).
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/maxOut
  ! SOURCE
  ! 
  integer,save :: maxOut=100 
  ! PURPOSE
  ! Maximal number of produced particles in one process.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/printPositions
  ! SOURCE
  ! 
  logical,save :: printPositions=.false. 
  ! PURPOSE
  ! Switch to turn of output of positions in real-pert collisions. 
  ! Produces statistical output.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/storeRho0Info
  ! SOURCE
  !
  logical, save :: storeRho0Info = .false.
  ! PURPOSE
  ! Flag whether in a rho0 decay the particle numbers of the 
  ! resulting charged pions are stored or not.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/storeRho0InfoOnlyDifr
  ! SOURCE
  !
  logical, save :: storeRho0InfoOnlyDifr = .false.
  ! PURPOSE
  ! Flag, whether the flag storeRho0Info is valid for all decays or
  ! only for rho0, which are marked to be diffractive.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/DoJustAbsorptive
  ! SOURCE
  !
  logical, save :: DoJustAbsorptive = .false.
  ! PURPOSE
  ! If this flag is true, then:
  ! for perturbative simulations all final state particles in a collision are set to zero;
  ! for real simulations %event index of incoming hadron is changed in the case of collision,
  ! but actual collision is not simulated. 
  ! This is a way to mimick Glauber like calculations.
  ! NOTES
  ! The "absorption" is done with sigmaTot, not just by sigmaInEl.
  !*****************************************************************************


  !*****************************************************************************
  !****g* collisionTerm/annihilate
  ! SOURCE
  !
  logical, save :: annihilate = .false.
  ! PURPOSE
  ! If this flag is true, then an annihilation of the antibaryons with 
  ! the closest baryons will be simulated (by hand) starting from 
  ! annihilationTime.
  !*****************************************************************************


  !*****************************************************************************
  !****g* collisionTerm/annihilationTime
  ! SOURCE
  !
  real, save :: annihilationTime = 1000.
  ! PURPOSE
  ! Time moment (in fm/c) when the annihilation will be started.
  ! NOTES
  ! This flag has an influence only when annihilate = .true.
  ! Before annihilationTime all the collision processes are not
  ! activated. They start to act (if the corresponding switches 
  ! oneBodyProcesses,twoBodyProcesses etc. are .true.) only after
  ! annihilationTime.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/justDeleteDelta
  ! SOURCE
  !
  logical, save :: justDeleteDelta = .false.
  ! PURPOSE
  ! Deletes final-state products in Delta N N -> NNN and Delta N -> N N.
  ! NOTE: Only for testing and comparing with the old Effenberger code!
  ! DO NOT USE OTHERWISE: Violates energy conservation!
  ! This switch is meant to simulate the treatment of the Delta in the old code.
  ! Only implemented for perturbative runs.
  !*****************************************************************************

  !*****************************************************************************
  !****g* collisionTerm/debug
  ! SOURCE
  logical,save  :: debug=.false.
  ! PURPOSE
  ! Switch for debug information.
  !*****************************************************************************

  logical, save :: debugSaveInfo=.false.
  logical, save :: initFlag=.true.

!!$  type(histogram), save :: hXStot,hXSelast

  Public :: collideMain, finalCheck, ForceDecays


contains


  !*****************************************************************************
  !****s* collisionTerm/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "collisionTerm".
  !*****************************************************************************
  subroutine readInput

    use output, only: Write_ReadingInput

    integer :: ios

    !***************************************************************************
    !****n* collisionTerm/collisionterm
    ! PURPOSE
    ! Namelist for collisionTerm includes:
    ! * oneBodyProcesses
    ! * twoBodyProcesses
    ! * threeBodyProcesses
    ! * twoPlusOneBodyProcesses
    ! * twoBodyProcessesRealReal
    ! * twoBodyProcessesRealPert
    ! * oneBodyAdditional
    ! * energyCheck
    ! * maxOut
    ! * debug
    ! * printPositions
    ! * useStatistics
    ! * noNucNuc
    ! * storeRho0Info
    ! * storeRho0InfoOnlyDifr
    ! * DoJustAbsorptive
    ! * annihilate
    ! * annihilationTime
    ! * justDeleteDelta
    !***************************************************************************
    NAMELIST /collisionTerm/ oneBodyProcesses, twoBodyProcesses, threeBodyProcesses, twoPlusOneBodyProcesses, &
                             twoBodyProcessesRealReal, twoBodyProcessesRealPert, oneBodyAdditional, &
                             energyCheck, maxOut, debug, printPositions, useStatistics, noNucNuc, &
                             storeRho0Info, storeRho0InfoOnlyDifr, DoJustAbsorptive, &
                             annihilate, annihilationTime, justDeleteDelta

    rewind(5)
    call Write_ReadingInput('collisionTerm',0)
    read(5,nml=collisionTerm,iostat=ios)

    if (ios>0) then
       write(*,*)
       write(*,*) '>>>>> parameter "minimumEnergy" moved to the new namelist "insertion"'
    endif
    call Write_ReadingInput('collisionTerm',0,ios)
    write(*,*) '  Precision of energy Check in GeV=',energyCheck

    write(*,*) '  Performing 1-body-induced reactions :', oneBodyProcesses

    if (oneBodyAdditional)  Write(*,*) '  + Perform additional Jetset Decays.'

    if (twoBodyProcesses) then
       write(*,*) '  Performing 2-body-induced reactions :', twoBodyProcesses, &
            '  ...RealReal,RealPert=',twoBodyProcessesRealReal,twoBodyProcessesRealPert
    else
       write(*,*) '  Performing 2-body-induced reactions :', twoBodyProcesses
    end if
    write(*,*) '  Performing 3-body-induced reactions :', threeBodyProcesses
    write(*,*) '  Performing 2+1 body collisions :',  twoPlusOneBodyProcesses
    write(*,*) '  Maximal number of final state particles :', maxout
    write(*,*) '  Statistical output?' ,useStatistics
    write(*,*) '  Switching off perturbative NN collisions:',noNucNuc
    write(*,*) '  Store pion numbers in rho0 decay:        ',storeRho0Info
    write(*,*) '                -"-      only diffractive: ',storeRho0InfoOnlyDifr

    if (DoJustAbsorptive) then 
       write(*,*)
       write(*,*) '  ATTENTION: absorptive treatment !!!!'
       write(*,*)
    endif


    if (annihilate) then 
       write(*,*)
       write(*,*) '  ATTENTION: annihilation of the antibaryons by hand !!!!'
       write(*,*) '  annihilation time:  ', annihilationTime
       write(*,*)
    endif

    if (justDeleteDelta) then 
       write(*,*)
       write(*,*) '  ATTENTION: Delta N N -> NNN and Delta N -> NN: Finalstate not populated !!!!'
       write(*,*) '  ATTENTION: ONLY FOR TESTING,  NOT FOR REGULAR USE !!!!'
       write(*,*)
    endif

    call Write_ReadingInput('collisionTerm',1)
    initFlag = .false.

!!$    call CreateHist(hXStot,   'XS tot',   0.,10.,0.05)
!!$    call CreateHist(hXSelast, 'XS elast', 0.,10.,0.05)

  end subroutine readInput


  !*****************************************************************************
  !****s* collisionTerm/ForceDecays
  ! NAME
  ! subroutine ForceDecays (teilchenPert, teilchenReal, time)
  ! PURPOSE
  ! This routine forces all instable particles to decay in the last timestep.
  ! We loop over the particle vectors several times in order to allow for
  ! multi-step decay chains and to make sure that everything has decayed.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: teilchenReal -- real particles 
  ! * type(particle), dimension(:,:) :: teilchenPert -- perturbative particles 
  ! * real                           :: time
  !*****************************************************************************
  subroutine ForceDecays (teilchenPert, teilchenReal, time)
    use particleDefinition
    use output, only: DoPr
    use AddDecay, only: PerformAddDecay
    use Dilepton_Analysis, only: Dilep_Decays
    use inputGeneral, only: delta_t, eventType
    use eventTypes, only: HeavyIon, hadron

    type(particle), dimension(:,:) :: teilchenReal, teilchenPert
    real, intent(in) :: time
    integer :: i,l
    logical :: dec

    !       call WriteHist(hXStot,121,file="CollTerm.XS.tot.dat",DoAve=.true.)
    !       call WriteHist(hXSElast,121,file="CollTerm.XS.elast.dat",DoAve=.true.)
         
    If (oneBodyProcesses) then 
      ! At the end of all time steps the resonances are forced to decay:
      Do i=1,100
        ! Evaluate final decays
        if (DoPr(2)) Write(*,*) 'No usual Collision Term; Just particle decays: ', i
        dec = oneBody(teilchenPert,teilchenReal,time,.true.)

        ! final Dilepton Analysis during/after forced decays
        if (dec) then
          l = 1
        else
          l = 2
        end if
        if (eventType==HeavyIon .or. eventType==hadron) then
          call Dilep_Decays (time+i*delta_t, teilchenReal, l)
        else
          call Dilep_Decays (time+i*delta_t, teilchenPert, l)
        end if
        
        if (.not. dec) exit  ! no further decays have occurred => exit the loop
      end do

      if (oneBodyAdditional) then
        if (DoPr(2)) Write(*,*) 'Perform additional Jetset Decays.'
        call PerformAddDecay(teilchenReal,time)
        call PerformAddDecay(teilchenPert,time)
      end if
    end if

  end subroutine


  !*****************************************************************************
  !****s* collisionTerm/collideMain
  ! NAME
  ! subroutine collideMain (teilchenPert, teilchenReal, time)
  ! PURPOSE
  ! Evaluates the collision term for a given real and perturbative particle vector. 
  ! There are reactions induced by a single particle (decays), 
  ! two particles and even three particles.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: teilchenReal -- real particles 
  ! * type(particle), dimension(:,:) :: teilchenPert -- perturbative particles 
  ! * real                           :: time
  !*****************************************************************************
  subroutine collideMain (teilchenPert, teilchenReal, time)
    use particleDefinition
    use output, only: DoPr, paragraph
    use deuterium_Pl, only : deuterium_PertOrigin!,deuterium_pertOrigin_flag
    use inputGeneral, only :localEnsemble
    use CollHistory, only: CollHist_SetSize
    use AntibaryonWidth, only: DoAnnihilation
    use CallStack, only: Traceback

    type(particle), intent(inout), dimension(:,:) :: teilchenPert, teilchenReal
    real, intent(in) :: time
    
    logical :: dec

    ! Initialize at first call
    If (initFlag) call ReadInput

    ! For security:
    teilchenReal%perturbative=.false.
    teilchenPert%perturbative=.true.

    If (size(teilchenPert,dim=1) /= size(teilchenReal,dim=1)) then
      write(*,*) 'Number of ensembles in real and perturbative particle vectors do not fit' 
      write(*,*) 'Real :' , size(teilchenReal,dim=1)
      write(*,*) 'Perturbative :' , size(teilchenPert,dim=1)
      call Traceback('Critical Error in collideMain! Stop!')
    end if
    
    if (annihilate) call DoAnnihilation(teilchenReal,time,annihilationTime)

    ! prepare for extra collision history
    call CollHist_SetSize(ubound(teilchenPert))

    ! Evaluate decays
    If (oneBodyProcesses) then
      if (DoPr(2)) Write(*,paragraph)   "1-Body Processes"
      dec = oneBody (teilchenPert, teilchenReal, time, .false.)
    end if

    ! For security:
    teilchenReal%perturbative=.false.
    teilchenPert%perturbative=.true.

    ! Evaluate two-body collisions
    If (twoBodyProcesses)  then
      if (DoPr(2)) Write(*,paragraph)   "2-Body Processes"
      If (localEnsemble) then
        call twoBody_local (teilchenPert, teilchenReal, time)
      else
        call twoBody (teilchenPert, teilchenReal, time)
       end If
    end if

    ! For security:
    teilchenReal%perturbative=.false.
    teilchenPert%perturbative=.true.

    ! Evaluate three-body interactions
    If (threeBodyProcesses)  then
      if (DoPr(2)) Write(*,paragraph)   "3-Body Processes"
      call threeBody (teilchenPert, teilchenReal, time)
    end if
    
    ! For security:
    teilchenReal%perturbative=.false.
    teilchenPert%perturbative=.true.
    
    nullify(deuterium_pertOrigin)

  end subroutine collideMain
  
  

  !*****************************************************************************
  !****s* collisionTerm/oneBody
  ! NAME
  ! logical function oneBody (teilchenPert, teilchenReal, time, ForceFlag) 
  ! PURPOSE
  ! Administrates the 1-body processes (i.e. decays).
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: teilchenPert -- perturbative particles
  ! * type(particle),dimension(:,:) :: teilchenReal -- real particles 
  ! * real    :: time      -- actual time step
  ! * logical :: ForceFlag --
  !   .true. = let particles decay (decay probability=1), 
  !   do not consider the width anymore for the decay probability.
  !   This is useful at the end of a run.
  !
  ! OUTPUT
  ! * teilchenPert and teilchenReal are changed.
  ! * The return value indicates whether any decays have occurred.
  !*****************************************************************************
  logical function oneBody (teilchenPert, teilchenReal, time, ForceFlag) 

    use master_1Body, only : decayParticle
    use inputGeneral, only : fullEnsemble, numEnsembles
    use collisionNumbering, only: ReportEventNumber
    use particleDefinition
    use output, only : WriteParticleVector
    use IdTable, only: isMeson, isBaryon
    use deuterium_PL, only : deuterium_pertOrigin,deuterium_pertOrigin_flag
    use pauliBlockingModule, only: checkPauli
    use twoBodyStatistics, only : rate
    use insertion, only: particlePropagated, setIntoVector
    use CallStack, only: Traceback
    
    type(particle), intent(inout), dimension(:,:) :: teilchenPert, teilchenReal
    real, intent(in) :: time
    logical, intent(in) :: ForceFlag 
    
    integer :: index,ensemble
    type(particle),target         :: resonance   ! particle which shall decay
    type(particle),dimension(1:3) :: finalState  ! final state of particles
    logical :: Flag, setFlag, pauliFlag

    oneBody = .false.

    If(debug) write(*,*) 'real decays'
    call RealDecay
    If(debug) write(*,*) 'perturbative decays'
    call PertDecay
    If(debug) write(*,*) 'after perturbative decays'

  contains

    !************************************
    ! real particles decaying
    !************************************
    subroutine RealDecay
      use idTable, only : photon

      logical :: flag1
      integer :: i

      ensemble_loop : Do ensemble=1,numEnsembles             !loops over ensemble
         index_loop : Do index=1,size(teilchenReal,dim=2)

            If (teilchenReal(ensemble,index)%ID < 0) exit index_loop !!! ID==-1
            If (teilchenReal(ensemble,index)%ID == 0) cycle index_loop
            If (teilchenReal(ensemble,index)%ID.eq.photon) cycle index_loop

            decay_loop: Do

            resonance=teilchenReal(ensemble,index)
            deuterium_pertOrigin=>resonance
            deuterium_pertOrigin_flag=11

            call decayParticle(resonance,finalState,flag,pauliFlag,ForceFlag,time)
            If ((.not. PauliFlag) .and. (.not. ForceFlag)) then
               ! Note : a) PauliFlag is set if the Pauli Blocking was already considered in the decay width
               !        b) If particles are forced to decay then we forget about Pauli blocking
               if (.not. flag) cycle index_loop
               call checkPauli(finalState,teilchenReal,flag)
            end if
            if (.not. flag) cycle index_loop
            
            flag = finalCheck ((/resonance/), finalState)
            if (.not. flag) cycle index_loop

            oneBody = .true.  ! at this point we know that the decay is happening

            call rate((/resonance/),finalState,time)  ! Compute various collision rates

            ! (2) Create final state particles:
            flag1=.false.
            do i=1,3
              if (finalState(i)%id <= 0) exit
              If (particlePropagated(finalState(i))) then
                 if (.not.flag1 .and. ((isBaryon(resonance%ID) .and. isBaryon(finalState(i)%ID) .and. &
                                        (resonance%antiparticle.eqv.finalState(i)%antiparticle)) .or. &
                                       (isMeson(resonance%ID) .and. isMeson(finalState(i)%ID)))) then
                    call setNumber(finalState(i))
                    teilchenReal(ensemble,index)=finalState(i)
                    flag1=.true.
                 else
                     ! Find empty space in the particle vector: 
                     If (fullensemble) then
                        call setIntoVector(finalState(i:i),teilchenReal,setFlag)
                     else
                        call setIntoVector(finalState(i:i),teilchenReal(ensemble:ensemble,:),setFlag)
                     end if
                     ! (3) Check that setting into real particle
                     !     vector worked out :
                     If (.not. setFlag) then
                        write(*,*) 'Real particle vector too small!'
                        write(*,*) size(finalState),&
                        lBound(teilchenReal), uBound(teilchenReal)
                        write(*,*) 'Dumping real particles to files "RealParticles_stop_*!"'
                        call WriteParticleVector("RealParticles_stop",teilchenReal)
                        write(*,*) 'Dumping perturbative particles to files "PertParticles_stop_*!"'
                        call WriteParticleVector("PertParticles_stop",teilchenPert)
                        call Traceback('collision Term, one-body')
                     end if
                 end if
              end If
            end do

            call ReportEventNumber((/resonance/),finalState, finalState(1)%event,time,11)

            ! (4) Check whether leading particle exists:
            if(.not.flag1) then            
!!$               write(*,*)'In onebody: no outgoing leading particle found !!!'
!!$               write(*,*) resonance(1)%ID,finalState%ID
               ! A.L.: This corresponds to the meson decay to photons.
               ! Thus, we must destroy the decaying particle:
               teilchenReal(ensemble,index)%ID=0               
               exit decay_loop
            end if

            ! (5) Decay also the outgoing particle if it is unstable 
            !     at the end of time evolution, otherwise go to the next 
            !     particle:
            if(.not.ForceFlag) exit decay_loop

            End do decay_loop ! Loop over possible decay of outgoing resonances

         End do index_loop ! Index of particle
      End do ensemble_loop ! Ensemble of particle
      nullify(deuterium_pertOrigin)

    end subroutine RealDecay

    !************************************
    ! perturbative decays
    !************************************
    subroutine PertDecay
      use idTable, only: photon, pion, rho, nucleon
      use statistics, only : saveInfo
      use PIL_rho0Dec, only: PIL_rho0Dec_PUT
      use CollHistory, only: CollHist_UpdateHist

      integer :: k
      logical :: r!,rr
      integer, dimension(2,3) :: posOut

      ensemble_loop : Do ensemble=1,numEnsembles             !loop over first particle=real particle
         index_loop : Do index=1,size(teilchenPert,dim=2)

            If (teilchenPert(ensemble,index)%ID < 0) exit index_loop !!! ID==-1
            If (teilchenPert(ensemble,index)%ID == 0) cycle index_loop
            If (teilchenPert(ensemble,index)%ID.eq.photon) cycle index_loop

            resonance=teilchenPert(ensemble,index)
            deuterium_pertOrigin=>resonance
            deuterium_pertOrigin_flag=1

            call decayParticle(resonance,finalState,flag,pauliFlag,ForceFlag,time)
            If((.not.PauliFlag).and.(.not.ForceFlag)) then   
               ! Note : a) PauliFlag is set if the Pauli Blocking was already considered in the decay width
               !        b) If particles are forced to decay then we forget about Pauli blocking
               if (.not.flag) cycle index_loop
               call checkPauli(finalState,teilchenReal,flag)
            end if
            if (.not.flag) cycle index_loop

            flag = finalCheck ((/resonance/), finalState)
            if (.not.flag) cycle index_loop

            oneBody = .true.  ! at this point we know that the decay is happening
           
            ! (1a) some statistics

            do k=lbound(finalstate,dim=1),ubound(finalstate,dim=1)
               if(finalState(k)%ID.le.0) cycle
               call setNumber(finalState(k))
               ! Statistical output
               if (printPositions .and. (finalState(k)%id==pion .or. finalState(k)%id==nucleon)) &
                 call positionPrinter (finalstate(k))
               if (useStatistics)  then
                  call saveInfo(finalState(k)%number,ensemble,finalstate(k)%position, resonance%ID,0,0)
                  if(debugSaveInfo) write(98,*)finalState(k)%number,ensemble,finalstate(k)%position, resonance%ID
               end if
            end do

            ! (1b) Store information for a rho0 Decay
            if (storeRho0Info .and. resonance%ID==rho .and. resonance%charge==0) then
               r = .true.
               !if (storeRho0InfoOnlyDifr)  rr = PIL_rhoDiffractive_GET(resonance(1)%number,r)

               if (r) then
                  call PIL_rho0Dec_PUT (finalState(1)%number, finalState(2)%number, &
                                        resonance%lastCollisionTime, resonance%productionTime, resonance%formationTime)
                  call PIL_rho0Dec_PUT (finalState(2)%number, finalState(1)%number, &
                                        resonance%lastCollisionTime, resonance%productionTime, resonance%formationTime)
               endif
            endif

            ! (2) Eliminate incoming particle
            teilchenPert(ensemble,index)%Id = 0
            ! (3) Create final state particles
            posOut = 0
            If (particlePropagated(finalState(1))) then
               teilchenPert(ensemble,index) = finalState(1)
               posOut(1:2,1) = (/ensemble,index/)
            end if
            If (fullensemble) then
               call setIntoVector (finalState(2:), teilchenPert, setFlag, .true., positions=posOut(:,2:))
            else
               call setIntoVector (finalState(2:), teilchenPert(ensemble:ensemble,:), setFlag, .true., positions=posOut(:,2:))
               posOut(1,2:) = ensemble
            end if

            call ReportEventNumber ((/resonance/), finalState, finalState(1)%event, time, 12)

            call CollHist_UpdateHist ((/resonance/), finalState, (/ensemble,index/), posOut, resonance%perweight)

            ! (4) Check that setting into perturbative particle vector worked out 
            If (.not.setFlag) then
               write(*,*) 'Perturbative particle vector too small!'
               write(*,*) size(finalState),  lBound(teilchenPert), uBound(teilchenPert)
               write(*,*) 'Dumping real particles to files "realParticles_stop_*!"'
               call WriteParticleVector("RealParticles_stop",teilchenReal)
               write(*,*) 'Dumping perturbative particles to files "PertParticles_stop_*!"'
               call WriteParticleVector("PertParticles_stop",teilchenPert)
               call Traceback('collision Term, one-body')
            end if

            If (debug) write(*,*)  'After Set into vector'

         End do index_loop ! Index of particle
      End do ensemble_loop! Ensemble of  particle
      nullify(deuterium_pertOrigin)

    end subroutine PertDecay

  end function oneBody


  !*****************************************************************************
  !****s* collisionTerm/twoBody
  ! NAME
  ! subroutine twoBody (teilchenPert, teilchenReal, time)
  ! PURPOSE
  ! Administrates the 2-body processes.
  ! The collision criteria is based upon the Kodama criteria.
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: teilchenPert -- perturbative particles
  ! * type(particle),dimension(:,:) :: teilchenReal -- real particles 
  ! * real                          :: time         -- actual time step
  !
  ! OUTPUT
  ! * teilchenPert and teilchenReal are changed
  !*****************************************************************************
  subroutine twoBody (teilchenPert, teilchenReal, time)
    
    use master_2Body, only: collide_2body
    use twoBodyStatistics, only : sqrts_distribution, rate
    use inputGeneral, only : fullEnsemble,eventtype,numEnsembles
    use eventtypes, only: HeavyIon
    use particleDefinition
    use output, only : WriteParticleVector, timeMeasurement
    use IdTable, only: isMeson, isBaryon
    use masterNBody, only: NBody_Analysis, check_for_Nbody, generate_3body_collision
    use XsectionRatios, only : accept_event
    use deuterium_PL, only : deuterium_pertOrigin,deuterium_pertOrigin_flag,deuteriumPL_ensemble
    use pauliBlockingModule, only: checkPauli
    use initHadron, only : particleId,antiparticle,particleCharge
    use insertion, only: particlePropagated, setIntoVector
    use CallStack, only: Traceback
    
    type(particle), intent(inout), target, dimension(:,:) :: teilchenPert, teilchenReal
    real, intent(in) :: time

    integer :: index1,index2,ensemble1,ensemble2,index1_start,index2_end,ensemble2_start,ensemble2_end  ! Loop variables
    type(particle), dimension(1:maxOut) :: finalState       ! final state of particles in 2-body collisions
    type(particle), dimension(1:maxOut) :: finalState_3Body ! final state of particles in 3-body collisions
    logical :: Flag, Flag_3Body, setFlag, HiEnergyFlag, HiEnergyFlag_3Body, flag1, flag2, flag3
    integer :: size_pert_dim2, size_real_dim2
    integer :: HiEnergyType, HiEnergyType_3Body
    integer :: i,number
    integer :: justdeletedelta_count=0,justdeletedelta_countO=0

    size_real_dim2=size(teilchenReal,dim=2)
    size_pert_dim2=size(teilchenPert,dim=2)

    if (twoBodyProcessesRealReal) call realReal
    if (twoBodyProcessesRealPert) call realPert
    
  contains

    !************************************
    ! real<->real collisions
    !************************************
    subroutine realReal

      use idTable, only: photon, nucleon, Lambda, SigmaResonance, Xi
      use collisionNumbering, only: check_justCollided, real_numbering, real_firstnumbering, ReportEventNumber

      type(particle), pointer :: part1, part2, part3
      type(particle) :: partIn(1:2)

      !***  Variables needed for many-body collisions  *************************
      integer :: n_found ! number of particles found around colliding pair
      integer, parameter :: Nmax=400 ! max possible value of n_found
      integer :: ind_min ! index of the most close particle to the collision point
      integer, dimension(1:Nmax) :: ind_found ! indices of found particles
      real :: coll_time  ! collision time instant
      real :: sigma,gamma ! cross section (mbarn) and gamma-factor in CM frame of colliding particles
      !*************************************************************************
      logical :: numbersAlreadySet
      integer :: imax
      logical :: pauliIsIncluded
      
      deuterium_pertOrigin_flag=99
      nullify(deuterium_pertOrigin)

      If(debug) write(*,*) 'real<->real'
      If(debug) call timeMeasurement(.true.)


      ensemble1_loop : do ensemble1=1,numEnsembles             !loop over first ensemble

         If (fullEnsemble) then
           ensemble2_End=numEnsembles
         else
           ensemble2_End=ensemble1
         end if

         ensemble2_loop : do ensemble2=ensemble1,ensemble2_End  ! loop over second ensemble
 
            if(ensemble2==ensemble1) then
              index1_start=2
            else
              index1_start=1
            end if
            
            deuteriumPL_ensemble=ensemble2
            
            index1_loop : do index1=index1_start,size_real_dim2

              If (teilchenReal(ensemble1,index1)%ID < 0) exit index1_loop !!! ID==-1
              If (teilchenReal(ensemble1,index1)%ID <= 0) cycle index1_loop
              If (teilchenReal(ensemble1,index1)%ID.eq.photon) cycle index1_loop

              part1 => teilchenReal(ensemble1,index1)
              partIn(1) = part1

              If (ensemble2.eq.ensemble1) then 
                index2_end=index1-1
              else
                index2_end=size_real_dim2
              end if

              if(DoJustAbsorptive .and. part1%event(1).ge.1000000) cycle index1_loop

              index2_loop : do index2=1,index2_end

                If (teilchenReal(ensemble2,index2)%ID < 0) exit index2_loop !!! ID==-1
                If (teilchenReal(ensemble2,index2)%ID <= 0) cycle index2_loop
                If (teilchenReal(ensemble2,index2)%ID.eq.photon) cycle index2_loop

                part2 => teilchenReal(ensemble2,index2)
                partIn(2) = part2

                if(DoJustAbsorptive .and. part2%event(1).ge.1000000) cycle index2_loop

                ! Check that they didn't collide just before:
                If (check_justCollided(part1,part2)) cycle index2_loop

                call resetNumberGuess
           
                call collide_2body ((/part1,part2/), finalState, time, flag, HiEnergyFlag, HiEnergyType, &
                                    collision_time=coll_time, PauliIncluded_out=pauliIsIncluded)

                coll_time=time+coll_time

                if (.not.flag) cycle index2_loop

                if (annihilate .and. HiEnergyType.eq.-3) cycle index2_loop

                if (finalState(3)%ID==0 .and. isElastic((/part1,part2/),finalState(1:2))) then
                  ! elastic
                  call sqrts_distribution((/part1,part2/),4)
                else
                  ! inelastic
                  call sqrts_distribution((/part1,part2/),5)
                end if

                flag_3Body=.false.
                part3=>teilchenReal(ensemble1,1)  ! This is just dummy setting to avoid stop when 
                                                  ! twoPlusOneBodyProcesses=.false.

              ! Check for possible presence of another particles nearby
              ! (currently only for heavy ion collisions, for the parallel 
              !  ensemble mode and only for the baryons):
                If (twoPlusOneBodyProcesses .and. eventtype==HeavyIon .and. .not.fullEnsemble &
                    .and. (isBaryon(part1%ID) .or. isBaryon(part2%ID))) then

                  if(debug) write(*,*)'before check_for_Nbody'

                  call check_for_Nbody((/index1,index2/),time,coll_time,&
                  & teilchenReal(ensemble1,:),n_found,ind_found,ind_min,sigma,gamma)

                  if(debug) write(*,*)'after check_for_Nbody:', n_found 

                  call Nbody_analysis((/index1,index2/),time,&
                  & teilchenReal(ensemble1,:),n_found,ind_found,sigma,gamma)

                  if(ind_min.ne.0) then
                    ! Fill histograms for stat. analysis:
                    call sqrts_distribution((/part1,part2/),2)
                    part3=>teilchenReal(ensemble1,ind_min)

                    call resetNumberGuess

                    call generate_3body_collision((/part1,part2,part3/),&
                    & finalState_3Body,flag_3Body,time,HiEnergyFlag_3Body,HiEnergyType_3Body)
                    if( .not.flag_3Body ) cycle index2_loop
                  else
                    call sqrts_distribution((/part1,part2/),1)                    
                  end if
             
                end If  ! End of many-body check

                if(flag_3Body) then
                  finalState=finalState_3Body
                  HiEnergyFlag=HiEnergyFlag_3Body
                  HiEnergyType=HiEnergyType_3Body
                end if

                if(debug) write(*,*) 'flag_3Body:', flag_3Body

                imax=0
                do i=1,maxout
                  if(finalState(i)%id <= 0) exit
                  imax=i
                  if(debug) then
                    write(*,*) 'Id, antiparticle ?, charge, number:',&
                               &finalState(i)%Id,finalState(i)%antiparticle,&
                               &finalState(i)%charge,finalState(i)%number
                  end if
                end do

                if (debug .and. flag_3Body) then
                  write(*,*)' 3-d particle:'
                  write(*,*) part3%Id, part3%antiparticle, part3%charge, part3%mass
                  write(*,*) part3%momentum
                  write(*,*) part3%velocity
                  write(*,*)' last final particle:'
                  write(*,*) finalState(imax)%Id, finalState(imax)%antiparticle, finalState(imax)%charge, finalState(imax)%mass
                  write(*,*) finalState(imax)%momentum
                  write(*,*) finalState(imax)%velocity
                end if

                if (.not. flag_3Body) then
                   ! Accept or reject the event randomly, using the in-medium cross section ratios:
                   if (.not.accept_event((/part1,part2/),finalState)) cycle index2_loop
                end if

                if (.not. pauliIsIncluded) call checkPauli(finalState,teilchenReal,flag)

                if (.not. flag) cycle index2_loop

                if (flag_3Body) then
                  if (debug) write(*,*)' In realreal, HiEnergyFlag=', HiEnergyFlag
                  flag = finalCheck ((/part1,part2,part3/), finalState, HiEnergyFlag)
                else
                  flag = finalCheck ((/part1,part2/), finalState, HiEnergyFlag)
                end if

                if (.not. flag) cycle index2_loop

                call rate((/part1,part2/),finalState,time)  ! Compute various collision rates

              ! (1) Label event by eventNumber
                number=real_numbering()
                if (DoJustAbsorptive) then
                   if (part1%id==particleId .and. (part1%antiparticle.eqv.antiparticle) .and. &
                       part1%charge==particleCharge) then
                     part1%event=number 
                     cycle index1_loop
                   else if (part2%id==particleId .and. (part2%antiparticle.eqv.antiparticle) .and. &
                            part2%charge==particleCharge) then
                     part2%event=number 
                     cycle index2_loop
                   end if
                end if
                finalState%event(1)=number
                finalState%event(2)=number

                !*** Removing Lambda/Sigma^0 history change if they were also in initial state:
                if ((part1%id==Lambda .or. (part1%id==SigmaResonance .and. part1%charge==0)) .and. &
                    .not.part1%antiparticle .and. part2%id==nucleon .and. .not.part2%antiparticle) then
                   do i=1,imax
                      if (finalState(i)%id==Lambda .or. (finalState(i)%id==SigmaResonance .and. finalState(i)%charge==0) .and. &
                          .not. finalState(i)%antiparticle) finalState(i)%history=part1%history
                   end do
                else if ((part2%id==Lambda .or. (part2%id==SigmaResonance .and. part2%charge==0)) .and. &
                         .not.part2%antiparticle .and. part1%id==nucleon .and. .not.part1%antiparticle) then
                   do i=1,imax
                      if (finalState(i)%id==Lambda .or. (finalState(i)%id==SigmaResonance .and. finalState(i)%charge==0) .and. &
                          .not. finalState(i)%antiparticle) finalState(i)%history=part2%history
                   end do
                else if (part1%id==Xi      .and. .not.part1%antiparticle .and. &
                         part2%id==nucleon .and. .not.part2%antiparticle) then
                   do i=1,imax
                      if (finalState(i)%id==Xi .and. .not. finalState(i)%antiparticle) finalState(i)%history=part1%history
                   end do
                else if (part2%id==Xi      .and. .not.part2%antiparticle .and. &
                         part1%id==nucleon .and. .not.part1%antiparticle) then
                   do i=1,imax
                      if(finalState(i)%id==Xi .and. .not. finalState(i)%antiparticle) finalState(i)%history=part2%history
                   end do
                end if

                finalState%perweight = min(part1%perweight, part2%perweight)
                number = min(part1%firstEvent, part2%firstEvent)
                if (number==0) number = real_firstnumbering()
                finalState%firstEvent = number
                finalState%perturbative = .false.
                
                ! (2) Create final state particles
                NumbersAlreadySet = AcceptGuessedNumbers()

                flag1=.false.
                flag2=.false.
                flag3=.false.
                do i=1,imax
                    If (particlePropagated(finalState(i))) then
                      if (.not.flag1 .and. ((isBaryon(part1%ID) .and. isBaryon(finalState(i)%ID) .and. &
                                             (part1%antiparticle.eqv.finalState(i)%antiparticle)) .or. &
                                            (isMeson(part1%ID).and.isMeson(finalState(i)%ID)))) then
                        if (.not.NumbersAlreadySet .or. finalState(i)%number==0) call setNumber(finalState(i))
                        part1 = finalState(i)
                        flag1 = .true.
                      else if (.not.flag2 .and. ((isBaryon(part2%ID) .and. isBaryon(finalState(i)%ID) .and. &
                                                  (part2%antiparticle.eqv.finalState(i)%antiparticle)) .or. &
                                                 (isMeson(part2%ID) .and. isMeson(finalState(i)%ID)))) then
                        if (.not.NumbersAlreadySet .or. finalState(i)%number==0) call setNumber(finalState(i))
                        part2 = finalState(i)
                        flag2 = .true.
                      else if (flag_3Body .and. .not.flag3 .and. ((isBaryon(part3%ID) .and. isBaryon(finalState(i)%ID) .and. &
                                                                   (part3%antiparticle.eqv.finalState(i)%antiparticle)) .or. &
                                                                  (isMeson(part3%ID) .and. isMeson(finalState(i)%ID)))) then
                        if (.not.NumbersAlreadySet .or. finalState(i)%number==0) call setNumber(finalState(i))
                        part3 = finalState(i)
                        flag3 = .true.
                      else
                      ! Find empty space in the particle vector: 
                        If (fullensemble) then
                          call setIntoVector (finalState(i:i), teilchenReal, setFlag, NumbersAlreadySet)
                        else
                          call setIntoVector (finalState(i:i), teilchenReal(ensemble1:ensemble1,:), setFlag, NumbersAlreadySet)
                        end if
                      ! (3) Check that setting into real particle
                      !     vector worked out :
                        If (.not.setFlag) then
                          write(*,*) 'Real particle vector too small!'
                          write(*,*) size(finalState),&
                          lBound(teilchenReal), uBound(teilchenReal)
                          write(*,*) 'Dumping real particles to files "RealParticles_stop_*!"'
                          call WriteParticleVector("RealParticles_stop",teilchenReal)
                          write(*,*) 'Dumping perturbative particles to files "PertParticles_stop_*!"'
                          call WriteParticleVector("PertParticles_stop",teilchenPert)
                          call Traceback('collision Term, two-body')
                        end if
                      end if
                    end If
                end do

                ! For statistical purposes and also for possible calculation of the time step
                ! also a 3-body collision is treated here as a 2-body collision:
                call ReportEventNumber (partIn, finalState, finalState(1)%event, time, 211, HiEnergyType)


                ! (4) destroy incoming particles if there are
                ! no corresponding final leading particles:
                if(.not.flag2) part2%Id=0
                if(flag_3Body .and. .not.flag3) part3%Id=0
                if(.not.flag1) then
                  part1%Id=0
                  cycle index1_loop
                endif

              end do index2_loop ! Index2 of second particle
            end do index1_loop ! Index1 of first particle
         end do ensemble2_loop ! Ensemble2 of second particle
      end do ensemble1_loop ! Ensemble1 of first particle


      If (twoPlusOneBodyProcesses .and. eventtype==HeavyIon .and. .not.fullEnsemble) &
        call Nbody_analysis((/index1,index2/),time,teilchenReal(ensemble1,:),n_found,ind_found,sigma,gamma,.true.)


      deuterium_pertOrigin_flag=0

    end subroutine realReal
    
    !************************************
    ! real<->pert collisions
    !************************************
    subroutine realPert
      use IdTable, only: photon, nucleon, Delta, pion
      use statistics, only : saveInfo
      use CollHistory, only: CollHist_UpdateHist
      use collisionNumbering, only: check_justCollided, pert_numbering, pert_firstnumbering, ReportEventNumber

      type(particle),pointer :: pair1, pair2
      real :: EventWeight
      logical :: numbersAlreadySet
      integer :: number2
      logical :: pauliIsIncluded
      integer, dimension(2,maxOut) :: posOut
      type(particle) :: pairS2

      If(debug) write(*,*) 'real<->pert'
      If(debug) call timeMeasurement(.false.)
      If(debug) call timeMeasurement(.true.)
      
      ensemble1_loop : do ensemble1=1,numEnsembles   ! loop over first particle = real particle
         If (fullEnsemble) then
            ensemble2_Start=1
            ensemble2_End=numEnsembles
         else
            ensemble2_Start=ensemble1
            ensemble2_End=ensemble1
         end if
         if(debug) write(*,*)  'Loop über',ensemble1, ensemble2_Start,   ensemble2_End

         index1_loop: do index1=1,size_real_dim2
            If (teilchenReal(ensemble1,index1)%ID < 0) exit index1_loop
            If (teilchenReal(ensemble1,index1)%ID <= 0) cycle index1_loop
            If (teilchenReal(ensemble1,index1)%ID.eq.photon) cycle index1_loop
            pair1=>teilchenReal(ensemble1,index1)
            NULLIFY(deuterium_pertOrigin)
            deuterium_pertOrigin=>teilchenReal(ensemble1,index1)
            deuterium_pertOrigin_flag=2
            ensemble2_loop : do ensemble2=ensemble2_Start,ensemble2_End  
               ! loop over second particle = perturbative
               index2_loop : do index2=1,size_Pert_dim2
                  If (teilchenPert(ensemble2,index2)%ID < 0) exit index2_loop
                  If (teilchenPert(ensemble2,index2)%ID <= 0) cycle index2_loop
                  If (teilchenPert(ensemble2,index2)%ID.eq.photon) cycle index2_loop
                  pair2 => teilchenPert(ensemble2,index2)
                  pairS2 = pair2
                  
                  ! Check that they didn't collide just before
                  If (check_justCollided(pair1,pair2)) cycle index2_loop

                  ! No Nuk-Nuk collisions
                  If (noNucNuc .and. (pair1%ID==nucleon .and. pair2%ID==nucleon)) cycle index2_loop

                  call resetNumberGuess

                  call collide_2body ((/pair1,pair2/), finalState, time, flag, HiEnergyFlag, HiEnergyType, &
                                      EventWeight, PauliIncluded_out=pauliIsIncluded)
                  if (.not.flag) cycle index2_loop
!                  If(debug) write(*,*) 'final states. IDs=',finalState%ID,'  charges=',finalState%Charge
!                  If(debug) write(*,*) 'HiEnergy: ',HiEnergyFlag,HiEnergyType
                  
                  flag = finalCheck ((/pair1,pair2/), finalState, HiEnergyFlag) 
                  if (.not.flag) cycle index2_loop
                  
                  if (.not.pauliIsIncluded) call checkPauli(finalState,teilchenReal,flag)
                  if (.not.flag) cycle index2_loop

!                  write(*,*) 'colliding:',pair1%number,pair2%number,HiEnergyType
                  
                  NumbersAlreadySet = AcceptGuessedNumbers()

                  if (DoJustAbsorptive) finalState%ID = 0

                  !simulate Effenberger treatment of Deltas
                  if(justDeleteDelta) then
                     if(((pair1%ID.eq.nucleon) .and. (pair2%ID.eq.delta)).or.((pair2%ID.eq.nucleon) .and. (pair1%ID.eq.delta))) then
                        !check that we have N N in finalstate
                        justdeletedelta_count=0
                        justdeletedelta_countO=0
                        do i=1,maxout
                           if(finalState(i)%id <= 0) exit
                           if(finalstate(i)%id.eq.1) justdeletedelta_count=justdeletedelta_count+1
                           if(finalstate(i)%id.ne.1) justdeletedelta_countO=justdeletedelta_countO+1
                        end do
                        if(justdeletedelta_count.eq.2.and.justdeletedelta_countO.eq.0) finalstate%id=0
                     end if
                  end if

                  
                  ! (1) Label event by eventNumber, such that we can track the perturbative particle to its 
                  ! real scattering partner.
                  
                  number = pert_numbering(pair1)
                  if (pair2%firstEvent==0) then
                     number2=pert_firstnumbering(pair1,pair2)
                  else
                     number2=pair2%firstEvent
                  endif

                  do i=1,maxout
                     if(finalState(i)%id <= 0) exit
                     
                     finalState(i)%event = number
                     finalState(i)%lastCollisionTime = time
                     finalState(i)%perweight = pair2%perweight * EventWeight
                     finalState(i)%firstEvent = number2
                     finalState(i)%perturbative = .true.
                     if (.not.NumbersAlreadySet) call setNumber(finalState(i))

                     !###
                     ! Writing statistical information
                     if (useStatistics .and. finalState(i)%ID==pion) then
                       call saveInfo(finalState(i)%number,ensemble2,finalstate(i)%position, pair1%ID,pair2%ID,0)
                       if (debugSaveInfo) write(99,*) finalState(i)%number, ensemble2,finalstate(i)%position, pair1%ID,pair2%ID,0
                     end if
                     if (printPositions .and. (finalState(i)%id==pion .or. finalState(i)%id==nucleon)) &
                       call positionPrinter (finalstate(i))
                     !####
                  end do
                  
                  ! (2) Eliminate incoming perturbative particles
                  teilchenPert(ensemble2,index2)%Id=0

                  ! (3) Create final state particles
                  posOut = 0
                  If (particlePropagated(finalState(1))) then
                     teilchenPert(ensemble2,index2) = finalState(1)
                     posOut(1:2,1) = (/ensemble2,index2/)
                  end if

                  If (fullensemble) then
                    call setIntoVector(finalState(2:),teilchenPert,setFlag,.true.,positions=posOut(:,2:))
                  else
                    call setIntoVector(finalState(2:),teilchenPert(ensemble1:ensemble1,:),setFlag,.true.,positions=posOut(:,2:))
                    posOut(1,2:) = ensemble1
                  end if

                  call ReportEventNumber ((/pair1,pairS2/), finalState, finalstate(1)%event, time, 212, HiEnergyType)

                  call CollHist_UpdateHist ((/pair1,pairS2/), finalState, (/ensemble2,index2/), posOut, finalState(1)%perweight)

                  ! (4) Check that setting into perturbative particle vector worked out 
                  If (.not.setFlag) then
                     write(*,*) 'Perturbative particle vector too small!'
                     write(*,*) size(finalState),  lBound(teilchenPert), uBound(teilchenPert)
                     write(*,*) 'Dumping real particles to files "RealParticles_stop_*!"'
                     call WriteParticleVector("RealParticles_stop",teilchenReal)
                     write(*,*) 'Dumping perturbative particles to files "PertParticles_stop_*!"'
                     call WriteParticleVector("PertParticles_stop",teilchenPert)
                     call Traceback('collision Term, two-body')
                  end if
               end do index2_loop ! Index2 of second particle
            end do ensemble2_loop ! Ensemble2 of second particle
         end do index1_loop! Index1 of first particle
      end do ensemble1_loop! Ensemble1 of first particle
      if (printPositions) call positionPrinter()

      If(debug) call timeMeasurement(.false.)
      nullify(deuterium_pertOrigin)
      deuterium_pertOrigin_flag=999


    end subroutine realPert
    
  end subroutine twoBody


  !*****************************************************************************
  !****s* collisionTerm/twoBody_local
  ! NAME
  ! subroutine twoBody_local (teilchenPert, teilchenReal, time)
  !
  ! PURPOSE
  ! * Administrates the 2-body processes.
  ! * The collision criteria is based upon the local collision criteria. 
  !   Therefore the volume elements (module VolumeElements) must be 
  !   initialized. Only scatterings within one volume element are allowed.
  !
  ! NOTES 
  ! For real-real collisions there is still the Kodama procedure implemented. 
  ! This shall be changed in the future.
  !
  ! INPUTS
  ! * type(particle),intent(inOUT),target,dimension(:,:) :: teilchenReal -- real particles 
  ! * type(particle),intent(inOUT),target,dimension(:,:) :: teilchenPert -- perturbative particles 
  ! * real, intent(in) :: time
  !
  ! OUTPUT
  ! * teilchenPert and teilchenReal are changed
  !*****************************************************************************
  subroutine twoBody_local (teilchenPert, teilchenReal, time)
    
    use master_2Body, only: collide_2body
    use inputGeneral, only : fullEnsemble,eventtype,numEnsembles
    use eventtypes, only: HeavyIon
    use particleDefinition
    use output, only : WriteParticleVector, timeMeasurement
    use IdTable, only: nucleon, delta, pion, photon, isMeson, isBaryon
    use masterNBody, only: NBody_Analysis, check_for_Nbody, generate_3body_collision
    use VolumeElements, only: VolumeElements_CLEAR, VolumeElements_SETUP_Pert, VolumeElements_SETUP_Real
    use pauliBlockingModule, only: checkPauli
    use insertion, only: particlePropagated, setIntoVector
    use CallStack, only: Traceback

    type(particle), intent(inout), dimension(:,:) :: teilchenPert
    type(particle), intent(inout), dimension(:,:), target :: teilchenReal
    real, intent(in) :: time

    integer :: index1,index2,ensemble1,ensemble2,index1_start,index2_end,ensemble2_end!,ensemble2_start
    type(particle), dimension(1:maxOut) :: finalState ! final state of particles
    logical :: Flag, setFlag, HiEnergyFlag, flag1, flag2
    integer :: size_real_dim2!,size_pert_dim2
    integer :: HiEnergyType,i,number
    real :: XS
    integer, save :: justdeletedelta_count=0, justdeletedelta_countO=0

    size_real_dim2=size(teilchenReal,dim=2)
    !size_pert_dim2=size(teilchenPert,dim=2)

    if (twoBodyProcessesRealReal) call realReal

    !=== Setup VolumeElements

    call VolumeElements_CLEAR
    call VolumeElements_SETUP_Pert(teilchenPert)
    call VolumeElements_SETUP_Real(teilchenReal)
!    call VolumeElements_Statistics

    if (twoBodyProcessesRealPert) call realPert
    
  contains

    !************************************
    ! real<->real collisions
    !************************************

    subroutine realReal

      use collisionNumbering, only: check_justCollided, real_numbering, real_firstnumbering, ReportEventNumber

      type(particle), pointer :: part1, part2, part3
      type(particle) :: partIn(1:2)

      !***  Variables needed for many-body collisions  ***************************
      integer :: n_found ! number of particles found around colliding part
      integer, parameter :: Nmax=60 ! max possible value of n_found
      integer :: ind_min ! index of the most close particle to the collision point
      integer, dimension(1:Nmax) :: ind_found ! indices of found particles
      real :: coll_time  ! collision time instant
      real :: sigma,gamma ! cross section (mbarn) and gamma-factor in CM frame of colliding particles
      !****************************************************************************
      type(particle), dimension(1:maxOut) :: finalState_3Body ! final state of particles in 3-body collisions
      logical :: Flag_3Body, HiEnergyFlag_3Body, flag3, numbersAlreadySet, pauliIsIncluded
      integer :: HiEnergyType_3Body, imax

      If(debug) write(*,*) 'real<->real'
      If(debug) call timeMeasurement(.true.)

      ensemble1_loop : do ensemble1=1,numEnsembles             !loop over first ensemble

         If (fullEnsemble) then
           ensemble2_End=numEnsembles
         else
           ensemble2_End=ensemble1
         end if

         ensemble2_loop : do ensemble2=ensemble1,ensemble2_End  ! loop over second ensemble
 
            if(ensemble2==ensemble1) then
              index1_start=2
            else
              index1_start=1
            end if
            
            index1_loop : do index1=index1_start,size_real_dim2

              If (teilchenReal(ensemble1,index1)%ID < 0) exit index1_loop !!! ID==-1
              If (teilchenReal(ensemble1,index1)%ID <= 0) cycle index1_loop
              If (teilchenReal(ensemble1,index1)%ID.eq.photon) cycle index1_loop
            
              part1 => teilchenReal(ensemble1,index1)
              partIn(1) = part1

              If (ensemble2.eq.ensemble1) then 
                index2_end=index1-1
              else
                index2_end=size_real_dim2
              end if

              index2_loop : do index2=1,index2_end

                If (teilchenReal(ensemble2,index2)%ID < 0) exit index2_loop !!! ID==-1
                If (teilchenReal(ensemble2,index2)%ID <= 0) cycle index2_loop
                If (teilchenReal(ensemble2,index2)%ID.eq.photon) cycle index2_loop

                part2 => teilchenReal(ensemble2,index2)
                partIn(2) = part2

                ! Check that they didn't collide just before:
                If (check_justCollided(part1,part2)) cycle index2_loop

                if(debug) then
                  write(*,*)'before collide_2body:'
                  write(*,*)'Ids:', part1%Id,part2%Id
                  write(*,*)'Antiparticles ?:', part1%antiparticle,part2%antiparticle
                  write(*,*)'charges:',part1%charge,part2%charge
                  write(*,*)'masses:',part1%mass,part2%mass
                  write(*,*)'numbers:', part1%number,part2%number
                  write(*,*)'events:',part1%event,part2%event
                  write(*,*) '4-mom 1-st part:',part1%momentum
                  write(*,*) 'velo 1-st part:',part1%velocity
                  write(*,*) '4-mom 2-nd part:',part2%momentum
                  write(*,*) 'velo 2-nd part:',part2%velocity
                end if
           
                call resetNumberGuess

                call collide_2body((/part1,part2/), finalState, time, flag, &
                     HiEnergyFlag, HiEnergyType, &
                     collision_time=coll_time, PauliIncluded_out=pauliIsIncluded)

                coll_time=time+coll_time

                if (.not.flag) cycle index2_loop

                flag_3Body=.false.
                part3=>teilchenReal(ensemble1,1)  ! This is just dummy setting to avoid stop when 
                                                  ! twoPlusOneBodyProcesses=.false.

              ! Check for possible presence of another particles nearby
              ! (currently only for heavy ion collisions, for the parallel 
              !  ensemble mode and only for the baryons):
                If (twoPlusOneBodyProcesses .and. eventtype==HeavyIon .and. .not.fullEnsemble &
                    .and. isBaryon(part1%ID) .and. isBaryon(part2%ID)) then

                  if(debug) write(*,*)'before check_for_Nbody'

                  call check_for_Nbody((/index1,index2/),time,coll_time,&
                  & teilchenReal(ensemble1,:),n_found,ind_found,ind_min,sigma,gamma)

                  if(debug) write(*,*)'after check_for_Nbody:', n_found 

                  call Nbody_analysis((/index1,index2/),time,&
                  & teilchenReal(ensemble1,:),n_found,ind_found,sigma,gamma)

                  if(ind_min.ne.0) then
                    part3=>teilchenReal(ensemble1,ind_min)

                    call resetNumberGuess

                    call generate_3body_collision((/part1,part2,part3/),&
                    & finalState_3Body,flag_3Body,time,HiEnergyFlag_3Body,HiEnergyType_3Body)
                  end if
             
                end If  ! End of many-body check

                if(flag_3Body) then
                  finalState=finalState_3Body
                  HiEnergyFlag=HiEnergyFlag_3Body
                  HiEnergyType=HiEnergyType_3Body
                end if

                if(debug) write(*,*) 'flag_3Body:', flag_3Body

                do i=1,maxout
                  if(finalState(i)%id <= 0) exit
                  imax=i
                  if(debug) then
                    write(*,*) 'Id, antiparticle ?, charge, number:',&
                               &finalState(i)%Id,finalState(i)%antiparticle,&
                               &finalState(i)%charge,finalState(i)%number
                  end if
                end do

                if (debug .and. flag_3Body) then
                  write(*,*)' 3-d particle:'
                  write(*,*) part3%Id,part3%antiparticle,part3%charge,part3%mass
                  write(*,*) part3%momentum
                  write(*,*) part3%velocity
                  write(*,*)' last final particle:'
                  write(*,*) finalState(imax)%Id,finalState(imax)%antiparticle,finalState(imax)%charge,finalState(imax)%mass
                  write(*,*) finalState(imax)%momentum
                  write(*,*) finalState(imax)%velocity
                end if

                if (.not. pauliIsIncluded) call checkPauli(finalState,teilchenReal,flag)

                if (.not. flag) cycle index2_loop

                if (flag_3Body) then
                  if(debug) write(*,*)' In realreal, HiEnergyFlag=', HiEnergyFlag
                  flag = finalCheck ((/part1,part2,part3/), finalState, HiEnergyFlag)
                else
                  flag = finalCheck ((/part1,part2/), finalState, HiEnergyFlag)
                end if

                if (.not. flag) cycle index2_loop

                ! (1) Label event by eventNumber
                number=real_numbering()
                finalState%event(1)=number
                finalState%event(2)=number

                finalState%lastCollisionTime = time
                finalState%perturbative = .false.
                finalState%perweight = min(Part1%perweight, Part2%perweight)
                number = min(part1%firstEvent,part2%firstEvent)
                if (number==0) number = real_firstnumbering()
                finalState%firstEvent=number
      
                ! (2) Create final state particles
                NumbersAlreadySet = AcceptGuessedNumbers()

                flag1=.false.
                flag2=.false.
                flag3=.false.
                do i=1,imax
                    If (particlePropagated(finalState(i))) then
                      if (.not.flag1 .and. ((isBaryon(part1%ID) .and. isBaryon(finalState(i)%ID) .and. &
                                             (part1%antiparticle.eqv.finalState(i)%antiparticle)) .or. &
                                            (isMeson(part1%ID) .and. isMeson(finalState(i)%ID)))) then
                        if (.not.NumbersAlreadySet .or. finalState(i)%number==0) call setNumber(finalState(i))
                        part1 = finalState(i)
                        flag1 = .true.
                      else if (.not.flag2 .and. ((isBaryon(part2%ID) .and. isBaryon(finalState(i)%ID) .and. &
                                                  (part2%antiparticle.eqv.finalState(i)%antiparticle)) .or. &
                                                 (isMeson(part2%ID) .and. isMeson(finalState(i)%ID)))) then
                        if (.not.NumbersAlreadySet .or. finalState(i)%number==0) call setNumber(finalState(i))
                        part2 = finalState(i)
                        flag2 = .true.
                      else if (flag_3Body .and. .not.flag3 .and. ((isBaryon(part3%ID) .and. isBaryon(finalState(i)%ID) .and. &
                                                                   (part3%antiparticle.eqv.finalState(i)%antiparticle)) .or. &
                                                                  (isMeson(part3%ID) .and. isMeson(finalState(i)%ID)))) then
                        if (.not.NumbersAlreadySet .or. finalState(i)%number==0) call setNumber(finalState(i))
                        part3 = finalState(i)
                        flag3 = .true.
                      else
                      ! Find empty space in the particle vector: 
                        If (fullensemble) then
                          call setIntoVector(finalState(i:i), teilchenReal, setFlag, NumbersAlreadySet)
                        else
                          call setIntoVector(finalState(i:i), teilchenReal(ensemble1:ensemble1,:), setFlag, NumbersAlreadySet)
                        end if
                      ! (3) Check that setting into real particle
                      !     vector worked out :
                        If (.not. setFlag) then
                          write(*,*) 'Real particle vector too small!'
                          write(*,*) size(finalState),&
                          lBound(teilchenReal), uBound(teilchenReal)
                          write(*,*) 'Dumping real particles to files "RealParticles_stop_*!"'
                          call WriteParticleVector("RealParticles_stop",teilchenReal)
                          write(*,*) 'Dumping perturbative particles to files "PertParticles_stop_*!"'
                          call WriteParticleVector("PertParticles_stop",teilchenPert)
                          call Traceback('collision Term, two-body')
                        end if
                      end if
                    end If
                end do

              ! For statistical purposes and also for possible calculation of the time step
              ! also a 3-body collision is treated here as a 2-body collision:
                call ReportEventNumber (partIn, finalState, finalState(1)%event, time, 211, HiEnergyType)


              ! (4) destroy incoming particles if there are
              ! no corresponding final leading particles:
                if(.not.flag2) part2%Id=0
                if(flag_3Body .and. .not.flag3) part3%Id=0
                if(.not.flag1) then
                  part1%Id=0
                  cycle index1_loop
                endif

              end do index2_loop ! Index2 of second particle
            end do index1_loop ! Index1 of first particle
         end do ensemble2_loop ! Ensemble2 of second particle
      end do ensemble1_loop ! Ensemble1 of first particle

      If (twoPlusOneBodyProcesses .and. eventtype==HeavyIon .and. .not.fullEnsemble) then
         call Nbody_analysis((/index1,index2/),time,&
              teilchenReal(ensemble1,:),n_found,ind_found,sigma,gamma,.true.)
      end If
    end subroutine realReal
    
    !************************************
    ! real<->pert collisions
    !************************************
    subroutine realPert
      use statistics, only : saveInfo
      use VolumeElements, only: VolumeElements_InitGetPart, VolumeElements_GetPart
      use CollHistory, only: CollHist_UpdateHist
      use collisionNumbering, only: check_justCollided, pert_numbering, pert_firstnumbering, ReportEventNumber

      type(particle),pointer :: pair1, pair2
      real :: EventWeight
      logical :: numbersAlreadySet, pauliIsIncluded
      integer :: nRealPart, number2
      !real :: dummy
      integer, dimension(2,maxOut) :: posOut
      type(particle) :: pairS2

      If (debug) then
        write(*,*) 'real<->pert'
        call timeMeasurement(.false.)
        call timeMeasurement(.true.)
      end if

      call VolumeElements_InitGetPart

      do 
         if (.not. VolumeElements_GetPart(pair1, pair2, nRealPart,ensemble2,index2)) exit

         ! Check that they didn't collide just before
         If (check_justCollided(pair1,pair2)) cycle

         ! No Nuk-Nuk collisions
         If (noNucNuc .and. (pair1%ID==nucleon) .and. (pair2%ID==nucleon)) cycle

         ! No Photon-Particle Collisions
         if ((pair1%ID==photon) .or. (pair2%ID==photon)) cycle


         call resetNumberGuess

         call collide_2body((/pair1,pair2/), finalState, time, flag, &
              HiEnergyFlag, HiEnergyType, &
              EventWeight, weightLocal=nRealPart, sigTot_out=XS, PauliIncluded_out=pauliIsIncluded)
         ! weightLocal=Number of real particles in volume, integer

         if (.not.flag) cycle
!         If(debug) write(*,*) 'final states. IDs=',finalState%ID,'  charges=',finalState%Charge
!         If(debug) write(*,*) 'HiEnergy: ',HiEnergyFlag,HiEnergyType
                  
         flag = finalCheck ((/pair1,pair2/), finalState, HiEnergyFlag) 
         if (.not.flag) cycle
                  
         if(.not.pauliIsIncluded) call checkPauli(finalState,teilchenReal,flag)

!!$         if (pair2%ID.eq.101) then
!!$            dummy = sqrtS((/pair1,pair2/))
!!$
!!$            if (flag) then
!!$               call AddHist(hXStot,dummy,1.0,1.0)
!!$               if (HiEnergyType.eq.-1) call AddHist(hXSelast,dummy,1.0,1.0)
!!$            else
!!$               call AddHist(hXStot,dummy,1.0,0.0)
!!$               if (HiEnergyType.eq.-1) call AddHist(hXSelast,dummy,1.0,0.0)
!!$            endif
!!$!            write(*,*) pair2%charge,pair1%charge, HiEnergyType,flag,dummy
!!$         endif

         if (.not.flag) cycle

!         write(*,*) 'colliding:',pair1%number,pair2%number,HiEnergyType
                  
         NumbersAlreadySet = AcceptGuessedNumbers()

         if (DoJustAbsorptive) finalState%ID = 0

         !simulate Effenberger treatment of Deltas
         if(justDeleteDelta) then
            if(((pair1%ID.eq.nucleon) .and. (pair2%ID.eq.delta)).or.((pair2%ID.eq.nucleon) .and. (pair1%ID.eq.delta))) then
               !check that we have N N in finalstate
               justdeletedelta_count=0
               justdeletedelta_countO=0
               do i=1,maxout
                  if(finalState(i)%id <= 0) exit
                  if(finalstate(i)%id.eq.1) justdeletedelta_count=justdeletedelta_count+1
                  if(finalstate(i)%id.ne.1) justdeletedelta_countO=justdeletedelta_countO+1
               end do
               if(justdeletedelta_count.eq.2.and.justdeletedelta_countO.eq.0) finalstate%id=0
            end if
         end if
                  
         ! (1) Label event by eventNumber, such that we can track the perturbative particle to its 
         ! real scattering partner.
                  
         number=pert_numbering(pair1)
         if (pair2%firstEvent==0) then
            number2=pert_firstnumbering(pair1,pair2)
         else
            number2=pair2%firstEvent
         endif

         do i=1,maxout
            if(finalState(i)%id <= 0) exit
            
            finalState(i)%event = number
            finalState(i)%lastCollisionTime = time
            finalState(i)%perweight = pair2%perweight * EventWeight
            finalState(i)%firstEvent = number2
            finalState(i)%perturbative = .true.
            if (.not.NumbersAlreadySet) call setNumber(finalState(i))
            
            !###
            ! Writing statistical information
            if (useStatistics .and. finalState(i)%ID==pion) then
              call saveInfo(finalState(i)%number,ensemble2,finalstate(i)%position, pair1%ID,pair2%ID,0)
              if(debugSaveInfo) write(99,*) finalState(i)%number, ensemble2,finalstate(i)%position, pair1%ID,pair2%ID,0
            end if
            if (printPositions .and. (finalState(i)%id==pion .or. finalState(i)%id==nucleon)) &
              call positionPrinter(finalstate(i))
            !####
         end do
                  
         ! (2) Eliminate incoming perturbative particles
         pairS2 = pair2
         pair2%Id = 0

         ! (3) Create final state particles
         posOut = 0
         If (particlePropagated(finalState(1))) then
            pair2 = finalState(1)
            posOut(1:2,1) = (/ensemble2,index2/)
         end if

         call setIntoVector(finalState(2:),teilchenPert,setFlag,.true.,positions=posOut(:,2:))

         call ReportEventNumber ((/pair1,pairS2/), finalState, finalstate(1)%event, time, 212, HiEnergyType)

         call CollHist_UpdateHist ((/pair1,pairS2/), finalState, (/ensemble2,index2/), posOut, finalState(1)%perweight)

         ! (4) Check that setting into perturbative particle vector worked out 
         If (.not.setFlag) then
            write(*,*) 'Perturbative particle vector too small!'
            write(*,*) size(finalState),  lBound(teilchenPert), uBound(teilchenPert)
            write(*,*) 'Dumping real particles to files "RealParticles_stop_*!"'
            call WriteParticleVector("RealParticles_stop",teilchenReal)
            write(*,*) 'Dumping perturbative particles to files "PertParticles_stop_*!"'
            call WriteParticleVector("PertParticles_stop",teilchenPert)
            call Traceback('collision Term, two-body')
         end if

      end do
      
      if (printPositions) call positionPrinter()
      
      If (debug) call timeMeasurement(.false.)
    end subroutine realPert
    
  end subroutine twoBody_local


  !*****************************************************************************
  !****s* collisionTerm/positionPrinter
  ! NAME
  ! subroutine positionPrinter (p)
  ! PURPOSE
  ! This routine saves the position of a given particle into a 2D Histogram as a function
  ! of cylinder coordinates z and rho. Each new entry is added to the previous ones.
  ! The histograms are meant to be used as tables of all production points.
  !
  ! If no input is given, the histogrames are printed to file.
  !
  ! INPUTS
  ! * type(particle), intent(in), optional :: p
  !
  ! OUTPUT
  ! * file 'ProdPlaces_pionPlus.dat'    -- All saved positions of pi^+ 
  ! * file 'ProdPlaces_pionNull.dat'    -- All saved positions of pi^0
  ! * file 'ProdPlaces_pionMinus.dat'   -- All saved positions of pi^-
  ! * file 'ProdPlaces_nucleon.dat'     -- All saved positions of nucleons
  !*****************************************************************************
  subroutine positionPrinter (p)
    use particleDefinition
    use hist2Df90
    use idTable, only: pion, nucleon

    type(particle), intent(in), optional :: p

    ! store position information in terms of z and the zylinder coordinate rho :
    type(histogram2D),save :: positions_pi(-1:1),positions_nuc

    logical, save :: initFlagge=.true.
    real, dimension(1:3) :: r
    real :: rho, weight

    if (initFlagge) then
       call createHist2d(positions_pi(1),  'ProdPlaces_pion_plus',  (/-8.,0./),(/8.,8./),(/0.25,0.25/))
       call createHist2d(positions_pi(0),  'ProdPlaces_pion_null',  (/-8.,0./),(/8.,8./),(/0.25,0.25/))
       call createHist2d(positions_pi(-1), 'ProdPlaces_pion_minus', (/-8.,0./),(/8.,8./),(/0.25,0.25/))
       call createHist2d(positions_nuc,    'ProdPlaces_nucleon',    (/-8.,0./),(/8.,8./),(/0.25,0.25/))
       initFlagge=.false.
    end if

    if (present(p)) then
       r = p%position
       rho = sqrt(r(1)**2+r(2)**2)
       weight = p%perweight
       select case (p%ID)
       case (pion)
         call AddHist2D(positions_pi(p%charge),(/r(3),rho/),weight)
       case (nucleon)
         call AddHist2D(positions_nuc,(/r(3),rho/),weight)
       end select
    else
       Open(101,file='ProdPlaces_pionPlus.dat')
       call WriteHist2D_Gnuplot(positions_pi(1),101)
       close (101)
       Open(101,file='ProdPlaces_pionNull.dat')
       call WriteHist2D_Gnuplot(positions_pi(0),101)
       close (101)
       Open(101,file='ProdPlaces_pionMinus.dat')
       call WriteHist2D_Gnuplot(positions_pi(-1),101)
       close (101)
       Open(101,file='ProdPlaces_nucleon.dat')
       call WriteHist2D_Gnuplot(positions_nuc,101)
       close (101)
    end if

  end subroutine positionPrinter


  !*****************************************************************************
  !****s* collisionTerm/threeBody
  ! NAME
  ! subroutine threeBody (teilchenPert, teilchenReal, time)
  ! PURPOSE
  ! Administrates the 3-body processes.
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: teilchenPert -- perturbative particles
  ! * type(particle),dimension(:,:) :: teilchenReal -- real particles 
  ! * real                          :: time         -- actual time step
  !
  ! OUTPUT
  ! * teilchenPert and teilchenReal are changed
  !*****************************************************************************
  subroutine threeBody (teilchenPert, teilchenReal, time)

    use IDTable, only: pion, nucleon, Delta
    use inputGeneral, only : fullEnsemble,localEnsemble
    use particleDefinition
    use master_3Body, only : GetRadiusNukSearch, nukSearch, make_3Body_Collision
    use collisionNumbering, only : pert_numbering,real_numbering, ReportEventNumber,real_firstnumbering
    use output, only : WriteParticleVector
    use VolumeElements, only : VolumeElements_NukSearch,VolumeElements_CLEAR,VolumeElements_SETUP_PERT, &
                               VolumeElements_SETUP_REAL !,VolumeElements_Statistics  
    use CollHistory, only: CollHist_UpdateHist
    use pauliBlockingModule, only: checkPauli
    use insertion, only: particlePropagated, setIntoVector
    use CallStack, only: Traceback
    
    type(particle), intent(inout), dimension(:,:) :: teilchenPert, teilchenReal
    real, intent(in) :: time
    
    INTEGER :: i,j,number,ii
    logical :: successFlag
    logical :: nukSearch_VE=.true.
    type(particle), pointer:: proton1, proton2                 ! Closest protons
    type(particle), pointer:: neutron1, neutron2               ! Closest neutrons
    type(particle), pointer:: scatterPartner1, scatterPartner2 ! Particles which one is scattering with
    type(particle), dimension(1:3) :: FinalState
   ! logical :: justDeleteDelta=.false. ! If Delta shall be just deleted - 
   !                                   ! therefore energy conservation is violated.
    integer, dimension(2,maxOut) :: posOut
    type(particle) :: pairS2
    
    if (localEnsemble) then
       !=== Setup VolumeElements
       call VolumeElements_CLEAR
       call VolumeElements_SETUP_Pert(teilchenPert)
       call VolumeElements_SETUP_Real(teilchenReal)
!       call VolumeElements_Statistics
    end if

    !**************************************************
    ! perturbative Particle with two real ones
    !**************************************************

    Do i=lBound(teilchenPert,dim=1),uBound(teilchenPert,dim=1)
       Do j=lBound(teilchenPert,dim=2),uBound(teilchenPert,dim=2)
          If( teilchenPert(i,j)%ID.ne.pion.and.teilchenPert(i,j)%ID.ne.delta) cycle
          
          IF (Associated(proton1))  Nullify(proton1)
          IF (Associated(proton2))  Nullify(proton2)
          IF (Associated(neutron1)) Nullify(neutron1)
          IF (Associated(neutron2)) Nullify(neutron2)
          
          ! Find closest nucleons
          If(fullensemble) then
             If(nukSearch_VE.and.localEnsemble) then
                call VolumeElements_NukSearch(teilchenPert(i,j), GetRadiusNukSearch(),&
                     & proton1,proton2,neutron1,neutron2, successFlag)
             else
                call NukSearch(teilchenPert(i,j),teilchenReal,&
                     & proton1,proton2,neutron1,neutron2, successFlag)
             end if
          else
             call NukSearch(teilchenPert(i,j),teilchenReal(i:i,:),&
                  & proton1,proton2,neutron1,neutron2, successFlag)
          end if
          
          If (.not.successFlag) cycle
          call make_3Body_Collision(teilchenPert(i,j),proton1,proton2,neutron1,neutron2, & 
                                    scatterPartner1,scatterPartner2,finalstate,successFlag)
          
          If (.not.successFlag) cycle
          
          If (teilchenPert(i,j)%ID==Delta .and. JustDeleteDelta) then 
             ! just delete the delta, no final state
             finalState%event(1)=pert_numbering(scatterPartner1)
             finalState%event(2)=pert_numbering(scatterPartner2)
             call ReportEventNumber((/teilchenPert(i,j),scatterPartner1,scatterPartner2/), &
                                    finalState, finalState(1)%event,time,3112)
             ! (2) Eliminate incoming perturbative particles
             teilchenPert(i,j)%Id=0
             cycle
          end if
          
          successflag = finalCheck ((/teilchenPert(i,j),scatterPartner1,scatterPartner2/), finalState, .false.) 
          If (.not. successFlag) cycle
          If (teilchenPert(i,j)%ID/=Delta) then
             call checkPauli(finalState,teilchenReal,successflag)
          else
             ! Pauli Blocking is already included in the delta "decay width", therefore we should not count it twice!
             successFlag = .true.
          end if
          if (.not. successflag) cycle 
          ! (1) Label event by eventNumber, such that we can track the perturbative particle to its 
          ! real scattering partner.
          finalState%event(1)=pert_numbering(scatterPartner1)
          finalState%event(2)=pert_numbering(scatterPartner2)

          finalState%lastCollisionTime = time
          finalState%perweight = teilchenPert(i,j)%perweight
          finalState%firstEvent = teilchenPert(i,j)%firstEvent
          finalState%perturbative = .true.
          
          do ii=1,size(finalState)
            If (particlePropagated(finalState(ii))) call setNumber(finalState(ii))
          enddo

          call ReportEventNumber((/teilchenPert(i,j),scatterPartner1,scatterPartner2/), finalState, finalState(1)%event, time, 3112)

          ! (2) Eliminate incoming perturbative particles
          pairS2 = teilchenPert(i,j)
          teilchenPert(i,j)%Id=0
          ! (3) Create final state particles
          posOut = 0
          If (particlePropagated(finalState(1))) then
             teilchenPert(i,j) = finalState(1)
             posOut(1:2,1) = (/i,j/)
          end if
          If (fullensemble) then
             call setIntoVector(finalState(2:), teilchenPert, successFlag, .true., positions=posOut(:,2:))
          else
             call setIntoVector(finalState(2:), teilchenPert(i:i,:), successFlag, .true., positions=posOut(:,2:))
             posOut(1,2:) = i
          end if

          call CollHist_UpdateHist((/pairS2,scatterPartner1,scatterPartner2/), finalState, (/i,j/), posOut, finalState(1)%perweight)

          ! (4) Check that setting into perturbative particle vector worked out 
          If (.not.successFlag) then
             write(*,*) 'Perturbative particle vector too small!'
             write(*,*) size(finalState),  lBound(teilchenPert), uBound(teilchenPert)
             write(*,*) 'Dumping real particles to files "RealParticles_stop_*!"'
             call WriteParticleVector("RealParticles_stop",teilchenReal)
             write(*,*) 'Dumping perturbative particles to files "PertParticles_stop_*!"'
             call WriteParticleVector("PertParticles_stop",teilchenPert)
             call Traceback('collision Term, three-body')
          end if
       end do
    end do
    
    !**************************************************************************
    ! real with two real ones
    !**************************************************************************

    Do i=lBound(teilchenReal,dim=1), uBound(teilchenReal,dim=1)
       Do j=lBound(teilchenReal,dim=2), uBound(teilchenReal,dim=2)
          If (teilchenReal(i,j)%ID/=pion .and. teilchenReal(i,j)%ID/=delta) cycle
          
          IF (Associated(proton1))  Nullify(proton1)
          IF (Associated(proton2))  Nullify(proton2)
          IF (Associated(neutron1)) Nullify(neutron1)
          IF (Associated(neutron2)) Nullify(neutron2)
          
          ! Find closest nucleons
          If (fullensemble) then
            call NukSearch (teilchenReal(i,j), teilchenReal, proton1, proton2, neutron1, neutron2, successFlag)
          else
            call NukSearch (teilchenReal(i,j), teilchenReal(i:i,:), proton1, proton2, neutron1, neutron2, successFlag)
          end if

          If (.not. successFlag) cycle
          call make_3Body_Collision (teilchenReal(i,j), proton1, proton2, neutron1, neutron2, & 
                                     scatterPartner1, scatterPartner2, finalstate, successFlag)
          If (.not. successFlag) cycle
          
          If (teilchenReal(i,j)%ID==Delta .and. JustDeleteDelta) then 
             number = real_numbering()
             call ReportEventNumber((/teilchenReal(i,j),scatterPartner1,scatterPartner2/), &
                                    finalState, (/number,number/), time, 3111)
             ! convert the delta to a nucleon:
             teilchenReal(i,j)%Id = nucleon
             ! dont care about total charge:
             teilchenReal(i,j)%charge = (teilchenReal(i,j)%charge + 1)/2
             teilchenReal(i,j)%event = number
             cycle
          end if

          successflag = finalCheck ((/teilchenReal(i,j),scatterPartner1,scatterPartner2/), finalState, .false.) 
          If (.not. successFlag) cycle
          If (teilchenReal(i,j)%ID/=Delta) then
             call checkPauli(finalState,teilchenReal,successflag)
          else
             ! Pauli Blocking is already included in the delta "decay width", therefore we should not count it twice!
             successFlag = .true.
          end if
          if (.not. successflag) cycle 

          ! (1) Label event by eventNumber
          number = real_numbering()
          finalState%event(1) = number
          finalState%event(2) = number
          
          finalState%lastCollisionTime = time
          finalState%perturbative = .false.
          finalState%perweight = min(scatterPartner1%perweight, scatterPartner2%perweight, teilchenReal(i,j)%perweight)
          number = min(scatterPartner1%firstEvent, scatterPartner2%firstEvent, teilchenReal(i,j)%firstEvent)
          if (number==0) number = real_firstnumbering()
          finalState%firstEvent = number
        
          call ReportEventNumber((/teilchenReal(i,j),scatterPartner1,scatterPartner2/), finalState, finalState(1)%event, time, 3111)

          ! (2) Eliminate incoming particles
          teilchenReal(i,j)%Id = 0
          scatterPartner1%Id = 0
          scatterPartner2%Id = 0

          ! (3) Create final state particles
          If (particlePropagated(finalState(1))) then
             call setNumber(finalState(1))
             scatterPartner1 = finalState(1)
          end if
          If (particlePropagated(finalState(2))) then
             call setNumber(finalState(2))
             scatterPartner2 = finalState(2)
          end if
          If (particlePropagated(finalState(3))) then
             call setNumber(finalState(3))
             teilchenReal(i,j) = finalState(3)
          end if
       end Do
    end Do

  end subroutine threeBody
 
  
  !*****************************************************************************
  !****f* collisionTerm/finalCheck
  ! NAME
  ! function finalCheck (teilchenIn, teilchenOut, HiEnergyFlagge, woher) result(flag)
  ! PURPOSE
  ! Checks the final state of a collision for the conservation of all 
  ! quantum numbers,  also including momentum and energy conservation. 
  !
  ! For HiEnergy events we do not check charge and momentum conservation, 
  ! this MUST be done separately. 
  ! The reason for this is, that some particles which are produced by 
  ! Pythia/Fritiof can not be propagated by BUU and therefore do not show up in 
  ! the final state vector "teilchenOut"
  ! ("unknown particles wont be propagated").
  !
  ! INPUTS
  ! * type(particle),dimension(:)  :: teilchenIn  -- Incoming particles
  ! * type(particle),dimension(:)  :: teilchenOut -- Outgoing particles
  ! * logical, optional            :: HiEnergyFlag -- 
  !   .true. if it was a HiEnergy event.
  !   if .true. then energy conservation is not checked
  !   and code does not stop if charge conservation is violated
  ! * character(*), optional      :: woher -- ...
  !
  ! OUTPUT
  ! * logical :: flag -- .true. if quantum numbers are conserved
  !*****************************************************************************
  function finalCheck (teilchenIn, teilchenOut, HiEnergyFlag, woher) result(flag)

    use particleDefinition
    use IdTable, only: isBaryon, isHadron
    use particleProperties, only: hadron, validCharge
    use callStack, only: traceBack
    use output, only: WriteParticle

    type(particle), intent(in), dimension(:) :: teilchenIn, teilchenOut
    logical, intent(in), optional :: HiEnergyFlag
    character(*), intent(in), optional :: woher
    logical :: flag

    integer :: totalCharge_in, totalCharge_out, baryon_number_In, baryon_number_Out, strangeness_In, strangeness_Out, k
    real, dimension(0:3) :: momentum_In, momentum_Out
    real :: m_miss2  ! ,s_In, s_Out

    flag=.false.

    ! Initialize at first call
    If (initFlag) call ReadInput

    ! Check charge of each particle
    Do k=lbound(teilchenOut,dim=1),ubound(teilchenOut,dim=1)
       If(.not.validCharge(teilchenOut(k))) then
          iF(present(woher)) write(*,*) 'Problem in ',trim(woher),' !!!!!!'
          Write(*,*) 'Charge of particle is not valid in finalCheck'
          call WriteParticle(6,1,k,teilchenOut(k))
          Write(*,*) 'Charge=',teilchenOut(k)%charge
          Write(*,*) 'Id=',teilchenOut(k)%id
          Write(*,*) 'Antiparticle=',teilchenOut(k)%antiparticle
          Write(*,*) 'Momentum=',teilchenOut(k)%momentum
          Write(*,*)
          If (Present(HiEnergyFlag)) then
             If (HiEnergyFlag) then
                write(*,*) 'It was a HiEnergy event!'
                call PYLIST(2)
             else
                write(*,*) 'It was a low-energy event!'
             end if
          else
             write(*,*) 'It was a decay-event!'
          end if
          Write(*,*)
          call printTeilchen (teilchenIn, .true.)
          call printTeilchen (teilchenOut, .false.)
          call Traceback('Severe problem. Code stops!!!!')
       end if
    End do

    if (present(HiEnergyFlag)) then
      if (HiEnergyFlag) then
        flag = .true.
        return
      end if
    end if

    ! Check baryon number conservation:
    baryon_number_In = 0
    do k=lbound(teilchenIn,dim=1),ubound(teilchenIn,dim=1)
      if (isBaryon(teilchenIn(k)%Id)) then
        if (.not.teilchenIn(k)%antiparticle) then
          baryon_number_In = baryon_number_In + 1
        else
          baryon_number_In = baryon_number_In - 1
        end if
      end if
    end do
    baryon_number_Out = 0
    do k=lbound(teilchenOut,dim=1),ubound(teilchenOut,dim=1)
      if (isBaryon(teilchenOut(k)%Id)) then
        if (.not.teilchenOut(k)%antiparticle) then
          baryon_number_Out = baryon_number_Out + 1
        else
          baryon_number_Out = baryon_number_Out - 1
        end if
      end if
    end do
    if (baryon_number_In /= baryon_number_Out) then
      if (Present(woher)) then 
        write(*,*)'Problems in '//trim(woher)//' : Baryon number conservation'
      else
        write(*,*)'Problems in collisionTerm: Baryon number conservation'
      end if
      write(*,*)'Baryon number in: ', baryon_number_In
      write(*,*)'Baryon number out:', baryon_number_Out
      call printTeilchen (teilchenIn, .true.)
      call printTeilchen (teilchenOut, .false.)
      call Traceback('stop in finalCheck, baryon number check')
    end if
    
    ! Check strangeness conservation:
    strangeness_In = 0
    do k=lbound(teilchenIn,dim=1),ubound(teilchenIn,dim=1)
      if (.not. isHadron(teilchenIn(k)%Id)) cycle
      if (.not. teilchenIn(k)%antiparticle) then
        strangeness_In = strangeness_In + hadron(teilchenIn(k)%Id)%strangeness
      else
        strangeness_In = strangeness_In - hadron(teilchenIn(k)%Id)%strangeness
      end if
    end do
    strangeness_Out = 0
    do k=lbound(teilchenOut,dim=1),ubound(teilchenOut,dim=1)
      if (.not. isHadron(teilchenOut(k)%Id)) cycle
      if (.not. teilchenOut(k)%antiparticle) then
        strangeness_Out = strangeness_Out + hadron(teilchenOut(k)%Id)%strangeness
      else
        strangeness_Out = strangeness_Out - hadron(teilchenOut(k)%Id)%strangeness
      end if
    end do
    if (strangeness_In /= strangeness_Out) then
      if (Present(woher)) then 
        write(*,*)'Problems in '//trim(woher)//' : Strangeness conservation'
      else
        write(*,*)'Problems in collisionTerm: Strangeness conservation'
      end if
      write(*,*)'Strangeness in: ', strangeness_In
      write(*,*)'Strangeness out:', strangeness_Out
      call printTeilchen (teilchenIn, .true.)
      call printTeilchen (teilchenOut, .false.)
      !  Check conservation of srts
      Do k=0,3
        momentum_IN(k) =Sum(teilchenIn%momentum(k))
        momentum_Out(k)=Sum(teilchenOut%momentum(k))
      End do
      m_miss2 = (momentum_IN(0)-momentum_Out(0))**2 - (momentum_IN(1)-momentum_Out(1))**2 &
               -(momentum_IN(2)-momentum_Out(2))**2 - (momentum_IN(3)-momentum_Out(3))**2
      write(*,*)' (Missing mass)^2 :', m_miss2
      call Traceback('stop in finalCheck, strangeness check')
    end if 


    !  Check total charge
    totalCharge_In  =   Sum(teilchenIn%charge)
    totalCharge_Out =   Sum(teilchenOut%charge)
    if (totalCharge_In /= totalCharge_Out) then
      if(Present(woher)) then 
        write(*,*)'Problems in '//trim(woher)//' : Charge conservation'
      else
        write(*,*)'Problems in collisionTerm: Charge conservation'
      end if
      write(*,*)'Charge in: ', totalCharge_In
      write(*,*)'Charge out:', totalCharge_Out
      call printTeilchen (teilchenIn, .true.)
      call printTeilchen (teilchenOut, .false.)
      call Traceback('stop in finalCheck, charge check')
    end if

!       if(debug) Print *, 'CollisionTerm: SqrtS=',sqrtS(teilchenIN(1),teilchenIN(2))

    !  Check conservation of srts
    Do k=0,3
      momentum_IN(k) =Sum(teilchenIn%momentum(k))
      momentum_Out(k)=Sum(teilchenOut%momentum(k))
    End do

    !s_In=momentum_IN(0)**2-Dot_Product(momentum_IN(1:3), momentum_IN(1:3))
    !s_Out=momentum_Out(0)**2-Dot_Product(momentum_Out(1:3), momentum_Out(1:3))

    Do k=0,3
      if(abs(momentum_In(k)-momentum_Out(k)).gt.energyCheck) then
        if (Present(woher)) then 
          write(*,*)'Problems in '//trim(woher)//' : Energy/momentum conservation'
        else
          write(*,*)'Problems in collisionTerm: Energy/momentum conservation'
        end if
        write(*,'(A,4G13.5)')'Momentum in: ', Momentum_In
        write(*,'(A,4G13.5)')'Momentum out:', Momentum_Out
        write(*,*)

        call printTeilchen (teilchenIn, .true.)
        call printTeilchen (teilchenOut, .false.)
        call Traceback('stop in finalCheck, momentum check')
      end if
    end do

    flag=.true.
    
  contains
  
    subroutine printTeilchen (teilchen, in)
      use mediumModule, only: mediumAt
      use output, only: writeParticle_debug
      type(particle), intent(in), dimension(:) :: teilchen
      logical, intent(in) :: in
      integer :: i
      Do i=lBound(teilchen,dim=1) , uBound(teilchen,dim=1)
        if (teilchen(i)%ID==0) cycle
        if (in) then
          write(*,*) 'Incoming Particle number #', i
        else
          write(*,*) 'Outgoing Particle number #', i
        end if
        call writeParticle_debug(teilchen(i),mediumAt(teilchen(i)%position))
        write(*,*)
      end do
      write(*,*)
      write(*,*) '**************************************************'
    end subroutine printTeilchen

  end function finalCheck


end module collisionTerm
