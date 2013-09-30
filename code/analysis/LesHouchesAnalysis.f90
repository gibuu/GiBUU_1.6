!***************************************************************************
!****m* /LesHouchesAnalysis
! NAME
! module LesHouchesAnalysis
!
! PURPOSE
! This module provides routines for writing the outgoing particle
! vector in the standard format according the "Les Houches Event Files".
!
! cf. http://arxiv.org/pdf/hep-ph/0609017v1
!
! INPUTS
! (none)
!
! NOTES
! These analysis routines are independent of the specific initialization
! and should work in principle for all event types.
!
!***************************************************************************

!***************************************************************************
!****n* LesHouchesAnalysis/LesHouches
! PURPOSE
! Namelist for LesHouchesAnalysis includes:
! * LesHouchesFinalParticles_Pert
! * LesHouchesFinalParticles_Real
!***************************************************************************

module LesHouchesAnalysis

  implicit none

  PRIVATE

  !*************************************************************************
  !****g* LesHouchesAnalysis/LesHouchesFinalParticles_Pert
  ! SOURCE
  !
  logical,save :: LesHouchesFinalParticles_Pert = .false.
  !
  ! PURPOSE
  ! flag: Print perturbative particle vector according LesHouches standard
  !*************************************************************************

  !*************************************************************************
  !****g* LesHouchesAnalysis/LesHouchesFinalParticles_Real
  ! SOURCE
  !
  logical,save :: LesHouchesFinalParticles_Real = .false.
  !
  ! PURPOSE
  ! flag: Print real particle vector according LesHouches standard
  !*************************************************************************

  logical, save :: init = .true.

  public :: DoLesHouchesAnalysis

contains

  !*************************************************************************
  !****s* LesHouchesAnalysis/initInput
  ! NAME
  ! subroutine initInput
  !
  ! PURPOSE
  ! read namelist
  !*************************************************************************
  subroutine initInput()
    use output
    
    NAMELIST /LesHouches/  LesHouchesFinalParticles_Pert, LesHouchesFinalParticles_Real

    integer :: ios

    call Write_ReadingInput('LesHouches',0)
    rewind(5)
    read(5,nml=LesHouches,iostat=ios)

    call Write_ReadingInput('LesHouches',0,ios)
    
    write(*,*) 'LesHouches output of final particles (real,pert): ',&
         &LesHouchesFinalParticles_Real,LesHouchesFinalParticles_Pert
    call Write_ReadingInput('LesHouches',1)

    init = .false.

  end subroutine initInput

  !*************************************************************************
  !****s* LesHouchesAnalysis/DoLesHouchesAnalysis
  ! NAME
  ! subroutine DoLesHouchesAnalysis
  !
  ! PURPOSE
  ! Do the actual writing out, if desired as in namelist
  !*************************************************************************
  subroutine DoLesHouchesAnalysis (realPart, pertPart)
    use particleDefinition

    type(particle), intent(in), dimension(:,:) :: realPart, pertPart

    integer, save :: nCall = 0
    character*10 :: BUF

    if (init) call initInput

    nCall = nCall+1

    !*****************************************************************************
    !****o* LesHouchesAnalysis/LesHouches.Pert.xml
    ! NAME
    ! file LesHouches.Pert.xml
    ! PURPOSE
    ! Contains all perturbative particles of a given run in Les Hoches format.
    ! Can be enabled by the switch LesHouchesFinalParticles_Pert.
    ! For each subsequent run a separate file will be produced:
    !  * LesHouches.Pert.00000001.xml
    !  * LesHouches.Pert.00000002.xml
    !  * etc
    !*****************************************************************************
    if (LesHouchesFinalParticles_Pert) then
       write(BUF,'(i8.8)') nCall
       call LesHouchesFileOpen(721,'LesHouches.Pert.'//trim(BUF)//'.xml')
       call LesHouchesWriteEvents(721,pertPart)
       call LesHouchesFileClose(721)
    end if

    !*****************************************************************************
    !****o* LesHouchesAnalysis/LesHouches.Real.xml
    ! NAME
    ! file LesHouches.Real.xml
    ! PURPOSE
    ! Contains all real particles of a given run in Les Hoches format.
    ! Can be enabled by the switch LesHouchesFinalParticles_Real.
    ! For each subsequent run a separate file will be produced:
    !  * LesHouches.Real.00000001.xml
    !  * LesHouches.Real.00000002.xml
    !  * etc
    !*****************************************************************************
    if (LesHouchesFinalParticles_Real) then
       write(BUF,'(i8.8)') nCall
       call LesHouchesFileOpen(722,'LesHouches.Real.'//trim(BUF)//'.xml')
       call LesHouchesWriteEvents(722,realPart)
       call LesHouchesFileClose(722)
    end if

  end subroutine DoLesHouchesAnalysis


  !*************************************************************************
  !****s* LesHouchesAnalysis/LesHouchesFileOpen
  ! NAME
  ! subroutine LesHouchesFileOpen (iFile, fName)
  ! 
  ! PURPOSE
  ! open a file for output event information according the 
  ! "Les Houches Event Files" standard.
  !*************************************************************************
  subroutine LesHouchesFileOpen (iFile, fName)
    integer, intent(in)       :: iFile
    character*(*), intent(in) :: fName

    open(iFile, file=fName, status='unknown')
    rewind(iFile)

    write(iFile,'(A)') '<LesHouchesEvents version="1.0">'
    write(iFile,'(A)') '<!--'
    write(iFile,'(A)') '   File generated by GiBUU.'
    write(iFile,*)
    write(iFile,'(A)') '   NOTES:'
    write(iFile,'(A)') '   * compulsory init information not applicable;'
    write(iFile,'(A)') '     Line may produce errors in some XML parsers.'
    write(iFile,'(A)') '   * There may be no closing tag in the last line.'
    write(iFile,'(A)') '   * In the case of HiLepton events, at the end of every event'
    write(iFile,'(A)') '     block, an additional line is given as'
    write(iFile,'(A)') '     " # 14 nu Q2 eps phiLepton Eventtype"'
    write(iFile,'(A)') '     (14 is the magic number of HiLepton)'
    write(iFile,'(A)') '     (cf. Electron_origin.f90 for a list of Eventtypes)'
    write(iFile,'(A)') '   * In the case of neutrino events, at the end of every event'
    write(iFile,'(A)') '     block, an additional line is given as'
    write(iFile,'(A)') '     " # 5 Eventtype Weight momLepIn(0:3) momLepOut(0:3)"'
    write(iFile,'(A)') '     (5 is the magic number of neutrino induced events)'
    write(iFile,'(A)') '     (Eventtype: 1=QE, 2-31=res ID, 32,33=1pi, 34=DIS, 35,36=2p2h, 37=2pi)'
    write(iFile,'(A)') '-->' 

    write(iFile,'(A)') '<header>'
    write(iFile,'(A)') '     <!-- individual XML tags may follow -->'
    write(iFile,'(A)') '</header>'

    write(iFile,'(A)') '<init>'
    write(iFile,'(1P,2I8,2E14.6,6I6)') 0,0, 0.0,0.0, 0,0,0,0,0,0
    write(iFile,'(A)') '</init>'

  end subroutine LesHouchesFileOpen


  !*************************************************************************
  !****s* LesHouchesAnalysis/LesHouchesFileClose
  ! NAME
  !   subroutine LesHouchesFileClose (iFile)
  ! 
  ! PURPOSE
  ! open a file for output event information according the 
  ! "Les Houches Event Files" standard.
  !*************************************************************************
  subroutine LesHouchesFileClose (iFile)
    integer, intent(in)                 :: iFile

    write(iFile,'(A)') '</LesHouchesEvents>'
    close(iFile)
  end subroutine LesHouchesFileClose


  !*************************************************************************
  !****is* LesHouchesAnalysis/LesHouchesWriteEvents
  ! NAME
  ! subroutine LesHouchesWriteEvents (iFile, Parts)
  !
  ! PURPOSE
  ! Do the actual printout
  !
  ! NOTES
  ! We have to sort the particles according their "firstevent" field.
  ! Therefore we allocate an array of "tParticleList". Unfortunately we can
  ! not use the "firstevent" entry directly as array index, since this
  ! is not starting with 1 and continously increasing for all kind
  ! of eventtypes. Therefore we (ab)use the module "PILIndex", which 
  ! implements methods of "indexing". (We do not use the possibility of
  ! reallocating as provided by the module "PILIndex".)
  !*************************************************************************
  subroutine LesHouchesWriteEvents (iFile, Parts)
    use particleDefinition
    use particlePointerListDefinition
    use particlePointerList, only: ParticleList_INIT, ParticleList_APPEND, ParticleList_CLEAR
    use PILIndex, only: tIndexList, PILIndex_PUT
    use ID_translation, only: KFfromBUU

    integer,        intent(in)                         :: iFile
    type(particle), intent(in), dimension(:,:), target :: Parts

    type(tIndexList), save :: IndexList
    type(tParticleList), allocatable, save :: ValueList(:)

    integer :: i,iEns,iPart, iFE,iiFE, KF

    type(particle), pointer :: pPart
    type(tParticleListNode),Pointer  :: pNode

    integer :: NUP,IDPRUP
    real :: XWGTUP,SCALUP,AQEDUP,AQCDUP

    character(len=15), parameter :: f1 = '(1P,2I6,4E14.6)'
    character(len=22), parameter :: f2 = '(1P,I8,5I5,5E18.10,A6)'

    IDPRUP = 0
    XWGTUP = 0.0 ! weight of event
    SCALUP = 0.0
    AQEDUP = 0.0
    AQCDUP = 0.0

    !
    ! Clean up the arrays:
    !
    
    if (allocated(ValueList)) then
       do i=1,size(ValueList)
          call ParticleList_CLEAR(ValueList(i))
       end do
    endif
    IndexList%nEntry = 0

    !
    ! Loop over all particles and group them according their first event:
    !

    do iEns = 1,size(Parts,dim=1)
       PartLoop:do iPart = 1,size(Parts,dim=2)
          pPart => Parts(iEns,iPart)
          if (pPart%ID <= 0) cycle PartLoop

          iFE = pPart%firstEvent
          if (iFE.eq.0) cycle PartLoop ! particle did not interact !
          iiFE = PILIndex_PUT(IndexList, iFE, 'LesHouches')
          if (iiFE>0) then
             call ParticleList_APPEND(ValueList(iiFE), pPart)
          else
             call ValuelistAllocate()
             call ParticleList_APPEND(ValueList(-iiFE), pPart)
          endif
 
       end do PartLoop
    end do

    !
    ! Loop over all events and write them to file:
    !
    if (.not.allocated(ValueList)) return

    do iiFE=1,size(ValueList)
       NUP = ValueList(iiFE)%nEntries ! number of particles
       if (NUP .eq. 0) cycle

       pNode => ValueList(iiFE)%first
       XWGTUP = pNode%V%perweight ! weight of event
       iFE = pNode%V%firstEvent

       write(iFile,'(A)') '<event>'
       write(iFile,f1) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP

       do 
          if (.not. ASSOCIATED(pNode)) exit

          pPart => pNode%V

          KF = KFfromBUU (pPart)

          write(iFile,f2) KF, 0, 0,0, 0,0, &
                          pPart%momentum(1:3), pPart%momentum(0), &
                          sqrts(pPart), '0. 9.'

          pNode => pNode%next
       end do

       call WriteAdditionalInfo(iFile,iFE)

       write(iFile,'(A)') '</event>'
    end do

    contains

      !*********************************************************************
      !****is* LesHouchesWriteEvents/ValueListAllocate
      ! NAME
      ! subroutine ValueListAllocate
      !
      ! PURPOSE
      ! Do the allocation stuff for the Particle Info List
      !*********************************************************************
      subroutine ValueListAllocate

        integer :: n0, n1,i
        type(tParticleList),allocatable :: L0(:)

        n1 = size(IndexList%PartNumber) ! new size

        if (.not.allocated(ValueList)) then
           allocate(ValueList(n1))
           do i=1,n1
              call ParticleList_INIT(ValueList(i))
           enddo
           return
        endif

        n0 = size(ValueList)            ! old size

        allocate(L0(n0))
        do i=1,n0
           L0(i)%first => ValueList(i)%first
           L0(i)%last  => ValueList(i)%last
           L0(i)%nEntries = ValueList(i)%nEntries
        end do
        deallocate(ValueList)
        allocate(ValueList(n1))
        do i=1,n0
           ValueList(i)%first => L0(i)%first
           ValueList(i)%last  => L0(i)%last
           ValueList(i)%nEntries = L0(i)%nEntries
        end do
        do i=n0+1,n1
           call ParticleList_INIT(ValueList(i))
        end do
        deallocate(L0)

      end subroutine ValueListAllocate

      !*********************************************************************
      !****is* LesHouchesWriteEvents/WriteAdditionalInfo
      ! NAME
      ! subroutine WriteAdditionalInfo(iFile,iFE)
      !
      ! PURPOSE
      ! Write additional info according the event.
      !
      ! This routine tries to find additional information about the event.
      ! It tries routines for different event types, which only return
      ! some information, if it was really stored.
      !
      ! it uses:
      ! * module "EventInfo_HiLep", EventInfo_HiLep_Get:
      !   If it finds something, then a line
      !     # 14 nu Q2 eps phiLepton Eventtype
      !   is added (14 is the magic number of "High energetic lepton induced) 
      ! * module "neutrinoProdInfo", NeutrinoProdInfo_Get:
      !   If it finds something, then a line
      !     # 5 Eventtype Weight  momLepIn(0:3) momLepOut(0:3)
      !   is added (5 is the magic naumber for neutrino events)
      ! 
      !*********************************************************************

      subroutine WriteAdditionalInfo(iFile,iFE)
        use EventInfo_HiLep, only: EventInfo_HiLep_Get
        use neutrinoProdInfo, only: NeutrinoProdInfo_Get

        integer, intent(in) :: iFile,iFE
        real :: weight,nu,Q2,eps,phiL
        integer :: EventType
        real,dimension(0:3) :: momLepIn, momLepOut, momBos

        if (EventInfo_HiLep_Get(0,iFE,Weight,nu,Q2,eps,EventType,phi_Lepton=phiL)) then
           write(iFile,'(A,1P,4e13.4,0P,I8)') '# 14 ',nu,Q2,eps,phiL,EventType
        else if(NeutrinoProdInfo_Get(iFE,EventType,Weight,momLepIn,momLepOut,momBos)) then
           write(iFile,'(A,I5,1P,e18.10,1P,2(" ",4e18.10),0P)') &
                '# 5 ', EventType, Weight, momLepIn, momLepOut
        end if

      end subroutine WriteAdditionalInfo

  end subroutine LesHouchesWriteEvents


end module LesHouchesAnalysis
