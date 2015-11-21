!*******************************************************************************
!****m* /LesHouchesAnalysis
! NAME
! module LesHouchesAnalysis
!
! PURPOSE
! This module provides routines for writing the outgoing particle
! vector in the standard format according to the "Les Houches Event Files".
!
! See also:
! * http://arxiv.org/abs/hep-ph/0609017
! * https://gibuu.hepforge.org/trac/wiki/LesHouches
!
! INPUTS
! (none)
!
! NOTES
! These analysis routines are independent of the specific initialization
! and should work in principle for all event types.
!*******************************************************************************
module LesHouchesAnalysis

  implicit none
  private

  !*****************************************************************************
  !****g* LesHouchesAnalysis/LesHouchesFinalParticles_Pert
  ! SOURCE
  !
  logical, save :: LesHouchesFinalParticles_Pert = .false.
  !
  ! PURPOSE
  ! Flag to print the perturbative particle vector according to the "Les Houches"
  ! standard, cf. http://arxiv.org/abs/hep-ph/0609017 .
  ! The output will be written to LesHouches.Pert.xml
  !*****************************************************************************

  !*****************************************************************************
  !****g* LesHouchesAnalysis/LesHouchesFinalParticles_Real
  ! SOURCE
  !
  logical, save :: LesHouchesFinalParticles_Real = .false.
  !
  ! PURPOSE
  ! Flag to print the real particle vector according to the "Les Houches"
  ! standard, cf. http://arxiv.org/abs/hep-ph/0609017 .
  ! The output will be written to LesHouches.Real.xml.
  !*****************************************************************************

  logical, save :: init = .true.

  public :: DoLesHouchesAnalysis

contains


  !*****************************************************************************
  !****s* LesHouchesAnalysis/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Read namelist 'LesHouches' from jobcard.
  !*****************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput

    !***************************************************************************
    !****n* LesHouchesAnalysis/LesHouches
    ! NAME
    ! namelist /LesHouches/
    ! PURPOSE
    ! Namelist for LesHouchesAnalysis includes:
    ! * LesHouchesFinalParticles_Pert
    ! * LesHouchesFinalParticles_Real
    !***************************************************************************
    namelist /LesHouches/ LesHouchesFinalParticles_Pert, LesHouchesFinalParticles_Real

    integer :: ios

    call Write_ReadingInput('LesHouches',0)
    rewind(5)
    read(5,nml=LesHouches,iostat=ios)

    call Write_ReadingInput('LesHouches',0,ios)
    
    write(*,*) 'LesHouches output of final particles (real,pert): ', &
               LesHouchesFinalParticles_Real, LesHouchesFinalParticles_Pert
    call Write_ReadingInput('LesHouches',1)

    init = .false.

  end subroutine initInput


  !*****************************************************************************
  !****s* LesHouchesAnalysis/DoLesHouchesAnalysis
  ! NAME
  ! subroutine DoLesHouchesAnalysis (realPart, pertPart)
  ! PURPOSE
  ! Do the actual writing out, if desired (as indicated in namelist).
  !*****************************************************************************
  subroutine DoLesHouchesAnalysis (realPart, pertPart)
    use particleDefinition

    type(particle), intent(in), dimension(:,:) :: realPart, pertPart

    integer, save :: nCall = 0
    character*10 :: BUF

    if (init) call initInput

    nCall = nCall+1

    !***************************************************************************
    !****o* LesHouchesAnalysis/LesHouches.Pert.xml
    ! NAME
    ! file LesHouches.Pert.xml
    ! PURPOSE
    ! Contains all perturbative particles of a given run in Les Hoches format.
    ! Can be enabled by the switch LesHouchesFinalParticles_Pert.
    ! For documentation of the file format see https://gibuu.hepforge.org/trac/wiki/LesHouches.
    ! For each subsequent run a separate file will be produced:
    !  * LesHouches.Pert.00000001.xml
    !  * LesHouches.Pert.00000002.xml
    !  * etc
    !***************************************************************************
    if (LesHouchesFinalParticles_Pert) then
       write(BUF,'(i8.8)') nCall
       call LesHouchesFileOpen(721,'LesHouches.Pert.'//trim(BUF)//'.xml')
       call LesHouchesWriteEvents_pert(721,pertPart)
       call LesHouchesFileClose(721)
    end if

    !***************************************************************************
    !****o* LesHouchesAnalysis/LesHouches.Real.xml
    ! NAME
    ! file LesHouches.Real.xml
    ! PURPOSE
    ! Contains all real particles of a given run in Les Hoches format.
    ! Can be enabled by the switch LesHouchesFinalParticles_Real.
    ! For documentation of the file format see https://gibuu.hepforge.org/trac/wiki/LesHouches.
    ! For each subsequent run a separate file will be produced:
    !  * LesHouches.Real.00000001.xml
    !  * LesHouches.Real.00000002.xml
    !  * etc
    !***************************************************************************
    if (LesHouchesFinalParticles_Real) then
       write(BUF,'(i8.8)') nCall
       call LesHouchesFileOpen(722,'LesHouches.Real.'//trim(BUF)//'.xml')
       call LesHouchesWriteEvents_real(722,realPart)
       call LesHouchesFileClose(722)
    end if

  end subroutine DoLesHouchesAnalysis


  !*****************************************************************************
  !****s* LesHouchesAnalysis/LesHouchesFileOpen
  ! NAME
  ! subroutine LesHouchesFileOpen (iFile, fName)
  ! PURPOSE
  ! Open a file for output event information according to the 
  ! "Les Houches Event Files" standard.
  !*****************************************************************************
  subroutine LesHouchesFileOpen (iFile, fName)
    integer, intent(in)       :: iFile
    character*(*), intent(in) :: fName

    open(iFile, file=fName, status='unknown')
    rewind(iFile)

    write(iFile,'(A)') '<LesHouchesEvents version="1.0">'
    write(iFile,'(A)') '<!-- File generated by GiBUU. For documentation see https://gibuu.hepforge.org/trac/wiki/LesHouches -->'

    write(iFile,'(A)') '<header>'
    write(iFile,'(A)') '     <!-- individual XML tags may follow -->'
    write(iFile,'(A)') '</header>'

    write(iFile,'(A)') '<init>'
    write(iFile,'(1P,2I8,2E14.6,6I6)') 0,0, 0.,0., 0,0,0,0,0,0
    write(iFile,'(A)') '</init>'

  end subroutine LesHouchesFileOpen


  !*****************************************************************************
  !****s* LesHouchesAnalysis/LesHouchesFileClose
  ! NAME
  ! subroutine LesHouchesFileClose (iFile)
  ! PURPOSE
  ! Close a file after outputting Les-Houches event information.
  !*****************************************************************************
  subroutine LesHouchesFileClose (iFile)
    integer, intent(in)                 :: iFile

    write(iFile,'(A)') '</LesHouchesEvents>'
    close(iFile)
  end subroutine LesHouchesFileClose


  !*****************************************************************************
  !****is* LesHouchesAnalysis/LesHouchesWriteEvents_real
  ! NAME
  ! subroutine LesHouchesWriteEvents_real (iFile, Parts)
  ! PURPOSE
  ! Do the actual printout for real particles.
  ! NOTES
  ! For the case of real particles, one event simply corresponds to one ensemble.
  !*****************************************************************************
  subroutine LesHouchesWriteEvents_real (iFile, Parts)
    use particleDefinition
    use ID_translation, only: KFfromBUU
    use IdTable, only: EOV, NOP

    integer,        intent(in)                         :: iFile
    type(particle), intent(in), dimension(:,:), target :: Parts
    
    integer :: iEns,iPart, KF
    integer, dimension(:), allocatable :: nParts

    integer :: NUP,IDPRUP
    real :: XWGTUP,SCALUP,AQEDUP,AQCDUP

    character(len=15), parameter :: f1 = '(1P,2I6,4E14.6)'
    character(len=22), parameter :: f2 = '(1P,I8,5I5,5E18.10,A6)'

    IDPRUP = 0
    XWGTUP = 1.0 ! weight of event
    SCALUP = 0.0
    AQEDUP = 0.0
    AQCDUP = 0.0

    allocate(nParts(1:size(Parts,dim=1)))
    nParts=0

    ! count particles per ensemble
    do iEns = 1,size(Parts,dim=1)
      do iPart = 1,size(Parts,dim=2)
        if (Parts(iEns,iPart)%ID==EOV) exit
        if (Parts(iEns,iPart)%ID==NOP) cycle
        nParts(iEns) = nParts(iEns) + 1
      end do
    end do

    ! Loop over all events and write them to file:
    do iEns = 1,size(Parts,dim=1)
       NUP = nParts(iEns) ! number of particles
       if (NUP == 0) cycle

       write(iFile,'(A)') '<event>'
       write(iFile,f1) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP

       do iPart = 1,size(Parts,dim=2)
         if (Parts(iEns,iPart)%ID==EOV) exit
         if (Parts(iEns,iPart)%ID==NOP) cycle

          KF = KFfromBUU (Parts(iEns,iPart))

          write(iFile,f2) KF, 0, 0,0, 0,0, &
                          Parts(iEns,iPart)%momentum(1:3), Parts(iEns,iPart)%momentum(0), &
                          sqrts(Parts(iEns,iPart)), '0. 9.'
       end do

       call WriteAdditionalInfo (iFile)

       write(iFile,'(A)') '</event>'
    end do

  end subroutine LesHouchesWriteEvents_real


  !*****************************************************************************
  !****is* LesHouchesAnalysis/LesHouchesWriteEvents_pert
  ! NAME
  ! subroutine LesHouchesWriteEvents_pert (iFile, Parts)
  ! PURPOSE
  ! Do the actual printout for perturbative particles.
  ! NOTES
  ! We have to sort the particles according their "firstevent" field.
  ! Therefore we allocate an array of "tParticleList". Unfortunately we can
  ! not use the "firstevent" entry directly as array index, since this
  ! is not starting with 1 and continously increasing for all kind
  ! of eventtypes. Therefore we (ab)use the module "PILIndex", which 
  ! implements methods of "indexing". (We do not use the possibility of
  ! reallocating as provided by the module "PILIndex".)
  !*****************************************************************************
  subroutine LesHouchesWriteEvents_pert (iFile, Parts)
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

       call WriteAdditionalInfo (iFile, iFE)

       write(iFile,'(A)') '</event>'
    end do

    contains

      !*************************************************************************
      !****is* LesHouchesWriteEvents_pert/ValueListAllocate
      ! NAME
      ! subroutine ValueListAllocate
      ! PURPOSE
      ! Do the allocation stuff for the Particle Info List.
      !*************************************************************************
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

  end subroutine LesHouchesWriteEvents_pert


  !*****************************************************************************
  !****s* LesHouchesAnalysis/WriteAdditionalInfo
  ! NAME
  ! subroutine WriteAdditionalInfo (iFile, iFE)
  ! PURPOSE
  ! Write additional info about the event, depending on eventtype.
  !
  ! This routine tries to find additional information about the event.
  ! It tries routines for different event types, which only return
  ! some information, if it was really stored.
  !
  ! The following cases are handled:
  ! * For eventtype "HiLep", the following line is added:
  !     # 14 nu Q2 eps phiLepton Eventtype
  !   (14 is the magic number of "HiLepton") 
  ! * For eventtype "neutrino", the following line is added:
  !     # 5 Eventtype Weight  momLepIn(0:3) momLepOut(0:3)
  !   (5 is the magic number for neutrino events)
  ! * For eventtype "heavyIon", the following line is added:
  !     # 1 b
  !   (1 is the magic number of "heavyIon", b is the impact parameter in fm)
  ! * For eventtype "hadron", the following line is added:
  !     # 300 b
  !   (300 is the magic number of "hadron", b is the impact parameter in fm)
  !*****************************************************************************
  subroutine WriteAdditionalInfo (iFile, iFE)
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use neutrinoProdInfo, only: NeutrinoProdInfo_Get
    use inputGeneral, only: eventType
    use eventtypes, only: hiLepton, neutrino, heavyIon, hadron
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b

    integer, intent(in) :: iFile
    integer, intent(in), optional :: iFE

    real :: weight,nu,Q2,eps,phiL
    integer :: evtType
    real,dimension(0:3) :: momLepIn, momLepOut, momBos, nuc_mom
    integer :: Chrg_Nuc

    select case (eventType)
    case (heavyIon)
      write(iFile,'(A,ES13.4)') '# 1 ', b_HI
    case (hadron)
      write(iFile,'(A,ES13.4)') '# 300 ', b_had
    case (neutrino)
      if (.not. present(iFE)) return
      if (NeutrinoProdInfo_Get (iFE,evtType,Weight,momLepIn,momLepOut,momBos,nuc_mom, Chrg_Nuc)) &
        write(iFile,'(A,I5,1P,e18.10,1P,3(" ",4e18.10),0P,A)') '# 5 ', evtType, Weight, momLepIn, momLepOut, nuc_mom
    case (hiLepton)
      if (.not. present(iFE)) return
      if (EventInfo_HiLep_Get (0,iFE,Weight,nu,Q2,eps,evtType,phi_Lepton=phiL)) &
        write(iFile,'(A,1P,4e13.4,0P,I8)') '# 14 ', nu, Q2, eps, phiL, evtType
    end select

  end subroutine WriteAdditionalInfo


end module LesHouchesAnalysis
