!***************************************************************************
!****m* /VolumeElements
! NAME
! module VolumeElements
!
! PURPOSE
! This module defines all stuff necessary to seperate the whole
! interaction volume in different "volume elements" (also called "VE-cells")
!
! This is necessary for the implementation of "local ensemble" runs.
!
! INPUTS
! (none)
!***************************************************************************
module VolumeElements

  use particleDefinition
  use particlePointerListDefinition
  use particlePointerList

  PRIVATE


  !*************************************************************************
  !****t* VolumeElements/tVolumeElements
  ! PURPOSE
  ! Discretize the whole possible space volume and hold in a 3D array
  ! Lists of particles, which are at the moment in a given coordinate cell
  !
  ! NOTES
  ! Die Arrays "zCoordFilled_xxx" sollen eine Abkürzung für die Loops
  ! über die z-Koordinate darstellen: Wenn für eine z-Koordinate
  ! für keine x- oder y-Zelle Einträge vorhanden sind, dann kann man
  ! diese z-Koord auch ganz schnell überspringen!
  !
  ! Denkbar wäre auch eine Verbesserung durch die Einführung von
  !   integer, dimension(3,2) :: iRange_Used
  ! wobei alle Schleifen statt über
  !   iRange(i,1)..iRange(i,2)
  ! über
  !   iRange_Used(i,1)..iRange_Used(i,2)
  ! laufen würden.
  ! Hierbei können aber nur zusammenhängende Bereiche benutzt werden,
  ! wodurch das benutzte Modell wiederum starken Auftrieb bekommt!!!
  !
  ! SOURCE
  !
  type tVolumeElements
     type(tParticleList), DIMENSION(:,:,:), ALLOCATABLE :: VE_real ! ~ 100 MB
     type(tParticleList), DIMENSION(:,:,:), ALLOCATABLE :: VE_pert ! ~ 100 MB
     real, dimension(3)      :: Delta
     real, dimension(3,2)    :: Range
     integer, dimension(3,2) :: iRange
     logical, dimension(:), ALLOCATABLE :: zCoordFilled_real
     logical, dimension(:), ALLOCATABLE :: zCoordFilled_pert
  end type tVolumeElements
  !*************************************************************************

  !*************************************************************************
  !****ig* VolumeElements/tVE
  ! PURPOSE
  ! The one and only instance of a tVolumeElements-object
  !
  ! SOURCE
  !
  type(tVolumeElements),save :: tVE
  !*************************************************************************


  integer, dimension(3) :: iPart
  type(tParticleListNode), POINTER :: pPert, pReal

  integer,save :: nEnsemble=0

  type(particle), POINTER, save :: p11,p12,p21


  PUBLIC:: VolumeElements_INIT, VolumeElements_CLEAR, cleanUp
  PUBLIC:: VolumeElements_SETUP_Pert,VolumeElements_SETUP_Real
  PUBLIC:: VolumeElements_Statistics
  PUBLIC:: VolumeElements_GetPart, VolumeElements_InitGetPart
  PUBLIC:: VolumeElements_boxSize
  PUBLIC:: VolumeElements_NukSearch

  logical, save :: initFlag=.true.


contains
  !*************************************************************************
  !****f* VolumeElements/VolumeElements_boxSize
  ! PURPOSE
  ! Returns size of box used for the local ensemble method, unit of fm^3
  ! OUTPUT
  ! (function value)
  !*************************************************************************

  real function VolumeElements_boxSize()

    if(initFlag) call VolumeElements_INIT

    VolumeElements_boxSize=tVE%Delta(1)*tVE%Delta(2)*tVE%Delta(3)
  end function VolumeElements_boxSize


  subroutine cleanUp
    implicit none
    call VolumeElements_CLEAR()
    if (allocated(tVE%VE_real)) deallocate(tVE%VE_real)
    if (allocated(tVE%VE_pert)) deallocate(tVE%VE_pert)
    if (allocated(tVE%zCoordFilled_real)) deallocate(tVE%zCoordFilled_real)
    if (allocated(tVE%zCoordFilled_pert)) deallocate(tVE%zCoordFilled_pert)
  end subroutine


  !*************************************************************************
  !****s* VolumeElements/VolumeElements_INIT
  ! NAME
  ! subroutine VolumeElements_INIT()
  !
  ! PURPOSE
  ! This routine initializes the tVE-instance.
  ! Initial sizes are set and all memory allocation is done.
  !
  ! NOTES
  ! The maximum volume size and also the volume elements size is still
  ! hard wired. maybe some more sophisticated init should be realized.
  !*************************************************************************
  subroutine VolumeElements_INIT()
    implicit none

    integer :: i,j,k

    initFlag=.false. ! to remember that routine has been called

!    tVE%Delta = (/0.25,0.25,0.25/) ! x-, y-, z-Binning
    tVE%Delta = (/0.50,0.50,0.50/) ! x-, y-, z-Binning
!    tVE%Delta = (/1.00,1.00,1.00/) ! x-, y-, z-Binning

    tVE%Range(1,1:2) = (/-40., 40./) ! x-Range
    tVE%Range(2,1:2) = (/-40., 40./) ! y-Range
    tVE%Range(3,1:2) = (/-40., 40./) ! z-Range

    do i=1,3
       tVE%iRange(i,1:2) = nint(tVE%Range(i,1:2)/tVE%Delta(i))
       tVE%iRange(i,1) = tVE%iRange(i,1)-1 ! for security
       tVE%iRange(i,2) = tVE%iRange(i,2)+1 ! for security
    end do

    write(*,'(79("#"))')
    write(*,*) ' tVE-Size   = ',&
         & (tVE%iRange(1,2)-tVE%iRange(1,1))* &
         & (tVE%iRange(2,2)-tVE%iRange(2,1))* &
         & (tVE%iRange(3,2)-tVE%iRange(3,1))
    write(*,'(79("#"))')

    ALLOCATE(tVE%VE_real(tVE%iRange(1,1):tVE%iRange(1,2),&
         & tVE%iRange(2,1):tVE%iRange(2,2),&
         & tVE%iRange(3,1):tVE%iRange(3,2)))

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_INIT(tVE%VE_real(i,j,k))
          end do
       end do
    end do

    ALLOCATE(tVE%zCoordFilled_real(tVE%iRange(3,1):tVE%iRange(3,2)))

    tVE%zCoordFilled_real = .FALSE.

    ALLOCATE(tVE%VE_pert(tVE%iRange(1,1):tVE%iRange(1,2),&
         & tVE%iRange(2,1):tVE%iRange(2,2),&
         & tVE%iRange(3,1):tVE%iRange(3,2)))

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_INIT(tVE%VE_pert(i,j,k))
          end do
       end do
    end do

    ALLOCATE(tVE%zCoordFilled_pert(tVE%iRange(3,1):tVE%iRange(3,2)))

    tVE%zCoordFilled_pert = .FALSE.

  end subroutine VolumeElements_INIT

  !*************************************************************************
  !****s* VolumeElements/VolumeElements_CLEAR
  ! NAME
  ! subroutine VolumeElements_CLEAR()
  !
  ! PURPOSE
  ! This routine resets the tVE-instance.
  !
  !*************************************************************************
  subroutine VolumeElements_CLEAR()
    implicit none

    integer :: i,j,k

!    write(*,*) 'VolumeElements_CLEAR: ...'

    if(initFlag) call VolumeElements_INIT

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_real(k)) CYCLE
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_CLEAR(tVE%VE_real(i,j,k))
          end do
       end do
    end do
    tVE%zCoordFilled_real = .FALSE.

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_pert(k)) CYCLE
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             call ParticleList_CLEAR(tVE%VE_pert(i,j,k))
          end do
       end do
    end do
    tVE%zCoordFilled_pert = .FALSE.

!    write(*,*) 'VolumeElements_CLEAR: ok'

  end subroutine VolumeElements_CLEAR

  !*************************************************************************
  !****s* VolumeElements/VolumeElements_SETUP_Pert
  ! NAME
  ! subroutine VolumeElements_SETUP_Pert(PartVec)
  !
  ! PURPOSE
  ! Build up the tVE structure of perturbative particles.
  !*************************************************************************

  subroutine VolumeElements_SETUP_Pert(PartVec)
    use output
    implicit none

    type(particle), TARGET:: PartVec(:,:)

    integer :: i,j,k, ii(3)
    type(particle), POINTER :: pPart
    integer :: nPart

!    write(*,*) 'VolumeElements_SETUP_Pert: Start'

    if(initFlag) then
       call VolumeElements_Init
       call VolumeElements_Clear
    end if

    nPart = 0

    nEnsemble = Size(PartVec,dim=1)

    ensemble_loop : do i=1,Size(PartVec,dim=1)
       index_loop : do j=1,Size(PartVec,dim=2)
          If (PartVec(i,j)%ID < 0) exit index_loop
          If (PartVec(i,j)%ID == 0) cycle index_loop

          do k=1,3
             ii(k) = nint(PartVec(i,j)%position(k)/tVE%Delta(k))

             if ((ii(k)<=tVE%iRange(k,1)).or.(ii(k)>=tVE%iRange(k,2))) then
!                write(*,*) '#### Particle not in volume:',&
!                     & k,PartVec(i,j)%position(k)
                cycle index_loop
             end if
          end do
          pPart => PartVec(i,j)

          call ParticleList_APPEND(tVE%VE_Pert(ii(1),ii(2),ii(3)), pPart)
          nPart = nPart + 1

          tVE%zCoordFilled_pert(ii(3)) = .TRUE.

       end do index_loop
    end do ensemble_loop

!    write(*,*) '...particles inserted: ',nPart

    p11 => PartVec(1,1)
    if (Size(PartVec,dim=1)>1) then
      p21 => PartVec(2,1)
    else
      p21 => PartVec(1,1) ! dummy value
    end if
    p12 => PartVec(1,2)

    if (DoPR(1)) write(*,*) 'VolumeElements_SETUP_Pert: particles inserted: ',nPart

  end subroutine VolumeElements_SETUP_Pert

  !*************************************************************************
  !****s* VolumeElements/VolumeElements_SETUP_Real
  ! NAME
  ! subroutine VolumeElements_SETUP_Real(PartVec)
  !
  ! PURPOSE
  ! Build up the tVE structure of real particles.
  !*************************************************************************

  subroutine VolumeElements_SETUP_Real(PartVec)
    use output
    implicit none

    type(particle), TARGET:: PartVec(:,:)

    integer :: i,j,k, ii(3)
    type(particle), POINTER :: pPart
    integer :: nPart

!    write(*,*) 'VolumeElements_SETUP_Real: Start'

    if(initFlag) then
       call VolumeElements_Init
       call VolumeElements_Clear
    end if


    nPart = 0

    ensemble_loop : do i=1,Size(PartVec,dim=1)
       index_loop : do j=1,Size(PartVec,dim=2)
          If (PartVec(i,j)%ID < 0) exit index_loop
          If (PartVec(i,j)%ID == 0) cycle index_loop

          do k=1,3
             ii(k) = nint(PartVec(i,j)%position(k)/tVE%Delta(k))

             if ((ii(k)<=tVE%iRange(k,1)).or.(ii(k)>=tVE%iRange(k,2))) then
!                write(*,*) '#### Particle not in volume:',&
!                     & k,PartVec(i,j)%position(k)
                cycle index_loop
             end if
          end do
          pPart => PartVec(i,j)

          call ParticleList_APPEND(tVE%VE_Real(ii(1),ii(2),ii(3)), pPart)
          nPart = nPart + 1

          tVE%zCoordFilled_real(ii(3)) = .TRUE.

       end do index_loop
    end do ensemble_loop

!    write(*,*) '...particles inserted: ',nPart

    if (DoPR(1)) write(*,*) 'VolumeElements_SETUP_Real: particles inserted: ',nPart

  end subroutine VolumeElements_SETUP_Real


  !*************************************************************************
  !****s* VolumeElements/VolumeElements_Statistics
  ! NAME
  ! subroutine VolumeElements_Statistics
  !
  ! PURPOSE
  ! This is a routine to produce some statistical informations about
  ! the elements in the tVE instance.
  !
  ! This routine is only for trial/documentational purposes.
  !*************************************************************************
  subroutine VolumeElements_Statistics
    implicit none

    integer :: i,j,k, ii, hh

    integer, dimension(-1:201) :: histPartR, histPartP

    integer, dimension(-1:201) :: histCollRR, histCollPR

    hh = (tVE%iRange(1,2)-tVE%iRange(1,1)) * (tVE%iRange(2,2)-tVE%iRange(2,1))

    histPartR = 0
    histPartP = 0
    histCollRR = 0
    histCollPR = 0


    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_real(k)) then
          histPartR(-1) = histPartR(-1) + hh
          CYCLE
       endif
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = min(tVE%VE_real(i,j,k)%nEntries, 201)
             histPartR(ii) = histPartR(ii)+1
          end do
       end do
    end do

!    write(311,'(210(i9," "))') histPartR

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_pert(k)) then
          histPartP(-1) = histPartP(-1) + hh
          CYCLE
       endif
       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = min(tVE%VE_pert(i,j,k)%nEntries, 201)
             histPartP(ii) = histPartP(ii)+1
          end do
       end do
    end do

!    write(312,'(210(i9," "))') histPartP

    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_real(k)) CYCLE

       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = (tVE%VE_real(i,j,k)%nEntries * (tVE%VE_real(i,j,k)%nEntries - 1) )/2
             ii = max(min(ii,201),0)
             histCollRR(ii) = histCollRR(ii)+1
          end do
       end do
    end do

!    write(313,'(210(i9," "))') histCollRR


    do k=tVE%iRange(3,1),tVE%iRange(3,2)
       if (.not.tVE%zCoordFilled_pert(k)) CYCLE
       if (.not.tVE%zCoordFilled_real(k)) CYCLE

       do j=tVE%iRange(2,1),tVE%iRange(2,2)
          do i=tVE%iRange(1,1),tVE%iRange(1,2)
             ii = min(tVE%VE_pert(i,j,k)%nEntries * tVE%VE_real(i,j,k)%nEntries, 201)
             histCollPR(ii) = histCollPR(ii)+1
          end do
       end do
    end do

!    write(314,'(210(i9," "))') histCollPR

!    write(*,*) 'VolumeElements_Statistics: End'


  end subroutine VolumeElements_Statistics


  !*************************************************************************
  !****if* VolumeElements/FindNextVE
  ! NAME
  ! logical function FindNextVE
  !
  ! PURPOSE
  ! This routine searches for the next volume element with non-vanishing
  ! number of real- and perturbative-test-particles.
  ! It remembers the values of the last call and finds really the "next"
  ! cell.
  !
  ! It proceeds via first increasing the x- and y-coordinates and then
  ! stepping to the next z-coordinate.
  !
  ! INPUTS
  ! * the static stored array iPart(1:3)
  !
  ! OUTPUT
  ! * the static stored array iPart(1:3) (changed!)
  ! * function value: .false. -> no more tVE-cells possible
  !*************************************************************************

  logical function FindNextVE()
    implicit none

    do
       iPart(1)=iPart(1)+1
       if (iPart(1) > tVE%iRange(1,2)) then
          iPart(1) = tVE%iRange(1,1)
          iPart(2) = iPart(2)+1
          if (iPart(2) > tVE%iRange(2,2)) then
             iPart(2) = tVE%iRange(2,1)
             do
                iPart(3) = iPart(3)+1
                if (iPart(3) > tVE%iRange(3,2)) then
                   iPart(3) = tVE%iRange(3,1)
                   FindNextVE = .false.
                   return
                endif
                if (.not.tVE%zCoordFilled_pert(iPart(3))) cycle
                if (.not.tVE%zCoordFilled_real(iPart(3))) cycle
                exit
             end do

          endif
       endif

       if (tVE%VE_pert(iPart(1),iPart(2),iPart(3))%nEntries == 0) cycle
       if (tVE%VE_real(iPart(1),iPart(2),iPart(3))%nEntries == 0) cycle

       pReal => tVE%VE_real(iPart(1),iPart(2),iPart(3))%first
       pPert => tVE%VE_pert(iPart(1),iPart(2),iPart(3))%first

       FindNextVE = .true.
       exit
    end do


  end function FindNextVE

  !*************************************************************************
  !****s* VolumeElements/VolumeElements_InitGetPart
  ! NAME
  ! subroutine VolumeElements_InitGetPart
  !
  ! PURPOSE
  ! initialize the "GetPart"-routines
  !
  ! NOTES
  ! this routine sets the x-,y-,z-indizes in such a way, that a call to
  ! "FindNextVE" will start at the very first cell.
  !*************************************************************************
  subroutine VolumeElements_InitGetPart
    implicit none

    NULLIFY(pPert,pReal)
    iPart(1:2) = tVE%iRange(1:2,2) ! really start with max value
    iPart(3) = tVE%iRange(3,1)-1
    if (.not.FindNextVE()) then
       write(*,*) 'VolumeElements_InitGetPart: warning'
!       stop
    endif

  end subroutine VolumeElements_InitGetPart


#ifdef f95
#define LOCcmd pointer
#elif defined nagfor

  integer function myLOC(x)
    use, intrinsic :: iso_c_binding
    implicit none
    type(particle), pointer :: x
    type(c_ptr) :: cptr
    integer(c_intptr_t) :: ptr
    cptr = C_LOC(x)
    ptr = transfer(cptr,ptr)
    myLOC = ptr
  end function

#define LOCcmd myLOC

#else
#define LOCcmd LOC
#endif


  !*************************************************************************
  !****f* VolumeElements/VolumeElements_GetPart
  ! NAME
  ! logical function VolumeElements_GetPart(Part1, Part2, nRealPart, iEns,iInd)
  !
  ! INPUTS
  !
  ! OUTPUT
  ! * type(particle), POINTER :: Part1, Part2 -- real and perturbative particle
  ! * integer :: nRealPart -- number of real particles in VE-cell
  ! * integer :: iEns,iInd -- coordinates of perturbative particle
  !
  ! PURPOSE
  ! This routine finds the next (possibly colliding?) pair
  ! of one perturbative particle and one real particle
  ! in the actual VE-cell.
  !
  ! If one stepped over all pert. particles in thi given VE-cell, the next
  ! cell s choosen (cf. FindNextVE) and everything goes on.
  !
  ! If no "next VE-cell" is possible any more, the routine returns
  ! .false. as failure-indicator. (Otherwise always .true. is returned.)
  !
  !
  ! NOTES
  ! * actually for a choosen VE-cell it returns the particle pairs
  !     (P_1,R_1), (P_2,R_2), ... (P_nPert, R_xxx)
  !   If there are more real particles than pert particles,
  !   "R_xxx" stands for "R_nPert". Otherwise, ie. if we have more
  !   pert particles than real particles, the loop restarts
  !   for the real particles
  !     (P_nReal,R_nReal), (P_nReal+1,R_1), ...
  !
  ! * random selection of access is not implemented yet.
  !   one should check whether ist really necessary !
  !   (maybe the insertion into different tVE-cells offers already
  !   enough random choice into the whole game)
  !
  !*************************************************************************

  logical function VolumeElements_GetPart(Part1, Part2, nRealPart, iEns,iInd)
    implicit none
    type(particle)         , POINTER :: Part1, Part2
    integer :: nRealPart
    integer, intent(OUT) :: iEns,iInd
    integer :: j2
    !logical :: flag

    if (.not.ASSOCIATED(pPert)) then
       if (.not.FindNextVE()) then
         VolumeElements_GetPart = .false.
         return
      endif
    endif

    if (.not.ASSOCIATED(pReal)) then
       pReal => tVE%VE_real(iPart(1),iPart(2),iPart(3))%first
    endif

    Part1 => pReal%V
    Part2 => pPert%V
    nRealPart = tVE%VE_real(iPart(1),iPart(2),iPart(3))%nEntries

    j2  =int( (LOCcmd(Part2)-LOCcmd(p11))/(LOCcmd(p12)-LOCcmd(p11)))
    if (LOCcmd(p21) == LOCcmd(p11)) then
      iEns = 1
    else
      iEns = int((LOCcmd(Part2)-LOCcmd(p11)-j2*(LOCcmd(p12)-LOCcmd(p11)))/(LOCcmd(p21)-LOCcmd(p11))) + 1
    end if
    iInd=j2+1

    pPert => pPert%next
    pReal => pReal%next

    VolumeElements_GetPart = .true.
    return
  end function VolumeElements_GetPart

  !*************************************************************************
  !****s* VolumeElements/VolumeElements_NukSearch
  ! NAME
  ! subroutine VolumeElements_NukSearch(partIn,RadiusNukSearch,proton1,proton2,neutron1,neutron2,FlagOK)
  !
  ! PURPOSE
  ! This routine searches for two protons and two neutrons given in the
  ! "volume elements particle vector array" "VE_real" in the vicinty
  ! of the particle given by "partIn".
  !
  ! NOTES
  ! we are looking for 2 protons and 2 neutrons, i.e. for 4 nucleons, while
  ! only 2 nucleons are necessary: possible combinations are p+p, p+n, n+n.
  !
  ! INPUTS
  ! * type(particle)          :: partIn
  ! * real                    :: RadiusNukSearch
  !
  ! OUTPUT
  ! * type(particle), pointer :: proton1,proton2   -- Closest protons
  ! * type(particle), pointer :: neutron1,neutron2 -- Closest neutrons
  ! * logical                 :: FlagOK
  !
  !*************************************************************************
  subroutine VolumeElements_NukSearch(partIn,RadiusNukSearch,proton1,proton2,neutron1,neutron2,FlagOK)
    use GridOrdering
    use IDTable, only : nucleon
    use collisionNumbering, only : check_justCollided
    use constants, only : rhoNull
    use dichteDefinition
    use densityModule

    implicit none

    type(particle), intent(in) :: partIn
    real , intent(in) :: RadiusNukSearch
    type(particle), pointer :: proton1, proton2
    type(particle), pointer :: neutron1, neutron2
    logical, intent(out) :: FlagOK

    type(dichte) :: dens
    real,    save :: radiusMax
    logical, save :: DoInit_GridOrdering = .true.

    integer,dimension(0:1) :: nFound
    integer :: iRadius, iiRadius, nRadiusMax
    integer, save :: nRadius

    integer :: ii0(3), ii(3), k!,i

    real, dimension(1:3) :: position,DeltaPos

    type(particle), pointer :: pPart
    type(tParticleListNode), pointer :: pNode

    FlagOK = .false.
    nFound = 0


    if (DoInit_GridOrdering) then
       call GridOrdering_INIT
       nRadius = ubound(nDistance,dim=1)-1

       DoInit_GridOrdering = .false.
    end if

    position = partIn%position
    dens = densityAt(position)

    if(dens%baryon(0).lt.5e-03) then  !If density too small, then no absorption
!       write(*,*) '#### VolumeElements_NukSearch: density too small:',&
!            & dens%baryon,'#####',position
       return
    end if

    radiusMax = min(RadiusNukSearch*(rhoNull/dens%baryon(0))**(1./3.),5.0)
    radiusMax = radiusMax/(float(nEnsemble))**(1./3.)

    nRadiusMax = (nint(radiusMax/minval(tVE%Delta(1:3))))**2 + 1

!    write(*,*) 'radiusMax...:',radiusMax,minval(tVE%Delta(1:3)),nRadiusMax,nRadius


    do k=1,3
       ii0(k) = nint(position(k)/tVE%Delta(k))

       if ((ii0(k)<=tVE%iRange(k,1)).or.(ii0(k)>=tVE%iRange(k,2))) then
          write(*,*) '#### VolumeElements_NukSearch: Particle not in volume:', &
               & k,partIn%position(:)
          return
       end if
    end do




    do iRadius=0,min(nRadius,nRadiusMax) ! = 0,1,2,3 ,...
       if (nDistance(iRadius,2) < 0) cycle ! for this radius no cells...

!       write(*,*) '#### iRadius = ',iRadius

!       if (iRadius>0) write(*,*) '#### VolumeElements_NukSearch: iRadius>0',iRadius


       call GridOrdering_RandomizeRadius(iRadius)

       ! loop over all possibilities to get iRadius=const:
       do iiRadius=nDistance(iRadius,1),nDistance(iRadius,2)

          ii(1:3) = ii0(1:3) + DeltaV(iDeltaV(iiRadius),1:3)
          pNode => tVE%VE_Real(ii(1),ii(2),ii(3))%first

          do
             if (.not.ASSOCIATED(pNode)) exit
             pPart => pNode%V
             pNode => pNode%next

             if (pPart%ID.ne.nucleon) cycle
             if (pPart%antiparticle) cycle
             if (check_justCollided(PartIn,pPart)) cycle
             if ((pPart%Charge < 0).or.(pPart%Charge>1)) cycle
             if (nFound(pPart%Charge) >= 2) cycle

             DeltaPos = position-pPart%position
             if (DOT_PRODUCT(DeltaPos,DeltaPos) > radiusMax**2) cycle



             select case (pPart%Charge)
             case (0)
                select case (nFound(0))
                case (0)
                   neutron1 => pPart
                case (1)
                   neutron2 => pPart
                case default
                   cycle
                end select
                nFound(0) = nFound(0)+1

             case (1)
                select case (nFound(1))
                case (0)
                   proton1 => pPart
                case (1)
                   proton2 => pPart
                case default
                   cycle
                end select
                nFound(1) = nFound(1)+1

             end select


             if (nFound(0)==2 .and. nFound(1)==2) then
                FlagOK = .TRUE.
                return
             endif


          end do



       end do



    end do





  end subroutine VolumeElements_NukSearch
end module VolumeElements
