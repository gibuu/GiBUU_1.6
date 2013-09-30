!*****************************************************************************
!****m* /PionBoxAnalysis
! NAME
! module BoxAnalysis
! PURPOSE
!*****************************************************************************
module PionBoxAnalysis

  use histf90
  use histMPf90
  use histMC

  implicit none

  Private

  PUBLIC :: DoPionBoxAnalysisTime

  type(histogram), save :: hMassRho
  type(histogramMC), save :: hMCMomPion
  type(histogramMP), save :: hMP_pSet2, hMP_pSet4
  real, dimension(:,:), allocatable, save :: arrMultSet2, arrMultSet4

contains
  !*************************************************************************
  !****s* PionBoxAnalysis/DoBoxAnalysisTime
  ! NAME
  ! subroutine DoPionBoxAnalysisTime(realPart,timestep)
  ! PURPOSE
  !*************************************************************************
  subroutine DoPionBoxAnalysisTime(realPart,timestep)

    use particleDefinition
    use output, only: Write_InitStatus, intTochar4, WriteParticleVector
    use history, only: history_getParents

    type(particle),dimension(:,:),intent(in),   target :: realPart
    integer, intent(in) :: timestep

    logical,save :: startFlag=.true.
    integer,save :: nHist

    integer :: nEns,nPart, i,j, iID, iCh
    real :: mulfak, mom, mass
    type(particle), POINTER :: pPart
    integer :: parents(1:3)

    if (startFlag) then
       call Write_InitStatus('PionBoxAnalysis',0)
       startFlag=.false.

       call CreateHist(hMassRho, "mass(rho)", 0.0, 1.2, 0.01)
       call CreateHistMC(hMCMomPion, "momentum(pion)", 0.0, 2.5, 0.02, 6)
       hMCMomPion%yDesc(1:6) = (/ "original  ",  &
            "rho       ", "sigma     ", "other dec ", &
            "pi pi     ", "other coll" /)

       call CreateHistMP(hMP_pSet2, "dN/p^2 dp", 0.0, 2.5, 0.05, 2)
       call CreateHistMP(hMP_pSet4, "dN/p^2 dp", 0.0, 2.5, 0.05, 4)

       nHist = Map2HistMP_getN(2)
       allocate( arrMultSet2(0:nHist, 2) )
       nHist = Map2HistMP_getN(4)
       allocate( arrMultSet4(0:nHist, 2) )

       open(123,file="PionBoxAnalysis_Mult_Set2.dat", status="unknown")
       call WriteHistMP_Names(hMP_pSet2,123)
       close(123)

       open(123,file="PionBoxAnalysis_Mult_Set4.dat", status="unknown")
       call WriteHistMP_Names(hMP_pSet4,123)
       close(123)

       call Write_InitStatus('PionBoxAnalysis',1)
    end if

    call ClearHistMP(hMP_pSet2)
    call ClearHistMP(hMP_pSet4)
    arrMultSet2 = 0.0
    arrMultSet4 = 0.0

    call ClearHistMC(hMCMomPion)

    nEns = size(realPart,dim=1)
    nPart = size(realPart,dim=2)

    mulfak = 1.0/nEns

    do i=1,nEns
       do j=1,nPart
          pPart => realPart(i,j)
          if(pPart%Id <  0) exit
          if(pPart%Id <= 0) cycle

          mom = absMom(pPart)

          select case(pPart%ID)
          case (101)
             parents = history_getParents(pPart%history)
!             write(*,*) parents
             if (parents(2) == 0) then
                select case(parents(1))
                case (0)
                   iCh = 1
                case (103)
                   iCh = 2
                case (104)
                   iCh = 3
                case default
                   iCh = 4
                end select
             else 
                if (parents(1)==101 .and. parents(2)==101) then
                   iCh = 5
                else
                   iCh = 6
                end if
             end if
             call AddHistMC(hMCMomPion, mom, iCh, 1.0/(mom**2))
          case (103)
             mass = sqrtS(pPart)
             call AddHist(hMassRho, mass, 1.0)
          end select

          call AddHistMP(hMP_pSet2, pPart, mom, 1.0/(mom**2), 1.0)
          call AddHistMP(hMP_pSet4, pPart, mom, 1.0/(mom**2), 1.0)

          arrMultSet2(0,1) = arrMultSet2(0,1) + 1.0
          arrMultSet4(0,1) = arrMultSet4(0,1) + 1.0

          iID = Map2HistMP_ID(pPart%ID,pPart%charge,pPart%antiparticle, 2)
          if (iID>0) then
             arrMultSet2(iID,1) = arrMultSet2(iID,1) + 1.0
          endif
          iID = Map2HistMP_ID(pPart%ID,pPart%charge,pPart%antiparticle, 4)
          if (iID>0) then
             arrMultSet4(iID,1) = arrMultSet4(iID,1) + 1.0
          endif
       end do
    end do

    if (mod(timestep,5)==1) then
       call WriteHistMP(hMP_pSet2, file='p_Set2_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=1)
       call WriteHistMP(hMP_pSet4, file='p_Set4_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=1)
       call WriteHist(hMassRho, file='massRho_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)

!       call WriteParticleVector('parts_'//intTochar(timestep),realPart)

       call WriteHistMC(hMCMomPion, file='MomPion_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
    end if

    open(123,file="PionBoxAnalysis_Mult_Set2.dat",status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
         & arrMultSet2(1:,1)*mulfak,arrMultSet2(0,1)*mulfak
    close(123)

    open(123,file="PionBoxAnalysis_Mult_Set4.dat",status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
         & arrMultSet4(1:,1)*mulfak,arrMultSet4(0,1)*mulfak
    close(123)
    

  end subroutine DoPionBoxAnalysisTime

end module PionBoxAnalysis
