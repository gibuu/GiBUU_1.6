!*****************************************************************************
!****m* /BoxAnalysis
! NAME
! module BoxAnalysis
! PURPOSE
!*****************************************************************************
module BoxAnalysis

  use histMPf90

  implicit none

  Private

  type(histogramMP), save :: hMP_p
  real, dimension(:,:), allocatable, save :: arrQ2_all, arrQ2_formed

  PUBLIC :: DoBoxAnalysisTime

contains
  !*************************************************************************
  !****s* BoxAnalysis/DoBoxAnalysisTime
  ! NAME
  ! subroutine DoBoxAnalysisTime(realPart,timestep)
  ! PURPOSE
  !*************************************************************************
  subroutine DoBoxAnalysisTime(realPart,timestep)
    use particleDefinition
    use output, only: Write_InitStatus, intTochar

    type(particle),dimension(:,:),intent(in),   target :: realPart
    integer, intent(in) :: timestep

    logical,save :: startFlag=.true.
    type(particle), POINTER :: pPart

    integer :: nEns,nPart, i,j, iID
    integer, save :: nHist
    real :: wForm, Q2, mulfak

    if (startFlag) then
       call Write_InitStatus('BoxAnalysis',0)
       startFlag=.false.

       call CreateHistMP(hMP_p, "dN/dp", 0.0, 2.5, 0.05, 2)

       nHist = Map2HistMP_getN(2)
       allocate(arrQ2_all(0:nHist,2))
       allocate(arrQ2_formed(0:nHist,2))

       open(123,file="Q2_all_formed.dat", status="unknown")
       call WriteHistMP_Names(hMP_p,123)
       close(123)
       open(123,file="N_all_formed.dat", status="unknown")
       call WriteHistMP_Names(hMP_p,123)
       close(123)
       
       call Write_InitStatus('BoxAnalysis',1)
    end if

    call ClearHistMP(hMP_p)
    arrQ2_all = 0.0
    arrQ2_formed = 0.0

    nEns = size(realPart,dim=1)
    nPart = size(realPart,dim=2)

    mulfak = 1.0/nEns

    do i=1,nEns
       do j=1,nPart
          pPart => realPart(i,j)
          if(pPart%Id <  0) exit
          if(pPart%Id <= 0) cycle

          if (pPart%in_Formation) then
             wForm = 0.0
          else
             wForm = 1.0
          end if

          call AddHistMP(hMP_p, pPart, absMom(pPart), 1.0, wForm)

          Q2 = 2*pPart%momentum(3)**2 - pPart%momentum(1)**2 - pPart%momentum(2)**2

          arrQ2_all(0,1) = arrQ2_all(0,1) + 1.0
          arrQ2_all(0,2) = arrQ2_all(0,2) + Q2

          arrQ2_formed(0,1) = arrQ2_formed(0,1) + wForm
          arrQ2_formed(0,2) = arrQ2_formed(0,2) + Q2*wForm

          iID = Map2HistMP_ID(pPart%ID,pPart%charge,pPart%antiparticle, 2)
          if (iID>0) then
             arrQ2_all(iID,1) = arrQ2_all(iID,1) + 1.0
             arrQ2_all(iID,2) = arrQ2_all(iID,2) + Q2

             arrQ2_formed(iID,1) = arrQ2_formed(iID,1) + wForm
             arrQ2_formed(iID,2) = arrQ2_formed(iID,2) + Q2*wForm
          end if

       end do
    end do

    if (mod(timestep,5)==1) then
       call WriteHistMP(hMP_p, file='p_all_'//intTochar(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=1)
       call WriteHistMP(hMP_p, file='p_formed_'//intTochar(timestep)//'.dat', add=1e-20, mul=mulfak, iColumn=3)
    end if

    
    open(123,file="N_all_formed.dat",  status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
         & arrQ2_all(1:nHist,1)*mulfak,arrQ2_all(0,1)*mulfak, &
         & arrQ2_formed(1:nHist,1)*mulfak,arrQ2_formed(0,1)*mulfak
    close(123)

    do i=0,nHist
       if (arrQ2_all(i,1)>0.0)     &
            & arrQ2_all(i,2)=arrQ2_all(i,2)/arrQ2_all(i,1)
       if (arrQ2_formed(i,1)>0.0)  &
            & arrQ2_formed(i,2)=arrQ2_formed(i,2)/arrQ2_formed(i,1)
    end do

    open(123,file="Q2_all_formed.dat", status="old",position='append')
    write(123,'(i11,1P,100E12.4,0P)') timestep, &
         & arrQ2_all(1:nHist,2),arrQ2_all(0,2),&
         & arrQ2_formed(1:nHist,2),arrQ2_formed(0,2)
    close(123)


  end subroutine DoBoxAnalysisTime

end module BoxAnalysis
