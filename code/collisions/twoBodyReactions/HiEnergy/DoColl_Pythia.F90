!***************************************************************************
!****m* /Coll_Pythia
! NAME
! module Coll_Pythia
! FUNCTION
! Implement a + b -> X processes/reactions done by PYTHIA
!
! NOTES
! There is a PreProcessor statement, which steers the weighting factor.
!***************************************************************************

#define PreWeight 2

module Coll_Pythia

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DoColl_Pythia

contains


  !*************************************************************************
  !****s* Coll_Pythia/DoColl_Pythia
  ! NAME
  ! subroutine DoColl_Pythia (inPart, outPart, flagOK, sqrtS, pcm, beta, DoReweight, weight)
  !
  ! PURPOSE
  ! Perform a collision of particles given in "inPart" with energy "sqrtS" and
  ! return outgoing particles in "outPart".
  ! "pcm" and "beta" are vectors used for Boost and Rotation of the event.
  ! if "flagOK" is false, no event happened, the output in "outPart" should 
  ! be neglected!
  ! You can call Pythia in the "Reweighting mode" by setting "DoReweight" to 
  ! true. Then "weight" has to be respected.
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: inPart   -- incoming particles
  ! * real                        :: sqrtS    -- energy of ollision
  ! * real, dimension(0:3)        :: pcm      -- boost-vector
  ! * real, dimension(1:3)        :: beta     -- boost vector
  ! * logical                     :: DoReweight -- Flag
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  ! outgoing particles
  ! * real                        :: weight   ! weight of event (default:1.0)
  ! * logical                     :: flagOK   ! event okay ?
  !
  ! NOTES
  ! in order to understand the meaning of "pcm" and "beta": 
  ! The (Pythia-)event is done in the restframe of the two given particles. 
  ! Then a call to PYROBO according 
  !       phi = atan2(pcm(2),pcm(1))
  !       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
  !       call PYROBO(1,N, theta,phi, beta(1),beta(2),beta(3))
  ! is performed in order to transform the system into the desired 
  ! (Lab-) system.
  !*************************************************************************
  subroutine DoColl_Pythia (inPart, outPart, flagOK, sqrtS, pcm, beta, DoReweight, weight)
    use particleDefinition, only: particle, IsSamePart
    use IdTable, only: nucleon, Delta
    use CollTools
    use twoBodyTools, only : IsChargeExchange 
    use hadronFormation, only: useJetSetVec
    use output, only: DoPR, WriteParticle
    use CollGetLeading, only: GetLeading_PY
    use ID_translation, only: KFfromBUU

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    real,                       intent(in)   :: sqrtS
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    logical,                    intent(in)   :: DoReweight
    real,                       intent(out)  :: weight
    logical,                    intent(out)  :: flagOK

    integer :: iTry, i
    real :: theta, phi
    integer :: ID1,ID2, IZ1,IZ2
    integer :: iz1c, iz2c        ! reduced charges
    integer :: id1c, id2c        ! reduced particle IDs
    integer :: kf1, kf2          ! KF-codes of incoming particles
    integer :: DeltaQ!,DeltaQ0    ! units of charges to be added
    logical :: EventIsAnti       ! particles converted to its anti ?
    
    type(particle) :: part1, part2 ! incoming particles

    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
    integer MDCY,MDME,KFDP
    double precision BRAT
    SAVE /PYDAT3/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/
    
    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/


    character*20 Buf,cBeam,cTarget

    integer :: iRewCl                     ! class of reweighting
    real ::  weight1(3,3),weight2(3,3)

    real :: EA,EB,pA

!---------------------------------------------------------------------

#if PreWeight == 2
!--- weight=2 (sqrt(s) = 10 ... 40 GeV) ------------------------------

    data (weight1(i,1),i=1,3) /-3.643d-3,-3.178d-4, 6.701d-5/ ! N+N 
    data (weight1(i,2),i=1,3) /-1.445d-2, 7.927d-4, 7.777d-5/ ! pi+N 
    data (weight1(i,3),i=1,3) /-1.471d-2, 8.118d-4, 7.768d-5/ ! K+N
    
    data (weight2(i,1),i=1,3) / 2.432d0 , 4.646d-2,-4.939d-4/ ! N+N 
    data (weight2(i,2),i=1,3) / 2.521d0 , 4.857d-2,-5.307d-4/ ! pi+N
    data (weight2(i,3),i=1,3) / 2.505d0 , 4.980d-2,-5.563d-4/ ! K+N

#elif PreWeight == 4
!--- weight=4 (sqrt(s) = 10 ... 40 GeV) ------------------------------


    data (weight1(i,1),i=1,3) /-3.683d-3,-3.134d-4, 6.691d-5/ ! N+N 
    data (weight1(i,2),i=1,3) /-1.449d-2, 7.963d-4, 7.769d-5/ ! pi+N
    data (weight1(i,3),i=1,3) /-1.511d-2, 8.476d-4, 7.704d-5/ ! K+N
                                                                    
    data (weight2(i,1),i=1,3) / 5.015d0 , 3.737d-1,-2.631d-3/ ! N+N 
    data (weight2(i,2),i=1,3) / 5.181d0 , 4.422d-1,-3.269d-3/ ! pi+N
    data (weight2(i,3),i=1,3) / 5.124d0 , 4.408d-1,-3.221d-3/ ! K+N

#elif PreWeight == 6
!--- weight=6 (sqrt(s) = 10 ... 40 GeV) ------------------------------

    data (weight1(i,1),i=1,3) /-3.922d-3,-2.908d-4, 6.660d-5/ ! N+N 
    data (weight1(i,2),i=1,3) /-1.482d-2, 8.302d-4, 7.755d-5/ ! pi+N
    data (weight1(i,3),i=1,3) /-1.386d-2, 7.198d-4, 7.996d-5/ ! K+N
                                                                    
    data (weight2(i,1),i=1,3) / 1.020d1 , 1.405d0 , 3.276d-2/ ! N+N 
    data (weight2(i,2),i=1,3) / 6.879d0 , 2.031d0 , 4.451d-2/ ! pi+N
    data (weight2(i,3),i=1,3) / 6.315d0 , 2.081d0 , 4.246d-2/ ! K+N
#else
!    write(*,*) 'PreWeight =',PreWeight
    write(*,*) 'wrong PreWeight'
    stop
#endif
!---------------------------------------------------------------------



    !...set default output

    outPart%ID = 0 ! reset outgoing particles
    weight = 1.0   ! default value
    flagOK = .FALSE.
    EventIsAnti = .FALSE.


    !...set incoming particles:
     if(inPart(1)%antiparticle) then ! (Shift antiparticle to second position)
        part1 = inPart(2)
        part2 = inPart(1)
     else
        part1 = inPart(1)
        part2 = inPart(2)
     endif

     if (part1%antiparticle) then
        if (DoPr(1)) write(*,*) 'DoColl_Pythia: converting anti+anti: '
        call WriteParticle(6,0,1,inPart(1))
        call WriteParticle(6,0,2,inPart(2))
        
        call ConvertToAnti(part1)
        call ConvertToAnti(part2)
        EventIsAnti = .TRUE.

        call WriteParticle(6,0,1,part1)
        call WriteParticle(6,0,2,part2)
        stop
     endif
        

    !...no collision for baryonic resonances Xi,Omega,Sigma_c,Lambda_c,... for sqrts<3:

    if (sqrtS.lt.3.0) then
       if (     ((part1%ID.lt.100).and.(part1%ID.gt.52)) &
            .or.((part2%ID.lt.100).and.(part2%ID.gt.52))) then
          if (DoPr(2)) write(*,'(A)') ' no PYTHIA collision for baryonic resonances Xi,Omega,Sigma_c,Lambda_c... for sqrts<3!'
          return
       endif
    endif
    
!!$    call WriteParticle(6,0,1,part1)
!!$    call WriteParticle(6,0,2,part2)
!!$    write(*,*) 'Sqrts:',sqrts
!!$    write(*,*) 'PCM:  ',PCM
!!$    write(*,*) 'Beta: ',BETA
    

100 continue

    !...Set ID etc

    ID1 = part1%ID
    ID2 = part2%ID
    if (part1%antiparticle) ID1 = -ID1
    if (part2%antiparticle) ID2 = -ID2
    IZ1 = part1%charge
    IZ2 = part2%charge

    !...reduce charge:      

    call ReduceCharge(ID1,IZ1,IZ1c)
    call ReduceCharge(ID2,IZ2,IZ2c)
    
    DeltaQ = (iz1+iz2)-(iz1c+iz2c)
    
    !...convert input particles:
    
    call ConvertInPartPythia(ID1,ID1c)
    call ConvertInPartPythia(ID2,ID2c)

    if (ID1c == 0) return
    if (ID2c == 0) return


    if (ID2c < -31) then
!!$       write(*,*) 'DoColl_Pythia: no "X+AntiBaryon_s" collisions !!!'
!!$       call WriteParticle(6,0,1,part1)
!!$       call WriteParticle(6,0,2,part2)
!!$       write(*,*) 'Trying to do as   "AntiX+Baryon_s"          0 !!!'
!!$       write(*,*)
       if (DoPr(1)) write(*,*) 'DoColl_Pythia: Retry "X+AntiBaryon_s" as ANTI!!!'

       if (EventIsAnti) then
          if (DoPr(2)) write(*,*) 'DoColl_Pythia: Event is already ANTI. QUIT!'
          return
       endif

       if ((ID1c .eq. 110).or.(ID1c .eq. 112)) then
          if (DoPr(1)) write(*,*) 'DoColl_Pythia: ANTI not possible for "K+AntiBaryon_s". QUIT'
          return
       endif


       if(inPart(1)%antiparticle) then
          call ConvertToAnti(inPart(1),part1) ! converted particle is anti,
          call ConvertToAnti(inPart(2),part2) ! should be stored at pos 2
       else
          call ConvertToAnti(inPart(2),part1)
          call ConvertToAnti(inPart(1),part2)
       endif
       EventIsAnti = .TRUE.
       goto 100! LabelReDo
    endif
       

    !...transpose BUU-code -> KF:
    
    KF1 = KFfromBUU (ID1c,IZ1c)
    KF2 = KFfromBUU (ID2c,IZ2c)
    
    !... if Reweighting: Get 'Class' (='pi+N','N+N')

    if (doReweight) then
       iRewCl = GetRewClass( (/ID1c,ID2c/) )
       if (iRewCl==0) then
          write(*,*) 'iRewCl=0 for ',ID1c,ID2c, sqrts
!          stop
          return
       endif
    endif

    !...set up PYTHIA/JETSET:

    call SetSomeDefaults_PY

    if (doReweight) then
       CKIN(3) = 1.5          ! min pThat 
       MSTP(142) = 1          ! use reweighted events
    else
       CKIN(3) = 0d0          ! min pThat (allow low pT events)
       MSTP(142) = 0          ! use reweighted events
    end if

    !... Initialize PYTHIA

    call PYNAME(KF1,Buf)
    write(cBeam,'(A)') Buf(1:10)
    call PYNAME(KF2,Buf)
    write(cTarget,'(A)') Buf(1:10)


    ! K^0 -> K_L^0:
    if (KF1==311) cBeam   = 'KL0'
    if (KF2==311) cTarget = 'KL0'


    if ((KF1==-311).or.(KF2==-311)) then

!!$       write(*,*) 'AUAAHHHH !!! (es gibt aber kein KL0bar!!!)'
!!$       if (KF1==-311) cBeam   = 'KL0'
!!$       if (KF2==-311) cTarget = 'KL0'

       if (EventIsAnti) then
          if (DoPr(2)) write(*,*) 'DoColl_Pythia: Event is already ANTI. QUIT!'
          return
       endif
       if (DoPr(1)) write(*,*) 'DoColl_Pythia: Retry "K0bar+X" as ANTI!!!'


       if (KF1==-311) then
          call ConvertToAnti(inPart(1),part1) ! converted baryon is anti,
          call ConvertToAnti(inPart(2),part2) ! should be stored at pos 2
       else
          call ConvertToAnti(inPart(1),part2)
          call ConvertToAnti(inPart(2),part1)
       endif

       EventIsAnti = .TRUE.
       goto 100! LabelReDo


    endif

!    write(*,*) '--',cBeam,cTarget, sqrts
!    call WriteParticle(6,0,1,inPart(1))
!    call WriteParticle(6,0,2,inPart(2))
   

    call CalcMomentum(sqrts,part1%mass,part2%mass, EA,EB,pA)

    P(1,1) = 0d0
    P(1,2) = 0d0
    P(1,3) = pA
    P(1,4) = EA
    P(1,5) = part1%mass

    P(2,1) = 0d0
    P(2,2) = 0d0
    P(2,3) = -pA
    P(2,4) = EB
    P(2,5) = part2%mass

    
!    call PYINIT('CMS', cBeam,cTarget, sqrts)
    call PYINIT('5MOM', cBeam,cTarget, sqrts)

    !... Start the Event Loop
    
    iTry = 0
    do
       outPart%ID = 0 ! reset outgoing particles
       
       iTry=iTry+1
       if(iTry.ge.100) then
          write(*,*)'DoColl_Pythia: itry=',iTry
          return
       endif

       if (useJetSetVec) call GetJetsetVecINIT
       
       !...Generate THE EVENT:
       
       call PYEVNT
!       call PYLIST(2)
!       stop

       if (MINT(51)==2) then
          N = 0
          ! print *,"failure in DoColl_Pythia!"
          return ! -> FAILURE
       endif

       if (useJetSetVec) then
          call GetJetsetVec(.TRUE.)
!          call PYLIST(2)
!          call GetJetSetVec_List(6,1,N)

          call GetJetsetVecCheckT(-1d-5)

          call GetJetsetVecPYEDIT
       endif
       
       call GetLeading_PY         ! find leading particles
       call PYEDIT(1)             ! clean up event list
       
       !...Rotate and boost the whole event to final system:
       
       phi = atan2(pcm(2),pcm(1))
       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
       call PYROBO(1,N, theta,phi, beta(1),beta(2),beta(3))

       if (useJetSetVec) then
          call GetJetsetVecPYROBO(theta,phi, beta(1),beta(2),beta(3))

!          call PYLIST(2)
!          call GetJetSetVec_List(6,1,N)
!          stop
       end if

       
       !...Copy Particles to ouput-vector
       
       call SetVectorFromPYJETS(outPart, real(max(1d0,VINT(52))) )

       ! If MSTI(1)==95 it was a low-pT-scattering and VINT(52)=0
       ! Otherwise VINT(52) holds something like
       !   Q2 = pT_hat^2 + (m_3^2+m_4^2)/2
       ! See PYTHIA manulal for details.
       !
       ! 1d0 is set as a (dummy) minimal value !!!!

       !...Correct Charges
       
!       DeltaQ0 = DeltaQ
       if (DeltaQ.ne.0) call CorrectChargeVector(outPart,DeltaQ)
       if (DeltaQ.ne.0) then
          if (DoPr(2)) then
             write(*,*) 'DoColl_Pythia: Charge correction failed. ReDo Event!!'
!             write(*,*) 'DeltaQ : ',DeltaQ0,DeltaQ
!             call PYLIST(2)
          endif

          cycle
       endif

       !...correct antiparticles if necessary 

       if (EventIsAnti) then
!      write(*,*) 'DoColl_Pythia: Event was done as ANTI'
          call ConvertToAnti(outPart)
       endif

       if (N.eq.2) then

          !...Test for elastic event

          if (IsSamePart(InPart(1),OutPart(1))&
               &.and.(IsSamePart(InPart(2),OutPart(2)))) cycle
          if (IsSamePart(InPart(1),OutPart(2))&
               &.and.(IsSamePart(InPart(2),OutPart(1)))) cycle


          !...Test for charge exchange event for antinucleon-nucleon, antinucleon-delta
          !...or nucleon-antidelta collisions

          if(InPart(1)%Id+InPart(2)%Id.le.nucleon+delta .and.&
            &(InPart(1)%antiParticle.neqv.InPart(2)%antiParticle) .and.&
            &IsChargeExchange(InPart(1),InPart(2),OutPart(1),OutPart(2))) cycle

       end if

       !...exit the loop
       exit
       
    end do
    
   
    !...Set weight:
   
!!$     ! PARI(1): sigma_tot in mb
!!$     ! PARI(2): unweighted: =PARI(1)/MSTI(5)
!!$     !          weighted:   =PARI(1)/sum(PARI(10))
!!$     ! PARI(10): compensating weight of event                               
!!$     ! MSTI(5): number of generated events
!!$     ! -> ave(PARI(10)) = PARI(1)/(PARI(2)*MSTI(5))
!!$      
!!$      weight = PARI(10) * MSTI(5)*PARI(2)/PARI(1) !=1.0 for unweighted events!
!!$
!!$ in the arrays: 
!!$   weight2: =1/ave(PARI(10)) 
!!$   weight1: = PARI(1)|(CKIN(3)=1.5) / PARI(1)|(CKIN(3)=0)
!!$ you can get these numbers by fitting the output of a routine like 
!!$ GetAveWeight2 (see below) 
   
   
   !... correct the weight:
   if (doReweight) then
      weight = PARI(10) * Polym( weight2(1:3,iRewCl), sqrts)
   endif

   if (CKIN(3).gt.0d0) then
      weight = weight * Polym( weight1(1:3,iRewCl), sqrts)
   endif

   flagOK = .TRUE.





!   call PYGIVE('MSTI(1)=')

   !  11 f + f' -> f + f' (QCD)
   !  12 f + fbar -> f' + fbar'
   !  13 f + fbar -> g + g     
   !  28 f + g -> f + g        
   !  53 g + g -> f + fbar     
   !  68 g + g -> g + g        
   !  95 Low-pT scattering     

!   call PYGIVE('MSTI(1)=')
!   call PYGIVE('VINT(52)=')
!   call PYLIST(2)
!   write(*,*) '===================='

  end subroutine DoColl_Pythia


  !*************************************************************************
  !****if* Coll_Pythia/Polym
  ! NAME
  ! function Polym(facs,x)
  !
  ! PURPOSE
  ! Calculate the value of a polynomial. 
  !
  ! INPUTS
  ! * real, dimension(:) :: facs -- coefficients
  ! * real               :: x    -- value
  !
  ! OUTPUT
  ! function value -- f(x)
  !
  ! NOTES
  ! This is very elementary. Check code in order not to duplicate routines!
  !*************************************************************************

  function Polym(facs,x)
    real                           :: Polym
    real, dimension(:), intent(in) :: facs
    real,               intent(in) :: x

    Polym = 0

    select case(size(facs))
    case(0)
       return
    case(1)
       Polym = facs(1)
    case(2)
       Polym = facs(1)+facs(2)*x
    case(3)
       Polym = facs(1)+(facs(2)+facs(3)*x)*x
       
    end select

  end function Polym


  !*************************************************************************
  !****is* Coll_Pythia/CalcMomentum
  ! NAME
  ! subroutine CalcMomentum(sqrts,mA,mB, EA,EB,pA)
  !
  ! PURPOSE
  ! calculate the momentum in a two particle system (fundamental kinematics)
  !
  ! INPUTS
  ! * sqrts,mA,mB, EA,EB -- cm energy, masses, energies
  !
  ! OUTPUT
  ! * pA -- momentum of particle A
  !
  ! NOTES
  ! This is very elementary. Check code in order not to duplicate routines!
  !*************************************************************************

  subroutine CalcMomentum(sqrts,mA,mB, EA,EB,pA)
    real, intent(in)  :: sqrts,mA,mB
    real, intent(out) :: EA,EB,pA

    real :: mA2,mB2,s

    s = sqrts**2
    mA2 = mA**2
    mB2 = mB**2
    
    EA = (s+(mA2-mB2))/(2*sqrts)
    EB = (s-(mA2-mB2))/(2*sqrts)
    pA = sqrt((s-mA2-mB2)**2-4*mA2*mB2)/(2*sqrts)

  end subroutine CalcMomentum


 !***************************************************************************

  integer function GetRewClass(ID)
    integer, intent(in), dimension(2) :: ID

    integer, dimension(2) :: II
    integer :: i

    GetRewClass = 0
    II = ID

    ! for simplicity: Lambda,Sigma,...reweighting like Nucleon 

    do i=1,2
       if (II(i)<56) II(i) = 1
    end do

    if ((II(1)==1).and.(II(2)==1)) then
       GetRewClass = 1
    else if (II(1)==1) then
       select case(II(2))
       case(101)
          GetRewClass = 2
       case(110,111)
          GetRewClass = 3
       end select
    else if (II(2)==1) then
       select case(II(1))
       case(101)
          GetRewClass = 2
       case(110,111)
          GetRewClass = 3
       end select
    else
    end if
    return
  end function GetRewClass

 !***************************************************************************

end module Coll_Pythia

!================================================================= 
! called by PYTHIA while using reweighted events

      SUBROUTINE PYEVWT(WTXS)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      COMMON/PYINT1/MINT(400),VINT(400)
      integer MINT
      double precision VINT
      SAVE /PYINT1/

!      write(*,*) 'PYEVWT called'


      pT2 = VINT(48)

#if PreWeight == 2
      WTXS = pT2    ! w2
#elif PreWeight == 4
      WTXS = pT2**2 ! w4
#elif PreWeight == 6
      WTXS = pT2**3 ! w6
#else
    write(*,*) 'PreWeight =',PreWeight
    stop
#endif
      return
      end

!================================================================= 


!!$c================================================================= 
!!$
!!$c as an example of how to get the weights for reweighted events:
!!$c the output (fort.85) has to be fitted with a polynomial (degree=2)
!!$
!!$      subroutine GetAveWeight2
!!$
!!$      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
!!$      integer MSTP,MSTI
!!$      double precision PARP,PARI
!!$      SAVE /PYPARS/
!!$
!!$      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
!!$      integer MSTU,MSTJ
!!$      double precision PARU,PARJ
!!$      SAVE /PYDAT1/
!!$
!!$      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
!!$      integer KCHG
!!$      double precision PMAS,PARF,VCKM
!!$      SAVE /PYDAT2/
!!$
!!$      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
!!$      integer MDCY,MDME,KFDP
!!$      double precision BRAT
!!$      SAVE /PYDAT3/ 
!!$
!!$      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
!!$      integer MSEL,MSELPD,MSUB,KFIN
!!$      double precision CKIN
!!$      SAVE /PYSUBS/
!!$
!!$      double precision results(4,0:100)
!!$
!!$      integer isqrts,nsqrts
!!$      double precision sqrts,sqrts1,sqrts2
!!$      integer i
!!$
!!$      parameter (sqrts1=10d0,sqrts2=40d0,nsqrts=30)
!!$
!!$      character*10 cPart1,cPart2
!!$      data cPart1,cPart2 /'pi+','p'/
!!$
!!$      integer iEV,nEV
!!$c      parameter (NEV=10000000) !
!!$c      parameter (NEV=1000000) ! 9m 12s (w=4), 4m 50s (w=0)
!!$      parameter (NEV=100000) ! 0m 46s
!!$c      parameter (NEV=10000) ! 
!!$c      parameter (NEV=100)
!!$
!!$      integer iH
!!$      double precision hhh
!!$
!!$     
!!$      MSTU(16) = 2              ! needed by GetLeading
!!$      PARP(2) = 1d0             ! lowest c.m. energy
!!$      CKIN(3) = 0d0             ! min pThat (allow low pT events)
!!$      PARP(111) = 0d0           ! otherwise problems at low sqrts
!!$      MSTP(111) = 0             ! master switch fragmentation/decay
!!$c      MSTP(122) = 0             ! switch init and max XS print-out
!!$
!!$      MSTJ(21) = 0              ! particle decay on/off
!!$      do i=556,560              ! switch off rho_0-decay-channels
!!$         MDME(i,1) = 0          ! (works only with v6.206)
!!$      enddo                     ! necessary, but should not be
!!$
!!$                                ! BRUTE FORCE !!!!
!!$      PMAS(102,1) = 0.138d0     ! set pi0-mass to BUU-value
!!$      PMAS(106,1) = 0.138d0     ! set pi+-mass to BUU-value
!!$
!!$      do isqrts=0,nsqrts
!!$
!!$         sqrts= sqrts1+(sqrts2-sqrts1)*isqrts/nsqrts
!!$
!!$         results(1,isqrts) = sqrts
!!$
!!$
!!$c************ ORIGINAL SETTING:
!!$
!!$         CKIN(3)=0d0            ! pThat_min
!!$         MSTP(142)=0            ! use reweighted events
!!$         CALL PYINIT('CMS',cPart1,cPart2,sqrts)
!!$         do iEV=1,NEV
!!$            CALL PYEVNT
!!$         enddo
!!$         results(2,isqrts) = PARI(1)
!!$
!!$c************ pThat_min > 0:
!!$
!!$         CKIN(3)=1.5d0          ! pThat_min
!!$         MSTP(142)=0            ! use reweighted events
!!$         CALL PYINIT('CMS',cPart1,cPart2,sqrts)
!!$         do iEV=1,NEV
!!$            CALL PYEVNT
!!$         enddo
!!$         results(3,isqrts) = PARI(1)
!!$
!!$c************ reweigthed:
!!$
!!$         CKIN(3)=1.5d0          ! pThat_min
!!$         MSTP(142)=1            ! use reweighted events
!!$         CALL PYINIT('CMS',cPart1,cPart2,sqrts)
!!$         do iEV=1,NEV
!!$            CALL PYEVNT
!!$         enddo
!!$         results(4,isqrts) = MSTI(5)*PARI(2)/PARI(1)
!!$        
!!$
!!$
!!$         write(*,*) results(1,isqrts),
!!$     $        results(3,isqrts)/results(2,isqrts),
!!$     $        results(4,isqrts)
!!$
!!$         write(85,*) results(1,isqrts),
!!$     $        results(3,isqrts)/results(2,isqrts),
!!$     $        results(4,isqrts)
!!$
!!$      enddo
!!$
!!$      stop
!!$
!!$      
!!$      end
!!$
!!$c================================================================= 

