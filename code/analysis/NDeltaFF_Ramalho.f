!******************************************************************************
!****m* /NDeltaFF_Ramalho
! NAME
! module NDeltaFF_Ramalho
! FUNCTION
! Calculate the Delta -> N gamma* transition form factor,
! as given in arXiv:1205.2575 [hep-ph], but with high-W extensions.
! "Model 3" (private communication).
! Author: Gilberto Ramalho.
!******************************************************************************
      module NDeltaFF_Ramalho

      private

      public :: NDeltaTL

      contains

! ============================================================
! Block 1   
!
! Routines: Init (initial parameters)
! 
! ============================================================



      ! ======================================================
      ! Subroutine INIT 
      ! Read parameters and initialize calculations
      ! ======================================================
      SUBROUTINE Init
      implicit none
      real*8 :: mpi,mN,mD,mrho,drho,dmx
      real*8 :: cp,cm,dp,dm
      real*8 :: kp,km,lamb  
      real*8 :: N20,beta1,beta2,N0
      real*8 :: alf1,alf2,NS
      real*8 :: lpiD,LLpiD

      COMMON /Masses/mpi,mN,mD,mrho,drho,dmx
      COMMON /current1/cp,cm,dp,dm
      COMMON /current2/kp,km,lamb  
      COMMON /NucleonWF/N0,beta1,beta2
      COMMON /DeltaWF/alf1,alf2,NS
      COMMON /PionCloud/lpiD,LLpiD


      ! Read parameters -----------------------------------------------
!       OPEN(unit=1,file="input",STATUS="old")
!       ! Masses
!       read (1,*) mpi,mN,mD,mrho,drho,dmx
!       ! Quark current
!       read (1,*) cp,cm,dp,dm
!       read (1,*) kp,km,lamb    
!       ! Nucleon wave function
!       read (1,*) N20,beta1,beta2
!       ! Delta wave function
!       read (1,*) alf1,alf2
!       read (1,*) lpiD,LLpiD
!       CLOSE(1)

      mpi = 0.138
      mN = 0.938
      mD = 1.232
      mrho = 0.775
      drho = 0.1491
      dmx = 0.5964

      cp = 4.1600
      cm = 1.1150
      dp = -0.6860
      dm = -0.6860

      kp = 1.6410
      km = 1.8240
      lamb = 1.2100

      N20 = 11.2692
      beta1 = 0.0490
      beta2 = 0.7170

      alf1 = 0.33660
      alf2 = 0.33660

      lpiD = 0.4410
      LLpiD = 1.530

!       write (*,99) mpi,mN,mD,mrho,drho,dmx
!       write (*,99) cp,cm,dp,dm
!       write (*,99) kp,km,lamb  
!       write (*,99) N20,beta1,beta2
!       write (*,99) alf1,alf2
!       write (*,99) lpiD,LLpiD
! 
!  99   format(F9.4,3X,F9.4,3X,F9.4,3X,F9.4,3X,F9.4,3X,F9.4)

      ! Convert LLpiD to natural units
      LLpiD=LLpiD/mN**2.0

      N0=sqrt(N20)

      ! Calculate Delta normalization constant
      CALL Normaliza(NS)

      END SUBROUTINE




! ============================================================
! Block 2   
!
! Routines: NDeltaTL
!            
! 
! ============================================================



      ! **********************************************
      ! NDeltaTL
      ! calculate NDelta FF in timelike region q2t > 0
      ! 
      ! input: q2t=q^2
      !        W= Delta mass
      !        GM2= abs(G_M*)**2
      ! **********************************************
      SUBROUTINE NDeltaTL(q2t,W,GM2)
      implicit none
      real*8 :: q2t,W,Q2,GMre,GMim,GM2
      real*8 :: alf1,alf2,NS,mpi,mN,mD0,mrho,drho
      real*8 :: GMa,GMb,GMpiA,GMpiB
      logical :: first = .true.

      COMMON /DeltaWF/alf1,alf2,NS
      COMMON /Masses/mpi,mN,mD0,mrho,drho

      if (first) then
        call Init
        first = .false.
      end if

      Q2=-q2t

      if ((q2t < 0d0) .or. (q2t > (W-mN)**2.0)) then 
         GM2 = 0.
         return
      end if    

      CALL DeltaTL3(Q2,NS,W,GMa,GMb)

      if (abs(Q2) > 1d-6) then
         CALL PionCloudR3(Q2,GMpiA,GMpiB)
      else 
         CALL PionCloudR3S(0d0,GMpiA,GMpiB)
      end if

      GMre=GMa+GMpiA
      GMim=GMb+GMpiB
      GM2=GMre**2.0 + GMim**2.0
     
      END SUBROUTINE




! ============================================================
! Block 3   
!
! Routines: DeltaFF0,PionCloudR (spacelike FF)
!           DeltaTL3,PionCloudR3,PionCloudR3S 
!           Normaliza,GAUSS 
! 
! ============================================================



      ! Time-like region
      ! gR=gR(Q2)

      SUBROUTINE DeltaTL3(Q2,NS,W,GMa,GMb)
      implicit none
      integer, parameter :: nmax=400
      real*8 :: mN,W
      real*8 :: q2vec,qvec,EN,omega
      integer :: np,Nk,Nz,ik,iz
      real*8 :: mpi,mD,mD0,mrho,drho,dmx
      real*8 :: Q2,N0,GMS,GM,GMa,GMb
      real*8 :: cp,cm,dp,dm
      real*8 :: kp,km,lamb    
      real*8 :: N2,beta1,beta2
      real*8 :: NS,alf1,alf2
      real*8 :: PI,mdel,mdel2,dm2,sm2,sm,dmp
      real*8 :: intfac,fac,GD3,F1m,F2m,jm
      real*8 :: Xk(nmax),WX(nmax),U(nmax),WU(nmax)
      real*8 :: cut,Q2m,qD,qD2
      real*8 :: k,dk,Es,z,dz,k2,z2,intS
      real*8 :: chiN,chiD,psiN,psiD
      real*8 :: chi1,psi1,chi2,psi2
      real*8 :: lpiD,LLpiD,mrho2,mR,Rfactor,gR
      real*8 :: Rfactor2A,Rfactor2B,mh2,gmh
      real*8::  Rfactor3A,Rfactor3B,Rfactor2
      real*8 :: RfactorA,RfactorB,F1mA,F1mB,F2mA,F2mB
      real*8 :: jmA,jmB,Q2R,mpi2


      COMMON /Masses/mpi,mN,mD0,mrho,drho,dmx
      COMMON /current1/cp,cm,dp,dm
      COMMON /current2/kp,km,lamb  
      COMMON /NucleonWF/N0,beta1,beta2
      COMMON /DeltaWF/alf1,alf2

      gmh=dmx
      gmh=4d0*(-Q2/(4d0*mN**2.0-Q2))**2.0*gmh

      ! Q2R threshold for gR(Q2)

      Q2R=.076176d0   ! = 4*mpi**2
      
      gR=0d0

      if ( Q2 < - Q2R) then
         gR=((-Q2R-Q2)/(-Q2R+mrho**2.0))**1.5*mrho/sqrt(-Q2)*drho
      end if   
      
      mD=W

      q2vec=(W**2.0+ mN**2.0+ Q2)**2.0/(4d0*W**2.0)-mN**2.0
      qvec=sqrt(q2vec)

      omega=(W**2.0-mN**2.0-Q2)/(2d0*W)

      En=sqrt(mN**2.0+q2vec)


      ! Constants  ----------------------------------------------------
      PI=4d0*ATAN(1d0)
      mdel=mD/mN        ! Delta mass (nucleon unities)
      mdel2=mdel*mdel

      dm2=mdel*mdel-1d0   ! md**2-mn**2 
      sm2=mdel*mdel+1d0   ! md**2+mn**2
      
      dmp=mdel-1d0
      sm=mdel+1d0

      ! Integral factor
      intfac=.5d0/(2d0*PI)**2.0
      GD3=3d0/(1d0+Q2/0.71d0)**2.0

      ! Integral coeficient
      fac=1d0/sqrt(3d0)

      ! Integration parameters ----------------------------------------

      Nk=200     ! k grid
      Nz=300     ! z grid

      CALL GAUSS(0d0,1d0,Xk,WX,Nk)
      CALL GAUSS(-1d0,1d0,U,WU,Nz)

      cut=0.5d0  

      Q2m=Q2/mN**2     ! Normalized momenta

      ! current functions

      mrho2=(mrho/mN)**2.0

      Rfactor=mrho**2/(mrho**2+Q2)

      RfactorA=mrho**2.0*(mrho**2+Q2)/
     &                 ( (mrho**2+Q2)**2.0 + mrho**2*gR**2)

      RfactorB=mrho**2.0*mrho*gR/
     &                 ( (mrho**2+Q2)**2.0 + mrho**2*gR**2)


      mh2=4d0*mN**2.0

      Rfactor2=mh2/((mh2+Q2)**2.0+mh2*gmh**2.0)


      Rfactor3A=(mh2+Q2)*Rfactor2
      Rfactor3B=sqrt(mh2)*gmh*Rfactor2
      
      Rfactor2A=Rfactor2**2.0*
     &          ( (mh2+Q2)**2.0 - Mh2*gmh**2.0)
      RFactor2B=Rfactor2**2.0*
     &          ( 2d0*(mh2+Q2)*sqrt(mh2)*gmh )


      ! F1m    
      F1m=lamb + (1d0-lamb)*Rfactor+
     &           cm*.25d0*Q2m/(1d0+.25d0*Q2m)**2.0



      ! F2m
      F2m=km*(   dm*Rfactor+
     &         (1d0-dm)/(1d0+.25d0*Q2m))
     
 

      ! Real part

      ! F1m    
      F1mA=lamb + (1d0-lamb)*RfactorA+
     &           cm*Rfactor2A*Q2/mh2

      ! F2m
      F2mA=km*(   dm*RfactorA+
     &         (1d0-dm)*Rfactor3A  ) 

      ! Imaginary part

      ! F1m    
      F1mB=lamb*0 + (1d0-lamb)*RfactorB+
     &           cm*Rfactor2B*Q2/mh2


      ! F2m
      F2mB=km*(   dm*RfactorB+
     &         (1d0-dm)*Rfactor3B )

      ! Initiate integration

      intS=0d0

      do ik=1,Nk
         k=cut*Xk(ik)/(1d0-Xk(ik))
         k2=k*k
         dk=cut/(1d0-Xk(ik))**2.0*WX(ik)
         Es=sqrt(1d0+k2)
         do iz=1,Nz
            z=U(iz)
            z2=z*z
            dz=WU(iz)

            ! Wave functions (Delta rest frame)
            chi1=2d0*Es
            psi1=1d0/((alf1-2d0+chi1)*(alf2-2d0+chi1)**2.0)

            chi2=2d0*(En*Es-qvec*k*z)/mN
            psi2=1d0/((beta1-2d0+chi2)*(beta2-2d0+chi2))

            intS=intS+psi1*psi2*k2*dk*dz/Es

         end do   ! End of iz cicle  
      end do      ! End of ik cicle

      intS=  N0*NS*intfac*intS

      jm=F1m+sm*F2m/2d0     ! F1m (normalized to 1)

      jmA=F1mA+sm*F2mA/2d0   
      jmB=F1mB+sm*F2mB/2d0 

      ! S-state Form Factors
      GMS=8d0*jm/(3d0*sm)*intS*fac   
  
      GMa=8d0*jmA/(3d0*sm)*intS*fac
      GMb=8d0*jmB/(3d0*sm)*intS*fac
 
      end subroutine





      SUBROUTINE PionCloudR3(Q2,GMpiA,GMpiB)
      implicit none
      real*8 :: mpi,mN,mD,mrho,drho,dmx
      real*8 :: Q2,GMpi,GD3,Q2m,GD3p,GMpiA,GMpiB
      real*8 :: lpiD,LLpiD
      real*8 :: gR,Rfactor,cut2,RfactorA,RfactorB,fpi
      real*8 :: pi,t1,t2,mpi2,mrho2
      real*8 :: a,b,gm,RfactorC

      COMMON /Masses/mpi,mN,mD,mrho,drho,dmx
      COMMON /PionCloud/lpiD,LLpiD

      pi=4d0*atan(1d0)
      
      Q2m=Q2/mN**2.0

      mpi2=mpi**2.0
      mrho2=mrho**2.0

      gR=drho
      gm=dmx

      gm=gm/mN

      gm=4d0*gm*(-Q2m/(LLpiD-Q2m))**2.0


      t1=Q2*log(-Q2/mpi2)*(gR/mpi)/pi
      t2=Q2*(gR/mpi)

      cut2=.71d0

      Rfactor=mrho2/( (mrho2+Q2 + t1)**2.0 + t2**2.0)

      RfactorA=Rfactor*(mrho2+Q2+t1)
      RfactorB=Rfactor*t2
      
      RfactorC=( LLpiD/((LLpiD+Q2m)**2+ LLpiD*gm**2) )**2

      fpi=lpiD*RfactorC

      a= (LLpiD+Q2m)**2.0 - LLpiD*gm**2.0
      b= 2d0*(LLpiD+Q2m)*sqrt(LLpiD)*gm

      GMpiA=3d0*fpi*(a*RfactorA - b*RfactorB)
      GMpiB=3d0*fpi*(b*RfactorA + a*RfactorB)

      END SUBROUTINE



      SUBROUTINE PionCloudR3S(Q2,GMpiA,GMpiB)
      implicit none
      real*8 :: mpi,mN,mD,mrho,drho,dmx
      real*8 :: Q2,GMpi,GD3,Q2m,GD3p,GMpiA,GMpiB
      real*8 :: lpiD,LLpiD
      real*8 :: gR,Rfactor,cut2,RfactorA,RfactorB,fpi
      real*8 :: pi,t1,t2,mpi2,mrho2
      real*8 :: a,b,gm,RfactorC

      COMMON /Masses/mpi,mN,mD,mrho,drho,dmx
      COMMON /PionCloud/lpiD,LLpiD

      pi=4d0*atan(1d0)
      
      mpi2=mpi**2.0
      mrho2=mrho**2.0

      gR=drho
      gm=dmx

      gm=gm/mN*0d0

      t1=Q2*log(-Q2/mpi2)*(gR/mpi)/pi
      t2=Q2*(gR/mpi)

      t1=0d0
      t2=0d0

      cut2=.71d0

      Rfactor=mrho2/( (mrho2+Q2 + t1)**2.0 + t2**2.0)

      RfactorA=Rfactor*(mrho2+Q2+t1)
      RfactorB=Rfactor*t2

      Q2m=Q2/mN**2.0

      RfactorC=( LLpiD/((LLpiD+Q2m)**2+ LLpiD*gm**2) )**2

      fpi=lpiD*RfactorC

      a= (LLpiD+Q2m)**2.0 - LLpiD*gm**2.0
      b= 2d0*(LLpiD+Q2m)*sqrt(LLpiD)*gm

      GMpiA=3d0*fpi*(a*RfactorA - b*RfactorB)
      GMpiB=3d0*fpi*(b*RfactorA + a*RfactorB)

      END SUBROUTINE



      !**************************************************************
      ! Routine DeltaFF0
      ! Input: Q2,NS,mD 
      !        
      ! Output: GM
      !
      ! Uses Delta rest frame coordenates
      ! SS transition (frame invariant): intS
      !**************************************************************
!      SUBROUTINE DeltaFF0(Q2,NS,mD,GM)
!      implicit none
!      integer, parameter :: nmax=400
!      integer :: np,Nk,Nz,ik,iz
!      real*8 :: mpi,mN,mD,mrho,mD0
!      real*8 :: Q2,N0,GMS,GM,GE,GC,fpi
!      real*8 :: cp,cm,dp,dm
!      real*8 :: kp,km,lamb    
!      real*8 :: N2,beta1,beta2
!      real*8 :: NS,alf1,alf2
!      real*8 :: PI,mdel,mdel2,dm2,sm2,sm,dmp
!      real*8 :: intfac,fac,GD3,F1m,F2m,jm
!      real*8 :: Xk(nmax),WX(nmax),U(nmax),WU(nmax)
!      real*8 :: cut,Q2m,qD,qD2,En
!      real*8 :: k,dk,Es,z,dz,k2,z2,intS
!      real*8 :: chiN,chiD,psiN,psiD
!      real*8 :: lpiD,LLpiD,mrho2
!
!      COMMON /Masses/mpi,mN,mD0,mrho
!      COMMON /current1/cp,cm,dp,dm
!      COMMON /current2/kp,km,lamb  
!      COMMON /NucleonWF/N0,beta1,beta2
!      COMMON /DeltaWF/alf1,alf2
!      COMMON /PionCloud/lpiD,LLpiD
!
!      
!
!      ! Constants  ----------------------------------------------------
!      PI=4d0*ATAN(1d0)
!      mdel=mD/mN        ! Delta mass (nucleon unities)
!      mdel2=mdel*mdel
!
!      dm2=mdel*mdel-1d0   ! md**2-mn**2 
!      sm2=mdel*mdel+1d0   ! md**2+mn**2
!      
!      dmp=mdel-1d0
!      sm=mdel+1d0
!
!      ! Integral factor
!      intfac=.5d0/(2d0*PI)**2.0
!      GD3=3d0/(1d0+Q2/0.71d0)**2.0
!
!      ! Integral coeficient
!      fac=1d0/sqrt(3d0)
!
!      ! Integration parameters ----------------------------------------
!
!      Nk=200     ! k grid
!      Nz=300     ! z grid
!
!      CALL GAUSS(0d0,1d0,Xk,WX,Nk)
!      CALL GAUSS(-1d0,1d0,U,WU,Nz)
!
!      cut=0.5d0  
!
!      Q2m=Q2/mN**2     ! Normalized momenta
!
!
!      ! Delta rest frame variables
!      qD2=.25d0*(sm*sm+Q2m)*(dmp*dmp+Q2m)/mdel2
!      qD=sqrt(qD2)
!
!      En=(mdel2+1d0+Q2m)/(2d0*mdel)
!
!      ! current functions
!
!      mrho2=(mrho/mN)**2.0
!      
!      ! F1m     (No pion cloud)
!      F1m=lamb + (1d0-lamb)/(1d0+Q2m/mrho2)+
!     &           cm*.25d0*Q2m/(1d0+.25d0*Q2m)**2.0
!
!      ! F2m
!      F2m=km*(   dm/(1d0+Q2m/mrho2)+
!     &         (1d0-dm)/(1d0+.25d0*Q2m))
!
!      ! Initiate integration
!
!      intS=0d0
!
!      do ik=1,Nk
!         k=cut*Xk(ik)/(1d0-Xk(ik))
!         k2=k*k
!         dk=cut/(1d0-Xk(ik))**2.0*WX(ik)
!         Es=sqrt(1d0+k2)
!         do iz=1,Nz
!            z=U(iz)
!            z2=z*z
!            dz=WU(iz)
!
!            ! Wave function parameters
!            ! Delta Rest Frame
!            chiN=2d0*En*Es+2d0*k*qD*z
!            chiD=2d0*Es
!              
!            ! Nucleon wave function
!            psiN=1d0/((beta1-2d0+chiN)*(beta2-2d0+chiN))
!
!            ! Delta wave functions
!            ! S-state
! 
!            psiD=1d0/((alf1-2d0+chiD)*(alf2-2d0+chiD)**2.0)
!
!            intS=intS+psiN*psiD*k2*dk*dz/Es
!
!
!         end do   ! End of iz cicle  
!      end do      ! End of ik cicle
!
!      intS=  N0*NS*intfac*intS
!
!      jm=F1m+sm*F2m/2d0     ! F1m (normalized to 1)
!
!      ! S-state Form Factors; no pion cloud
!      GMS=8d0*jm/(3d0*sm)*intS*fac   
!  
!      GM=GMS
!
!      END SUBROUTINE





!      SUBROUTINE PionCloudR(Q2,GMpi)
!      implicit none
!      real*8 :: mpi,mN,mD,mrho,drho
!      real*8 :: pi,Q2,GMpi,GD3,Q2m,GD3p
!      real*8 :: lpiD,LLpiD
!
!      COMMON /Masses/mpi,mN,mD,mrho,drho
!      COMMON /PionCloud/lpiD,LLpiD
!
!      pi=4d0*atan(1d0)
!
!      GD3=3d0/(1d0+Q2/0.71d0)**2.0
!
!      GD3p=3d0*mrho**2.0/
!     &  (mrho**2.0 + Q2 + (drho/mpi)*Q2/pi*log(Q2/mpi**2.0))
!
!      Q2m=Q2/mN**2.0
!
!      GMpi=lpiD*(LLpiD/(LLpiD+Q2m))**2*GD3p
!
!      END SUBROUTINE





      !**************************************************************
      ! Routine Normalize
      ! Evaluate normalization constants NS
      ! INPUT: alf1,alf2 (COMMON)
      !
      ! OUTPUT: NS
      !**************************************************************
      SUBROUTINE Normaliza(NS)
      implicit none
      integer, parameter :: nmax=400
      integer :: Nk,Nz,ik,iz
      real*8 :: mpi,mN,mD,mrho
      real*8 :: NS,alf1,alf2
      real*8 :: PI,mdel,mdel2,dm2,sm2,sm,dm
      real*8 :: intfac
      real*8 :: Xk(nmax),WX(nmax),U(nmax),WU(nmax)
      real*8 :: cut,qD,qD2,En
      real*8 :: k,dk,Es,z,dz,k2,z2
      real*8 :: intS,chiD,psiD

      COMMON /Masses/mpi,mN,mD,mrho 
      COMMON /DeltaWF/alf1,alf2



      ! Constants  ----------------------------------------------------
      PI=4d0*ATAN(1d0)
      mdel=mD/mN        ! Delta mass (nucleon unities)
      mdel2=mdel*mdel

      dm2=mdel2-1d0   ! md**2-mn**2 
      sm2=mdel2+1d0   ! md**2+mn**2
      
      dm=mdel-1d0
      sm=mdel+1d0

      ! Integral factor
      intfac=.5d0/(2d0*PI)**2.0
      !GD3=3d0/(1d0+Q2/0.71d0)**2.0

      ! Integration parameters ----------------------------------------

      Nk=200     ! k grid
      Nz=300     ! z grid

      CALL GAUSS(0d0,1d0,Xk,WX,Nk)
      CALL GAUSS(-1d0,1d0,U,WU,Nz)

      cut=0.5d0  

      ! Delta rest frame variables Q2=0
      qD2=dm2/(2d0*mdel)
      qD=sqrt(qD2)
      En=sm2/(2d0*mdel)
      

      ! Initiate integration

      intS=0d0

      do ik=1,Nk
         k=cut*Xk(ik)/(1d0-Xk(ik))
         k2=k*k
         dk=cut/(1d0-Xk(ik))**2.0*WX(ik)
         Es=sqrt(1d0+k2)
         do iz=1,Nz
            z=U(iz)
            z2=z*z
            dz=WU(iz)

            ! Wave function parameters
            ! Delta Rest Frame
            chiD=2d0*Es   
              
            ! Delta wave function initial state
            psiD=1d0/((alf1-2d0+chiD)*(alf2-2d0+chiD)**2.0)

            intS=intS+(psiD*psiD)*k2*dk*dz/Es


         end do   ! End of iz cicle  
      end do      ! End of ik cicle

      intS= intfac*intS

      NS=1d0/sqrt(intS)
      
      END SUBROUTINE






!********************************************************************
!     Gaussian method: calculate coordinates and weights
!********************************************************************
      SUBROUTINE GAUSS(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),W(N)
      PARAMETER (EPS=3.D-12)
      M=(N+1)/2  
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(3.14159265358979D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0   
          DO 11 J=1,N   
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END SUBROUTINE



      end module
