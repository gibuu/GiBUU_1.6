!***************************************************************************
!****m* /spline
! NAME
! module spline
! PURPOSE
! Includes routines for x-y -- cubic B-splines (i.e. the curve is 'smoother' than the data 
! points).
! 
!***************************************************************************
module spline

Private

public :: Bsplint2

contains


  !*********************************************************************
  !****if* spline/compareSize
  ! NAME
  ! logical function compareSize(x,y,z)
  ! PURPOSE
  ! check whether all three arrays have the same boundaries
  !*********************************************************************
!   logical function compareSize(x,y,z)
!     implicit none
!     real, intent(in) :: x(:),y(:),z(:)
!     compareSize=.false.
!     if (lbound(x,dim=1).ne.lbound(y,dim=1)) return
!     if (lbound(x,dim=1).ne.lbound(z,dim=1)) return
!     if (ubound(x,dim=1).ne.ubound(y,dim=1)) return
!     if (ubound(x,dim=1).ne.ubound(z,dim=1)) return
!     compareSize=.true.
!     return
!   end function compareSize


  !*********************************************************************
  !****s* spline/Bsplint2
  ! NAME
  ! subroutine Bsplint2(XA,YA,X,Y, Ind,Pref,SkipCalc)
  ! PURPOSE
  ! Returns an interpolated value for the data set (xa,ya)
  ! at position x.
  !
  ! Interpolation is done via cubic B-splines.
  ! INPUTS
  ! * real, intent(in) :: XA(:),YA(:) -- data set
  ! * real, intent(in) :: Y2A(:)      -- Derivative at data set points
  ! * real, intent(in) :: x           -- value where spline shall be evaluated
  ! * logical, intent(in) :: SkipCalc -- if true then use given Ind and Pref values [OPTIONAL] 
  ! RESULT
  ! * real,    intent(out):: y           -- y value of spline at x
  ! * integer, intent(inout):: Ind(-1:2)   -- indizes of values [OPTIONAL]
  ! * real,    intent(inout):: Pref(-1:2)  -- pre factors of values [OPTIONAL]
  ! NOTES 
  ! This "subroutine" should be a "function" !
  !
  ! No init as for "splint2" is necessary.
  !*********************************************************************
  subroutine Bsplint2(XA,YA,X,Y, Ind,Pref,SkipCalc)
    implicit none
    real, intent(in)  :: XA(:),YA(:)
    real, intent(in)  :: x
    real, intent(out) :: y
    integer, intent(inout), optional:: Ind(-1:2)
    real,    intent(inout), optional:: Pref(-1:2)
    logical, intent(in),    optional:: SkipCalc
    
    integer :: N, K, KLO, KHI
    real    :: H, t!, t1

    integer :: aInd(-1:2)
    real    :: aPre(-1:2)

    if (present(SkipCalc)) then
       if (SkipCalc) then
          if (.not.present(Ind) .or. .not.present(Pref)) then
             write(*,*) 'If SkipCalc, Ind and Pref has to be given. STOP!'
             stop
          endif
          aInd = Ind
          aPre = Pref
          goto 1000
       endif
    endif
    
    n=size(xa)
    KLO=1
    KHI=N
1   IF (KHI-KLO.GT.1) THEN
       K=(KHI+KLO)/2
       IF(XA(K).GT.X)THEN
          KHI=K
       ELSE
          KLO=K
       ENDIF
       GOTO 1
    ENDIF
    H=XA(KHI)-XA(KLO)
    t = (x-XA(KLO))/H
    
    aInd(-1) = klo-1
    aInd( 0) = klo
    aInd( 1) = klo+1
    aInd( 2) = klo+2

    
    if (klo .lt. 1) then
       write(*,*) "Error:",klo,lbound(xA)
       stop
    else if (klo .lt. 2) then
       aInd(-1) = aInd( 0)
!       write(*,*) 't=',t
       if(6.0*t .le. 1.0) then
!          write(*,*) 'linear (1)'

          aPre(-1) = 0
          aPre( 0) = ( 1.0-t )
          aPre( 1) = ( t )
          aPre( 2) = 0

          if (present(Ind)) Ind = aInd
          if (present(Pref)) Pref = aPre

          goto 1000

       else
!          write(*,*) 'solve (1)'
          call FindT(t, -3.,-3.,-1.+6.*t,  t)
       endif
    endif

    if (klo .gt. N-1) then
       write(*,*) "Error:",klo,ubound(xA)
       stop
    else if (klo .gt. N-2) then 
       aInd(2) = aInd(1)
       if(6.0*t .ge. 5.0) then
!          write(*,*) 'linear (2)'

          aPre(-1) = 0
          aPre( 0) = ( 1.0-t )
          aPre( 1) = ( t )
          aPre( 2) = 0

          if (present(Ind)) Ind = aInd
          if (present(Pref)) Pref = aPre

          goto 1000
       else
!          write(*,*) 'solve (2)'
          call FindT(t, 0.,-6.,6.*t, t)
       endif
    endif

    aPre(-1) =    -t**3 + 3.*t**2 - 3.*t + 1.
    aPre( 0) =  3.*t**3 - 6.*t**2        + 4.
    aPre( 1) = -3.*t**3 + 3.*t**2 + 3.*t + 1.
    aPre( 2) =     t**3

    aPre = aPre / 6.

    if (present(Ind)) Ind = aInd
    if (present(Pref)) Pref = aPre

1000 continue
    
    y = aPre(-1)*yA(aInd(-1))+aPre(0)*yA(aInd(0))+aPre(1)*yA(aInd(1))+aPre(2)*yA(aInd(2))



  contains
    !*********************************************************************
    !****is* Bsplint2/FindT
    ! NAME
    ! subroutine FindT(x,a2,a1,a0,t)
    ! PURPOSE
    ! Solve[t^3 + a2 t^2 + a1 t + a0 == 0, t] with the constraint
    ! 0 <= t <= 1
    !*********************************************************************

    subroutine FindT(x,a2,a1,a0,t)
    implicit none
    real, intent(in) :: x,a2,a1,a0
    real, intent(out):: t

    real :: r,q,D
    complex :: s1!,s2
    real :: z3!,z1,z2
    
!    write(*,*) 'T --> x:',x
    q = a1/3-a2**2/9
    r = (a1*a2-3*a0)/6-a2**3/27
    D = q**3+r**2

    if (D > 0) then
       write(*,*) 'D>0:',D,x
       stop
    endif

    s1 = (cmplx(r, sqrt(-D)))**(1./3.)
    !s2 = (cmplx(r,-sqrt(-D)))**(1./3.)

!    write(*,*) 's1,s2:',s1,s2

    !z1 = 2*REAL(s1)-a2/3
    !z2 = -REAL(s1)-sqrt(3.)*IMAG(s1)-a2/3
    z3 = -REAL(s1)+sqrt(3.)*AIMAG(s1)-a2/3

!    write(*,*) 'z1,z2,z3:',z1,z2,z3

    t = z3
!    write(*,*) 't=',t

  end subroutine FindT

  end subroutine Bsplint2

end module spline
