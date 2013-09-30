!***************************************************************
!****m* /errorFunction
! NAME
! module errorFunction
! PURPOSE
! Includes error function
! SOURCE
!***************************************************************

module errorFunction
  
  implicit none
  
contains
  

  !*******************************************************
  !****f* errorFunction/errorFunc
  ! NAME
  ! real function errorFunc(z)
  ! PURPOSE
  ! Evaluates Error function
  ! errf(z)=2/SQRT(pi) Integral(e^(r^2) dr) on (0,z)
  ! INPUTS
  ! real z, 0<z<4
  ! NOTES
  ! Saves the error function in steps of dz with maximal value nz*dz
  !*******************************************************
  real function errorFunc(z)
    use constants, only : pi
    use output, only: Write_InitStatus 
    
    real, intent(in) :: z
    logical, save :: initFlag=.true.
    real, parameter :: dz=0.1
    integer, parameter :: nz=40
    integer :: index
    real, save, dimension(0:nz) :: field

    If(initFlag) then
       call init !Evaluates an array of possible results for the error function
       initFlag=.false.
    end if

    index=int(z/dz)
    If (index.gt.(nz-1)) then
       Write (*,*) 'z out of bounds in errorfunc :',z
       write (*,*) 'stop program'
       stop
    end if
    
    !Linear extrapolation in between data points
    errorFunc=field(index)+(field(index+1)-field(index))/dz*(z-float(index)*dz)
    
    contains
    
    subroutine init

      use output, only: Write_InitStatus 

      integer k,m
      double precision :: int, r
      integer,parameter :: steps=2000
      double precision :: dr

      call Write_InitStatus('array for error function',0)
      write(*,*) ' Maximal input for errf(z): z=',(nz-1)*dz

      field(0)=0.
      Do k=1,nz  
         !Evaluate Integral Int(e^(-r^2) dr) on the interval (0,z)
         !z=k*dz
         dr=(float(k)*dz)/float(steps)
         Int=0
         Do m=0,steps-1 
            r=float(m)*dr  !r=0,...,(nz-1)*dz
            Int=Int+Exp(-r**2)*dr
         End do
         field(k)=2/sqrt(pi)*Int
      end do

!!$      open(15,file='errorFunction.dat')
!!$      Do k=0,nz
!!$         Write(15,*) k*dz,field(k)
!!$      End do
!!$      close(15)

      call Write_InitStatus('array for error function',1)
    end subroutine init
  end function errorFunc

end module errorFunction
