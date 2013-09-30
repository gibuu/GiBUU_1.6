!***************************************************************************
!****m* /CallStack
! NAME
! module CallStack
!
! PURPOSE
! implement a wrapper around the ifort routine TRACEBACKQQ, so it can be used
! also when compiling with other compilers
!***************************************************************************
module CallStack

  implicit none
  private
  
  public :: traceback

contains

  !*************************************************************************
  !****s* CallStack/TRACEBACK
  ! NAME
  ! subroutine TRACEBACK(string,user_exit_code)
  !
  ! PURPOSE
  ! write out the call stack of the program in case ifort is used for compiling.
  !
  ! INPUTS
  ! * character*(*), optional :: string -- header line to write
  ! * integer, optional :: user_exit_code -- code whether return or stop 
  !   program
  ! * By specifying a user exit code of -1, control returns to the calling program. Specifying a user exit 
  !   code with a positive value requests that specified value be returned to the operating system. The default 
  !   value is 0, which causes the application to abort execution.
  ! NOTES
  ! This is a wrapper for the IFORT routine TRACEBACKQQ.
  ! See also the documentation there.
  !*************************************************************************
  subroutine TRACEBACK(string,user_exit_code)
#ifdef ifort
    use IFCORE
#endif

    character*(*), intent(in), optional:: string
    integer, intent(in), optional :: user_exit_code

#ifdef ifort
!#warning "compiling with ifort"

    if (present(string)) then
       if (present(user_exit_code)) then
          call TRACEBACKQQ(string=string,user_exit_code=user_exit_code)
       else
          call TRACEBACKQQ(string=string)
       endif
    else
       if (present(user_exit_code)) then
          call TRACEBACKQQ(user_exit_code=user_exit_code)
       else
          call TRACEBACKQQ()
       endif
    endif

#else
!#warning "not compiling with ifort"
    if (present(string)) then
       write(*,'(A)') string
    endif
# ifdef __GFORTRAN__
    ! ABORT is a GNU extension and gives a backtrace with gfortran 4.7 and above
    call abort()
# else

    write(*,*) '--- no call stack trace possible ---'

    if (present(user_exit_code)) then
       if (user_exit_code.eq.-1) return
       stop 123
    end if
    stop
# endif
#endif

  end subroutine TRACEBACK

end module CallStack
