program testRandom
  implicit none
  !call test_rnCos
  call test_rn_trueFalse
end program testRandom



subroutine test_rn_trueFalse
  use random
  implicit none
  integer,dimension(1:2) :: hits
  integer :: i

  do i=1,10000000
     if(rn_trueFalse()) then
        hits(1)=hits(1)+1
     else
        hits(2)=hits(2)+1
     end if
  end do
  write(*,*) hits
end subroutine test_rn_trueFalse

subroutine test_rnCos
use histf90
use random
implicit none
type(histogram) :: h,cosh
integer :: i
real :: x,y
integer :: n=1000000


call createHist(H, "x",-4.,4.,0.01)
call createHist(cosH, "cos(x)",-4.,4.,0.01)
do i=1,n
   x=rnCos()
   y=cos(x)
   call AddHist(h, x,1.)
   call AddHist(cosh, y,1.)
end do
call writeHist(h,11,0.,1./float(n))
call writeHist(cosh,12,0.,1./float(n))
end subroutine test_rnCos
