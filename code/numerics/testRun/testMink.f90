
program testMink
  call test
end program testMink

subroutine test()
  use minkowski
  use matrix_Module
  implicit none
  integer :: i,j, mu, nu
  complex,dimension(1:2,1:2) :: sig
  complex,dimension(0:3,0:3) :: gam
  real,dimension(0:3) :: q,p
  character, dimension(20) :: formatGamma='(4("(",2F5.1,")"))'

  call printMatrices()

  ! Test Slashing:
  write(*,*)
  write(*,*) 'Test 1:'
  q=(/1.23,1.9,2.3,7.1/)
  write(*,'(A,4F15.4)') 'q=',q
  write(*,*) 'tr [slashed(q) slashed(q)] - 4 q^mu q^_mu=0'
  write(*,'(A, 4F15.4)')' = ', trace(MatMul(slashed(q),slashed(q)))-4*SP(q,q)

  write(*,*)
  write(*,*) 'Test 2:'
  p=(/21.23,17.9,7.3,7.121/)
  write(*,'(A, 4F15.4)') 'p=',q
  write(*,'(A, 4F15.4)') 'q=',q
  write(*,*) 'tr [slashed(q) slashed(p)] - 4 q^mu p^_mu=0'
  write(*,'(A,4F15.4)')' = ', trace(MatMul(slashed(q),slashed(p)))-4*SP(q,p)

  ! Test daggering
  write(*,*)
  write(*,*) 'Test 3: gamma(i)^dagger-gamma_0 gamma(i) gamma_0 = 0'
  do i=0,3
     gam=dagger(gamma(:,:,i))-MatrixMult(gamma0,gamma(:,:,i),gamma0)
     write(*,*) 'i=', i
     do j=0,3
        write(*,'(4("(",2F5.1,")"))') gam(j,:)
       end do
     write(*,*)
  end do

  ! Test sigma^mu nu
  write(*,*)
  write(*,*) 'Test 3: sigma^mu nu-i/2*(2 gamma^mu gamma^nu-2g^mu nu) = 0'
  do mu=0,3
  do nu=0,3
     write(*,*) '#mu,nu=',mu,nu
     gam=sigma4(mu,nu)-ii*(MatrixMult(gamma(:,:,mu),gamma(:,:,nu))-metricTensor(mu,nu)*unit4)
     call printMatrixF(gam)
     write(*,*)

  end do
  end do


     write(*,*) 's_z(initial)=-1/2'
     write(*,*)
      write(*,*) '(sigma4(1,0)+ii*sigma4(2,0))*(/0.,1.,0.,0. /)'
     write(*,*) MatMul(MatMul(sigma4(1,0)+ii*sigma4(2,0),gamma5),(/0.,1.,0.,0. /))
     write(*,*)
      write(*,*) '(sigma4(1,3)+ii*sigma4(2,3))*(/0.,1.,0.,0. /)'
     write(*,*) MatMul(MatMul(sigma4(1,3)+ii*sigma4(2,3),gamma5),(/0.,1.,0.,0. /))
     write(*,*)
      write(*,*) '(gamma(1)+ii*gamma(2))*(/0.,1.,0.,0. /)'
     write(*,*) MatMul(MatMul(gamma1+ii*gamma2,gamma5),(/0.,1.,0.,0. /))

     write(*,*)
     write(*,*) 's_z(initial)=1/2'
     write(*,*)
     write(*,*) 'MatMul(MatMul(sigma4(0,3),gamma(5)),(/1.,0.,0.,0. /))'
     write(*,*) MatMul(MatMul(sigma4(0,3),gamma5),(/1.,0.,0.,0. /))
     write(*,*)
     write(*,*) 'MatMul(MatMul(gamma(0),gamma(5)),(/1.,0.,0.,0. /))'
     write(*,*) MatMul(MatMul(gamma0,gamma5),(/1.,0.,0.,0. /))
     write(*,*)
     write(*,*) 'MatMul(MatMul(gamma(3),gamma(5)),(/1.,0.,0.,0. /))'
     write(*,*) MatMul(MatMul(gamma3,gamma5),(/1.,0.,0.,0. /))
     write(*,*)
     write(*,*)
     write(*,*) 'gamma(1)+i gamma(2)'
!     call printMatrixF(gamma(1)+ii*gamma(2))

 !    write(*,*) '10+i 20'
  !   call printMatrixF(sigma4(1,0)+ii*sigma4(2,0))
  !   write(*,*) '13+i 23'
   !  call printMatrixF(sigma4(1,3)+ii*sigma4(2,3))

contains


  subroutine printMatrixF(A)
    implicit none
    integer :: j
    complex, intent(in),  dimension(0:3,0:3) :: a 
    
    do j=0,3
       write(*,'(4("(",2F5.1,")"))') A(j,:)
    end do
  end subroutine printMatrixF
end subroutine test
