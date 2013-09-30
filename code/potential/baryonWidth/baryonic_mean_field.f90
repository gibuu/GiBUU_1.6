


module mean_field_common
  !****m* /mean_field_common
  ! NAME mean_field_common
  ! module 
  ! NOTES
  !  auxiliary module for module mean_field_baryon
  !
  !***
  save
  double precision :: hbc,m0      
  parameter( hbc=197.326d0,m0=938.21897d0)  

  double precision ::  rho, ueffn,ueffp, mnp,mnn, mn ,ueff
  double precision ::  fnn,fpn,fpp
  double precision ::  ueff_1,mn1,omega_f1
  double precision ::  omega_fp,omega_fn
  double precision ::  za, taun,taup,tau,kfp,kfn,tfac
  double precision :: iso,iso1


  integer ndim, l, li,kmin,kmax,lini, realflag
  parameter (ndim=5)  
  !******free parameters
  double precision omega_min,omega_max
  double precision domega,p_max,dp
  integer nstep_omega,nstep_p

  parameter (p_max=1.25d0)
  parameter (nstep_p=120)
  parameter (omega_min=-0.5d0,omega_max=0.5d0)
  parameter (nstep_omega=120)
  parameter (domega=(omega_max-omega_min)/nstep_omega)
  parameter (dp=p_max/nstep_p)

  double precision :: sigma_field(0:nstep_omega,0:nstep_p,2)
  double precision :: rsigma_field(0:nstep_omega,0:nstep_p,2)
  double precision :: rsigma_field1(0:nstep_omega,0:nstep_p,2)
  double precision :: strength,omega_f,p,omega,gamma_av,rsigma_av
  double precision :: min_sigma,omega_r,p_r     
!!$  double precision intas_field(0:nstep_p,0:nstep_p,0:nstep_omega)
end module mean_field_common

!***************************************************************************



module mean_field_skyrme1
  !****m* /mean_field_skyrme1
  ! NAME 
  ! module mean_field_skyrme1
  ! NOTES
  !  auxiliary module for module mean_field_baryon
  !
  !***
  double precision  t0,t1,t2,t3
  parameter(t0=-2490.23,t1=489.53,t2=-566.58,t3=13803.0)
  double precision  x0,x1,x2,x3
  parameter(x0=1.1318,x1=-0.8426,x2=-1.,x3=1.9219)
  double precision  ctbf
  parameter(ctbf=1./6.)
end module mean_field_skyrme1



module mean_field_baryon
  !****m* /
  ! NAME  mean_field_baryon
  ! module mean_field_baryon
  ! NOTES
  !  calculates the energy for the width of the baryon
  !
  !***

  PRIVATE
!!$  PUBLIC :: proton_omega_onshell_width
  PUBLIC :: omega_f_width
contains

  !***********************************************************************
!!$  double precision function intas(p2min,p2max,omega2)
!!$    use mean_field_common
!!$    implicit none      
!!$    double precision p2min,p2max,omega2
!!$    integer ip2min,ip2max,iomega
!!$
!!$    ip2min=min(nint(sqrt(p2min)/dp),nstep_p)
!!$    ip2max=min(nint(sqrt(p2max)/dp),nstep_p)
!!$    iomega=max(min(nint((omega2-omega_min)/domega),nstep_omega),0)
!!$    intas=intas_field(ip2min,ip2max,iomega)
!!$    return
!!$  end function intas

  double precision function ueff1(rho1,y1,tau1)   
    use mean_field_common
    use mean_field_skyrme1
    implicit none
    double precision rho1,y1,tau1
    ueff1=0.25*t0*(2*(2+x0)-2*(2*x0+1)*y1)*rho1        &
         +1./24.*t3*ctbf*rho1**(ctbf+1)*((2+x3)        & 
         -(2*x3+1)*(y1**2+(1-y1)**2))                  &  
         +1./12.*t3*rho1**(ctbf+1)*((2+x3)-(2*x3+1)*y1) & 
         +1./8.*(t1*(2+x1)+t2*(2+x2))*tau               &
         +1./8.*(t2*(2*x2+1)-t1*(2*x1+1))*tau1          


  end function ueff1
  !*************************************************************************


  subroutine mean_field_skyrme
    use mean_field_common     
    use mean_field_skyrme1
    use constants, only: pi
    implicit none
    !logical, save:: init=.true.
    double precision::yn
    tfac=hbc**2/(2*m0)

    !        za=0.4
    !        rho=0.16d0   

    kfp=(3*pi**2*za*rho)**(1./3.)
    kfn=(3*pi**2*(1-za)*rho)**(1./3.)     
    taup=3./5.*(3*pi**2)**(2./3.)*(za*rho)**(5./3.)
    taun=3./5.*(3*pi**2)**(2./3.)*((1-za)*rho)**(5./3.)
    tau=taup+taun
    ueffp=ueff1(rho,za,taup)/1000
    yn=1-za
    ueffn=ueff1(rho,yn,taun)/1000 
    mnp=mass(rho,za)/1000           
    mnn=mass(rho,yn)/1000
    omega_fp=(hbc/1000*kfp)**2/(2*mnp)+ueffp
    omega_fn=(hbc/1000*kfn)**2/(2*mnn)+ueffn

!20  format(4(e12.4,1x))
!!$ if (init) then
!!$      write(*,*)'rho=', rho,'[1/fm^3]'
!!$      write(*,*)'Z/A=', za
!!$      write(*,*)'kfp=',kfp, '[1/fm]'
!!$      write(*,*)'kfn=',kfn, '[1/fm]'
!!$      write(*,*)'mnn=',mnn, '[GeV]' 
!!$      write(*,*)'mnp=',mnp,'[GeV]' 
!!$      write(*,*)'ueffp=',ueffp,'[GeV]' 
!!$      write(*,*)'ueffn=',ueffn, '[GeV]' 
!!$      write(*,*)'omega_fp=',omega_fp, '[GeV]'
!!$      write(*,*)'omega_fn=',omega_fn, '[GeV]'
!!$  init=.false.
!!$ end if
  end subroutine mean_field_skyrme

  !*************************************************************************

  double precision function mass(rho1,y1)
    use mean_field_common
    use mean_field_skyrme1
    implicit none
    double precision rho1,y1
    mass=m0/(1.+0.25*m0/(hbc**2)*((t1*(2+x1)+t2*(2+x2))*rho1&
         +(t2*(2*x2+1)-t1*(2*x1+1))*y1*rho1))

  end function mass



  !***************************************************************************
!   double precision function eos(rho1)    
!     use mean_field_common
!     use mean_field_skyrme1
!     implicit none
!     double precision :: rho1
! !    double precision,external :: fm
!     eos=3d0/5d0*tfac*(3*pi**2/2.*rho1)**(2d0/3d0)*fm(5d0/3d0)       &
!          +1d0/8d0*t0*rho1*(2*(x0+2)-(2*x0+1)*fm(2d0))               &
!          +1d0/48d0*t3*rho1**(ctbf+1)*(2*(x3+2)-(2*x3+1)*fm(2d0))    &
!          +3d0/40d0*(3*pi**2/2d0)**(2d0/3d0)*rho1**(5d0/3d0)*((t1*(x1&
!          +2)+t2*(x2+2))*fm(5d0/3d0)                                 &
!          +1d0/2d0*(t2*(2*x2+1)-t1*(2*x1+1))*fm(8d0/3d0))
! 
!   end function eos


  !**************************************************************************
!   double precision function fm(m1)
!     use mean_field_common
!     use mean_field_skyrme1
!     implicit none
!     double precision :: m1   
!     fm=2**(m1-1)*(za**m1+(1-za)**m1)
!     return
!   end function fm
  !**************************************************************************

!!$  real function proton_omega_onshell_width(p,rho1,za1)
!!$
!!$    !****f* mean_field_baryon/proton_omega_onshell
!!$    ! NAME
!!$    ! function proton_omega_onshell_width(p,rho,za)
!!$    ! INPUTS
!!$    ! *  REAL, INTENT(IN) :: rho1   density of the particle rho_proton/neutron
!!$    ! *  REAL, INTENT(IN) :: p   
!!$    ! *  REAL, INTENT(IN) :: za
!!$    ! Result 
!!$    ! real
!!$    ! NOTES
!!$    ! Rho is the Density of the PARTICLE, rho_proton/neutron
!!$    !***
!!$
!!$    use mean_field_common , only : mnp, Ueffp,rho,za
!!$    use mean_field_skyrme1
!!$    implicit none
!!$    REAL, INTENT(IN) :: rho1   
!!$    REAL, INTENT(IN) :: p   
!!$    REAL, INTENT(IN) :: za1   
!!$    REAL :: za_old =0.
!!$    REAL :: rho_old=0.
!!$    if ((rho_old .NE. rho1 ) .OR. (za_old .NE. za1)) then
!!$       rho_old=rho1
!!$       rho=rho1
!!$       za_old=za1
!!$       za=za1
!!$       call mean_field_skyrme
!!$    end if
!!$
!!$    proton_omega_onshell_width= p**2/(2*mnp) + Ueffp 
!!$  end function proton_omega_onshell_width

  real function omega_f_width(rho1,za1)
    use mean_field_common  
    use mean_field_skyrme1
    implicit none
    REAL, INTENT(IN) :: rho1   
    REAL, INTENT(IN) :: za1   
    REAL, save  :: za_old =0.
    REAL, save  :: rho_old=0.
    if ((rho_old .NE. rho1 ) .OR. (za_old .NE. za1)) then
       rho_old=rho1
       rho=rho1
       za_old=za1
       za=za1
       call mean_field_skyrme
    end if
    omega_f_width=omega_fp
  end function omega_f_width


!!$  real function k_f(rho1,za1)
!!$    use mean_field_common  
!!$    use mean_field_skyrme1
!!$    implicit none
!!$    REAL, INTENT(IN) :: rho1   
!!$    REAL, INTENT(IN) :: za1   
!!$    REAL, save  :: za_old =0.
!!$    REAL, save  :: rho_old=0.
!!$    if ((rho_old .NE. rho1 ) .OR. (za_old .NE. za1)) then
!!$       rho_old=rho1
!!$       rho=rho1
!!$       za_old=za1
!!$       za=za1
!!$       call mean_field_skyrme
!!$    end if
!!$    k_f=(hbc/1000.)*kfp
!!$  end function k_f

end module mean_field_baryon

