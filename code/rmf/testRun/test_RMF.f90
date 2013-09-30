program test_RMF

use particleProperties, only: initParticleProperties
use RMF
use constants, only : pi, hbarc

implicit none

real :: rmass, mass_seed, rhobar, shift, rhos, energyDensity, press, U, p_fermi
real :: ScalarPotential, VectorPotential
integer :: i

call initParticleProperties

open(1,file='NL3.dat',status='unknown')
write(1,*)'# rhobar, fm^-3:    p_fermi, GeV:   rhos, fm^-3:   m^*, GeV:'
write(1,*)'# (cont.)  E/A-m_N, GeV:  p, GeV^-3:  S, GeV:  V, GeV:   U, GeV:'

! NL1:
! rmass = 0.938

! NL3:
rmass = 0.939

! NL2-Lang:
!rmass = 0.938

! NLZ2:
!rmass = 0.9389


mass_seed = rmass

do i=1,200

   rhobar = i*0.001

   p_fermi = hbarc*(1.5*pi**2*rhobar)**0.333333

   call walecka(rhobar,shift,&
!               &em0=mass_seed,&
               &rhoscalar=rhos,&
               &endens=energyDensity,pressure=press,&
               &S=ScalarPotential,V=VectorPotential,potential=U)

   write(1,'(9(e13.6,1x))') rhobar, p_fermi, rhos, rmass-shift,&
                           &energyDensity/rhobar-rmass, press,&
                           &ScalarPotential, VectorPotential, U

   mass_seed = rmass - shift

end do

end program test_RMF
