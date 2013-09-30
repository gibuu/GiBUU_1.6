program testNuc

  use particleDefinition
  use potentialModule, only: potential_LRF
  use IdTable, only: nucleon
  use particleProperties, only: initParticleProperties, hadron

  implicit none
  type(particle) :: teilchen
  real :: dens,mom
  integer :: index_mom, id
  integer, parameter :: densPoints=10  !number of points in density
  integer, parameter :: momPoints=150   !number of points in momentum

  call initParticleProperties

  id=nucleon
  
  open(400,file='NucleonPotential.dat')
  write(400,*) '# 1: momentum'
  write(400,*) '# 2: pot at 0.16 fm^-3'
  write(400,*) '# 3: pot at 0.08 fm^-3'
  write(400,*) '# 4: pot at 0.04 fm^-3'
  write(400,*) '# 4: pot at 0.01 fm^-3'
  dens=0.16
  do index_mom=0,momPoints
     mom=index_mom*3.0/float(momPoints-1)
     teilchen%Id=id
     teilchen%mass=hadron(id)%mass
     teilchen%momentum(1:3)=(/mom,0.,0./)
     teilchen%momentum(0)=FreeEnergy(teilchen)
     write(400,'(5E20.8)') mom, potential_LRF(teilchen,0.16)  &
          &                   , potential_LRF(teilchen,0.08)  &
          &                   , potential_LRF(teilchen,0.04)  &
          &                   , potential_LRF(teilchen,0.01)
  end do
  close(400)


end program testNuc
