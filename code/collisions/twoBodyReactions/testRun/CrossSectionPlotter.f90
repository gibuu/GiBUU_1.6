!***************************************************************************
!****m* /CrossSectionPlotter
! NAME
! program CrossSectionPlotter
!
! PURPOSE
! This is a standalone main routine. It is used to plot the cross sections
! used in the code, but is restriced to the total and elastic cross section (cf. GiBUU Homepage).
!
! INPUTS
! The Namelist "Plotter" in the Jobcard and additional (usual namelists).
!
! NOTES
! Particle 2 is the moving one, i.e. everything is in the rest frame of 
! particle 1.
!
!***************************************************************************

program CrossSectionPlotter
  use inputGeneral, only: readInputGeneral
  use version, only: printVersion
  use particleProperties, only: initParticleProperties, hadron, PartName
  use master_2Body, only: XsectionMaster, HiEnergyContrib
  use particleDefinition
  use mediumDefinition
  use preEventDefinition
  use energyCalc, only: energyDetermination
  use dichteDefinition
  use densitymodule, only: densityAt
  use RMF, only : getRMF_flag
  use twoBodyTools, only : sqrtS_free

  implicit none

  type(particle), dimension(1:2)  :: pair               ! incoming pair of particles
  type(preEvent), dimension(1:4)  :: finalState         ! produced final state
  real :: momentumLRF(0:3),betaToCF(1:3),rHiEnergy,srts,srtS_XS,srtS_vacuum,srtS_corr,mstar(2)
  real :: sigma_lo(2),sigma_hi(2),sigmaCEX,sigmaAnni,sigmaLbar,sigmaSbar,sigmaXiBar,sigmaJPsi
  integer :: i,q1,q2,id1,id2
  logical :: HiFlag,anti1,anti2
  type(medium) :: tMedium
  type(dichte) :: Dens
  real :: srts_max = 10.
!  real :: dens=0.168
!   real :: h

  call PrintVersion

  call readInputGeneral
  call initParticleProperties

  call readinput

!!$  tMedium%useMedium=.true.
!!$  tMedium%temperature=0.
!!$  tMedium%densityProton=dens/2.
!!$  tMedium%densityNeutron=dens/2.
!!$  tMedium%density=dens

  betaToCF = 0.

  pair%Id = (/id1,id2/)
  pair%charge = (/q1,q2/)
  pair%antiparticle = (/anti1,anti2/)
  pair%perturbative = .false.
  pair(1)%event = 1
  pair(2)%event = 2

  pair(1)%mass = hadron(id1)%mass
  pair(2)%mass = hadron(id2)%mass

  pair(1)%position = (/1.,0.,0./)
  pair(2)%position = (/0.,0.,0./)

  pair(1)%momentum(1:3) = 0.

  !pair(1)%momentum(0)=sqrt(pair(1)%mass**2+Dot_Product(pair(1)%momentum(1:3),pair(1)%momentum(1:3)))
  call energyDetermination (pair(1), betaToCF)
  call energyDetermination (pair(2), betaToCF)

  pair(1)%velocity = pair(1)%momentum(1:3)/pair(1)%momentum(0)


  Dens = densityAt(pair(1)%position)

  tMedium%useMedium     =.true.
  tMedium%temperature   =0.
  tMedium%density       =Dens%baryon(0)
  tMedium%densityProton =Dens%proton(0)
  tMedium%densityNeutron=Dens%neutron(0)


  write(*,*) '***********************'
  write(*,*) 'positions:'
  write(*,*) pair(1)%position
  write(*,*) pair(2)%position
  write(*,*) 'momenta:'
  write(*,*) pair(1)%momentum
  write(*,*) pair(2)%momentum
  write(*,*) 'sqrt(s)=',sqrts(pair(1),pair(2))
  write(*,*) 'Total momentum=',pair(1)%momentum+pair(2)%momentum
  write(*,*) '***********************'

  open(23,file='XS.dat')

  Do i=1,200000
     pair(2)%momentum(1:3)=(/0.,float(i)*0.01,0./)

!!$  Do i=1,1000
!!$     h = exp(i*0.001*(log(15.0)-log(0.1))+log(0.1))
!!$     pair(2)%momentum(1:3)=(/0.,h,0./)

     !pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))

     call energyDetermination(pair(2),betaToCF)
     pair(2)%velocity=pair(2)%momentum(1:3)/pair(2)%momentum(0)

!     srts=sqrts(pair(1),pair(2))

     srtS = sqrtS(pair,"generateFinalState, srtS")
     srtS_vacuum=sqrtS_free(pair)

     if (.not. getRMF_flag()) then
        srtS_XS = srtS_vacuum ! This is the srtS value the XS should be calculated with 
     else
        mstar(1) = sqrtS(pair(1),'generateFinalState, mstar(1)')
        mstar(2) = sqrtS(pair(2),'generateFinalState, mstar(2)')
        srtS_corr = srtS - mstar(1) - mstar(2) + pair(1)%mass + pair(2)%mass 
        srtS_XS = srtS_corr
     end if

     if (srts > srts_max) exit

     write(*,*) 'plab = ', absmom(pair(2)), ', srts = ', srts, srts_XS

     momentumLRF = pair(2)%momentum+pair(1)%momentum
     rHiEnergy = HiEnergyContrib(srts,pair%ID,pair%Antiparticle)
     sigma_lo = 0.
     sigma_hi = 0.

     if (rHiEnergy<1.0) &
        call XsectionMaster (srts_XS, pair, tMedium, momentumLRF, finalState, sigma_lo(1), sigma_lo(2), &
                             sigmaCEX, sigmaAnni, sigmaLbar, sigmaSbar, sigmaXiBar, sigmaJPsi, &
                             HiFlag, .true., ForceHiEnergy=.false.)

     if (rHiEnergy>0.0) &
        call XsectionMaster (srts_XS, pair, tMedium, momentumLRF, finalState, sigma_hi(1), sigma_hi(2), &
                             sigmaCEX, sigmaAnni, sigmaLbar, sigmaSbar, sigmaXiBar, sigmaJPsi, &
                             HiFlag, .true., ForceHiEnergy=.true.)


     write(23,'(10F12.4)') absmom(pair(2)), srts, (1.0-rHiEnergy)*sigma_lo+rHiEnergy*sigma_hi, sigma_lo, sigma_hi, &
                           rHiEnergy, kineticEnergy(pair(2))


  End do

  close(23)

contains

  subroutine readInput

    use output, only: Write_ReadingInput
    character(15) :: name(2)
    character(40) :: name1, name2

    NAMELIST /Plotter/ id1,id2,q1,q2,anti1,anti2,srts_max

    call Write_ReadingInput('Plotter',0)
    rewind(5)
    read(5,nml=Plotter)
    write(*,*) '  Id of first particle      : ',id1
    write(*,*) '  Id of second particle     : ',id2
    write(*,*) '  Charge of first particle  : ', q1
    write(*,*) '  Charge of second particle : ', q2
    write(*,*) '  antiparticle 1            : ', anti1
    write(*,*) '  antiparticle 2            : ', anti2
    write(*,*) '  srt(s)_max                : ', srts_max

    name(1) = PartName(id1,q1,anti1)
    name(2) = PartName(id2,q2,anti2)

    write(*,*)
    write(*,*) '  >> ',trim(name(2)),' --> ',trim(name(1)),' <<'
    write(*,*)

    call Write_ReadingInput('Plotter',1)

    !*****
    ! The following generates a html formated table of the collision scenario
    ! Better to use name(i) instead of this table ???
    !*****

    name1 = hadron(id1)%name
    name2 = hadron(id2)%name

    open(100,file="XS.html")
    write(100,*) '<table border="1" cellspacing="5" cellpadding="10">'
    write(100,*) '<tr>'
    write(100,*) '<th>  Scattering particles  </th> <th>',  name(1),'  </th> <th> ', name(2),'  </th> '
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<th>  </th> <th>',  name1,'  </th> <th> ', name2,'  </th> '
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<td>  Charge  </td> <td>',  q1,'  </td> <td> ', q2,'  </td> '
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<td>  Antiparticle  </td> <td>',  anti1,'  </td> <td> ', anti2,'  </td> '
    write(100,*) '</tr>'
    write(100,*) '</table>'
    close(100)

  end subroutine readInput


end program CrossSectionPlotter
