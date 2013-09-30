!***************************************************************************
!****m* /AnaEventDefinition
! NAME
! module AnaEventDefinition
! PURPOSE
! Here type(tAnaEvent) is defined. 
! Routines to work with this type are defined elsewhere, as e.g. in 
! module AnaEvent
!***************************************************************************
module AnaEventDefinition

  use particlePointerListDefinition, only: tParticleList

  implicit none

  !*************************************************************************
  !****t* AnaEventDefinition/tEvent
  ! NAME
  ! Type tEvent
  ! PURPOSE
  ! Tpye definition for events.
  !
  ! SOURCE
  !
  Type tAnaEvent
     type(tParticleList) :: particleList                 ! particles in the event
     integer,dimension(1:12,-2:2) :: numberParticles =0  ! Counters for stable particles (under strong decays)
  End Type tAnaEvent
  ! NOTES
  ! The field numberParticles includes the multiplicities of particles. 
  ! Antiparticles are not counted!
  !
  ! 1st Index (Particle ID):
  ! * 1=pion
  ! * 2=eta
  ! * 3=kaon
  ! * 4=kaonBar
  ! * 5=dMeson
  ! * 6=dBar
  ! * 7=Nucleon
  ! * 8=Lambda
  ! * 9=Sigma
  ! * 10=Xi
  ! * 11=Omega
  ! * 12=Any other
  ! 
  ! 2nd Index: 
  ! * Particle Charge (-2:2)
  !*************************************************************************

end module AnaEventDefinition
