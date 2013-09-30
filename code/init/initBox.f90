!*****************************************************************************
!****m* /initBox
! NAME
! module initBox
! PURPOSE
! Initializes nucleons for a calculation in a box of nuclear matter
!*****************************************************************************
module initBox

  use constants, only : rhoNull

  implicit none

  Private

  !*************************************************************************
  !****g* initBox/neutron_Density
  ! SOURCE
  !
  real, save :: neutron_Density= 0.084
  ! PURPOSE
  ! * neutron Density [fm^-3]
  !*************************************************************************

  !*************************************************************************
  !****g* initBox/proton_Density
  ! SOURCE
  !
  real, save :: proton_Density = 0.084
  ! PURPOSE
  ! * proton Density [fm^-3]
  !*************************************************************************

  !*************************************************************************
  !****g* initBox/energy_Density
  ! SOURCE
  !
  real, save :: energy_Density= 1.0
  ! PURPOSE
  ! * energy density [GeV fm^-3]
  !*************************************************************************

  !*************************************************************************
  !****g* initBox/fermiMotion
  ! SOURCE
  !
  logical, save :: fermiMotion=.true.
  ! PURPOSE
  ! * true=switch on fermi motion
  ! * false=switch off fermi motion
  !*************************************************************************


  Public :: initializeBox
  Public :: BoostToEps

contains

  !*************************************************************************
  !****s* initBox/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'initBox'.
  !*************************************************************************
  subroutine initInput

    use output, only: Write_ReadingInput
    use dichteDefinition
    use densityModule, only: get_densitySwitch,set_densityInput,densityAt
    use callstack, only: traceBack

    type(dichte) :: density

    !***********************************************************************
    !****n* initBox/initbox
    ! NAME 
    ! NAMELIST initBox
    ! PURPOSE
    ! Includes the input parameters:
    ! * proton_Density
    ! * neutron_Density
    ! * energy_Density
    ! * fermiMotion
    !***********************************************************************

    NAMELIST /initBox/ proton_Density, neutron_Density, energy_density, &
         & fermiMotion

    call Write_ReadingInput('initBox',0)
    rewind(5)
    read(5,nml=initBox)
    write(*,*) '  proton  Density   [1/fm^3] =',proton_Density
    write(*,*) '  neutron Density   [1/fm^3] =',neutron_Density
    write(*,*) '  Energy Density  [GeV/fm^3] =',energy_Density
    write(*,*) '  fermiMotion=', fermiMotion


    ! Taking care of the density routine
    select case(get_densitySwitch())
    case(1)
       ! ok
    case(3)
       write(*,'(A)') 'WARNING : Set parameters densityInput_XXX in '// &
            'densityModule to the values which are given here!!'
       call set_densityInput(proton_density,neutron_density)
    case default
       Write(*,*) 'Wrong densitySwitch for box-calculations:', get_densitySwitch()
       call TRACEBACK('Only 1 or 3 allowed!')
    end select

    density=densityAt((/0.,0.,0./))
    write(*,'(A,4F8.3)') '  proton  Density with densityAt=',density%proton
    write(*,'(A,4F8.3)') '  neutron Density with densityAt=',density%neutron

    call Write_ReadingInput('initBox',1)

  end subroutine initInput

  !*************************************************************************
  !****s* initBox/initializeBox
  ! NAME
  ! subroutine initializeBox(teilchen)
  ! PURPOSE
  ! Initialize nucleons in a box
  !*************************************************************************
  subroutine initializeBox(teilchen)
    use particleDefinition
    use IdTable, only: nucleon
    use random, only: rn
    use densityModule, only: gridsize
    use insertion, only: GarbageCollection
    use constants, only: mN
    use output, only: Write_InitStatus

    type(particle), dimension(:,:),intent(inOut) :: teilchen
    integer :: numberNeutrons, numberProtons
    integer :: i,k,index,offset, producedProtons



    ! Calculate the number of particles which are to be initialized

    call Write_InitStatus('box of nucleons',0)
    call initInput

    numberNeutrons=NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*neutron_density)
    numberProtons= NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*proton_density)
    write(*,*) ' Number of neutrons per ensemble=',numberNeutrons
    write(*,*) ' Number of protons  per ensemble=',numberProtons
    write(*,'(A,3F9.3)') '  Gridsize   =',gridsize(:)
    write(*,*) ' Size of box= (8*Gridsize) = ',8.*gridsize(1)*gridsize(2)*gridsize(3)
    write(*,*) ' Number Ensembles             =',size(teilchen(:,1))
    write(*,*) ' Number Particles per Ensemble=',size(teilchen(1,:))


    Do i=1,size(teilchen,dim=1)  !Loop over all ensembles
       producedProtons=0
       offset=0
       Do k=1,numberNeutrons+numberProtons   !Loop over all particles in the box
          Do !Search for empty space in particle vector
             index=k+offset
             If (teilchen(i,index)%ID > 0) then
                offset=offset+1
             else if (index.gt.size(teilchen,dim=2)) then
                Write (*,*) 'Real particle vector too small. Stop in initializeBox.'
                stop
             else
                exit
             end if
          end do
          call setToDefault(teilchen(i,index)) !set teilchen to its default values
          call setNumber(teilchen(i,index)) ! give each particle a unique number
          Teilchen(i,index)%event=1
          Teilchen(i,index)%ID=Nucleon
          Teilchen(i,index)%antiparticle=.false.
          Teilchen(i,index)%perturbative=.false.
          Teilchen(i,index)%productionTime=0.
          Teilchen(i,index)%mass=mN
          call chooseCharge
          Teilchen(i,index)%position(1)=(1.-2.*rn())*gridSize(1)
          Teilchen(i,index)%position(2)=(1.-2.*rn())*gridSize(2)
          Teilchen(i,index)%position(3)=(1.-2.*rn())*gridSize(3)
          call chooseMomentum
       End do
       If (producedProtons.ne.numberProtons) then
          Write(*,*) 'Problem in initPhaseSpace', producedProtons, numberProtons
       end if
    End do

    call GarbageCollection(teilchen)

    call Write_InitStatus('box of nucleons',1)

  contains

    !Choose randomly the charge of the nucleon
    subroutine chooseCharge
      real :: probabilityProton

      probabilityProton=float(numberProtons-producedProtons)/float(numberProtons+numberNeutrons-k+1)
      if(rn().le.probabilityProton) then
         Teilchen(i,index)%charge=1
         producedProtons=producedProtons+1
      else
         Teilchen(i,index)%charge=0
      end if
    end subroutine chooseCharge

    subroutine chooseMomentum
      use constants, only: mN, pi, hbarc

      real :: pFermi
      real, dimension(1:3) :: p
      real :: dens

      If(fermiMotion) then
         If (teilchen(i,index)%charge.eq.1) then
            dens=proton_density
         else
            dens=neutron_density
         End if
         pFermi=(3.*pi**2*dens)**(1./3.)*hbarc
         do   ! Monte Carlo distribution of momentum in sphere of radius pFermi
            p(1)=1.-2.*rn()
            p(2)=1.-2.*rn()
            p(3)=1.-2.*rn()
            if (sqrt(dot_product(p,p)).le.1.) exit
         end do
         p=p*pfermi

         Teilchen(i,index)%momentum(1:3)=p
      else
         Teilchen(i,index)%momentum(1:3)=0
      end if
      Teilchen(i,index)%momentum(0)=Sqrt(dot_product(p,p)+mN**2)
      !Assume vacuum dispersion relation:
      Teilchen(i,index)%velocity(1:3)=Teilchen(i,index)%momentum(1:3)/Teilchen(i,index)%momentum(0)
    end subroutine chooseMomentum


  end subroutine initializeBox


  !*************************************************************************
  !****s* initBox/BoostToEps
  ! NAME
  ! subroutine BoosToEps(teilchen)
  ! PURPOSE
  ! Boost the particles along z-axis such, that the energy density (=E/V) 
  ! corresponds to the input parameter
  !*************************************************************************
  subroutine BoostToEps(Teilchen)
    use particleDefinition
    use random, only : rn
    use densityModule, only : gridsize
    use constants, only: mN
    use lorentzTrafo, only: lorentz
    use densityModule, only : gridsize
    use output, only: Write_InitStatus

    type(particle), dimension(:,:),intent(inOut) :: teilchen
    integer :: i,j,nEns
    real :: betaCM, SumE
    real,dimension(1:3) :: betaV

    call Write_InitStatus('BoostToEps',0)

    betaCM = ( (proton_Density+neutron_Density)*mN / energy_density )**2
    if (betaCM >= 1.0) then
       write (*,*) 'betaCM invalid! stop!'
       stop
    end if
    betaCM = sqrt(1-betaCM)

    write(*,*) 'betaCM = ',abs(betaCM)
    write(*,*) ' ---> gamma = ', 1./sqrt(1-betaCM**2)

    SumE = 0.0
    nEns = size(teilchen,dim=1)
    Do i=1,nEns  !Loop over all ensembles
       Do j=1,size(teilchen,dim=2)  !Loop over all particles
          if (Teilchen(i,j)%ID <= 0) cycle
          if (rn() > 0.5) then
             betaV = (/0., 0.,  betaCM/)
          else
             betaV = (/0., 0., -betaCM/)
          end if
          call lorentz(betaV,Teilchen(i,j)%momentum, "BoostToEps")

          Teilchen(i,j)%event = -Teilchen(i,j)%number

          SumE = SumE + Teilchen(i,j)%momentum(0)

       end Do
    end Do

    write(*,*) 'E/V = ',SumE/(nEns* 8*gridsize(1)*gridsize(2)*gridsize(3))

    call Write_InitStatus('BoostToEps',1)

  end subroutine BoostToEps


end module initBox
