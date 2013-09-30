!***************************************************************************
!****m* /nucleus
! NAME
! module nucleus
!
! PURPOSE
! Includes some routines handling the type "tNucleus" given in module 
! "nucleusDefinition".
! 
! INPUTS
! Namelists "target" and "projectile".
!***************************************************************************
module nucleus

  use nucleusDefinition, only: tNucleus

  implicit none
  PRIVATE

  PUBLIC :: getProjectile, getTarget

  type(tNucleus),save,pointer :: targetNuc     => NULL()     ! target nucleus
  type(tNucleus),save,pointer :: projectileNuc => NULL()     ! projectile nucleus

contains

  !*************************************************************************
  !****f* nucleus/getTarget
  ! NAME
  ! function getTarget()
  ! PURPOSE
  ! * Returns an initialized target nucleus resting at 0. with velocity=0.
  ! INPUTS
  ! * NONE
  ! OUTPUT
  ! * type(nucleus) :: getTarget 
  !*************************************************************************
  function getTarget()
    use nucleusDefinition
    type(tnucleus),pointer :: getTarget
    if (.not. associated(targetNuc)) then
      allocate(targetNuc)
      call initTarget
    end if
    getTarget => targetNuc
  end function getTarget


  !*************************************************************************
  !****f* nucleus/getProjectile
  ! NAME
  ! function getProjectile()
  ! PURPOSE
  ! Returns an initialized projectile nucleus resting at 0. with velocity=0.
  ! INPUTS
  ! * NONE
  ! OUTPUT
  ! * type(nucleus) :: getProjectile 
  !*************************************************************************
  function getProjectile()
    use nucleusDefinition
    type(tnucleus),pointer :: getProjectile
    if (.not. associated(projectileNuc)) then
      allocate(ProjectileNuc)
      call initProjectile
    end if
    getProjectile => ProjectileNuc
  end function getProjectile


  !*************************************************************************
  !****s* nucleus/initProjectile
  ! NAME
  ! subroutine initProjectile
  ! PURPOSE
  ! Initializes the projectile nucleus resting in the origin
  ! according to Information in namelist "projectile" in jobcard.
  !*************************************************************************
  subroutine initProjectile
    use output, only: Write_ReadingInput
    use nucleusDefinition, only: WriteNucleusStaticDens

    !***********************************************************************
    !****g* initProjectile/Projectile_A
    ! SOURCE
    integer,save :: Projectile_A=40
    ! PURPOSE
    ! mass of projectile nucleus
    !***********************************************************************

    !***********************************************************************
    !****g* initProjectile/Projectile_Z
    ! SOURCE
    integer,save :: Projectile_Z=20
    ! PURPOSE
    ! charge of projectile nucleus
    !***********************************************************************

    !***********************************************************************
    !****g* initProjectile/fermiMotion
    ! SOURCE
    logical,save :: fermiMotion=.true.
    ! PURPOSE
    ! switch, whether particles feel Fermi motion or not
    !***********************************************************************

    !***********************************************************************
    !****g* initProjectile/densitySwitch_static
    ! SOURCE
    !
    integer,save :: densitySwitch_static=3
    ! PURPOSE
    ! This switch is important, because it decides, which static density is 
    ! used to set up the testparticles in the nuclei before the first
    ! time-step. 
    ! 
    ! Possible values:
    ! * 0 : density=0.0
    ! * 1 : Static density uses Woods-Saxon according to H. Lenske
    ! * 2 : Static density according to NPA 554  
    ! * 3 : Static density according to Horst Lenske, 
    !   implements different radii for neutrons and protons
    ! * 4 : Static density according oscillator shell model
    ! * 5 : Density distribution is a sphere with density according to the
    !       input value of "fermiMomentum_input".
    ! * 6 : Static Density based on LDA, implemented by Birger Steinmueller
    ! * 7 : Static Density based on LDA + Welke potential
    ! * 8 : Static Density prescription according Relativistic Thomas-Fermi
    !       (Valid only in RMF-mode)
    ! 
    ! Possible nuclei for the different prescriptions:
    ! * 1 : A > 2 (only A > 16 makes sense)
    ! * 2 : 
    ! * 3 : 
    !   6->C(12), 8->O(16),O(18), 13->Al(27), 20->Ca(40),Ca(44), 79->Au(197)
    !   82->Pb(208)
    ! * 4: 2->He(4), 4->Be(9), 5->B(11), 6->C(12), 8->O(16)    
    !***********************************************************************

    !***********************************************************************
    !****g* initProjectile/fermiMomentum_input
    ! SOURCE
    !
    real,save :: fermiMomentum_input=0.225
    ! PURPOSE
    ! Input value of the fermi momentum for densitySwitch_static=5.
    !***********************************************************************

    !***********************************************************************
    !****n* initProjectile/projectile
    ! NAME
    ! NAMELIST /projectile/
    ! 
    ! PURPOSE
    ! Includes the input parameters for the projectile nucleus:
    ! * Projectile_A
    ! * Projectile_Z
    ! * fermiMotion
    ! * densitySwitch_static
    ! * fermiMomentum_input
    !***********************************************************************
    NAMELIST /projectile/ Projectile_A, Projectile_Z, fermiMotion, &
                          densitySwitch_static, fermiMomentum_input

    integer :: ios

    if (.not.projectileNuc%DoInit) return

    call Write_ReadingInput('projectile',0)
    rewind(5)
    read(5,nml=projectile,iostat=ios)
    call Write_ReadingInput('projectile',0,ios)
    write(*,*) ' A=',Projectile_A,' Z=',Projectile_Z
    write(*,*) ' fermiMotion=', fermiMotion
    write(*,*) ' densitySwitch_static = ',densitySwitch_static
    if (densitySwitch_static==5) then
       write(*,*) ' -> A naive Fermi gas distribution is initialized. No surface - just a sphere!'
       write(*,*) '    Fermi momentum = ',fermiMomentum_input
       projectileNuc%fermiMomentum_input=fermiMomentum_input
    end if

    projectileNuc%mass   = Projectile_A
    projectileNuc%charge = Projectile_Z
    projectileNuc%fermiMotion=fermiMotion
    projectileNuc%densitySwitch_static=densitySwitch_static

    call InitNucleus(projectileNuc,fermiMomentum_input)

    projectileNuc%DoInit = .false.

    !*************************************************************************
    !****o* nucleus/DensTab_projectile.dat
    ! NAME
    ! file DensTab_projectile.dat
    ! PURPOSE
    ! Density tabulation of the projectile nucleus at initialization.
    !*************************************************************************
    call WriteNucleusStaticDens('DensTab_projectile.dat',targetNuc)

    call Write_ReadingInput('projectile',2)

  end subroutine initProjectile


  !*************************************************************************
  !****s* nucleus/initTarget
  ! NAME
  ! subroutine initTarget
  ! PURPOSE
  ! Initializes the target nucleus resting in the origin
  ! according to Information in namelist "target" in jobcard.
  !*************************************************************************
  subroutine initTarget
    use output, only: Write_ReadingInput
    use nucleusDefinition, only: WriteNucleusStaticDens

    !***********************************************************************
    !****g* initTarget/Target_A
    ! SOURCE
    integer,save :: Target_A=40
    ! PURPOSE
    ! mass of target nucleus
    !***********************************************************************

    !***********************************************************************
    !****g* initTarget/Target_Z
    ! SOURCE
    integer,save :: Target_Z=20
    ! PURPOSE
    ! charge of target nucleus
    !***********************************************************************

    !***********************************************************************
    !****g* initTarget/fermiMotion
    ! SOURCE
    logical,save :: fermiMotion=.true.
    ! PURPOSE
    ! switch, whether particles feel Fermi motion or not
    !***********************************************************************

    !***********************************************************************
    !****g* initTarget/densitySwitch_static
    ! SOURCE
    !
    integer,save :: densitySwitch_static=3  
    ! PURPOSE
    ! This switch is important, because it decides, which static density is 
    ! used to set up the testparticles in the nuclei before the first
    ! time-step. 
    ! 
    ! Possible values:
    ! * 0 : density=0.0
    ! * 1 : Static density uses Woods-Saxon according to H. Lenske
    ! * 2 : Static density according to NPA 554  
    ! * 3 : Static density according to Horst Lenske, 
    !   implements different radii for neutrons and protons
    ! * 4 : Static density according oscillator shell model
    ! * 5 : Density distribution is a sphere with density according to the
    !       input value of "fermiMomentum_input".
    ! * 6 : Static Density based on LDA, implemented by Birger Steinmueller
    ! * 7 : Static Density based on LDA + Welke potential
    ! * 8 : Static Density prescription according Relativistic Thomas-Fermi
    !       (Valid only in RMF-mode)
    !
    ! Possible nuclei for the different prescriptions:
    ! * 1 : A > 2 (only A > 16 makes sense)
    ! * 2 : Be (9), C(12), O(16,18), Al(27), Ca(40), Ca(44), Fe(56), Cu(63), As(75), Ce(142), Sn(112, 116,120,124), 
    !       Ta(181), Au(197), Pb(208)    see densityStatic.f90   subroutine denspar for more info
    ! * 3 : 
    !   6->C(12), 8->O(16),O(18), 13->Al(27), 20->Ca(40),Ca(44), 79->Au(197)
    !   82->Pb(208)
    ! * 4: 2->He(4), 4->Be(9), 5->B(11), 6->C(12), 8->O(16)    
    !***********************************************************************

    !***********************************************************************
    !****g* initTarget/fermiMomentum_input
    ! SOURCE
    !
    real,save :: fermiMomentum_input=0.225
    ! PURPOSE
    ! Input value of the fermi momentum for densitySwitch_static=5.
    !***********************************************************************

    !***********************************************************************
    !****g* initTarget/ReAdjustForConstBinding
    ! SOURCE
    !
    logical, save :: ReAdjustForConstBinding = .false.
    ! PURPOSE
    ! If this flag is set to true, we use the selected density distribution
    ! only for a preliminary step, where we calculate the baryonic potential
    ! as function of r (which depends on the density distribution). 
    ! From the condition, that the binding energy has to be constant, we 
    ! deduce the distribution of the fermi momentum and thus the 'new' 
    ! density distribution. 
    !
    ! The tabulated density distribution is replaced via the 'new' one
    ! and all behaviour is as usual.
    !***********************************************************************

    !***********************************************************************
    !****g* initTarget/ConstBinding
    ! SOURCE
    !
    real, save :: ConstBinding = -0.008
    ! PURPOSE
    ! if 'ReAdjustForConstBinding' equals true, we a trying to readjust
    ! the fermi momentum and the density such, we quarantee this value
    ! for the binding energy. 
    !***********************************************************************

    !***********************************************************************
    !****n* initTarget/target
    ! NAME
    ! NAMELIST /target/
    ! 
    ! PURPOSE
    ! Includes the input parameters for the target nucleus:
    ! * Target_A    ---  mass of nucleus
    ! * Target_Z    ---  charge of nucleus
    ! * fermiMotion ---  Fermi motion yes/no
    ! * densitySwitch_static
    ! * fermiMomentum_input
    ! * ReAdjustForConstBinding
    ! * ConstBinding
    !***********************************************************************
    NAMELIST /target/ Target_A, Target_Z, fermiMotion, densitySwitch_static, &
                      fermiMomentum_input, ReAdjustForConstBinding, ConstBinding

    integer :: ios

    if (.not.targetNuc%DoInit) return
    
    call Write_ReadingInput('target',0)
    rewind(5)
    read(5,nml=target,iostat=ios)
    call Write_ReadingInput('target',0,ios)
    write(*,*) ' A=',Target_A,' Z=',Target_Z
    write(*,*) ' fermiMotion = ', fermiMotion
    write(*,*) ' densitySwitch_static = ',densitySwitch_static
    if (densitySwitch_static==5) then
       write(*,*) ' -> A naive Fermi gas distribution is initialized. No surface - just a sphere!'
       write(*,*) '    Fermi momentum = ',fermiMomentum_input
       targetNuc%fermiMomentum_input=fermiMomentum_input
    end if

    targetNuc%mass   = Target_A
    targetNuc%charge = Target_Z
    targetNuc%fermiMotion=fermiMotion
    targetNuc%densitySwitch_static=densitySwitch_static

    if (ReAdjustForConstBinding) then
       write(*,*)
       write(*,*) 'We us the initial density distribution only'
       write(*,*) 'for getting some educated guess of the baryonic'
       write(*,*) 'potential. With this we get a density distribution'
       write(*,*) 'ensuring a constant binding energy. !!!!!!!!!!!'
       write(*,*) 'E_binding = ',ConstBinding

       if (ConstBinding>0.) then
          write(*,*) 'Binding energy must be negative. stop'
          stop
       endif
    end if

    targetNuc%ReAdjustForConstBinding = ReAdjustForConstBinding
    targetNuc%ConstBinding = ConstBinding
    

    call InitNucleus(targetNuc,fermiMomentum_input)

    targetNuc%DoInit = .false.

    !*************************************************************************
    !****o* nucleus/DensTab_target.dat
    ! NAME
    ! file DensTab_target.dat
    ! PURPOSE
    ! Density tabulation of the target nucleus at initialization.
    !*************************************************************************
    call WriteNucleusStaticDens('DensTab_target.dat',targetNuc)
    call Write_ReadingInput('target',2)

  end subroutine initTarget


  !*************************************************************************
  !****s* nucleus/initNucleus
  ! NAME
  ! subroutine initNucleus(Nuc)
  ! PURPOSE
  ! Initializes a nucleus resting in the frame of calculation at r=0.
  ! INPUTS
  ! * type(tNucleus), pointer :: Nuc
  ! * real, intent(in) :: fermiMomentum_input
  ! NOTES
  ! In 'Nuc', mass and charge have to be set as input variables.
  !*************************************************************************
  subroutine initNucleus (Nuc, fermiMomentum_input)
    use NucD, only: nuclfit
    use densityStatic, only: staticDensityInit
    use nucleusDefinition, only: WriteNucleus
    use constants, only: pi, hbarc

    type(tNucleus), pointer :: Nuc
    real,intent(in) :: fermiMomentum_input
    real :: radius,surface,density 

    if (Nuc%charge > Nuc%mass) then
      write(*,*) 'Wrong nucleus in the input:'
      write(*,*) 'A, Z: ', Nuc%mass, Nuc%charge
      write(*,*) 'STOP'
      stop
    end if       

    select case(nuc%densitySwitch_static)
    case(5)
       If (Nuc%mass/=2*Nuc%charge) then
          write(*,*)
          write(*,*) 'WARNING: Using neutron-proton-symmetric Fermi-gas model'
          write(*,*) '         with non-symmetric nucleus! (in nucleus/initNucleus)'
          write(*,*) '         densitySwitch_static=', nuc%densitySwitch_static, ', A=',Nuc%mass,' , Z=',Nuc%charge
          write(*,*)
       end if
       Nuc%surface=0.
       Nuc%density=(fermiMomentum_input/hbarc)**3/(3.*pi**2)*2.
       ! A=4/3 pi r^3 *density
       Nuc%radius=(float(nuc%mass)/(4./3.*pi*nuc%density))**(1./3.)
       write(*,*) ' Fermi-gas parameters: '
       write(*,*) ' # Radius of nucleus =' ,nuc%radius
       write(*,*) ' # Density of nucleus=' ,nuc%density

    case DEFAULT

! the following is just a dummy stuff. DO IT BETTER !!!

       If (Nuc%mass>=12) then
          call nuclfit(Nuc%mass,Nuc%charge,radius,surface,density)
          Nuc%radius=radius
          Nuc%surface=surface
          Nuc%density=density
       else
          write(*,*)
          write(*,*) 'WARNING: radius parameters not initialized !!!!!!!'
          write(*,*)
          if (Nuc%mass == 2) then
             Nuc%radius=2.
          else
             Nuc%radius=1.2*Nuc%mass**(1/3)
          end if
          write(*,*) '...setting  Nuc%radius = ',Nuc%radius
       end if
    end select

    call staticDensityInit(nuc)
    call WriteNucleus(Nuc)

  end subroutine initNucleus


end module nucleus
