!****m* /NNPION_NN
! NAME
! module NNPION_NN
! PURPOSE
! Implements the Matrix elements for N N Pion -> N N
!***
module NNPION_NN

  implicit none
  Private

  Public :: matrix_NN_NNPion,gamma_NNPion_NN,gamma_NNPion_NN_output

contains

  function gamma_NNPion_NN(srts,Epion,Enucleon,rhoProton,rhoNeutron,chargePion,chargeNucleon,OutPutFlag) Result(gamma)
    !****s* NNPION_NN/gamma_NNPion
    ! NAME
    ! function gamma_NNPion_NN(srts,Epion,Enucleon,rhoProton,rhoNeutron,chargePion,chargeNucleon,OutPutFlag) Result(gamma)
    ! PURPOSE
    ! Evaluates Gammas for pion absorption on a pair of nucleons. There exist the following possible channels:
    ! * 1: p p pi^0-> p p  ; which is equivalent to n n -> n n pi^0 
    ! * 2: n n pi^+-> n p  ; which is equivalent to  p n -> p p pi^-
    ! * 3: p n pi^0-> p n  ; 
    ! * 4: p n pi^+-> p p  ; which is equivalent to n n -> p n pi^-
    ! OUTPUT
    ! * real :: gamma ! absorption rate of pion on one specific pair on nucleons of the input
    ! INPUTS
    ! * real, intent(in) :: srts ! SQRT(s)
    ! * real, intent(in) :: Epion ! Energy of pion
    ! * real, intent(in),dimension(1:2) :: Enucleon !  Energies of nucleons
    ! * real, intent(in) :: rhoProton 
    ! * real, intent(in) :: rhoNeutron
    ! * integer, intent(in) :: chargePion
    ! * integer, dimension(1:2), intent(in) :: chargeNucleon
    ! * real,intent(out) :: gamma
    ! * logical, optional,intent(in) :: outputFlag
    !***

    use constants, only : pi, mN, mPi

    real, intent(in) :: srts ! SQRT(s)
    real, intent(in) :: Epion ! Energy of pion
    real, intent(in),dimension(1:2) :: Enucleon !  Energies of nucleons
    real, intent(in) :: rhoProton 
    real, intent(in) :: rhoNeutron
    integer, intent(in) :: chargePion
    integer, dimension(1:2), intent(in) :: chargeNucleon
    logical, optional,intent(in) :: outputFlag

    real :: gamma ! OUTPUT

    real, dimension(1:4) ::matrixelements,rhoTwo
    real :: pcm
    integer :: k
    logical, save,dimension(-1:1) :: openFlag=.true.
    character(40), dimension(-1:1) :: fileName

    call matrix_NN_NNPion(srts,matrixElements)

    rhoTwo=0.

    gamma=0.
    Select Case(chargePion)
    Case(-1)
       If(sum(chargeNucleon).eq.2) then           !pi^- p p -> pn
          rhoTwo(2)=rhoProton**2/2.
          k=2
       else If(sum(chargeNucleon).eq.1) then      ! pi^- pn -> nn
          rhoTwo(4)=rhoProton*rhoNeutron
          k=4
       else
          gamma=0
          return
       end if
    Case(0)
       If(sum(chargeNucleon).eq.2) then           ! pi^0 pp -> pp
          rhoTwo(1)=(rhoProton**2)/2.
          k=1
       else if(sum(chargeNucleon).eq.0) then      ! pi^0 nn -> nn
          rhoTwo(1)=(rhoNeutron**2)/2.
          k=1
       else if(sum(chargeNucleon).eq.1) then      ! pi^0 pn -> pn
          rhoTwo(3)=rhoProton*rhoNeutron
          k=3
       end if

    Case(1)
       If(sum(chargeNucleon).eq.0) then           ! pi^+ nn -> pn
          rhoTwo(2)=rhoNeutron**2/2.
          k=2
       else If(sum(chargeNucleon).eq.1) then      ! pi^+ pn -> pp
          rhoTwo(4)=rhoProton*rhoNeutron
          k=4
       else
          gamma=0
          return
       end if
    end select

    pcm=SQRT(Max(0.,srts**2/4.-mN**2))
    

    gamma=Min(100000.,gamma+matrixElements(k)*pcm/ srts/4. /pi*rhoTwo(k)/8./ePion/eNucleon(1)/eNucleon(2)*(0.197)**6*2.57672) 
    !factor 2.57672 because matrixs contains still the mb of cross section

     ! Make output
     IF(Present(outputFlag)) then
        If(OutPutFlag) then
           fileName(1)='gamma_PiPlus.dat'
           fileName(0)='gamma_PiNull.dat'
           fileName(-1)='gamma_PiMinus.dat'
           If(Openflag(chargePion)) then
              Open(11, File=fileName(chargePion))
              openFlag(chargePion)=.false.
           else
              Open(11, File=fileName(chargePion),position='Append')
           end if
           write(11,'(5F15.3,2I5)') ePion-mPi,eNucleon,srts,gamma,chargeNucleon
           close(11)
        end if
     end IF
   end function gamma_NNPion_NN

   !*********************************************************************************************


  subroutine gamma_NNPion_NN_output(srts,Epion,rhoProton, rhoNeutron,chargePion,gamma,OutPutFlag)
    !****s* NNPION_NN/gamma_NNPion_NN_output
    ! NAME
    ! subroutine gamma_NNPion_NN_output(srts,Epion,rhoProton, rhoNeutron,chargePion,gamma,OutPutFlag)
    ! NOTES
    ! Evaluates Gamma for the channels:
    ! * 1: p p pi^0-> p p  ; which is equivalent to n n -> n n pi^0 
    ! * 2: n n pi^+-> n p  ; which is equivalent to  p n -> p p pi^-
    ! * 3: p n pi^0-> p n  ; 
    ! * 4: p n pi^+-> p p  ; which is equivalent to n n -> p n pi^-
    ! The nucleons are considered to rest.
    ! This routine is mainly suited for output.
    ! OUTPUT
    ! * real, dimension(1:4),intent(out) :: gamma
    ! INPUTS
    ! * real, intent(in) :: srts ! SQRT(s)
    ! * real, intent(in) :: Epion ! Energy of pion
    ! * real, intent(in) :: rhoProton 
    ! * real, intent(in) :: rhoNeutron
    ! * integer, intent(in) :: chargePion
    ! * logical, optional,intent(in) :: outputFlag  -  if .true. then results are written to files "gamma_NNPiPlus_NN.dat","gamma_NNPiNull_NN.dat","gamma_NNPiMinus_NN.dat"
    !***

    use constants, only : pi, mN, mPi

    real, intent(in) :: srts ! SQRT(s)
    real, intent(in) :: Epion ! Energy of pion
    real, intent(in) :: rhoProton 
    real, intent(in) :: rhoNeutron
    integer, intent(in) :: chargePion
    real, dimension(1:4),intent(out) :: gamma
    logical, optional,intent(in) :: outputFlag

    real, dimension(1:4) ::matrixelements,rhoTwo
    real :: pcm
    integer :: k
    logical, save,dimension(-1:1) :: openFlag=.true.

    character(40), dimension(-1:1) :: fileName
    real, dimension(1:2) :: Enucleon !  Energies of nucleons


    call matrix_NN_NNPion(srts,matrixElements)
    
    eNucleon=mN  ! nucleons rest

    Select Case(chargePion)
       Case(-1)
       rhoTwo(1)=0. ! channel not open
       rhoTwo(2)=rhoProton**2/2.
       rhoTwo(3)=0. ! channel not open
       rhoTwo(4)=rhoProton*rhoNeutron
    
       Case(0)
       rhoTwo(1)=(rhoProton**2+rhoNeutron**2)/2.
       rhoTwo(2)=0. ! channel not open
       rhoTwo(3)=rhoProton*rhoNeutron
       rhoTwo(4)=0. ! channel not open
    
       Case(1)
       rhoTwo(1)=0. ! channel not open
       rhoTwo(2)=rhoNeutron**2/2.
       rhoTwo(3)=0. ! channel not open
       rhoTwo(4)=rhoProton*rhoNeutron
    end select

    pcm=SQRT(Max(0.,srts**2/4.-mN**2))
    
    do k=1,4
       gamma(k)=Min(100000.,matrixElements(k)*pcm/ srts/4.   /pi*rhoTwo(k)/8./ePion/eNucleon(1)/eNucleon(2)*(0.197)**6*2.57672) 
       !factor 2.57672 because matrixs contains still the mb of cross section
    end do

     ! Make output
     IF(Present(outputFlag)) then
        If(OutPutFlag) then
           fileName(1)='gamma_NNPiPlus_NN.dat'
           fileName(0)='gamma_NNPiNull_NN.dat'
           fileName(-1)='gamma_NNPiMinus_NN.dat'
           If(Openflag(chargePion)) then
              Open(11, File=fileName(chargePion))
              openFlag(chargePion)=.false.
              write(11,'(A)') '# Ekin of pion,mass of nucleon, energy of Nucleon, srts,gamma(1:4), sum(gamma)'
              write(11,'(A)') '# gamma(1): p p pi^0-> p p  +  n n pi^0 -> nn' 
              write(11,'(A)') '# gamma(2): n n pi^+-> n p  ; which is equivalent to  p p pi^- ->pn '
              write(11,'(A)') '# gamma(3): p n pi^0-> p n   '
              write(11,'(A)') '# gamma(4): p n pi^+-> p p  ; which is equivalent to p n pi^- -> nn'
           else
              Open(11, File=fileName(chargePion),position='Append')
           end if
           write(11,'(9F15.6)') ePion-mPi,eNucleon,srts,gamma, Sum(gamma)
           close(11)
        end if
     end IF
   end subroutine gamma_NNPion_NN_output

   !***********************************************************************

  subroutine matrix_NN_NNPion(srts,matrixElements,OutPutFlag)
    !****s* NNPION_NN/matrix_NN_NNPion
    ! NAME
    ! subroutine matrix_NN_NNPion(srts,matrixElements)
    ! PURPOSE
    ! Implements the Matrix elements for N N Pion -> N N in units of  "GeV**2 * mB"  for the channels:
    ! * 1: p p -> p p pi^0 ; which is equivalent to n n -> n n pi^0 
    ! * 2: p n -> n n pi^+ ; which is equivalent to  p n -> p p pi^-
    ! * 3: p n -> p n pi^0 ; 
    ! * 4: p p -> p n pi^+ ; which is equivalent to n n -> p n pi^-
    ! INPUTS
    ! * real, intent(in) ::srts
    ! OUTPUT
    ! * real, dimension(1:4), intent(out) ::matrixelements
    !***
    use barBar_barbarMes, only: NN_NNpi_direct
    use threeBodyPhaseSpace, only : Integrate_3bodyPS
    use constants, only : pi, mN, mPi

    real, intent(in) ::srts
    real, dimension(1:4), intent(out) ::matrixelements
    logical, optional, intent(in) ::outputFlag

    real, dimension (1:4) :: Sigma_NNPion
    real :: ps, pCM, factor
    integer :: k
    logical, save :: openFlag=.true.

    ! detailed balance: pi N N -> N N
    ! (see p.46 of Effenbergers diploma thesis)
    ! or page 56 of Buss' thesis
     ps = Integrate_3bodyPS (srts, mN, mN, mPi)
     !*ps=\int dm12 dm23 /s (!) 

     Sigma_NNPion = NN_NNpi_direct (srts)
     ! Field Sigma_NNPion:
     ! 1: p p -> p p pi^0 ; which is equivalent to n n -> n n pi^0 
     ! 2: p n -> n n pi^+ ; which is equivalent to  p n -> p p pi^-
     ! 3: p n -> p n pi^0 ; 
     ! 4: p p -> p n pi^+ ; which is equivalent to n n -> p n pi^-

     pcm=SQRT(Max(0.,srts**2/4.-mN**2))
     
     if(ps.gt.10**(-12)) then
        do k=1,4
           !*factor for identical particles
           factor=1.
           if(k.eq.2) factor=2.
           if(k.eq.4) factor=0.5
           matrixElements(k)=factor*64.*(2.*pi)**3*pCM*srts*Sigma_NNPion(k)/ps !units GeV**2 mB
        end do
     else
        do k=1,4
           matrixElements(k)=0
        end do
     end if
     
     ! Make output
     IF(Present(outputFlag)) then
        If(OutPutFlag) then
           If(Openflag) then
              Open(11, File='matrixelements_NN_NNPION.dat')
              openFlag=.false.
           else
              Open(11, File='matrixelements_NN_NNPION.dat',position='append')
           end if
           write(11,'(5F25.9)') srts, matrixElements
           close(11)
        end if
     end IF

   end subroutine matrix_NN_NNPion

end module NNPION_NN
