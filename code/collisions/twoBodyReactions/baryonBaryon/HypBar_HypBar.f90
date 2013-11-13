!*******************************************************************************
!****m* /hypNuc_hypNuc
! NOTES
! This module includes the calculation of hyperon-nucleon -> hyperon-nucleon
! cross sections. The properties of all final channels are calculated.
!*******************************************************************************
module hypNuc_hypNuc

  use preEventDefinition
  implicit none

  Private

  public :: hypNuc_hypNuc_Main
  public :: get_Channels_YN

  ! ID's & charges of final states.
  ! 1-index: Number of isospin channels
  ! 2-index: number of particle in 2-body final state (Baryon,Hyperon)
  Integer, dimension(1:3,1:2), SAVE :: IdsOut,ChargesOut
  ! YN->YN Xsection (in units of mb)

contains

  !*******************************************************************************
  !****f* hypNuc_hypNuc/hypNuc_hypNuc_Main
  ! NAME
  ! function hypNuc_hypNuc_Main (srts, teilchenIN) result (sigma_yn)
  ! PURPOSE
  ! * This routine prepares the final state for each isospin channel.
  ! INPUTS
  ! * type(particle), dimension (1:2), intent(in) :: teilchenIN  -- incoming particles
  ! * real,                            intent(in) :: srts        -- sqrt(s) [GeV]
  ! OUTPUT
  ! * see global variables of this module.
  ! NOTES 
  ! Included reactions:
  ! * Lambda Nucleon -> Lambda Nucleon
  ! * Lambda Nucleon -> Sigma  Nucleon
  ! * Sigma  Nucleon -> Sigma  Nucleon
  ! * Sigma  Nucleon -> Lambda Nucleon
  ! * Xi     Nucleon -> Xi     Nucleon
  ! * Xi     Nucleon -> Lambda Lambda
  ! * Xi     Nucleon -> Lambda Sigma
  !*******************************************************************************
  function hypNuc_hypNuc_Main (srts, teilchenIN) result (sigma_yn)
    use particleDefinition
    use particleProperties, only: hadron
    use IdTable, only: nucleon,Lambda,SigmaResonance,Xi
    use random, only: rn
    !-----------------------------------------------------------------------------
    real,           intent(in)                  :: srts
    type(particle), intent(in), dimension (1:2) :: teilchenIN
    real, dimension(1:3) :: sigma_yn
    !-----------------------------------------------------------------------------    
    integer :: inuc, ihyp
    real    :: M_Nuc,M_Hyp
    !-----------------------------------------------------------------------------
    ! initialize global variables:
    !-----------------------------------------------------------------------------
    sigma_yn(:) = 0.
    IdsOut(:,:) = 0
    ChargesOut(:,:) = 9999
    !-----------------------------------------------------------------------------
    ! Check input: only Lambda/Sigma/Xi + Nucleon scattering!
    !-----------------------------------------------------------------------------
    if ( hadron(teilchenIn(1)%ID)%strangeness==0) then
       inuc = 1
    else
       inuc = 2
    endif
    if (teilchenIN(inuc)%ID /= nucleon) return  ! No nucleon resonances!
    ihyp = teilchenIN(3-inuc)%ID
    if ( (ihyp > SigmaResonance) .and. (ihyp /= Xi) ) return
    !-----------------------------------------------------------------------------
    M_Nuc = teilchenIN(inuc)%mass
    M_Hyp = teilchenIN(3-inuc)%mass
    !-----------------------------------------------------------------------------
    ! elastic channel:
    IdsOut    (1,1:2) = teilchenIN(1:2)%ID
    ChargesOut(1,1:2) = teilchenIN(1:2)%charge
    !-----------------------------------------------------------------------------
    ! Prepare all isospin channels and get the PreEvent-vector:
    !-----------------------------------------------------------------------------
    Select Case (sum(teilchenIN(:)%ID))

       Case((Lambda+nucleon))

          ! Lambda proton:
          if (TeilchenIN(inuc)%charge==1) then

             ! Lambda proton-->Lambda proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,1)

             ! Lambda proton-->Sigma^+ neutron:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,2)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = SigmaResonance
             ChargesOut(2,inuc) = 0
             ChargesOut(2,3-inuc) = 1

             ! Lambda proton-->Sigma^0 proton:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,3)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 1
             ChargesOut(3,3-inuc) = 0

          ! Lambda neutron:
          else if (TeilchenIN(inuc)%charge==0) then

             ! Lambda neutron-->Lambda neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,4)

             ! Lambda neutron-->Sigma^- proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,5)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = SigmaResonance
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = -1

             ! Lambda neutron-->Sigma^0 neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,6)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 0

          endif

       Case((SigmaResonance+nucleon))

          ! Sigma^0 proton:
          if ( TeilchenIN(3-inuc)%charge==0 .and. TeilchenIN(inuc)%charge==1 ) then
             ! Sigma^0 proton-->Sigma^0 proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,7)

             ! Sigma^0 proton-->Lambda proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,8)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = Lambda
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = 0

             ! Sigma^0 proton-->Sigma^+ neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,9)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 1

          ! Sigma^0 neutron:
          else if ( TeilchenIN(3-inuc)%charge==0 .and. TeilchenIN(inuc)%charge==0 ) then
             ! Sigma^0 neutron-->Sigma^0 neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,10)

             ! Sigma^0 neutron-->Sigma^- proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,11)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = SigmaResonance
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = -1

             ! Sigma^0 neutron-->Lambda neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,12)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = Lambda
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 0

          ! Sigma^- proton:
          else if ( TeilchenIN(3-inuc)%charge==-1 .and. TeilchenIN(inuc)%charge==1 ) then
             ! Sigma^- proton-->Sigma^- proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,13)

             ! Sigma^- proton-->Lambda neutron:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,14)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = Lambda
             ChargesOut(2,inuc) = 0
             ChargesOut(2,3-inuc) = 0

             ! Sigma^- proton-->Sigma^0 neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,15)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 0

          ! Sigma^- neutron:
          else if ( TeilchenIN(3-inuc)%charge==-1 .and. TeilchenIN(inuc)%charge==0 ) then
             ! Sigma^- neutron-->Sigma^- neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,16)

          ! Sigma^+ proton:
          else if ( TeilchenIN(3-inuc)%charge==1 .and. TeilchenIN(inuc)%charge==1 ) then
             ! Sigma^+ proton-->Sigma^+ proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,17)

          ! Sigma^+ neutron:
          else if ( TeilchenIN(3-inuc)%charge==1 .and. TeilchenIN(inuc)%charge==0 ) then
             ! Sigma^+ neutron-->Sigma^+ neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,18)

             ! Sigma^+ neutron-->Lambda proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,19)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = Lambda
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = 0

             ! Sigma^+ neutron-->Sigma^0 proton:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,20)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 1
             ChargesOut(3,3-inuc) = 0

          endif

       Case((Xi+nucleon))
          ! Xi + Nucleon (I=0) --> Xi N & Lambda Lambda
          if ( ( TeilchenIN(3-inuc)%charge==-1 .and. TeilchenIN(inuc)%charge==1 ) .or. & 
               ( TeilchenIN(3-inuc)%charge==0  .and. TeilchenIN(inuc)%charge==0 ) ) then
             ! Xi + Nucleon (I=0) --> Xi + Nucleon 
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,21)

             ! Xi + Nucleon (I=0) --> Lambda Lambda
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,22)
             IdsOut(2,inuc) = Lambda
             IdsOut(2,3-inuc) = Lambda
             ChargesOut(2,inuc) = 0
             ChargesOut(2,3-inuc) = 0

          ! Xi + Nucleon (I=1) --> Xi N & Lambda Sigma
          else if ( TeilchenIN(3-inuc)%charge==0 .and. TeilchenIN(inuc)%charge==1 ) then
             ! Xi + Nucleon (I=1) --> Xi + Nucleon 
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,23)

             ! Xi + Nucleon (I=1) --> Lambda Sigma
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,24)
             if (rn() < 0.5) then
                IdsOut(2,inuc) = Lambda
                IdsOut(2,3-inuc) = SigmaResonance
                ChargesOut(2,inuc) = 0
                ChargesOut(2,3-inuc) = 1
             else
                IdsOut(2,3-inuc) = Lambda
                IdsOut(2,inuc) = SigmaResonance
                ChargesOut(2,3-inuc) = 0
                ChargesOut(2,inuc) = 1
             endif

             ! Note: I=-1 channel no theor. info; Xi^- neutron --> Sigma^- Lambda

          endif

       Case default 

          write(*,*) 'HypBar_HypBar_Main: undefined channel!', teilchenIN%ID
          STOP

       end Select

  end function hypNuc_hypNuc_Main


  !*******************************************************************************
  !****f* hypNuc_hypNuc/get_Channels_YN
  ! NAME
  ! function get_Channels_YN (InChan) result (ch)
  ! PURPOSE
  ! Returns channel infos (IDs and charges) as a PreEvent vector.
  ! INPUTS
  ! integer :: InChan -- Number of final channel
  !*******************************************************************************
  function get_Channels_YN (InChan) result (ch)
    integer, intent(in) :: InChan
    type(preEvent), dimension(1:2) :: ch
    
    ch(1:2)%ID = IdsOut(InChan,1:2)
    ch(1:2)%charge = ChargesOut(InChan,1:2)

  end function get_Channels_YN


  !*******************************************************************************
  !****f* hypNuc_hypNuc/xsectionYN
  ! NAME
  ! real function xsectionYN (srts, M_Nuc, M_Hyp, ichan)
  ! PURPOSE
  ! * This function calculates the hyperon-baryon -> hyperon baryon cross sections 
  !   for a particular isospin channel "ichan". 
  ! * The cross sections are given in mb.
  ! INPUTS
  ! * integer, intent(in) :: ichan  -- incoming particles
  ! * real,    intent(in) :: M_Nuc,M_Hyp -- masses of incoming particles
  ! * real,    intent(in) :: srts -- sqrt(s) [GeV]
  ! NOTES
  ! * Fits to experimentally known cross sections. 
  ! * Not for all channels the XS were experimentally known. The XS's for the unknown 
  !   channels were calculated using detailed balance or charge symmetry. 
  ! * For the important channels, where a Lambda is converted into a Sigma 
  !   (or a Sigma to a Lambda), the Xsections are known. 
  ! * References for exp. data on hyperon+nucleon scattering: 
  !   Landolt-Boernstein, New Series I/12b, p. 323: "Hyperon induced reactions"
  !   M.M. Nagels, T.A. Rijken, J.J. De Swart, Annals of Physics 79 (1973) 338.
  !   T.A. Rijken, Y. Yamamoto, Phys. Rev. C73 (2006) 044008.
  ! * Xi+N scattering cross sections: no exp. data available here; fits to theoretical 
  !   calculations: ESC04a model, 
  !   T.A. Rijken, Y. Yamamoto, nucl-th/0608074 (30.08.2006), Tables IX,X,XI & XII.
  !*******************************************************************************
  real function xsectionYN (srts, M_Nuc, M_Hyp, ichan)

    use particleProperties, only : hadron
    use IdTable, only : nucleon,Lambda, SigmaResonance
    use twoBodyTools, only : pCM
    use constants, only: mN

    real,    intent(in) :: srts, M_Nuc, M_Hyp
    integer, intent(in) :: ichan

    real :: plab,p_ab,p_cd,balFac,pmev

    if (ichan<1 .or. ichan>24) then
       write(*,*) 'hypNuc_hypNuc/xsection: ichan not defined! ichan = ',ichan
       STOP
    endif

    ! hyperon as beam-particle
    plab = sqrt( (srts**2-(m_hyp-m_nuc)**2)*(srts**2-(m_hyp+m_nuc)**2) )/(2.*m_nuc)
    pmev = plab*1000.

    p_ab = pCM (srts, M_Nuc, M_Hyp)  ! needed for detailed balance

    xsectionYN = 0.

    select case (ichan)

    ! Lambda/Sigma + Nucleon channels:

    case(1,4)

       ! Lambda + N --> Lambda + N   (parametrization by Oliver Arnold, following Rijken et al, PRC 59, 21-40, 1999)
       if (plab < 0.4) then
         xsectionYN = 203.56*exp(-14.47*plab**2) + 253.88*exp(-76.19*plab**2)
       else
          xsectionYN = 14.4*plab**(-0.12)
       endif

    case(2,3,5,6,8,12)

       ! Lambda + N --> Sigma0 + N  (ichan==3,6)
       !       if (plab < 0.436) then
       if (plab < 0.8) then
          xsectionYN = 30.*plab**(4.9)
       else if (plab > 0.8 .and. plab < 1.206) then
          xsectionYN = 5.7*plab**(-2.5) !orig.
       else if (plab > 1.206) then
          xsectionYN = min (19., 9.9456-10.254*plab+4.1135*plab**2)
       endif

       if (ichan==2 .or. ichan==5) then
         ! Lambda + Proton  --> Sigma+ + Neutron  (ichan==2)
         ! Lambda + Neutron --> Sigma- + Proton   (ichan==5)
         xsectionYN = 0.5*xsectionYN
       else if (ichan==8 .or. ichan==12) then
         ! Sigma0 + N --> Lambda + N
         p_cd = pCM(srts,hadron(lambda)%mass,mN)
         If (p_ab<1E-10) then
           write(*,*) 'WARNING: pInitial is zero in xsectionYN', p_ab
           balFac= 0.
         else  
           balFac= (p_cd/p_ab)**2        
         end if
         xsectionYN = xsectionYN*balFac ! detailed balance
       end if

    case(7,10,13,18)

       ! Sigma0 + N --> Sigma0 + N  (ichan==7,10)
       ! Sigma- + p --> Sigma- + p  (ichan==13)
       ! Sigma+ + n --> Sigma+ + n  (ichan==18)
       xsectionYN =  13.5*plab**(-1.25)

    case(9,11,15,20)

       ! Sigma- + Proton  --> Sigma0 + Neutron (ichan==15)
       ! Sigma+ + Neutron --> Sigma0 + Proton  (ichan==20)
       xsectionYN = 13.5*plab**(-1.25)

       if (ichan==9 .or. ichan==11) then
         ! Sigma0 + Proton  --> Sigma+ + Neutron (ichan==9)
         ! Sigma0 + Neutron --> Sigma- + Proton  (ichan==11)
         p_cd = pCM(srts,hadron(SigmaResonance)%mass,mN)
         If (p_ab<1E-10) then
           write(*,*) 'WARNING: pInitial is zero in xsectionYN', p_ab
           balFac= 0.
         else  
           balFac= (p_cd/p_ab)**2        
         end if
         xsectionYN = xsectionYN*balFac ! detailed balance
       end if

    case(14,19)

       ! Sigma- + Proton  --> Lambda + Neutron (ichan==14)
       ! Sigma+ + Neutron --> Lambda + Proton  (ichan==19)
       xsectionYN = 13.2*plab**(-1.18)

    case(16,17)

       ! Sigma- + Neutron --> Sigma- + Neutron (ichan==16)
       ! Sigma+ + Proton  --> Sigma+ + Proton  (ichan==17)
       xsectionYN = 38.*plab**(-0.62)

    ! Xi Nucleon channels:

    case(21)
       ! Xi N --> Xi N (I=0)
       xsectionYN = 17.3886*exp(-0.01*pmev)+2572.95*pmev**(-1.18381)

    case(22)
       ! Xi N --> Lambda Lambda (I=0)      
       xsectionYN = 417.996*exp(-0.00567813*pmev) - 340.722*pmev**(-0.359183) + 31.*pmev**0.02

    case(23)
       ! Xi N --> Xi N (I=1)
       xsectionYN = 328.083*exp(-0.00418763*pmev)+0.00602*pmev
             
    case(24)
       ! Xi N --> Sigma Lambda (I=1)
       if (pmev > 590.) xsectionYN = 4.626*(pmev/590.+1.)**(-0.604771)

    end select


! another parametrization for Xi N; however, not properly implemented. 
!    Case(2) 
!       if (ichan==21) then
!          xsectionYN = 60662.*pmev**(-1.5189)
!       else if (ichan==22) then
!          xsectionYN = min(300.,5.4451e+09*pmev**(-3.0405))
!       else if (ichan==23) then
!          xsectionYN = 16.91*exp(-0.00752*pmev)+pmev**0.2241
!       else if (ichan==24) then
!          if (pmev < 590.) then
!             xsectionYN = 0.0
!          else
!             xsectionYN = 6.746*(pmev/590.+1.)**(-1.02643)
!          endif
!       endif

  end function xsectionYN


end module hypNuc_hypNuc
