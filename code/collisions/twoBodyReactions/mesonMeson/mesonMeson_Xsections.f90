!***************************************************************************
!****m* /mesonMeson_Xsections
! NAME
! module mesonMeson_Xsections
! PURPOSE
! Provide parametrizations for (some) meson-meson cross sections
!***************************************************************************

module mesonMeson_Xsections

  implicit none
  Private

  PUBLIC :: sig_pipi2kkbar
  PUBLIC :: kkbar_cross
  PUBLIC :: kkbar_out
  PUBLIC :: kstarkbar_cross
  PUBLIC :: kstarkbar_out

  logical, save :: debug = .false.

contains

  !*************************************************************************
  !****f* mesonMeson_Xsections/sig_pipi2kkbar
  ! NAME
  ! real function sig_pipi2kkbar(srts,srts0)
  ! PURPOSE
  ! Cross section of pi pi --> K Kbar averaged over isospins of
  ! incoming pions and summed over isospins of outgoing pions
  !
  ! INPUTS
  ! * real:: srts  -- c.m. energy of incoming particles (GeV),
  ! * real:: srts0 -- threshold c.m. energy (GeV),
  ! OUTPUT
  ! * sig_pipi2kkbar -- cross section (mbarn)
  !*************************************************************************

  real function sig_pipi2kkbar(srts,srts0)
    real, intent(in) :: srts,srts0
    if(srts.le.srts0) then
       sig_pipi2kkbar = 0.
       return
    endif
    !   Parameterization from W. Cassing et al., NPA 614 (1997) 415:
    sig_pipi2kkbar = 2.7*(1.-(srts0/srts)**2)**0.76
  end function sig_pipi2kkbar


  !*************************************************************************
  !****s* mesonMeson_Xsections/kkbar_cross
  ! NAME
  ! subroutine kkbar_cross(izt,srts,pinitial2,const,msigbg)
  ! PURPOSE
  ! calculate the cross-sections for the KKbar incoming channel
  !
  ! INPUTS 
  ! * integer:: izt -- total charge of K and Kbar,
  ! * real   :: srts -- c.m. energy (GeV),
  ! * real   :: pinitial2 -- c.m. momentum squared of incoming particles (GeV/c)**2
  ! * real   :: const -- cross section of nonstrange+nonstrange --> K Kbar
  !   or --> Kstar Kbar, K Kbarstar (mbarn),  
  ! OUTPUT
  ! * real, dimension(21) :: msigbg(i), i=1,2,..,21 --- partial cross sections for different outgoing channels (mbarn)
  !*************************************************************************
  subroutine kkbar_cross(izt,srts,pinitial2,const,msigbg)
      
    use IdTable
    use particleProperties, only : hadron
    use constants, only: mPi, mK

    integer, intent(in) :: izt
    real, intent(in) :: srts,pinitial2,const
    real, dimension(21), intent(out) :: msigbg
    
    real :: srts0,sfree,pfinal2,isofac,spinfac
    
    if(debug) write(*,*) 'In kkbar_cross:', srts
    
    !     fritiof might produce lower mass kaons
    if(srts.lt.2.*mK) then 
       srts0 = 2.*0.493
    else
       srts0 = 2.*mK
    endif

    if(srts0.gt.srts) then
       write(*,*)'srts0 gt srts in mesmes',srts0,srts
       stop
    endif

    sfree = srts**2

    !     KKbar -> pipi:
    isofac = 0.5
    pfinal2 = sfree/4. - mPi**2
    msigbg(1) = isofac*9./4.*sig_pipi2kkbar(srts,srts0)*pfinal2/pinitial2
      
    !     KKbar -> rhorho:
    !     (same isofac as before but additional spinfactor)     
    if(srts.gt.2.*hadron(rho)%mass) then
       spinfac = 9.
       pfinal2 = sfree/4. - hadron(rho)%mass**2
       msigbg(2) = spinfac*isofac*9./4.*sig_pipi2kkbar(srts,srts0) &
            *pfinal2/pinitial2
    else
       msigbg(2) = 0.
    endif

    !     KKbar -> pi rho:
    !**** (this is possible only in p-wave, hence for simplicity put to zero)
    msigbg(3) = 0.

    !     KKbar -> pi eta:
    if(izt.eq.0) then
       isofac = 0.5
    else
       isofac = 1.
    endif
    if(srts.gt.mPi+hadron(eta)%mass) then
       pfinal2 = (sfree + mPi**2 - hadron(eta)%mass**2)**2 &
            /(4.*sfree) - mPi**2
       msigbg(4) = isofac*const*pfinal2/pinitial2
    else
       msigbg(4) = 0.
    endif

    !     KKbar -> pi sigma, because of parity:
    msigbg(5) = 0.
    

    !     KKbar -> pi omega, because of p-wave:
    msigbg(6) = 0.
    
    !     KKbar -> pi etap:
    if(srts.gt.mPi+hadron(etaPrime)%mass) then
       pfinal2 = (sfree + mPi**2 - hadron(etaPrime)%mass**2)**2 &
            /(4.*sfree) - mPi**2
       msigbg(7) = isofac*const*pfinal2/pinitial2
    else
       msigbg(7) = 0.
    endif
    
    !     KKbar -> eta eta:
    if(srts.gt.2.*hadron(eta)%mass.and.izt.eq.0) then
       pfinal2 = sfree/4. - hadron(eta)%mass**2 
       !****    additional factor 0.5 is due to identical particles in fin. st.:
       msigbg(8) = 0.5*isofac*const*pfinal2/pinitial2
    else
       msigbg(8) = 0.
    endif

    !     KKbar -> eta rho (p-wave):
    msigbg(9) = 0.
    

    !     KKbar -> eta sigma, because of parity:
    msigbg(10) = 0.

    !     KKbar -> eta omega (p-wave):
    msigbg(11) = 0.


    !     KKbar -> eta etap:
    if(srts.gt.hadron(eta)%mass+hadron(etaPrime)%mass.and.izt.eq.0) then
       pfinal2 = (sfree + hadron(eta)%mass**2 - hadron(etaPrime)%mass**2)**2 &
            /(4.*sfree) - hadron(eta)%mass**2
       msigbg(12) = isofac*const*pfinal2/pinitial2
    else
       msigbg(12) = 0.
    endif
    
    !     KKbar -> rho sigma (p-wave):
    msigbg(13) = 0.
    
    !     KKbar -> rho omega:
    if(srts.gt.hadron(rho)%mass+hadron(omegaMeson)%mass) then
       spinfac = 9.
       pfinal2 = &
            (sfree + hadron(rho)%mass**2 - hadron(omegaMeson)%mass**2)**2 &
            /(4.*sfree) - hadron(rho)%mass**2
       msigbg(14) = spinfac*isofac*const*pfinal2/pinitial2
    else
       msigbg(14) = 0.
    endif
    
    !     KKbar -> rho etap (p-wave):
    msigbg(15) = 0.
    
    !     KKbar -> sigma sigma:
    if(srts.gt.2.*hadron(sigmaMeson)%mass.and.izt.eq.0) then
       pfinal2 = sfree/4. - hadron(sigmaMeson)%mass**2
       ! factor 0.5 is due to identical particles in fin. st.
       msigbg(16) = 0.5*isofac*const*pfinal2/pinitial2
    else
       msigbg(16) = 0.
    endif
    
    !     KKbar -> sigma omega (p-wave):
    msigbg(17) = 0.
    
    !     KKbar -> sigma etap (parity):
    msigbg(18) = 0.
    
    !     KKbar -> omega omega:
    if(srts.gt.2.*hadron(omegaMeson)%mass.and.izt.eq.0) then
       spinfac = 9.
       pfinal2 = sfree/4. - hadron(omegaMeson)%mass**2
       ! factor 0.5 is due to identical particles in fin. st.
       msigbg(19) = 0.5*spinfac*isofac*const*pfinal2/pinitial2
    else 
       msigbg(19) = 0.
    endif
    
    !     KKbar -> omega etap (p-wave):
    msigbg(20) = 0.
    
    !     KKbar -> etap etap:
    if(srts.gt.2.*hadron(etaPrime)%mass.and.izt.eq.0) then
       pfinal2 = sfree/4. - hadron(etaPrime)%mass**2
       !****    factor 0.5 is due to identical particles in fin. st.
       msigbg(21) = 0.5*isofac*const*pfinal2/pinitial2
    else
       msigbg(21) = 0.
    endif
    
  end subroutine kkbar_cross
  
  !*************************************************************************
  !****s* mesonMeson_Xsections/kkbar_out
  ! NAME
  ! subroutine kkbar_out(izt,msigbg,sigbgt,teilchenOut)
  ! PURPOSE
  ! choose outgoing state for the K Kbar annihilation
  ! INPUTS 
  ! * integer :: izt --- total charge of K and Kbar,
  ! * real, dimension(:) ::msigbg(i), i=1,2,...,21 --- partial cross sections
  !   for different outgoing channels (mbarn),
  ! * real ::sigbgt --- total background cross section (mbarn)
  ! OUTPUT 
  ! * type(preEvent),dimension(21) :: teilchenOut   ---   outgoing particles
  !*************************************************************************
  subroutine kkbar_out(izt,msigbg,sigbgt,teilchenOut)

    use IdTable
    use random, only : rn
    use preEventDefinition

    integer, intent(in) :: izt
    real, dimension(21), intent(in) :: msigbg
    real, intent(in) :: sigbgt
    type(preEvent),dimension(:), intent(out) :: teilchenOut

    logical :: flag
    real :: msig,x,xx
    integer :: mch,n1,n2
    
    !     determine outgoing channel
    flag = .true.
    mch = 0
    msig = 0.
    x = rn()
    do while (flag)
       mch = mch + 1
       if(mch.gt.21) then
          write(*,*)'problems in mesmes (kkbar_out)'
          stop
       end if
       msig = msig + msigbg(mch)
       if(x.le.msig/sigbgt) flag = .false.
    end do
    
    select case(mch)
    case( 1)
       teilchenOut(1:2)% ID = (/pion, pion/)
    case( 2)
       teilchenOut(1:2)% ID = (/rho, rho/)
    case( 3)
       teilchenOut(1:2)% ID = (/pion, rho/)
    case( 4)
       teilchenOut(1:2)% ID = (/pion, eta/)
    case( 5)
       teilchenOut(1:2)% ID = (/pion, sigmaMeson/)
    case( 6)
       teilchenOut(1:2)% ID = (/pion, omegaMeson/)
    case( 7)
       teilchenOut(1:2)% ID = (/pion, etaPrime/)
    case( 8)
       teilchenOut(1:2)% ID = (/eta, eta/)
    case( 9)
       teilchenOut(1:2)% ID = (/rho, eta/)
    case(10)
       teilchenOut(1:2)% ID = (/eta, sigmaMeson/)
    case(11)
       teilchenOut(1:2)% ID = (/eta, omegaMeson/)
    case(12)
       teilchenOut(1:2)% ID = (/eta, etaPrime/)
    case(13)
       teilchenOut(1:2)% ID = (/rho, sigmaMeson/)
    case(14)
       teilchenOut(1:2)% ID = (/rho, omegaMeson/)
    case(15)
       teilchenOut(1:2)% ID = (/rho, etaPrime/)
    case(16)
       teilchenOut(1:2)% ID = (/sigmaMeson, sigmaMeson/)
    case(17)
       teilchenOut(1:2)% ID = (/sigmaMeson, omegaMeson/)
    case(18)
       teilchenOut(1:2)% ID = (/sigmaMeson, etaPrime/)
    case(19)
       teilchenOut(1:2)% ID = (/omegaMeson, omegaMeson/)
    case(20)
       teilchenOut(1:2)% ID = (/omegaMeson, etaPrime/)
    case(21)
       teilchenOut(1:2)% ID = (/etaPrime, etaPrime/)
    end select


    !     Determine charges:

    if(mch.le.3) then

       !       outgoing pipi, rhorho, pirho:
       !       (for pirho this, actually, will not happen)

       !       in order to distribute charge equally : random decision
       x = rn()
       if(x.lt.0.5) then
          n1 = 1
          n2 = 2
       else
          n2 = 1
          n1 = 2
       endif
       
       if(izt.eq.0) then
          xx = rn()
          if(xx.lt.5./6.) then
             teilchenOut(n1)%charge = 1
             teilchenOut(n2)%charge = -1
          else
             teilchenOut(n1)%charge = 0
             teilchenOut(n2)%charge = 0
          endif
       else
          teilchenOut(n1)%charge = izt
          teilchenOut(n2)%charge = 0           
       endif
       
    else
       
       !       in all other cases:
       teilchenOut(1)%charge = izt
       teilchenOut(2)%charge = 0
       
    endif
    
  end subroutine kkbar_out
      
  !*************************************************************************
  !****s* mesonMeson_Xsections/kstarkbar_cross
  ! NAME
  ! subroutine kstarkbar_cross(izt,srts,pinitial2,const,msigbg)
  ! PURPOSE
  ! calculate cross-sections for Kstar Kbar or K Kstarbar incoming channels
  !
  ! INPUTS 
  ! * integer :: izt --- total charge of Kstar and Kbar,
  ! * real :: srts --- c.m. energy (GeV),
  ! * real :: pinitial2 --- c.m. momentum squared of incoming particles (GeV/c)**2,
  ! * real :: const --- cross section of nonstrange+nonstrange --> K Kbar
  !   or --> Kstar Kbar, K Kbarstar (mbarn),  
  ! OUTPUT 
  ! * real, dimension(:) :: msigbg(i), i=1,2,..,8 --- partial cross sections
  !   for different outgoing channels (mbarn)
  !*************************************************************************
  subroutine kstarkbar_cross(izt,srts,pinitial2,const,msigbg)

    use IdTable
    use particleProperties, only : hadron
    use constants, only: mPi, mK

    integer, intent(in) :: izt
    real, intent(in) :: srts,pinitial2,const
    real, dimension(:), intent(out) :: msigbg

    real :: srts0,sfree,pfinal2,isofac!,spinfac

    if(debug) write(*,*) 'In kstarkbar_cross:', srts

    !     fritiof might produce lower mass kaons
    if(srts.lt.2.*mK) then 
       srts0 = 2.*0.493
    else
       srts0 = 2.*mK
    endif

    if(srts0.gt.srts) then
       write(*,*)'srts0 gt srts in mesmes',srts0,srts
       stop
    endif

    sfree = srts**2

    !     Kbar Kstar->pi rho:
    if(srts.gt.mPi+hadron(rho)%mass) then
       isofac = 1.
       pfinal2 = (sfree + mPi**2 - hadron(rho)%mass**2)**2 &
            /(4.*sfree) - mPi**2
       msigbg(1) = isofac*9./4.*sig_pipi2kkbar(srts,srts0)*pfinal2/pinitial2
    else
       msigbg(1) = 0.
    endif

    !     Kbar Kstar->pi omega:
    if(izt.eq.0) then
       isofac = 0.25
    else
       isofac = 0.5
    endif
    if(srts.gt.mPi+hadron(omegaMeson)%mass) then
       pfinal2 = (sfree + mPi**2 - hadron(omegaMeson)%mass**2)**2&
            /(4.*sfree) - mPi**2
       msigbg(2) = isofac*const*pfinal2/pinitial2
    else
       msigbg(2) = 0.
    endif

    !     Kbar Kstar->eta rho:
    if(srts.gt.hadron(eta)%mass+hadron(rho)%mass) then
       pfinal2 = (sfree + hadron(eta)%mass**2 - hadron(rho)%mass**2)**2 &
            /(4.*sfree) - hadron(eta)%mass**2
       msigbg(3) = isofac*const*pfinal2/pinitial2
    else
       msigbg(3) = 0.
    endif

    !     Kbar Kstar->eta omega:
    if(izt.eq.0.and.srts.gt.hadron(eta)%mass+hadron(omegaMeson)%mass) then
       pfinal2 = (sfree + hadron(eta)%mass**2 - hadron(omegaMeson)%mass**2)**2&
            /(4.*sfree) - hadron(eta)%mass**2
       msigbg(4) = isofac*const*pfinal2/pinitial2
    else
       msigbg(4) = 0.
    endif

    !     Kbar Kstar->rho sigma (parity):
    msigbg(5) = 0.

    !     Kbar Kstar->rho etap:
    if(srts.gt.hadron(rho)%mass+hadron(etaPrime)%mass) then
       pfinal2 = (sfree + hadron(rho)%mass**2 - hadron(etaPrime)%mass**2)**2 &
            /(4.*sfree) - hadron(rho)%mass**2
       msigbg(6) = isofac*const*pfinal2/pinitial2
    else
       msigbg(6) = 0.
    endif

    !     Kbar Kstar->sigma omega (parity):
    msigbg(7) = 0.

    !     Kbar Kstar->etap omega:
    if(izt.eq.0.and.srts.gt.hadron(etaPrime)%mass+hadron(omegaMeson)%mass) then
       pfinal2 = (sfree + hadron(etaPrime)%mass**2 &
            - hadron(omegaMeson)%mass**2)**2/(4.*sfree) &
            - hadron(etaPrime)%mass**2
       msigbg(8) = isofac*const*pfinal2/pinitial2
    else
       msigbg(8) = 0.
    endif

  end subroutine kstarkbar_cross
  
  !*************************************************************************
  !****s* mesonMeson_Xsections/kstarkbar_out
  ! NAME 
  ! subroutine kstarkbar_out(izt,msigbg,sigbgt,teilchenOut)
  ! PURPOSE
  ! choose outgoing state for the Kstar Kbar or K Kstarbar annihilation
  ! INPUTS
  ! INPUT: 
  ! * integer ::izt --- total charge of incoming particles,
  ! * real, dimension(:) :: msigbg(i), i=1,2,...,8 --- partial cross sections
  !   for different outgoing channels (mbarn),
  ! * real :: sigbgt --- total background cross section (mbarn)
  ! OUTPUT
  ! * type(preEvent),dimension(:) :: teilchenOut   ---   outgoing particles
  !*************************************************************************
  subroutine kstarkbar_out(izt,msigbg,sigbgt,teilchenOut)

    use IdTable
    use random, only : rn
    use preEventDefinition

    integer, intent(in) :: izt
    real, dimension(:), intent(in) :: msigbg
    real, intent(in) :: sigbgt
    type(preEvent),dimension(:), intent(out) :: teilchenOut     ! produced particles

    logical :: flag
    real :: msig,x,xx
    integer :: mch,n1,n2

    !     determine outgoing channel
    flag = .true.
    mch = 0
    msig = 0.
    x = rn()
    do while (flag)
       mch = mch + 1
       if(mch.gt.8) then
          write(*,*)'problems in mesmes (kstarkbar_out)'
          stop
       end if
       msig = msig + msigbg(mch)
       if(x.le.msig/sigbgt) flag = .false.
    end do

    select case (mch)
    case( 1)
       teilchenOut(1:2)% ID = (/pion, rho/)
    case( 2)
       teilchenOut(1:2)% ID = (/pion, omegaMeson/)
    case( 3)
       teilchenOut(1:2)% ID = (/rho, eta/)
    case( 4)
       teilchenOut(1:2)% ID = (/eta, omegaMeson/)
    case( 5)
       teilchenOut(1:2)% ID = (/rho, sigmaMeson/)
    case( 6)
       teilchenOut(1:2)% ID = (/rho, etaPrime/)
    case( 7)
       teilchenOut(1:2)% ID = (/sigmaMeson, omegaMeson/)
    case( 8)
       teilchenOut(1:2)% ID = (/etaPrime, omegaMeson/)
    end select

    !     Determine charges:

    if(mch.eq.1) then

       !       outgoing pi rho:

       !       in order to distribute charge equally : random decision
       x = rn()
       if(x.lt.0.5) then
          n1 = 1
          n2 = 2
       else
          n2 = 1
          n1 = 2
       endif

       if(izt.eq.0) then
          xx = rn()
          if(xx.lt.5./6.) then
             teilchenOut(n1)%charge = 1
             teilchenOut(n2)%charge = -1
          else
             teilchenOut(n1)%charge = 0
             teilchenOut(n2)%charge = 0
          endif
       else
          teilchenOut(n1)%charge = izt
          teilchenOut(n2)%charge = 0           
       endif

    else

       !       in all other cases:
       teilchenOut(1)%charge = izt
       teilchenOut(2)%charge = 0

    endif

  end subroutine kstarkbar_out

end module mesonMeson_Xsections

