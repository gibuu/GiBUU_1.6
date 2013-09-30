      subroutine SetSwitchPythiaHermes(inFlag)
      IMPLICIT NONE
      logical inFlag

      common /DataSwitchPythiaHermes/ UseHermes
      logical UseHermes
      save /DataSwitchPythiaHermes/

      data UseHermes /.FALSE./

      UseHermes = inFlag

c$$$      if (UseHermes) then
c$$$         write(*,*) '## PYTHIA with HERMES tuning...'
c$$$      else
c$$$         write(*,*) '## PYTHIA w/o  HERMES tuning...'
c$$$      endif
      end


      logical function GetSwitchPythiaHermes()
      IMPLICIT NONE

      common /DataSwitchPythiaHermes/ UseHermes
      logical UseHermes
      save /DataSwitchPythiaHermes/

c$$$      write(*,*) '####### GetSwitchPythiaHermes:',UseHermes

      GetSwitchPythiaHermes = UseHermes

      return
      end


      SUBROUTINE PYDIFF
      IMPLICIT NONE

      common /DataSwitchPythiaHermes/ UseHermes
      logical UseHermes
      save /DataSwitchPythiaHermes/

      if (UseHermes) then
         call PYDIFF_hermes
      else
         call PYDIFF_orig
      endif
      end


      SUBROUTINE PYGAGA(IGAGA,WTGAGA)
      IMPLICIT NONE
      integer iGAGA
      double precision WTGAGA

      common /DataSwitchPythiaHermes/ UseHermes
      logical UseHermes
      save /DataSwitchPythiaHermes/

      if (UseHermes) then
         call PYGAGA_hermes(IGAGA,WTGAGA)
      else
         call PYGAGA_modif(IGAGA,WTGAGA) ! modif !
      endif
      end


      SUBROUTINE PYSIGH(NCHN,SIGS)
      IMPLICIT NONE
      integer NCHN
      double precision SIGS

      common /DataSwitchPythiaHermes/ UseHermes
      logical UseHermes
      save /DataSwitchPythiaHermes/

      if (UseHermes) then
         call PYSIGH_hermes(NCHN,SIGS)
      else
         call PYSIGH_modif(NCHN,SIGS)
      endif
      end


      SUBROUTINE PYXTOT
      IMPLICIT NONE

      common /DataSwitchPythiaHermes/ UseHermes
      logical UseHermes
      save /DataSwitchPythiaHermes/

      if (UseHermes) then
         call PYXTOT_hermes
      else
         call PYXTOT_orig
      endif
      end

