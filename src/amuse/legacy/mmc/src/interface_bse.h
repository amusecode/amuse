C------------------------------------------------------------
C
C       interface_bse.h 
C
C       Common for stellar evolution related variables. 
C
C------------------------------------------------------------
C
c      include "common.h"
      integer nmax
c      parameter (nmax=700000)
      parameter (nmax=2500000)
c      parameter (nmax=11000)
*
      integer kstar(0:nmax)
*
      integer kstarb(0:nmax)
      integer iprim(0:nmax),isec(0:nmax)
C
      real*8 age(0:nmax),epoch(0:nmax)
      real*8 ms_lifetime(0:nmax),standard_timestep(0:nmax)
      real*8 zmass0(0:nmax),zmass(0:nmax),zmassc(0:nmax)
      real*8 zinit(0:nmax),yinit(0:nmax)
      real*8 zcurr(0:nmax),ycurr(0:nmax)
      real*8 spin(0:nmax),radius(0:nmax),zlum(0:nmax)
      real*8 radc(0:nmax),renv(0:nmax),menv(0:nmax)
      real*8 rlof(0:nmax)
      real*8 xstar(3,0:nmax),vstar(3,0:nmax)
*
      real*8 sep(0:nmax),tb(0:nmax),ecc(0:nmax)
*
      character*20 label(0:15),labelb(6)

      data label / 'main_sequence (0)', 'main_sequence (1)', 
     &             'hertzsprung_gap (2)', 'giant_branch (3)',
     &             'core_helium (4)', 'first_agb (5)', 
     &             'second_agb (6)', 'helium_ms (7)', 
     &             'helium_hg (8)', 'helium_gb (9)', 
     &             'he_white_dwarf (10)', 'co_white_dwarf (11)',
     &             'one_white_dwarf (12)', 'neutron_star (13)', 
     &             'black_hole (14)', 'no_remnant (15)' /
      data labelb / 'detached', 'semi_detached', 
     &              'contact', 'merged', 
     &              'obliterated', 'disrupted' /
C
      common /stellar/ age,epoch,ms_lifetime,standard_timestep,
     &                 zmass0,zmass,zmassc,zinit,yinit,
     &                 zcurr,ycurr,spin,radius,radc,zlum,
     &                 menv,renv,rlof,xstar,vstar,kstar
      common /binary/ sep,tb,ecc,kstarb,iprim,isec
C
C------------------------------------------------------------
