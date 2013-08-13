      subroutine proto_star(ZMBAR,RBAR,mass1,mass2,ECC,SEMI)
*
*       Pre-mainsequence binary evolution (Kroupa MN 277, 1491, 1995).
*       --------------------------------------------------------------
*
      implicit none
      real*8          mass1,mass2,ecc,semi,period
      real*8          ecc_initial,period_initial
      real*8          qnew,qold,mtot,Ro,mtot_initial
      real*8          R_periastron,alpha,beta
      real*8          ZMBAR,RBAR,au,Rsun
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
* astr. unit, solar radius, all in AU (1pc=206259.591AU)
      parameter(au=206259.591D0,Rsun=4.6523D-3)
*
* At this stage we have the masses, eccentricity and period of each binary
* at "birth", i.e. prior to circularisation and "feeding". Now evolve these
* to very, very roughly take into account complete circularisation,
* partial circularisation and "feeding". Do this if option KZ(41)=1:
* (i.e. mass-exchange at proto-stellar time):
*
* Define the best model: (alpha==lambda, beta==chi).
      alpha = 28.D0
      beta = 0.75D0
*
* in Msun:
      mtot = (mass1+mass2)*ZMBAR
      mtot_initial = mtot
* in AU:
      semi = semi*RBAR*au
* in years:
      period = semi*semi*semi/mtot
      period = DSQRT(period)
      ecc_initial = ecc
      period_initial = period
*
* 1) Circularisation and evolution of orbit as a function of periastron.
* Note that the algorithm used here leads to circularised orbits for
* logP<=1 approximately!! (if beta=1.5,alpha=35 approximately)
      Ro = alpha *Rsun
      R_periastron =  semi*(1.D0-ecc)
      alpha = -1.D0*(Ro/R_periastron)**beta
      if (ecc.GT.0.D0) then
         ecc = DEXP(alpha + DLOG(ecc))
      else
         ecc = ecc_initial
      end if
*
* 2) Change mass-ratio towards unity as a function of initial periastron.
      qold = mass1/mass2
      if (qold.GT.1.D0) qold = 1.D0/qold
      alpha = -1.D0*alpha
      if (alpha.GT.1.D0) then
         qnew = 1.D0
      else
         qnew = qold + (1.D0-qold) * alpha
      end if
*
* Set new masses in model units (remembering q=m1/m2<1) if mass is conserved.
*      mtot = mtot/ZMBAR
*      mass1 = mtot/(qnew+1.D0)
*      mass2 = mtot_initial/ZMBAR - mass1
*
* Keep the mass of primary fixed and adjust mass of secondary. Note that this
* algorithm leads to a gain in mass of the binary, hence of whole cluster.
*
C Added 20.06.96 write statements: added 3.09.2002: write to feeding.dat.
*       write(55,'(2(4F10.3))')
*    &    mass1*ZMBAR,mass2*ZMBAR,ecc_initial,
*    &    LOG10(period_initial*365.25),
*    &    DMAX1(mass1,mass2)*ZMBAR,qnew*mass1*ZMBAR,ecc,
*    &    LOG10(period*365.25) 
*       write(6,*)
        if(rank.eq.0)then
        write(6,*)' FEEDING in binpop_pk.f'
        write(6,'(a,2F8.3)')' old masses [Msun]:', 
     +                        mass1*ZMBAR,mass2*ZMBAR 
        end if
        mass1 = DMAX1(mass1,mass2)
        mass2 = qnew*mass1
        if(rank.eq.0)
     +  write(6,'(a,2F8.3)')' new masses [Msun]:', 
     +                        mass1*ZMBAR,mass2*ZMBAR 
C End added bit.
*
* In Msun:
        mtot = (mass1+mass2)*ZMBAR
*
C This below is wrong as in ecc formula above constant Rperi was assumed!
c* Duquennoy et al. 1992 in "Binaries as tracers of stellar evolution":
c      period = period_initial * DEXP((57.D0/14.D0) *
c     & (ecc*ecc - ecc_initial*ecc_initial))
C This below is correct:
       period = period_initial*((1.D0-ecc_initial)/(1.D0-ecc))**1.5D0
       period = period * DSQRT(mtot_initial/mtot)
*
* Form new semi-major axis and convert back to model units.
      semi = mtot * period*period
      semi = semi**(1.D0/3.D0)
      semi = semi/(RBAR*au)
*
      return
      end
