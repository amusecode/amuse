************************************************************************
*
*      INTERFACE_BSE.F 
*
*      STELLAR AND BINARY EVOLUTION MODULE 
*      [comprised of interface_bse.f, interface_bse.h, const_bse.h, 
*       & libstr.a]
*      Version 0.0 March 20, 2003 by JH
c     Version 3.1 June 14, 2004 by DCH
*
*      Requires libstr.a - BSE library of stellar and binary evolution 
*      functions/prescriptions described in the papers: 
*
*      "Comprehensive Analytic Formulae for Stellar Evolution as a 
*       Function of Mass and Metallicity", 
*       Hurley J.R., Pols O.R. and Tout C.A. 2000, MNRAS, 315, 543. 
*       (the SSE paper). 
*
*      "Evolution of Binary Stars and the Effect of Tides on Binary Populations", 
*       Hurley J.R., Tout C.A. and Pols O.R. 2002, MNRAS, 329, 897
*       (the BSE paper). 
*
*      [note that libstr.a was compiled with either g77 or fort 
*       - the Compaq Fortran compiler - on Linux (or Mac OSX) 
*       and using the -O2 option in all cases]
*      [contact jhurley@amnh.org regarding any problems with the module]
*
*      The subroutines evStar and evBinary provide the link between the 
*      main program and the SSE/BSE library (or module). These routines 
*      evolve the star (or binary) forward by a user specified interval. 
*      The library routines return a recommended update timestep based on 
*      the evolution stage of the star, and the parameters dmmax (maximum 
*      allowed change in mass) and drmax (maximum allowed change in radius). 
*      The main program may or may not utilise this timestep. 
*
*      Input options for the SSE/BSE library are set in the subroutine initPar 
*      and these may be modified by the user (see the routine for an 
*      explanation of these options). The common blocks that convey these 
*      options to the library are declared in const_bse.h. 
*      Single stars are initialized by calling initStar and binaries by 
*      calling initBinary (which in turn initializes two single stars). 
*
*      The following quantities are available for star i: 
*
*          Age     - time to which star has been evolved (Myr)
*          MStime  - main sequence lifetime (Myr)
*          Tstep   - recommended time between updates (Myr)
*          Epoch   - effective zero-age of star (Myr)(age = Age - Epoch)
*          Mass    - current stellar mass (Msun)
*          Mass0   - initial stellar mass (Msun)
*          Radius  - current stellar radius (Rsun)
*          Lum     - current stellar luminosity (Lsun) 
*          Y0      - initial helium mass fraction (0->1)(Y = 0.24 + 2*Z)
*          Z0      - initial metallicity mass fraction (0->1)(0.02 is solar)
*          Y       - current helium mass fraction (0->1)
*          Z       - current metallicity mass fraction (0->1)
*          Mc      - core mass (Msun)
*          Rc      - core radius (Rsun)
*          Menv    - mass of the convective envelope (Msun)
*          Renv    - radius of the convective envelope (Msun)
*          Spin    - stellar rotation (1/yr)
*          Rl      - Roche-lobe radius (Rsun) 
*          Pos     - x,y,z co-ordinates (?)
*          Vel     - velocity in x,y,z (km/s)
*          Type    - index of stellar type (0->15, e.g. MS star = 1)
*
*      and these additional quantities are available for binaries: 
*
*          Btype   - index of binary type (1->6, e.g. detached = 1)
*          Iprim   - array index of the primary star (0->nmax)
*          Isec    - array index of the primary star (0->nmax)
*          Mprim   - mass of the primary (primary has greater radius/rl)
*          Msec    - mass of the secondary 
*          Semi    - semi-major axis of the orbit (Rsun) 
*          Ecc     - eccentricity of the orbit 
*          Tb      - orbital period (day)
*
*      Note that quantity X may be obtained by a call to the subroutine getX 
*      and may be set by a call to setX. The arrays that store these quantities 
*      are declared in interface_bse.h (where the user may choose to alter the 
*      size of these arrays, i.e. nmax). 
*
*      Additional subroutines included in this interface: 
*
*         ssupdatetime - provides next update time for star i 
*         bsupdatetime - calls ssupdatetime using index of primary star 
*         binaryexists - determines if binary remains bound based on Btype 
*         getLabel     - text label associated with index of stellar type 
*         getLabelb    - text label associated with index of binary type 
*         printstar    - formatted output for a single star 
*         printbinary  - formatted output for a binary 
* 
************************************************************************
*
       SUBROUTINE initPar
       implicit none
*
* Input options for the SSE/BSE library functions 
* (see BSE paper for details on most of the options). 
*
       include "const_bse.h"
       include "interface_bse.h"
       integer i,iflagns,iflagbh
       common /fflags/ iflagns,iflagbh
*
* Random number seed used to determine supernovae velocity kicks. 
* [NOTE THAT THIS SHOULD IDEALLY BE COMMUNICATED FROM MAIN PROGRAM]
       idum = -999
*
* Flag for use in common-envelope routine (not currently utilised).
       ceflag = 0
*
* Flag to activate tidal circularisation (0=off; 1=on).
       tflag = 1
*
* Flag to choose which WD IFMR to use (0=SSE; 1=HPE, 1995, MN, 272, 800).
       ifflag = 0
*
* Flag to determine NS/BH mass (0=SSE; 1=Belczynski et al. 2002, ApJ, 572, 407).
       nsflag = 1
*
* Flag to choose WD cooling track (0=simple-Mestel; 1="modified-Mestel").
       wdflag = 1
*
* Flag to allow kicks at birth for BHs (0=no; 1=yes).
       bhflag = iflagbh
*
* Maximum NS mass (set at 1.8 Msun for nsflag=0; 3.0 Msun for nsflag=1).
       mxns = 3.d0
*
* Reimers mass-loss coefficient (neta*4x10^-13; =0.5 normally).
       neta = 0.5d0
*
* Binary enhanced mass-loss parameter (inactive for single stars).
       bwind = 0.d0
*
* Common-envelope efficiency parameter (1.0).
       alpha1 = 1.d0
*
* Binding energy factor for common-envelope evolution (0.5).
       lambda = 0.5d0
*
* Dispersion in the Maxwellian for the SN kick speed (190 km/s).
       if (iflagns.eq.0) then
          sigma = 0.d0
       else
          sigma = 190.d0
       endif
       print*,'initpar: iflagns,sigma',iflagns,sigma
*

* Wind velocity parameter (proportional to V_wind**2; 0.125).
       beta = 0.125d0
*
* Wind accretion efficiency factor for momentum transfer (1.0).
       xi = 1.d0
*
* Bondi-Hoyle accretion factor (1.5).
       acc2 = 1.5d0
*
* Fraction of matter retained in a nova eruption (0.001).
       epsnov = 0.001d0
*
* Eddington limit factor for mass transfer (1.0; set large to allow super-Edd).
       eddfac = 1.d0
*
* Parameter to determine angular momentum change during RLOF mass-loss 
* (>0: lost material carries with it a fraction gamma of orbital angular momentum;
*  -1: material carries with it specific angular momentum of the primary; 
*  -2: material is lost from system as if a wind from the secondary).
       gamma = -1.d0
*
* Parameters to determine the timestep based on the evolution phase, 
* i.e. dt = pts*t_phase (pts1=MS,HG,HeMS; pts2=GB,CHeB,HeGB; pts3=AGB).
       pts1 = 0.05d0
       pts2 = 0.01d0
       pts3 = 0.02d0
*
* Maximum allowed fractional change in mass per timestep.
       dmmax = 0.05d0
*
* Maximum allowed fractional change in radius per timestep.
       drmax = 0.1d0
*
       aursun = 214.95d0
       yeardy = 365.24d0
*
* Set the collision matrix (ktype). 
*
       CALL instar
*
* Initialize the stellar index array to negative value. 
*
       do i = 0,nmax
          CALL setType(i,-1)
       enddo
*
       RETURN
       END
***
       SUBROUTINE printStar(idstar)
       integer idstar,kw,iend
       real*8 mass,rad
       character*20 label
*
       CALL getType(idstar,kw)
       CALL getLabel(kw,label,iend)
       CALL getMass(idstar,mass)
       mass = MIN(mass,999.d0)
       CALL getRadius(idstar,rad)
       rad = MIN(rad,9999.d0)
       WRITE(6,500)label(1:iend),mass,rad
  500  FORMAT(a,',  mass = ',f12.8,' radius = ',f13.8)
*
       RETURN
       END
***
       SUBROUTINE printBinary(idbin)
       implicit none
       integer idbin,idstar,ip,is
       integer k,kw,iend1,iend2,iend3
       real*8 time,epoch1,epoch2
       real*8 m1,r1,m2,r2,a,ecc,v(3),vd
       character*20 label1,label2,label3
*
       CALL getIprim(idbin,idstar)
       ip = idstar
       CALL getType(idstar,kw)
       CALL getLabel(kw,label1,iend1)
       CALL getMass(idstar,m1)
       m1 = MIN(m1,999.d0)
       CALL getRadius(idstar,r1)
       r1 = MIN(r1,9999.d0)
*
       CALL getAge(idstar,time)
       CALL getVel(idstar,v)
       call getepoch(idstar,epoch1)	
       vd = 0.d0
       do k = 1,3
          vd = vd + v(k)**2
       enddo
       if(vd.gt.0.d0) vd = SQRT(vd)
*
       CALL getIsec(idbin,idstar)
       is = idstar
       CALL getType(idstar,kw)
       CALL getLabel(kw,label2,iend2)
       CALL getMass(idstar,m2)
       m2 = MIN(m2,999.d0)
       CALL getRadius(idstar,r2)
       call getepoch(idstar,epoch2)
       r2 = MIN(r2,9999.d0)
*
       CALL getBtype(idbin,kw)
       CALL getLabelb(kw,label3,iend3)
       CALL getSemi(idbin,a)
       a = MIN(a,9.999d+07)
       CALL getEcc(idbin,ecc)
       ecc = MIN(ecc,99.d0)
       vd = MIN(vd,99999.d0)
*
       WRITE(6,600)time
       WRITE(6,601)label3(1:iend3),label1(1:iend1),label2(1:iend2)
cThis is arranged so that the component of lowest id is printed first
       if (ip.lt.is) then
          WRITE(6,602)m1,r1
          WRITE(6,603)m2,r2
       elseif (ip.gt.is) then
          WRITE(6,602)m2,r2
          WRITE(6,603)m1,r1
       else
          write (6,*) 'ip = is 1, stopping'
          stop
       endif
       WRITE(6,604)a,ecc,vd
       print*,'epochs ',epoch1,epoch2
*
  600  FORMAT(' status at time = ',f12.4)
  601  FORMAT(a,' (',a,', ',a,')')
  602  FORMAT(' M = ',f12.8,' R = ',f13.8)
  603  FORMAT(' m = ',f12.8,' r = ',f13.8)
  604  FORMAT(' a = ',f14.4,' e = ',f7.3,' v = ',f8.1)
*
       RETURN
       END
***
       SUBROUTINE initStar(in,idstar,m_init)
       implicit none
       integer in,idstar,k,kw
       real*8 m_init,x(3),v(3),zini
       common /zset/ zini
*
       if(in.eq.1)then
          CALL initPar
          in = 2
       endif
*
       CALL getType(idstar,kw)
       if(kw.ge.0)then
          WRITE(6,*)' WARNING: id already in use ',idstar
       endif
*
       CALL setType(idstar,1)
       CALL setAge(idstar,0.d0)
       CALL setEpoch(idstar,0.d0)
       CALL setMStime(idstar,1.0d+10)
       CALL setTstep(idstar,0.d0)
       CALL setMass0(idstar,m_init)
       CALL setMass(idstar,m_init)
       CALL setMc(idstar,0.d0)
       CALL setMenv(idstar,0.d0)
       CALL setRadius(idstar,m_init)
       CALL setRc(idstar,0.d0)
       CALL setRenv(idstar,0.d0)
       CALL setLum(idstar,m_init)
       CALL setSpin(idstar,0.d0)
       CALL setRl(idstar,1.0d+10)
*
* Assume solar abundance for now. 
* (Note that chemical evolution is not yet accounted for).
*
       CALL setZ0(idstar,zini)
       CALL setY0(idstar,0.28d0)
       CALL setZ(idstar,zini)
       CALL setY(idstar,0.28d0)
*
       if(idstar.lt.10) print*,'idstar,zini =',idstar,zini
*                            
c       CALL setZ0(idstar,0.002d0)
c       CALL setY0(idstar,0.28d0)
c       CALL setZ(idstar,0.002d0)
c       CALL setY(idstar,0.28d0)
*
       do k = 1,3
          x(k) = 0.d0
          v(k) = 0.d0
       enddo
       CALL setPos(idstar,x)
       CALL setVel(idstar,v)
*
* Obtain initial values for stellar parameters. 
*
       CALL evStar(idstar,0.d0)
*
       RETURN
       END
***
       SUBROUTINE initBinary(in,idbin,a,e,m1,m2)
       implicit none
       integer in,idbin,idstar
       real*8 a,e,tb0,m1,m2,mx,aursun,yeardy
       COMMON /PARAMS/ aursun,yeardy
*
       if(m2.gt.m1)then
          mx = m1
          m1 = m2
          m2 = mx
       endif
*
       idstar = idbin
       CALL initStar(in,idstar,m1)
       CALL setIprim(idbin,idstar)
       idstar = idbin + 1
       CALL initStar(in,idstar,m2)
       CALL setIsec(idbin,idstar)
*
       CALL setSemi(idbin,a)
       CALL setEcc(idbin,e)
       tb0 = (a/aursun)*SQRT(a/(aursun*(m1+m2)))*yeardy
       CALL setTb(idbin,tb0)
       CALL setBtype(idbin,1)
*
       CALL evBinary(idbin,0.d0)
*
       RETURN
       END
***
       SUBROUTINE getLabel(kw,labelx,i2)
       implicit none
       include "interface_bse.h"
       integer kw,i1,i2
       character*20 labelx
       labelx = label(kw)
       CALL strip(labelx,i1,i2)
*     i2 = i2 + 3
*     if(kw.ge.10) i2 = i2 + 1
       RETURN
       END
***
       SUBROUTINE getLabelb(kw,labelx,i2)
       implicit none
       include "interface_bse.h"
       integer kw,i1,i2
       character*20 labelx
       labelx = labelb(kw)
       CALL strip(labelx,i1,i2)
       RETURN
       END
***
       SUBROUTINE strip(word,i1,i2)
       implicit none
       integer i,i1,i2,n
       parameter(n=80)
       character*(n) word
*
       do i = 1,n
          i1 = i
          if(word(i:i).ne.' ') goto 1
       enddo
  1    do i = i1,n
          i2 = i
          if(word(i:i).eq.' ') goto 2
       enddo
  2    i2 = i2 - 1
*
       RETURN
       END
***
       SUBROUTINE ssupdatetime(idstar,time)
       implicit none
       integer idstar
       real*8 time,age,dt
       CALL getAge(idstar,age)
       CALL getTstep(idstar,dt)
       time = age + dt
       RETURN
       END
***
       SUBROUTINE bsupdatetime(idbin,time)
       implicit none
       integer idbin,idstar
       real*8 time
       CALL getIprim(idbin,idstar)
       CALL ssupdatetime(idstar,time)
       RETURN
       END
***
       SUBROUTINE binaryexists(idbin,iexist)
       implicit none
       integer idbin,iexist,kw
       CALL getBtype(idbin,kw)
       iexist = 1
       if(kw.gt.3) iexist = 0
       RETURN
       END
***
       SUBROUTINE setAge(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       age(idstar) = time
       RETURN
       END
***
       SUBROUTINE getAge(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = age(idstar)
       RETURN
       END
***
       SUBROUTINE setMStime(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       ms_lifetime(idstar) = time
       RETURN
       END
***
       SUBROUTINE getMStime(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = ms_lifetime(idstar)
       RETURN
       END
***
       SUBROUTINE setTstep(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       standard_timestep(idstar) = time
       RETURN
       END
***
       SUBROUTINE getTstep(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = standard_timestep(idstar)
       RETURN
       END
***
       SUBROUTINE setEpoch(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       epoch(idstar) = time
       RETURN
       END
***
       SUBROUTINE getEpoch(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = epoch(idstar)
       RETURN
       END
***
       SUBROUTINE setMass(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       zmass(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMass(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = zmass(idstar)
       RETURN
       END
***
       SUBROUTINE setMass0(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       zmass0(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMass0(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = zmass0(idstar)
       RETURN
       END
***
       SUBROUTINE setRadius(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       radius(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRadius(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = radius(idstar)
       RETURN
       END
***
       SUBROUTINE setLum(idstar,lum)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 lum
       zlum(idstar) = lum
       RETURN
       END
***
       SUBROUTINE getLum(idstar,lum)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 lum
       lum = zlum(idstar)
       RETURN
       END
***
       SUBROUTINE setY0(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       yinit(idstar) = y
       RETURN
       END
***
       SUBROUTINE getY0(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       y = yinit(idstar)
       RETURN
       END
***
       SUBROUTINE setZ0(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       zinit(idstar) = z
       RETURN
       END
***
       SUBROUTINE getZ0(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       z = zinit(idstar)
       RETURN
       END
***
       SUBROUTINE setY(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       ycurr(idstar) = y
       RETURN
       END
***
       SUBROUTINE getY(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       y = ycurr(idstar)
       RETURN
       END
***
       SUBROUTINE setZ(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       zcurr(idstar) = z
       RETURN
       END
***
       SUBROUTINE getZ(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       z = zcurr(idstar)
       RETURN
       END
***
       SUBROUTINE setMc(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       zmassc(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMc(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = zmassc(idstar)
       RETURN
       END
***
       SUBROUTINE setRc(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       radc(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRc(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = radc(idstar)
       RETURN
       END
***
       SUBROUTINE setMenv(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       menv(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMenv(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = menv(idstar)
       RETURN
       END
***
       SUBROUTINE setRenv(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       renv(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRenv(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = renv(idstar)
       RETURN
       END
***
       SUBROUTINE setSpin(idstar,sp)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 sp
       spin(idstar) = sp
       RETURN
       END
***
       SUBROUTINE getSpin(idstar,sp)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 sp
       sp = spin(idstar)
       RETURN
       END
***
       SUBROUTINE setRl(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rlof(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRl(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = rlof(idstar)
       RETURN
       END
***
       SUBROUTINE setPos(idstar,x)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 x(3)
       do k = 1,3
          xstar(k,idstar) = x(k)
       enddo
       RETURN
       END
***
       SUBROUTINE getPos(idstar,x)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 x(3)
       do k = 1,3
          x(k) = xstar(k,idstar)
       enddo
       RETURN
       END
***
       SUBROUTINE setVel(idstar,v)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 v(3)
       do k = 1,3
          vstar(k,idstar) = v(k)
       enddo
       RETURN
       END
***
       SUBROUTINE getVel(idstar,v)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 v(3)
       do k = 1,3
          v(k) = vstar(k,idstar)
       enddo
       RETURN
       END
***
       SUBROUTINE setType(idstar,ktype1)
       implicit none
       include "interface_bse.h"
       integer idstar,ktype1
       kstar(idstar) = ktype1
       RETURN
       END
***
       SUBROUTINE getType(idstar,ktype1)
       implicit none
       include "interface_bse.h"
       integer idstar,ktype1
       ktype1 = kstar(idstar)
       RETURN
       END
***
       SUBROUTINE setBtype(idbin,ktype1)
       implicit none
       include "interface_bse.h"
       integer idbin,ktype1
       kstarb(idbin) = ktype1
       RETURN
       END
***
       SUBROUTINE getBtype(idbin,ktype1)
       implicit none
       include "interface_bse.h"
       integer idbin,ktype1
       ktype1 = kstarb(idbin)
       RETURN
       END
***
       SUBROUTINE setIprim(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       iprim(idbin) = idstar
       RETURN
       END
***
       SUBROUTINE getIprim(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       idstar = iprim(idbin)
       RETURN
       END
***
       SUBROUTINE setIsec(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       isec(idbin) = idstar
       RETURN
       END
***
       SUBROUTINE getIsec(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       idstar = isec(idbin)
       RETURN
       END
***
       SUBROUTINE getMprim(idbin,mass)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       real*8 mass
       idstar = iprim(idbin)
       mass = zmass(idstar)
       RETURN
       END
***
       SUBROUTINE getMsec(idbin,mass)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       real*8 mass
       idstar = isec(idbin)
       mass = zmass(idstar)
       RETURN
       END
***
       SUBROUTINE setSemi(idbin,arsun)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 arsun
       sep(idbin) = arsun
       RETURN
       END
***
       SUBROUTINE getSemi(idbin,arsun)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 arsun
       arsun = sep(idbin)
       RETURN
       END
***
       SUBROUTINE setEcc(idbin,e)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 e
       ecc(idbin) = e
       RETURN
       END
***
       SUBROUTINE getEcc(idbin,e)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 e
       e = ecc(idbin)
       RETURN
       END
***
       SUBROUTINE setTb(idbin,tbday)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 tbday
       tb(idbin) = tbday
       RETURN
       END
***
       SUBROUTINE getTb(idbin,tbday)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 tbday
       tbday = tb(idbin)
       RETURN
       END
***
       SUBROUTINE evStar(idstar,time)
       implicit none
       integer idstar,kw
       real*8 time,tphys,tphysf,dmmax,drmax
       COMMON /TSTEPC/ dmmax,drmax
       real*8 z,zpars(20)
       real*8 mass,mt,mc,me
       real*8 rad,rc,re
       real*8 epch,tm,lum,ospin,vs(3)
*
       CALL getAge(idstar,tphys)
       tphysf = time
       if(tphysf.ge.tphys.or.tphys.le.1.0d-10)then
*
          CALL getZ0(idstar,z)
          CALL zcnsts(z,zpars)
*
          CALL getType(idstar,kw)
          CALL getEpoch(idstar,epch)
          CALL getMStime(idstar,tm)
          CALL getMass0(idstar,mass)
          CALL getMass(idstar,mt)
          CALL getMc(idstar,mc)
          CALL getMenv(idstar,me)
          CALL getRadius(idstar,rad)
          CALL getRc(idstar,rc)
          CALL getRenv(idstar,re)
          CALL getLum(idstar,lum)
          CALL getSpin(idstar,ospin)
          CALL getVel(idstar,vs)
c          print*,rad, 'before evolvlb'
*
          CALL evolv1b(kw,mass,mt,rad,lum,mc,rc,me,re,ospin,
     &                epch,tm,tphys,tphysf,z,zpars,dmmax,drmax,vs)
*
c         if (kw.ne.13) then
c          print*,rad, 'after evolvlb'
          CALL setType(idstar,kw)
          CALL setAge(idstar,tphys)
          CALL setEpoch(idstar,epch)
          CALL setMStime(idstar,tm)
          CALL setTstep(idstar,tphysf-tphys)
          CALL setMass0(idstar,mass)
          CALL setMass(idstar,mt)
          CALL setMc(idstar,mc)
          CALL setMenv(idstar,me)
          CALL setRadius(idstar,rad)
          CALL setRc(idstar,rc)
          CALL setRenv(idstar,re)
          CALL setLum(idstar,lum)
          CALL setSpin(idstar,ospin)
          CALL setVel(idstar,vs)
c         else
          if (kw.eq.13) then
              write (6,*) 'interface ',idstar,kw,tphys,epch,tm,tphysf,
     &            mass,
     &           mt,mc,me,rad,rc,re,lum,ospin,vs
          endif
*
       endif
*
       RETURN
       END
***
       SUBROUTINE evBinary(idbin,time)
       implicit none
       integer idbin,idstar,id1,id2
       integer kw1,kw2,kw(2),k
       real*8 time,tphys,tphysf,dmmax,drmax,aursun,yeardy
       COMMON /TSTEPC/ dmmax,drmax
       COMMON /PARAMS/ aursun,yeardy
       real*8 z,zpars(20)
       real*8 mass(2),mt(2),mc(2),me(2)
       real*8 rad(2),rc(2),re(2),rol(2)
       real*8 epch(2),tm(2),lum(2),ospin(2),dmdt(2)
       real*8 semi,tb,ecc,vs(3)
       logical iprint
*
       CALL getIprim(idbin,idstar)
       CALL getAge(idstar,tphys)
       tphysf = time
       iprint = .false.
c      iprint = .true.
*
       if(tphysf.ge.tphys.or.tphys.le.1.0d-10)then
*
          CALL getZ0(idstar,z)
          CALL zcnsts(z,zpars)
          CALL getVel(idstar,vs)
*
          id1 = idstar
          do k = 1,2
             if(k.eq.2) CALL getIsec(idbin,idstar)
             CALL getType(idstar,kw(k))
             CALL getEpoch(idstar,epch(k))
             CALL getMStime(idstar,tm(k))
             CALL getMass0(idstar,mass(k))
             CALL getMass(idstar,mt(k))
             CALL getMc(idstar,mc(k))
             CALL getMenv(idstar,me(k))
             CALL getRadius(idstar,rad(k))
             CALL getRc(idstar,rc(k))
             CALL getRenv(idstar,re(k))
             CALL getLum(idstar,lum(k))
             CALL getSpin(idstar,ospin(k))
             CALL getRl(idstar,rol(k))
          enddo
          id2 = idstar
          kw1 = kw(1)
          kw2 = kw(2)
*
          CALL getSemi(idbin,semi)
          CALL getEcc(idbin,ecc)
          tb = (semi/aursun)*SQRT(semi/(aursun*(mt(1)+mt(2))))*yeardy
*
          CALL evolv2b(kw,mass,mt,rad,lum,mc,rc,me,re,ospin,rol,
     &               dmdt,epch,tm,tphys,tphysf,z,zpars,dmmax,drmax,
     &                tb,ecc,vs)
*
c          print*,'evbin period',tb
          if(rad(2)/rol(2).gt.rad(1)/rol(1))then
             CALL setIprim(idbin,id2)
             CALL setIsec(idbin,id1)
             iprint = .true.
          elseif(kw1.ne.kw(1).or.kw2.ne.kw(2))then
             iprint = .true.
          endif
*
          idstar = id1
          do k = 1,2
             if(k.eq.2) idstar = id2
             CALL setType(idstar,kw(k))
             CALL setAge(idstar,tphys)
             CALL setEpoch(idstar,epch(k))
             CALL setMStime(idstar,tm(k))
             CALL setTstep(idstar,tphysf-tphys)
             CALL setMass0(idstar,mass(k))
             CALL setMass(idstar,mt(k))
             CALL setMc(idstar,mc(k))
             CALL setMenv(idstar,me(k))
             CALL setRadius(idstar,rad(k))
             CALL setRc(idstar,rc(k))
             CALL setRenv(idstar,re(k))
             CALL setLum(idstar,lum(k))
             CALL setSpin(idstar,ospin(k))
             CALL setRl(idstar,rol(k))
             CALL setVel(idstar,vs)
          enddo
*
* Set binary type and semi-major axis. 
*
          CALL getBtype(idbin,id1)
          if(tb.le.0.d0)then
             if(kw(1).lt.15.and.kw(2).lt.15)then
                id2 = 6
             elseif(kw(1).lt.15.or.kw(2).lt.15)then
                id2 = 4
             else
                id2 = 5
             endif
             semi = -1.d0
          else
             if(rad(1)/rol(1).lt.1.d0.and.rad(2)/rol(2).lt.1.d0)then
                id2 = 1
             else
                if(rad(1)/rol(1).ge.1.d0.and.rad(2)/rol(2).ge.1.d0)then
                   id2 = 3
                else
                   id2 = 2
                endif
             endif
             semi = aursun*((mt(1)+mt(2))*(tb/yeardy)**2)**(1.d0/3.d0)
          endif
          CALL setBtype(idbin,id2)
          CALL setSemi(idbin,semi)
          CALL setEcc(idbin,ecc)
          CALL setTb(idbin,tb)
          if(id2.ne.id1) iprint = .true.
*
          if(iprint)then
             CALL printBinary(idbin)
             call print_roche_data(idbin)
          endif
*
       endif
*
       RETURN
       END
*
************************************************************************
cAdditions for v3.1:
cAlso, an extra call to print_roche_data was also added to evBinary
       subroutine init_binaries (in, id, aRsun, e, m1, m2)
       implicit none
       integer in,id
       double precision aRsun,e,m1,m2,a
       character*200 datadir
       common /AMUSE/ datadir

C        if (in.ne.1) then
C           write (6,*) 'in != 1 in init_binaries not implemented,',
C      &        ' stopping'
C           stop
C        endif
       a = aRsun
       call initBinary(in,id,a,e,m1,m2)
       open (7,file=trim(datadir)//'/BSE.data')
       write (7,*) '      TIME    M1     M2  KW1 KW2    SEP    ECC R1/',
     &     'ROL1 R2/ROL2  TYPE'
       return
       end

       subroutine out_binary(id)
       implicit none
       integer id
       call printBinary(id)
       call print_roche_data(id)
       return
       end

       subroutine get_bs_updatetime(id, updt)
       implicit none
       double precision updt
       integer id
       call bsupdatetime(id,updt)
       return
       end

       subroutine morph_binary(idb, id3, a_factor, e, outcome)
       implicit none
       integer idb,id3,outcome
       double precision a_factor,e
       if (outcome.eq.0) then
          a_factor = 1
          call getEcc(idb,e)
          return
       endif
       if (outcome.eq.1) then
          call exchange_hiid(idb,id3,a_factor,e)
          return
       endif
       if (outcome.eq.2) then
          call exchange_loid(idb,id3,a_factor,e)
          return
       endif
       write (6,*) 'outcome ',outcome,
     &     ' unexpected in morph_binary, stopping'
       stop
       end

       subroutine out_star(ids)
       implicit none
       integer ids,kw,iend2
       character*20 label2
       call getType(ids,kw)
       CALL getLabel(kw,label2,iend2)
       write (6,*) 'type (',label2,'), '
       return
       end

       subroutine ev_binary(id,tMyr)
       implicit none
       integer id
       real*8 tMyr,updatetime
c 10    continue
c       call bsupdatetime(id,updatetime)
c       if (updatetime.ge.tMyr) then
          call evBinary(id,tMyr)
c          print*,'binary',id,tMyr
c       else
c          call evBinary(id,updatetime)
c          print*,'binary',id,updatetime
c          goto 10
c       endif
       return
       end

       subroutine binary_exists(id, iexist)
       implicit none
       integer id,iexist
       call binaryexists(id,iexist)
       return
       end

       subroutine get_sma(id, aRsun)
       implicit none
       integer id
       real*8 aRsun
       call getSemi(id,aRsun)
       return
       end

       subroutine get_ecc(id, e)
       implicit none
       integer id
       real*8 e
       call getEcc(id,e)
       return
       end

       subroutine get_loid_mass(id, m1)
       implicit none
       integer id,ip,is
       real*8 m1
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.lt.is) then
          call getMprim(id,m1)
       elseif (ip.gt.is) then
          call getMsec(id,m1)
       else
          write (6,*) 'ip = is 2, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine get_hiid_mass(id, m2)
       implicit none
       integer id,ip,is
       real*8 m2
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.gt.is) then
          call getMprim(id,m2)
       elseif (ip.lt.is) then
          call getMsec(id,m2)
       else
          write (6,*) 'ip = is 3, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine out_scatter(id,id3,tMyr,write_text)
       implicit none
       integer id,id3
       real*8 tMyr
       logical write_text
       return
       end

       subroutine init_stars(in, id3, m3)
       implicit none
       integer in,id3
       real*8 m3
       call initStar(in,id3,m3)
       return
       end

       subroutine ev_star(id3, tmyr)
       implicit none
       integer id3
       real*8 tmyr
       call evStar(id3,tmyr)
       return
       end

       subroutine get_mass(id3, m3)
       implicit none
       integer id3
       real*8 m3
       call getMass(id3,m3)
       return
       end

       subroutine get_radius(id3, r3)
       implicit none
       integer id3
       real*8 r3
       call getRadius(id3,r3)
       return
       end

       subroutine exchange_loid(idb,id3,a_factor,e)
       implicit none
       integer idb,id3,ip,is
       real*8 a_factor,e
       CALL getIprim(idb,ip)
       CALL getIsec(idb,is)
       if (ip.lt.is) then
          call setIprim(idb,id3)
       elseif (ip.gt.is) then
          call setIsec(idb,id3)
       else
          write (6,*) 'ip = is 4, stopping...'
          stop
       endif
       call update_binary_parameters(idb,a_factor,e)
       return
       end

       subroutine exchange_hiid(idb,id3,a_factor,e)
       implicit none
       integer idb,id3,ip,is
       real*8 a_factor,e
       CALL getIprim(idb,ip)
       CALL getIsec(idb,is)
       if (ip.lt.is) then
          call setIsec(idb,id3)
       elseif (ip.gt.is) then
          call setIprim(idb,id3)
       else
          write (6,*) 'ip = is 5, stopping...'
          stop
       endif
       call update_binary_parameters(idb,a_factor,e)
       return
       end

       subroutine update_binary_parameters(idb,a_factor,e)
       implicit none
       integer idb
       real*8 a_factor,e,arsun
       call getSemi(idb,arsun)
       arsun = arsun*a_factor
       call setSemi(idb,arsun)
       call setEcc(idb,e)
       return
       end

       subroutine print_roche_data(idbin)
       implicit none
       integer idbin,idstar
       integer k,kw,iend1,iend2,iend3
       real*8 time
       real*8 m1,r1,m2,r2,a,ecc,v(3),vd
       real*8 rl1,rl2
       integer kw1,kw2,id1,id2,id10,id20,kw10,kw20
       character*20 label1,label2,label3,type,label30
       logical first,swap
       data first/.true./
       save id10,id20,kw10,kw20,swap,label30
*
       CALL getIprim(idbin,idstar)
       id1 = idstar
       CALL getType(idstar,kw)
       kw1 = kw
       call getRl(idstar,rl1)
       CALL getLabel(kw,label1,iend1)
       CALL getMass(idstar,m1)
       m1 = MIN(m1,999.d0)
       CALL getRadius(idstar,r1)
       r1 = MIN(r1,9999.d0)
*
       CALL getAge(idstar,time)
       CALL getVel(idstar,v)
       vd = 0.d0
       do k = 1,3
          vd = vd + v(k)**2
       enddo
       if(vd.gt.0.d0) vd = SQRT(vd)
*
       CALL getIsec(idbin,idstar)
       id2 = idstar
       CALL getType(idstar,kw)
       kw2 = kw
       call getRl(idstar,rl2)
       CALL getLabel(kw,label2,iend2)
       CALL getMass(idstar,m2)
       m2 = MIN(m2,999.d0)
       CALL getRadius(idstar,r2)
       r2 = MIN(r2,9999.d0)
*
       CALL getBtype(idbin,kw)
       CALL getLabelb(kw,label3,iend3)
       CALL getSemi(idbin,a)
       a = MIN(a,9.999d+07)
       CALL getEcc(idbin,ecc)
       ecc = MIN(ecc,99.d0)
       vd = MIN(vd,99999.d0)
*
       if (first) then
          first = .false.
          type = 'INITIAL'
          swap = .false.
       else
          if (id1.eq.id10) then
             if (id2.eq.id20) then
cBoth components still in place; could be kw change or change of binary type 
c(e.g. no longer detached) or end of run
                if (kw1.ne.kw10.or.kw2.ne.kw20) then
                   type = 'KW CHNGE'
                else
                   if (label3.ne.label30) then
                      type = label3
                   else
                      type = 'MAX TIME'
                   endif
                endif
             else
cStar 2 must have been exchanged
                type = 'EXCHANGE'
             endif
          else
cStar 1 no longer in place
             if (id1.eq.id20) then
c..but star 1 is still in the binary; swop printing order
                swap = .not.swap
                if (id2.eq.id10) then
cStars 1 and 2 are interchanged
                   if (kw1.ne.kw20.or.kw2.ne.kw10) then
                      type = 'KW CHNGE'
                   else
                      if (label3.ne.label30) then
                         type = label3
                      else
cThis should not happen
                         type = 'ERROR 1'
                      endif
                   endif
                else
cStar 1 has swapped, but star 2 is no longer in the system
                   type = 'EXCHANGE'
                endif
             else
cStar 1 is no longer in the sytem...
                if (id2.eq.id20) then
c...but star 2 is still in place
                   type = 'EXCHANGE'
                else if (id2.eq.id10) then
c...but star 2 has been swapped
                   type = 'EXCHANGE'
                   swap = .not.swap
                else
cNeither star surviving
                   type = 'ERROR 2'
                endif
             endif
          endif
       endif
       if (.not.swap) then
          WRITE(7,100)time,m1,m2,kw1,kw2,a,ecc,r1/rl1,r2/rl2,type,label3
       else
          WRITE(7,100)time,m2,m1,kw2,kw1,a,ecc,r2/rl2,r1/rl1,type,label3
       endif
  100  format (f11.4,2f7.3,2i3,f10.3,f6.2,2f7.3,1x,a10,a20)
*
       id10 = id1
       id20 = id2
       kw10 = kw1
       kw20 = kw2
       label30 = label3
c
       RETURN
       END
c===========================================================================
cNew wrappers added Edinburgh 1/5/6
cThe following two corrected 13/11/6 following Mirek's email of 12/11/6

       subroutine get_loid(id,k)
       implicit none
       integer id,k
       integer ip,is
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.lt.is) then
          k = ip
       elseif (ip.gt.is) then
          k = is
       else
          write (6,*) 'ip = is 2, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine get_hiid(id,k)
       implicit none
       integer id,k
       integer ip,is
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.lt.is) then
          k = is
       elseif (ip.gt.is) then
          k = ip
       else
          write (6,*) 'ip = is 3, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine get_ss_type(i,k)
       implicit none
       integer i,k
       call getType(i,k)
       return
       end

       subroutine get_bs_type(i,k)
       implicit none
       integer i,k
       call getBtype(i,k)
       return
       end

       subroutine get_ss_updatetime(i,t)
       implicit none
       integer i
       double precision t
       call ssupdatetime(i,t)
       return
       end
c===========================================================================
cNew subroutine added Edinburgh 4/5/6
cModelled on initBinary
       subroutine create_binary(idbin,id1,id2,sma,e,time)
       implicit none
c      include "interface_bse.h"
       double precision mx,m1,m2,sma,t1,t2,tb0,e
       integer idbin,idstar,id1,id2,idx,ix1,ix2
       integer kw1,kw2,kw(2),k
       real*8 time,tphys,tphysf,dmmax,drmax,aursun,yeardy
       COMMON /TSTEPC/ dmmax,drmax
       COMMON /PARAMS/ aursun,yeardy
       real*8 z,zpars(20)
       real*8 mass(2),mt(2),mc(2),me(2)
       real*8 rad(2),rc(2),re(2),rol(2)
       real*8 epch(2),tm(2),lum(2),ospin(2),dmdt(2)
       real*8 semi,tb,ecc,vs(3)
       logical iprint
c      
       call getAge(id1,t1)
       call getAge(id2,t2)
cPut the star of higher mass first
       call getMass(id1,m1)
       call getMass(id2,m2)
cCheck that both stars have been evolved to the same time
       if (t1.ne.t2) then
          print*,'create_binary: stars have different ages',id1,id2,  
     &            idbin,t1,t2,m1,m2,sma,e,time
c          stop
       endif
                                
       if(m2.gt.m1)then
          mx = m1
          m1 = m2
          m2 = mx
          idx = id1
          id1 = id2
          id2 = idx
       endif
*
       print*,'id1,id2,m1,m2,idbin =',id1,id2,m1,m2,idbin
       CALL setIprim(idbin,id1)
       CALL setIsec(idbin,id2)
       CALL setSemi(idbin,sma)
       CALL setEcc(idbin,e)
       tb0 = (sma/aursun)*SQRT(sma/(aursun*(m1+m2)))*yeardy
       CALL setTb(idbin,tb0)
cIt would be nice to distinguish primordial from three-body binaries, but still....
       CALL setBtype(idbin,1)
cNow we have to set all the parameters that would be needed in a call to evBinary
       tphys = t1
       tphysf = t1
       idstar = id1
       CALL getZ0(idstar,z)
       CALL zcnsts(z,zpars)
       CALL getVel(idstar,vs)
*
       iprint = .false.
       do k = 1,2
          if(k.eq.2) CALL getIsec(idbin,idstar)
          CALL getType(idstar,kw(k))
          CALL getEpoch(idstar,epch(k))
          CALL getMStime(idstar,tm(k))
          CALL getMass0(idstar,mass(k))
          CALL getMass(idstar,mt(k))
          CALL getMc(idstar,mc(k))
          CALL getMenv(idstar,me(k))
          CALL getRadius(idstar,rad(k))
          CALL getRc(idstar,rc(k))
          CALL getRenv(idstar,re(k))
          CALL getLum(idstar,lum(k))
          CALL getSpin(idstar,ospin(k))
          CALL getRl(idstar,rol(k))
       enddo
       kw1 = kw(1)
       kw2 = kw(2)
*
       CALL getSemi(idbin,semi)
       CALL getEcc(idbin,e)
       tb = (semi/aursun)*SQRT(semi/(aursun*(mt(1)+mt(2))))*yeardy
*
       CALL evolv2b(kw,mass,mt,rad,lum,mc,rc,me,re,ospin,rol,
     &    dmdt,epch,tm,tphys,tphysf,z,zpars,dmmax,drmax,
     &     tb,e,vs)
*
       if(rad(2)/rol(2).gt.rad(1)/rol(1))then
          CALL setIprim(idbin,id2)
          CALL setIsec(idbin,id1)
          iprint = .true.
       elseif(kw1.ne.kw(1).or.kw2.ne.kw(2))then
          iprint = .true.
       endif
*
       idstar = id1
       do k = 1,2
          if(k.eq.2) idstar = id2
          CALL setType(idstar,kw(k))
          CALL setAge(idstar,tphys)
          CALL setEpoch(idstar,epch(k))
          CALL setMStime(idstar,tm(k))
          CALL setTstep(idstar,tphysf-tphys)
          CALL setMass0(idstar,mass(k))
          CALL setMass(idstar,mt(k))
          CALL setMc(idstar,mc(k))
          CALL setMenv(idstar,me(k))
          CALL setRadius(idstar,rad(k))
          CALL setRc(idstar,rc(k))
          CALL setRenv(idstar,re(k))
          CALL setLum(idstar,lum(k))
          CALL setSpin(idstar,ospin(k))
          CALL setRl(idstar,rol(k))
          CALL setVel(idstar,vs)
       enddo
*
* Set binary type and semi-major axis. 
*
       CALL getBtype(idbin,ix1)
       if(tb.le.0.d0)then
          if(kw(1).lt.15.and.kw(2).lt.15)then
             ix2 = 6
          elseif(kw(1).lt.15.or.kw(2).lt.15)then
             ix2 = 4
          else
             ix2 = 5
          endif
          semi = -1.d0
       else
          if(rad(1)/rol(1).lt.1.d0.and.rad(2)/rol(2).lt.1.d0)then
             ix2 = 1
          else
             if(rad(1)/rol(1).ge.1.d0.and.rad(2)/rol(2).ge.1.d0)then
                ix2 = 3
             else
                ix2 = 2
             endif
          endif
          semi = aursun*((mt(1)+mt(2))*(tb/yeardy)**2)**(1.d0/3.d0)
       endif
       CALL setBtype(idbin,ix2)
       CALL setSemi(idbin,semi)
       CALL setEcc(idbin,e)
       CALL setTb(idbin,tb)
       if(id2.ne.id1) iprint = .true.
*     
       if(iprint)then
          CALL printBinary(idbin)
          call print_roche_data(idbin)
       endif
*
       return
       end
*
*
      SUBROUTINE collStars(in,id1,id2,time)
      implicit none
*
* Developed by J. Hurley (May 1, 2006) for use 
* with the SSE/BSE McScatter interface. 
*
* Determines collision product of Star 1 and Star 2 
* based on Hurley et al. 2002, MNRAS, 329, 897. 
* Assumes that product is to be placed in position 
* occupied by Star 1. 
* Any cases that potentially lead to common-envelope 
* (subgiant/giant + any other star) are treated here 
* as mergers without mass-loss. If instead it is 
* determined that the stars form a close binary prior 
* to contact, then the binary can be evolved with the 
* routine evBinary which contains the common-envelope 
* treatment.  
*
      include "const_bse.h"
*
      integer oldkw
      integer in,id1,id2,idstar
      integer kw1,kw2,kwstar
      real*8 time
      real*8 m1,m2,mstar
      real*8 m01,m02,m0star
      real*8 mc1,mc2,mcstar
      real*8 aj1,aj2,epch,ajstar,f1,f2
      real*8 tms1,tms2,tms
      real*8 xs(3),vs(3)
      real*8 z,zpars(20),tscls(20),lums(10),gb(10),tn
      logical tzo
*
cSee email from jhurley 15/x/8
c      tzo = .true.
      tzo = .false.
*
*     tn is not deffined so is set to 0
*
      tn = 0.d0      
*
      CALL getType(id1,kw1)
      CALL getType(id2,kw2)
      print*,'id1,id2,kw1,kw2,time = ',id1,id2,kw1,kw2,time
      if(kw1.lt.0) kw1 = abs(kw1)
      if(kw2.lt.0) kw2 = abs(kw2)
      kwstar = ktype(kw1,kw2)
      print*,'id1,id2,kw1,kw2,kw,time = ',id1,id2,kw1,kw2,kwstar,time
      if(kwstar.ge.113) kwstar = kwstar - 100
*
      CALL getMass(id1,m1)
      CALL getMass(id2,m2)
      mstar = m1 + m2
      CALL getMass0(id1,m01)
      CALL getMass0(id2,m02)
      m0star = mstar
*
      CALL getMc(id1,mc1)
      CALL getMc(id2,mc2)
      mcstar = mc1 + mc2
*
      CALL getEpoch(id1,epch)
      aj1 = time - epch
      CALL getEpoch(id2,epch)
      aj2 = time - epch
      ajstar = 0.d0
      CALL getMStime(id1,tms1)
      CALL getMStime(id2,tms2)
      print*,'kw1,kw2,kwstar,m1,m2,m01,m02,mc1,mc2,aj1,aj2,time= ',
     &        kw1,kw2,kwstar,m1,m2,m01,m02,mc1,mc2,aj1,aj2,time
      print*,'mstar,m0star,mcstar = ',mstar,m0star,mcstar
*
      CALL getPos(id1,xs)
      CALL getVel(id1,vs)
*
      idstar = id1
      CALL setType(idstar,-1)
      print*,' before initS '
      call flush(6)
      CALL initStar(in,idstar,mstar)
      print*,' after initS '
      call flush(6)
      CALL getZ0(idstar,z)
      print*,'z of merger',z
      CALL zcnsts(z,zpars)
      CALL getMStime(idstar,tms)
      print*,'z, tms',z,tms
*
      if(kwstar.eq.1)then
         if(kw1.eq.7) kwstar = 7
         ajstar = 0.1d0*tms*(aj1*m1/tms1+aj2*m2/tms2)/mstar
      elseif(kwstar.eq.4)then
         if(kw1.le.1)then
            mcstar = m2
            ajstar = aj2/tms2
         else
            mcstar = m1
            ajstar = aj1/tms1
         endif
         print*,'a- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
         CALL gntage(mcstar,mstar,kwstar,zpars,m0star,ajstar)
         print*,'b- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
      elseif(kwstar.eq.7)then
         if(kw1.lt.10)then
            ajstar = (tms/mstar)*(aj1*m1/tms1)
         else
            ajstar = (tms/mstar)*(aj2*m2/tms2)
         endif
      elseif(kwstar.le.9)then
         print*,'c- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
         CALL gntage(mcstar,mstar,kwstar,zpars,m0star,ajstar)
         print*,'d- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
      elseif(kwstar.le.12)then
         if(kwstar.lt.12.and.mstar.ge.1.44)then
            mstar = 0.d0
            kwstar = 15
         endif
      elseif(kwstar.eq.13.or.kwstar.eq.14)then
*       Absorb all mass into NS/BH unless taking the 
*       unstable Thorne-Zytkow object option. 
         if(tzo)then
            mstar = 0.d0
            if(kw1.ge.10) mstar = mstar + m1
            if(kw2.ge.10) mstar = mstar + m2
            m0star = mstar
            mcstar = mstar
            print*,'TZO kw1,kw2,mstar,m1,m2=',kw1,kw2,mstar,m1,m2
         endif
         if(kwstar.eq.13.and.mstar.gt.mxns) kwstar = 14
      elseif(kwstar.eq.15)then
         mstar = 0.d0
      elseif(kwstar.gt.100)then
*       Common envelope cases.  
*       In the absence of CE treatment assume that 
*       merger proceeds with no mass lost from the system. 
         kwstar = kwstar - 100
         if(kwstar.eq.4.or.kwstar.eq.7)then
            f1 = mc1
            if(kw1.le.3.or.kw1.eq.10) f1 = 0.d0
            if(kw1.eq.7) f1 = m1*aj1/tms1
            if(kw1.eq.4.or.kw1.eq.5)then
        print*,'e- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               CALL star(kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars)
        print*,'f- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               f1 = f1*(aj1 - tscls(2))/(tscls(13) - tscls(2))
            endif
            f2 = mc2
            if(kw2.le.3.or.kw2.eq.10) f2 = 0.d0
            if(kw2.eq.7) f2 = m2*aj2/tms2
            if(kw2.eq.4.or.kw2.eq.5)then
        print*,'g- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               CALL star(kw2,m02,m2,tms2,tn,tscls,lums,gb,zpars)
        print*,'h- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               f2 = f2*(aj2 - tscls(2))/(tscls(13) - tscls(2))
            endif
            ajstar = (f1+f2)/mcstar
         endif
         if(kwstar.eq.7)then
            ajstar = ajstar*tms
cccc         elseif(kwstar.eq.7)then
         else
         print*,'i- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
            oldkw = kwstar
            CALL gntage(mcstar,mstar,kwstar,zpars,m0star,ajstar)
         print*,'j- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
cThe following two lines are really meant to deal with black hole outcomes
cThe previous version gave NaN stellar radii
            print*,'collision routine kw,ajstar',kwstar,ajstar,
     &        'changed to original value',oldkw,'and 0.d0'
c            kwstar = oldkw
c            ajstar = 0.d0
         endif
         print*,'kwstar,mstar,m0star,mcstar,ajstar = ',
     &           kwstar,mstar,m0star,mcstar,ajstar
      else
*       This should not be reached.
        print*,' *************** '
        kwstar = 1
      endif
*
      call getz0(idstar,z)
c      print*,z
      call getz(idstar,z)
c      print*,z
      epch = time - ajstar
      ajstar = MAX(0.d0,time-1.0d-14)
      CALL setAge(idstar,time)       
      CALL setEpoch(idstar,epch)
      CALL setType(idstar,kwstar)
      CALL setMass0(idstar,m0star)
      CALL setMass(idstar,mstar)
      CALL setMc(idstar,mcstar)
      CALL setPos(idstar,xs)
      CALL setVel(idstar,vs)
*
      call getz0(idstar,z)
      print*,z
      print*,'before evSTAR kwstar,time,epoch,mstar,m0star,mcstar =',
     &        kwstar,time,epch,mstar,m0star,mcstar
      CALL evStar(idstar,time)
      call getz0(idstar,z)
      print*,z
      call getMass(idstar,mstar)
      call getMass0(idstar,m0star)
      call getMc(idstar,mcstar)
      call getEpoch(idstar,epch)
      call getType(idstar,kwstar)
      print*,'after evStar kwstar,time,epoch,mstar,m0star,mcstar =',
     &        kwstar,time,epch,mstar,m0star,mcstar
           
      CALL printStar(idstar)
*
      RETURN
      END
*
***
*
* Developed by J. Hurley (August 29, 2006) for use
* with the SSE/BSE McScatter interface.
*
*   TURN  - an initial guess at the turn-off mass (set large if unsure)
*   TPHYS - current age in Myr
*   ZPARS - the metallicity parameter array set by an earlier call
*            to ZCNSTS.
*
***
      SUBROUTINE MTURN(TURN,TPHYS,ZPARS)
*
*
*       Current MS turn-off mass.
*       --------------------------------------
*
      IMPLICIT NONE
      INTEGER  I,II,IMAX
      PARAMETER(IMAX=30)
      REAL*8 TURN,TPHYS,ZPARS(20)
      REAL*8 TM,TURN2,DM,FMID,TACC
      PARAMETER(TACC=0.001D0)
      REAL*8 THOOKF,TBGBF
      EXTERNAL THOOKF,TBGBF
*
      TURN2 = 100.D0
      TM = MAX(ZPARS(8),THOOKF(TURN2))*TBGBF(TURN2)
      IF(TM.GT.TPHYS)THEN
         TURN = TURN2
         GOTO 40
      ENDIF
*
      II = 0
 25   TM = MAX(ZPARS(8),THOOKF(TURN))*TBGBF(TURN)
      IF(TM.GT.TPHYS)THEN
         IF(TPHYS.LE.0.D0.OR.TURN.GT.98.D0) GOTO 40 
         TURN = 2.D0*TURN
         II = II + 1
         GOTO 25
      ENDIF
      TURN2 = TURN
      DM = TURN
      DO 30 , I = 1,IMAX
         DM = 0.5D0*DM
         TURN = TURN2 - DM
         TM = MAX(ZPARS(8),THOOKF(TURN))*TBGBF(TURN)
         FMID = TM - TPHYS
         IF(FMID.LT.0.0) TURN2 = TURN
         IF(DM.LT.TACC.OR.ABS(FMID).LT.1.0D-14) GOTO 40
         IF(I.EQ.IMAX)THEN
            GOTO 40
         ENDIF
 30   CONTINUE
 40   CONTINUE
*
      RETURN
*
      END
***
        



