      subroutine input
*
*
*       parameter input.
*       ----------------
*
      include 'common.h'
*
      integer i,ixx1,ixx2,ixx3
*
      real*8 xxx4
*
*
*              open input file and read initial parameters
*
*
*      open(1,file='src/mont.run')
      open(1,file=trim(datadir)//'/static/mont.run')
*
      read(1,*) ixx1,ixx2,istart,ncor,nmin,ixx3,nzonc,nminzo,ntwo,
     &          imodel,iprint,ib3f,iexch,tcrit,tcomp,qe,alphal,
     &          alphah,brakem,body1,bodyn,fracb,amin,amax,qvir,rbar,
     &          xxx4,w0,bmin0,bmax,tau0,gamma,xtid,rplum,dttp,dtte,
     &          dtte0,tcrevo,xtau,ytau,ybmin,zini,ikroupa,iflagns,
     &          iflagbh,nitesc
*
      read(1,*) (flagr(i),i=1,nlagra)
*
      close(1)
      write(6,*)ixx1,ixx2,istart,ncor,nmin,ixx3,nzonc,nminzo,ntwo,
     &          imodel,iprint,ib3f,iexch,tcrit,tcomp,qe,alphal,
     &          alphah,brakem,body1,bodyn,fracb,amin,amax,qvir,rbar,
     &          xxx4,w0,bmin0,bmax,tau0,gamma,xtid,rplum,dttp,dtte,
     &          dtte0,tcrevo,xtau,ytau,ybmin,zini,ikroupa,iflagns,
     &          iflagbh,nitesc
      call flush(6)
*
      bmin = bmin0
      bmax = 2.d0*bmin0 
      if(istart.eq.1) then
        irun = ixx1
        nt = ixx2
        print*, "NT in input function is:"
        print*, nt
        nt0 = nt
        iseed = ixx1q
        nz0 = ixx3
        zmbar = xxx4
        nt00 = nt0
      endif
*
*
*     ------------------------------------------------------------------------
*
*                  INPUT PARAMETERS
*                  ----------------
*
*     irun    -  initial sequence of random numbers
*     nt      -  total number of objects (stars and binaries) at T=0
*                ns - number of single stars, nb - number of binaries
*                (nt = ns + nb), nss - number of stars (nss = nt + nb)
*     istart  -  1 - initial model,    .ne.1 - restart
*     ncor    -  number of stars to calculate the central parameters
*     nmin    -  minimum number of stars to calculate the central
*                parameters
*     nz0     -  number of stars in each zone at T=0
*     nzonc   -  minimum number of zones in the core
*     nminzo  -  minimum number of stars in a zone
*     ntwo    -  maximum index of 2
*     imodel  -  initial model: 1- uniform & isotropic, 2- Plummer,
*                3- King, 4 - M67
*     iprint  -  0- full diagnostic information, 1- diagnostic info.
*                suppressed
*     ib3f    -  1 - Spitzer's, 2 - Heggie's formula for three-body binary
*                interaction with field stars, 3 - use Pmax for interaction
*                probability  4 - three- and four-body numerical integration
*     iexch   -  0 - no exchange in any interactions, 1 - exchange only in 
*                binary field star interacions, 2 - exchange in all
*                interactions (binary - field and binary - binary)
*     tcrit   -  termination time in units of the crossing time
*     tcomp   -  maximum computing time in hours
*     qe      -  energy tolerance
*     alphal  -  power-law index for initial mass function for masses
*                smaller than breake mass: -1 - equal mass case
*     alphah  -  power-law index for initial mass function for masses
*                greater than breake mass. If alphal=alphah the IMF does
*                not have a break
*     brakem  -  the mass in which the IMF is broken. If brakem is smaller
*                than the minimum mass (bodyn) than the break mass is as for 
*                the Kroupa mass function (brakem = 0.5 Mo) 
*     body1   -  maximum particle mass before scaling (solar mass)
*     bodyn   -  minimum particle mass before scaling (solar mass)
*     fracb   -  primordial binary fraction by number. nb = fracb*nt,
*                ns = (1 - fracb)*nt, nss = (1 + fracb)*nt
*                fracb > 0 - primordial binaries
*                fracb = 0 - only dynamical binaries
*     amin    -  minimum semi-major axis of binaries (in sollar units)
*                = 0 then amin = 2*(R1+R2), > 0 then amin = amin
*     amax    -  maximum semi-major axis of binaries (in sollar units)
*     qvir    -  virial ratio  (qvir = 0.5 for equilibrium)
*     rbar    -  tidal radius in pc, halfmass radius in pc for isolated
*                cluster. No scaling - rbar = 1
*     zmbar   -  total mass of the cluster in sollar mass, 
*                no scaling zmbar = 1
*     w0      -  king model parameter
*     bmin    -  minimum value of sin(beta^2/2)
*     bmax    -  maximum value of sin(beta^2/2)
*     tau0    -  time step for a complite cluster model
*     gamma   -  parameter in the Coulomb logarithm (standard value = 0.11)
*     xtid    -  coeficient in the front of cluster tidal energy:
*                                                            -xtid*smt/rtid
*     rplum   -  for M67 rtid = rplum*rsplum (rsplum - scale radius for
*                plummer model)
*     dttp    -  time step (Myr) for profile output
*     dtte    -  time step (Myr) for mloss call for all objects
*                greater then tcrevo. For tphys less then tcrevo time step
*                is eqiual to dtte0
*     dtte0   -  time step (Myr) for mloss call for all objects for tphys
*                less then tcrevo. For tphys greater then tcrevo time step
*                is eqiual to dtte
*     tcrevo  -  critical time for which time step for mloss call changes from
*                dtte0 to dtte
*     xtau    -  call mloss for a particlular object when
*                (uptime(im1) - olduptime(im1))/tau/tscale < xtau
*     ytau    -  multiplication of tau0 (tau = ytau*tau0) after time
*                greater than tcrevo
*     ybmin   -  multiplication of bmin0 (bmin = ybmin*bmin0) after time
*                greater than tcrevo
*     zini    -  initial metalicity (solar z = 0.02, globular clusters
*                M4 - z = 0.002,  NGC6397 - z = 0.0002)
*     ikroupa -  0 - the initial binary parameters are picked up
*                according Kroupa's eigenevolution and feeding algorithm
*                (Kroupa 1995, MNRAS 277, 1507)
*                1 - the initial binary parameters are picked as for M67
*                model (Hurley et al. 2005)
*     iflagns -  0 - no SN natal kiks for NS formation, 1 - SN natal kicks
*                only for single NS formation, 2 - SN natal kick for single
*                NS formation and NS formation in binaries
*     iflagbh -  0 - no SN natal kiks for BH formation, 1 - SN natal kicks
*                only for single BH formation, 2 - SN natal kick for single
*                BH formation and BH formation in binaries
*     nitesc  -  0 - no iteration of the tidal radius and induced mass loss
*                due to stellar evolution, 1 - iteration of the tidal radius
*                and induced mass loss due to stellar evolution
*
*     ------------------------------------------------------------------------
*
*
      print*, "nt @ end input fct",nt

      return
*
      end
*
*
*
*
