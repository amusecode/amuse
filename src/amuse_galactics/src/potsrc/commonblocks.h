      include 'parameters.h'
      common /dfparameters/ npsi,nint,
     +	     tableE(10000),dfsersic(10000),dfnfw(10000),
     +     denspsibulge(10000),denspsihalo(10000)
      common /potconstants/ apot(20,0:nmax), fr(20,0:nmax), 
     +     dr, nr, lmax, potcor
      common /gparameters/  chalo, v0, a, 
     +                      nnn, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /newparameters/ drtrunchalo, cusp
      common /bulgepars/ stream
      common /legendre/ plcon(0:40)
      common /moreconstants/ v02, v03, rdisk2, diskconst
      common /energytable/ psic, psid
      common /flags/ idiskflag, ibulgeflag, ihaloflag, ibhflag
      common /sersicparameters/ Re,butt,Rho0,emax,ppp
      common /diskblackhole/ idiskbhflag
      parameter(nrmax=1000)
      common /diskpars/ sigr0, disksr, nrdisk
      common /splines/ rr(0:nrmax),fdrat(0:nrmax),drat2(0:nrmax),
     +              fszrat(0:nrmax), szrat2(0:nrmax), nrspl
      common /scalings/ rscale, vscale, mscale, Ro
      real psizh,psi00
      common /cutoffs/ fcut_halo, fcut_bulge

      real adens(20,0:nmax),s1(0:nmax),s2(0:nmax)
      real rref(10),zref(10),pref(10),oldpref(10)
      real fr2(20,0:nmax)

      real minlog,potcor,psi0,outdisk,potcor1,psic,psid,mscale
      real nnn,ppp


