c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MERCURY6_1.FOR    (ErikSoft   3 May 2002)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Mercury is a general-purpose N-body integration package for problems in
c celestial mechanics.
c
c------------------------------------------------------------------------------
c This package contains some subroutines taken from the Swift integration 
c package by H.F.Levison and M.J.Duncan (1994) Icarus, vol 108, pp18.
c Routines taken from Swift have names beginning with `drift' or `orbel'.
c
c The standard symplectic (MVS) algorithm is described in J.Widsom and
c M.Holman (1991) Astronomical Journal, vol 102, pp1528.
c
c The hybrid symplectic algorithm is described in J.E.Chambers (1999)
c Monthly Notices of the RAS, vol 304, pp793.
c
c RADAU is described in E.Everhart (1985) in ``The Dynamics of Comets:
c Their Origin and Evolution'' p185-202, eds. A.Carusi & G.B.Valsecchi,
c pub. Reidel.
c
c The Bulirsch-Stoer algorithms are described in W.H.Press et al. (1992)
c ``Numerical Recipes in Fortran'', pub. Cambridge.
c------------------------------------------------------------------------------
c
c Variables:
c ---------
c  M      = mass (in solar masses)
c  XH     = coordinates (x,y,z) with respect to the central body (AU)
c  VH     = velocities (vx,vy,vz) with respect to the central body (AU/day)
c  S      = spin angular momentum (solar masses AU^2/day)
c  RHO    = physical density (g/cm^3)
c  RCEH   = close-encounter limit (Hill radii)
c  STAT   = status (0 => alive, <>0 => to be removed)
c  ID     = name of the object (8 characters)
c  CE     = close encounter status
c  NGF    = (1-3) cometary non-gravitational (jet) force parameters
c   "     =  (4)  beta parameter for radiation pressure and P-R drag
c  EPOCH  = epoch of orbit (days)
c  NBOD  = current number of bodies (INCLUDING the central object)
c  NBIG  =    "       "    " big bodies (ones that perturb everything else)
c  TIME  = current epoch (days)
c  TOUT  = time of next output evaluation
c  TDUMP = time of next data dump
c  TFUN  = time of next periodic effect (e.g. next check for ejections)
c  H     = current integration timestep (days)
c  EN(1) = initial energy of the system
c  " (2) = current    "    "  "    "
c  " (3) = energy change due to collisions, ejections etc.
c  AM(1,2,3) = as above but for angular momentum
c
c Integration Parameters :
c ----------------------
c  ALGOR = 1  ->  Mixed-variable symplectic
c          2  ->  Bulirsch-Stoer integrator
c          3  ->         "           "      (conservative systems only)
c          4  ->  RA15 `radau' integrator
c          10 ->  Hybrid MVS/BS (democratic-heliocentric coords)
c          11 ->  Close-binary hybrid (close-binary coords)
c          12 ->  Wide-binary hybrid (wide-binary coords)
c
c TSTART = epoch of first required output (days)
c TSTOP  =   "      final required output ( "  )
c DTOUT  = data output interval           ( "  )
c DTDUMP = data-dump interval             ( "  )
c DTFUN  = interval for other periodic effects (e.g. check for ejections)
c  H0    = initial integration timestep (days)
c  TOL   = Integrator tolerance parameter (approx. error per timestep)
c  RMAX  = heliocentric distance at which objects are considered ejected (AU)
c  RCEN  = radius of central body (AU)
c  JCEN(1,2,3) = J2,J4,J6 for central body (units of RCEN^i for Ji)
c
c Options:
c  OPT(1) = close-encounter option (0=stop after an encounter, 1=continue)
c  OPT(2) = collision option (0=no collisions, 1=merge, 2=merge+fragment)
c  OPT(3) = time style (0=days 1=Greg.date 2/3=days/years w/respect to start)
c  OPT(4) = o/p precision (1,2,3 = 4,9,15 significant figures)
c  OPT(5) = < Not used at present >
c  OPT(6) = < Not used at present >
c  OPT(7) = apply post-Newtonian correction? (0=no, 1=yes)
c  OPT(8) = apply user-defined force routine mfo_user? (0=no, 1=yes)
c
c File variables :
c --------------
c  OUTFILE  (1) = osculating coordinates/velocities and masses
c     "     (2) = close encounter details
c     "     (3) = information file
c  DUMPFILE (1) = Big-body data
c     "     (2) = Small-body data
c     "     (3) = integration parameters
c     "     (4) = restart file
c
c Flags :
c -----
c  NGFLAG = do any bodies experience non-grav. forces?
c                            ( 0 = no non-grav forces)
c                              1 = cometary jets only
c                              2 = radiation pressure/P-R drag only
c                              3 = both
c  OPFLAG = integration mode (-2 = synchronising epochs)
c                             -1 = integrating towards start epoch
c                              0 = main integration, normal output
c                              1 = main integration, full output
c
c------------------------------------------------------------------------------
c
      implicit none
      include 'mercury.inc'
c
      integer j,algor,nbod,nbig,opt(8),stat(NMAX),lmem(NMESS)
      integer opflag,ngflag,ndump,nfun
      real*8 m(NMAX),xh(3,NMAX),vh(3,NMAX),s(3,NMAX),rho(NMAX)
      real*8 rceh(NMAX),epoch(NMAX),ngf(4,NMAX),rmax,rcen,jcen(3)
      real*8 cefac,time,tstart,tstop,dtout,h0,tol,en(3),am(3)
      character*8 id(NMAX)
      character*80 outfile(3), dumpfile(4), mem(NMESS)
      external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy
      external mco_dh2h,mco_h2dh
      external mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden
c
      data opt/0,1,1,2,0,1,0,0/
c
c------------------------------------------------------------------------------
c
c Get initial conditions and integration parameters
      call mio_in (time,tstart,tstop,dtout,algor,h0,tol,rmax,rcen,jcen,
     %  en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,id,
     %  epoch,ngf,opt,opflag,ngflag,outfile,dumpfile,lmem,mem)
c
c If this is a new integration, integrate all the objects to a common epoch.
      if (opflag.eq.-2) then
  20    open (23,file=outfile(3),status='old',position='append',err=20)
        write (23,'(/,a)') mem(55)(1:lmem(55))
        write (*,'(a)') mem(55)(1:lmem(55))
        call mxx_sync (time,tstart,h0,tol,jcen,nbod,nbig,m,xh,vh,s,rho,
     %    rceh,stat,id,epoch,ngf,opt,ngflag)
        write (23,'(/,a,/)') mem(56)(1:lmem(56))
        write (*,'(a)') mem(56)(1:lmem(56))
        opflag = -1
        close (23)
      end if
c
c Main integration
      if (algor.eq.1) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_mvs,mco_h2mvs,mco_mvs2h)
c
      if (algor.eq.9) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_mvs,mco_iden,mco_iden)
c
      if (algor.eq.2) call mal_hvar (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_bs1)
c
      if (algor.eq.3) call mal_hvar (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_bs2)
c
      if (algor.eq.4) call mal_hvar (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_ra15)
c
      if (algor.eq.10) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_hy,mco_h2dh,mco_dh2h)
c
c Do a final data dump
      do j = 2, nbod
        epoch(j) = time
      end do
      call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,
     %  rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,
     %  id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
c
c Calculate and record the overall change in energy and ang. momentum
  50  open  (23, file=outfile(3), status='old', position='append',
     %  err=50)
      write (23,'(/,a)') mem(57)(1:lmem(57))
      call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
c
      write (23,231) mem(58)(1:lmem(58)), 
     %  abs((en(2) + en(3) - en(1)) / en(1))
      write (23,232) mem(59)(1:lmem(59)), 
     %  abs((am(2) + am(3) - am(1)) / am(1))
      write (23,231) mem(60)(1:lmem(60)), abs(en(3) / en(1))
      write (23,232) mem(61)(1:lmem(61)), abs(am(3) / am(1))
      close (23)
      write (*,'(a)') mem(57)(1:lmem(57))
c
c------------------------------------------------------------------------------
c
 231  format (/,a,1p1e12.5)
 232  format (a,1p1e12.5)
      stop
      end
