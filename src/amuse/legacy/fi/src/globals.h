      use memMod
      implicit none

      !                user changeable/compiled parameters
      !

      LOGICAL sortpart
      PARAMETER(sortpart=.TRUE.)
      ! whether to sort particles
            
      CHARACTER*15 metallicity
      PARAMETER(metallicity='solar')
      ! metallicity, can be:
      ! scaled from solar abundances:
      ! (Grevesse & Sauval, 1998, Space Sci. Rev. 85, 161)
      !   '2.5-solar','solar','0.2-solar','0.1-solar','0.02-solar','0.01-solar'
      ! solar, depleted Si and Fe (1/30):
      !    'FeSi-depl'
      ! solar, depleted C (1/30):
      !    'C-depl'	 
      
      LOGICAL ionCool,simpleeff,cosmrayheat,H2cooling
      PARAMETER(H2cooling=.FALSE.,ionCool=.TRUE.)
      PARAMETER(simpleeff=.FALSE.,cosmrayheat=.TRUE.)
      !  ISM model flags:
      ! H2cooling: use (approx.) molecular H cooling
      ! IonCool: use eq. ionization  cooling (otherwise use x_e=0.1)
      ! simpleeff: use constant grain heating efficiency (only exp((T/2.e4)**2) falloff) 
      ! cosmrayheat: use cosmic ray heating
      INTEGER nbodsmax,nsphmax
      PARAMETER(nbodsmax=300000,nsphmax=200000)
      ! hard limits on maximum number of particles: all & sph
      
      LOGICAL periodic
      PARAMETER (periodic=.FALSE.)
      ! periodic boundary conditions
            
      LOGICAL usepm
      INTEGER npm
      PARAMETER(usepm=.false.,npm=128)
      ! PM gravity parameters	
      
      REAL salpetertime,bh_rad_eff,bh_fdbk_eff
      INTEGER bh_method
      PARAMETER(salpetertime=450.,bh_rad_eff=0.1,bh_method=1)
      PARAMETER(bh_fdbk_eff=0.01)
      ! BH model parameters
      ! salpetertime in Myr; bh_method:
      !  1=bondi
      !  2=divv
      !  3=bondi, edd. limited
      !  4=divv, edd. limited
      
      REAL rhomin,nsupernova,esupernova
      PARAMETER(rhomin=1.e-6)
      PARAMETER(nsupernova=0.00886,esupernova=1.e51)
      ! # supernova/ Msun, kanonical energy/supernova in erg

      !    simultation parameters - organized into thematic common blocks        
      !
      COMMON/io_param/ input(128),output(128),firstsnap,stepout,steplog, &
     &        verbosity,datadir,inputfile,outputfile,halofile
      COMMON/sim_param/ pboxsize,usesph,radiate,starform,cosmo
      COMMON/units_param/ unitl_in_kpc,unitm_in_msun
      COMMON/step_param/ dtime,tstepcrit,tstpcr2,freev,freea,freevexp,   &
     &  freeaexp,nsteps,max_tbin,minppbin,sqrttstp,acc_tstp,             &
     &  freetstp
      COMMON/grav_param/ bh_tol,eps,gdgtol,nn_tol,targetnn,usequad,      &
     &  directsum,selfgrav,fixthalo,adaptive_eps,gdgop
      COMMON/sph_param/ epsgas,gamma,alpha,beta,epssph,courant,removgas, &
     &  consthsm,nsmtol,nsmooth,smoothinput,consph,sphinit,uentropy,     &
     &  isotherm,eps_is_h,hupdatemethod,sph_visc
      COMMON/rad_param/ graineff,crionrate,heat_par1,heat_par2,cool_par, &
     &  optdepth
      COMMON/cosmo_param/ comove
      COMMON/sf_param/ tcollfac,masscrit,sfeff,tbubble,sne_eff,tsnbeg,   &
     &  rhomax,sfmode,feedback
      ! rhomax( must be: < (tcollfac/dtmin)**2/pi/4)      
      
      ! 
      !
      REAL year,amu,kboltz,lightspeed,solarmass,kpc,mhydrogen
      PARAMETER(year=3.15576e7,amu=1.66054e-24,kboltz=1.38065e-16)
      PARAMETER(lightspeed=2.9979e10,solarmass=1.989e33)
      PARAMETER(kpc=3.086e21,mhydrogen=1.0079*amu)
      
      REAL pi,sqrtpi,piinv
      PARAMETER(pi=3.14159265358979323846)
      PARAMETER(piinv=0.31830988618379067154)
      PARAMETER(sqrtpi=1.77245385090551602730)
   
      !
      !
      INTEGER ncells,ninterp,ndim,nsubcell,nbodcell,nbods1
      
      PARAMETER(ncells=3*nbodsmax/5,ninterp=5000)
      PARAMETER(ndim=3,nsubcell=2**ndim)
      PARAMETER(nbodcell=nbodsmax+ncells,nbods1=nbodsmax+1)
      
      INTEGER initialseed
      PARAMETER(initialseed=7654321)
            
      INTEGER upars,ulog,ubodsin,ubodsout,uboddump
      CHARACTER*16 parsfile,logfile,dumpfile
      
      PARAMETER(upars=10,ulog=11,ubodsin=12,ubodsout=13,uboddump=14)
      PARAMETER(parsfile='runinfo',logfile='simlog',dumpfile='dumpfile')
      
      INTEGER uvfilenr,uvinterp,uvedge
      CHARACTER*16 uvfile
      
      PARAMETER(uvfile='uvpol10_10',uvfilenr=22,uvinterp=100,uvedge=10)
      
      INTEGER sfrfilenr,bhfilenr,erasenr
      PARAMETER(sfrfilenr=23,bhfilenr=24,erasenr=25)
      
      INTEGER root,nbodies,incells,nttot,ntmin,ntmax,ntavg,nsteps,       &
     &  stepout,steplog,incellsg,targetnn,nstot,nsmin,nsmax,nsavg,       &
     &  nttotfuv,ntminfuv,ntmaxfuv,ntavgfuv
      LOGICAL usequad,usesph,fixthalo,selfgrav,adaptive_eps, directsum,  &
     &  isotherm
      REAL bh_tol,eps,rsize,rmin,tnow,tpos,dtime,tiny,                   &
     &  etot,mtot,ektot,amvec,cmpos,cmvel,eptot,mstar,                   &
     &  mgas,snheat,nn_tol,esofttot,enerror
      
      INTEGER verbosity            

      INTEGER rnseed,nrndtable
      PARAMETER(nrndtable=2048)
      REAL rndtable
      COMMON/rndcom/ rndtable(nrndtable),rnseed
      
      INTEGER nsnap
      COMMON/paramcom/ nbodies,nsnap
      COMMON/cellcom/ rmin(ndim),rsize,incells,incellsg
      COMMON/pointers/ root
      COMMON/forcecom/ nttot,ntmin,ntmax,ntavg,nttotfuv,ntminfuv,        &
     &  ntmaxfuv,ntavgfuv
      COMMON/softcom/ nstot,nsmin,nsmax,nsavg
      COMMON/timecom/ tnow,tpos
      COMMON/misccom/ tiny
      COMMON/enrgycom/ mtot,etot,ektot,eptot,mstar,mgas,snheat,esofttot, &
     &  enerror,amvec(ndim),cmpos(ndim),cmvel(ndim)
      
           
      INTEGER firstsnap
      CHARACTER*16 outputfile
      CHARACTER*30 halofile
      CHARACTER*16 inputfile
      CHARACTER*200 datadir
      
      INTEGER output,input
      
      INTEGER npactive,minppbin,nsphact,max_tbin,active_bin
      REAL etol,tsteppos
      
      COMMON/actcom/ npactive,nsphact
      COMMON/stepcom/ etol,tsteppos,active_bin
      
      INTEGER syncflag, entropyflag
      COMMON/statecom/ syncflag, entropyflag
      ! syncflag=1 -> synced
      ! entropyflag=0 -> internal energy in ethermal/entropy
      
      CHARACTER*2 symmetry
      CHARACTER*4 sph_visc
      INTEGER nsph,nsmooth,nstar
      LOGICAL sphinit,starform, uentropy,eps_is_h
      REAL deldr2i,alpha,beta,epssph,gamma,gamma1,courant,wsmooth,       &
     &  dwsmooth,teth,ethtot,epsgas,consthsm,                            &
     &  nsmtol,massres,tstarform,tsnfeedback
            
      COMMON/ethcom/ teth,ethtot
      COMMON/symmcom/ symmetry
      COMMON/sphparam/ gamma1,nstar,nsph
      COMMON/smoocom/ massres
      COMMON/interpoc/ deldr2i
      COMMON/skerncom/ wsmooth(0:1+ninterp),dwsmooth(0:1+ninterp)
      
      INTEGER nntot,nnmin,nnmax,nnavg
      
      COMMON/neighcom/ nntot,nnmin,nnmax,nnavg
      
      LOGICAL cosmo,comove
      REAL pboxsize,hboxsize
      
      COMMON/misccom/ hboxsize
            
      REAL poshalo,masshalo
      COMMON/halocom/ poshalo(3),masshalo
            
      LOGICAL radiate
      REAL heat_par1,heat_par2,cool_par,eradiate,trad,                   &
     &  meanmwt,mhboltz,fhydrogn,unitm_in_msun,unitl_in_kpc,                       &
     &  mumhkbol,mumhkgam,mumhkg1,optdepth,graineff,                     &
     &  efuvheat,eradcool,mcold,mluke,mwarm,mhot,estar,heatconst,        &
     &  densconst,timescale,lengthscale,velscale,flxscale,crionrate,     &
     &  massscale
      
      COMMON/radcom/ eradiate,trad,meanmwt,                              &
     &  mhboltz,fhydrogn,mumhkbol,mumhkgam,                              &
     &  mumhkg1,efuvheat,eradcool,mcold,mluke,mwarm,mhot,estar,          &
     &  heatconst,densconst,timescale,lengthscale,velscale,flxscale,     &
     &  massscale
      
      CHARACTER*4 feedback
      INTEGER totptag
      
      REAL tcollfac,masscrit,sfeff,removgas,tsnbeg,tbubble,              &
     &  sne_eff,snenergy
      
      COMMON/snsfcom1/ tstarform,tsnfeedback,snenergy,totptag
            
      REAL h2time
      REAL smoothuv
      
      COMMON/horh2/ h2time
      COMMON/comgrav/ smoothuv(0:uvinterp)
      
      ! blackhole stuff
      REAL tbh
      INTEGER nbh
      COMMON/bhcom/ tbh,nbh
      
      CHARACTER*4 hupdatemethod
      CHARACTER*10 sfmode
      
      LOGICAL gdgop,sqrttstp,acc_tstp,freetstp,consph,smoothinput
      
      REAL gdgtol,tstepcrit,tstpcr2
      REAL freev,freea,freevexp,freeaexp,rhomax
      
      INTEGER ntreemin
      PARAMETER(ntreemin=100)
      
      ! pm stuff
      REAL tpm,rcut,rcut2,pmsoft,pmpot,pmacc,pmdr
      COMMON/pmcom/ tpm,rcut,rcut2,pmsoft,pmdr,pmpot(0:ninterp+1),       &
     &  pmacc(0:ninterp+1)
      
      REAL searchmagic
      PARAMETER(searchmagic=0.4)
      INTEGER searchn,searchreuse,reuseflag,ncalls,nsearches
      REAL searchpos,searchh,searchdelta,searchacc4
      COMMON/searchcom/ searchpos(3),searchh,searchdelta,searchacc4,     &
     &  searchn,searchreuse,reuseflag,ncalls,nsearches
!$omp   THREADPRIVATE(/searchcom/)
      
      
      INTEGER peanovalue(8,24)
      DATA peanovalue/ 1, 8, 4, 5, 2, 7, 3, 6, 1, 2, 4, 3, 8, 7, 5, 6,   &
     &  7, 8, 6, 5, 2, 1, 3, 4, 3, 6, 4, 5, 2, 7, 1, 8, 1, 4, 8, 5, 2,   &
     &  3, 7, 6, 5, 8, 4, 1, 6, 7, 3, 2, 1, 2, 8, 7, 4, 3, 5, 6, 3, 2,   &
     &  4, 1, 6, 7, 5, 8, 7, 2, 6, 3, 8, 1, 5, 4, 5, 6, 4, 3, 8, 7, 1,   &
     &  2, 7, 8, 2, 1, 6, 5, 3, 4, 7, 6, 8, 5, 2, 3, 1, 4, 3, 4, 6, 5,   &
     &  2, 1, 7, 8, 1, 4, 2, 3, 8, 5, 7, 6, 5, 4, 8, 1, 6, 3, 7, 2, 5,   &
     &  8, 6, 7, 4, 1, 3, 2, 1, 8, 2, 7, 4, 5, 3, 6, 7, 2, 8, 1, 6, 3,   &
     &  5, 4, 5, 6, 8, 7, 4, 3, 1, 2, 3, 2, 6, 7, 4, 1, 5, 8, 7, 6, 2,   &
     &  3, 8, 5, 1, 4, 5, 4, 6, 3, 8, 1, 7, 2, 3, 4, 2, 1, 6, 5, 7, 8,   &
     &  3, 6, 2, 7, 4, 5, 1, 8 /
      
      INTEGER pshapefromi(8,24)
      DATA pshapefromi/ 2, 3, 4, 4, 5, 6, 5, 6, 1, 7, 8, 7, 9, 10, 8,    &
     &  10, 11, 1, 11, 12, 13, 9, 13, 12, 12, 8, 1, 1, 12, 8, 10, 13,    &
     &  14, 13, 12, 13, 1, 1, 15, 15, 10, 16, 10, 8, 1, 1, 15, 15, 17,   &
     &  2, 18, 19, 20, 2, 20, 19, 18, 18, 2, 6, 4, 4, 2, 20, 21, 20,     &
     &  21, 20, 2, 3, 18, 18, 6, 2, 6, 19, 22, 2, 4, 19, 3, 17, 23, 18,  &
     &  3, 21, 23, 21, 18, 18, 5, 3, 4, 4, 21, 3, 3, 5, 23, 5, 3, 22,    &
     &  23, 4, 5, 23, 17, 17, 21, 23, 22, 22, 22, 22, 19, 23, 5, 6, 5,   &
     &  6, 19, 6, 17, 17, 19, 20, 22, 22, 7, 11, 14, 16, 24, 24, 14,     &
     &  16, 12, 8, 7, 11, 12, 8, 9, 9, 16, 7, 15, 7, 16, 10, 24, 10, 9,  &
     &  9, 24, 24, 7, 16, 7, 8, 9, 9, 24, 24, 14, 11, 12, 11, 15, 15,    &
     &  14, 16, 10, 13, 14, 16, 11, 14, 11, 15, 13, 14, 13, 24, 21, 20,  &
     &  21, 20, 17, 17, 19, 23 /
      
      INTEGER pshapefromp(8,24)
      DATA pshapefromp/ 2, 5, 5, 4, 4, 6, 6, 3, 1, 7, 7, 8, 8, 10, 10,   &
     &  9, 9, 13, 13, 12, 12, 11, 11, 1, 10, 12, 12, 1, 1, 8, 8, 13,     &
     &  14, 1, 1, 13, 13, 15, 15, 12, 8, 15, 15, 10, 10, 1, 1, 16, 17,   &
     &  2, 2, 20, 20, 19, 19, 18, 6, 18, 18, 2, 2, 4, 4, 20, 3, 20, 20,  &
     &  18, 18, 21, 21, 2, 4, 19, 19, 6, 6, 2, 2, 22, 18, 23, 23, 21,    &
     &  21, 3, 3, 17, 21, 4, 4, 3, 3, 18, 18, 5, 22, 3, 3, 5, 5, 23,     &
     &  23, 4, 5, 17, 17, 23, 23, 22, 22, 21, 23, 6, 6, 22, 22, 5, 5,    &
     &  19, 20, 22, 22, 19, 19, 17, 17, 6, 7, 14, 14, 24, 24, 16, 16,    &
     &  11, 11, 8, 8, 9, 9, 12, 12, 7, 24, 10, 10, 16, 16, 7, 7, 15,     &
     &  16, 9, 9, 7, 7, 24, 24, 8, 12, 24, 24, 11, 11, 9, 9, 14, 13,     &
     &  16, 16, 15, 15, 14, 14, 10, 15, 11, 11, 14, 14, 13, 13, 24, 19,  &
     &  21, 21, 17, 17, 20, 20, 23 /
      
      INTEGER indexfromp(8,24)
      DATA indexfromp/ 1, 5, 7, 3, 4, 8, 6, 2, 1, 2, 4, 3, 7, 8, 6, 5,   &
     &  6, 5, 7, 8, 4, 3, 1, 2, 7, 5, 1, 3, 4, 2, 6, 8, 1, 5, 6, 2, 4,   &
     &  8, 7, 3, 4, 8, 7, 3, 1, 5, 6, 2, 1, 2, 6, 5, 7, 8, 4, 3, 4, 2,   &
     &  1, 3, 7, 5, 6, 8, 6, 2, 4, 8, 7, 3, 1, 5, 7, 8, 4, 3, 1, 2, 6,   &
     &  5, 4, 3, 7, 8, 6, 5, 1, 2, 7, 5, 6, 8, 4, 2, 1, 3, 6, 5, 1, 2,   &
     &  4, 3, 7, 8, 1, 3, 4, 2, 6, 8, 7, 5, 4, 8, 6, 2, 1, 5, 7, 3, 6,   &
     &  8, 7, 5, 1, 3, 4, 2, 1, 3, 7, 5, 6, 8, 4, 2, 4, 2, 6, 8, 7, 5,   &
     &  1, 3, 7, 8, 6, 5, 1, 2, 4, 3, 6, 2, 1, 5, 7, 3, 4, 8, 7, 3, 4,   &
     &  8, 6, 2, 1, 5, 6, 8, 4, 2, 1, 3, 7, 5, 4, 3, 1, 2, 6, 5, 7, 8,   &
     &  7, 3, 1, 5, 6, 2, 4, 8 /
