*       COMMON.H
*       --------
*
*
*       MONT-CAR commons
*       ----------------
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical isnan
*
      include 'params.h'

      CHARACTER*200 datadir
*
*
*      REAL*8  BODY,R,VR,VT,U,RO,VRO,VTO,UO,BIN2,BIN3,BINP,BIN4,SMTO,
*     &        BIN2IN,BIN3IN,BIN4IN,BINPIN,SMT,ETOT,ESCSTA,ENRAD,
      REAL*8  BODY,R,VR,VT,U,RO,VRO,VTO,UO,BIN,SMTO,TSCALE,DTTE0,XTAU,
     &        BININ,SMT,ETOT,ESCSTA,ENRAD,EHB3B3,EHBI3P,TCREVO,
     &        ESCBI2,ESCBI3,ESCBI4,ESCBIP,ESCBB2,ESCBB3,ESCBB4,ESCBBP,
     &        EHMLEV,ESCB2S,ESCB3S,ESCB4S,ESCBPS,EHBIN2,EHBIN3,TTE,
     &        EHBIN4,EHBINP,ECBIN2,ECBIN3,ECBIN4,ECBINP,SLOSEV,SLOSCO,
     &        ECCOLL,EHCOLL,PBIN2,PBIN3,PBIN4,PBIN2S,PBIN3S,PBIN4S,
     &        PBINPS,PB2B2,PB2B3,PB2B4,PB2BP,PB3B3,PB3B4,PB3BP,PB4B4,
     &        PB4BP,PBPBP,PCOLL,BMIN,BMAX,TAU0,TCRIT,TCOMP,QE,ALPHAL,
     &        ALPHAH,BRAKEM,BODY1,BODYM,BODYN,QVIR,RBAR,ZMBAR,W0,ZKIN,
     &        POT,RSCALE,TCR,TIME,ONE2,ONE3,ONE4,ONE5,ONE6,ONE9,ONE12,
     &        PI,TWOPI,SLOSES,FLAGR,SMC,VC,ROC,RC,RTID,TAU,ERROR,CPU,
     &        RN,OLDTIME,DTTP,SLOB3S,SLOB3B,XESCP,XESCT,ERELB3,ERB3IN,
     &        XMIN,XMAX,DTTE,XMAXZ,VRR,XGMIN,XGMAX,TTP,FRACB,AMIN,AMAX,
     &        UPTIME,RPLUM,TSCALE0,TIMET,XTID,BMIN0,YTAU,YBMIN,EKICKT,
     &        ZINI,TOLM,SESCRT,EKICKTBS,EKICKTBD,EKICKTB2
*
      REAL*8  X,XDOT,RLAG,V2RL,V2TL,ANI,RTIDKG,RN1,RN2,VRN1,VRN2,
     &        VTN1,VTN2,ESCENE,ESCMAS,ESCDIS,ESCANG,RSUNTOPC,ENEKIN,
     &        SINFI,COSFI,FI,RMGLOB,GAMMA,ESCTIM,ENEPOT,VRP,VTP

      INTEGER  NZS,NZST,LTWO,NT0,NT,NESCST,NESCB2,NTO,IKIND,IINTE3,
     &         NESCB3,NESCB4,NESCBP,NESB2S,NESB3S,NESB4S,NESBPS,
     &         NMLOEV,NMLOCO,NCOLL,NECB2,NECB3,NECB4,NECBP,NEHB2,
     &         NEHB3,NEHB4,NEHBP,NBIN2,NBIN3,NBIN4,NBINP,IRUN,IB3F,	
     &         ISTART,NCOR,NMIN,NZ0,NZONC,NMINZO,NTWO,IPRINT,IMODEL,
     &         INAME,NC,IZON,NSUZON,INAMEO,NZON0,IVNEWG,IVRR,ISEED,
     &         NB3FIN,NB3B3,NDIST3,NDIST4,NSS0,NAMES,NAMEB,NWHICH,
     &         IBSM,IBSC,IBS3,IBS4,NMERGE,IBSTRA,NESCM,NDISTE,NT00,
     &         NBINAR,IKICKT,MNSBH,IEXCH,NEXCHANG,INEXCH,NEXCHANG2,
     &         INEXCH2,IFLAGNS,IFLAGBH,NTNEW,NKICK,IOBT,NESCRT,
     &         IKROUPA,IKICKTBS,IKICKTBD,NTSN1,NTSN2,NTSNB,NITESC
*
      REAL*4 TIMEOLD
c      INTEGER  IENAME,IESC,IESC1,IENAM1,IEKIND,IDUM2,IY,IV
      INTEGER  IENAME,IESC,IESC1,IENAM1,IEKIND,ITIME
*
      INTEGER  IDUM2,IY,IV
*
*
*       ----------------------------------------------------------------------
*
*                        VARIABLES
*                        ---------
*
*       TIME    -  evolution time in N-body units
*       TIMET   -  evolution time at the end of the last step (Myr)
*       UPTIME  -  evolution time of single stars and binaries
*       OLDTIME -  time of the last evolution step
*       TSCALE  -  scale time to the physical units (Myr) -
*                  TSCALE = TSCALE0*ln(gammaN0)/ln(gammaN)
*       TSCALE0 -  initial scale time to the physical units (Myr)
*       TCRIT   -  termination time in units of the crossing time
*       TCOMP   -  maximum computing time in minutes
*       TIMEOLD -  the runtime time from the previous time step
*       ITIME   -  number of 24h passed from the beginning of the run
*       QE      -  energy tolerance
*       ALPHAL  -  power-law index for initial mass function for masses
*                  smaller than breake mass: -1 - equal mass case
*       ALPHAH  -  power-law index for initial mass function for masses
*                  greater than breake mass. If alphal=alphah the IMF does
*                  not have a break
*       BRAKEM  -  the mass in which the IMF is broken. If brakem is smaller
*                  than the minimum mass (bodyn) than the break mass is as
*                  for the Kroupa mass function (brakem = 0.5 Mo)
*       BODY1   -  maximum particle mass before scaling (solar mass)
*       BODYM   -  avarage particle mass in scaled units
*       BODYN   -  minimum particle mass before scaling (solar mass)
*       FRACB   -  primordial binary fraction by number. nb = fracb*nt,
*                  ns = (1 - fracb)*nt, nss = (1 + fracb)*nt
*                  fracb > 0 - primordial binaries
*                  fracb = 0 - only dynamical binaries
*       AMIN    -  minimum semi-major axis of binaries (in sollar units)
*       AMAX  -  maximum semi-major axis of binaries (in sollar units)
*       QVUR    -  virial ratio  (qvir = 0.5 for equilibrium)
*       RBAR    -  mean cluster radius in pc (set = 0 for isolated cluster)
*       ZMBAR   -  mean mass in solar units (set > 0 for isolated cluster)
*       W0      -  king model parameter
*       ZKIN    -  total kinetic energy
*       POT     -  total potential energy
*       RSCALE  -  scale radius
*       TCR     -  crrosing time
*       FLAGR   -  fraction of the total mass (definition of lagrangian radii)
*       RTID    -  tidal radius
*       BMIN    -  minimum value of sin(beta^2/2), at any time
*       BMIN0   -  BMIN at time = 0
*       BMAX    -  maximum value of sin(beta^2/2)
*       XTID    -  coeficient in front of cluster tidal energy
*                                                         -xtid*smt/rtid
*       TAU0    -  time-step for a complite cluster model at T=0
*       TAU     -  time-step for a complite cluster model at any time
*       ERROR   -  total error of energy scaled by kinetic energy
*       CPU     -  CPU time for curent run
*       BODY    -  masses of single bodies (s.b.) and centre of masses
*                  of multiple configurations (c.m.)
*       GAMMA   -  parameter in Coulomb logarithm (standard value = 0.11)
*       RPLUM   -  for M67 rtid = rplum*rsplum (rsplum - scale radius for
*                  plummer model)
*       DTTP    -  time step (Myr) for profile output
*       DTTE    -  time step (Myr) for mloss call for all objects for tphys
*                  greater then tcrevo. For tphys less then tcrevo time step
*                  is eqiual to dtte0
*       DTTE0   -  time step (Myr) for mloss call for all objects for tphys
*                  less then tcrevo. For tphys greater then tcrevo time step
*                  is eqiual to dtte
*       TCREVO  -  critical time for which time step for mloss call changes 
*                  from dtte0 to dtte
*       XTAU    -  call mloss for a particlular object when
*                  (uptime(im1) - olduptime(im1))/tau/tscale < xtau
*       YTAU    -  multiplication of TAU0 (TAU = YTAU*TAU0) after time
*                  greater than TCREVO
*       YBMIN   -  multiplication of BMIN0 (BMIN = YBMIN*BMIN0) after time
*                  greater than TCREV0
*       ZINI    -  initial metalicity (solar z = 0.02)
*       IKROUPA -  0 - the initial binary parameters are picked up
*                  according Kroupa's eigenevolution and feeding algorithm
*                  (Kroupa 1995, MNRAS 277, 1507)
*                  1 - the initial binary parameters are picked as for M67
*                  model (Hurley et al. 2005)
*       NITESC  -  0 - no iteration of the tidal radius and induced mass loss
*                  due to stellar evolution, 1 - iteration of the tidal radius
*                  and induced mass loss due to stellar evolution
*       IFLAGNS -  0 - no SN natal kiks for NS, 1 - SN natal kicks for NS
*       IFLAGBH -  0 - no SN natal kiks for BH, 1 - SN natal kicks for BH
*       TOLM    -  minimum mass for which the object mass is set to 0.0
*                  tolm = 1.0d-5/zmbar
*       SMC     -  core mass
*       RC      -  core radius
*       VC      -  central velocity
*       ROC     -  central density
*       NC      -  number of stars in the core       
*       R       -  radial positions of s.b. and c.m.
*       RN      -  new R
*       RO      -  old R
*       VR      -  radial velocities of s.b. and c.m.
*       VRO     -  old VR
*       VT      -  tangential velocities of s.b. and c.m.
*       VTO     -  old VT
*       U       -  potentials of s.b. and c.m.
*       UO      -  old U
*       AZ      -  coefficients in formulae for sin(beta^2/2)
*       SMT     -  cluster total mass
*       ETOT    -  cluster total binding energy
*       RMGLOB  -  minimum value for minimum distance from the cluster centre
*                  equal to smt/nt
*       XESCP   -  position of o star which for the first time has binding
*                  energy greater than zero
*       XESCT   -  the same as for XESCP but for time
*       VRR     -  keeps for each star the energy error due to potential
*                  adjustment
*       BIN     -  binary parameters
*                  1 -  mass of the first component
*                  2 -  mass of the second component
*                  3 -  semi-major axis
*                  4 -  eccentricity
*                  5 -  binding energy
*                  6 -  time of live
*                  7 -  number of c.m. in the global list
*                  8 -  position of the binary
*       BININ  -  information about binaries interactions
*                  1  - position
*                  2  - energy change
*                  3  - mass 1 + mass 2
*                  4  - mass 3 + mass 4
*                  5  - time
*                  6  - binding energy
*                  7  - 0 - formation, 1 - single,  2 - binary, 3 - merger
*                       4 - binary-single merger, 5 - binary-binary merger
*                       6 - dsruption evolution, 7 - disruption binary-binary
*                       8 - disruption binary-single, 9 - 3b binary exchange
*                       10 - binary-binary exchange, 13 - merger after
*                            binary exchange, 15 - merger after binary-binary
*                            exchange
*                       -1 - binary evolution only for binarym.dat printout
*       EKICKT  -  kinetic energy of NS/BH after supernove explosion and 
*                  natal kick
*       EKICKTBS-  kinetic energy of binary CM after supernove explosion and
*                  natal kick
*       EKICKTBD-  kinetic energy of binary components after supernove 
*                  explosion binary disruption and natal kick 
*       EKICKTB2-  kinetic energy of second binary component after supernove
*                  explosion binary disruption and natal kick
*       ESCSTA  -  energy of escapers      
*       ESCBI2  -  energy of two-body escapers
*       ESCBI3  -  energy of three-body escapers
*       ESCBI4  -  energy of four body escapers
*       ESCBIP  -  energy of primprodial escapers
*       ESCBB2  -  internal binding energy of two-body escapers
*       ESCBB3  -  internal binding energy of three-body escapers
*       ESCBB4  -  internal binding energy of four-body escapers
*       ESCBBP  -  internal binding energy of primorodial escapers
*       ESCB2S  -  energy of star escapers - two-body binaries
*       ESCB3S  -  energy of star escapers - three-body binaries
*       ESCB4S  -  energy of star escapers - four-body binaries
*       ESCBPS  -  energy of star escapers - primorodial binaries
*       EHBIN2  -  heating by two-body binaries
*       EHBIN3  -  heating by three-body binaries
*       EHBI3P  -  initial binding energy of primordial binaries
*       EHB3B3  -  heating by three-body--three-body interactions
*       EHBIN4  -  heating by four-body binaries
*       EHBINP  -  heating by primorodial binaries
*       ECBIN2  -  cooling by two-body binaries
*       ECBIN2  -  cooling by three-body binaries
*       ECBIN2  -  cooling by four-body binaries
*       ECBIN2  -  cooling by primorodial binaries
*       EHMLEV  -  heating by mass loss due to stellar evolution
*       ERELB3  -  energy of three-body binary escapers due to relaxation
*       ERB3IN  -  internal binding energy of three-body escapers due to
*                  relaxation
*       ENEPOT  -  energy put to the system due to collisions, formation
*                   or dissruption of binaries
*       ENEKIN  -  kinetic energu changes due to collisions
*       ECCOLL  -  cooling due to stellar collisions
*       EHCOLL  -  heating due to stellar collisions
*       SLOSES  -  mass loss due to the relaxation process
*       SLOB3S  -  mass loss - stars in interaction between three-body
*                  binaries and field stars
*       SLOB3B  -  mass loss - binaries in interaction between three-body
*                  binaries and field stars
*       SLOSEV  -  mass loss due to stellar evolutin
*       SLOSCO  -  mass loss due to stellar collisions
*       SESCRT  -  mass loss due to adjustment ot the tidal radius 
*       PBIN2   -  cumulative probability of two-body binary formation
*       PBIN3   -  cumulative probability of three-body binary formation
*       PBIN4   -  cumulative probability of four-body binary formation
*       PBIN2S  -  cumulative probability of two-body star interaction
*       PBIN3S  -  cumulative probability of three-body star interaction
*       PBIN4S  -  cumulative probability of four-body star interaction
*       PBINPS  -  cumulative probability of primorodial star interaction
*       PB2B2   -  cumulative probability of two-body two-body interaction
*       PB2B3   -  cumulative probability of two-body three-body interaction
*       PB2B4   -  cumulative probability of two-body four-body interaction
*       PB2BP   -  cumulative probability of two-body primorodial interaction
*       PB3B3   -  cumulative probability of three-body three-body interaction
*       PB3B4   -  cumulative probability of three-body four-body interaction
*       PB3BP   -  cumulative probability of three-body primorodial interaction
*       PB4B4   -  cumulative probability of four-body four-body interaction
*       PB4BP   -  cumulative probability of four-body primorodial interaction
*       PBPBP   -  cumulative probability of primorodial primorodial interaction
*       PB2BP   -  cumulative probability of two-body primorodial interaction
*       PCOLL   -  cumulative probability of collisions
*
*       IKIND   -  type of stars: 1 - single star, 2 - binary, 3 - binary
*                  mergger, 4 - single stars collision merger 
*       IINTE3  -  numerate interactions with field stars of each three-body
*                  binary
*       IRUN    -  initial parameter for random number generator
*       ISEED   -  descriptor of a run
*       NT00    -  total number of objects (stars and binaries) at T=0
*                  ns - number of single stars, nb - number of binaries
*                  (nt = ns + nb), nss - number of stars (nss = nt + nb)
*       NT0     -  total number of objects - initial plus new objects after
*                  binary disruption
*       NSS0    -  initial number of single stars  (T=0)
*       NT      -  total number of objects at any time
*       NTNEW   -  number of objects just before the next overall time step
*       ISTART  -  1 - initial model,    .ne.1 - restart
*       NCOR    -  number of stars to calculate the central parameters
*       NMIN    -  minimum number of stars to calculate the central parameters
*       NZ0     -  number of stars in each zone at T=0
*       IZON    -  = 1 then do not change number of zones and stars in them
*       NSUZON  -  number of super-zones in the cluster
*       NZONC   -  minimum number of zones in the core
*       NMINZO  -  minimum number of stars in a zone
*       NTWO    -  maximum index of two 
*       IMODEL  -  initial model: 1- uniform & isotropic, 2- Plummer, 3- King
*       IPRINT  -  0- full diagnostic info., 1- diagnostic info. suppressed 
*       IB3F    -  1 - Spitzer's, 2 - Heggie's formula for three-body binary
*                  interaction with field stars, 3 - use Pmax for interaction
*                  probability  4 - three- and four-body numerical integration
*       IEXCH   -  0 - no exchange in any interactions, 1 - exchange only in 
*                  binary field star interactions, 2 - exchange in all
*                  interactions (binary - field and binary - binary)
*      NEXCHANG -  Total number of exchanges in 3b interactions
*     NEXCHANG2 -  Total number of exchanges in 4b interactions
*       NZS     -  number of the last star in each zone
*       NZST    -  number of the last star in each super-zone
*       NZSTE   -  number of the last star in each super-zone reduced by escapers
*       LTWO    -  index for each super-zone
*       INAME   -  mames of particles (after sort)
*       INAMEO  -  old INAME
*       NAMES   -  name of single star (initial name)
*       NAMEB   -  name of primordial and formed binaries (initial name)
*       NKICK   -  flag for object (NS/BH) which received natal kick because
*                  of supernove explosion (0 - nokick, 1 - kick) 
*       INEXCH  -  number of 3b star exchanges for a particular binary
*       INEXCH2 -  number of 4b star exchanges for a particular binary
*       NWHICH  -  returns the binary possition in the BIN aray in the 
*                  case of binary formation
*       NBINAR  -  returns the object name from the position in the BIN
*                  aray
*       IBSTRA  -  keep info about blue stargglers and chanels of their
*                  formation: 1 - mergers, 2 - collisions, 3 - 3-body
*                  interactions, 4 - 4-body interactions, 0 - not blue
*                  straggler 
*       IB2F    -  1 - Spitzer's, 2 - Heggie's formula for three-body binary
*                  interaction with field stars, 3 - use Pmax for interaction
*                  probability  4 - three- and four-body numerical integration
*       NLAGRA  -  number of Lagrangian radii
*       NMERGE  -  total number of binary mergers
*       NESCST  -  total number of star escapers
*       NESCM   -  total number of massless objects removed from the system
*       NESCB2  -  total number of two-body binary escapers
*       NESCB3  -  total number of three-body bimary escapers
*       NESCB4  -  total number of four-body escapers
*       NESCBP  -  total number of primorodial binaries escapers
*       NESB2S  -  total number of star escapers from two-body binaries
*       NESB3S  -  total number of star escapers from three-body binaries
*       NESB4S  -  total number of star escapers from four-body binaries
*       NESBPS  -  total number of star escapers from primorodial binaries
*       NESCRT  -  total number of escapers due to adjustment ot the 
*                  tidal radius
*       NMLOEV  -  total number of mass loss due to stellar evolution
*       NMLOCO  -  total number of mass loss due to stellar collisions
*       NCOLL   -  total number of stellar collisions
*       NECB2   -  total number of cooling events for two-body binaries
*       NEHB2   -  total number of heating events for two-body binaries
*       NECB3   -  total number of cooling events for three-body binaries
*       NEHB3   -  total number of heating events for three-body binaries
*       NECB4   -  total number of cooling events for four-body binaries
*       NEHB4   -  total number of heating events for four-body binaries
*       NECBP   -  total number of cooling events for primorodial binaries
*       NEHBP   -  total number of heating events for primorodial binaries
*       NBIN2   -  total number of two-body binaries
*       NBIN3   -  total number of three-body binaries
*       NB3FIN  -  total number of interactions between three-body binaries
*                  and field stars
*       NB3B3   -  total number of interactions between three-body binaries
*       NBIN4   -  total number of four-body binaries
*       NBINP   -  total number of primorodial binaries
*       NZON0   -  total number of zones at T=0
*       NDIST3  -  total number of binaries destroyed in 3-body interactions
*       NDIST4  -  total number of binaries destroyed in 4-body interactions
*       NDISTE  -  total number of binaries destroyed due to binary evolution
*       IBSM    -  total number of blue stragglers formed in binary mergers
*       IBSC    -  total number of blue stragglers formed in collisions
*       IBS3    -  total number of blue stragglers formed in 3-body 
*                  interactions
*       IBS4    -  total number of blue stragglers formed in 3-body 
*                  interactions
*       IKICKT  -  total number of supenove explosions - NS/BH natal kicks
*       IKICKTBS-  total number of supenove explosions in binaries - binary
*                  survived NS/BH natal kick
*       IKICKTBD-  total number of supenove explosions in binaries - binary   
*                  not survived NS/BH natal kick
*       MNSBH   -  number of NS/BH escaped from the system
*       IOBT    -  total number of obliterated object
*       NTSN1   -  total number of first binary component removed from the 
*                  system after SN natal kick
*       NTSN2   -  total number of second binary component removed from the 
*                  system after SN natal kick
*       NTSNB   -  total number of binaries removed from the system after SN
*                  natal kick
*
*          OTHER COMMONS
*
*       X       -  3D positions of stars at T=0
*       XDOT    -  3D velocities of stars at T=0
*       RLAG    -  lagrangian radii
*       V2RL    -  radial velocities ^2 in lagrangian shells
*       V2TL    -  tangential velocities ^2 in lagrangian shells
*       ANI     -  anisotropy in lagrangian shells A=2-v2tl/v2rl
*       RN1     -  new position of the first interacting star
*       RN2     -  new position of the second interacting star
*       VRN1    -  new radial velocity of the first interacting star
*       VRT1    -  new tangential velocity of the first interacting star
*       VRN2    -  the same as VRN1 but for second star
*       VRT2    -  the same as VTN1 but for second star
*       RTIDKG  -  tidal radius from the king model
*       ESCENE  -  energy of escapers
*       ESCDIS  -  place of escape event
*       ESCTIM  -  time of escape event
*       ESCMAS  -  mass of escapers
*       ESCANG  -  angular momentum of escapers
*       ENRAD   -  energy error due to problems with the new velocity
*       SMTO    -  total mass of the old model
*       SINFI   -  sin(fi)
*       COSFI   -  cos(fi)
*       FI      -  angle between tangential velocities of interacting stars
*       XMIN    -  pericentre distance distances for each particle
*       XMAX    -  apocentre distance distances for each particle
*       XMAXZ   -  the maximum radius for paricles in the curent zone
*       XGMIN   -  minimum value of q - determination of a new position
*       XGMAX   -  maximum value of q - determination of a new position
*       TTP     -  time for profile output (in 10^6 years)
*       TTE     -  time for mloss call for all objects (in 10^6 years)
*       VRP     -  radial velocity computed in newpos and used in timepot
*       VTP     -  tangential velocity computed in newpos and used in timepot
*
*       IENAME  -  names of escaping stars in each cycle
*       IESC    -  number of escaping stars in each cycle
*       IENAM1  -  names of escaping stars in each model
*       IEKIND  -  type of escaping stars
*       IESC1   -  number of escaping stars in each model
*       IVNEWG  -  number of events when 'vnew' is < 0
*       IVRR    -  total number of relaxed stars from time = 0
*       NTO     -  total number of star in the old model
*       IDUM2   -  parameter for random number function
*       IY      -  parameter for random number function
*       IV      -  parameters for random number function  (32)
*
*       ----------------------------------------------------------------------
*
*
      COMMON /PARAM/ BMIN,BMAX,TAU0,TCRIT,TCOMP,QE,ALPHAL,ALPHAH,
     &               BRAKEM,BODY1,BODYM,BODYN,QVIR,RBAR,ZMBAR,W0,
     &               TSCALE,TTP,FRACB,AMIN,AMAX,RPLUM,DTTP,DTTE,TTE,
     &               DTTE0,TCREVO,XTAU,TSCALE0,TIMET,XTID,TOLM
*
      COMMON /COEFI/ ONE2,ONE3,ONE4,ONE5,ONE6,ONE9,ONE12,PI,TWOPI,
     &               RSUNTOPC,FLAGR(NLAGRA)
*
      COMMON /BODY/ BODY(NMAX),VR(NMAX),VT(NMAX),U(NMAX),
     &              XESCP(NMAX),XESCT(NMAX),VRR(NMAX),R(NMAX)
*
      COMMON /BODYO/ RO(NMAX),VRO(NMAX),VTO(NMAX),UO(NMAX),RN(NMAX)
*
*      COMMON /BINAR/ BIN2(NBMAX2,7),BIN3(NBMAX3,7),BINP(NBMAXP,7),
*     &               BIN4(NBMAX4,15)
      COMMON /BINAR/ BIN(NBMAX3,8)
*
*      COMMON /BININF/ BIN2IN(1000*NBMAX2,6),BIN3IN(1000*NBMAX3,6),
*     &                BIN4IN(1000*NBMAX4,6),BINPIN(1000*NBMAXP,6)
      COMMON /BININF/ BININ(50*NBMAX3,7)
*
      COMMON /SYSTEM/ SMT,ETOT,ESCSTA,ESCBI2,ESCBI3,ESCBI4,ESCBIP,
     &                ESCBB2,ESCBB3,ESCBB4,ESCBBP,EHMLEV,SLOSES,
     &                ESCB2S,ESCB3S,ESCB4S,ESCBPS,EHBIN2,EHB3B3,
     &                EHBIN3,EHBIN4,EHBINP,ECBIN2,ECBIN3,ECBIN4,
     &                ECBINP,SLOSEV,SLOSCO,ECCOLL,EHCOLL,ZKIN,POT,
     &                RSCALE,TCR,TIME,RTID,TAU,ERROR,CPU,ENRAD,
     &                RMGLOB,GAMMA,SLOB3S,SLOB3B,ERELB3,ERB3IN,
     &                ENEPOT,EHBI3P,BMIN0,YTAU,YBMIN,EKICKT,ENEKIN,
     &                SESCRT,EKICKTBS,EKICKTBD,EKICKTB2
*
      COMMON /PROBA/  PBIN2,PBIN3,PBIN4,PBIN2S,PBIN3S,PBIN4S,PBINPS,
     &                PB2B2,PB2B3,PB2B4,PB2BP,PB3B3,PB3B4,PB3BP,
     &                PB4B4,PB4BP,PBPBP,PCOLL
*
      COMMON /COREP/  SMC,RC,ROC,VC,NC
*
      COMMON /IPARAM/ IRUN,NT0,NT,ISTART,NCOR,NMIN,NZ0,NZONC,NMINZO,
     &                NTWO,IPRINT,IMODEL,ISEED,IB3F,NSS0,NT00,IEXCH,
     &                NTNEW,IKROUPA,NITESC
*
      COMMON /IBODY/ NZS(NZONMA),NZST(NSUPZO),LTWO(NSUPZO),IKIND(NMAX),
     &               NAMES(NMAX),NAMEB(NMAX),NWHICH(NMAX),IBSTRA(NMAX),
     &               IINTE3(NBMAX3),NBINAR(NMAX),NKICK(NMAX),IZON,
     &               NSUZON,NZON0

*
      COMMON /ISYSTE/ NESCST,NESCB2,NESCB3,NESCB4,NESCBP,NESB2S,
     &                NESB3S,NESB4S,NESBPS,NMLOEV,NMLOCO,NCOLL,NECB2,
     &                NECB3,NECB4,NECBP,NEHB2,NEHB3,NEHB4,NEHBP,NBIN2,
     &                NBIN3,NBIN4,NBINP,IOBT,INAME(NMAX),INAMEO(NMAX),
     &                IVNEWG,IVRR,NB3FIN,NB3B3,NDIST3,NDIST4,NDISTE,
     &                IBSM,IBSC,IBS3,IBS4,NMERGE,NESCM,IKICKT,MNSBH,
     &                NEXCHANG,INEXCH(NMAX),NEXCHANG2,INEXCH2(NMAX),
     &                NESCRT,IKICKTBS,IKICKTBD,NTSN1,NTSN2,NTSNB
*
*       OTHER COMMONS
*
*
      COMMON /POSVEL/ X(NMAX,3),XDOT(NMAX,3)
*
      COMMON /LAGRAN/ RLAG(NLAGRA),V2RL(NLAGRA),V2TL(NLAGRA),
     &                ANI(NLAGRA)
*
      COMMON /KINGM/  RTIDKG
*
      COMMON /UPTIME/ UPTIME(NMAX),OLDTIME(NMAX)
*    
      COMMON /TWOPOS/ RN1,RN2
*
      COMMON /TWOVEL/ VRN1,VRN2,VTN1,VTN2
*
      COMMON /ESCAP/  ESCENE(100000),ESCANG(100000),ESCDIS(100000),
     &                ESCMAS(100000),ESCTIM(100000),IENAME(100000),
     &                IENAM1(100000),IEKIND(100000),IESC,IESC1
*
      COMMON /OLDPOT/ SMTO,NTO
*
      COMMON /ANGLEFI/ SINFI,COSFI,FI
*
      COMMON /INTEGRAL/ XMIN(NMAX),XMAX(NMAX),XMAXZ(NMAX),
     &                  XGMIN(NMAX),XGMAX(NMAX)
*
      COMMON /RANDX/ IDUM2,IY,IV(32)
*
      COMMON /TIMEP/ VRP(NMAX),VTP(NMAX)
*
      COMMON /ZSET/ ZINI
*
      COMMON /FFLAGS/ IFLAGNS,IFLAGBH
*
      COMMON /RUNTIME/ TIMEOLD
*
      COMMON /IRUNTIME/ ITIME            
*
*
*
*
     COMMON /AMUSE/ datadir
