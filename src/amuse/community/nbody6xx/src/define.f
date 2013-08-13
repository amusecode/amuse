      SUBROUTINE DEFINE
*
*
*       Definition of input parameters, options & counters.
*       ---------------------------------------------------
*
*
*       Input parameters
*       ****************
*
*       ---------------------------------------------------------------------
***
* NBODY6:
*
*       KSTART  Control index (1: new run; >1: restart; 3, 4, 5: new params).
*       TCOMP   Maximum CPU time in minutes (saved in CPU).
*       TCRITp  Termination time in Myrs.
*       isernb  Max size of sequential irr blocks on parallel machine 
*               for single CPU dummy
*       iserreg as isernb for reg blocks
*               for single CPU dummy
***
* INPUT:
*
*       N       Number of objects (N_s + 2*N_b; singles + 3*NBIN0 < NMAX).
*       NFIX    Output frequency of data save or binaries (options 3 & 6).
*       NCRIT   Final particle number (alternative termination criterion).
*       NRAND   Random number sequence skip.
*       NNBMAX  Maximum number of neighbours (< LMAX - 5).
*       NRUN    Run identification index.
*
*       ETAI    Time-step parameter for irregular force polynomial.
*       ETAR    Time-step parameter for regular force polynomial.
*       RS0     Initial radius of neighbour sphere (N-body units).
*       DTADJ   Time interval for parameter adjustment (N-body units).
*       DELTAT  Output time interval (N-body units).
*       TCRIT   Termination time (N-body units).
*       QE      Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1).
*       RBAR    Virial cluster radius in pc (set = 1 for isolated cluster).
*       ZMBAR   Mean mass in solar units (=1.0 if 0; final depends on #20).
*
*       KZ(J)   Non-zero options for alternative paths (see table).
*
*       DTMIN   Time-step criterion for regularization search.
*       RMIN    Distance criterion for regularization search.
*       ETAU    Regularized time-step parameter (6.28/ETAU steps/orbit).
*       ECLOSE  Binding energy per unit mass for hard binary (positive).
*       GMIN    Relative two-body perturbation for unperturbed motion.
*       GMAX    Secondary termination parameter for soft KS binaries.
***
* INPUT: if (kz(4).gt.0)
*
*       DELTAS  Output interval for binary search (in TCR; suppressed).
*       ORBITS  Minimum periods for binary output (level 1).
*       GPRINT  Perturbation thresholds for binary output (9 levels).
***
* DATA:
*
*       ALPHAS  Power-law index for initial mass function (used if #20 < 2).
*       BODY1   Maximum particle mass before scaling (KZ(20): solar mass).
*       BODYN   Minimum particle mass before scaling.
*       NBIN0   Number of primordial binaries (for IMF2 with KZ(20) > 1).
*       NHI0    Primordial hierarchies (may be needed in IMF if > 0).
*       ZMET    Metal abundance (in range 0.03 - 0.0001).
*       EPOCH0  Evolutionary epoch (in 10**6 yrs; NB! < 0 for PM evolution).
*       DTPLOT  Plotting interval for HRDIAG (N-body units; >= DELTAT).
***
* SETUP: if (kz(5).eq.2)
*
*       APO     Separation of two Plummer models (SEMI = APO/(1 + ECC).
*       ECC     Eccentricity of two-body orbit (ECC < 0.999).
*       N2      Membership of second Plummer model (N2 <= N).
*       SCALE   Second scale factor (>= 0.2 for limiting minimum size).
*
*        if (kz(5).eq.3)
*
*       APO     Separation between the perturber and Sun.
*       ECC     Eccentricity of orbit (=1 for parabolic encounter).
*       DMIN    Minimum distance of approach (pericentre).
*       SCALE   Perturber mass scale factor (=1 for Msun).
*
*        if (kz(5).eq.4)
*
*       SEMI    Semi-major axis (slightly modified; ignore if ECC > 1).
*       ECC     Eccentricity (ECC > 1: NAME = 1 & 2 free-floating).
*       M1      Mass of first member (in units of mean mass).
*       M2      Mass of second member (rescaled total mass = 1).
*
*        if (kz(5).ge.6) & (kz(24).lt.0)
*
*       ZMH     Mass of single BH (in N-body units).
*       RCUT    Radial cutoff in Zhao cusp distribution (MNRAS, 278, 488).
***
* SCALE:
*
*       Q       Virial ratio (Q = 0.5 for equilibrium).
*       VXROT   XY-velocity scaling factor (> 0 for solid-body rotation).
*       VZROT   Z-velocity scaling factor (not used if VXROT = 0).
*       RTIDE   Unscaled tidal radius (#14 >= 2; otherwise copied to RSPH2).
***
* XTRNL0: if (kz(14).eq.2)
*
*       GMG     Point-mass galaxy (solar masses, linearized circular orbit).
*       RG0     Central distance (in kpc).
*
*         if (kz(14).eq.3)
*
*       GMG     Point-mass galaxy (solar masses).
*       DISK    Mass of Miyamoto disk (solar masses).
*       A       Softening length in Miyamoto potential (in kpc).
*       B       Vertical softening length (kpc).
*       VCIRC   Galactic circular velocity (km/sec) at RCIRC (=0: no halo).
*       RCIRC   Central distance for VCIRC with logarithmic potential (kpc).
*       RG      Initial position; DISK+VCIRC=0, VG(3)=0: A(1+E)=RG(1), E=RG(2).
*       VG      Initial cluster velocity vector (km/sec).
*
*         if (kz(14).eq.3.or.kz(14).eq.4)
*
*       MP      Total mass of Plummer sphere (in scaled units).
*       AP      Plummer scale factor (N-body units; square saved in AP2).
*       MPDOT   Decay time for gas expulsion (MP = MP0/(1 + MPDOT*(T-TD)).
*       TDELAY  Delay time for starting gas expulsion (T > TDELAY).
***
* BINPOP: if (kz(8).eq.1.or.kz(8).gt.2) (in nb6++ PK binpop used, RSp)
*         note that NBIN, Number of initial binaries = NBIN0, not read here.
*--------
*       SEMI    Initial semi-major axis (= 0 for range of energies).
*       ECC     Initial eccentricity (for BINPOP_4NEW)
*               <=1 AND >=0 for one particular fixed ecc. for all systems
*               < 0 for thermal distribution,
*               =20 for uniform distribution,
*               =30 for f(e)=0.1765/(e*e)
*               =40 for general f(e)=a*e^b, e0<=e<=1 with a=(1+b)/(1-e0^(1+b))
*                   e0 and b must be defined in binpop routine
*       RATIO   Mass ratio M1/(M1 + M2); (= 1.0: M1 = M2 = <M>).
*       NBGR    Number of binaries in fixed energy groups.
*       REDUCE  Reduction factor in semi-major axis for each group.
*       RANGE   Energy range for uniform logarithmic distribution.
*       NSKIP   Binary frequency of mass spectrum (starting from body #1).
*       IDORM   Indicator for dormant binaries (>0: merged components).
*       ICIRC   Eigenevolution & period distribution (RANGE: minimum period).
***
* HIPOP: if (kz(8).gt.0.and.kz(18).gt.1)
*
*       SEMI    Max semi-major axis in model units (all equal if RANGE = 0).
*       ECC     Initial eccentricity (< 0 for thermal distribution).
*       RATIO   Mass ratio (= 1.0: M1 = M2; random in [0.5-0.9]).
*       RANGE   Range in SEMI for uniform logarithmic distribution (> 0).
*       ICIRC   Circularization & collision check (not implemented yet).
***
* INTIDE: if (kz(27).gt.0)  (currently suppressed)
*
*       RSTAR   Size of typical star in A.U.
*       IMS     # idealized main-sequence stars.
*       IEV     # idealized evolved stars.
*       RMS     Scale factor for main-sequence radii (>0: fudge factor).
*       REV     Scale factor for evolved radii (initial size RSTAR).
***
* CLOUD0: if (kz(13).gt.0)
*
*       NCL     Number of interstellar clouds.
*       RB2     Radius of cloud boundary in pc (square is saved).
*       VCL     Mean cloud velocity in km/sec.
*       SIGMA   Velocity dispersion (#13 > 1: Gaussian).
*       CLM     Individual cloud masses in solar masses (maximum MCL).
*       RCL2    Half-mass radii of clouds in pc (square is saved).
*       ---------------------------------------------------------------------
*
*
*       Options KZ(J)
*       *************
*
*       ---------------------------------------------------------------------
*       1  COMMON save unit 1 (=1: 'touch STOP'; =2: every 100*NMAX steps).
*       2  COMMON save unit 2 (=1: at output; =2: restart if DE/E > 5*QE).
*       3  Basic data unit 3 at output time (unformatted, frequency NFIX;
*             =1/2: standard /and tail; =3: tail only; >3: cluster + tail).
*       4  Binary diagnostics on unit 4 (# threshold levels = KZ(4) < 10).
*                                       (currently suppressed in ksint.f.)
*       5  Initial conditions (#22 =0; =0: uniform & isotropic sphere);
*                =1: Plummer; =2: two Plummer models in orbit, extra input;
*                =3: massive perturber and planetesimal disk, extra input;
*                =4: massive initial binary, extra input: A, E, M1, M2;
*                =5: Jaffe model;
*               >=6: Zhao BH cusp model, extra input if #24 < 0: ZMH, RCUT.
*       6  Significant & regularized binaries at main output (=1, 2, 3 & 4).
*       7  Lagrangian radii (>0: RSCALE; =2, 3, 4: output units 6, 7, 12);
*                >=2: half-mass radii of 50% mass, also 1% heavies, unit 6;
*                >=2: Lagrangian radii for two mass groups on unit 31 & 32;
*                >=2: harmonic radii for three mass groups on unit 6;
*                 =5: density, rms velocity & mean mass on unit 26, 27 & 36;
*                 =6: pairwise values of mean mass and radii on unit 28.
*       8  Primordial binaries (=1 & >=3; >0: BINOUT; >2: BINDAT; >3: HIDAT;
*                               =4: Kroupa 1995 period distribution;
*                               >4: standard setup using RANGE & SEMI0).
*       9  Individual bodies on unit 6 at main output (MIN(5**KZ9,NTOT)).
*      10  Diagnostic KS output (>0: begin KS; >1: end; >=3: each step).
*      11  Not used in NBODY6++.
*      12  HR diagnostics of evolving stars (interval DTPLOT).
*      13  Interstellar clouds (=1: constant velocity; >1: Gaussian).
*      14  External force (=1: standard tidal field; =2: point-mass galaxy;
*               =3: point-mass + disk + halo + Plummer; =4: Plummer model).
*      15  Triple, quad, chain (#30 > 0) or merger search (>1: full output).
*      16  Updating of regularization parameters (>0: RMIN, DTMIN & ECLOSE);
*                  >1: RMIN expression based on core radius (experimental);
*                  >2: modify RMIN for GPERT > 0.05 or < 0.002 in chain.
*      17  Modification of ETAI, ETAR (>=1) and ETAU (>1) by tolerance QE.
*      18  Hierarchical systems (=1: diagnostics; =2: primordial; =3: both).
*      19  Mass loss (=1: old supernova scheme; =3: Eggleton, Tout & Hurley;
*                                               >3: extra diagnostics).
*      20  Initial mass function (=0: Salpeter type using ALPHAS; =1: Scalo;
*              =2, 4, 6: Kroupa; =3, 5: Eggleton; > 1: primordial binaries).
*      21  Extra output (>0: MODEL #, TCOMP, DMIN, AMIN; >1: NESC by JACOBI).
*      22  Initial m,r,v on #10 (=1: output; >=2: input; >2 & < 6: no scaling).
*                           (=4&8: starlab, 3&7 tree, 2&6 nbody input format)
*      23  Escaper removal (>1: diagnostics in file ESC with V_inf in km/s);
*                           >=3: initialization & integration of tidal tail.
*      24  Initial conditions for subsystem (M,X,V routine SCALE; KZ(24)= #);
*                           <0: ZMH & RCUT (N-body units) Zhao model (#5>=6).
*      25  Velocity kicks for white dwarfs (=1: type 11 & 12; >1: all WDs).
*      25  Partial reflection of KS binary orbit (GAMMA < GMIN; suppressed).
*      26  Slow-down of two-body motion (>=1: KS; >=2: chain; =3: rectify).
*      27  Tidal effects (=1: sequential; =2: chaos; =3: GR energy loss);
*                         =-1: collision detector, no coalescence, #13 < 0.
*      28  GR radiation for NS & BH binaries (with #19 = 3; choice of #27);
*                         =3: and #5 >= 6: input of ZMH = 1/SQRT(2*N).
*      29  Boundary reflection for hot system (suppressed).
*      30  Multiple regularization (=1: all; >1: BEGIN/END; >2: each step);
*                                =-1: CHAIN only; =-2: TRIPLE & QUAD only. 
*      31  Centre of mass correction after energy check.
*      32  Increase of output intervals (based on single particle energy).
*      33  Histograms at main output (>=1: STEP; =2: STEPR, NBHIST & BINARY).
*      34  Roche-lobe overflow (=1: Roche & Synch; =2: Roche & BSE synch).
*      35  Time offset (global time from TTOT = TIME + TOFF; offset = 100).
*      36  Step reduction for hierarchical systems (suppressed).
*      37  Neighbour additions in CHECKL (>0: high-velocity; >1: all types).
*      38  Force polynomial corrections (=0: standard, no corrections;
*                                =1: all gains & losses included;
*                                =2: small FREG change skipped;
*                                =3: fast neighbour loss only).
*      39  No unique density centre (skips velocity modification of RS(I)).
*      40  Neighbour number control (=1: increase if <NNB>  <  NNBMAX/2);
*                     >=2: fine-tuning at NNBMAX/5; =3: reduction of NNBMAX.
*      41-50  Currently free.
*       ---------------------------------------------------------------------
*
* NBODY6: Restart from fort.1
*
*       KSTART TCOMP (KSTART = 2, 3, 4, 5)
*
*       DTADJ DELTAT TADJ TNEXT TCRIT QE J KZ(J) (if > 0 & KSTART = 3 or 5).
*       
*       ETAI ETAR ETAU DTMIN RMIN NCRIT (if > 0 & KSTART = 4 or 5).
*
*       ---------------------------------------------------------------------
*
*       Output counters
*       ***************
*
*       ---------------------------------------------------------------------
*       NSTEPI  Irregular integration steps.
*       NSTEPR  Regular integration steps.
*       NSTEPU  Regularized integration steps.
*       NNPRED  Coordinate & velocity predictions of all particles.
*       NBPRED  Coordinate & velocity prediction of neighbours (NNB counted).
*       NBCORR  Force polynomial corrections.
*       NBFULL  Too many neighbours with standard criterion.
*       NBVOID  No neighbours inside 1.26 times the basic sphere radius.
*       NICONV  Irregular step reduction (force convergence test).
*       NBSMIN  Retained neighbours inside 2*RS (STEP < SMIN).
*       NLSMIN  Small step neighbours selected from other neighbour lists.
*       NBDIS   Second component of recent KS pair added as neighbour (#18).
*       NBDIS2  Second component of old KS pair added as neighbour (#18 > 1).
*       NCMDER  C.m. values for force derivatives of KS component.
*       NBDER   Large F3DOT corrections not included in D3 & D3R.
*       NFAST   Fast particles included in LISTV (option 18).
*       NBFAST  Fast particles included in neighbour list (option 18).
*       NBLOCK  Number of blocks (block-step version).
*       NMDOT   Mass loss events (option 19).
*       NBSTAT  Diagnostic data on binary interactions (option 4).
*       NKSTRY  Two-body regularization attempts.
*       NKSREG  Total KS regularizations.
*       NEWKS   Enforced KS regularization using wider criterion (~8 > 0).
*       NKSHYP  Hyperbolic KS regularizations.
*       NKSPER  Unperturbed KS binary orbits.
*       NPRECT  Initialization of NKSPER after exceeding 2*10**9.
*       NKSREF  Partial reflections of KS binary (option 25; suppressed).
*       NKSMOD  Slow KS motion restarts (option 26).
*       NTTRY   Search for triple, quad & chain regularization or mergers.
*       NTRIP   Three-body regularizations (option 15).
*       NQUAD   Four-body regularizations (option 15).
*       NCHAIN  Chain regularizations (options 15 & 30).
*       NMERG   Mergers of stable triples or quadruples (option 15).
*       NEWHI   New hierarchical systems (counted by routine HIARCH).
*       NSTEPT  Triple regularization integration steps (option 15).
*       NSTEPQ  Quadruple regularization integration steps (option 15).
*       NSTEPC  Chain regularization steps (# DIFSY calls).
*       NDISS   Tidal dissipations at pericentre (option 27).
*       NTIDE   Tidal captures from hyperbolic motion (option 27).
*       NSYNC   Number of synchronous binaries (option 27).
*       NCOLL   Stellar collisions (option 27).
*       NSESC   Escaped single particles (option 23).
*       NBESC   Escaped binaries (option 23).
*       NMESC   Escaped mergers (options 15 & 23).
*       NRG     Red giants.
*       NHE     Helium stars.
*       NRS     Red supergiants.
*       NNH     Naked Helium stars.
*       NWD     White dwarfs.
*       NSN     Neutron stars.
*       NBH     Black holes.
*       NBS     Blue stragglers.
*       ---------------------------------------------------------------------
*
*
*       Stellar evolution types
*       ***********************
*
*       ---------------------------------------------------------------------
*       0       Low main sequence (M < 0.7).
*       1       Main sequence.
*       2       Hertzsprung gap (HG).
*       3       Red giant.
*       4       Core Helium burning.
*       5       First AGB.
*       6       Second AGB.
*       7       Helium main sequence.
*       8       Helium HG.
*       9       Helium GB.
*      10       Helium white dwarf.
*      11       Carbon-Oxygen white dwarf.
*      12       Oxygen-Neon white dwarf.
*      13       Neutron star.
*      14       Black hole.
*      15       Massless supernova remnant.
*      19       Circularizing binary (c.m. value).
*      20       Circularized binary.
*      21       First Roche stage (inactive).
*      22       Second Roche stage.
*       ---------------------------------------------------------------------
*
      RETURN
*
      END
