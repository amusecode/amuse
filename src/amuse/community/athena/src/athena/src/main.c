#include "copyright.h"
#define MAIN_C
/*============================================================================*/
/*! \file main.c
 *  \brief Athena main program file.
 *
 * //////////////////////////// ATHENA Main Program \\\\\\\\\\\\\\\\\\\\\\\\\\\
 *
 *  Athena - C version developed by JM Stone, TA Gardiner, & PJ Teuben.
 *
 *  Significant additional contributions from X. Bai, K. Beckwith, N. Lemaster,
 *  I. Parrish, & A. Skinner.  See also the F90 version developed by JF Hawley
 *  & JB Simon.
 *
 *  History:
 * - v1.0 [Feb 2003] - 1D adiabatic and isothermal MHD
 * - v1.1 [Sep 2003] - bug fixes in eigensystems
 * - v2.0 [Dec 2004] - 2D adiabatic and isothermal MHD
 * - v3.0 [Feb 2007] - 3D adiabatic and isothermal MHD with MPI
 * - v3.1 [Jan 2008] - multiple species, self-gravity
 * - v3.2 [Sep 2009] - viscosity, resistivity, conduction, particles, special
 *                     relativity, cylindrical coordinates
 * - v4.0 [Jul 2010] - static mesh refinement with MPI
 *
 * See the GNU General Public License for usage restrictions. 
 *									        
 * PRIVATE FUNCTION PROTOTYPES:
 * - change_rundir() - creates and outputs data to new directory
 * - usage()         - outputs help message and terminates execution	      */
/*============================================================================*/
static char *athena_version = "version 4.0 - 01-Jul-2010";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   change_rundir - creates and outputs data to new directory
 *   usage         - outputs help message and terminates execution
 *============================================================================*/

static void change_rundir(const char *name);
static void usage(const char *prog);

/* Maximum number of mkdir() and chdir() file operations that will be executed
 * at once in the change_rundir() function when running in parallel, passed to
 * baton_start() and baton_end().
 */ 
#define MAX_FILE_OP 256

/*----------------------------------------------------------------------------*/
/*! \fn int main(int argc, char *argv[]) 
 *  \brief Athena main program  
 *
 * Steps in main:
 * - 1 - check for command line options and respond
 * - 2 - read input file and parse command line for changes
 * - 3 - set up diagnostic log files
 * - 4 - initialize Mesh, Domains, and Grids
 * - 5 - set initial conditions
 * - 6 - set boundary condition function pointers, and use to set BCs
 * - 7 - set function pointers for desired algorithms and physics
 * - 8 - write initial conditions to output file(s)
 * - 9 - main loop
 * - 10 - finish by writing final output(s), diagnostics, and free memory     */

int main(int argc, char *argv[])
{
  MeshS Mesh;             /* the entire mesh hierarchy, see athena.h */
  VDFun_t Integrate;      /* function pointer to integrator, set at runtime */
#ifdef SELF_GRAVITY
  VDFun_t SelfGrav;      /* function pointer to self-gravity, set at runtime */
#endif
  int nl,nd;
  char *definput = "athinput";  /* default input filename */
  char *athinput = definput;
  int ires=0;             /* restart flag, set to 1 if -r argument on cmdline */
  char *res_file = NULL;  /* restart filename */
  char *rundir = NULL;    /* directory to which output is directed */
  int nstep_start=0;      /* number of cycles already completed on restart */
  char *name = NULL;
  char level_dir[80];     /* names of directories for levels > 0 */
  FILE *fp;               /* file pointer for data outputs */
  int nflag=0;            /* set to 1 if -n argument is given on command line */
  int i,nlim;             /* cycle index and limit */
  Real tlim;              /* time limit (in code units) */

  int out_level, err_level, lazy; /* diagnostic output & error log levels */
  int iflush, nflush;             /* flush buffers every iflush cycles */

  int iquit=0;  /* quit signal sent to ath_sig_act, our system signal handler */

/* local variables used for timing and performance measures */

  time_t start, stop;
  int have_time = time(&start);  /* Is current calendar time (UTC) available? */
  int zones;
  double cpu_time, zcs;
  long clk_tck = sysconf(_SC_CLK_TCK);
  struct tms tbuf;
  clock_t time0,time1, have_times;
  struct timeval tvs, tve;
  Real dt_done;

#ifdef MPI_PARALLEL
  char *pc, *suffix, new_name[MAXLEN];
  int len, h, m, s, err, use_wtlim=0;
  double wtend;
  if(MPI_SUCCESS != MPI_Init(&argc, &argv))
    ath_error("[main]: Error on calling MPI_Init\n");
#endif /* MPI_PARALLEL */

/*----------------------------------------------------------------------------*/
/* Steps in main:
 *  1 - check for command line options and respond
 *  2 - read input file and parse command line for changes
 *  3 - set up diagnostic log files
 *  4 - initialize Mesh, Domains, and Grids
 *  5 - set initial conditions
 *  6 - set boundary condition function pointers, and use to set BCs
 *  7 - set function pointers for desired algorithms and physics
 *  8 - write initial conditions to output file(s)
 *  9 - main loop
 *  10 - finish by writing final output(s), diagnostics, and free memory
 */

/*--- Step 1. ----------------------------------------------------------------*/
/* Check for command line options and respond.  See comments in usage()
 * for description of options.  */

  for (i=1; i<argc; i++) {
/* If argv[i] is a 2 character string of the form "-?" then: */
    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'i':                      /* -i <file>   */
	athinput = argv[++i];
	break;
      case 'r':                      /* -r <file>   */
	ires = 1;
	res_file = argv[++i];
/* If input file is not set on command line, use the restart file */
	if(athinput == definput) athinput = res_file;
	break;
      case 'd':                      /* -d <directory>   */
	rundir = argv[++i];
	break;
      case 'n':                      /* -n */
	nflag = 1;
	break;
      case 'h':                      /* -h */
	usage(argv[0]);
	break;
      case 'c':                      /* -c */
	show_config();
	exit(0);
	break;
#ifdef MPI_PARALLEL
      case 't':                      /* -t hh:mm:ss */
	use_wtlim = 1; /* Logical to use a wall time limit */
	sscanf(argv[++i],"%d:%d:%d",&h,&m,&s);
	wtend = MPI_Wtime() + s + 60*(m + 60*h);
	printf("Wall time limit: %d hrs, %d min, %d sec\n",h,m,s);
	break;
#else
      default:
	usage(argv[0]);
	break;
#endif /* MPI_PARALLEL */
      }
    }
  }

/*--- Step 2. ----------------------------------------------------------------*/
/* Read input file and parse command line.  For MPI_PARALLEL jobs, parent reads
 * input file and distributes information to children  */

#ifdef MPI_PARALLEL
/* Get proc id (rank) in MPI_COMM_WORLD, store as global variable */

  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &myID_Comm_world))
    ath_error("[main]: Error on calling MPI_Comm_rank\n");

/* Only rank=0 processor reads input parameter file, parses command line,
 * broadcasts the contents of the (updated) parameter file to the children. */

  if(myID_Comm_world == 0){
    par_open(athinput);   /* for restarts, default is athinput=resfile */ 
    par_cmdline(argc,argv);
  }
  par_dist_mpi(myID_Comm_world,MPI_COMM_WORLD);

/* Modify the problem_id name in the <job> block to include information about
 * processor ids, so that all output filenames constructed from this name will
 * include myID_Comm_world as an identifier.  Only child processes modify
 * name, rank=0 process does not have myID_Comm_world in the filename */

  if(myID_Comm_world != 0){
    name = par_gets("job","problem_id");
    sprintf(new_name,"%s-id%d",name,myID_Comm_world);
    free(name);
    par_sets("job","problem_id",new_name,NULL);
  }

  show_config_par(); /* Add the configure block to the parameter database */

/* Share the restart flag with the children */

  if(MPI_SUCCESS != MPI_Bcast(&ires, 1, MPI_INT, 0, MPI_COMM_WORLD))
    ath_error("[main]: Error on calling MPI_Bcast\n");

/* rank=0 needs to send the restart file name to the children.  This requires 
 * sending the length of the restart filename string, the string, and then
 * having each child add my_id to the name so it opens the appropriate file */

/* Parent finds length of restart filename */

  if(ires){ 
    if(myID_Comm_world == 0)
      len = 1 + (int)strlen(res_file);

/* Share this length with the children */

    if(MPI_SUCCESS != MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD))
      ath_error("[main]: Error on calling MPI_Bcast\n");

    if(len + 10 > MAXLEN)
      ath_error("[main]: Restart filename length = %d is too large\n",len);

/* Share the restart filename with the children */

    if(myID_Comm_world == 0) strcpy(new_name, res_file);
    if(MPI_SUCCESS != MPI_Bcast(new_name, len, MPI_CHAR, 0, MPI_COMM_WORLD))
      ath_error("[main]: Error on calling MPI_Bcast\n");

/* Assume the restart file name is of the form
 *  [/some/dir/]basename.0000.rst and search for the periods in the name. */

    pc = &(new_name[len - 5]);
    if(*pc != '.'){
      ath_error("[main]: Bad Restart filename: %s\n",new_name);
    }

    do{ /* Position the char pointer at the first period */
      pc--;
      if(pc == new_name)
	ath_error("[main]: Bad Restart filename: %s\n",new_name);
    }while(*pc != '.');

/* Only children add myID_Comm_world to the filename */

    if(myID_Comm_world == 0) {
      strcpy(new_name, res_file);
    } else {       
      suffix = ath_strdup(pc);
      sprintf(pc,"-id%d%s",myID_Comm_world,suffix);
      free(suffix);
      res_file = new_name;
    }
  }

/* Quit MPI_PARALLEL job if code was run with -n option. */

  if(nflag){          
    par_dump(0,stdout);   
    par_close();
    MPI_Finalize();
    return 0;
  }

#else
/* For serial (single processor) job, there is only one process to open and
 * read input file  */

  myID_Comm_world = 0;
  par_open(athinput);   /* opens AND reads */
  par_cmdline(argc,argv);
  show_config_par();   /* Add the configure block to the parameter database */

/* Quit non-MPI_PARALLEL job if code was run with -n option. */

  if(nflag){
    par_dump(0,stdout);
    par_close();
    return 0;
  }
#endif /* MPI_PARALLEL */

/*--- Step 3. ----------------------------------------------------------------*/
/* set up the simulation log files */

/* Open <problem_id>.out and <problem_id>.err files if file_open=1 in the 
 * <log> block of the input file.  Otherwise, diagnositic output will go to
 * stdout and stderr streams. */

  if(par_geti_def("log","file_open",0)){
    iflush = par_geti_def("log","iflush",0);
    name = par_gets("job","problem_id");
    lazy = par_geti_def("log","lazy",1);
    /* On restart we use mode "a", otherwise we use mode "w". */
    ath_log_open(name, lazy, (ires ? "a" : "w"));
    free(name);  name = NULL;
  }
  else{
    iflush = par_geti_def("log","iflush",1);
  }
  iflush = iflush > 0 ? iflush : 0; /* make iflush non-negative */

/* Set the ath_log output and error logging levels */
  out_level = par_geti_def("log","out_level",0);
  err_level = par_geti_def("log","err_level",0);
#ifdef MPI_PARALLEL
    if(myID_Comm_world > 0){   /* Children may use different log levels */
    out_level = par_geti_def("log","child_out_level",-1);
    err_level = par_geti_def("log","child_err_level",-1);
  }
#endif /* MPI_PARALLEL */
  ath_log_set_level(out_level, err_level);

  if(have_time > 0) /* current calendar time (UTC) is available */
    ath_pout(0,"Simulation started on %s\n",ctime(&start));

/*--- Step 4. ----------------------------------------------------------------*/
/* Initialize nested mesh hierarchy. */

  init_mesh(&Mesh);
  init_grid(&Mesh);
#ifdef PARTICLES
  init_particle(&Mesh);
#endif

/*--- Step 5. ----------------------------------------------------------------*/
/* Set initial conditions, either by reading from restart or calling problem
 * generator.  But first start by setting variables in <time> block (these
 * control execution time), and reading EOS parameters from <problem> block.  */

  CourNo = par_getd("time","cour_no");
  nlim = par_geti_def("time","nlim",-1);
  tlim = par_getd("time","tlim");

#ifdef ISOTHERMAL
  Iso_csound = par_getd("problem","iso_csound");
  Iso_csound2 = Iso_csound*Iso_csound;
#else
  Gamma = par_getd("problem","gamma");
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;
#endif
/* initialize gravity constants <0, selfg_init will test these values below to
 * ensure user has set values in problem generator */
#ifdef SELF_GRAVITY
  grav_mean_rho = -1.0;
  four_pi_G = -1.0;
#endif

  if(ires) {
    restart_grids(res_file, &Mesh);  /*  Restart */
    nstep_start = Mesh.nstep;
  } else {                           /* New problem */
    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL) problem(&(Mesh.Domain[nl][nd]));
      }
    }
  }

/* restrict initial solution so grid hierarchy is consistent */
#ifdef STATIC_MESH_REFINEMENT
  SMR_init(&Mesh);
  RestrictCorrect(&Mesh);
#endif

/* Initialize the first nstep value to flush the output and error logs. */
  nflush = nstep_start + iflush;

/*--- Step 6. ----------------------------------------------------------------*/
/* set boundary value function pointers using BC flags in <grid> blocks, then
 * set boundary conditions for initial conditions.  With SMR, this includes
 * a prolongation step to set ghost zones at internal fine/coarse boundaries  */

  bvals_mhd_init(&Mesh);

#ifdef SELF_GRAVITY
  bvals_grav_init(&Mesh);
#endif
#if defined(SHEARING_BOX) || (defined(FARGO) && defined(CYLINDRICAL))
  bvals_shear_init(&Mesh);
#endif
#ifdef PARTICLES
  bvals_particle_init(&Mesh);
#endif

  for (nl=0; nl<(Mesh.NLevels); nl++){ 
    for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
      if (Mesh.Domain[nl][nd].Grid != NULL){
        bvals_mhd(&(Mesh.Domain[nl][nd]));
#ifdef PARTICLES
        bvals_particle(&(Mesh.Domain[nl][nd]));
#ifdef FEEDBACK
        exchange_feedback_init(&(Mesh.Domain[nl][nd]));
#endif
#endif
      }
    }
  }

/* Now that BC set, prolongate solution into child Grid GZ with SMR */
#ifdef STATIC_MESH_REFINEMENT
  Prolongate(&Mesh);
#endif

/* For new runs, set initial timestep */

  if(ires == 0) new_dt(&Mesh);

/*--- Step 7. ----------------------------------------------------------------*/
/* Set function pointers for integrator; self-gravity (based on dimensions)
 * Initialize gravitational potential for new runs
 * Allocate temporary arrays */

  init_output(&Mesh); 
  lr_states_init(&Mesh);
  Integrate = integrate_init(&Mesh);
#ifdef SELF_GRAVITY
  SelfGrav = selfg_init(&Mesh);
  for (nl=0; nl<(Mesh.NLevels); nl++){ 
    for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
      if (Mesh.Domain[nl][nd].Grid != NULL){
        (*SelfGrav)(&(Mesh.Domain[nl][nd]));
        bvals_grav(&(Mesh.Domain[nl][nd]));
      }
    }
  }
#endif
#if defined(RESISTIVITY) || defined(VISCOSITY) || defined(THERMAL_CONDUCTION)
  integrate_diff_init(&Mesh);
#endif

/*--- Step 8. ----------------------------------------------------------------*/
/* Setup complete, output initial conditions */

  if(out_level >= 0){
    fp = athout_fp();
    par_dump(0,fp);      /* Dump a copy of the parsed information to athout */
  }
  change_rundir(rundir); /* Change to run directory */
  ath_sig_init();        /* Install a signal handler */
  for (nl=1; nl<(Mesh.NLevels); nl++){
    sprintf(level_dir,"lev%d",nl);
    mkdir(level_dir, 0775); /* Create directories for levels > 0 */
  }

  gettimeofday(&tvs,NULL);
  if((have_times = times(&tbuf)) > 0)
    time0 = tbuf.tms_utime + tbuf.tms_stime;
  else
    time0 = clock();

/* Force output of everything (by passing last argument of data_output = 1) */

  if (ires==0) data_output(&Mesh, 1);

  ath_pout(0,"\nSetup complete, entering main loop...\n\n");
  ath_pout(0,"cycle=%i time=%e next dt=%e\n",Mesh.nstep, Mesh.time, Mesh.dt);

/*--- Step 9. ----------------------------------------------------------------*/
/* START OF MAIN INTEGRATION LOOP ==============================================
 * Steps are: (a) Check for data ouput
 *            (b) Add explicit diffusion with operator splitting
 *            (c) Integrate all Grids over Mesh hierarchy
 *            (d) Restrict solution and correct fine/coarse boundaries
 *            (e) Userwork
 *            (f) Self-gravity
 *            (g) Update time, set new timestep
 *            (h) Set boundary values
 *            (i) check for stopping criteria
 */

  while (Mesh.time < tlim && (nlim < 0 || Mesh.nstep < nlim)) {

/*--- Step 9a. ---------------------------------------------------------------*/
/* Only write output's with t_out>t (last argument of data_output = 0) */

    data_output(&Mesh, 0);

/*--- Step 9b. ---------------------------------------------------------------*/
/* operator-split explicit diffusion: thermal conduction, viscosity, resistivity
 * Done first since CFL constraint is applied which may change dt  */

#if defined(RESISTIVITY) || defined(VISCOSITY) || defined(THERMAL_CONDUCTION)
    integrate_diff(&Mesh);
    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL){
          bvals_mhd(&(Mesh.Domain[nl][nd]));
        }
      }
    }
#endif

/*--- Step 9c. ---------------------------------------------------------------*/
/* Loop over all Domains and call Integrator */

    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL){
          (*Integrate)(&(Mesh.Domain[nl][nd]));

#ifdef FARGO
          Fargo(&(Mesh.Domain[nl][nd]));
#ifdef PARTICLES
          advect_particles(&level0_Grid, &level0_Domain);
#endif
#endif /* FARGO */
        }
      }
    }

/*--- Step 9d. ---------------------------------------------------------------*/
/* With SMR, restrict solution from Child --> Parent grids  */

#ifdef STATIC_MESH_REFINEMENT
    RestrictCorrect(&Mesh);
#endif

/*--- Step 9e. ---------------------------------------------------------------*/
/* User work (defined in problem()) */

    Userwork_in_loop(&Mesh);

/*--- Step 9f. ---------------------------------------------------------------*/
/* Compute gravitational potential using new density, and add second-order
 * correction to fluxes for accelerations due to self-gravity. */

#ifdef SELF_GRAVITY
    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL){
          (*SelfGrav)(&(Mesh.Domain[nl][nd]));
          bvals_grav(&(Mesh.Domain[nl][nd]));
          selfg_fc(&(Mesh.Domain[nl][nd]));
        }
      }
    }
#endif

/*--- Step 9g. ---------------------------------------------------------------*/
/* Update Mesh time, and time in all Grid's.  Compute new dt */

    Mesh.nstep++;
    Mesh.time += Mesh.dt;
    for (nl=0; nl<(Mesh.NLevels); nl++){
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){
        if (Mesh.Domain[nl][nd].Grid != NULL){
          Mesh.Domain[nl][nd].Grid->time = Mesh.time;
        }
      }
    }

    dt_done = Mesh.dt;
    new_dt(&Mesh);

/*--- Step 9h. ---------------------------------------------------------------*/
/* Boundary values must be set after time is updated for t-dependent BCs.
 * With SMR, ghost zones at internal fine/coarse boundaries set by Prolongate */

    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL){
          bvals_mhd(&(Mesh.Domain[nl][nd]));
#ifdef PARTICLES
          bvals_particle(&level0_Grid, &level0_Domain);
#endif
        }
      }
    }

#ifdef STATIC_MESH_REFINEMENT
    Prolongate(&Mesh);
#endif

/*--- Step 9i. ---------------------------------------------------------------*/
/* Force quit if wall time limit reached.  Check signals from system */

#ifdef MPI_PARALLEL
    if(use_wtlim && (MPI_Wtime() > wtend))
      iquit = 103; /* an arbitrary, unused signal number */
#endif /* MPI_PARALLEL */
    if(ath_sig_act(&iquit) != 0) break;

/* Print diagnostic message, flush message buffers, and continue... */

    ath_pout(0,"cycle=%i time=%e next dt=%e last dt=%e\n",
	     Mesh.nstep,Mesh.time,Mesh.dt,dt_done);

    if(nflush == Mesh.nstep){
      ath_flush_out();
      ath_flush_err();
      nflush += iflush;
    }
  } /* END OF MAIN INTEGRATION LOOP ==========================================*/

/*--- Step 10. ---------------------------------------------------------------*/
/* Finish up by computing zc/sec, dumping data, and deallocate memory */

/* Print diagnostic message as to why run terminated */

  if (Mesh.nstep == nlim)
    ath_pout(0,"\nterminating on cycle limit\n");
#ifdef MPI_PARALLEL
  else if(use_wtlim && iquit == 103)
    ath_pout(0,"\nterminating on wall-time limit\n");
#endif /* MPI_PARALLEL */
  else
    ath_pout(0,"\nterminating on time limit\n");

/* Get time used */

  gettimeofday(&tve,NULL);
  if(have_times > 0) {
    times(&tbuf);
    time1 = tbuf.tms_utime + tbuf.tms_stime;
    cpu_time = (time1 > time0 ? (double)(time1 - time0) : 1.0)/
      (double)clk_tck;
  } else {
    time1 = clock();
    cpu_time = (time1 > time0 ? (double)(time1 - time0) : 1.0)/
      (double)CLOCKS_PER_SEC;
  }

/* Calculate and print the zone-cycles/cpu-second on this processor */

  zones = 0;
  for (nl=0; nl<(Mesh.NLevels); nl++){ 
  for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
    if (Mesh.Domain[nl][nd].Grid != NULL) {
      zones += (Mesh.Domain[nl][nd].Grid->Nx[0])*
               (Mesh.Domain[nl][nd].Grid->Nx[1])*
               (Mesh.Domain[nl][nd].Grid->Nx[2]);
    }
  }}
  zcs = (double)zones*(double)((Mesh.nstep) - nstep_start)/cpu_time;

  ath_pout(0,"  tlim= %e   nlim= %i\n",tlim,nlim);
  ath_pout(0,"  time= %e  cycle= %i\n",Mesh.time,Mesh.nstep);
  ath_pout(0,"\nzone-cycles/cpu-second = %e\n",zcs);

/* Calculate and print the zone-cycles/wall-second on this processor */

  cpu_time = (double)(tve.tv_sec - tvs.tv_sec) +
    1.0e-6*(double)(tve.tv_usec - tvs.tv_usec);
  zcs = (double)zones*(double)((Mesh.nstep) - nstep_start)/cpu_time;
  ath_pout(0,"\nelapsed wall time = %e sec.\n",cpu_time);
  ath_pout(0,"\nzone-cycles/wall-second = %e\n",zcs);

/* Calculate and print total zone-cycles/wall-second on all processors */
#ifdef MPI_PARALLEL
  zones = 0;
  for (nl=0; nl<(Mesh.NLevels); nl++){ 
  for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
    zones += (Mesh.Domain[nl][nd].Nx[0])*
               (Mesh.Domain[nl][nd].Nx[1])*
               (Mesh.Domain[nl][nd].Nx[2]);
  }}
  zcs = (double)zones*(double)(Mesh.nstep - nstep_start)/cpu_time;
  ath_pout(0,"\ntotal zone-cycles/wall-second = %e\n",zcs);
#endif /* MPI_PARALLEL */

/* complete any final User work */

  Userwork_after_loop(&Mesh);

/* Final output everything (last argument of data_output = 1) */

  data_output(&Mesh, 1);

/* Free all memory */

  lr_states_destruct();
  integrate_destruct();
  data_output_destruct();
#ifdef PARTICLES
  particle_destruct(&level0_Grid);
  bvals_particle_destruct(&level0_Grid);
#endif
#if defined(SHEARING_BOX) || (defined(FARGO) && defined(CYLINDRICAL))
  bvals_shear_destruct();
#endif
#if defined(RESISTIVITY) || defined(VISCOSITY) || defined(THERMAL_CONDUCTION)
  integrate_diff_destruct();
#endif
  par_close();

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif /* MPI_PARALLEL */

  if(time(&stop)>0) /* current calendar time (UTC) is available */
    ath_pout(0,"\nSimulation terminated on %s",ctime(&stop));

  ath_log_close(); /* close the simulation log files */

  return EXIT_SUCCESS;
}

/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*! \fn void change_rundir(const char *name) 
 *  \brief Change run directory;  create it if it does not exist yet
 */

void change_rundir(const char *name)
{
#ifdef MPI_PARALLEL

  int err=0;
  int rerr, gerr, my_id, status;
  char mydir[80];

  status = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  if(status != MPI_SUCCESS)
    ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",status);

  if(name != NULL && *name != '\0'){

    if(my_id == 0)
      mkdir(name, 0775); /* May return an error, e.g. the directory exists */

    MPI_Barrier(MPI_COMM_WORLD); /* Wait for rank 0 to mkdir() */

    baton_start(MAX_FILE_OP, ch_rundir0_tag);

    if(chdir(name)){
      ath_perr(-1,"[change_rundir]: Cannot change directory to \"%s\"\n",name);
      err = 1;
    }

    baton_stop(MAX_FILE_OP, ch_rundir0_tag);

    /* Did anyone fail to make and change to the run directory? */
    rerr = MPI_Allreduce(&err, &gerr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(rerr) ath_perr(-1,"[change_rundir]: MPI_Allreduce error = %d\n",rerr);

    if(rerr || gerr){
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(EXIT_FAILURE);
    }
  }

/* Next, change to the local run directory */

  sprintf(mydir, "id%d", my_id);

  baton_start(MAX_FILE_OP, ch_rundir1_tag);

  mkdir(mydir, 0775); /* May return an error, e.g. the directory exists */
  if(chdir(mydir)){
    ath_perr(-1,"[change_rundir]: Cannot change directory to \"%s\"\n",mydir);
    err = 1;
  }

  baton_stop(MAX_FILE_OP, ch_rundir1_tag);

  /* Did anyone fail to make and change to the local run directory? */
  rerr = MPI_Allreduce(&err, &gerr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(rerr) ath_perr(-1,"[change_rundir]: MPI_Allreduce error = %d\n",rerr);

  if(rerr || gerr){
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(EXIT_FAILURE);
  }

#else /* Serial job */

  if(name == NULL || *name == '\0') return;

  mkdir(name, 0775); /* May return an error, e.g. the directory exists */
  if(chdir(name))
    ath_error("[change_rundir]: Cannot change directory to \"%s\"\n",name);

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void usage(const char *prog)
 *  \brief Outputs help
 *
 *    athena_version is hardwired at beginning of this file
 *    CONFIGURE_DATE is macro set when configure script runs
 */

static void usage(const char *prog)
{
  ath_perr(-1,"Athena %s\n",athena_version);
  ath_perr(-1,"  Last configure: %s\n",CONFIGURE_DATE);
  ath_perr(-1,"\nUsage: %s [options] [block/par=value ...]\n",prog);
  ath_perr(-1,"\nOptions:\n");
  ath_perr(-1,"  -i <file>       Alternate input file [athinput]\n");
  ath_perr(-1,"  -r <file>       Restart a simulation with this file\n");
  ath_perr(-1,"  -d <directory>  Alternate run dir [current dir]\n");
  ath_perr(-1,"  -h              This Help, and configuration settings\n");
  ath_perr(-1,"  -n              Parse input, but don't run program\n");
  ath_perr(-1,"  -c              Show Configuration details and quit\n");
  ath_perr(-1,"  -t hh:mm:ss     With MPI, wall time limit for final output\n");
  show_config();
  exit(0);
}
