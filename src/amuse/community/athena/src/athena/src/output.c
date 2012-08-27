#include "copyright.h"
/*============================================================================*/
/*! \file output.c
 *  \brief Controls output of data.
 *
 * PURPOSE: Controls output of data.  Output is divided into three types:
 * - 1. dump_*(): ALL variables are written in * format over WHOLE grid
 * - 2. output_*(): ONE variable is written in * format with various options
 * - 3. restarts: special form of a dump, includes extra data
 *   The number and types of outputs are all controlled by <ouputN> blocks in
 *   the input files, parsed by the functions in par.c.  
 *
 * TOTAL NUMBER of outputs is controlled by 'maxout' in <job> block in input
 *   file.  Only the first 'maxout' <outputN> blocks are processed, where
 *   N < maxout.  If N > maxout, that <outputN> block is ignored.
 *
 * OPTIONS available in an <outputN> block are:
 * - out       = cons,prim,d,M1,M2,M3,E,B1c,B2c,B3c,ME,V1,V2,V3,P,S,cs2,G
 * - out_fmt   = bin,hst,tab,rst,vtk,pdf,pgm,ppm
 * - dat_fmt   = format string used to write tabular output (e.g. %12.5e)
 * - dt        = problem time between outputs
 * - time      = time of next output (useful for restarts)
 * - id        = any string
 * - dmin/dmax = max/min applied to all outputs
 * - palette   = rainbow,jh_colors,idl1,idl2,step8,step32,heat
 * - x1,x2,x3  = range over which data is averaged or sliced; see parse_slice()
 * - usr_expr_flag = 1 for user-defined expression (defined in problem.c)
 * - level,domain = integer indices of level and domain to be output with SMR
 *   
 * EXAMPLE of an <outputN> block for a VTK dump:
 * - <output1>
 * - out_fmt = vtk
 * - out_dt  = 0.1
 *
 * EXAMPLE of an <outputN> block for a ppm image of a x1-x2 slice with data
 * averaged over 0.5-10 in x3 in ppm format:
 * - <output5>
 * - out_fmt = ppm
 * - dt      = 100.0
 * - out     = d
 * - id      = d
 * - x3      = 0.5:10.0
 * - dmin    = 0.25
 * - dmax    = 2.9
 * - palette = rainbow
 *
 * EXAMPLE of an <outputN> block for restarts:
 * - <ouput3>
 * - out_fmt = rst
 * - out_dt  = 1.0
 *
 * CONTROL of output proceeds as follows:
 *  -init_output(): called by main(), parses the first maxout output blocks.
 *     The info in each block is stored in an element of a global array of
 *     "Output_s" structures, including a pointer to the appropriate output
 *     function.
 *  -data_output(): called in main loop, compares integration time with time 
 *     for output for each element in Output array, and calls output functions.
 *    
 *   To add permanently a new type of output X, write a new function output_X, 
 *   modify init_output() to set the output function pointer when out_fmt=X
 *   in the input file (see below for examples of pgm, ppm, etc.)
 *
 *   See Users Manual to add a problem-specific user-defined output function in
 *   the problem definition file.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - init_output() -
 * - data_output() -
 * - data_output_destruct()
 * - OutData1,2,3()   -
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - expr_*()
 * - get_expr()
 * - free_output()
 * - parse_slice()
 * - getRGB()
 *
 * VARIABLE TYPE AND STRUCTURE DEFINITIONS: none
 *============================================================================*/

#ifndef __POWERPC__
#define HAVE_DLFCN
#endif

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_DLFCN
#include <dlfcn.h>
#endif

#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "palette.h"
#include "prototypes.h"
#include "particles/prototypes.h"

#define MAXOUT_DEFAULT     10

static int out_count = 0;           /* Number of elements in the OutArray */
static OutputS *OutArray = NULL;    /* Array of Output modes */
static OutputS rst_out;             /* Restart Output */
static int rst_flag = 0;            /* (0,1) -> Restart Outputs are (off,on) */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   expr_*
 *   get_expr
 *   free_output
 *   parse_slice
 *   getRGB
 *============================================================================*/

Real expr_d  (const GridS *pG, const int i, const int j, const int k);
Real expr_M1 (const GridS *pG, const int i, const int j, const int k);
Real expr_M2 (const GridS *pG, const int i, const int j, const int k);
Real expr_M3 (const GridS *pG, const int i, const int j, const int k);
Real expr_E  (const GridS *pG, const int i, const int j, const int k);
Real expr_B1c(const GridS *pG, const int i, const int j, const int k);
Real expr_B2c(const GridS *pG, const int i, const int j, const int k);
Real expr_B3c(const GridS *pG, const int i, const int j, const int k);
Real expr_ME (const GridS *pG, const int i, const int j, const int k);
Real expr_V1 (const GridS *pG, const int i, const int j, const int k);
Real expr_V2 (const GridS *pG, const int i, const int j, const int k);
Real expr_V3 (const GridS *pG, const int i, const int j, const int k);
Real expr_P  (const GridS *pG, const int i, const int j, const int k);
Real expr_cs2(const GridS *pG, const int i, const int j, const int k);
Real expr_S  (const GridS *pG, const int i, const int j, const int k);
#ifdef SPECIAL_RELATIVITY
Real expr_G  (const GridS *pG, const int i, const int j, const int k);
#endif
#ifdef PARTICLES
/*! \fn Real expr_dpar (const GridS *pG, const int i, const int j, const int k)
 *  \brief Particle density. */
extern Real expr_dpar (const GridS *pG, const int i, const int j, const int k);
/*! \fn Real expr_M1par(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Particle 1-momentum */
extern Real expr_M1par(const GridS *pG, const int i, const int j, const int k);
/*! \fn Real expr_M2par(const GridS *pG, const int i, const int j, const int k)
 *  \brief Particle 2-momentum */
extern Real expr_M2par(const GridS *pG, const int i, const int j, const int k);
/*! \fn Real expr_M3par(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Particle 3-momentum */
extern Real expr_M3par(const GridS *pG, const int i, const int j, const int k);
/*! \fn Real expr_V1par(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Particle 1-velocity */
extern Real expr_V1par(const GridS *pG, const int i, const int j, const int k);
/*! \fn Real expr_V2par(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Particle 2-velocity */
extern Real expr_V2par(const GridS *pG, const int i, const int j, const int k);
/*! \fn Real expr_V3par(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Particle 3-velocity */
extern Real expr_V3par(const GridS *pG, const int i, const int j, const int k);
int check_particle_binning(char *out);
#endif
static ConsFun_t getexpr(const int n, const char *expr);
static void free_output(OutputS *pout);
static void parse_slice(char *block, char *axname, Real *l, Real *u, int *flag);
float *getRGB(char *name);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void init_output(MeshS *pM)
 *  \brief Initializes data output. */

void init_output(MeshS *pM)
{
  int i,j,outn,maxout,nl,nd;
  char block[80], *fmt, defid[10];
  OutputS new_out;
  int usr_expr_flag;

  maxout = par_geti_def("job","maxout",MAXOUT_DEFAULT);

/* allocate output array */

  if((OutArray = malloc(maxout*sizeof(OutputS))) == NULL){
    ath_error("[init_output]: Error allocating output array\n");
  }

/*--- loop over maxout output blocks, reading parameters into a temporary -----*
 *--- OutputS called new_out --------------------------------------------------*/

  for (outn=1; outn<=maxout; outn++) {

    sprintf(block,"output%d",outn);

/* An output format or output name is required.
 * If neither is present we write an error message and move on. */
    if((par_exist(block,"out_fmt") == 0) && (par_exist(block,"name") == 0)){
      ath_perr(-1,"[init_output]: neither %s/out_fmt, nor %s/name exist\n",
	       block, block);
      continue;
    }

/* Zero (NULL) all members of the temporary OutputS structure "new_out" */
    memset(&new_out,0,sizeof(OutputS));

/* The next output time and number */
    new_out.t   = par_getd_def(block,"time",pM->time);
    new_out.num = par_geti_def(block,"num",0);

    new_out.dt  = par_getd(block,"dt");
    new_out.n   = outn;

/* level and domain number can be specified with SMR  */
    nl = new_out.nlevel = par_geti_def(block,"level",-1);
    nd = new_out.ndomain = par_geti_def(block,"domain",-1);

    if (par_exist(block,"dat_fmt")) new_out.dat_fmt = par_gets(block,"dat_fmt");

/* set id in output filename to input string if present, otherwise use "outN"
 * as default, where N is output number */
    sprintf(defid,"out%d",outn);
    new_out.id = par_gets_def(block,"id",defid);

    if(par_exist(block,"out_fmt")) 
      fmt = new_out.out_fmt = par_gets(block,"out_fmt");

/* out:     controls what variable can be output (all, prim, or any of expr_*)
 * out_fmt: controls format of output (single variable) or dump (all cons/prim)
 * if "out" doesn't exist, we assume 'cons' variables are meant to be dumped */

    new_out.out = par_gets_def(block,"out","cons");

#ifdef PARTICLES
    /* check input for particle binning (=1, default) or not (=0) */
    new_out.out_pargrid = par_geti_def(block,"pargrid",
                                       check_particle_binning(new_out.out));
    if ((new_out.out_pargrid < 0) || (new_out.out_pargrid >1)) {
      ath_perr(-1,"[init_output]: %s/pargrid must be 0 or 1\n",
	       block, block);
      continue;
    }

/* set particle property selection function. By default, will select all the
 * particles. Used only when particle output is called, otherwise useless. */
    if(par_exist(block,"par_prop")) {
      new_out.par_prop = get_usr_par_prop(par_gets(block,"par_prop"));
      if (new_out.par_prop == NULL) {
        ath_pout(0,"[init_output]: Particle selection function not found! \
Now use the default one.\n");
        new_out.par_prop = property_all;
      }
    }
    else
      new_out.par_prop = property_all;
#endif

/* First handle data dumps of all CONSERVED variables (out=cons) */

    if(strcmp(new_out.out,"cons") == 0){
/* check for valid data dump: dump format = {bin, hst, tab, rst, vtk} */
      if(par_exist(block,"name")){
	/* The output function is user defined - get its name */
	char *name = par_gets(block,"name");
	/* Get a pointer to the output function via its name */
	new_out.out_fun = get_usr_out_fun(name);
	if(new_out.out_fun == NULL){
	  free_output(&new_out);
	  ath_error("Unsupported output named %s in %s/out_fmt=%s\n",
		    name,block,fmt);
	}
	free(name);  name = NULL;
	goto add_it;
      }
      else if (strcmp(fmt,"bin")==0){
	new_out.out_fun = dump_binary;
#ifdef PARTICLES
        new_out.out_pargrid = 1; /* bin particles */
#endif
#ifdef PARTICLES
        new_out.out_pargrid = 1; /* bin particles */
#endif
	goto add_it;
      }
      else if (strcmp(fmt,"hst")==0){
	new_out.out_fun = dump_history;
	goto add_it;
      }
#ifdef PARTICLES
      else if (strcmp(fmt,"phst")==0){
        new_out.fun = dump_particle_history;/* do not bin particles (default) */
        goto add_it;
      }
#endif
      else if (strcmp(fmt,"tab")==0){
	new_out.out_fun = dump_tab_cons;
#ifdef PARTICLES
        new_out.out_pargrid = 1; /* bin particles */
#endif
	goto add_it;
      }
      else if (strcmp(fmt,"rst")==0){
	new_out.res_fun = dump_restart;
        rst_flag = 1;
        rst_out = new_out;
	ath_pout(0,"Added out%d\n",outn);
	continue;
      }
      else if (strcmp(fmt,"vtk")==0){
	new_out.out_fun = dump_vtk;
#ifdef PARTICLES
        new_out.out_pargrid = 1; /* bin particles */
#endif
	goto add_it;
      }
#ifdef PARTICLES
      else if (strcmp(fmt,"lis")==0){ /* dump particle list */
	new_out.fun = dump_particle_binary; /* do not bin particles (default) */
	goto add_it;
      }
#endif
      else{    /* Unknown data dump (fatal error) */
	ath_error("Unsupported dump mode for %s/out_fmt=%s for out=cons\n",
          block,fmt);
      }
    }

/* Next handle data dumps of all PRIMITIVE variables (out=prim) */

    if(strcmp(new_out.out,"prim") == 0){
/* check for valid data dump: dump format = {bin, tab, vtk} */
      if(par_exist(block,"name")){
        /* The output function is user defined - get its name */
        char *name = par_gets(block,"name");
        /* Get a pointer to the output function via its name */
        new_out.out_fun = get_usr_out_fun(name);
        if(new_out.out_fun == NULL){
          free_output(&new_out);
          ath_error("Unsupported output named %s in %s/out_fmt=%s\n",
                    name,block,fmt);
        }
        free(name);  name = NULL;
        goto add_it;
      }
      else if (strcmp(fmt,"bin")==0){
        new_out.out_fun = dump_binary;
        goto add_it;
      }
      else if (strcmp(fmt,"tab")==0){
        new_out.out_fun = dump_tab_prim;
#ifdef PARTICLES
        new_out.out_pargrid = 1; /* bin particles */
#endif
        goto add_it;
      }
      else if (strcmp(fmt,"vtk")==0){
        new_out.out_fun = dump_vtk;
        goto add_it;
      }
      else{    /* Unknown data dump (fatal error) */
        ath_error("Unsupported dump mode for %s/out_fmt=%s for out=prim\n",
          block,fmt);
      }
    }

/* Now handle data outputs (ouput of SINGLE variable).  There are lots more
 * options for outputs than dumps.  Need to choose variable, format, size
 * of domain to be output, scaling to min/max (if necessary),...    */

/* Is this a user defined expression? This allows the user to output any
 * problem-specific quantity using the formats and options supported here.
 * new_out.out must point to an expression defined in the user's problem.c */

    if(par_exist(block,"usr_expr_flag"))
      usr_expr_flag = par_geti(block,"usr_expr_flag");
    else
      usr_expr_flag = 0;

/* Get the expression function pointer */
    if(usr_expr_flag)
      new_out.expr = get_usr_expr(new_out.out);
    else
      new_out.expr = getexpr(outn, new_out.out);

    if (new_out.expr == NULL) {
      ath_perr(-1,"Could not parse expression %s, skipping it\n",
	      new_out.out);
      free_output(&new_out);
      continue;
    }

/* x1, x2, x3:  parse coordinate range for slicing and averaging */
    new_out.ndim = 1;
    for (i=1; i<3; i++) if (pM->Nx[i]>1) new_out.ndim++;

    new_out.x1l = pM->RootMinX[0];
    new_out.x1u = pM->RootMaxX[0];
    new_out.reduce_x1 = 0;
    parse_slice(block,"x1",&new_out.x1l,&new_out.x1u,&new_out.reduce_x1);
    if (new_out.reduce_x1 != 0) new_out.ndim--;

    new_out.x2l = pM->RootMinX[1];
    new_out.x2u = pM->RootMaxX[1];
    new_out.reduce_x2 = 0;
    parse_slice(block,"x2",&new_out.x2l,&new_out.x2u,&new_out.reduce_x2);
    if (pM->Nx[1] > 1 && new_out.reduce_x2 != 0) new_out.ndim--;

    new_out.x3l = pM->RootMinX[2];
    new_out.x3u = pM->RootMaxX[2];
    new_out.reduce_x3 = 0;
    parse_slice(block,"x3",&new_out.x3l,&new_out.x3u,&new_out.reduce_x3);
    if (pM->Nx[2] > 1 && new_out.reduce_x3 != 0) new_out.ndim--;

    if (new_out.ndim <= 0) ath_error("Too many slices specified in %s\n",block);

/* dmin/dmax & sdmin/sdmax */
    if(par_exist(block,"dmin") != 0){ /* Use a fixed minimum scale? */
      new_out.sdmin = 1;
      new_out.dmin = par_getd(block,"dmin");
    }
    new_out.gmin = (HUGE_NUMBER);

    if(par_exist(block,"dmax") != 0){ /* Use a fixed maximum scale? */
      new_out.sdmax = 1;
      new_out.dmax = par_getd(block,"dmax");
    }
    new_out.gmax = -1.0*(HUGE_NUMBER);

/* palette: default is rainbow */
    if (strcmp(fmt,"ppm") == 0) {
      new_out.palette = par_gets_def(block,"palette","rainbow");

      new_out.rgb = getRGB(new_out.palette);
      if ( (new_out.der = (float *) malloc(3*256*sizeof(float))) == NULL) {
	free_output(&new_out);
	ath_error("[init_output]: malloc returned a NULL pointer\n");
      }
      for(j=0; j<3; j++)    /* compute derivates to speed up interpolations */
	for (i=0; i<255; i++)
	  new_out.der[3*i+j] = new_out.rgb[3*(i+1)+j] - new_out.rgb[3*i+j];
    }

/* check for valid data output option (output of single variables)
 *  output format = {pdf, pgm, ppm, tab, vtk}.  Note for pdf and tab
 *  outputs we also get the format for the print statements.
 */

    if(par_exist(block,"name")){
      /* The output function is user defined - get its name */
      char *name = par_gets(block,"name");
      /* Get a pointer to the output function via its name */
      new_out.out_fun = get_usr_out_fun(name);
      if(new_out.out_fun == NULL){
	free_output(&new_out);
	ath_error("Unsupported output named %s in %s/out_fmt=%s\n",
		  name,block,fmt);
      }
      free(name);  name = NULL;
    }
    else if (strcmp(fmt,"pdf")==0)
      new_out.out_fun = output_pdf;
    else if (strcmp(fmt,"pgm")==0)
      new_out.out_fun = output_pgm;
    else if (strcmp(fmt,"ppm")==0)
      new_out.out_fun = output_ppm;
    else if (strcmp(fmt,"vtk")==0)
      new_out.out_fun = output_vtk;
    else if (strcmp(fmt,"tab")==0)
      new_out.out_fun = output_tab;
    else {
/* unknown output format is fatal */
      free_output(&new_out);
      ath_error("Unsupported %s/out_fmt=%s\n",block,fmt);
    }

  add_it:

/* Now copy data in "new_out" into OutArray structure, and increment index. */
    
    ath_pout(1,"OUTPUT: %d %d %s %s [%g : %g]\n",
             new_out.n, new_out.ndim, new_out.out_fmt,
             new_out.out, new_out.dmin, new_out.dmax); /* DEBUG */

    OutArray[out_count] = new_out;
    out_count++;
    ath_pout(0,"Added out%d\n",outn);

  } /*---------------------- end loop over output blocks ----------------------*/

}

/*----------------------------------------------------------------------------*/
/*! \fn void data_output(MeshS *pM, const int flag)
 *  \brief Called by main(), tests whether time for output, and calls
 *   appropriate output functions.  
 *
 *   Setting the input argument flag=1 forces a
 *   write of all output's.  If the input argument flag=0, then only those
 *   output's whose next output time has passed will be written.        */

void data_output(MeshS *pM, const int flag)
{
  int n;
  int dump_flag[MAXOUT_DEFAULT+1];
  char block[80];

/* Loop over all elements in output array
 * set dump flag to input argument, check whether time for output */

  for (n=0; n<out_count; n++) {
    dump_flag[n] = flag;
    if (pM->time >= OutArray[n].t) {
      OutArray[n].t += OutArray[n].dt;
      dump_flag[n] = 1;
    }
  }

/* Now check for restart dump, and make restart if dump_flag != 0 */

  if(rst_flag){
    dump_flag[out_count] = flag;
    if(pM->time >= rst_out.t){
      rst_out.t += rst_out.dt;
      dump_flag[out_count] = 1;
    }

    if(dump_flag[out_count] != 0){
/* Update the output numbers and times in the output blocks */
      for(n=0; n<out_count; n++){
/* User enrolled outputs have outn < 0 */
	if(OutArray[n].n > 0){
	  sprintf(block,"output%d",OutArray[n].n);
          if (dump_flag[n] != 0) {
/* About to write this output, so increase the output
 * number given in the restart file */
	    par_seti(block,"num","%d",OutArray[n].num+1,"Next Output Number");
          } else {
	    par_seti(block,"num","%d",OutArray[n].num,"Next Output Number");
          }
	  par_setd(block,"time","%.15e",OutArray[n].t,"Next Output Time");
	}
      }
/* Now do the same for the restart output block */
      sprintf(block,"output%d",rst_out.n);
      par_seti(block,"num","%d",rst_out.num+1,"Next Output Number");
      par_setd(block,"time","%.15e",rst_out.t,"Next Output Time");

/* Write the restart file */
      (*(rst_out.res_fun))(pM,&(rst_out));

      rst_out.num++;
    }
  }

/* Loop over all elements in output array, if dump_flag != 0, make output */

  for (n=0; n<out_count; n++) {
    if(dump_flag[n] != 0) {

#ifdef PARTICLES
      if (OutArray[n].out_pargrid == 1)      /* binned particles are output */
        particle_to_grid(pM, OutArray[n].par_prop);
#endif
      (*OutArray[n].out_fun)(pM,&(OutArray[n]));

      OutArray[n].num++;

    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void data_output_destruct(void) 
 *  \brief Free all memory associated with Output, called by
 *   main() at end of run */

void data_output_destruct(void)
{
  int i;
  double global_min, global_max;
#ifdef MPI_PARALLEL
  int ierr;
#endif

  for (i=0; i<out_count; i++) {

/* print the global min/max computed over the calculation */

    if (OutArray[i].out != NULL){
      if((strcmp(OutArray[i].out,"cons") != 0) &&
         (strcmp(OutArray[i].out,"prim") != 0)){
/* get global min/max with MPI calculation */
#ifdef MPI_PARALLEL
        ierr = MPI_Allreduce(&OutArray[i].gmin, &global_min, 1, MPI_DOUBLE,
          MPI_MIN, MPI_COMM_WORLD);
        ierr = MPI_Allreduce(&OutArray[i].gmax, &global_max, 1, MPI_DOUBLE,
          MPI_MAX, MPI_COMM_WORLD);
#else
        global_min = OutArray[i].gmin;
        global_max = OutArray[i].gmax;
#endif
	ath_pout(0,"Global min/max for %s: %g %g\n",OutArray[i].out,
		 global_min, global_max);
      }

      free(OutArray[i].out);
    }
    if (OutArray[i].out_fmt != NULL) free(OutArray[i].out_fmt);
    if (OutArray[i].dat_fmt != NULL) free(OutArray[i].dat_fmt);
    if (OutArray[i].id      != NULL) free(OutArray[i].id);
  }

  if(rst_flag){
    if (rst_out.out     != NULL) free(rst_out.out);
    if (rst_out.out_fmt != NULL) free(rst_out.out_fmt);
    if (rst_out.dat_fmt != NULL) free(rst_out.dat_fmt);
    if (rst_out.id      != NULL) free(rst_out.id);
  }

  if (OutArray != NULL) {
    free(OutArray);
    OutArray = NULL;
    out_count = 0;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn Real ***OutData3(GridS *pgrid, OutputS *pout, int *Nx1, int *Nx2, 
 *                       int *Nx3)
 *  \brief Creates 3D array of output data with dimensions equal to Grid
 * using output expression (function pointer) stored in Output structure.
 *
 * Dimensions of array created also returned in arguments. */

Real ***OutData3(GridS *pgrid, OutputS *pout, int *Nx1, int *Nx2, int *Nx3)
{
  Real ***data;
  int i,j,k,il,jl,kl,iu,ju,ku;

  if (pout->ndim != 3) ath_error("[OutData3] <output%d> %s is %d-D, not 3-D\n",
    pout->n,pout->out, pout->ndim);

#ifdef WRITE_GHOST_CELLS
  if(pgrid->Nx[0] > 1){
    il = pgrid->is - nghost;
    iu = pgrid->ie + nghost;
  } else{
    il = pgrid->is;
    iu = pgrid->ie;
  }

  if(pgrid->Nx[1] > 1){
    jl = pgrid->js - nghost;
    ju = pgrid->je + nghost;
  } else{
    jl = pgrid->js;
    ju = pgrid->je;
  }

  if(pgrid->Nx[2] > 1){
    kl = pgrid->ks - nghost;
    ku = pgrid->ke + nghost;
  } else{
    kl = pgrid->ks;
    ku = pgrid->ke;
  }
#else
  il = pgrid->is;
  iu = pgrid->ie;
  jl = pgrid->js;
  ju = pgrid->je;
  kl = pgrid->ks;
  ku = pgrid->ke;
#endif
  *Nx1 = iu-il+1;
  *Nx2 = ju-jl+1;
  *Nx3 = ku-kl+1;

  data = (Real***) calloc_3d_array(*Nx3,*Nx2,*Nx1,sizeof(Real));
  if (data == NULL) ath_error("[OutData3] Error creating 3D data array\n");
  for (k=0; k<*Nx3; k++)
    for (j=0; j<*Nx2; j++)
      for (i=0; i<*Nx1; i++)
        data[k][j][i] = (*pout->expr)(pgrid,i+il,j+jl,k+kl);
  return data;
}

/*----------------------------------------------------------------------------*/
/*! \fn Real **OutData2(GridS *pgrid, OutputS *pout, int *Nx1, int *Nx2)
 *  \brief Creates 2D array of output data with two dimensions equal to Grid
 * and one dimension reduced according to range stored in x1l/x1u, etc.  
 *
 * Data is computed using output expression (function pointer) stored in Output
 * structure.  If slice range lies outside of coordinate range in Grid, the
 * NULL pointer is returned.  Dimensions of array created are also returned in
 * arguments */

Real **OutData2(GridS *pgrid, OutputS *pout, int *Nx1, int *Nx2)
{
  Real **data;
  Real factor,x1fc,x2fc,x3fc;
  int Nx3;
  int i,j,k,il,jl,kl,iu,ju,ku;
  int istart,iend,jstart,jend,kstart,kend;

  if (pout->ndim != 2) ath_error("[OutData2] <output%d> %s is %d-D, not 2-D\n",
    pout->n,pout->out, pout->ndim);

#ifdef WRITE_GHOST_CELLS
  if(pgrid->Nx[0] > 1){
    il = pgrid->is - nghost;
    iu = pgrid->ie + nghost;
  } else{
    il = pgrid->is;
    iu = pgrid->ie;
  }

  if(pgrid->Nx[1] > 1){
    jl = pgrid->js - nghost;
    ju = pgrid->je + nghost;
  } else{
    jl = pgrid->js;
    ju = pgrid->je;
  }

  if(pgrid->Nx[2] > 1){
    kl = pgrid->ks - nghost;
    ku = pgrid->ke + nghost;
  } else{
    kl = pgrid->ks;
    ku = pgrid->ke;
  }
#else
  il = pgrid->is;
  iu = pgrid->ie;
  jl = pgrid->js;
  ju = pgrid->je;
  kl = pgrid->ks;
  ku = pgrid->ke;
#endif
  *Nx1 = iu-il+1;
  *Nx2 = ju-jl+1;
  Nx3 = ku-kl+1;

/* data is already 2D in 2D simulations */

  if (pgrid->Nx[2] == 1) {
    data = (Real**) calloc_2d_array(*Nx2,*Nx1,sizeof(Real));
    if (data == NULL) ath_error("[OutData2] Error creating 2D data array\n");
    for (j=0; j<*Nx2; j++) {
      for (i=0; i<*Nx1; i++) {
	data[j][i] = (*pout->expr)(pgrid,i+il,j+jl,kl);
      }
    }
    return data;
  }

/* Slice 3D data into 2D arrays according to reduce_x* flags */
	  
/* Nx3,Nx2,Nx1 -> Nx2,Nx1 */
  if (pout->reduce_x3 != 0) {
    if (pout->x3u < pgrid->MinX[2] || pout->x3l >= pgrid->MaxX[2]) return NULL;

    /* find k indices of slice range */
    k=kl+1;
    fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
    while (pout->x3l >= x3fc) {
      k++;
      fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
    }
    kstart = k-1;

    k=ku;
    fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
    while (pout->x3u < x3fc) {
      k--;
      fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
    }
    kend = k;

    /* allocate array and compute data */
    data = (Real**) calloc_2d_array(*Nx2,*Nx1,sizeof(Real));
    if (data == NULL) ath_error("[OutData2] Error creating 2D data array\n");
    factor = 1.0/(kend - kstart + 1);
    for (j=0; j<*Nx2; j++) {
      for (i=0; i<*Nx1; i++) {
	data[j][i] = 0.0;
	for (k=kstart; k<=kend; k++)
	  data[j][i] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
	data[j][i] *= factor;
      }
    }

/* Nx3,Nx2,Nx1 -> Nx3,Nx1 */
  } else if (pout->reduce_x2 != 0) {
    if (pout->x2u < pgrid->MinX[1] || pout->x2l >= pgrid->MaxX[1]) return NULL;

    /* find j indices of slice range */
    j=jl+1;
    fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x2l >= x2fc) {
      j++;
      fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    }
    jstart = j-1;

    j=ju;
    fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x2u < x2fc) {
      j--;
      fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    }
    jend = j;

    /* allocate array and compute data */
    data = (Real**) calloc_2d_array(Nx3,*Nx1,sizeof(Real));
    if (data == NULL) ath_error("[OutData2] Error creating 2D data array\n");
    factor = 1.0/(jend - jstart + 1);
    for (k=0; k<Nx3; k++) {
      for (i=0; i<*Nx1; i++) {
	data[k][i] = 0.0;
	for (j=jstart; j<=jend; j++)
	  data[k][i] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
	data[k][i] *= factor;
      }
    }
    *Nx2 = Nx3; /* return second dimension of array created */

/* Nx3,Nx2,Nx1 -> Nx3,Nx2 */
  } else if (pout->reduce_x1 != 0) {
    if (pout->x1u < pgrid->MinX[0] || pout->x1l >= pgrid->MaxX[0]) return NULL;

    /* find i indices of slice range */
    i=il+1;
    fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x1l >= x1fc) {
      i++;
      fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    }
    istart = i-1;

    i=iu;
    fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x1u < x1fc) {
      i--;
      fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    }
    iend = i;

    /* allocate array and compute data */

    data = (Real**) calloc_2d_array(Nx3,*Nx2,sizeof(Real));
    if (data == NULL) ath_error("[OutData2] Error creating 2D data array\n");
    factor = 1.0/(iend - istart + 1);
    for (k=0; k<Nx3; k++) {
      for (j=0; j<*Nx2; j++) {
	data[k][j] = 0.0;
	for (i=istart; i<=iend; i++)
	  data[k][j] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
	data[k][j] *= factor;
      }
    }
    *Nx1 = *Nx2;
    *Nx2 = Nx3; /* return dimensions of array created */
  } else
    ath_perr(-1,"[OutData2]: Should not reach here\n");
  return data;
}

/*----------------------------------------------------------------------------*/
/*! \fn Real *OutData1(GridS *pgrid, OutputS *pout, int *Nx1)
 *  \brief Creates 1D array of output data with one dimensions equal to Grid
 * and two dimensions reduced according to range stored in x1l/x1u, etc.  
 *
 * Data is computed using output expression (function pointer) stored in Output 
 * structure.  If slice range lies outside of coordinate range in Grid, the
 * NULL pointer is returned.  Dimension of array created is also returned in
 * arguments. */

Real *OutData1(GridS *pgrid, OutputS *pout, int *Nx1)
{
  Real *data;
  Real factor,x1fc,x2fc,x3fc;
  int Nx2, Nx3;
  int i,j,k,il,jl,kl,iu,ju,ku;
  int istart,iend,jstart,jend,kstart,kend;

  if (pout->ndim != 1) ath_error("[OutData1] <output%d> %s is %d-D, not 1-D\n",
    pout->n,pout->out, pout->ndim);

#ifdef WRITE_GHOST_CELLS
  if(pgrid->Nx[0] > 1){
    il = pgrid->is - nghost;
    iu = pgrid->ie + nghost;
  } else{
    il = pgrid->is;
    iu = pgrid->ie;
  }

  if(pgrid->Nx[1] > 1){
    jl = pgrid->js - nghost;
    ju = pgrid->je + nghost;
  } else{
    jl = pgrid->js;
    ju = pgrid->je;
  }

  if(pgrid->Nx[2] > 1){
    kl = pgrid->ks - nghost;
    ku = pgrid->ke + nghost;
  } else{
    kl = pgrid->ks;
    ku = pgrid->ke;
  }
#else
  il = pgrid->is;
  iu = pgrid->ie;
  jl = pgrid->js;
  ju = pgrid->je;
  kl = pgrid->ks;
  ku = pgrid->ke;
#endif
  *Nx1 = iu-il+1;
  Nx2 = ju-jl+1;
  Nx3 = ku-kl+1;

/* data is already 1D in 1D simulations */

  if (pgrid->Nx[1] == 1) {
    data = (Real*) calloc_1d_array(*Nx1,sizeof(Real));
    if (data == NULL) ath_error("[OutData1] Error creating 1D data array\n");
    for (i=0; i<*Nx1; i++) {
      data[i] = (*pout->expr)(pgrid,i+il,jl,kl);
    }
    return data;
  }

/* Slice 2D and 3D data into 2D arrays according to reduce_x* flags */

/* Nx3,Nx2,Nx1 -> Nx1 */
  if (pout->reduce_x1 == 0) {
    if (pout->x3u < pgrid->MinX[2] || pout->x3l >= pgrid->MaxX[2] ||
        pout->x2u < pgrid->MinX[1] || pout->x2l >= pgrid->MaxX[1]) return NULL;

    /* find k indices of slice range */
    if (pgrid->Nx[2] == 1){
      kstart=kl;
      kend=ku;
    } else {
      k=kl+1;
      fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      while (pout->x3l >= x3fc) {
        k++;
        fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      }
      kstart = k-1;

      k=ku;
      fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      while (pout->x3u < x3fc) {
        k--;
        fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      }
      kend = k;
    }

    /* find j indices of slice range */
    j=jl+1;
    fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x2l >= x2fc) {
      j++;
      fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    }
    jstart = j-1;

    j=ju;
    fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x2u < x2fc) {
      j--;
      fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    }
    jend = j;

    /* allocate array and compute data */
    data = (Real*) calloc_1d_array(*Nx1,sizeof(Real));
    factor = 1.0/(kend - kstart + 1)/(jend - jstart + 1);
    for (i=0; i<*Nx1; i++) {
      data[i] = 0.0;
      for (k=kstart; k<=kend; k++)
	for (j=jstart; j<=jend; j++)
	  data[i] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
      data[i] *= factor;
    }

/* Nx3,Nx2,Nx1 -> Nx2 */
  } else if (pout->reduce_x2 == 0) {
    if (pout->x3u < pgrid->MinX[2] || pout->x3l >= pgrid->MaxX[2] ||
        pout->x1u < pgrid->MinX[0] || pout->x1l >= pgrid->MaxX[0]) return NULL;

    /* find k indices of slice range */
    if (pgrid->Nx[2] == 1){
      kstart=kl;
      kend=ku;
    } else {
      k=kl+1;
      fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      while (pout->x3l >= x3fc) {
        k++;
        fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      }
      kstart = k-1;

      k=ku;
      fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      while (pout->x3u < x3fc) {
        k--;
        fc_pos(pgrid,il,jl,k,&x1fc,&x2fc,&x3fc);
      }
      kend = k;
    }

    /* find i indices of slice range */
    i=il+1;
    fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x1l >= x1fc) {
      i++;
      fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    }
    istart = i-1;

    i=iu;
    fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x1u < x1fc) {
      i--;
      fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    }
    iend = i;

    /* allocate array and compute data */
    data = (Real*) calloc_1d_array(Nx2,sizeof(Real));
    factor = 1.0/(kend - kstart + 1)/(iend - istart + 1);
    for (j=0; j<Nx2; j++) {
      data[j] = 0.0;
      for (k=kstart; k<=kend; k++)
	for (i=istart; i<=iend; i++)
	  data[j] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
      data[j] *= factor;
    }
    *Nx1 = Nx2; /* return dimensions of array created */

/* Nx3,Nx2,Nx1 -> Nx3. Data must be 3D in this case. */
  } else if (pout->reduce_x3 == 0) {
    if (pout->x2u < pgrid->MinX[1] || pout->x2l >= pgrid->MaxX[1] ||
        pout->x1u < pgrid->MinX[0] || pout->x1l >= pgrid->MaxX[0]) return NULL;

    /* find j indices of slice range */
    j=jl+1;
    fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x2l >= x2fc) {
      j++;
      fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    }
    jstart = j-1;

    j=ju;
    fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x2u < x2fc) {
      j--;
      fc_pos(pgrid,il,j,kl,&x1fc,&x2fc,&x3fc);
    }
    jend = j;

    /* find i indices of slice range */
    i=il+1;
    fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x1l >= x1fc) {
      i++;
      fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    }
    istart = i-1;

    i=iu;
    fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    while (pout->x1u < x1fc) {
      i--;
      fc_pos(pgrid,i,jl,kl,&x1fc,&x2fc,&x3fc);
    }
    iend = i;

    /* allocate array and compute data */
    data = (Real*) calloc_1d_array(Nx3,sizeof(Real));
    factor = 1.0/(jend - jstart + 1)/(iend - istart + 1);
    for (k=0; k<Nx3; k++) {
      data[k] = 0.0;
      for (j=jstart; j<=jend; j++)
	for (i=istart; i<=iend; i++)
	  data[k] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
      data[k] *= factor;
    }
    *Nx1 = Nx3; /* return dimensions of array created */
  } else {
    ath_perr(-1,"[OutData1]: Should not reach here\n");
  }

  return data;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/* expr_*: where * are the conserved variables d,M1,M2,M3,E */
/*! \fn Real expr_d(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Density */
Real expr_d(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].d;
}
/*! \fn Real expr_M1(const GridS *pG, const int i, const int j, const int k) 
 *  \brief 1-component of momentum */ 
Real expr_M1(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M1;
}
/*! \fn Real expr_M2(const GridS *pG, const int i, const int j, const int k) 
 *  \brief 2-component of momentum */
Real expr_M2(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M2;
}
/*! \fn Real expr_M3(const GridS *pG, const int i, const int j, const int k) 
 *  \brief 3-component of momentum */
Real expr_M3(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M3;
}
#ifndef BAROTROPIC
/*! \fn Real expr_E(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Total energy */
Real expr_E(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].E;
}
#endif

/*--------------------------------------------------------------------------- */
/* expr_*: where * are magnetic field variables: B1c, B2c, B3c, B^2 */

#ifdef MHD
/*! \fn Real expr_B1c(const GridS *pG, const int i, const int j, const int k) 
 *  \brief 1-component of cell-centered B-field */
Real expr_B1c(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].B1c;
}
/*! \fn Real expr_B2c(const GridS *pG, const int i, const int j, const int k) 
 *  \brief 2-component of cell-centered B-field */
Real expr_B2c(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].B2c;
}
/*! \fn Real expr_B3c(const GridS *pG, const int i, const int j, const int k)  
 *  \brief 3-component of cell-centered B-field */
Real expr_B3c(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].B3c;
}
/*! \fn Real expr_ME(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Magnetic field energy */
Real expr_ME(const GridS *pG, const int i, const int j, const int k) {
  return 0.5*(pG->U[k][j][i].B1c*pG->U[k][j][i].B1c + 
	      pG->U[k][j][i].B2c*pG->U[k][j][i].B2c + 
	      pG->U[k][j][i].B3c*pG->U[k][j][i].B3c);
}
#endif

/*--------------------------------------------------------------------------- */
/* expr_*: where * are the primitive variables */

/*! \fn Real expr_V1(const GridS *pG, const int i, const int j, const int k)  
 *  \brief 1-velocity */
Real expr_V1(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M1/pG->U[k][j][i].d;
}
/*! \fn Real expr_V2(const GridS *pG, const int i, const int j, const int k)  
 *  \brief 2-velocity */
Real expr_V2(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M2/pG->U[k][j][i].d;
}
/*! \fn Real expr_V3(const GridS *pG, const int i, const int j, const int k)  
 *  \brief 3-velocity */
Real expr_V3(const GridS *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M3/pG->U[k][j][i].d;
}

/*! \fn Real expr_P(const GridS *pG, const int i, const int j, const int k) 
 *  \brief Pressure */
Real expr_P(const GridS *pG, const int i, const int j, const int k) {
#ifdef ISOTHERMAL
  return  pG->U[k][j][i].d*Iso_csound2;
#else
  ConsS *gp = &(pG->U[k][j][i]);
  return Gamma_1*(gp->E 
#ifdef MHD
		  - 0.5*(gp->B1c*gp->B1c + gp->B2c*gp->B2c + gp->B3c*gp->B3c)
#endif /* MHD */
		  - 0.5*(gp->M1*gp->M1 + gp->M2*gp->M2 + gp->M3*gp->M3)/gp->d);
#endif /* ISOTHERMAL */
}

/*--------------------------------------------------------------------------- */
/*! \fn Real expr_cs2(const GridS *pG, const int i, const int j, const int k)
 *  \brief Sound speed squared  */

#ifdef ADIABATIC
Real expr_cs2(const GridS *pG, const int i, const int j, const int k)
{
  ConsS *gp = &(pG->U[k][j][i]);
  return (Gamma*Gamma_1*(gp->E 
#ifdef MHD
	  - 0.5*(gp->B1c*gp->B1c + gp->B2c*gp->B2c + gp->B3c*gp->B3c)
#endif /* MHD */
	  - 0.5*(gp->M1*gp->M1 + gp->M2*gp->M2 + gp->M3*gp->M3)/gp->d)/gp->d);
}
#endif /* ADIABATIC */

/*--------------------------------------------------------------------------- */
/*! \fn Real expr_S(const GridS *pG, const int i, const int j, const int k)
 *  \brief entropy = P/d^{Gamma}  */

#ifdef ADIABATIC
Real expr_S(const GridS *pG, const int i, const int j, const int k)
{
  ConsS *gp = &(pG->U[k][j][i]);
  Real P = Gamma_1*(gp->E 
#ifdef MHD
		   - 0.5*(gp->B1c*gp->B1c + gp->B2c*gp->B2c + gp->B3c*gp->B3c)
#endif /* MHD */
		   - 0.5*(gp->M1*gp->M1 + gp->M2*gp->M2 + gp->M3*gp->M3)/gp->d);
  return P/pow((double)gp->d, (double)Gamma);
}
#endif /* ADIABATIC */

/*--------------------------------------------------------------------------- */
/*! \fn Real expr_G(const GridS *pG, const int i, const int j, const int k)
 *  \brief gamma = 1/sqrt(1-v^2)  */

#ifdef SPECIAL_RELATIVITY
Real expr_G(const GridS *pG, const int i, const int j, const int k)
{
  PrimS W;
  W = Cons_to_Prim(&(pG->U[k][j][i]));
  return 1.0/sqrt(1.0 - (SQR(W.V1)+SQR(W.V2)+SQR(W.V3)));
}
#endif /* SPECIAL_RELATIVITY */

/*---------------------------------------------------------------------------_*/
/*! \fn int check_particle_binning(char *out)
 *  \brief Check if particle binning is need */
#ifdef PARTICLES
int check_particle_binning(char *out)
{/* 1: need binning; 0: no */
  if (strcmp(out,"dpar")==0)
    return  1;
  else if (strcmp(out,"M1par")==0)
    return  1;
  else if (strcmp(out,"M2par")==0)
    return  1;
  else if (strcmp(out,"M3par")==0)
    return  1;
  else if (strcmp(out,"V1par")==0)
    return  1;
  else if (strcmp(out,"V2par")==0)
    return  1;
  else if (strcmp(out,"V3par")==0)
    return  1;
  else
    return 0;
}
#endif /* PARTICLES */

/*--------------------------------------------------------------------------- */
/*! \fn static ConsFun_t getexpr(const int n, const char *expr)
 *  \brief Return a function pointer for a simple expression - no parsing.
 *
 *   For a user defined expression, get_usr_expr() in problem.c is used.  */

static ConsFun_t getexpr(const int n, const char *expr)
{
  char ename[32];

  sprintf(ename,"expr_out%d",n);

  if (strcmp(expr,"d")==0)
    return expr_d;
  else if (strcmp(expr,"M1")==0)
    return expr_M1;
  else if (strcmp(expr,"M2")==0)
    return expr_M2;
  else if (strcmp(expr,"M3")==0)
    return expr_M3;
#ifndef BAROTROPIC
  else if (strcmp(expr,"E")==0)
    return expr_E;
#endif /* BAROTROPIC */
#ifdef MHD
  else if (strcmp(expr,"B1c")==0)
    return expr_B1c;
  else if (strcmp(expr,"B2c")==0)
    return expr_B2c;
  else if (strcmp(expr,"B3c")==0)
    return expr_B3c;
  else if (strcmp(expr,"ME")==0)
    return expr_ME;
#endif
  else if (strcmp(expr,"V1")==0)
    return expr_V1;
  else if (strcmp(expr,"V2")==0)
    return expr_V2;
  else if (strcmp(expr,"V3")==0)
    return expr_V3;
  else if (strcmp(expr,"P")==0)
    return expr_P;
#ifdef ADIABATIC
  else if (strcmp(expr,"cs2")==0)
    return  expr_cs2;
#endif /* ADIABATIC */
#ifdef ADIABATIC
  else if (strcmp(expr,"S")==0)
    return  expr_S;
#endif /* ADIABATIC */
#ifdef SPECIAL_RELATIVITY
  else if (strcmp(expr,"G")==0)
    return  expr_G;
#endif /* SPECIAL_RELATIVITY */
#ifdef PARTICLES
  else if (strcmp(expr,"dpar")==0)
    return  expr_dpar;
  else if (strcmp(expr,"M1par")==0)
    return  expr_M1par;
  else if (strcmp(expr,"M2par")==0)
    return  expr_M2par;
  else if (strcmp(expr,"M3par")==0)
    return  expr_M3par;
  else if (strcmp(expr,"V1par")==0)
    return  expr_V1par;
  else if (strcmp(expr,"V2par")==0)
    return  expr_V2par;
  else if (strcmp(expr,"V3par")==0)
    return  expr_V3par;
#endif
  else {
    ath_perr(-1,"Unknown data expression\n");
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn static void free_output(OutputS *pOut)
 *  \brief free memory associated with Output structure.  
 *
 *   Only used when
 *   error occurs in adding a new output; this function frees memory and returns
 *   control to calling function */

static void free_output(OutputS *pOut)
{
  if(pOut->out     != NULL) free(pOut->out);
  if(pOut->out_fmt != NULL) free(pOut->out_fmt);
  if(pOut->dat_fmt != NULL) free(pOut->dat_fmt);
  if(pOut->id      != NULL) free(pOut->id);
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void parse_slice(char *block, char *axname, Real *l, Real *u, 
 *			        int *flag)
 *  \brief Sets the lower and upper bounds of a slice along an axis, 
 *   using values of x1, x2 or x3 in the <output> block.  
 *
 *   These are used to
 *   slice the data for outputs, averaged between l and u.  Valid formats are:
 *   -   x1 = 5e3         both l and u set to 5.0e3
 *   -   x1 = 5.3:10e4    l set to 5.3, u set to 1.0e5
 *   -   x1 = :           l set to RootMinX, u set to RootMaxX
 *   -   x1 = 5:          l set to 5.0, u set to RootMaxX
 *   -   x1 = :10         l set to RootMinX, u set to 10.0
 *   If values for x1,x2,x3 are not set in the <output> block, then l and u
 *   are not changed (default values should be set in calling function).
 *
 *   Note that data is always reduced along the directions specified by x1/2/3.
 *   It is not possible to create a smaller 3D array by specifying ranges for
 *   all three of x1,x2 and x3 at once.  Instead, this would reduce the data to
 *   a single point (not allowed).
 *
 *   This function only parses the input text to extract values of l and u,
 *   the actual slicing and averaging is done by OutData1,2,3().
 */

static void parse_slice(char *block, char *axname, Real *l, Real *u, int *flag)
{
  char *expr, *cp;

  if (par_exist(block,axname)) {
    expr = par_gets(block,axname);
    cp = strchr(expr, ':');
    if (cp) {             /* either ':'  or 'lower:upper'  */
      *cp++ = 0;
      while (*cp && isspace(*cp)) cp++;
      if (*cp)
	*u = atof(cp);
      cp = expr;
      while (*cp && isspace(*cp)) cp++;
      if (*cp)
	*l = atof(cp);
    } else {               /* single slice  */
      *l = *u = atof(expr);
    }
    if (*l > *u) {
      ath_error("[parse_slice]: lower slice limit %d > upper %d in %s\n",
      *l,*u,expr);
    }
    free(expr);
    *flag = 1;
  }

}

/*----------------------------------------------------------------------------*/
/*! \fn float *getRGB(char *name)
 *  \brief function for accessing palettes stored stored in structure RGB.
 *
 *   Compares argument with strings (names) of palettes in RGB, and returns 
 *   pointer to first element of matching palette.  */

float *getRGB(char *name)
{
  int i;

  for (i=0; rgb[i].name && rgb[i].rgb; i++) {
    if (strcmp(name,rgb[i].name) == 0)
      return rgb[i].rgb;
  }

/* failed to find a matching palette: print them all and exit...  */

  ath_perr(-1,"Fatal error: could not find palette=%s, valid names are:\n",
    name);
  for (i=0; rgb[i].name && rgb[i].rgb; i++)
    ath_perr(-1,"%s ",rgb[i].name);
  ath_perr(-1,"\n");
  exit(EXIT_FAILURE);

  return NULL; /* Never reached, avoids compiler warnings */
}
