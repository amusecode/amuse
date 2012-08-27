#include "../copyright.h"
/*============================================================================*/
/*! \file dump_particle_history.c
 *  \brief Functions to write dumps of particle "history" variables in a
 *   formatted table.
 *
 * PURPOSE: Functions to write dumps of particle "history" variables in a
 *   formatted table.  "History" dumps are scalars (usually volume averages)
 *   written periodically, giving the time-evolution of that quantitity.
 *
 *   The particle history dump has two sets. The first set is global, particle
 *   type independent quantities, including the following:
 *   - scal[0] = time
 *   - scal[1] = maximum particle density
 *   - scal[2] = energy dissipation rate from the drag
 *   - scal[3] = maximum stiffness parameter
 *   - scal[4] = particle mass
 *   - scal[5] = particle x1 momentum
 *   - scal[6] = particle x2 momentum
 *   - scal[7] = particle x3 momentum
 *   - scal[8] = particle x1 kinetic energy
 *   - scal[9] = particle x2 kinetic energy
 *   - scal[10] = particle x3 kinetic energy
 *
 *   The second set is particle type dependent quantities, which contains
 *   - array[0] = particle x1 average position
 *   - array[1] = particle x2 average position
 *   - array[2] = particle x3 average position
 *   - array[3] = particle x1 average velocity
 *   - array[4] = particle x2 average velocity
 *   - array[5] = particle x3 average velocity
 *   - array[6] = particle x1 position variation
 *   - array[7] = particle x2 position variation
 *   - array[8] = particle x3 position variation
 *   - array[9] = particle x1 velocity dispersion
 *   - array[10]= particle x2 velocity dispersion
 *   - array[11]= particle x3 velocity dispersion
 *
 * More variables can be hardwired by increasing NSCAL=number of variables, and
 * adding calculation of desired quantities below.
 *
 * Alternatively, up to MAX_USR_H_COUNT new history variables can be added using
 * dump_parhistory_enroll() in the problem generator.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - dump_particle_history()      - Writes variables as formatted table
 * - dump_parhistory_enroll()     - Adds new user-defined history variables   */
/*============================================================================*/

#include <stdio.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"

#ifdef PARTICLES /* endif at the end of the file */

/* Maximum Number of default history dump columns. */
#define NSCAL 11
#define NARAY 12

/* Maximum number of history dump columns that the user routine can add. */
#define MAX_USR_SCAL 10
#define MAX_USR_ARAY 10

/* Array of strings / labels for the user added history column. */
static char *usr_label_scal[MAX_USR_SCAL];
static char *usr_label_aray[MAX_USR_ARAY];

/* Array of history dump function pointers for user added history columns. */
static Parfun_t phst_scalfun[MAX_USR_SCAL];
static Parfun_t phst_arayfun[MAX_USR_ARAY];

static int usr_scal_cnt = 0; /* User History Counter <= MAX_USR_SCAL */
static int usr_aray_cnt = 0; /* User History Counter <= MAX_USR_ARAY */

extern Real expr_dpar(const Grid *pG, const int i, const int j, const int k);

/*============================================================================*/
/*----------------------------- Public Functions -----------------------------*/

/*----------------------------------------------------------------------------*/
/*! \fn void dump_particle_history(Grid *pGrid, Domain *pD, Output *pOut)
 *  \brief  Writes particle variables as formatted table */
void dump_particle_history(Grid *pGrid, Domain *pD, Output *pOut)
{
  FILE *fid;
  int i,j,k,n,prp,mhst;
  long p,vol_rat,*npar;
  int tot_scal_cnt,tot_aray_cnt;
  Real scal[NSCAL+MAX_USR_SCAL],**array,rho,dvol;
  char fmt[20];
  Grain *gr;
#ifdef MPI_PARALLEL
  Real my_scal[NSCAL+MAX_USR_SCAL],*sendbuf,*recvbuf;
  int err;
  long *my_npar;

  my_npar = (long*)calloc_1d_array(pGrid->partypes, sizeof(long));
  sendbuf = (Real*)calloc_1d_array( (NARAY+MAX_USR_ARAY)*pGrid->partypes,
                                                             sizeof(Real));
  recvbuf = (Real*)calloc_1d_array( (NARAY+MAX_USR_ARAY)*pGrid->partypes,
                                                             sizeof(Real));
#endif

  npar  = (long*)calloc_1d_array(pGrid->partypes, sizeof(long));
  array = (Real**)calloc_2d_array(NARAY+MAX_USR_ARAY, pGrid->partypes,
                                                             sizeof(Real));

  tot_scal_cnt = 11 + usr_scal_cnt;
  tot_aray_cnt = 12 + usr_aray_cnt;

  particle_to_grid(pGrid, pD, property_all);

/* Add a white space to the format */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%13.5e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

  for (i=1; i<tot_scal_cnt; i++) {
    scal[i] = 0.0;
  }

/*--------------------- Compute scalar history variables ---------------------*/
  scal[0] = pGrid->time;

  /* Maximum density and energy dissipation rate */
  scal[1] = 0.0;
  scal[2] = 0.0;
  scal[3] = 0.0;
  for (k=pGrid->ks; k<=pGrid->ke; k++)
  for (j=pGrid->js; j<=pGrid->je; j++)
  for (i=pGrid->is; i<=pGrid->ie; i++)
  {
    scal[1] = MAX(scal[1],expr_dpar(pGrid,i,j,k));
#ifdef FEEDBACK
    scal[2] = MAX(scal[2],pGrid->Coup[k][j][i].FBstiff);
    scal[3]+= pGrid->Coup[k][j][i].Eloss;
#endif
  }

  /* particle mass, momentum and kinetic energy */
  for(p=0; p<pGrid->nparticle; p++) {
    gr = &(pGrid->particle[p]);
    if (gr->pos == 1) /* grid particle */
    {
#ifdef FEEDBACK
      rho = pGrid->grproperty[gr->property].m;/* contribution to total mass */
#else
      rho = 1.0;                              /* contribution to total number */
#endif
      mhst = 4;
      scal[mhst] += rho;
      mhst++;
      scal[mhst] += rho*gr->v1;
      mhst++;
      scal[mhst] += rho*gr->v2;
      mhst++;
      scal[mhst] += rho*gr->v3;
      mhst++;
      scal[mhst] += 0.5*rho*SQR(gr->v1);
      mhst++;
      scal[mhst] += 0.5*rho*SQR(gr->v2);
      mhst++;
      scal[mhst] += 0.5*rho*SQR(gr->v3);

      /* Calculate the user defined history variables */
      for(n=0; n<usr_scal_cnt; n++){
        mhst++;
        scal[mhst] += (*phst_scalfun[n])(pGrid, gr);
      }
    }
  }

#ifdef MPI_PARALLEL
  for (i=1; i<tot_scal_cnt; i++)
    my_scal[i] = scal[i];

  err = MPI_Reduce(&(my_scal[1]),&(scal[1]),2,
                                 MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  err = MPI_Reduce(&(my_scal[3]),&(scal[3]),8+usr_scal_cnt,
                                 MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif

  /* Average the sums */
  vol_rat = (pD->ide - pD->ids + 1)*(pD->jde - pD->jds + 1)*
            (pD->kde - pD->kds + 1);

  if(pGrid->my_id == 0){ /* I'm the parent */

    dvol = 1.0/(double)vol_rat;

    for (i=3; i<tot_scal_cnt; i++)
      scal[i] *= dvol;
  }

/* Compute particle type dependent history variables */

  mhst = 0;

  for (i=0; i<pGrid->partypes; i++)
  {
    npar[i] = 0;

    for (j=0; j<tot_aray_cnt; j++)
      array[j][i] = 0.0;
  }

  /* average position and velocity */
  for (p=0; p<pGrid->nparticle; p++)
  {
    gr = &(pGrid->particle[p]);
    if (gr->pos == 1)
    {
      prp = gr->property;
      npar[prp] += 1;
      array[0][prp] += gr->x1;
      array[1][prp] += gr->x2;
      array[2][prp] += gr->x3;
      array[3][prp] += gr->v1;
      array[4][prp] += gr->v2;
      array[5][prp] += gr->v3;
    }
  }

#ifdef MPI_PARALLEL
  for (i=0; i<pGrid->partypes; i++)
    my_npar[i] = npar[i];

  for (j=0; j<6; j++)
  {
    n = j*pGrid->partypes;
    for (i=0; i<pGrid->partypes; i++)
      sendbuf[n+i] = array[j][i];
  }

  err = MPI_Allreduce(my_npar, npar,  pGrid->partypes,
                                      MPI_LONG,  MPI_SUM,MPI_COMM_WORLD);
  err = MPI_Allreduce(sendbuf,recvbuf,6*pGrid->partypes,
                                      MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  for (j=0; j<6; j++)
  {
    n = j*pGrid->partypes;
    for (i=0; i<pGrid->partypes; i++)
      array[j][i] = recvbuf[n+i];
  }
#endif

  for (j=0; j<6; j++) {
  for (i=0; i<pGrid->partypes; i++) {
    array[j][i] = array[j][i]/MAX(1,npar[i]);
  }}

  /* position scatter, velocity dispersion, and user defined work */
  for (p=0; p<pGrid->nparticle; p++)
  {
    gr = &(pGrid->particle[p]);
    if (gr->pos == 1)
    {
      prp = gr->property;
      array[6][prp] += SQR(gr->x1-array[0][prp]);
      array[7][prp] += SQR(gr->x2-array[1][prp]);
      array[8][prp] += SQR(gr->x3-array[2][prp]);
      array[9][prp] += SQR(gr->v1-array[3][prp]);
      array[10][prp]+= SQR(gr->v2-array[4][prp]);
      array[11][prp]+= SQR(gr->v3-array[5][prp]);
      mhst = 11;
      for (n=0; n<usr_aray_cnt; n++){
        mhst++;
        array[mhst][prp] += (*phst_arayfun[n])(pGrid, gr);
      }
    }
  }

#ifdef MPI_PARALLEL
  for (j=6; j<tot_aray_cnt; j++)
  {
    n = (j-6)*pGrid->partypes;
    for (i=0; i<pGrid->partypes; i++)
      sendbuf[n+i] = array[j][i];
  }

  err = MPI_Reduce(sendbuf,recvbuf,(tot_aray_cnt-6)*pGrid->partypes,
                                   MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (pGrid->my_id == 0)
  {
    for (j=6; j<tot_aray_cnt; j++)
    {
      n = (j-6)*pGrid->partypes;
      for (i=0; i<pGrid->partypes; i++)
        array[j][i] = recvbuf[n+i];
    }

  }
#endif

  if (pGrid->my_id == 0) {

    for (i=0; i<pGrid->partypes; i++) {
      for (j=6; j<12; j++)
        array[j][i] = sqrt(array[j][i]/MAX(1,npar[i]-1));

      for (j=12; j<tot_aray_cnt; j++)
        array[j][i] = array[j][i]/MAX(1,npar[i]);
    }

  }

/*------------------------ Output particle history file ----------------------*/

  if (pGrid->my_id == 0) {

#ifdef MPI_PARALLEL
    fid = ath_fopen("../",pGrid->outfilename,0,0,NULL,"phst","a");
#else
    fid = ath_fopen(NULL,pGrid->outfilename,0,0,NULL,"phst","a");
#endif
    if(fid == NULL){
      ath_perr(-1,"[dump_particle_history]: Unable to open the history file\n");
      return;
    }

    /* Write out column headers */
    if(pOut->num == 0){

      fprintf(fid,"#Global scalars:\n");
      mhst = 1;
      fprintf(fid,"#  [%i]=time    ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=d_max   ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=stiffmax",mhst);
      mhst++;
      fprintf(fid,"  [%i]=Edot    ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=mass    ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x1 Mom. ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x2 Mom. ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x3 Mom. ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x1-KE   ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x2-KE   ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x3-KE  ",mhst);
      for(n=0; n<usr_scal_cnt; n++){
        mhst++;
        fprintf(fid,"  [%i]=%s",mhst,usr_label_scal[n]);
      }

      fprintf(fid,"\n#Particle type dependent scalars:\n");
      mhst = 1;
      fprintf(fid,"#  [%i]=x1 avg  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x2 avg  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x3 avg  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=v1 avg  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=v2 avg  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=v3 avg  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x1 var  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x2 var  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=x3 var  ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=v1 var ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=v2 var ",mhst);
      mhst++;
      fprintf(fid,"  [%i]=v3 var ",mhst);
      for(n=0; n<usr_aray_cnt; n++){
        mhst++;
        fprintf(fid,"  [%i]=%s",mhst,usr_label_aray[n]);
      }
      fprintf(fid,"\n#\n");
    }

    /* Write data */
    for (i=0; i<tot_scal_cnt; i++) {
      fprintf(fid,fmt,scal[i]);
    }
    fprintf(fid,"\n");

    for (j=0; j<pGrid->partypes; j++) {
//      fprintf(fid,"%2d:",j);
      for (i=0; i<tot_aray_cnt; i++)
        fprintf(fid,fmt,array[i][j]);
      fprintf(fid,"\n");
   }
   fprintf(fid,"\n");

    fclose(fid);
  }

  free(npar);
  free_2d_array(array);
#ifdef MPI_PARALLEL
  free(my_npar);
  free(sendbuf);
  free(recvbuf);
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void dump_parhistory_enroll(const Parfun_t pfun, const char *label,
 *                                               const int  set)
 *  \brief Set: global variable (0) or particle type dependent variable (1)
 */
void dump_parhistory_enroll(const Parfun_t pfun, const char *label,
                                                 const int  set)
{
  if (set == 0) /* global variable */
  {
    if(usr_scal_cnt >= MAX_USR_SCAL)
      ath_error("[par_history_enroll]: MAX_USR_SCAL = %d exceeded\n",
                MAX_USR_SCAL);

    /* Copy the label string */
    if((usr_label_scal[usr_scal_cnt] = ath_strdup(label)) == NULL)
      ath_error("[par_history_enroll]: Error on sim_strdup(\"%s\")\n",label);

    /* Store the function pointer */
    phst_scalfun[usr_scal_cnt] = pfun;
    usr_scal_cnt++;
  }
  else /* particle type dependent variable */
  {
    if(usr_aray_cnt >= MAX_USR_ARAY)
      ath_error("[par_history_enroll]: MAX_USR_ARAY = %d exceeded\n",
                MAX_USR_ARAY);

    /* Copy the label string */
    if((usr_label_aray[usr_aray_cnt] = ath_strdup(label)) == NULL)
      ath_error("[par_history_enroll]: Error on sim_strdup(\"%s\")\n",label);

    /* Store the function pointer */
    phst_arayfun[usr_aray_cnt] = pfun;
    usr_aray_cnt++;
  }

  return;
}

#undef NSCAL
#undef NARAY

#undef MAX_USR_SCAL
#undef MAX_USR_ARAY

#endif /* PARTICLES */

