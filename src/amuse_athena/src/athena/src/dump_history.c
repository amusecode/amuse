#include "copyright.h"
/*============================================================================*/
/*! \file dump_history.c
 *  \brief Functions to write dumps of scalar "history" variables in a
 *  formatted table.
 *
 * PURPOSE: Functions to write dumps of scalar "history" variables in a
 *   formatted table.  "History" variables are mostly VOLUME AVERAGES, e.g.
 *   - S = \Sum_{ijk} q[k][j][i]*dx1[i]*dx2[j]*dx3[k] /
 *             \Sum_{ijk} dx1[i]*dx2[j]*dx3[k],
 *
 *   although other quantities like time and dt are also output.  To compute
 *   TOTAL, VOLUME INTEGRATED quantities, just multiply by the Domain volume,
 *   which is output in the first line of the history file.  History data is
 *   written periodically, giving the time-evolution of that quantitity.
 *   Default (hardwired) values are:
 *   - scal[0] = time
 *   - scal[1] = dt
 *   - scal[2] = mass
 *   - scal[3] = total energy
 *   - scal[4] = d*v1
 *   - scal[5] = d*v2
 *   - scal[6] = d*v3
 *   - scal[7] = 0.5*d*v1**2
 *   - scal[8] = 0.5*d*v2**2
 *   - scal[9] = 0.5*d*v3**2
 *   - scal[10] = 0.5*b1**2
 *   - scal[11] = 0.5*b2**2
 *   - scal[12] = 0.5*b3**2
 *   - scal[13] = d*Phi
 *   - scal[14+NSCALARS] = passively advected scalars
 *   - scal[15+NSCALARS] = angular momentum
 *
 * More variables can be hardwired by increasing NSCAL=number of variables, and
 * adding calculation of desired quantities below.
 *
 * Alternatively, up to MAX_USR_H_COUNT new history variables can be added using
 * dump_history_enroll() in the problem generator.
 *
 * With SMR, data is averaged over each Domain separately, and dumped to
 * separate files with the level and domain number encoded in the filename.
 * Dumps are always made for all levels and domains, and are written in lev#
 * directories of the root (rank=0) process FOR THAT DOMAIN COMMUNICATOR.
 *
 * With SR, the default (hardwired) variables are:
 *   - scal[0] = time
 *   - scal[1] = dt
 *   - scal[2] = mass
 *   - scal[3] = total energy
 *   - scal[4] = x1-mom
 *   - scal[5] = x2-mom
 *   - scal[6] = x3-mom
 *   - scal[7] = (U^t)**2
 *   - scal[8] = (U^x)**2
 *   - scal[9] = (U^y)**2
 *   - scal[10] = (U^z)**2
 *   - scal[11] = (b^t)**2
 *   - scal[12] = (b^x)**2
 *   - scal[13] = (b^y)**2
 *   - scal[14] = (b^z)**2
 *   - scal[15] = |b|**2
 *   - scal[16] = T^00_EM
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - dump_history()        - Writes variables as formatted table
 * - dump_history_enroll() - Adds new user-defined history variables	      */
/*============================================================================*/

#include <stdio.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* Maximum Number of default history dump columns. */
#define NSCAL 14

/* Maximum number of history dump columns that the user routine can add. */
#define MAX_USR_H_COUNT 30

/* Array of strings / labels for the user added history column. */
static char *usr_label[MAX_USR_H_COUNT];

/* Array of history dump function pointers for user added history columns. */
static ConsFun_t phst_fun[MAX_USR_H_COUNT];

static int usr_hst_cnt = 0; /* User History Counter <= MAX_USR_H_COUNT */

/*----------------------------------------------------------------------------*/
/*! \fn void dump_history(MeshS *pM, OutputS *pOut)
 *  \brief Function to write dumps of scalar "history" variables in a
 *  formatted table. */

void dump_history(MeshS *pM, OutputS *pOut)
{
  GridS *pG;
  DomainS *pD;
  int i,j,k,is,ie,js,je,ks,ke,nl,nd;
  double dVol, scal[NSCAL + NSCALARS + MAX_USR_H_COUNT], d1;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL,*pdir=NULL,fmt[80];
  char levstr[8],domstr[8],dirstr[20];
  int n, total_hst_cnt, mhst, myID_Comm_Domain=1;
#ifdef MPI_PARALLEL
  double my_scal[NSCAL + NSCALARS + MAX_USR_H_COUNT]; /* My Volume averages */
  int ierr;
#endif
#ifdef CYLINDRICAL
  Real x1,x2,x3;
#endif
#ifdef SPECIAL_RELATIVITY
  PrimS W;
  Real g, g2, g_2;
  Real bx, by, bz, vB, b2, Bmag2;
#endif


  total_hst_cnt = 9 + NSCALARS + usr_hst_cnt;
#ifdef ADIABATIC
  total_hst_cnt++;
#endif
#ifdef MHD
  total_hst_cnt += 3;
#endif
#ifdef SELF_GRAVITY
  total_hst_cnt += 1;
#endif
#ifdef CYLINDRICAL
  total_hst_cnt++;  /* for angular momentum */
#endif
#ifdef SPECIAL_RELATIVITY
  total_hst_cnt = 12 + usr_hst_cnt;
#ifdef MHD
   total_hst_cnt += 6;
#endif
#endif


/* Add a white space to the format */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%14.6e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

/* store time and dt in first two elements of output vector */

  scal[0] = pM->time;
  scal[1] = pM->dt;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){

    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        is = pG->is, ie = pG->ie;
        js = pG->js, je = pG->je;
        ks = pG->ks, ke = pG->ke;

        for (i=2; i<total_hst_cnt; i++) {
          scal[i] = 0.0;
        }
 
/* Compute history variables */

        for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
              dVol = 1.0; 
              if (pG->dx1 > 0.0) dVol *= pG->dx1;
              if (pG->dx2 > 0.0) dVol *= pG->dx2;
              if (pG->dx3 > 0.0) dVol *= pG->dx3;
#ifndef SPECIAL_RELATIVITY
#ifdef CYLINDRICAL
              cc_pos(pG,i,j,k,&x1,&x2,&x3);
              dVol *= x1;
#endif

              mhst = 2;
              scal[mhst] += dVol*pG->U[k][j][i].d;
              d1 = 1.0/pG->U[k][j][i].d;
#ifndef BAROTROPIC
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].E;
#endif
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].M1;
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].M2;
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].M3;
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].M1)*d1;
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].M2)*d1;
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].M3)*d1;
#ifdef MHD
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B1c);
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B2c);
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B3c);
#endif
#ifdef SELF_GRAVITY
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].d*pG->Phi[k][j][i];
#endif
#if (NSCALARS > 0)
              for(n=0; n<NSCALARS; n++){
                mhst++;
                scal[mhst] += dVol*pG->U[k][j][i].s[n];
              }
#endif

#ifdef CYLINDRICAL
              mhst++;
              scal[mhst] += dVol*(x1*pG->U[k][j][i].M2);
#endif

#else /* SPECIAL_RELATIVITY */

              W = Cons_to_Prim (&(pG->U[k][j][i]));
        
              /* calculate gamma */
              g   = pG->U[k][j][i].d/W.d;
              g2  = SQR(g);
              g_2 = 1.0/g2;

              mhst = 2;
              scal[mhst] += dVol*pG->U[k][j][i].d;
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].E;
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].M1;
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].M2;
              mhst++;
              scal[mhst] += dVol*pG->U[k][j][i].M3;

              mhst++;
              scal[mhst] += dVol*SQR(g);
              mhst++;
              scal[mhst] += dVol*SQR(g*W.V1);
              mhst++;
              scal[mhst] += dVol*SQR(g*W.V2);
              mhst++;
              scal[mhst] += dVol*SQR(g*W.V3);

	      mhst++;
	      scal[mhst] += dVol*W.P;

#ifdef MHD

              vB = W.V1*pG->U[k][j][i].B1c + W.V2*W.B2c + W.V3*W.B3c;
              Bmag2 = SQR(pG->U[k][j][i].B1c) + SQR(W.B2c) + SQR(W.B3c);
        
              bx = g*(pG->U[k][j][i].B1c*g_2 + vB*W.V1);
              by = g*(W.B2c*g_2 + vB*W.V2);
              bz = g*(W.B3c*g_2 + vB*W.V3);
        
              b2 = Bmag2*g_2 + vB*vB;

              mhst++;
              scal[mhst] += dVol*(g*vB*g*vB);
              mhst++;
              scal[mhst] += dVol*bx*bx;
              mhst++;
              scal[mhst] += dVol*by*by;
              mhst++;
              scal[mhst] += dVol*bz*bz;
              mhst++;
              scal[mhst] += dVol*b2;
              mhst++;
              scal[mhst] += dVol*(Bmag2*(1.0 - 0.5*g_2) - SQR(vB) / 2.0);

#endif /* MHD */

#endif  /* SPECIAL_RELATIVITY */

/* Calculate the user defined history variables */
              for(n=0; n<usr_hst_cnt; n++){
                mhst++;
                scal[mhst] += dVol*(*phst_fun[n])(pG, i, j, k);
              }
            }
          }
        }

/* Compute the sum over all Grids in Domain */

#ifdef MPI_PARALLEL 
        for(i=2; i<total_hst_cnt; i++){
          my_scal[i] = scal[i];
        }
        ierr = MPI_Reduce(&(my_scal[2]), &(scal[2]), (total_hst_cnt - 2),
          MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#endif

/* Only the parent (rank=0) process computes the average and writes output.
 * For single-processor jobs, myID_Comm_world is always zero. */

#ifdef MPI_PARALLEL
        ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
#endif
        if((myID_Comm_Domain==0) || (myID_Comm_world==0)){  /* I'm the parent */

/* Compute volume averages */

          dVol = pD->MaxX[0] - pD->MinX[0];
#ifdef CYLINDRICAL
          dVol = (pD->MaxX[0]*pD->MaxX[0]) - (pD->MinX[0]*pD->MinX[0]);
#endif
          if (pD->Nx[1] > 1) dVol *= (pD->MaxX[1] - pD->MinX[1]);
          if (pD->Nx[2] > 1) dVol *= (pD->MaxX[2] - pD->MinX[2]);
          for(i=2; i<total_hst_cnt; i++){
            scal[i] /= dVol;
          }

/* Create filename and open file.  History files are always written in lev#
 * directories of root process (rank=0 in MPI_COMM_WORLD) */
#ifdef MPI_PARALLEL
          if (nl>0) {
            plev = &levstr[0];
            sprintf(plev,"lev%d",nl);
            pdir = &dirstr[0];
            sprintf(pdir,"../id0/lev%d",nl);
          }
#else
          if (nl>0) {
            plev = &levstr[0];
            sprintf(plev,"lev%d",nl);
            pdir = &dirstr[0];
            sprintf(pdir,"lev%d",nl);
          }
#endif

          if (nd>0) {
            pdom = &domstr[0];
            sprintf(pdom,"dom%d",nd);
          }

          fname = ath_fname(pdir,pM->outfilename,plev,pdom,0,0,NULL,"hst");
          if(fname == NULL){
            ath_perr(-1,"[dump_history]: Unable to create history filename\n");
          }
          pfile = fopen(fname,"a");
          if(pfile == NULL){
            ath_perr(-1,"[dump_history]: Unable to open the history file\n");
          }

/* Write out column headers, but only for first dump */

          mhst = 0;
          if(pOut->num == 0){
            fprintf(pfile,
         "# Athena history dump for level=%i domain=%i volume=%e\n",nl,nd,dVol);
            mhst++;
            fprintf(pfile,"#   [%i]=time   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=dt      ",mhst);
#ifndef SPECIAL_RELATIVITY
            mhst++;
            fprintf(pfile,"   [%i]=mass    ",mhst);
#ifdef ADIABATIC
            mhst++;
            fprintf(pfile,"   [%i]=total E ",mhst);
#endif
            mhst++;
            fprintf(pfile,"   [%i]=x1 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1-KE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2-KE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3-KE   ",mhst);
#ifdef MHD
            mhst++;
            fprintf(pfile,"   [%i]=x1-ME   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2-ME   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3-ME   ",mhst);
#endif
#ifdef SELF_GRAVITY
            mhst++;
            fprintf(pfile,"   [%i]=grav PE ",mhst);
#endif
#if (NSCALARS > 0)
            for(n=0; n<NSCALARS; n++){
              mhst++;
              fprintf(pfile,"  [%i]=scalar %i",mhst,n);
            }
#endif

#ifdef CYLINDRICAL
            mhst++;
            fprintf(pfile,"   [%i]=Ang.Mom.",mhst);
#endif

#else /* SPECIAL_RELATIVITY */
            mhst++;
            fprintf(pfile,"   [%i]=mass    ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=total E ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3 Mom." ,mhst);
            mhst++;
            fprintf(pfile,"   [%i]=Gamma   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1-KE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2-KE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3-KE  " ,mhst);
            mhst++;
            fprintf(pfile,"   [%i]=Press  " ,mhst);
#ifdef MHD
            mhst++;
            fprintf(pfile,"   [%i]=x0-ME  " ,mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1-ME  " ,mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2-ME  " ,mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3-ME  " ,mhst);
            mhst++;
            fprintf(pfile,"   [%i]=bsq    " ,mhst);
            mhst++;
            fprintf(pfile,"   [%i]=T^00_EM" ,mhst);
#endif
#endif /* SPECIAL_RELATIVITY */

            for(n=0; n<usr_hst_cnt; n++){
              mhst++;
              fprintf(pfile,"  [%i]=%s",mhst,usr_label[n]);
            }
            fprintf(pfile,"\n#\n");
          }

/* Write out data, and close file */

          for (i=0; i<total_hst_cnt; i++) {
            fprintf(pfile,fmt,scal[i]);
          }
          fprintf(pfile,"\n");
          fclose(pfile);
  
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void dump_history_enroll(const ConsFun_t pfun, const char *label)
 *  \brief Adds new user-defined history variables	      */

void dump_history_enroll(const ConsFun_t pfun, const char *label){

  if(usr_hst_cnt >= MAX_USR_H_COUNT)
    ath_error("[dump_history_enroll]: MAX_USR_H_COUNT = %d exceeded\n",
	      MAX_USR_H_COUNT);

/* Copy the label string */
  if((usr_label[usr_hst_cnt] = ath_strdup(label)) == NULL)
    ath_error("[dump_history_enroll]: Error on sim_strdup(\"%s\")\n",label);

/* Store the function pointer */
  phst_fun[usr_hst_cnt] = pfun;

  usr_hst_cnt++;

  return;

}

#undef NSCAL
#undef MAX_USR_H_COUNT
