#include <stdio.h>
#include <string.h>
#include "checkpoint.h"
#include "init.h"
#include "userFunc.h"

/**************************************************************************/
/* Read and write checkpoint routines                                     */
/**************************************************************************/

iostatus readCheckpoint(
			/* Checkpoint file name and number */
			const char *checkname,
			unsigned long *checknum,
			/* Control parameters */
			bool *eos_func, bool *massSrc_func,
			bool *intEnSrc_func,
			/* Output control parameters */
			unsigned long *nUserOut, 
			/* Simulation state parameters */
			double *t, double *dt,
			/* Diagnostic outputs */
			unsigned long *nStep,
			unsigned long *nIter,
			unsigned long *nFail,
			/* Grid and workspace */
			grid **grd, wksp **w,
			/* Simulation outputs */
			unsigned long *nOut,
			double **tOut, double **colOut,
			double **presOut, double **eIntOut,
			double **mBndOut, double **eBndOut,
			double **mSrcOut, double **eSrcOut,
			double **userOut,
			/* Control parameters */
			void *params,
			const bool userReadCheckpoint,
			const unsigned long verbosity) {

  unsigned long nr;
  FILE *fp;
  double *eInt_tmp, *mSrc_tmp, *eSrc_tmp, *user_tmp;

  /* Open file; bail if opening fails */
  if (!(fp = fopen(checkname, "r"))) {
    fprintf(stderr, "vader: error: unable to open checkpoint file %s!\n",
	    checkname);
    return IO_ERROR;
  }

  /* Read the checkpoint sequence number */
  fread(checknum, sizeof(unsigned long), 1, fp);
   
  /* Read control parameters */
  fread(eos_func, sizeof(bool), 1, fp);    
  fread(massSrc_func, sizeof(bool), 1, fp);   
  fread(intEnSrc_func, sizeof(bool), 1, fp);
  fread(nUserOut, sizeof(unsigned long), 1, fp);

  /* Read step info */
  fread(t, sizeof(double), 1, fp);
  fread(dt, sizeof(double), 1, fp);
  fread(nStep, sizeof(unsigned long), 1, fp);
  fread(nIter, sizeof(unsigned long), 1, fp);
  fread(nFail, sizeof(unsigned long), 1, fp);

  /* Read size of grid */
  nr = 0;
  fread(&nr, sizeof(unsigned long), 1, fp);

  /* Allocate memory to hold the grid and workspace */
  if (nr > 0) {
    *grd = gridAlloc(nr);
    *w = wkspAlloc(nr);
  } else {
    fclose(fp);
    fprintf(stderr, "vader: error: bad checkpoint file %s, grid size = %lu!\n",
	    checkname, nr);
    return IO_ERROR;
  }

  /* Read grid data */
  fread(&((*grd)->linear), sizeof(bool), 1, fp);    
  fread((*grd)->r_g, sizeof(double), (*grd)->nr+2, fp); 
  fread((*grd)->r_h, sizeof(double), (*grd)->nr+1, fp);
  fread((*grd)->dr_g, sizeof(double), (*grd)->nr+2, fp);
  fread((*grd)->area, sizeof(double), (*grd)->nr, fp);
  fread((*grd)->vphi_g, sizeof(double), (*grd)->nr+2, fp);
  fread((*grd)->vphi_h, sizeof(double), (*grd)->nr+1, fp);
  fread((*grd)->beta_g, sizeof(double), (*grd)->nr+2, fp);
  fread((*grd)->beta_h, sizeof(double), (*grd)->nr+1, fp);
  fread((*grd)->psiEff_g, sizeof(double), (*grd)->nr+2, fp);
  fread((*grd)->psiEff_h, sizeof(double), (*grd)->nr+1, fp);
  fread((*grd)->g_h, sizeof(double), (*grd)->nr+1, fp);

  /* Allocate memory to hold output data stored in checkpoint */
  fread(nOut, sizeof(unsigned long), 1, fp);
  if (!outputAlloc(*nOut, *eos_func, *massSrc_func, *intEnSrc_func,
		   *nUserOut, *grd, tOut, colOut, presOut,
		   eIntOut, mBndOut, eBndOut, mSrcOut,
		   eSrcOut, userOut)) {
    fclose(fp);
    gridFree(*grd);
    wkspFree(*w);
    grd = NULL;
    w = NULL;
    fprintf(stderr,
	    "vader: error: unable to allocate memory to read checkpoint %s!\n",
	    checkname);
    return ALLOCATION_ERROR;
  }

  /* Read checkpoint data */
  fread(*tOut, sizeof(double), *nOut, fp);
  fread(*colOut, sizeof(double), (*grd)->nr*(*nOut), fp);
  fread(*presOut, sizeof(double), (*grd)->nr*(*nOut), fp);
  if (*eos_func)
    fread(*eIntOut, sizeof(double), (*grd)->nr*(*nOut), fp);
  fread(*mBndOut, sizeof(double), 2*(*nOut), fp);
  fread(*eBndOut, sizeof(double), 2*(*nOut), fp);
  if (*massSrc_func)
    fread(*mSrcOut, sizeof(double), (*grd)->nr*(*nOut), fp);
  if (*massSrc_func || *intEnSrc_func)
    fread(*eSrcOut, sizeof(double), (*grd)->nr*(*nOut), fp);

  /* Read user outputs */
  if (*nUserOut > 0)
    fread(*userOut, sizeof(double), (*grd)->nr*(*nOut)*(*nUserOut), fp);

  /* Read the Anderson acceleration contraint vector */
#if AA_M > 0
  gsl_vector_fread(fp, (*w)->constraint);
#endif

  /* Call the user-provided reading routine if requested */
  if (userReadCheckpoint) {
    if (*eos_func) eInt_tmp = *eIntOut;
    else eInt_tmp = NULL;
    if (*massSrc_func) mSrc_tmp = *mSrcOut;
    else mSrc_tmp = NULL;
    if (*massSrc_func || *intEnSrc_func) eSrc_tmp = *eSrcOut;
    else eSrc_tmp = NULL;
    if (*nUserOut > 0) user_tmp = *userOut;
    else user_tmp = NULL;
    userCheckRead(fp, *grd, *nOut, *tOut, *colOut, *presOut,
		  eInt_tmp, *mBndOut, *eBndOut, mSrc_tmp, eSrc_tmp,
		  *nUserOut, user_tmp, params);
  }

  /* Make sure there are no errors, and bail if there are */
  if (ferror(fp)) {
    fclose(fp);
    gridFree(*grd);
    wkspFree(*w);
    outputFree(tOut, colOut, presOut, eIntOut, mBndOut, eBndOut,
	       mSrcOut, eSrcOut, userOut);
    *grd = NULL;
    *w = NULL;
    fprintf(stderr,
	    "vader: error: checkpoint file %s appears to be damaged!\n",
	    checkname);
    return IO_ERROR;
  }
  
  /* Close file, print status, return */
  fclose(fp);
  if (verbosity > 0)
    printf("Successfully read checkpoint file %s\n", checkname);

  return GOOD_IO;
}



iostatus saveCheckpoint(
			/* Checkpoint file name and number */
			const char *checkname,
			const unsigned long checknum,
			/* Control parameters */
			const bool eos_func,
			const bool massSrc_func,
			const bool intEnSrc_func,
			/* Output control parameters */
			const unsigned long nUserOut,
			/* Simulation state parameters */
			const double t,
			const double dt,
			/* Diagnostic outputs */
			const unsigned long nStep,
			const unsigned long nIter,
			const unsigned long nFail,
			/* The grid */
			const grid *grd,
			/* Simulation outputs */
			const unsigned long nOut,
			const double *tOut,
			const double *colOut,
			const double *presOut, 
			const double *eIntOut,
			const double *mBndOut,
			const double *eBndOut,
			const double *mSrcOut,
			const double *eSrcOut,
			const double *userOut,
#if AA_M > 0
			const gsl_vector *constraint,
#endif
			/* Control parameters */
			const void *params,
			const bool userWriteCheckpoint,
			const unsigned long verbosity) {
  char *fname;
  FILE *fp;
  
  /* Construct the name of the checkpoint file */
  fname = calloc(strlen(checkname) + 13, sizeof(char));
  sprintf(fname, "%s_%05lu.vader", checkname, checknum);

  /* Open the file; bail if we fail */
  if (!(fp = fopen(fname, "w"))) {
    fprintf(stderr, "vader: warning: unable to open checkpoint file %s!\n",
	    fname);
    return IO_ERROR;
  }
  
  /* Write the checkpoint file sequence number */
  fwrite(&checknum, sizeof(unsigned long), 1, fp);

  /* Write out the control parameters that affect what's in the
     checkpoint */
  fwrite(&eos_func, sizeof(bool), 1, fp);
  fwrite(&massSrc_func, sizeof(bool), 1, fp);
  fwrite(&intEnSrc_func, sizeof(bool), 1, fp);
  fwrite(&nUserOut, sizeof(unsigned long), 1, fp);

  /* Write out the step number and related information */
  fwrite(&t, sizeof(double), 1, fp);
  fwrite(&dt, sizeof(double), 1, fp);
  fwrite(&nStep, sizeof(unsigned long), 1, fp);
  fwrite(&nIter, sizeof(unsigned long), 1, fp);
  fwrite(&nFail, sizeof(unsigned long), 1, fp);

  /* Write out the grid */
  fwrite(&(grd->nr), sizeof(unsigned long), 1, fp);
  fwrite(&(grd->linear), sizeof(bool), 1, fp);
  fwrite(grd->r_g, sizeof(double), grd->nr+2, fp);
  fwrite(grd->r_h, sizeof(double), grd->nr+1, fp);
  fwrite(grd->dr_g, sizeof(double), grd->nr+2, fp);
  fwrite(grd->area, sizeof(double), grd->nr, fp);
  fwrite(grd->vphi_g, sizeof(double), grd->nr+2, fp);
  fwrite(grd->vphi_h, sizeof(double), grd->nr+1, fp);
  fwrite(grd->beta_g, sizeof(double), grd->nr+2, fp);
  fwrite(grd->beta_h, sizeof(double), grd->nr+1, fp);
  fwrite(grd->psiEff_g, sizeof(double), grd->nr+2, fp);
  fwrite(grd->psiEff_h, sizeof(double), grd->nr+1, fp);
  fwrite(grd->g_h, sizeof(double), grd->nr+1, fp);

  /* Write the simulation outputs */
  fwrite(&nOut, sizeof(unsigned long), 1, fp);
  fwrite(tOut, sizeof(double), nOut, fp);
  
  /* Write out the output holders; note that we write all of them,
     which results in some redundancy if we have multiple
     checkpoints, but this is worth because having the full history
     in each checkpoint makes life much easier, and the data aren't
     going to be that large most of the time */
  fwrite(colOut, sizeof(double), grd->nr*nOut, fp);
  fwrite(presOut, sizeof(double), grd->nr*nOut, fp);
  if (eos_func)
    fwrite(eIntOut, sizeof(double), grd->nr*nOut, fp);
  fwrite(mBndOut, sizeof(double), 2*nOut, fp);
  fwrite(eBndOut, sizeof(double), 2*nOut, fp);
  if (massSrc_func)
    fwrite(mSrcOut, sizeof(double), grd->nr*nOut, fp);
  if (massSrc_func || intEnSrc_func)
    fwrite(eSrcOut, sizeof(double), grd->nr*nOut, fp);

  /* Write out user outputs */
  if (nUserOut > 0)
    fwrite(userOut, sizeof(double), grd->nr*nOut*nUserOut, fp);

  /* Write the constraint vector if using Anderson acceleration */
#if AA_M > 0
  gsl_vector_fwrite(fp, constraint);
#endif

  /* Call the user-provided writing routine if requested */
  if (userWriteCheckpoint) {
    userCheckWrite(fp, grd, nOut, tOut, colOut, presOut, eIntOut,
		   mBndOut, eBndOut, mSrcOut, eSrcOut,
		   nUserOut, userOut, params);
  }
  
  /* Close */
  fclose(fp);

  /* Print status if verbose */
  if (verbosity > 0)
    printf("Wrote checkpoint file %s\n", fname);

  /* Return success */
  return GOOD_IO;
}
