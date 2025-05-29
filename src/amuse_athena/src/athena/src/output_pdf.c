#include "copyright.h"
/*============================================================================*/
/*! \file output_pdf.c
 *  \brief Outputs Probability Distribution Functions of selected variables
 *   in formatted tabular form.  
 *
 * PURPOSE: Outputs Probability Distribution Functions of selected variables
 *   in formatted tabular form.  Fully MPI enabled, which requires passing
 *   lots of global sums and means (only the parent process produces output).
 *   With SMR, dumps are made for all levels and domains, unless nlevel and
 *   ndomain are specified in <output> block.

 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - output_pdf() - output PDFs
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

static int size_dat=0; /* Number of elements in the data[] array */
static double *data=NULL; /* Computed data array: data[size_dat] */

static int size_pdf=0; /* Number of elements in the pdf[] array */
static int *pdf=NULL; /* (non-normalized) PDF */
#ifdef MPI_PARALLEL
static int *cd_pdf=NULL; /* (non-normalized) complete Domain PDF */
#endif /* MPI_PARALLEL */

static char def_fmt[]="%21.15e"; /* A default tabular dump data format */

/*----------------------------------------------------------------------------*/
/*! \fn void output_pdf(MeshS *pM, OutputS *pOut)
 *  \brief Outputs PDFs. */

void output_pdf(MeshS *pM, OutputS *pOut)
{
  GridS *pG;
  FILE *pfile;
  char fmt[80];
  char fid[80]; /* File "id" for the statistics table */
  char *fname,*plev=NULL,*pdom=NULL,*pdir=NULL;
  char levstr[8],domstr[8],dirstr[8];
  int nl,nd,i,j,k,is,ie,js,je,ks,ke;
  int n, data_cnt;
  double dmin, dmax, delta, dpdf, dat, scl;
  double mean=0.0, var=0.0; /* mean and variance of the distribution */
  double adev=0.0, sdev=0.0; /* average & standard deviation */
  double skew=0.0, kurt=0.0; /* skewness and kurtosis of the distribution */
  double r, s, ep=0.0; /* Temp. variables for calculating the variance, etc. */
#ifdef MPI_PARALLEL
  DomainS *pD;
  int ierr, cd_data_cnt, myID_Comm_Domain;
#endif /* MPI_PARALLEL */

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pG = pM->Domain[nl][nd].Grid;

#ifdef MPI_PARALLEL
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
        cd_data_cnt = (pD->Nx[0])*(pD->Nx[1])*(pD->Nx[2]);
#endif
        is = pG->is, ie = pG->ie;
        js = pG->js, je = pG->je;
        ks = pG->ks, ke = pG->ke;
        data_cnt = (ie - is + 1)*(je - js + 1)*(ke - ks + 1);

/* Are the requisite arrays allocated? */
        if(data == NULL){
          size_dat = data_cnt;
          data = (double *)calloc(size_dat,sizeof(double));
          if(data == NULL)
            ath_error("[output_pdf]: Failed to allocate data array\n");

/* This choice for size_pdf represents a balance between
 * resolution in the PDF and "shot noise" in the data binning. */
#ifdef MPI_PARALLEL
          size_pdf = (int)sqrt((double)cd_data_cnt);
#else /* MPI_PARALLEL */
          size_pdf = (int)sqrt((double)size_dat);
#endif /* MPI_PARALLEL */
          pdf = (int *)calloc(size_pdf,sizeof(int));
          if(pdf == NULL)
            ath_error("[output_pdf]: Failed to allocate pdf array\n");

#ifdef MPI_PARALLEL
          if(myID_Comm_Domain == 0){ /* I'm the parent */
            cd_pdf = (int *)calloc(size_pdf,sizeof(int));
            if(cd_pdf == NULL)
              ath_error("[output_pdf]: Failed to allocate cd_pdf array\n");
          }
#endif /* MPI_PARALLEL */
        }


/* Initialize dmin, dmax */
        dmin = dmax = (*pOut->expr)(pG,is,js,ks);

/* Fill the data array */
        n=0;
        for(k = ks; k<=ke; k++){
          for(j = js; j<=je; j++){
            for(i = is; i<=ie; i++){
              data[n] = (double)(*pOut->expr)(pG,i,j,k);
              dmin = data[n] < dmin ? data[n] : dmin;
              dmax = data[n] > dmax ? data[n] : dmax;
              mean += data[n];
              n++;
            }
          }
        }

#ifdef MPI_PARALLEL

        dat = dmin;
        ierr = MPI_Allreduce(&dat,&dmin,1,MPI_DOUBLE,MPI_MIN,pD->Comm_Domain);

        dat = dmax;
        ierr = MPI_Allreduce(&dat,&dmax,1,MPI_DOUBLE,MPI_MAX,pD->Comm_Domain);

        dat = mean;
        ierr = MPI_Allreduce(&dat,&mean,1,MPI_DOUBLE,MPI_SUM,pD->Comm_Domain);

        mean /= (double)cd_data_cnt; /* Complete the calc. of the mean */

#else /* MPI_PARALLEL */

        mean /= (double)data_cnt; /* Complete the calc. of the mean */

#endif /* MPI_PARALLEL */

        if(data_cnt > 1){
/* Calculate the variance, etc. with the corrected 2-pass formula */
          for(n=0; n<data_cnt; n++){
            s = data[n] - mean;
            adev += fabs(s);
            ep += s;
            var += (r = s*s);
            skew += (r *= s);
            kurt += (r *= s);
          }

#ifdef MPI_PARALLEL
          dat = ep;
          ierr = MPI_Allreduce(&dat,&ep,1,MPI_DOUBLE,MPI_SUM,pD->Comm_Domain);

          dat = var;
          ierr = MPI_Allreduce(&dat,&var,1,MPI_DOUBLE,MPI_SUM,pD->Comm_Domain);

          dat = skew;
          ierr = MPI_Allreduce(&dat,&skew,1,MPI_DOUBLE,MPI_SUM,pD->Comm_Domain);

          dat = kurt;
          ierr = MPI_Allreduce(&dat,&kurt,1,MPI_DOUBLE,MPI_SUM,pD->Comm_Domain);

          adev /= (double)cd_data_cnt;
          var = (var - ep*ep/(double)cd_data_cnt)/(double)(cd_data_cnt-1);
          sdev = sqrt(var);
          if(sdev > 0.0){
            skew /= var*sdev*cd_data_cnt;
            kurt = kurt/(var*var*cd_data_cnt) - 3.0;
          }
#else /* MPI_PARALLEL */
          adev /= (double)data_cnt;
          var = (var - ep*ep/(double)data_cnt)/(double)(data_cnt-1);
          sdev = sqrt(var);
          if(sdev > 0.0){
            skew /= var*sdev*data_cnt;
            kurt = kurt/(var*var*data_cnt) - 3.0;
          }
#endif /* MPI_PARALLEL */
        }

/* Store the global maximum and minimum of the quantity */
        pOut->gmin = dmin < pOut->gmin ? dmin : pOut->gmin;
        pOut->gmax = dmax > pOut->gmax ? dmax : pOut->gmax;

/* Compute the pdf directly using sampling. Define size_pdf bins, each of equal
 * size, and fill them with the number of cells whose data value falls in the
 * range spanned by the bin. */

        if(dmax - dmin > 0.0){
/* Initialize pdf[] to zero */
          for(n=0; n<size_pdf; n++) pdf[n] = 0;
/* Calculate the number of cells whose data falls in each bin */
          scl = (double)size_pdf/(dmax - dmin);
          for(n=0; n<data_cnt; n++){
            i = (int)(scl*(data[n] - dmin));
            i = i < size_pdf ? i : size_pdf - 1;
            pdf[i]++;
          }
        }

#ifdef MPI_PARALLEL

/* Sum up the pdf in the array cd_pdf */
        ierr=MPI_Reduce(pdf,cd_pdf,size_pdf,MPI_INT,MPI_SUM,0,pD->Comm_Domain);

#endif /* MPI_PARALLEL */

#ifdef MPI_PARALLEL
/* For parallel calculations, only the parent writes the output. */
        if(myID_Comm_Domain != 0) return;
#endif /* MPI_PARALLEL */

/* Create filename and open file.  pdf files are always written in lev#
 * directories of root (rank=0) process. */
#ifdef MPI_PARALLEL
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
          pdir = &dirstr[0];
          sprintf(pdir,"../lev%d",nl);
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

        fname = ath_fname(pdir,pM->outfilename,plev,pdom,num_digit,
          pOut->num,pOut->id,"prb");
        if(fname == NULL){
          ath_perr(-1,"[output_pdf]: Unable to create filename\n");
        }
        pfile = fopen(fname,"w");
        if(pfile == NULL){
          ath_perr(-1,"[output_pdf]: Unable to open pdf file\n");
        }

/* Write out some extra information in a header */
        fprintf(pfile,"# Time = %21.15e\n",pG->time);
        fprintf(pfile,"# expr = \"%s\"\n",pOut->out);
        fprintf(pfile,"# Nbin = %d\n",((dmax - dmin) > 0.0 ? size_pdf : 1));
        fprintf(pfile,"# dmin = %21.15e\n",dmin); 
        fprintf(pfile,"# dmax = %21.15e\n",dmax); 
        fprintf(pfile,"# mean = %21.15e\n",mean); 
        fprintf(pfile,"# variance = %21.15e\n",var); 
        fprintf(pfile,"# std. dev. = %21.15e\n",sdev); 
        fprintf(pfile,"# avg. dev. = %21.15e\n",adev); 
        fprintf(pfile,"# skewness = %21.15e\n",skew); 
        fprintf(pfile,"# kurtosis = %21.15e\n#\n",kurt); 

/* Add a white space to the format */
        if(pOut->dat_fmt == NULL)
          sprintf(fmt,"%s  %s\n",def_fmt, def_fmt);
        else
          sprintf(fmt,"%s  %s\n",pOut->dat_fmt,pOut->dat_fmt);

/* write out the normalized Proabability Distribution Function */
        if(dmax - dmin > 0.0){
          delta = (dmax - dmin)/(double)(size_pdf);
#ifdef MPI_PARALLEL
          scl = (double)size_pdf/(double)(cd_data_cnt*(dmax - dmin));
#else
          scl = (double)size_pdf/(double)(data_cnt*(dmax - dmin));
#endif /* MPI_PARALLEL */
          for(n=0; n<size_pdf; n++){
/* Calculate the normalized Prob. Dist. Fun. */
            dat = dmin + (n + 0.5)*delta;
#ifdef MPI_PARALLEL
            dpdf = (double)(cd_pdf[n])*scl;
#else
            dpdf = (double)(pdf[n])*scl;
#endif /* MPI_PARALLEL */
            fprintf(pfile, fmt, dat, dpdf);
          }
        }
        else
          fprintf(pfile,fmt,dmax,1.0);

        fclose(pfile);

/* Also write a history type file on the statistics */
        sprintf(fid,"prb_stat.%s",pOut->id);

        fname = ath_fname(pdir,pM->outfilename,plev,pdom,0,0,fid,"tab");
        if(fname == NULL){
          ath_perr(-1,"[output_pdf]: Unable to create stats filename\n");
        }
        pfile = fopen(fname,"a");
        if(pfile == NULL){
          ath_perr(-1,"[output_pdf]: Unable to open stats file\n");
        }

        if(pOut->num == 0){
          fprintf(pfile,"# expr = \"%s\"\n#\n",pOut->out);
          fprintf(pfile,"# time  dmin  dmax  mean  variance  \"std. dev.\"  ");
          fprintf(pfile,"\"avg. dev.\"  skewness  kurtosis\n#\n");
        }

/* Add a white space to the format */
        if(pOut->dat_fmt == NULL) sprintf(fmt," %s",def_fmt);
        else                      sprintf(fmt," %s",pOut->dat_fmt);

        fprintf(pfile,"%21.15e",pG->time);

/* write out the table of statistics */
        fprintf(pfile,fmt,dmin);
        fprintf(pfile,fmt,dmax);
        fprintf(pfile,fmt,mean);
        fprintf(pfile,fmt,var);
        fprintf(pfile,fmt,sdev);
        fprintf(pfile,fmt,adev);
        fprintf(pfile,fmt,skew);
        fprintf(pfile,fmt,kurt);
        fprintf(pfile,"\n");

        fclose(pfile);
      }}
    }
  }

  return;
}
