/* HOP_INPUT.C, Daniel Eisenstein, 1997 */
/* Based on a paper by Daniel Eisenstein & Piet Hut, 
"HOP: A New Group-Finding Algorithm for N-body Simulations."
See the included documentation or view it at 
http://www.sns.ias.edu/~eisenste/hop/hop_doc.html */

/* Version 1.0 (12/15/97) -- Original Release */

/* The routine ReadSimulationFile() is just a wrapper for whatever
routine you need to get your simulation format read from disk and
put in the KD format.  I've included three examples, ReadSimple(),
ReadASCII(), and ReadTPM(), but you can do what you like. */

/* Since you will probably need to write a new version of this,
here's what has to happen:

You are reading from the file fp (which might be stdin and need not be
opened or closed) and putting data into the structure kd.

1) Set kd->nActive:

    kd->nActive = The number of particles you intend to run the algorithm on.

2) Initialize space for kd->p, an array to hold the information for all the
particles.  Do this by:

    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));

3) Read in the particle information and put it in the kd->p array.
Note that the array is zero-offset, i.e. it runs from 0 to nActive-1.

You only need to set the position and mass of the particles.  The positions
are at:
	kd->p[j].r[0-2] -- The position information of particle j along
				Cartesian directions 0, 1, and 2.

The masses depend on whether you are compiling with the DIFFERENT_MASSES
flag.  If yes, then kd->p[j].fMass holds the mass of particle j.  If no,
then all the particles are assumed to have the same mass, equal to kd->fMass

The mass and length scales you chose are up to you, but remember that
all your density thresholds will be in those units. 

That's it! */

/* I've included two routines f77read() and f77write() at the bottom
of this file, in case you want to read or write FORTRAN's "unformatted" 
output. */

/* If you only want to run on a subset of the particles, you need to 
adjudicate that within this subroutine and make sure that kd->nActive
and kd->p contain only the particles that you want to include in the 
grouping. */

/* The following variables in the structure KD aren't used by the HOP
program.  You only need to set them if you want to use them for your
own custom purposes (e.g. changing the output to the TIPSY format):
    kd-> nDark, nGas, nStar, nParticles, fTime, bDark, bGas, bStar */

/* Sorry, but I haven't tested ReadASCII or ReadSimple since my files
aren't in that format.  They look ok by eye, but that's been known
to fail :-).  In any case, the point is to give an example of what 
has to be done. */

/* ================================================================ */
/* ====================== ReadSimulationFile() =================== */
/* ================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kd.h"


#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_INT(x)   (*(unsigned int *)&(x)   = SWAP_4(*(unsigned int *)&(x)))
#define FIX_FLOAT(x) FIX_INT(x)


int ReadSimulationFile(KD kd, char *inputfile)
{
  int ReadRamses(KD kd,char *inputfile);

  ReadRamses(kd, inputfile);

  return;
}

/* ================================================================ */
/* =================== An Example: ReadSimple() =================== */
/* ================================================================ */

/* Read the following simple format -- unformatted FORTRAN output 
written by the following statements:

	int*4 n_particles
	real*4 pos_x(n_particles), pos_y(n_particles), pos_z(n_particles)
	write(*) n_particles
	write(*) (pos_x(j),j=1,n_particles)
	write(*) (pos_y(j),j=1,n_particles)
	write(*) (pos_z(j),j=1,n_particles)

and that all particles have equal masses, chosen to be 1/n_particles */

int ReadSimple(KD kd,FILE *fp)
{
    int f77read(FILE *f, void *p, int maxbytes);
    int j;
    int header[100]; 	/* Assuming that the first FORTRAN block 
			is smaller than this */
    float *readlist;

    /* First, find out how many particles are involved */
    f77read(fp,header,400);  
    kd->nActive = header[0];  /* The number of particles */

    /* Allocate space to hold their positions */
    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
    assert(kd->p != NULL);

    /* Allocate temporary space to hold the data array */
    readlist = (float *)malloc(4*kd->nActive);

    /* Now read the X positions and transcribe them */
    f77read(fp,readlist,kd->nActive*4);  
    for (j=0; j<kd->nActive; j++)
	kd->p[j].r[0] = readlist[j];	/* Note zero-offset! */

    /* Repeat for Y and Z */
    f77read(fp,readlist,kd->nActive*4);  
    for (j=0; j<kd->nActive; j++)
	kd->p[j].r[1] = readlist[j];	
    f77read(fp,readlist,kd->nActive*4);  
    for (j=0; j<kd->nActive; j++)
	kd->p[j].r[2] = readlist[j];	

    /* Assume the particle mass is 1/kd->nActive */
#ifdef DIFFERENT_MASSES
    for (j=0;j<kd->nActive;j++) kd->p[j].fMass= 1.0/kd->nActive;
#else
    kd->fMass = 1.0/kd->nActive;	
#endif

    /* Give up the temp space */
    free(readlist);
    return kd->nActive;
}

/* ================================================================ */
/* =================== An Example: ReadASCII() =================== */
/* ================================================================ */

/* Read the following format -- an ASCII file with each particle's
information listed line by line:

	Line 1:		N_particles
	Line 2 to N+1:	n X Y Z Mass	

where n is the number of the particle, (X, Y, Z) is the position vector,
and Mass is the mass. */

int ReadASCII(KD kd,FILE *fp)
{
    int j, npart, dummy;
    float pos[3], mass;
    char line[200];	/* We'll read the file line-by-line */
    void f77error(char *s);  /* Report and die */

#ifndef DIFFERENT_MASSES
    /* Our format calls for individual masses, yet we have nowhere to put 
	them and no logical fallback. See the Makefile to compile with
	-DDIFFERENT_MASSES */
    fprintf(stderr,"Don't know what to do with masses.");
    exit(1);
#endif

    /* First, find out how many particles are involved */
    if (fgets(line,200,fp)==NULL) f77error("Unexpected EOF in first line.");
    if (sscanf(line,"%d",&npart)!=1) f77error("Couldn't parse first line.");
    kd->nActive = npart;  /* The number of particles */

    /* Allocate space to hold their positions */
    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
    assert(kd->p != NULL);

    for (j=0;j<npart;j++) {	/* Look at each line */
	if (fgets(line,200,fp)==NULL) f77error("Unexpected EOF.");
	if (sscanf(line,"%d %f %f %f %f", &dummy, pos, pos+1, pos+2, &mass)!=5){
		fprintf(stderr,"Couldn't parse line %d.\n",j+1);
		exit(1);
	}
	/* I won't compare dummy to anything, although this could be a check */
	kd->p[j].r[0] = pos[0];
	kd->p[j].r[1] = pos[1];
	kd->p[j].r[2] = pos[2];
#ifdef DIFFERENT_MASSES 	
	kd->p[j].fMass = mass;
#endif
    }
    return kd->nActive;
}

/* ================================================================ */
/* ====================== An Example: ReadTPM() =================== */
/* ================================================================ */

/* We need to read from Guohong Xu's TPM format */
/* To give info to the user: INFORM("info"); */
#define INFORM(string) printf(string); fflush(stdout)

int ReadTPM(KD kd,FILE *fp)
{
    int f77read(FILE *f, void *p, int maxbytes);
    int header[100];
    int i, j, bl, blocksize;
    float *readlist, masspart;

    f77read(fp,header,400);  /* All the cosmological information */
    f77read(fp,header,8);	/* The particle and block count */
    kd->nActive = header[0];  /* The number of particles */

    /* We're going to use all the particles */
    /* We won't set the following variables; they aren't used unless
    you want to *output* in tipsy format, which isn't my convention: */
    /* kd-> nDark, nGas, nStar, nParticles, fTime, bDark, bGas, bStar */

    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
    assert(kd->p != NULL);

    /* This file format has divided the particle information into 
    header[1] sets.  Each set contains px[], py[], pz[], vx[], vy[],
    and vz[] vectors for the given fraction of the particles. */

    blocksize = header[0]/header[1];
    readlist = (float *)malloc((size_t)(4*blocksize)); 
    assert(readlist!=NULL);
    /* readlist is zero offset for particles */
    
    printf("nActive = %d, blocksize = %d\n", 
	kd->nActive, blocksize);
    INFORM("Reading particles..");
    for (bl=0;bl<header[1];bl++) {
	f77read(fp,readlist,blocksize*4);  /* position_x */
	for (j=0,i=bl*blocksize; j<blocksize; j++,i++)
	    kd->p[i].r[0] = readlist[j];
	f77read(fp,readlist,blocksize*4);  /* position_y */
	for (j=0,i=bl*blocksize; j<blocksize; j++,i++)
	    kd->p[i].r[1] = readlist[j];
	f77read(fp,readlist,blocksize*4);  /* position_z */
	for (j=0,i=bl*blocksize; j<blocksize; j++,i++)
	    kd->p[i].r[2] = readlist[j];
	f77read(fp,readlist,blocksize*4);  /* velocity_x */
	f77read(fp,readlist,blocksize*4);  /* velocity_y */
	f77read(fp,readlist,blocksize*4);  /* velocity_z */
	INFORM(".");
    }
    free(readlist); 

    masspart = 1.0/kd->nActive;	/* All particles have the same mass,
			    chosen so that the average density is 1. */

#ifdef DIFFERENT_MASSES
    for (i=0;i<kd->nActive;i++) kd->p[i].fMass=masspart;
#else
    kd->fMass = masspart;	
#endif

    INFORM("Done!\n");
    return kd->nActive;
}

/* ================================================================ */
/* ===================== Some FORTRAN utilities =================== */
/* ================================================================ */

void f77error(char *s)
{
    fprintf(stderr,"%s\n",s); exit(1);
}

int f77read(FILE *f, void *p, int maxbytes)
/* Read a FORTRAN style block from the given file */
/* maxbytes is the amount of space (in bytes) the pointer p points to */
/* Space must be allocated to read the whole block into p */
/* Return amount read, scream if there's a problem */
/* Reading is done ZERO-OFFSET */
{
    int size, size2;
    if (fread(&size,4,1,f)!=1) 
        f77error("f77read(): Error reading begin delimiter.");
    if (size>maxbytes) 
        f77error("f77read(): Block delimiter exceeds size of storage.");
    if (size<maxbytes) 
        fprintf(stderr,"f77read(): Block size is smaller than size of storage.");
    if (fread(p,1,size,f)!=size) f77error("f77read(): Error reading data.");
    if (fread(&size2,4,1,f)!=1) 
        f77error("f77read(): Error reading end delimiter.");
    if (size!=size2) 
	f77error("f77read(): Delimiters do not match.");
    return size;
}

/*  For completeness.... */
int f77write(FILE *f, void *p, int len)
/* len is number of bytes to be written from p[0..len-1] */
/* Return 0 if successful, 1 if not */
{
    if (fwrite(&len,4,1,f)!=1) return 1;
    if (fwrite(p,1,len,f)!=len) return 1;
    if (fwrite(&len,4,1,f)!=1) return 1;
    return 0;
}
/* ************************************************* */
/* My Routine*/


struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;


struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P;


int ReadRamses(KD kd,char *inputfile)
{
  int dummy,npart,npart2,i,j,icpu,ncpu,ndim,ncurr,nstar_tot,nsink,nstar;
  int localseed[4];
  double mstar_tot,mstar_lost,m_tot;
  float *temp;
  double *temp_dbl;
  double *x_dbl,*y_dbl,*z_dbl,*a_dbl,*m_dbl;
  int *temp_int;
  int *temp_id;
  float pmax;
  FILE *fp,*fascii;
  char *tempi;
  char currfile[80];

  //**********************************************
  // First Read: Finding the number of Ramses files
  //**********************************************

  icpu=1;
  tempi=inputfile;
  strncat(tempi,"%05d",19);

  sprintf(currfile,tempi,icpu);

  fp=fopen(currfile,"r");
  printf("opening %s\n",currfile);
  
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&ncpu, sizeof(ncpu), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&ndim, sizeof(ndim), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&npart, sizeof(npart), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&localseed[0], sizeof(int), 4, fp);
  fread(&dummy, sizeof(dummy), 1, fp);      
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&nstar, sizeof(nstar), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);

  printf("Number of Ramses files = %d \n",ncpu);
  /*  printf("ndim = %d \n",ndim);
  printf("npart= %d \n",npart);
  printf("localseed: %d \n",localseed[0]);
  printf("localseed: %d \n",localseed[1]);
  printf("localseed: %d \n",localseed[2]);
  printf("localseed: %d \n",localseed[3]);
  printf("nstar= %d \n",nstar); */
  if(nstar>0)printf("Found stars in Ramses files. \n");
  fclose(fp);


  //**********************************************
  // Second Read: Finding the number particles
  //**********************************************

  npart=0;
  
  for(icpu=1;icpu<=ncpu;icpu++)
    {
      sprintf(currfile,tempi,icpu);
      
      fp=fopen(currfile,"r");
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&ncpu, sizeof(ncpu), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&ndim, sizeof(ndim), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&npart2, sizeof(npart2), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fclose(fp);
      npart=npart+npart2;
    }

  //**********************************************
  // if star are present, we need to remove them
  //**********************************************
  if(nstar>0){

    printf("Maximum number of particles = %d \n",npart);
    printf("Removing stars from particle list. \n");

    ncurr=0;
    for(icpu=1;icpu<=ncpu;icpu++)
      {
	sprintf(currfile,tempi,icpu);
	
	fp=fopen(currfile,"r");
	//	printf("opening %s\n",currfile);
	
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&ncpu, sizeof(ncpu), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&ndim, sizeof(ndim), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&npart2, sizeof(npart2), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&localseed[0], sizeof(int), 4, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&nstar_tot, sizeof(nstar_tot), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&mstar_tot, sizeof(mstar_tot), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&mstar_lost, sizeof(mstar_lost), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&nsink, sizeof(nsink), 1, fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	temp_int=(int *) malloc(npart2*sizeof(int));
	temp_id=(int *) malloc(npart2*sizeof(int));
	temp_dbl=(double *)malloc(npart2*sizeof(double));

	// Skipping positions
	for(i=0;i<=ndim-1;i++)
	  {
	    fread(&dummy, sizeof(dummy), 1, fp);
	    fread(&temp_dbl[0],sizeof(double),npart2,fp);
	    fread(&dummy, sizeof(dummy), 1, fp);
	  }
	
	// Skipping velocities
	for(i=0;i<=ndim-1;i++)
	  {
	    fread(&dummy, sizeof(dummy), 1, fp);
	    fread(&temp_dbl[0],sizeof(double),npart2,fp);
	    fread(&dummy, sizeof(dummy), 1, fp);
	  }
	
	//Skipping masses
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&temp_dbl[0],sizeof(double),npart2,fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	//Skipping identity
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&temp_id[0],sizeof(int),npart2,fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	//Skipping level
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&temp_int[0],sizeof(int),npart2,fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	//Reading age
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&temp_dbl[0],sizeof(double),npart2,fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	for(j=0;j<=npart2-1;j++)
	  {
	    if(temp_dbl[j]==0. && temp_id[j]>0)ncurr=ncurr+1;
	  }
	
	free(temp_dbl);
	free(temp_int);
	free(temp_id);

	fclose(fp);
      }
    
    npart=ncurr;

  }

  printf("Actual number of DM particles = %d \n",npart);

  kd->nActive = npart;
  kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
  printf("Memory Allocated\n");
  
  printf("Reading Positions\n");

  ncurr=0;
  for(icpu=1;icpu<=ncpu;icpu++)
    {
      sprintf(currfile,tempi,icpu);
      
      fp=fopen(currfile,"r");
      //      printf("opening %s\n",currfile);
      
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&ncpu, sizeof(ncpu), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&ndim, sizeof(ndim), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&npart2, sizeof(npart2), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&localseed[0], sizeof(int), 4, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&nstar_tot, sizeof(nstar_tot), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&mstar_tot, sizeof(mstar_tot), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&mstar_lost, sizeof(mstar_lost), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&nsink, sizeof(nsink), 1, fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      temp_dbl=(double *)malloc(npart2*sizeof(double));
      x_dbl=(double *)malloc(npart2*sizeof(double));
      y_dbl=(double *)malloc(npart2*sizeof(double));
      z_dbl=(double *)malloc(npart2*sizeof(double));
      m_dbl=(double *)malloc(npart2*sizeof(double));

      if(nstar>0){
	temp_int=(int *)malloc(npart2*sizeof(int));
	temp_id=(int *)malloc(npart2*sizeof(int));
	a_dbl=(double *)malloc(npart2*sizeof(double));
      }

      // Read x position
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&x_dbl[0],sizeof(double),npart2,fp);
      fread(&dummy, sizeof(dummy), 1, fp);

      // Read y position
      fread(&dummy, sizeof (dummy), 1, fp);
      fread(&y_dbl[0],sizeof(double),npart2,fp);
      fread(&dummy, sizeof(dummy), 1, fp);

      // Read z position
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&z_dbl[0],sizeof(double),npart2,fp);
      fread(&dummy, sizeof(dummy), 1, fp);

      // Skipping velocities
      for(i=0;i<=ndim-1;i++)
	{
	  fread(&dummy, sizeof(dummy), 1, fp);
	  fread(&temp_dbl[0],sizeof(double),npart2,fp);
	  fread(&dummy, sizeof(dummy), 1, fp);
	}
      
      //Reading Mass
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&m_dbl[0],sizeof(double),npart2,fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      //If stars are present, remove them
      if(nstar>0){
	//Skipping identity
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&temp_id[0],sizeof(int),npart2,fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	//Skipping level
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&temp_int[0],sizeof(int),npart2,fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	//Reading age
	fread(&dummy, sizeof(dummy), 1, fp);
	fread(&a_dbl[0],sizeof(double),npart2,fp);
	fread(&dummy, sizeof(dummy), 1, fp);
	
	for(j=0;j<=npart2-1;j++)
	  {
	    if(a_dbl[j]==0.&&temp_id[j]>0){
	      kd->p[ncurr].r[0]=x_dbl[j];
	      kd->p[ncurr].r[1]=y_dbl[j];
	      kd->p[ncurr].r[2]=z_dbl[j];
#ifdef DIFFERENT_MASSES
	      kd->p[ncurr].fMass=m_dbl[j];
#else
	      kd->fMass = m_dbl[j];
#endif
	      ncurr=ncurr+1;
	    }
	  }
      }
      
      //If no star are present, just copy particles
      if(nstar==0){
	for(j=0;j<=npart2-1;j++)
	  {
	    kd->p[ncurr].r[0]=x_dbl[j];
	    kd->p[ncurr].r[1]=y_dbl[j];
	    kd->p[ncurr].r[2]=z_dbl[j];
#ifdef DIFFERENT_MASSES
	    kd->p[ncurr].fMass=m_dbl[j];
#else
	    kd->fMass = m_dbl[j];
#endif
	    ncurr=ncurr+1;
	  }

      }

      free(temp_dbl);
      free(x_dbl);
      free(y_dbl);
      free(z_dbl);
      free(m_dbl);
      if(nstar>0){
	free(a_dbl);
	free(temp_int);
	free(temp_id);
      }

      fclose(fp);
      
    }

#ifdef DIFFERENT_MASSES
  m_tot=0.;
  
  for(j=0;j<=npart-1;j++){
    m_tot=m_tot+kd->p[j].fMass;
    
  }
  printf("Total Mass  =%f\n", m_tot);
  printf("Renormalizing particle masses to get Mtot=1\n");
  
  for(j=0;j<=npart-1;j++){
    kd->p[j].fMass=kd->p[j].fMass/m_tot;   
  }
  
#else
  printf("Total Mass =%f\n", kd->fMass*(npart));
  printf("Renormalizing particle masses to get Mtot=1\n");
  kd->fMass=1./float(npart);
#endif  

  return kd->nActive;	

}
/************************************************** */


int ReadASCII2(KD kd,FILE *fp)
{
  int j,npart,dummy;
  float pos[3], mass;

  fscanf(fp,"%d\n",&npart);  
  fscanf(fp,"%f\n",&mass);

  kd->nActive = npart;
  kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));

  for(j=0;j<npart;j++)
    {
      fscanf(fp,"%d %f %f %f \n",&dummy,&pos[0],&pos[1],&pos[2]);
      kd->p[j].r[0] = pos[0];
      kd->p[j].r[1] = pos[1];
      kd->p[j].r[2] = pos[2];
    }
  kd->fMass=1.e10*mass;
 /* kd->fMass=1000000.*mass;*/
  printf("%f\n",kd->fMass);
  return kd->nActive;
}

