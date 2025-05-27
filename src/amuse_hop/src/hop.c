/* HOP.C -- Daniel Eisenstein, 1997
Based on a paper by Daniel Eisenstein & Piet Hut, 
"HOP: A New Group-Finding Algorithm for N-body Simulations."
See the included documentation or view it at 
http://www.sns.ias.edu/~eisenste/hop/hop_doc.html */

/* main() was customized from that of SMOOTH v2.0.1, by Joachim Stadel
   and the NASA HPCC ESS at the University of Washington Dept. of Astronomy. */
/* PrepareKD() is a code fragment from the same program. */
/* ssort() is a C-translation of the Slatec FORTRAN routine of
   R.E. Jones and J.A. Wisniewski (SNLA) */
/* The other routines were written by DJE. */

/* Version 1.0 (12/15/97) -- Original Release */
/* Version 1.1 (04/02/01) -- No changes to this file.
			     Tiny bug fix in regroup.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "kd.h"
#include "smooth.h"
/* To give info to the user: INFORM("info"); */
#define INFORM(string) printf(string); fflush(stdout)

#ifndef HUGE
#define HUGE 3.40282e+38
#endif

int ReadSimulationFile(KD, FILE *);

void smDensityTH(SMX smx,int pi,int nSmooth,int *pList,float *fList);

void smHop(SMX smx,int pi,int nSmooth,int *pList,float *fList);
void FindGroups(SMX smx);
void SortGroups(SMX smx);

void MergeGroupsHash(SMX smx);
void smMergeHash(SMX smx,int pi,int nSmooth,int *pList,float *fList);
void ReSizeSMX(SMX smx, int nSmooth);

void PrepareKD(KD kd);
void binOutHop(SMX smx, FILE *fp);
void binOutDensity(SMX smx, FILE *fp);
void binInDensity(SMX smx, FILE *fp);
void outGroupMerge(SMX smx, FILE *fp);

void usage(void)
{
    fprintf(stderr,"USAGE:\n");
    fprintf(stderr,"hop [-den <file>] Name of file from which to input densities.\n");
    fprintf(stderr,"    If omitted, calculate the densities: \n");
    fprintf(stderr,"    [-nd <int>] N_dens -- # of particles to smooth for density, default 64.\n");
    fprintf(stderr,"    [-tophat or -th] Use tophat kernal for density, else cubic spline.\n");
    fprintf(stderr,"    [-g] Use gather-only for spline kernal, else gather-scatter.\n\n");

    fprintf(stderr,"    [-nh <int>] N_hop -- # of particles over which to look for density max, default 16.\n");
    fprintf(stderr,"    [-dt <float>] Density below which we don't assign to any group, default none.\n");
    fprintf(stderr,"    [-nm <int>] N_merge -- # of particles to catalog group boundaries, default 4.\n\n");

    fprintf(stderr,"    [-in <file>] Simulation particle data file.\n");
    fprintf(stderr,"    	     If omitted, read from stdin.\n");
    fprintf(stderr,"    [-p <float>] Periodicity of the simulation box, default infinite.\n\n");

    fprintf(stderr,"    [-o/-out <file>] Root of output file names.\n");
    fprintf(stderr,"    [-nodensity] [-nohop] [-nogbound] Turn off output of these files.\n");
    fprintf(stderr,"    [-densityonly] Compute the densities and quit.\n");
    fprintf(stderr,"    [-b <int>] Performance tuning in kd-tree search.\n");
    exit(1);
}


void oldmain(int argc,char **argv)
{
	KD kd;
	SMX smx;
	int nBucket,nSmooth,i,j;
	FILE *fp, *fpb;
	char ach[80],achFile[80], *inputfile, *densfile;
	float fPeriod[3];
	int bDensity,bGroup,bSym,bMerge,nDens,nHop,nMerge,bTopHat;
	float fDensThresh;
	
  nBucket = 16;
	nSmooth = 64;
	nDens = 64;
	nHop = -1;
	fDensThresh = -1.0;
	bDensity = 3;
	bGroup = 3;
	bMerge = 3;
	bSym = 1;
	bTopHat = 0;
	strcpy(achFile,"output_hop");
	inputfile = NULL;
	i = 1;
	for (j=0;j<3;++j) fPeriod[j] = HUGE;
	nMerge = 4;

	if (argc==1) usage();	/* No arguments is a request for help */
	while (i < argc) {
		if (!strcmp(argv[i],"-b")) {
			++i;
			if (i >= argc) usage();
			nBucket = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-nd")) {
			++i;
			if (i >= argc) usage();
			nDens = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-nh")) {
			++i; if (i>=argc) usage();
			nHop = atoi(argv[i]); ++i;
			}
		else if (!strcmp(argv[i],"-tophat")||!strcmp(argv[i],"-th")) {
			++i; bTopHat = 1;
			}
		else if (!strcmp(argv[i],"-dt")) {
			++i; if (i>=argc) usage();
			sscanf(argv[i],"%f", &fDensThresh); ++i;
			}
		else if (!strcmp(argv[i],"-nm")) {
			++i; if (i>=argc) usage();
			nMerge = atoi(argv[i]); ++i;
			}
		else if (!strcmp(argv[i],"-in")) {
			++i;
			if (i >= argc) usage();
			inputfile = argv[i];
			++i;
			}
		else if (!strcmp(argv[i],"-den")) {
			++i; if (i >= argc) usage();
			densfile = argv[i];
			bDensity = 0; ++i;
			}
		else if (!strcmp(argv[i],"-o") || !strcmp(argv[i],"-out")) {
			++i;
			if (i >= argc) usage();
			strcpy(achFile,argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-g")) {
			bSym = 0;
			++i;
			}
		else if (!strcmp(argv[i],"-p")) {
			++i;
			if (i >= argc) usage();
			/* Setting all directions to have the same periodicity.
			You could change this if you want; for example, see
			the SMOOTH program from the HPCC. */
			fPeriod[0] = atof(argv[i]); 
			fPeriod[1] = atof(argv[i]);
			fPeriod[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-nodensity")) {
			bDensity &= 1;	/* The 2-bit must be off */
			++i;
			}
		else if (!strcmp(argv[i],"-nohop")) {
			bGroup &= 1;	/* The 2-bit must be off */
			++i;
			}
		else if (!strcmp(argv[i],"-densityonly")) {
			bGroup = 0; ++i;
			bMerge = 0;  /* Merging without Grouping isn't coded */
			}
		else if (!strcmp(argv[i],"-nogbound")) {
			bMerge = 0; ++i; /* If we're not writing it, we needn't
					waste time computing it */
			}
		else usage();
		}
	if (nHop<0) nHop=nDens;
	if (bDensity==0) nSmooth = nHop+1;
	    else nSmooth = nDens+1;
	/* When smSmooth() is asked for nSmooth particles, it seems to
	generally return nSmooth-1 particles, including primary itself.
	Hence, when we want nDens or nSmooth particles (including the
	primary) we should ask for one more.  By convention, I've chosen
	to have nMerge reflect the number *not* including the primary,
	so in this case we need to ask for two more! */

	assert(!bMerge || nMerge<nHop-1);
	    /* This is generally the case of interest, and there's an
	    optimization to be had.  I haven't the coded the other case */
	assert(!bMerge || bGroup);
	    /* Can't run Merging without Grouping */
	if (!bMerge) nMerge = 9999999;	
		/* If we're not Merging, we don't want smHop() to set 
		up for Merging.  This kludge avoids that */

	kdInit(&kd,nBucket);

	if (inputfile==NULL) fp = stdin;
	    else if ((fp=fopen(inputfile,"r"))==NULL) {
		fprintf(stderr,"Input file not found.\n"); exit(1);
	    }
	/* ReadSimulationFile(kd,fp); */ // LHarms 30-10-2010
	PrepareKD(kd);
	if (fp!=NULL) fclose(fp);

	smInit(&smx,kd,nSmooth,fPeriod);
	smx->nHop = nHop;
	smx->nDens = nDens;
	smx->nMerge = nMerge;
	smx->nGroups = 0;
	smx->fDensThresh = fDensThresh;

	if (bDensity==0) {
	    INFORM("Reading Densities...\n");
	    fp=fopen(densfile,"r");
	    assert(fp!=NULL);
	    binInDensity(smx,fp);
	}
	INFORM("Building Tree...\n");
	kdBuildTree(kd);

	if (bDensity) {
	    INFORM("Finding Densities...\n");
	    if (bTopHat) smSmooth(smx,smDensityTH);
	    else if (bSym) smSmooth(smx,smDensitySym);
	    else smSmooth(smx,smDensity);
	}  /* Else, we've read them */
	if (bGroup) {
	     INFORM("Finding Densest Neighbors...\n");
	     if (bDensity && nHop<nSmooth) smReSmooth(smx,smHop);
	     else {
		if (nHop>=nSmooth) {
		    nSmooth = nHop+1;
		    ReSizeSMX(smx,nSmooth);
		}
		smSmooth(smx,smHop);
	    }
	}

	INFORM("Grouping...\n");
	if (bGroup) FindGroups(smx);
	if (bGroup) SortGroups(smx);

	if (bMerge) {
	    INFORM("Merging Groups...\n");
	    MergeGroupsHash(smx);
	}

	kdOrder(kd);
	INFORM("Writing Output...\n");
	if (bDensity&2) {
	    strcpy(ach,achFile); strcat(ach,".den");
	    fp = fopen(ach,"w"); assert(fp != NULL);
	    binOutDensity(smx,fp);
	    fclose(fp);
	}

	if (bMerge&2) {
	    strcpy(ach,achFile); strcat(ach,".gbound");
	    fp = fopen(ach,"w"); assert(fp != NULL);
	    smx->nSmooth=nSmooth; /* Restore this for output */
	    outGroupMerge(smx,fp);  
	    fclose(fp);
	}
	if (bMerge) free(smx->hash);

	if (bGroup&2) {
	    strcpy(ach,achFile); strcat(ach,".hop");
	    fp = fopen(ach,"w"); assert(fp != NULL);
	    binOutHop(smx,fp);
	    fclose(fp);
	}
	if (bGroup) {free(smx->densestingroup); free(smx->nmembers);}
	smFinish(smx);
	kdFinish(kd);
	INFORM("All Done!");
	return;
}


/* ============================================================= */
/* ===================== New Density Routine =================== */
/* ============================================================= */

void smDensityTH(SMX smx,int pi,int nSmooth,int *pList,float *fList)
/* Find density only using top-hat kernal. */
{
#ifdef DIFFERENT_MASSES
    int j;
    float totalmass;
    for (j=0,totalmass=0.0; j<nSmooth; j++)
	totalmass += smx->kd->p[pList[j]].fMass;
    smx->kd->p[pi].fDensity = totalmass*0.75*M_1_PI/
	smx->pfBall2[pi]/sqrt(smx->pfBall2[pi]);
#else
    /* This case is simple: the total mass is nSmooth times the mass
	per particle */
    smx->kd->p[pi].fDensity = (nSmooth)*smx->kd->fMass*0.75*M_1_PI/
	smx->pfBall2[pi]/sqrt(smx->pfBall2[pi]);
#endif
    return;
}

/* ============================================================= */
/* ================== Hop to Neighbors/Form Groups ============= */
/* ============================================================= */

void smHop(SMX smx,int pi,int nSmooth,int *pList,float *fList)
/* Look at the nHop nearest particles and find the one with the
highest density.  Store its ID number in iHop as -1-ID (to make it negative) */
/* nSmooth tends to be the expected value (smx->nSmooth-1) but can vary plus
or minus 1 (at least). */
/* If Merging is turned off, nMerge should be huge so as to avoid extra
sorting below */
{
    int i,max, search, didsort;
    float maxden;
    void ssort(float X[], int Y[], int N, int KFLAG);

    /* If the density is less than the threshold requirement, then assign 0 */
    if (smx->kd->p[pi].fDensity<smx->fDensThresh) {
	smx->kd->p[pi].iHop = 0;
	return;
    }

    /* If we have exactly the right number, then it doesn't matter how
    we search them.  Otherwise, we need to sort first. */
    /* We can destroy pList and fList if we want. fList holds the square radii*/

    search = smx->nHop;
    if (smx->nHop>nSmooth) search = nSmooth;	/* Very rare exception */

    if (smx->nHop<nSmooth || smx->nMerge+2<nSmooth) {
	ssort(fList-1, pList-1, nSmooth, 2);
	didsort = 1;
    } else didsort = 0;

    maxden = 0.0;
    for (i=0;i<search;++i) {
	if (smx->kd->p[pList[i]].fDensity>maxden) {
	    max = i;
	    maxden = smx->kd->p[pList[i]].fDensity;
	}
    }
    
    smx->kd->p[pi].iHop = -1-pList[max];

    /* If a sort was done, then we can save time in the Merge step by
    recording the new Ball radius. */
    /* Place the new radius in between the two boundary because occasionally
    the floating point math comes out strange when comparing two floats */
    if (didsort && smx->nMerge+2<nSmooth) 
	smx->pfBall2[pi] = 0.5*(fList[smx->nMerge+1]+fList[smx->nMerge]);
    return;
}

/* ----------------------------------------------------------------- */

void FindGroups(SMX smx)
/* Number the maxima. Trace each particle uphill to a maximum. */
/* The local maxima were stored as in iHop as -1-ID (negative numbers);
now we will store group number as positive numbers (1..n) in the same spot */
/* Zero is left as an error condition */
/* The particles MUST be left in the BuildTree order, for that is how the
iHop tracing is done */
/* Allocate space for densestingroup, from 0 to nGroups 
(inclusive) and store the particle number of the maxima, which is the 
densest particle in the group.  Ditto for nmembers[], the number of
particles in the group */
{
    int j, ng, current, localmax;
    PARTICLE *p;

    smx->nGroups = 0;
    /* First look for maxima, where particle ID = iHop.  Number the groups */
    for (j=0, p=smx->kd->p;j<smx->kd->nActive;j++,p++) 
	if (p->iHop == -1-j) {  /* Was p->iOrder */
	    /* Yes, it's a maximum */
	    smx->nGroups++;
	    /* p->iHop = smx->nGroups; */
	}

    /* Now catalog the maxima, before numbering the groups */
    smx->densestingroup = (int *)malloc((size_t)((smx->nGroups+1)*4));
    assert(smx->densestingroup!=NULL);
    smx->nmembers = (int *)malloc((size_t)(4*(smx->nGroups+1)));
    assert(smx->nmembers!=NULL);

    ng = 0;
    for (j=0,p=smx->kd->p;j<smx->kd->nActive;j++,p++)
	if (p->iHop== -1-j) {  
	    /* It's a maximum */
	    ng++;
	    smx->densestingroup[ng] = p->iOrder;
	    p->iHop = ng;
	}
    
    /* Now take the remaining particles and trace up to a maximum */
    for (j=0,p=smx->kd->p;j<smx->kd->nActive;j++,p++) {
	if (p->iHop>=0) continue;	/* It's a maximum or an error */
	localmax = -1-p->iHop;
	while (smx->kd->p[localmax].iHop<0)
	    localmax = -1-smx->kd->p[localmax].iHop;
	ng = smx->kd->p[localmax].iHop;
	/* Now assign this group number to the whole lineage! */
	/* Note that errors (iHop=0) will propagate */
	current = -1-p->iHop;
	p->iHop = ng;
	while (smx->kd->p[current].iHop<0) {
	    localmax = -1-smx->kd->p[current].iHop;
	    smx->kd->p[current].iHop = ng;
	    current = localmax;
	}
    }
    return;
}

/* ----------------------------------------------------------------- */

void SortGroups(SMX smx)
/* Renumber the groups in order of number of members. */
/* Move the group numbering from unit offset to zero offset */
/* Move error condition (group=0) to (group=-1) */
/* Store the number of members and most dense particle in each group */
{
    int j, *indx, *irank, *ip;
    PARTICLE *p;
    void make_rank_table(int n, int *ivect, int *rank);

    indx = (int *)malloc((size_t)(4*(smx->nGroups+1)));
    assert(indx!=NULL);
    irank = (int *)malloc((size_t)(4*(smx->nGroups+1)));
    assert(irank!=NULL);

    /* Count the number of members in each group */
    for (j=0;j<=smx->nGroups;j++) smx->nmembers[j]=0;
    for (j=0,p=smx->kd->p;j<smx->kd->nActive;j++,p++)
	smx->nmembers[p->iHop]++;

    make_rank_table(smx->nGroups, smx->nmembers,irank);

    for (j=1;j<=smx->nGroups;j++) irank[j] = smx->nGroups-irank[j];
    irank[0] = -1;	/* Move old 0's to new -1's */
    /* irank[j] is the new group number of group j: zero-offset, ordered
    large to small */

    /* Relabel all the particles */
    for (j=0,p=smx->kd->p;j<smx->kd->nActive;j++,p++) 
	p->iHop = irank[p->iHop];	

    /* Sort the nmembers and densestingroup lists to reflect the new ordering */
    /* Use indx as a temp space */
    for (j=1;j<=smx->nGroups;j++) indx[irank[j]]=smx->densestingroup[j];
    ip = smx->densestingroup;
    smx->densestingroup = indx;
    indx = ip;
    /* Number of error (old_group=0) is in nmembers[0]; move it to [nGroups] */
    for (j=1;j<=smx->nGroups;j++) indx[irank[j]]=smx->nmembers[j];
    indx[smx->nGroups]=smx->nmembers[0];
    free(smx->nmembers);
    smx->nmembers = indx;

    free(irank);
    /* Note that the memory allocated to indx is now used by smx->densestingroup
    and so it should not be free'd.  */
    return;
}

/* ================================================================== */
/* ========================== Group Merging ========================= */
/* ================================================================== */

void MergeGroupsHash(SMX smx)
/* We're going to look through the particles looking for boundary particles,
defined as particles with close neighbors of a different group, and store
the most dense boundary point (average of the two points) */
/* The matrix of boundaries is stored in a hash table */
/* SortGroups() should be called previous to this, so that all the
particles are in the assumed group numbering, i.e. 0 to ngroup-1, with
-1 being unattached. The tags are not altered */ 
/* In smHop, if nMerge+2 was smaller than nSmooth, we set the new radius
for searching.  If not, we left the old radius alone.  Either way, we're
ready to go. */
{
    int j, k, g, next, newgroup;
    
    ReSizeSMX(smx, smx->nMerge+2);	/* Alter the smoothing scale on smx */
    smx->nHashLength = smx->nGroups*10+1;
    smx->hash = (Boundary *)malloc(smx->nHashLength*sizeof(Boundary));
    assert(smx->hash!=NULL);
    for (j=0;j<smx->nHashLength;j++) {
	smx->hash[j].nGroup1 = -1;
	smx->hash[j].nGroup2 = -1;
	smx->hash[j].fDensity = -1.0;
    } /* Mark the slot as unused */

    smReSmooth(smx,smMergeHash);	/* Record all the boundary particles */
    return;
}

/* ----------------------------------------------------------------- */

void smMergeHash(SMX smx,int pi,int nSmooth,int *pList,float *fList)
/* Look at the list for groups which are not that of the particle */
/* If found, and if density is high enough, then mark it as a boundary */
/* by recording it in the hash table */
{
    int j,group;
    float averdensity;
    int g1,g2,count;
    unsigned long hashpoint;
    Boundary *hp;
    int search;
    void ssort(float X[], int Y[], int N, int KFLAG);

    group = smx->kd->p[pi].iHop;
    if (group==(-1)) return;	/* This particle isn't in a group */

    /* It seems that BallGather doesn't always get the right number of
    particles....*/
    search = nSmooth;
    if (nSmooth>smx->nMerge+1) {
	ssort(fList-1,pList-1,nSmooth,2);
	search = smx->nMerge+1;
    }
    for (j=0;j<search;j++) {
	g2=smx->kd->p[pList[j]].iHop;
	if (g2==-1 || g2==group) continue;  /* Same group or unassigned */
	/* It's in a different group; we need to connect the two */
	if (group<g2) g1=group;
	    else {g1=g2; g2=group;}
	averdensity = 0.5*(smx->kd->p[pi].fDensity + 
			smx->kd->p[pList[j]].fDensity);
	hashpoint = (g1+1)*g2;  /* Avoid multiplying by 0 */
	hashpoint = hashpoint % smx->nHashLength;
	hp = smx->hash+hashpoint;
	count = 0;
	for (;;) {
	    if (hp->nGroup1==(-1)) { 	/* Empty slot */
		hp->nGroup1 = g1;
		hp->nGroup2 = g2;
		hp->fDensity = averdensity;
		break;	
	    }
	    if (hp->nGroup1==g1 && hp->nGroup2==g2) {	
		/* We've seen this pair of groups before */
		if (hp->fDensity > averdensity) break;
		else {
		    hp->fDensity = averdensity;
		    break;
		}
	    }
	    /* Else, this slot was full, go to the next one */
	    hp++; 
	    if (hp>=smx->hash+smx->nHashLength) hp = smx->hash;
	    if (++count>1000) {
		fprintf(stderr,"Hash Table is too full.\n");
		exit(1);
	    }
	}
	/* Look at the next particle */
    }
    return;
}


/* ----------------------------------------------------------------- */

void ReSizeSMX(SMX smx, int nSmooth)
/* Set a new smoothing length, resizing the arrays which depend on this,
but leaving the particle information intact. */
/* However, because we won't always have resized pfBall2 (the search 
radius) correctly, we won't reduce the size of the fList and pList
arrays */
{
    PQ_STATIC;
    if (nSmooth>smx->nSmooth) {	/* We're increasing the size */
	smx->nListSize = nSmooth+RESMOOTH_SAFE;
	free(smx->fList);
	smx->fList = (float *)malloc(smx->nListSize*sizeof(float));
	assert(smx->fList != NULL);
	free(smx->pList);
	smx->pList = (int *)malloc(smx->nListSize*sizeof(int));
	assert(smx->pList != NULL);
    }
    smx->nSmooth=nSmooth;
    free(smx->pq);
    smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
    assert(smx->pq != NULL);
    PQ_INIT(smx->pq,nSmooth);
    return;
}

/* ===================================================================== */
/* ===================== Input/Output, Binary and ASCII ================ */
/* ===================================================================== */

void PrepareKD(KD kd)
/* This labels all the particles and finds the min/max of the positions */
/* It used to appear in kd.c within kdReadTipsy(), but it seems so general
that I'll spare the user the trouble of including it in any custom input
routines */
{
    BND bnd;
    int i, j;

    /* Label the particles, so that we can restore the order at the end */
    for (i=0;i<kd->nActive;i++) kd->p[i].iOrder=i;

    /*
     ** Calculate Bounds.
     */
    for (j=0;j<3;++j) {
        bnd.fMin[j] = kd->p[0].r[j];
        bnd.fMax[j] = kd->p[0].r[j];
        }
    for (i=1;i<kd->nActive;++i) {
        for (j=0;j<3;++j) {
            if (bnd.fMin[j] > kd->p[i].r[j])
                bnd.fMin[j] = kd->p[i].r[j];
            else if (bnd.fMax[j] < kd->p[i].r[j])
                bnd.fMax[j] = kd->p[i].r[j];
            }
        }
    kd->bnd = bnd;
    return;
}

void binOutHop(SMX smx, FILE *fp)
/* Write Group tag for each particle.  Particles should be ordered. */
/* Binary file: nActive, nGroups, list of Groups */
{
    int j,dummy;

    dummy=smx->kd->nActive; fwrite(&dummy, 4, 1, fp);
    dummy=smx->nGroups; fwrite(&dummy, 4, 1, fp);
    for (j=0;j<smx->kd->nActive;j++) fwrite(&(smx->kd->p[j].iHop),4,1,fp);
    return;
}

/* ----------------------------------------------------------------- */

void binOutDensity(SMX smx, FILE *fp)
/* Write the density for each particle.  Particles should be ordered. */
/* Binary file: nActive, list of Densities */
{
    int j,dummy;

    dummy=smx->kd->nActive; fwrite(&dummy, 4, 1, fp);
    for (j=0;j<smx->kd->nActive;j++) fwrite(&(smx->kd->p[j].fDensity),4,1,fp);
    return;
}

/* ----------------------------------------------------------------- */

void binInDensity(SMX smx, FILE *fp)
/* Read in the densities for each particle.  Particles should be ordered. */
/* Binary file: nActive, list of Densities */
{
    int j, dummy;

    if (fread(&dummy,4,1,fp)!=1) {
	fprintf(stderr,"Format of density file seems wrong.\n"); exit(1);
    }
    assert(dummy==smx->kd->nActive);
    for (j=0;j<smx->kd->nActive;j++) 
	if (fread(&(smx->kd->p[j].fDensity),4,1,fp)!=1) {
	    fprintf(stderr,"Error reading density file.\n");
	    exit(1);
	}
    return;
}

/* ----------------------------------------------------------------- */

void outGroupMerge(SMX smx, FILE *fp)
/* Write an ASCII file with information on the groups and group merging */
/* Start the boundary list with the only ### line */
/* Groups should be ordered before calling this (else densities will be wrong)*/
{
    int j, den;
    Boundary *hp;

    fprintf(fp,"%d\n", smx->nGroups);
    fprintf(fp,"# Number of Particles: %d\n", smx->kd->nActive);
    fprintf(fp,"# Number of Groups: %d\n", smx->nGroups);
    fprintf(fp,"# Number of Particles not in a Group: %d\n",
	    smx->nmembers[smx->nGroups]);
    fprintf(fp,"# nSmooth = %d, nHop = %d, nMerge = %d, fDensThresh = %6.3f.\n",
	    smx->nSmooth, smx->nHop, smx->nMerge, smx->fDensThresh);
    fprintf(fp,"# 1) Group ID.\n");
    fprintf(fp,"# 2) Number of Particles in Group.\n");
    fprintf(fp,"# 3) Densest Particle ID (zero-offset).\n");
    fprintf(fp,"# 4-6) Position of that Particle.\n");
    fprintf(fp,"# 7) Density at that Particle.\n#\n");
    for (j=0;j<smx->nGroups;j++) {
	den = smx->densestingroup[j];
	fprintf(fp,"%4d %5d %7d %6.4f %6.4f %6.4f %6.2f\n",
	    j, smx->nmembers[j], den, smx->kd->p[den].r[0], 
	    smx->kd->p[den].r[1], smx->kd->p[den].r[2], 
	    smx->kd->p[den].fDensity);
    }
    fprintf(fp,"### Begin list of boundaries. Group, Group, Average Density.\n");
    for (j=0, hp=smx->hash;j<smx->nHashLength; j++,hp++) 
	if (hp->nGroup1>=0)
	    fprintf(fp,"%5d %5d %6.2f\n",
		hp->nGroup1, hp->nGroup2, hp->fDensity);
    return;
}


/* ================================================================== */
/* ======================= Sorting ================================== */
/* ================================================================== */

typedef struct index_struct {
    float value;
    int index;
} *ptrindex;

int cmp_index(const void *a, const void *b) 
{
    if ( ((ptrindex)a)->value<((ptrindex)b)->value) return -1;
    else if ( ((ptrindex)a)->value>((ptrindex)b)->value) return 1;
    else return 0;
}

void make_rank_table(int n, int *ivect, int *rank)
/* Given a vector of integers ivect[1..n], construct a rank table rank[1..n]
so that rank[j] contains the ordering of element j, with rank[j]=n indicating
that the jth element was the highest, and rank[j]=1 indicating that it 
was the lowest.  Storage for rank[] should be declared externally */
/* I don't think this routine is particularly fast, but it's a 
miniscule fraction of the runtime */
{
    int j;
    ptrindex sortvect;

    sortvect = (ptrindex)malloc(n*sizeof(struct index_struct));
    for (j=0;j<n;j++) sortvect[j].value = (float)ivect[j+1];
    for (j=0;j<n;j++) sortvect[j].index = j+1;  /* Label them prior to sort */
    qsort(sortvect,n,sizeof(struct index_struct),cmp_index);
    /* Now sortvect is in order (smallest to largest) */
    for (j=0;j<n;j++) rank[sortvect[j].index]=j+1;
    free(sortvect);
    return;
}

/* DJE -- This is a C-translation of the Slatec FORTRAN routine ssort(). */
/* I have kept the variable names and program flow unchanged; hence
all the goto's and capitalized names.... */
/* The input vectors must be UNIT-OFFSET */

/* Here is the original Slatec header:

***BEGIN PROLOGUE  SSORT
***PURPOSE  Sort an array and optionally make the same interchanges in
            an auxiliary array.  The array may be sorted in increasing
            or decreasing order.  A slightly modified QUICKSORT
            algorithm is used.
***LIBRARY   SLATEC
***CATEGORY  N6A2B
***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
***AUTHOR  Jones, R. E., (SNLA)
           Wisniewski, J. A., (SNLA)
***DESCRIPTION

   SSORT sorts array X and optionally makes the same interchanges in
   array Y.  The array X may be sorted in increasing order or
   decreasing order.  A slightly modified quicksort algorithm is used.

   Description of Parameters
      X - array of values to be sorted   (usually abscissas)
      Y - array to be (optionally) carried along
      N - number of values in array X to be sorted
      KFLAG - control parameter
            =  2  means sort X in increasing order and carry Y along.
            =  1  means sort X in increasing order (ignoring Y)
            = -1  means sort X in decreasing order (ignoring Y)
            = -2  means sort X in decreasing order and carry Y along.

***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
                 for sorting with minimal storage, Communications of
                 the ACM, 12, 3 (1969), pp. 185-187.
***ROUTINES CALLED  XERMSG
***REVISION HISTORY  (YYMMDD)
   761101  DATE WRITTEN
   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
   890531  Changed all specific intrinsics to generic.  (WRB)
   890831  Modified array declarations.  (WRB)
   891009  Removed unreferenced statement labels.  (WRB)
   891024  Changed category.  (WRB)
   891024  REVISION DATE from Version 3.2
   891214  Prologue converted to Version 4.0 format.  (BAB)
   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
   920501  Reformatted the REFERENCES section.  (DWL, WRB)
   920519  Clarified error messages.  (DWL)
   920801  Declarations section rebuilt and code restructured to use
           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
***END PROLOGUE  SSORT
*/

#define TYPEOFY int 	/* DJE--To make the variable type of Y customizable */
			/* because it has to be changed in two places...*/

void ssort(float X[], TYPEOFY Y[], int N, int KFLAG)
/* Note that the second array is an int array.   If you want otherwise,
alter the type above */
{
     /* .. Local Scalars .. */
      float R, T, TT;
      TYPEOFY TTY, TY;
      int I, IJ, J, K, KK, L, M, NN;
     /* .. Local Arrays .. */
      int IL[31], IU[31];  /* DJE--These were 21, but I suspect that this
		limits N to 2^21.  Since memory is cheap, I'll set this a
		little higher */

      /* ***FIRST EXECUTABLE STATEMENT  SSORT */
      NN = N;
      if (NN < 1) {
	 fprintf(stderr,"The number of values to be sorted is not positive.");
	 abort();
      }

      KK = abs(KFLAG);
      if (KK != 1 && KK != 2) {
	 fprintf(stderr,"The sort control parameter, K, is not 2, 1, -1, or -2.");
	 abort();
      }

     /* Alter array X to get decreasing order if needed */

      if (KFLAG <= -1) 
	 for (I=1; I<=NN; I++)
            X[I] = -X[I];
      

      if (KK == 2) goto line100;

     /* Sort X only */

      M = 1;
      I = 1;
      J = NN;
      R = 0.375E0;

line20: if (I == J) goto line60;
      if (R <= 0.5898437E0) 
         R = R+3.90625E-2;
      else R = R-0.21875E0;
      

line30: K = I;

     /* Select a central element of the array and save it in location T */

      IJ = I + (int)((J-I)*R);
      T = X[IJ];

     /* If first element of array is greater than T, interchange with T */

      if (X[I] > T) {
         X[IJ] = X[I];
         X[I] = T;
         T = X[IJ];
      }
      L = J;

     /* If last element of array is less than than T, interchange with T */

      if (X[J] < T) {
         X[IJ] = X[J];
         X[J] = T;
         T = X[IJ];

        /* If first element of array is greater than T, interchange with T */

         if (X[I] > T) {
            X[IJ] = X[I];
            X[I] = T;
            T = X[IJ];
         }
      }

     /* Find an element in the second half of the array which is smaller */
     /* than T */

line40: L = L-1;
      if (X[L] > T) goto line40;

     /* Find an element in the first half of the array which is greater */
     /* than T */

line50: K = K+1;
      if (X[K] < T) goto line50;

     /* Interchange these elements */

      if (K <= L) {
         TT = X[L];
         X[L] = X[K];
         X[K] = TT;
         goto line40;
      }

     /* Save upper and lower subscripts of the array yet to be sorted */

      if (L-I > J-K) {
         IL[M] = I;
         IU[M] = L;
         I = K;
         M = M+1;
      } else {
         IL[M] = K;
         IU[M] = J;
         J = L;
         M = M+1;
      }
      goto line70;

     /* Begin again on another portion of the unsorted array */

line60: M = M-1;
      if (M == 0) goto line190;
      I = IL[M];
      J = IU[M];

line70: if (J-I >= 1) goto line30;
      if (I == 1) goto line20;
      I = I-1;

line80: I = I+1;
      if (I == J) goto line60;
      T = X[I+1];
      if (X[I] <= T) goto line80;
      K = I;

line90: X[K+1] = X[K];
      K = K-1;
      if (T < X[K]) goto line90;
      X[K+1] = T;
      goto line80;

     /* Sort X and carry Y along */

line100: M = 1;
      I = 1;
      J = NN;
      R = 0.375E0;

line110: if (I == J) goto line150;
      if (R <= 0.5898437E0) 
         R = R+3.90625E-2;
      else R = R-0.21875E0;

line120: K = I;

     /* Select a central element of the array and save it in location T */

      IJ = I + (int)((J-I)*R);
      T = X[IJ];
      TY = Y[IJ];

     /* If first element of array is greater than T, interchange with T */

      if (X[I] > T) {
         X[IJ] = X[I];
         X[I] = T;
         T = X[IJ];
         Y[IJ] = Y[I];
         Y[I] = TY;
         TY = Y[IJ];
      }
      L = J;
;
     /* If last element of array is less than T, interchange with T */

      if (X[J] < T) {
         X[IJ] = X[J];
         X[J] = T;
         T = X[IJ];
         Y[IJ] = Y[J];
         Y[J] = TY;
         TY = Y[IJ];

        /* If first element of array is greater than T, interchange with T */

         if (X[I] > T) {
            X[IJ] = X[I];
            X[I] = T;
            T = X[IJ];
            Y[IJ] = Y[I];
            Y[I] = TY;
            TY = Y[IJ];
         }
      }

     /* Find an element in the second half of the array which is smaller */
     /* than T */

line130: L = L-1;
      if (X[L] > T) goto line130;

     /* Find an element in the first half of the array which is greater */
     /* than T */

line140: K = K+1;
      if (X[K] < T) goto line140;

     /* Interchange these elements */

      if (K <= L) {
         TT = X[L];
         X[L] = X[K];
         X[K] = TT;
         TTY = Y[L];
         Y[L] = Y[K];
         Y[K] = TTY;
         goto line130;
      }

     /* Save upper and lower subscripts of the array yet to be sorted */

      if (L-I > J-K) {
         IL[M] = I;
         IU[M] = L;
         I = K;
         M = M+1;
      } else {
         IL[M] = K;
         IU[M] = J;
         J = L;
         M = M+1;
      }
      goto line160;

     /* Begin again on another portion of the unsorted array */

line150: M = M-1;
      if (M == 0) goto line190;
      I = IL[M];
      J = IU[M];

line160: if (J-I >= 1) goto line120;
      if (I == 1) goto line110;
      I = I-1;

line170: I = I+1;
      if (I == J) goto line150;
      T = X[I+1];
      TY = Y[I+1];
      if (X[I] <= T) goto line170;
      K = I;

line180: X[K+1] = X[K];
      Y[K+1] = Y[K];
      K = K-1;
      if (T < X[K]) goto line180;
      X[K+1] = T;
      Y[K+1] = TY;
      goto line170;

     /* Clean up */

line190: if (KFLAG <= -1) 
	 for (I=1; I<=NN; I++)
            X[I] = -X[I];
      
     return;
}
