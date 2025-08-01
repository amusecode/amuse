/* REGROUP.C, Daniel Eisenstein, 1997 */
/* Based on a paper by Daniel Eisenstein & Piet Hut, 
"HOP: A New Group-Finding Algorithm for N-body Simulations."
See the included documentation or view it at 
http://www.sns.ias.edu/~eisenste/hop/hop_doc.html */

/* Version 1.0 (12/15/97) -- Original Release */
/* Version 1.1 (04/02/02) -- Changed the type of variable 'block' in
			     densitycut() function from float to int */

#include "slice.h"
#include <string.h>
#include <stdio.h>
#include <limits.h>

#define MINDENS (-FLT_MAX/3.0)
/* This is the most negative density that can be accomodated.  Note
that MINDENS*2.0 is referenced in the code and so must be properly
represented by the machine.  There's no reason for this to be close to
the actual minimum of the density. */

#define INFORM(pstr) printf(pstr); fflush(stdout)
/* Used for messages, e.g. INFORM("Doing this"); */

#define UNBOUND -2      /* The tag marker for unbound particles */

/* ----------------------------------------------------------------------- */
/* The following structures track all the information about the groups */

typedef struct groupstruct {
    int npart;          /* Number of particles in the group */
    int npartcum;       /* Cumulative number of particles */
    int nread;          /* Number read so far, also a utility field */
    double compos[3], comvel[3];/* Lists of group CoM position and velocities */
    double comtemp[3];  /* Temporary CoM position */
    int idmerge;        /* The new group ID # after merging.  Not necessarily
                                unique! */
    int rootgroup;	/* The fully traced group id */
} Group;  /* Type Group is defined */

typedef struct groupliststruct {
    int npart;          /* Number of particles in the simulation */
    int ngroups;        /* Number of groups in list */
    int nnewgroups;     /* Number of groups after relabeling */
    int npartingroups;  /* Number of particles in groups */
    Group *list;        /* List of groups, zero-offset */
} Grouplist; /* Type Grouplist is defined */

/* ----------------------------------------------------------------------- */
/* Prototypes */
void initgrouplist(Grouplist *g);
void readtags(Slice *s, Grouplist *g, char *fname);
void densitycut(Slice *s, char *fname, float densthresh);
void writegmerge(Slice *s, Grouplist *gl, char *fname, float pt, float mt);
void readgmerge(Slice *s, Grouplist *gl, char *fname);
void merge_groups_boundaries(Slice *s, Grouplist *gl, char *fname,
	float peakdensthresh, float saddledensthresh, float densthresh);
void translatetags(Slice *s, Grouplist *gl);
void writetags(Slice *s, Grouplist *gl, char *fname);
void writetagsf77(Slice *s, Grouplist *gl, char *fname);
void count_membership(Slice *s, Grouplist *g);
void sort_groups(Slice *s, Grouplist *gl, int mingroupsize, char *fname);

/* ----------------------------------------------------------------------- */
/* We use the following structure to handle the user interface: */

typedef struct controlstruct {
    char *tagname;	/* Input file for group tags */
    char *densname;	/* Input file for density file */
    char *gmergename;	/* Input file for group boundary specifications, OR
				input file for group merging data */
    char *outsizename;	/* Output file for size output*/
    char *outtagname;	/* Output file for group tags */
    char *outgmergename;	/* Output file for group merging */

    int qdenscut;	/* =1 if we're making a density cut, =0 otherwise */
    float densthresh;	/* The outer density threshold (delta_outer)*/

    int qgbound;	/* =1 if we are to read the boundaries file and 
				determine the merging.*/
    float peak_thresh;	/* Density threshold for peak (delta_peak) */
    float saddle_thresh;  /* Density threshold for merging (delta_saddle) */
    int qgmerge_given;	/* =1 if we are to use a group translation from file */

    int mingroupsize;	/* The minimum group size we follow */
    int qoutput;	/* =1 if we are to write the tags */
    int qf77;		/* =1 if binary output if in f77 format */
    int qpipe;		/* =1 if we are to write the output tags to stdout */
    int qsort;		/* =1 if we are to sort */

    /* The following aren't used in the present version, but I included
    them in case the user wants to customize the program: */
    char *dataname;	/* Input file for particle data */
    int qunbind;	/* =1 if we are to unbind at all */
} Controls;	/* Type Controls is defined */

/* ====================================================================== */
/* ===================== User Interface ================================= */
/* ====================================================================== */

void usage()
{
    fprintf(stderr, "USAGE for regroup():\n");
    fprintf(stderr,"regroup [-hop <file_name>] -- The .hop file: input group memberships\n");
    fprintf(stderr,"        [-den <file_name>] -- The .den file: density at each particle\n");
    fprintf(stderr,"        [-gbound <file_name>] -- The .gbound file: the boundaries between groups\n");
    fprintf(stderr,"        [-root <file_name>] -- Set the above 3 to root + default suffixes.\n\n");
    fprintf(stderr,"        [-douter <float>] -- Density required for a particle to be in a group.\n");
    fprintf(stderr,"        [-nodens] -- Turn off the density cut.\n\n");
    fprintf(stderr,"        [-dpeak <float>] -- Density required at group center\n");
    fprintf(stderr,"        [-dsaddle <float>] -- Density required at boundary\n");
    fprintf(stderr,"        [-gmerge <file_name>] -- Read the merging from an existing file.\n");
    fprintf(stderr,"        [-nomerge] -- Turn off group merging.\n\n");
    fprintf(stderr,"        [-mingroup <int>] -- Minimum group size.\n");
    fprintf(stderr,"        [-nosort] -- Don't sort the groups by size order.\n\n");
    fprintf(stderr,"        [-o/-out <file_name>] -- Output file root for .gmerge, .size, .tag files\n");
    fprintf(stderr,"        [-otag/-outtag <file_name>] -- Give a special name for the tag output.\n");
    fprintf(stderr,"        [-notagoutput] -- Don't output the tags.\n");
    fprintf(stderr,"        [-pipe] -- Write the tag file to stdout.\n");
    fprintf(stderr,"        [-pipequiet] -- As -pipe, but don't write .gmerge or .size files\n");
    fprintf(stderr,"        [-f77] -- Write tag output in f77 format.\n");
    exit(1); return;
}

void parsecommandline(int argc, char *argv[], Controls *c)
{
    int narg, qmerge;
    char *outname, *rootname;
    narg = 1;
    rootname = c->dataname = c->densname = c->gmergename = c->tagname = 
	outname = c->outsizename = c->outtagname = c->outgmergename = NULL;
    c->qdenscut = -1;
    qmerge = 1;
    c->qgmerge_given = 0;

    c->qunbind = 0; 
    c->qoutput = 1;
    c->qsort = 1;
    c->qpipe = 0;
    c->qf77 = 0; 

    c->mingroupsize = -1;
    if (2.0*MINDENS>=MINDENS || MINDENS>=0) 
	myerror("MINDENS seems to be illegal.");  
	/* Need MINDENS<0 and 2*MINDENS to be machine-representable */
    c->densthresh = 2.0*MINDENS;
    c->saddle_thresh = 2.0*MINDENS; 	
    c->peak_thresh = 2.0*MINDENS; 

    if (argc==1) usage();	/* No arguments is a request for help */
    while (narg<argc) {
	if (!strcmp(argv[narg],"-in")) {
	    myerror("Particle data is not needed unless program is customized.");
	    narg++; if (narg>=argc) usage();
	    c->dataname = argv[narg]; narg++;
	}
	if (!strcmp(argv[narg],"-root")) {
	    narg++; if (narg>=argc) usage();
	    rootname = argv[narg];
	    narg++;
	}
	else if (!strcmp(argv[narg],"-hop")) {
	    narg++; if (narg>=argc) usage();
	    c->tagname = argv[narg]; narg++;
	}
	else if (!strcmp(argv[narg],"-dens")||!strcmp(argv[narg],"-den")) {
	    narg++; if (narg>=argc) usage();
	    c->densname = argv[narg]; narg++;
	}
	else if (!strcmp(argv[narg],"-gbound")) {
	    narg++; if (narg>=argc) usage();
	    if (c->gmergename!=NULL) 
		myerror("You shouldn't use both -gbound and -gmerge.");
	    c->gmergename = argv[narg]; narg++;
	    c->qgmerge_given = 0; 
	}
	else if (!strcmp(argv[narg],"-gmerge")) {
	    narg++; if (narg>=argc) usage();
	    if (c->gmergename!=NULL) 
		myerror("You shouldn't use both -gbound and -gmerge.");
	    c->gmergename = argv[narg]; narg++;
	    c->qgmerge_given = 1; 
	}
	else if (!strcmp(argv[narg],"-o") || !strcmp(argv[narg],"-out")) {
	    narg++; if (narg>=argc) usage();
	    outname = argv[narg]; narg++;
	}
	else if (!strcmp(argv[narg],"-otag")||!strcmp(argv[narg],"-outtag")) {
	    narg++; if (narg>=argc) usage();
	    c->outtagname = argv[narg]; narg++;
	}
	else if (!strcmp(argv[narg],"-douter")) {
	    narg++; if (narg>=argc) usage();
	    sscanf(argv[narg],"%f",&(c->densthresh)); narg++;
	    if (c->qdenscut!=0) c->qdenscut=1;
	}
	else if (!strcmp(argv[narg],"-dsaddle")) {
	    narg++; if (narg>=argc) usage();
	    sscanf(argv[narg],"%f",&(c->saddle_thresh)); narg++;
	}
	else if (!strcmp(argv[narg],"-dpeak")) {
	    narg++; if (narg>=argc) usage();
	    sscanf(argv[narg],"%f",&(c->peak_thresh)); narg++;
	}
	else if (!strcmp(argv[narg],"-mingroup")) {
	    narg++; if (narg>=argc) usage();
	    sscanf(argv[narg],"%d",&(c->mingroupsize)); narg++;
	}
	else if (!strcmp(argv[narg],"-nodens")) {
	    c->qdenscut= 0; narg++;
	}
	else if (!strcmp(argv[narg],"-nomerge")) {
	    c->qgmerge_given = qmerge = 0; narg++; 
	}
	else if (!strcmp(argv[narg],"-nosort")) {
	    c->qsort= 0; narg++;
	}
	else if (!strcmp(argv[narg],"-notagoutput")) {
	    c->qoutput=0; narg++;
	}
	else if (!strcmp(argv[narg],"-nounbind")) {
	    myerror("The program isn't customized to use this feature.");
	    c->qunbind=0; narg++;
	}
	else if (!strcmp(argv[narg],"-unbind")) {
	    myerror("The program isn't customized to use this feature.");
	    c->qunbind=1; narg++;
	}
	else if (!strcmp(argv[narg],"-pipe")) {
	    c->qpipe = 1; narg++;
	}
	else if (!strcmp(argv[narg],"-pipequiet")) {
	    c->qpipe = -1; narg++;
	}
	else if (!strcmp(argv[narg],"-f77")) {
	    c->qf77 = 1; narg++;
	}
	else usage();
    }

    /* Get the input files ready */
    if (c->qdenscut==-1) {
	/* Neither -douter nor -nodens was chosen. */
	mywarn("Outer density threshold left unspecified.  Skipping this cut.");
	c->qdenscut = 0;
    } else if (c->qdenscut==1) {
	/* We have a chosen density.  Need to figure out the density file. */
	if (c->densname==NULL) {
	    if (rootname==NULL)
		myerror("No density file name or root has been specified.");
	    c->densname = (char *)malloc(80);
	    strcpy(c->densname,rootname); strcat(c->densname, ".den");
	}
    } else c->densname = NULL;	/* We have no reason to read it */

    if (c->tagname==NULL) {
	if (rootname==NULL)
	    myerror("No .hop file name or root has been specified.");
	c->tagname = (char *)malloc(80);
	strcpy(c->tagname,rootname); strcat(c->tagname, ".hop");
    }

    if (qmerge==1) {
	if (c->qgmerge_given==0) {  
	    /* We need to have a .gbound file */
	    c->qgbound = 1;
	    if (c->saddle_thresh<MINDENS || c->peak_thresh<MINDENS)
		myerror("-dsaddle and -dpeak need to be specified.");
	    if (c->gmergename==NULL) {
		if (rootname==NULL)
		    myerror("No .gbound file name or root has been specified.");
		c->gmergename = (char *)malloc(80);
		strcpy(c->gmergename,rootname); 
		strcat(c->gmergename, ".gbound");
	    }
	} else c->qgbound = 0;    /* We know c->mergename is ready to go */
    } else c->gmergename = NULL;  /* No reason to read it */

    /* Get the output files ready */
    /* If a default name wasn't given, we'll assume zregroup */
    if (outname==NULL) {
	outname = (char *)malloc(20);
	strcpy(outname,"zregroup");
    }
    /* Primary tag output: */
    if (c->qoutput) {  /* Need to figure out where we're sending the output */
	if (c->qpipe&&c->outtagname)
	    myerror("Conflicting instructions--gave specific output name and told to pipe.");
	if (c->qpipe>0) mywarn("Writing tags to stdout.");
	if (c->qpipe) c->outtagname = NULL;  /* Our signal to send to stdout */
	else if (c->outtagname==NULL) {
	    c->outtagname = (char *)malloc(80);
	    strcpy(c->outtagname, outname);
	    strcat(c->outtagname, ".tag");
	} /* Otherwise the name was set by the user */
    } else {
	/* We're not outputing tags */
	if (c->qpipe) myerror("Conflicting instructions--told to pipe and not to output.");
    }

    if (c->qsort) {
	if (c->qpipe>=0) {	/* The user didn't specify quiet */
	    c->outsizename = (char *)malloc(80);
	    strcpy(c->outsizename, outname);
	    strcat(c->outsizename, ".size");
	}
    }

    if (c->qpipe>=0) {	/* The user didn't specify quiet */
	c->outgmergename = (char *)malloc(80);
	strcpy(c->outgmergename, outname);
	strcat(c->outgmergename, ".gmerge");
    }

    if (c->mingroupsize >= 0 && !c->qsort)
	myerror("Imposition of a certain group size occurs within the sort routine.");
    if (c->qsort && c->mingroupsize < 0) {
	mywarn("No minimum group size specified.  Assuming 10 particles.");
	c->mingroupsize = 10;
    }

    if (c->densthresh<MINDENS) c->densthresh=MINDENS;
	/* This is our default--a very negative number */

    return;
}

/* ====================================================================== */
/* ============================== MAIN() ================================ */
/* ====================================================================== */

void main(int argc, char *argv[])
{
    Grouplist gl;
    Slice *s;
    FILE *f;
    Controls c;

    parsecommandline(argc, argv, &c);

    initgrouplist(&gl);
    s=newslice();

    /* We need to read the tag file and perhaps perform a density cut */
    readtags(s,&gl,c.tagname);
    if (c.qdenscut) densitycut(s,c.densname,c.densthresh);

    /* Next do the merging of the input groups */
    if (c.qgbound) {
	/*  We're going to read a .gbound file and merge groups */
	merge_groups_boundaries(s,&gl,c.gmergename,
		c.peak_thresh, c.saddle_thresh, c.densthresh);
	/* Renumber the groups from large to small; remove any tiny ones */
	if (c.qsort) sort_groups(s, &gl, c.mingroupsize, c.outsizename);
	writegmerge(s, &gl, c.outgmergename, c.peak_thresh, c.saddle_thresh);
	translatetags(s,&gl);
    }
    else if (c.qgmerge_given) {
	/* We're going to read a .gmerge file and merge groups as it says */
 	readgmerge(s, &gl, c.gmergename);
 	translatetags(s, &gl);
    } /* Else we'll use the tags as given by the original .hop file */

    /* If one wants to manipulate the groups any more, this is a good 
    place to do it.  For example, you might want to remove unbound particles:
	if (c.qunbind) {
	    get_particle_data(s, &gl, c.dataname);
	    unbind_particles(s, &gl, c.mingroupsize);
	}
    */

    /* Write the output */
    if (c.qoutput) {
	if (c.qf77) writetagsf77(s, &gl, c.outtagname);
	else writetags(s, &gl, c.outtagname);
    }
    
    free_slice(s);
    return;
}

/* ================================================================= */
/* =================== Initialization Routines ===================== */
/* ================================================================= */

void initgrouplist(Grouplist *g)
/* Just make sure this stuff is zero */
{
    g->list = NULL;
    g->npartingroups = g->npart = g->ngroups = 0; g->nnewgroups = 0;
    return;
}

void readtags(Slice *s, Grouplist *g, char *fname)
/* Read the tag file named fname into s->ntag[] */
/* Groups need not be sorted, but must be numbered from 0 to ngroups-1 */
{
    FILE *f;

    if ((f=fopen(fname,"r"))==NULL) myerror("Input tag file not found.");
    if (fread(&(g->npart),4,1,f)!=1) myerror("Tag file read error.");
    if (fread(&(g->ngroups),4,1,f)!=1) myerror("Tag file read error.");
    fprintf(stderr,"Number of particles: %d.   Number of groups: %d.\n",
	g->npart, g->ngroups);

    s->numpart = g->npart;
    s->numlist = g->npart;
    s->ntag = ivector(1,s->numlist);
    fread(s->ntag+1, 4, s->numlist, f);	/* Read in all the tags */
    fclose(f);
    return;
}

/* ========================== Density Cut ======================== */

#define MAXBLOCK 65536 		/* Read the file 256k at a time */

void densitycut(Slice *s, char *fname, float densthresh)
/* Read the density file and change the tag on any particle with density
less than densthresh to -1, thus removing them from groups */
/* This will leave some groups with no particles, which is fine */
/* We read the file in segments, so as to reduce memory consumption */
{
    FILE *f;
    int j, numread, npart, block;	/* block was a float by mistake */
    float density[MAXBLOCK];

    if ((f=fopen(fname,"r"))==NULL)
	myerror("Density file not found.");
    npart = 0; fread(&npart,4,1,f);
    if (npart!=s->numpart) 
	mywarn("Density file doesn't match slice description.");

    numread = 0;
    block = MAXBLOCK;	/* Start off big */
    while (numread<npart) {
	if (npart-numread<block) block = npart-numread;
	if (fread(density,4,block,f)!=block)
	    myerror("Read error in density file.");
	for (j=1;j<=block;j++) 
	    if (density[j-1]<densthresh)	/* density is zero-offset */
		s->ntag[numread+j]=(-1);	/* s->ntag is unit-offset */
	numread+=block;
    }
    fclose(f);
    return;
}

/* ====================== Read/Write .gmerge files ======================= */
/* The gmerge file is just a map from the old (pre-regroup) group numbers 
to the new (post-regroup) group numbers.  Of course, there are more "old"
groups than "new" groups, since the point of regroup() is to merge groups. */

void writegmerge(Slice *s, Grouplist *gl, char *fname, float pt, float mt)
/* Write the translation between old groups and new groups, ASCII */
{
    FILE *f;
    int j;
    Group *gr;

    if (fname==NULL) return; /* We've been told not to write anything */

    if ((f=fopen(fname,"w"))==NULL) myerror("Can't open gmerge file for write.");
    fprintf(f,"%d\n%d\n%d\n", gl->npart, gl->ngroups, gl->nnewgroups);
    fprintf(f,"%f\n%f\n", pt, mt);
    for (j=0,gr=gl->list;j<gl->ngroups;j++,gr++) 
	 fprintf(f,"%d %d\n", j, gr->idmerge);
    fclose(f);
    return;
}

void readgmerge(Slice *s, Grouplist *gl, char *fname)
/* Read the translation between old groups and new groups, ASCII */
/* Also, set up gl->list for translation */
{
    FILE *f;
    int j, dummy;
    Group *gr;
    float pt, mt;

    if ((f=fopen(fname,"r"))==NULL) myerror("Can't open gmerge read file.");
    if (fscanf(f,"%d\n%d\n%d\n", &(gl->npart), &(gl->ngroups), 
	&(gl->nnewgroups))!=3) myerror("Error in header of gmerge file.");
    if (gl->npart!=s->numpart) myerror("Number of particles in gmerge file doesn't match that of tags file.");
    fscanf(f,"%f %f\n", &pt, &mt);

    if (gl->list!=NULL) free(gl->list);
    gl->list = (Group *)malloc((size_t)(gl->ngroups *sizeof(Group)));
    if (gl->list==NULL) myerror("Error in allocating gl->list.");

    for (j=0,gr=gl->list; j<gl->ngroups; j++,gr++) {
	if (fscanf(f,"%d %d\n", &dummy, &(gr->idmerge))!=2 || dummy!=j)
		myerror("Error in reading gmerge file.");
	gr->npart = -1;	/* We're not setting this */
    }
    fclose(f);
    return;
}

/* ====================== GROUP MERGING BY BOUNDARIES ================ */

void merge_groups_boundaries(Slice *s, Grouplist *gl, char *mergename,
	float peakdensthresh, float saddledensthresh, float densthresh) 
/* Read in the gmerge file and decide which groups are to be merged.
Groups are numbered 0 to ngroups-1.  Groups with boundaries greater
than saddledensthresh are merged.  Groups with maximum densities
less than peakdensthresh are merged to the group with
maxdensity above peakdensthresh with which it shares the highest
density border. */
/* Only groups with maximum densities above peakdensthresh can be group
centers. */
/* Allocate space for the grouplist and store the merging results in
the idmerge field. */
/* I think this will work even if saddledensthresh<densthresh */
{
    int j, k, g1, g2, ngroups, dummy[3];
    Group *gr;
    float *gdensity, *densestbound, fdum[3], dens;
    int *densestboundgroup, changes;
    char line[80], *tempfilename; 
    FILE *fp;
    FILE *boundfp;

    tempfilename = tmpnam(NULL);
    if (densthresh<MINDENS) densthresh=MINDENS;	
	/* Have a 2*MINDENS condition below... */
    if ((fp=fopen(mergename,"r"))==NULL) myerror("Error opening gmerge file.");
    if (fgets(line,80,fp)==NULL) myerror("Unexpected EOF in gmerge.");
    if (sscanf(line, "%d", &ngroups)!=1) myerror("Error in gmerge format.");
    gdensity = vector(0,ngroups-1);
    densestbound = vector(0,ngroups-1);
    densestboundgroup = ivector(0,ngroups-1);

    for (j=0;;) {
	if (fgets(line, 80, fp)==NULL) myerror("Unexpected EOF in gmerge.");
	if (line[0]=='#' && line[1]=='#' && line[2]=='#') break;
			/* This is the flag for the end of the group list */
	if (line[0]=='#') continue;	/* A comment line */
	if (j>=ngroups) myerror("Group list did not end where expected.");
	if (sscanf(line, "%d %d %d %f %f %f %f", dummy, dummy+1, dummy+2,
		fdum, fdum+1, fdum+2, gdensity+j)!=7)
		myerror("Error parsing gmerge line.");
	if (dummy[0]!=j) myerror("Group list format is improper.");
	if (gdensity[j]<0.0) myerror("Negative density read.");
	j++;
    }

    /* Now allocate the grouplist */
    gl->ngroups = ngroups;
    if (gl->list!=NULL) free(gl->list);
    gl->list = (Group *)malloc((size_t)(gl->ngroups *sizeof(Group)));
    if (gl->list==NULL) myerror("Error in allocating gl->list.");
    for (j=0,gr=gl->list;j<gl->ngroups;j++,gr++) {
	/* If group is too underdense, it cannot be a group center */
	if (gdensity[j]<peakdensthresh) gr->idmerge=(-1);
	    else gr->idmerge = j;
	gr->npart = -1;	/* Not doing anything with this */
	densestbound[j] = 2.0*MINDENS;	/* Initialize */
	densestboundgroup[j] = -1;	/* Initialize */
    }

    /* Now step through the list of boundaries */
    /* If a boundary is between two groups with max densities above 
    peakdensthresh and if the boundary is above saddledensthresh, then
    merge the groups (keeping the lower number of the two). */
    /* If one of the groups is above peakdensthresh and the other is
    below, and if the boundary density is higher than any seen previously 
    for the lower density group, then record this information */
    /* If neither group is above peakdensthresh, skip the boundary */

    if ((boundfp=fopen(tempfilename,"w"))==NULL)
	myerror("Error opening scratch file");

    while (fgets(line,80,fp)!=NULL) {
	if (sscanf(line,"%d %d %f", &g1, &g2, &dens)!=3 || g1<0 || g2<0 || 
		dens<0.0 || g1>=ngroups || g2>=ngroups) 
		myerror("Error reading boundary.");
	if (gdensity[g1]<peakdensthresh && gdensity[g2]<peakdensthresh) {
	    if (gdensity[g1]>densthresh && gdensity[g2]>densthresh &&
		    dens>densthresh) 
		fprintf(boundfp,"%d %d %f\n", g1, g2, dens);
	    continue;  	/* group isn't dense enough */
	}
	if (gdensity[g1]>=peakdensthresh && gdensity[g2]>=peakdensthresh)
	    if (dens<saddledensthresh) continue;
		/* Boundary is not dense enough to merge */
	    else {	/* Groups should be merged */
		/* Trace each group to its root */
		while (g1!=gl->list[g1].idmerge)
		    g1=gl->list[g1].idmerge;
		while (g2!=gl->list[g2].idmerge)
		    g2=gl->list[g2].idmerge;
		if (g1<g2) gl->list[g2].idmerge=g1;
		else gl->list[g1].idmerge=g2;
		continue;	/* Go to the next boundary */
	    }
	/* Else one is above peakdensthresh, the other below.   */
	/* Make the high one g1 */
	if (gdensity[g1]<gdensity[g2]) {
	    dummy[0] = g2; g2=g1; g1=dummy[0];
	}
	if (dens>densestbound[g2]) {
	    /* It's the densest boundary yet */
	    densestbound[g2] = dens;
	    densestboundgroup[g2] = g1;
	}
    } /* Get the next boundary line */
    fclose(fp);
    fclose(boundfp);

    /* Now the fringe groups are connected to the proper group 
    (>peakdensthresh) with the largest boundary.  But we want to look
    through the boundaries between fringe groups to propagate this
    along.  Connections are only as good as their smallest boundary */
    /* Keep the density of the connection in densestbound, and the
    proper group it leads to in densestboundgroup */
    do {
	if ((boundfp=fopen(tempfilename,"r"))==NULL)
		myerror("Error opening scratch file for read");
	changes = 0;
	while (fgets(line,80,boundfp)!=NULL) {
	    if (sscanf(line,"%d %d %f", &g1, &g2, &dens)!=3)
		myerror("Error reading boundary.");
	    /* If the density of this boundary and the densestbound of
	    the other group is higher than a group's densestbound, then
	    replace it. */
	    /* Make g1 have the higher densestbound */
	    if (densestbound[g2]>densestbound[g1]) {
		dummy[0] = g2; g2=g1; g1=dummy[0];
	    }
	    if (dens>densestbound[g2]&&densestbound[g1]>densestbound[g2]) {
		changes++;
		if (dens<densestbound[g1]) densestbound[g2]=dens;
		    else densestbound[g2]=densestbound[g1];
		densestboundgroup[g2] = densestboundgroup[g1];
	    }
	}
	fclose(boundfp);
    } while (changes);

    /* Now connect the low-density groups to their densest boundaries */
    /* But only if the boundary exceeds densthresh! */
    for (j=0;j<gl->ngroups;j++) 
	if (densestbound[j]>=densthresh)
	    gl->list[j].idmerge = densestboundgroup[j];
    
    /* Now we want to number the newly merged groups */
    /* The center groups are given negative numbers <-1 */
    for (j=0,gl->nnewgroups=0; j<gl->ngroups; j++)
	if (gl->list[j].idmerge==j) 
	    gl->list[j].idmerge = -2-(gl->nnewgroups++);

    /* Now trace each group through until a negative number is reached */
    for (j=0; j<gl->ngroups; j++) {
	if (gl->list[j].idmerge<0) continue;
	g1 = j;
	while ((g1=gl->list[g1].idmerge)>=0);
	g2 = j;
	do gl->list[g2].idmerge = g1;
	    while ((g2=gl->list[g2].idmerge)>=0);
    }

    /* Finally, renumber the groups 0..N-1 */
    for (j=0,gr=gl->list;j<gl->ngroups;j++,gr++)
	gr->idmerge = -2-gr->idmerge;	/* Keep -1 -> -1 */
    
    /* And delete the tempfile */
    remove(tempfilename);
    free_vector(gdensity,0,ngroups-1);
    free_vector(densestbound,0,ngroups-1);
    free_ivector(densestboundgroup,0,ngroups-1);
    return;
}

/* ======================================================================= */
/* =============== Update the tags and write them out ==================== */
/* ======================================================================= */

void translatetags(Slice *s, Grouplist *gl)
/* Alter s->ntag to have the new groups.  Reset gl so as to reflect the
new number of groups. */
{
    int j;

    for (j=1;j<=s->numlist;j++) 
	if (s->ntag[j]>=0)
	    s->ntag[j] = gl->list[s->ntag[j]].idmerge;
	/* Otherwise, translate the unbound particles */
	else if (s->ntag[j]<-1)
	    s->ntag[j] = UNBOUND - gl->list[UNBOUND-s->ntag[j]].idmerge;
    free(gl->list);
    gl->list = NULL;
    gl->ngroups = gl->nnewgroups;
    return;
}

void writetags(Slice *s, Grouplist *gl, char *fname)
/* Write s->ntag to file */
/* If fname==NULL, write to stdout */
{
    FILE *f;

    if (fname!=NULL) {
	if ((f=fopen(fname,"w"))==NULL) myerror("Error opening new tag file.");
    } else f=stdout;
    fwrite(&(s->numpart),4,1,f);
    fwrite(&(gl->ngroups),4,1,f);
    fwrite(s->ntag+1,4,s->numlist,f);
    fclose(f);
    return;
}

void writetagsf77(Slice *s, Grouplist *gl, char *fname)
/* Write s->ntag to file */
/* If fname==NULL, write to stdout */
/* Use a format readable for FORTRAN unformatted read commands */
{
    FILE *f;
    int dummy;
    if (fname!=NULL) {
	if ((f=fopen(fname,"w"))==NULL) myerror("Error opening new tag file.");
    } else f=stdout;
    dummy = 8; fwrite(&dummy,4,1,f);
    fwrite(&(s->numpart),4,1,f);
    fwrite(&(gl->ngroups),4,1,f);
    fwrite(&dummy,4,1,f);
    dummy = s->numlist*4; fwrite(&dummy,4,1,f);
    fwrite(s->ntag+1,4,s->numlist,f);
    fwrite(&dummy,4,1,f);
    fclose(f);
    return;
}

/* ====================================================================== */
/* ========================== Sorting the Groups ======================== */
/* ====================================================================== */

void sort_groups(Slice *s, Grouplist *gl, int mingroupsize, char *fname)
/* Sort the groups, as labeled by the idmerge field not their original 
number, from largest to smallest.  Alter the idmerge field to this new
numbering, setting any below mingroupsize to -1. */
/* If fname!=NULL, write a little output file listing the group sizes */
{
    FILE *f;
    int j,k, *order, partingroup, igr, *newnum, nmergedgroups;
    float *gsize;
    Group *c;
    void make_index_table(int n, float *fvect, int *index);

    nmergedgroups = gl->nnewgroups;
    gsize = vector(0,nmergedgroups-1);
    order = ivector(1,nmergedgroups);
    newnum = ivector(0,nmergedgroups-1);

    /* First we need to find the number of particles in each group */
    for (j=0,c=gl->list;j<gl->ngroups;j++,c++) c->npart=0; 

    for (j=1;j<=s->numlist;j++) {	/* Look through all the particles */
	igr = s->ntag[j];
	if (igr>=0)
	    if (igr<gl->ngroups) gl->list[igr].npart++;
	    else myerror("Group tag is out of bounds.");
    }
    /* Now combine these to find the number in the new groups */
    for (j=0;j<nmergedgroups;j++) gsize[j]=0;
    for (j=0,c=gl->list;j<gl->ngroups;j++,c++) 
	if (c->idmerge>=0 && c->idmerge<nmergedgroups)
	    gsize[c->idmerge]+=c->npart;
	else if (c->idmerge>=nmergedgroups)
	    myerror("Group idmerge is out of bounds.");
	    
    make_index_table(nmergedgroups, gsize-1, order);
    /* But remember that order[] thinks that gsize is unit-offset */
    for (j=nmergedgroups,k=0;j>0; j--,k++) 
	if (gsize[order[j]-1]>mingroupsize-0.5) newnum[order[j]-1]=k;
	else break; 	/* All of the rest are too small */
	
    gl->nnewgroups = k;
    for (;j>0;j--) newnum[order[j]-1]=(-1);
    /* Newnum[] holds the new sorted number for merged group j */

    /* Now assign sorted group numbers to idmerge */
    partingroup = 0;
    for (j=0,c=gl->list;j<gl->ngroups;j++,c++) 
	if (c->idmerge>=0)
	    if ((c->idmerge = newnum[c->idmerge])>=0) 
		partingroup+=c->npart;

    /* Output the .size file, if inputed name isn't NULL */
    if (fname!=NULL) {
	f = fopen(fname,"w");
	fprintf(f,"%d\n%d\n%d\n", s->numpart, partingroup, gl->nnewgroups);
	for (j=0;j<gl->nnewgroups;j++)
	    fprintf(f,"%d %d\n", j, (int)gsize[order[nmergedgroups-j]-1]);
    }
    fclose(f);
    free_ivector(order,1,nmergedgroups);
    free_vector(gsize,0,nmergedgroups-1);
    free_ivector(newnum,0,nmergedgroups-1);
    return;
}

/* ======================== Sorting ============================ */

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

void make_index_table(int n, float *fvect, int *index)
/* Given a vector of floats fvect[1..n], construct a index table index[1..n]
so that index[j] contains the ID number of the jth lowest element.
Storage for index[] should be declared externally */
/* This isn't fast, but it takes a tiny fraction of the runtime */
{
    int j;
    ptrindex sortvect;

    sortvect = (ptrindex)malloc(n*sizeof(struct index_struct));
    for (j=0;j<n;j++) sortvect[j].value = fvect[j+1];
    for (j=0;j<n;j++) sortvect[j].index = j+1;  /* Label them prior to sort */
    qsort(sortvect,n,sizeof(struct index_struct),cmp_index);
    /* Now sortvect is in order (smallest to largest) */
    for (j=0;j<n;j++) index[j+1] = sortvect[j].index;
    free(sortvect);
    return;
}
