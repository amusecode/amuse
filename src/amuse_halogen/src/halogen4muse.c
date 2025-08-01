/* halogen.c
**
** Program for generating spherical structures of the alpha-beta-gamma model family
**
** written by Marcel Zemp
**
** This is a stripped-down version of HALOGEN for the MUSE project
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <assert.h>
#include <time.h>
#include "definitions.h"
#include "functions.h"
#include "routines.h"

int main(int argc, char **argv) {

    /*
    ** Variables
    */

    INT i, j, index;
    INT outputgridr, outputgriddf;
    DOUBLE randomseed;
    DOUBLE t0, t1, t2, t3, t4, t5, t6, t7;
    PARTICLE *bh;
    SI *halo;
    GI *gi;
    CHAR FILENAME[STRINGSIZE], INPUTNAME[STRINGSIZE];
    FILE *file;

    randomseed = time(NULL);

    t0 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);

    initialise_fixed_structures(&bh, &halo, &gi);
    set_standard_values_for_parameters(&outputgridr, &outputgriddf, gi, halo);
    sprintf(INPUTNAME,"none");

    initialise_black_hole(bh);
    initialise_parameters(halo);

    /*
    ** Read in and calculate model parameters
    */

    i = 1;
    while (i < argc) {
	/*
	** Halo parameters
	*/
	if (strcmp(argv[i],"-a") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->alpha = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-b") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->beta = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-c") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->gamma = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-M") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->M = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rs") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->rs = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rcutoff") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->rcutoff = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-N") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->N0 = (INT) (atof(argv[i]));
	    i++;
	    }
	/*
	** Black hole parameters
	*/
	else if (strcmp(argv[i],"-MBH") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    bh->mass = atof(argv[i]);
	    i++;
	    }
	/*
	** Model name
	*/
	else if (strcmp(argv[i],"-name") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    sprintf(INPUTNAME,"%s",argv[i]);
	    i++;
	    }
	/*
	** Output parameters
	*/
	else if (strcmp(argv[i],"-ogr") == 0) {
	    outputgridr = 1;
	    i++;
	    }
	else if (strcmp(argv[i],"-ogdf") == 0) {
	    outputgriddf = 1;
	    i++;
	    }
	/*
	** Special parameters
	*/
	else if (strcmp(argv[i],"-randomseed") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    randomseed = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-dorvirexact") == 0) {
            gi->dorvirexact = 1;
            i++;
            }
	/*
	** Help or failure
	*/
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
	    usage();
	    }
	else {
	    usage();
	    }
	}

    fprintf(stderr,"Checking parameters, calculating halo properties and initialising grid in r... \n");
    srand48(randomseed);

    /*
    ** Check main input parameters
    */

    if (strcmp(INPUTNAME,"none") == 0) {
	fprintf(stderr,"You have not set a name for the output model.\n");
	usage();
	}
    if ((NGRIDR-1) % (NGRIDDF-1) != 0) {
	fprintf(stderr,"Bad choice of NGRIDR and NGRIDDF!\n");
	fprintf(stderr,"These numbers have to fulfill the condition (NGRIDR-1) mod (NGRIDDF-1) == 0.\n");
	usage();
	}

    check_main_parameters(halo);
    calculate_parameters(gi,halo);

    initialise_all_grids(gi, bh, halo, outputgridr, outputgriddf, t0, &t1, FILENAME, INPUTNAME);

    /*
    ** Initialise structure
    */

    initialise_structure(gi,halo);

    /*
    ** Set particle positions
    */

    t2 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle positions... \n",t2-t1);

    set_positions(halo);

    /*
    ** Set particle velocities
    */

    t3 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle velocities... \n",t3-t2);
	
    set_velocities(gi,halo);

    /*
    ** Set remaining attributes
    */

    t4 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting remaining particle attributes... \n",t4-t3);

    set_attributes(gi,halo);

    /*
    ** Calculate a few things and do center of mass correction
    */

    t5 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nCalculating a few things, doing mass scaling and correct center of mass position and velocity... \n",t5-t4);

    double_particles(halo);
    calculate_stuff(gi,bh,halo);

    /*
    ** Write Output
    */

    t6 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nWriting Output... \n",t6-t5);

    sprintf(FILENAME,"%s.IC.ascii",INPUTNAME);
    file = fopen(FILENAME,"w");

    index = 1;
    if (bh->mass != 0) {
	assert(fprintf(file,OFI1" ",index) > 0);
        assert(fprintf(file,OFD5" ",bh->mass) > 0);
	for (j = 0; j < 3; j++) {
	    assert(fprintf(file,OFD6" ",bh->r[j+1]) > 0);
            }
        for (j = 0; j < 3; j++) {
            assert(fprintf(file,OFD6" ",bh->v[j+1]) > 0);
            }
	assert(fprintf(file,"\n") > 0);
	index++;
        }
    for (i = 0; i < halo->N; i++) {
        assert(fprintf(file,OFI1" ",index) > 0);
        assert(fprintf(file,OFD5" ",halo->p[i].mass) > 0);
        for (j = 0; j < 3; j++) {
            assert(fprintf(file,OFD6" ",halo->p[i].r[j+1]) > 0);
            }
        for (j = 0; j < 3; j++) {
            assert(fprintf(file,OFD6" ",halo->p[i].v[j+1]) > 0);
            }
	assert(fprintf(file,"\n") > 0);
        index++;
	}
    fclose(file);
    
    /* 
    ** Print some output in file
    */

    sprintf(FILENAME,"%s.out",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    fprintf(file,"HALOGEN4MUSE by Marcel Zemp / Version 30. August 2007\n\n");
    fprintf(file,"Command line\n\n");
    for (i = 0; i < argc; i++) fprintf(file,"%s ",argv[i]);
    fprintf(file,"\n\n");
    fprintf(file,"Model properties\n\n");
    fprintf(file,"alpha = "OFD1"\n",halo->sp->alpha);
    fprintf(file,"beta  = "OFD1"\n",halo->sp->beta);
    fprintf(file,"gamma = "OFD1"\n",halo->sp->gamma);
    fprintf(file,"rho0  = "OFD3" MU LU^-3\n",halo->sp->rho0);
    fprintf(file,"\n");
    if (bh->mass > 0) {
	fprintf(file,"MBH = "OFD3" MU\n",bh->mass);
	fprintf(file,"\n");
	}
    fprintf(file,"rs      = "OFD3" LU\n",halo->sp->rs);
    fprintf(file,"rhalf   = "OFD3" LU\n",halo->sp->rhalf);
    fprintf(file,"rvir    = "OFD3" LU\n",halo->sp->rvir);
    fprintf(file,"rinner  = "OFD3" LU\n",gi->rinner);
    fprintf(file,"router  = "OFD3" LU\n",gi->router);
    if (halo->sp->rcutoff != SBI) {
        fprintf(file,"rcutoff = "OFD3" LU\n",halo->sp->rcutoff);
        fprintf(file,"rdecay  = "OFD3" LU\n",halo->sp->rdecay);
        }
    fprintf(file,"\n");
    fprintf(file,"M(rs)      = "OFD3" MU\n",Menc(halo->sp->rs,gi));
    fprintf(file,"M(rhalf)   = "OFD3" MU\n",Menc(halo->sp->rhalf,gi));
    fprintf(file,"M(rvir)    = "OFD3" MU\n",Menc(halo->sp->rvir,gi));
    fprintf(file,"M(rinner)  = "OFD3" MU\n",Menc(gi->rinner,gi));
    fprintf(file,"M(router)  = "OFD3" MU\n",Menc(gi->router,gi));
    if (halo->sp->rcutoff != SBI) {
	fprintf(file,"M(rcutoff) = "OFD3" MU\n",Menc(halo->sp->rcutoff,gi));
	}
    fprintf(file,"\n");
    fprintf(file,"Sampling properties\n\n");
    fprintf(file,"|Cr| = "OFD3" MU       Cr = ("OFD4", "OFD4", "OFD4") MU\n",gi->stuff->Cr[0],gi->stuff->Cr[1],gi->stuff->Cr[2],gi->stuff->Cr[3]);
    fprintf(file,"|Cv| = "OFD3" MU TU^-1 Cv = ("OFD4", "OFD4", "OFD4") MU TU^-1\n",gi->stuff->Cv[0],gi->stuff->Cv[1],gi->stuff->Cv[2],gi->stuff->Cv[3]);
    fprintf(file,"\n");
    fprintf(file,"Etot = "OFD4" MU LU^2 TU^-2\n",gi->stuff->Etot);
    fprintf(file,"Ekin = "OFD4" MU LU^2 TU^-2\n",gi->stuff->Ekin);
    fprintf(file,"Epot = "OFD4" MU LU^2 TU^-2\n",gi->stuff->Epot);
    fprintf(file,"Rvir = |2*Ekin/Epot| = %g\n",fabs(2*gi->stuff->Ekin/gi->stuff->Epot));
    fprintf(file,"\n");
    fprintf(file,"Ntot                = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->N,gi->stuff->N);
    fprintf(file,"rimp                = "OFD3" LU\n",halo->rimp);
    fprintf(file,"r1                  = "OFD3" LU\n",halo->r1);
    fprintf(file,"r100                = "OFD3" LU\n",halo->r100);
    fprintf(file,"Mtheo               = "OFD3" MU\n",bh->mass + halo->sp->M);
    fprintf(file,"Msamp               = "OFD3" MU\n",gi->stuff->Mp);
    fprintf(file,"(Msamp-Mtheo)/Mtheo = "OFD3"\n",gi->stuff->Mp/(bh->mass + halo->sp->M)-1.0);
    fprintf(file,"Random seed         = "OFD3"\n",randomseed);
    fprintf(file,"\n");
    fprintf(file,"Times for individual steps\n\n");
    fprintf(file,"Calculation of halo properties and initialisation of grid in r: "OFD1" seconds.\n",t1-t0);
    fprintf(file,"Initialisation of grid for distribution function: "OFD1" seconds.\n",t2-t1);
    fprintf(file,"Setting particle positions: "OFD1" seconds\n",t3-t2);
    fprintf(file,"Setting particle velocities: "OFD1" seconds\n",t4-t3);
    fprintf(file,"Setting remaining particle attributes: "OFD1" seconds\n",t5-t4);
    fprintf(file,"Calculating a few things and correct center of mass: "OFD1" seconds\n",t6-t5);
    t7 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(file,"Writing output: "OFD1" seconds\n",t7-t6);
    fprintf(file,"Total time: "OFD1" seconds\n",t7-t0);
    fclose(file);
   
    fprintf(stderr,"Done in "OFD1" seconds\nTotal time needed was "OFD1" seconds\n",t7-t6,t7-t0);

    free(bh);
    free(halo);
    free(gi);
    exit(0);

    } /* end of main function */


