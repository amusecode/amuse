/* 
** routines.c 
**
** Routines for HALOGEN
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc/malloc.h>
#include <assert.h>
#include <time.h>
#include "definitions.h"
#include "functions.h"
#include "routines.h"

/*
** Initialise fixed structures
*/

void initialise_fixed_structures(PARTICLE **bh, SI **halo, GI **gi) {
    *bh = malloc(sizeof(PARTICLE));
    assert(*bh != NULL);

    *halo = malloc(sizeof(SI));
    assert(*halo != NULL);
    (*halo)->sp = malloc(sizeof(SP));
    assert((*halo)->sp != NULL);
    (*halo)->griddf = malloc(sizeof(GRIDDF));
    assert((*halo)->griddf != NULL);

    *gi = malloc(sizeof(GI));
    assert(*gi != NULL);
    (*gi)->stuff = malloc(sizeof(STUFF));
    assert((*gi)->stuff != NULL);
    (*gi)->gridr = malloc(sizeof(GRIDR));
    assert((*gi)->gridr != NULL);
}

/*
** Set standard values for parameters
*/

void set_standard_values_for_parameters(INT *outputgridr, INT *outputgriddf, GI *gi, SI *halo) {
    *outputgridr = 0;
    *outputgriddf = 0;
    gi->dorvirexact = 0;
    sprintf(halo->systemname, "halo");
}

/*
** Routine for initialising systems
*/

void initialise_parameters(SI *si) {

    si->sp->alpha = -1;
    si->sp->beta = -1;
    si->sp->gamma = -1;
    si->sp->delta = -1;
    si->sp->M = 1;
    si->sp->rs = 1;
    si->sp->rhalf = -1;
    si->sp->rcutoff = -1;
    si->sp->rho0 = -1;
    si->sp->rvir = -1;
    si->sp->rdecay = -1;
    si->N0 = -1;
    si->rimp = SBI;
    si->r1 = -1;
    si->r100 = -1;
}

/*
** Routine for initialising black hole
*/

void initialise_black_hole(PARTICLE *p) {

    INT i;

    for (i = 0; i < 4; i++) {
        p->r[i] = 0;
        p->v[i] = 0;
    }
    p->mass = 0;
}

/*
** Routine for checking parameters of system
*/

void check_main_parameters(const SI *si) {

    if (si->sp->alpha == -1) {
        fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
        fprintf(stderr,"You have not set a value for alpha.\n");
        usage();
    }
    if (si->sp->beta == -1) {
        fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
        fprintf(stderr,"You have not set a value for beta.\n");
        usage();
    }
    if (si->sp->gamma == -1) {
        fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
        fprintf(stderr,"You have not set a value for gamma.\n");
        usage();
    }
    if (si->sp->gamma >= 3) {
        fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
        fprintf(stderr,"You have chosen gamma = "OFD1".\n",si->sp->gamma);
        fprintf(stderr,"This means your cumulative mass function is diverging at the centre.\n");
        fprintf(stderr,"Use a smaller value for gamma.\n");
        usage();
    }
    if (si->sp->M == -1) {        
        fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
        fprintf(stderr,"You have not set a value for the mass M.\n");
        usage();
    }
    if (si->N0 == -1) {
        fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
        fprintf(stderr,"You have not set a value for the number of paricles N0.\n");
        usage();
    }
}

/*
** Routine for calculating system parameters
*/

void calculate_parameters(const GI *gi, SI *si) {

    DOUBLE IM, IMcutoff;

    if (si->sp->beta > 3) {
        /*
        ** Finite mass models
        */
        if (si->sp->rs == -1) {
            fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
            fprintf(stderr,"For finite mass models you have to set a value for the scale radius rs.\n");
            usage();
        }
        if (si->sp->rcutoff != -1) {
            fprintf(stderr,"Warning for %s!\n",si->systemname);
            fprintf(stderr,"For finite mass models the cutoff radius rcutoff is not needed!\n");
            fprintf(stderr,"Hence, your input for the cutoff radius rcutoff (= "OFD1" kpc) was ignored.\n",si->sp->rcutoff);
        }
        si->sp->rdecay = 0;
        si->sp->delta = 0;
        IM = exp(lgamma((si->sp->beta-3)/si->sp->alpha));
        IM *= exp(lgamma((3-si->sp->gamma)/si->sp->alpha));
        IM /= (si->sp->alpha*exp(lgamma((si->sp->beta-si->sp->gamma)/si->sp->alpha)));
        si->sp->rho0 = si->sp->M/(4*M_PI*(si->sp->rs*si->sp->rs*si->sp->rs)*IM);
        si->sp->rcutoff = SBI;
    }
    else {
        /*
        ** Cutoff models
        */
        if ((si->sp->rs == -1) || (si->sp->rcutoff == -1)) {
            fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
            fprintf(stderr,"Specify values for the scale radius rs and cutoff radius rcutoff for models with cutoff.\n");
            usage();
        }
        si->sp->rdecay = CutoffFac*si->sp->rcutoff;
        si->sp->delta = si->sp->rcutoff/si->sp->rdecay + dlrhodlr(si->sp->rcutoff,si);
        IM = pow(1e-6,3-si->sp->gamma)/(3-si->sp->gamma); /* approximate inner integral */
        IM += integral(integrandIM,1e-6,si->sp->rcutoff/si->sp->rs,si);
        IMcutoff = 1/tau(si->sp->rcutoff,si);
        IMcutoff *= 1/(si->sp->rs*si->sp->rs*si->sp->rs);
        IMcutoff *= integral(integrandIMcutoff,si->sp->rcutoff,si->sp->rcutoff+2000*si->sp->rdecay,si);
        IM += IMcutoff;
        si->sp->rho0 = si->sp->M/(4*M_PI*(si->sp->rs*si->sp->rs*si->sp->rs)*IM);
    }
}

/* 
** Routine for initialising gridr 
*/

void initialise_gridr(GI *gi, PARTICLE *bh, SI *halo) {

    INT i;
    DOUBLE dlogr, logr, r, r3;
    DOUBLE Mencr, MencHalor, DeltaMencHalor;
    DOUBLE rhor, rhoHalor;
    DOUBLE rhoencr, rhoencHalor;
    DOUBLE Potr, Potoutr;
    GRIDR *gridr;
    SP *hp;

    gridr = gi->gridr;
    hp = halo->sp;
    dlogr = (log(gi->router)-log(gi->rinner))/(NGRIDR-1);
    i = 0;
    logr = log(gi->rinner);
    r = exp(logr);
    r3 = r*r*r;
    rhoHalor = rho(r,halo);
    rhor = rhoHalor;
    DeltaMencHalor = 4*M_PI*rhor*r3/(3-hp->gamma); /* approximate inner gridpoint by analytical calculation */
    MencHalor = DeltaMencHalor;
    Mencr = bh->mass + DeltaMencHalor;
    rhoencHalor = MencHalor/(4*M_PI*r3/3.0);
    rhoencr = Mencr/(4*M_PI*r3/3.0);
    gridr->r[i] = r;
    gridr->logr[i] = logr;
    halo->logr[i] = logr;
    gridr->rho[i] = rhor;
    gridr->logrho[i] = log(rhor);
    gridr->rhoHalo[i] = rhoHalor;
    gridr->logrhoHalo[i] = log(rhoHalor);
    gridr->rhoenc[i] = rhoencr;
    gridr->logrhoenc[i] = log(rhoencr);
    gridr->rhoencHalo[i] = rhoencHalor;
    gridr->logrhoencHalo[i] = log(rhoencHalor);
    halo->logrhoenc[i] = log(rhoencHalor);
    gridr->Menc[i] = Mencr;
    gridr->logMenc[i] = log(Mencr);
    gridr->MencHalo[i] = MencHalor;
    gridr->logMencHalo[i] = log(MencHalor);
    halo->logMenc[i] = log(MencHalor);
    gridr->eqrvcmax[i] = Mencr - 4*M_PI*rhor*r3;
    for (i = 1; i < NGRIDR; i++) {
        logr = log(gi->rinner) + i*dlogr;
        r = exp(logr);
        r3 = r*r*r;
        rhoHalor = rho(r,halo);
        rhor = rhoHalor;
        DeltaMencHalor =  integral(integrandMenc,gridr->r[i-1],r,halo);
        MencHalor = gridr->MencHalo[i-1] + DeltaMencHalor;
        Mencr = gridr->Menc[i-1] + DeltaMencHalor;
        rhoencHalor = MencHalor/(4*M_PI*r3/3.0);
        rhoencr = Mencr/(4*M_PI*r3/3.0);
        gridr->r[i] = r;
        gridr->logr[i] = logr;
        halo->logr[i] = logr;
        gridr->rho[i] = rhor;
        gridr->logrho[i] = log(rhor);
        gridr->rhoHalo[i] = rhoHalor;
        gridr->logrhoHalo[i] = log(rhoHalor);
        gridr->rhoenc[i] = rhoencr;
        gridr->logrhoenc[i] = log(rhoencr);
        gridr->rhoencHalo[i] = rhoencHalor;
        gridr->logrhoencHalo[i] = log(rhoencHalor);
        halo->logrhoenc[i] = log(rhoencHalor);
        gridr->Menc[i] = Mencr;
        gridr->logMenc[i] = log(Mencr);
        gridr->MencHalo[i] = MencHalor;
        gridr->logMencHalo[i] = log(MencHalor);
        halo->logMenc[i] = log(MencHalor);
        gridr->eqrvcmax[i] = Mencr - 4*M_PI*rhor*r3;
    }
    i = NGRIDR-1;
    Potoutr = 0; /* approximate outer gridpoint by 0 */
    Potr = (-1)*G*(gridr->Menc[i]/gridr->r[i]+Potoutr);
    gridr->Pot[i] = Potr;
    gridr->logPot[i] = log(-Potr);
    gridr->Potoutr[i] = Potoutr;
    for (i = (NGRIDR-2); i >= 0; i--) {
        Potoutr = integral(integrandPot,gridr->r[i],gridr->r[i+1],halo)+gridr->Potoutr[i+1];
        Potr = (-1)*G*(gridr->Menc[i]/gridr->r[i]+Potoutr);
        gridr->Pot[i] = Potr;
        gridr->logPot[i] = log(-Potr);
        gridr->Potoutr[i] = Potoutr;
    }
}

/*
** Routine for initialising griddf
*/

void initialise_griddf(const GI *gi, SI *si) {

    INT i, j, dj;
    DOUBLE fE;
    GRIDR *gridr;
    GRIDDF *griddf;

    gridr = gi->gridr;
    griddf = si->griddf;
    dj = (NGRIDR-1) / (NGRIDDF-1);
    for (i = (NGRIDDF-1), j = (NGRIDR-1); i >= 0; i--, j = j-dj) {
        fE = integraldf(j,gi,si);
        if (fE < 0) {
            fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
            fprintf(stderr,"You have chosen parameters that lead to a negative and hence unphysical distribution function.\n");
            usage();
            break;
        }
        griddf->r[i] = gridr->r[j];
        griddf->logr[i] = gridr->logr[j];
        griddf->E[i] = gridr->Pot[j];
        griddf->logE[i] = gridr->logPot[j];
        griddf->fE[i] = fE;
        griddf->logfE[i] = log(fE);
    }
}

/*
** Initialise all grids
*/

void initialise_all_grids(GI *gi, PARTICLE *bh, SI *halo, 
        INT outputgridr, INT outputgriddf, DOUBLE t0, DOUBLE *t1, 
        CHAR *FILENAME, CHAR *INPUTNAME) {
    
    FILE *file;
    
    /*
    ** Initialise gridr
    */

    gi->rinner = FACTORRINNER*halo->sp->rs;
    gi->router = halo->sp->rs;
    while (rho(halo->sp->rs,halo)/rho(gi->router,halo) < FACTORROUTER) {
        gi->router = gi->router*10;
    }

    initialise_gridr(gi,bh,halo);
    
    if (outputgridr == 1) {
        sprintf(FILENAME,"%s.gridr.dat",INPUTNAME);
        file = fopen(FILENAME,"w");
        assert(file != NULL);
        write_gridr(gi->gridr,file);
        fclose(file);
    }

    /*
    ** Initialise griddf
    */

    *t1 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nInitialising grid for distribution function... \n",*t1-t0);

    initialise_griddf(gi,halo);
    
    if (outputgriddf == 1) {
        sprintf(FILENAME,"%s.griddf.halo.dat",INPUTNAME);
        file = fopen(FILENAME,"w");
        assert(file != NULL);
        write_griddf(halo,file);
        fclose(file);
    }
}

/*
** Routine for initialising structure
*/

void initialise_structure(const GI *gi, SI *si) {

    si->N = si->N0/2;
    si->p = malloc(si->N*sizeof(PARTICLE));
    si->Mmin = exp(lininterpolate(NGRIDR,si->logr,si->logMenc,log(gi->rinner)));
    si->Mmax = exp(lininterpolate(NGRIDR,si->logr,si->logMenc,log(gi->router)));
    si->mass = si->sp->M/(2.0*si->N);
}

/*
** Routine for setting position of particles
*/

void set_positions(SI *si) {
    
    INT i, N;
    DOUBLE Mrand, logMrand, Mmin, Mmax;
    DOUBLE rrand, logrrand;
    DOUBLE theta, phi;
    PARTICLE *p;

    N = si->N;
    p = si->p;
    Mmin = si->Mmin;
    Mmax = si->Mmax;
    for (i = 0; i < N; i++) {
        Mrand = Mmin + rand01()*(Mmax - Mmin);
        logMrand = log(Mrand);
        logrrand = lininterpolate(NGRIDR,si->logMenc,si->logr,logMrand);
        rrand = exp(logrrand);
        theta = acos(2.0*rand01() - 1.0);
        phi = rand01()*2.0*M_PI;
        p[i].r[0] = rrand;
        p[i].r[1] = rrand*sin(theta)*cos(phi);
        p[i].r[2] = rrand*sin(theta)*sin(phi);
        p[i].r[3] = rrand*cos(theta);
    }
}

/*
** Routine for setting velocities of particles 
*/

void set_velocities(const GI *gi, SI *si) {
    
    INT i, N;
    DOUBLE r, Erand, Potr;
    DOUBLE fEmax, fErand, fEcheck;
    DOUBLE vesc, vrand;
    DOUBLE theta, phi;
    GRIDR *gridr;
    GRIDDF *griddf;
    PARTICLE *p;

    gridr = gi->gridr;
    griddf = si->griddf;
    N = si->N;
    p = si->p;
    for (i = 0; i < N; i++) {
        r = p[i].r[0];
        Potr = Pot(r,gi);
        vesc = vescape(r,gi);
        fEmax = f2(r,si);
        vrand = 0;
        Erand = 0;
        fErand = 0;
        fEcheck = 1;
        while (fEcheck > fErand) {
            vrand = pow(rand01(),1.0/3.0)*vesc;
            Erand = 0.5*vrand*vrand + Potr;
            fErand = f1(Erand,si);
            fEcheck = rand01()*fEmax;
        }
        theta = acos(2.0*rand01() - 1.0);
        phi = rand01()*2.0*M_PI;
        p[i].v[0] = vrand;
        p[i].v[1] = vrand*sin(theta)*cos(phi);
        p[i].v[2] = vrand*sin(theta)*sin(phi);
        p[i].v[3] = vrand*cos(theta);
    }
}

/*
** Routine for setting the remaining attributes
*/

void set_attributes(const GI *gi, SI *si) {

    INT i, N;
    DOUBLE mass;
    GRIDR *gridr;
    PARTICLE *p;
    
    gridr = gi->gridr;
    N = si->N;
    p = si->p;
    mass = si->mass;
    for (i = 0; i < N; i++) {
        p[i].mass = mass;
    }
}

/*
** Routine for doubling particles with mirror halo
*/

void double_particles(SI *si) {

    INT i, N;
    PARTICLE *p;

    N = si->N;
    p = si->p;
    p = realloc(p,2*N*sizeof(PARTICLE));
    if (N > 0) {
        assert(p != NULL);
    }
    for (i = 0; i < N; i++) {
        p[N+i].r[0] = p[i].r[0];
        p[N+i].r[1] = -p[i].r[1];
        p[N+i].r[2] = -p[i].r[2];
        p[N+i].r[3] = -p[i].r[3];
        p[N+i].v[0] = p[i].v[0];
        p[N+i].v[1] = -p[i].v[1];
        p[N+i].v[2] = -p[i].v[2];
        p[N+i].v[3] = -p[i].v[3];
        p[N+i].mass = p[i].mass; 
    }
    si->N = 2*N;
    si->p = p;
}

/*
** Routine for calculating center of mass position and velocity, and angular momentum
*/

void calculate_stuff(GI *gi, PARTICLE *bh, SI *halo) {

    INT i, j, N;
    DOUBLE mass, x, y, z, sumrvir;
    STUFF *stuff;
    PARTICLE *p;
    
    sumrvir = 0;
    stuff = gi->stuff;
    stuff->N = 0;
    stuff->Mp = 0;
    stuff->Ekin = 0;
    stuff->Epot = 0;
    for(i = 0; i < 4; i++) {
        stuff->Cr[i] = 0;
        stuff->Cv[i] = 0;
    }
    N = halo->N;
    p = halo->p;
    stuff->N += halo->N;
    if (bh->mass > 0) {
        stuff->N++;
        stuff->Mp += bh->mass;
    }
    for (i = 0; i < N; i++) {
        if(p[i].r[0] < halo->rimp) {
            halo->rimp = p[i].r[0];
        }
        mass = p[i].mass;
        stuff->Mp += mass;
        stuff->Cr[1] += mass*p[i].r[1];
        stuff->Cr[2] += mass*p[i].r[2];
        stuff->Cr[3] += mass*p[i].r[3];
        stuff->Cv[1] += mass*p[i].v[1];
        stuff->Cv[2] += mass*p[i].v[2];
        stuff->Cv[3] += mass*p[i].v[3];
        stuff->Ekin += mass*p[i].v[0]*p[i].v[0]/2.0;
        stuff->Epot += mass*Pot(p[i].r[0],gi);
        if (gi->dorvirexact == 1) {
            for (j = i+1; j < N; j++) {
                x = p[i].r[1] - p[j].r[1];
                y = p[i].r[2] - p[j].r[2];
                z = p[i].r[3] - p[j].r[3];
                sumrvir += p[i].mass*p[j].mass/sqrt(x*x + y*y +z*z);
            }
        }
    }        
    for(i = 1; i < 4; i++) {
        stuff->Cr[i] = stuff->Cr[i]/stuff->Mp;
        stuff->Cv[i] = stuff->Cv[i]/stuff->Mp;
    }
    stuff->Cr[0] = sqrt(stuff->Cr[1]*stuff->Cr[1]+stuff->Cr[2]*stuff->Cr[2]+stuff->Cr[3]*stuff->Cr[3]);
    stuff->Cv[0] = sqrt(stuff->Cv[1]*stuff->Cv[1]+stuff->Cv[2]*stuff->Cv[2]+stuff->Cv[3]*stuff->Cv[3]);
    stuff->Epot = stuff->Epot/2.0;
    stuff->Etot = stuff->Ekin + stuff->Epot;
    halo->sp->rhalf = exp(lininterpolate(NGRIDR,gi->gridr->Menc,gi->gridr->logr,halo->sp->M/2.0));
    halo->r1 = pow(((3.0-halo->sp->gamma)*halo->mass)/(4*M_PI*halo->sp->rho0*pow(halo->sp->rs,halo->sp->gamma)),1.0/(3.0-halo->sp->gamma));
    halo->r100 = pow(((3.0-halo->sp->gamma)*100*halo->mass)/(4*M_PI*halo->sp->rho0*pow(halo->sp->rs,halo->sp->gamma)),1.0/(3.0-halo->sp->gamma));
    if (gi->dorvirexact == 1) {
        halo->sp->rvir = (stuff->Mp*stuff->Mp)/(2.0*sumrvir);
    }
    else {
        sumrvir = (-1)*stuff->Epot/G;
        halo->sp->rvir = (stuff->Mp*stuff->Mp)/(2.0*sumrvir);
    }
}

/*
** Routines for writing out grids
*/

void write_gridr(GRIDR *gridr, FILE *file) {

    INT i;

    for (i = 0; i < NGRIDR; i++) {
        fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
                gridr->r[i],gridr->logr[i],gridr->rho[i],gridr->logrho[i],gridr->rhoenc[i],gridr->logrhoenc[i],
                gridr->Menc[i],gridr->logMenc[i],gridr->Pot[i],gridr->logPot[i]);
    }
}

void write_griddf(SI *si, FILE *file) {

    INT i;
    GRIDDF *griddf;

    griddf = si->griddf;
    for (i = 0; i < NGRIDDF; i++) {
        fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
                griddf->r[i],griddf->logr[i],griddf->E[i],griddf->logE[i],griddf->fE[i],griddf->logfE[i]);
    }
}

/*
** Usage description
*/

#ifdef AMUSE_FLAG
extern int commit_parameters_result;

void usage() {
    commit_parameters_result = -2;
}
#else
void usage() {

    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments.\n\n");
    fprintf(stderr,"-a <value>          : alpha parameter in density profile\n");
    fprintf(stderr,"-b <value>          : beta parameter in density profile\n");
    fprintf(stderr,"-c <value>          : gamma parameter in density profile\n");
    fprintf(stderr,"-rs <value>         : scale radius (default rs = 1)\n");
    fprintf(stderr,"-rcutoff <value>    : cutoff radius for cutoff models\n");
    fprintf(stderr,"-N <value>          : total number of particles\n");
    fprintf(stderr,"-M <value>          : total mass (default M = 1)\n");
    fprintf(stderr,"-MBH <value>        : mass of black hole (default MBH = 0)\n");
    fprintf(stderr,"-name <value>       : name of the output file\n");
    fprintf(stderr,"-ogr                : set this flag for outputting grid in r\n");
    fprintf(stderr,"-ogdf               : set this flag for outputting grid for distribution function\n");
    fprintf(stderr,"-randomseed <value> : set this flag for setting a value for a random seed (default: random value)\n");
    fprintf(stderr,"-dorvirexact        : set this flag for calculating rvir exactly via N^2 sum - Warning: time consuming for large N!\n");
    fprintf(stderr,"\n");
    exit(1);
}
#endif
