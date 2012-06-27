/* 
** functions.c
**
** Functions for HALOGEN4MUSE
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "definitions.h"
#include "functions.h"


/*
** Functions for integration via Simpson's rule
** Iteration until precision (TOL) is reached 
*/

DOUBLE integral(DOUBLE (*function)(DOUBLE,const SI (*)), DOUBLE a, DOUBLE b, const SI *si) {

    INT i;
    DOUBLE T[3], SA, SB;

    T[0] = trapez(function,a,b,0,si);
    T[1] = trapez(function,a,b,1,si);
    T[2] = trapez(function,a,b,2,si);
    SA = (4.0*T[1]-T[0])/3.0;
    SB = (4.0*T[2]-T[1])/3.0;
    for(i = 3; i < (NINTMAX+1); i++) {
	T[0] = T[1];
	T[1] = T[2];
	T[2] = trapez(function,a,b,i,si);
	SA = SB;
	SB = (4.0*T[2]-T[1])/3.0;
	if(i > NINTMIN) {
	    if((fabs(SB-SA) < TOL*fabs(SB)) || (SA == 0.0 && SB == 0.0)) {
		return SB;
		}
	    }
	}
    fprintf(stderr,"Warning!\n");
    fprintf(stderr,"Too many steps in function integral!\n");
    fprintf(stderr,"a = "OFD3"\n",a);
    fprintf(stderr,"b = "OFD3"\n",b);
    fprintf(stderr,"sum = "OFD3"\n",SB);
    fprintf(stderr,"Sum not converged within tolerance of %e\n",TOL);
    fprintf(stderr,"Abort tolerance was %e\n",fabs((SA-SB)/SB));
    return SB;
    }

DOUBLE trapez(DOUBLE (*function)(DOUBLE, const SI(*)), DOUBLE a, DOUBLE b, INT n, const SI *si) {

    INT i, j, N[n+1];
    DOUBLE deltax, sum, sumN;

    if (n == 0) {
	return (0.5*(b-a)*((*function)(a,si)+(*function)(b,si)));
	}
    else {
	for (i = 0; i < (n+1); i++) {
	    N[i] = pow(2,i);
	    }
	sum = 0.5*((*function)(a,si)+(*function)(b,si));
	for(i = 1; i < (n+1); i++) {
	    deltax = (b-a)/N[i];
	    sumN = 0;
	    for (j = 1; j < (N[i-1]+1); j++) {
		sumN += (*function)(a+(2*j-1)*deltax,si);
		}
	    sum = sum + sumN;
	    }
	return (sum*(b-a)/N[n]);
	}
    }

DOUBLE integraldf(INT j, const GI *gi, const SI *si) {
    
    INT i;
    DOUBLE T[3], SA, SB;

    T[0] = trapezdf(j,0,gi,si);
    T[1] = trapezdf(j,1,gi,si);
    T[2] = trapezdf(j,2,gi,si);
    SA = (4.0*T[1]-T[0])/3.0;
    SB = (4.0*T[2]-T[1])/3.0;
    for(i = 3; i < (NINTMAXDF+1); i++) {
	T[0] = T[1];
	T[1] = T[2];
	T[2] = trapezdf(j,i,gi,si);
	SA = SB;
	SB = (4.0*T[2]-T[1])/3.0;
	if(i > NINTMINDF) {
	    if((fabs(SB-SA) < TOLDF*fabs(SB)) || (SA == 0.0 && SB == 0.0)) {
		return SB;
		}
	    }
	}
    fprintf(stderr,"Warning!\n");
    fprintf(stderr,"Too many steps in function integraldf!\n");
    fprintf(stderr,"E = "OFD3" kpc^2 Gyr^-2\n",gi->gridr->Pot[j]);
    fprintf(stderr,"r = "OFD3" kpc\n",gi->gridr->r[j]);
    fprintf(stderr,"f(E) = "OFD3" MU Gyr^3 kpc^-6\n",SB);
    fprintf(stderr,"f(E) not converged within tolerance of %e\n",TOLDF);
    fprintf(stderr,"Abort tolerance was %e\n",fabs((SA-SB)/SB));
    if (SB != SB) {
	fprintf(stderr,"f(E) value set to %e MU Gyr^3 kpc^-6\n",DFFAILUREMAX);
	return DFFAILUREMAX;
	}
    else {
	return SB;
	}
    }

DOUBLE trapezdf(INT k, INT n, const GI *gi, const SI *si) {

    INT i, j, N[n+1];
    DOUBLE deltax, sum, sumN;
    DOUBLE x, xlower, xupper;
    DOUBLE r, rlower, rupper;
    DOUBLE Phi, Philower, Phiupper;
    GRIDR *gridr;

    gridr = gi->gridr;
    Philower = gridr->Pot[NGRIDR-1];
    Phiupper = gridr->Pot[k];
    xlower = asin(sqrt(Philower/Phiupper));
    xupper = M_PI/2;
    rlower = gridr->r[NGRIDR-1];
    rupper = gridr->r[k];
    if (n == 0) {
	return ((sqrt(-Phiupper/2)/(M_PI*M_PI))*0.5*(xupper-xlower)*(d2rhodPhi2(rlower,gi,si)*sin(xlower)+d2rhodPhi2(rupper,gi,si)*sin(xupper)));
	}
    else {
	for (i = 0; i < (n+1); i++) {
	    N[i] = pow(2,i);
	    }
	sum = 0.5*(d2rhodPhi2(rlower,gi,si)*sin(xlower)+d2rhodPhi2(rupper,gi,si)*sin(xupper));
	for(i = 1; i < (n+1); i++) {
	    deltax = (xupper-xlower)/N[i];
	    sumN = 0;
	    for (j = 1; j < (N[i-1]+1); j++) {
		x = xlower+(2*j-1)*deltax;
		Phi =  Phiupper*sin(x)*sin(x);
		r = exp(lininterpolate(NGRIDR,gridr->logPot,gridr->logr,log(-Phi))); 
		sumN += d2rhodPhi2(r,gi,si)*sin(x);
		}
	    sum = sum + sumN;
	    }
	return ((sqrt(-Phiupper/2)/(M_PI*M_PI))*sum*(xupper-xlower)/N[n]);
	}
    }

/*
** Function for locating x on a monotonic increasing or decreasing grid 
*/

INT locate(INT n, const DOUBLE *grid, DOUBLE x) {

    INT jl, jm, ju, ascend;

    jl = -1;
    ju = n;
    ascend = (grid[n-1] >= grid[0]);
    while ((ju-jl) > 1) {
	jm = (ju+jl)/2;
	if ((x >= grid[jm]) == ascend) {
	    jl = jm;
	    }
	else {
	    ju = jm;
	    }
	}
    if (x == grid[0]) {
	return 0;
	}
    else if (x == grid[n-1]) {
	return (n-2);
	}
    else {
	return jl;
	}
    }

/*
** Function for linear interpolation between two grid points 
*/

DOUBLE lininterpolate(INT n, const DOUBLE *gridx, const DOUBLE *gridy, DOUBLE x) {

    INT i; 
    DOUBLE dgx, dgy, dx, m, y, RE;

    i = locate(n,gridx,x);
    if (i < 0) {
	fprintf(stderr,"Warning!\n");
	fprintf(stderr,"x = "OFD3" was below range of array! Index i = "OFI1" / n = "OFI1".\n",x,i,n);
	fprintf(stderr,"Array between "OFD3" and "OFD3".\n",gridx[0],gridx[n-1]);
	RE = fabs((x-gridx[0])/gridx[0]);
	if (RE < TOLLININT) {
	    fprintf(stderr,"Relative error (= "OFD1") at lower boundary was within tolerance of %e\n",RE,TOLLININT);
	    fprintf(stderr,"Index set to i = 0\n");
	    i = 0;
	    }
	}
    else if (i > n-2) {
	fprintf(stderr,"Warning!\n");
	fprintf(stderr,"x = "OFD3" was above range of array! Index i = "OFI1" / n = "OFI1".\n",x,i,n);
	fprintf(stderr,"Array between "OFD3" and "OFD3".\n",gridx[0],gridx[n-1]);
	RE = fabs((x-gridx[n-1])/gridx[n-1]);
	if (RE < TOLLININT) {
	    fprintf(stderr,"Relative error (= "OFD1") at upper boundary was within tolerance of %e\n",RE,TOLLININT);
	    fprintf(stderr,"Index set to i = n-2 = "OFI1"\n",n-2);
	    i = n-2;
	    }
	}
    dgy = gridy[i+1] - gridy[i];
    dgx = gridx[i+1] - gridx[i];
    m = dgy/dgx;
    dx = x - gridx[i];
    y = gridy[i] + dx*m;

    return y;
    }

/*
** Function for generating random numbers between [0,1] 
*/

DOUBLE rand01() {

    return drand48();
//    return ( ((DOUBLE) rand()) / ((DOUBLE) RAND_MAX) );
    }

/* 
** alpha-beta-gamma density function with exponential cutoff 
** except for finite mass models 
*/

DOUBLE rho(DOUBLE r, const SI *si) {

    DOUBLE fac1, fac2;
    SP *sp;

    sp = si->sp;
    if (sp->beta > 3) {
	/*
	** Finite mass models
	*/
	return (sp->rho0/tau(r,si));
	}
    else {
	/*
	** Cutoff models
	*/
	if (r <= sp->rcutoff) {
	    return (sp->rho0/tau(r,si));
	    }
	else {
	    fac1 = pow((r/sp->rcutoff),sp->delta);
	    fac2 = exp(-(r-sp->rcutoff)/sp->rdecay);
	    return (sp->rho0/tau(sp->rcutoff,si)*fac1*fac2);
	    }
	}
    }

/* 
** Derivative of density drho/dr 
*/

DOUBLE drhodr(DOUBLE r, const SI *si) {

    SP *sp;

    sp = si->sp;
    if (sp->beta > 3) {
	/*
	** Finite mass models
	*/
	return (-rho(r,si)*eta(r,si));
	}
    else {
	/*
	** Cutoff models
	*/
	if (r <= sp->rcutoff) {
	    return (-rho(r,si)*eta(r,si));
	    }
	else {
	    return (rho(r,si)*(sp->delta/r-1/sp->rdecay));
	    }
	}
    }

/* 
** Derivative of density d^2rho/dr^2 
*/

DOUBLE d2rhodr2(DOUBLE r, const SI *si) {

    SP *sp;

    sp = si->sp;
    if (sp->beta > 3) {
	/*
	** Finite mass models
	*/
	return (rho(r,si)*(pow(eta(r,si),2)-detadr(r,si)));
	}
    else {
	/*
	** Cutoff models
	*/
	if (r <= sp->rcutoff) {
	    return (rho(r,si)*(pow(eta(r,si),2)-detadr(r,si)));
	    }
	else {
	    return (rho(r,si)*(pow((sp->delta/r-1/sp->rdecay),2)-sp->delta/(r*r)));
	    }
	}
    }

/* 
** Derivative of density d^2rho/dPhi^2 
*/

DOUBLE d2rhodPhi2(DOUBLE r, const GI *gi, const SI *si) {

    DOUBLE Mencr;
    DOUBLE fac1, fac2;
    SP *sp;

    sp = si->sp;
    Mencr = Menc(r,gi);
    fac1 = r*r/(G*G*Mencr*Mencr);
    fac2 = 2*r-4*M_PI*r*r*r*r*rho(r,si)/Mencr;
    return (fac1*(r*r*d2rhodr2(r,si)+fac2*drhodr(r,si)));
    }

/* 
** Derivative of density dlrho/dlr 
*/

DOUBLE dlrhodlr(DOUBLE r, const SI *si) {

    return (r/rho(r,si)*drhodr(r,si));
    }

/*
** Auxiliary functions 
*/

DOUBLE eta(DOUBLE r, const SI *si) {

    DOUBLE fac1, fac2, fac3;
    SP *sp;

    sp = si->sp;
    fac1 = (sp->beta-sp->gamma)/sp->rs;
    fac2 = pow((r/sp->rs),(sp->alpha-1));
    fac3 = 1+pow((r/sp->rs),sp->alpha);
    return ((sp->gamma/r)+(fac1*fac2/fac3));
    }

DOUBLE detadr(DOUBLE r, const SI *si) {

    DOUBLE fac1, fac2, fac3, fac4;
    SP *sp;

    sp = si->sp;
    fac1 = (sp->beta-sp->gamma)/(sp->rs*sp->rs);
    fac2 = pow((r/sp->rs),(sp->alpha-2));
    fac3 = (sp->alpha-1)*(1+pow((r/sp->rs),sp->alpha))-sp->alpha*pow((r/sp->rs),sp->alpha);
    fac4 = pow((1+pow((r/sp->rs),sp->alpha)),2);
    return (-sp->gamma/(r*r)+fac1*fac2*fac3/fac4);
    }

DOUBLE tau(DOUBLE r, const SI *si) {
    
    DOUBLE exp1, exp2, exp3;
    DOUBLE fac1, fac2, fac3;
    SP *sp;

    sp = si->sp;
    exp1 = sp->gamma;
    exp2 = sp->alpha;
    exp3 = (sp->beta-sp->gamma)/sp->alpha;
    fac1 = pow(r/sp->rs,exp1);
    fac2 = 1+pow(r/sp->rs,exp2);
    fac3 = pow(fac2,exp3);
    return (fac1*fac3);
    }

/*
** Integrand of IM integral 
*/

DOUBLE integrandIM(DOUBLE r, const SI *si) {
    
    DOUBLE exp1, exp2, exp3;
    DOUBLE fac1, fac2, fac3;
    SP *sp;

    sp = si->sp;
    exp1 = 2-sp->gamma;
    exp2 = sp->alpha;
    exp3 = (sp->beta-sp->gamma)/sp->alpha;
    fac1 = pow(r,exp1);
    fac2 = 1+pow(r,exp2);
    fac3 = pow(fac2,exp3);
    return (fac1/fac3);
    }

/*
** Integrand of IMcutoff integral
*/

DOUBLE integrandIMcutoff(DOUBLE r, const SI *si) {

    DOUBLE fac1, fac2, fac3;
    SP *sp;

    sp = si->sp;
    fac1 = r*r;
    fac2 = pow(r/sp->rcutoff,sp->delta);
    fac3 = exp(-(r-sp->rcutoff)/sp->rdecay);
    return (fac1*fac2*fac3);
    }

/*
** Integrand of enclosed mass integral 
*/

DOUBLE integrandMenc(DOUBLE r, const SI *si) {

    return (4*M_PI*r*r*(rho(r,si)));
    }

/*
** Integrand of outer potential integral 
*/

DOUBLE integrandPot(DOUBLE r, const SI *si) {

    return (4*M_PI*r*(rho(r,si)));
    }

/* 
** Function for calculating total enclosed mass within radius r 
*/

DOUBLE Menc(DOUBLE r, const GI *gi) {

    return (exp(lininterpolate(NGRIDR,gi->gridr->logr,gi->gridr->logMenc,log(r))));
    }

/*
** Function for calculating potential at radius r 
*/

DOUBLE Pot(DOUBLE r, const GI *gi) {

    return (-exp(lininterpolate(NGRIDR,gi->gridr->logr,gi->gridr->logPot,log(r))));
    }

/*
** Function for calculating the escape velocity at position r 
*/

DOUBLE vescape(DOUBLE r, const GI *gi) {

    return (sqrt(2.0*fabs(Pot(r,gi)-gi->gridr->Pot[NGRIDR-1])));
    }

/*
** Function for calculating the dynamical time at position r 
*/

DOUBLE Tdyn(DOUBLE r, const GI *gi) {

    return (2*M_PI*sqrt((r*r*r)/(G*Menc(r,gi))));
    }

/*
** Function for calculating the value of the distribution function at energy E 
*/

DOUBLE f1(DOUBLE E, const SI *si) {

    return (exp(lininterpolate(NGRIDDF,si->griddf->logE,si->griddf->logfE,log(-E))));
    }

/*
** Function for calculating the value of the distribution function at position r 
*/

DOUBLE f2(DOUBLE r, const SI *si) {

    return (exp(lininterpolate(NGRIDDF,si->griddf->logr,si->griddf->logfE,log(r))));
    }
