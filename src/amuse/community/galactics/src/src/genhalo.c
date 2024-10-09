#include "main.h"

// random skip factor, must be larger than the number of random draws per particle
#define SKIP  100000

int main(argc,argv)
int argc;
char **argv;
{
	int i, j, k, nobj=10000;
  long long lseed=0;
	int seed= -123;
	int nrtrunc;
	float rtrunc;
	float rhoran, massapp;
	float x, y, z, vx, vy, v, v2, R, vmax, vmax2;
	float phi, cph, sph, cth, sth, vR, vp, vz;
	float E, rad, rad2;
	float f0, frand, fmax, psi;
	float fr, fz;
	float dfhalo_(), halodens_(), pot_();
	float ran1();
  void ran_seed();
	float dr, rhomax1;
	float t, mass;
	float u1, v1, u1max, v1max;
	float xcm, ycm, zcm, vxcm, vycm, vzcm;
	float zip = 0.0;
	float rhomax, rhomin, rhotst;
	float rhocur, rcyl;
	float stream=0.5;
	float psid,psic,nnn;
	float dfmax;
	float b, c, d, int_const, sig2, aeff;
	float v02, vb2, coef1, coef2;
	int icofm=1;
	int jj;
	char harmfile[80];
	float fcut_halo, fcut_bulge;
	float kinetic, potential;

	strcpy(harmfile,"dbh.dat");
	fquery("Enter streaming fraction",&stream);
	iquery("Enter the number of particles",&nobj);
	llquery("Enter negative integer seed",&lseed);
	iquery("Center particles (0=no, 1=yes)",&icofm);
	readmassrad(); /* reads in mass and truncation radius of components */
	readharmfile_(harmfile,&gparam,&cparam,&bparam);
	readdenspsihalo_();

	chalo = gparam.chalo; v0 = gparam.v0; a = gparam.a; 
	psi0 = gparam.psi0;

	psic = cparam.psic; psid = cparam.psid;

	//	fcut_halo = bparam.fcut_halo; fcut_bulge = bparam.fcut_bulge;

	fcut_halo = dfhalo_(&psic);
	//	fcut_halo = 0.;

	rtrunc = haloedge;
	mass = halomass/nobj;
	sig2 = psi0;
	u1max = rtrunc;
	v1max = 0.5*M_PI;

	r = (phase *) calloc(nobj,sizeof(phase));

/* `Find maximum of rho*r^2 */

	nrtrunc = 10000;
	dr = rtrunc/nrtrunc;
	rhomax1 = 0.0;
  	for(i=1; i<nrtrunc; i++) {
  	  float z, rcyl, rhocur; 
  	  rcyl = i*dr; 
  	  z = 0; 
  	  rhocur = halodens_(&rcyl,&z); 
  	  rhocur *= (rcyl*rcyl); 
  	  if( rhocur > rhomax1  )
  	    rhomax1 = rhocur; 
  	} 
  	j = 0; 
  	for(i=1; i<nrtrunc; i++) { 
  	  rcyl = 0; 
  	  z = i*dr; 
  	  rhocur = halodens_(&rcyl,&z);
  	  rhocur *= (z*z); 
  	  if( rhocur > rhomax1  )
  	    rhomax1 = rhocur; 
  	} 
  	fprintf(stderr,"rhomax1 = %g\n",rhomax1); 

  	rhomax = 1.5*rhomax1; 
  	rhomin = 1.0e-10*rhomax; 

	kinetic = potential = 0.;

  	fprintf(stderr,"Calculating halo positions and velocities\n");
	if(lseed<0) lseed=-lseed;
  ran_seed(SKIP*lseed);
  for(i=0; i<nobj;) {
restart:
	  u1 = u1max*ran1(&seed);
	  v1 = 2.0*v1max*(ran1(&seed) - 0.5);
	  R = u1;
	  z = R*tan(v1);
	  
	  rhotst = halodens_(&R,&z);
	  rhotst = rhotst*(R*R + z*z);
	  j++;
	  if( rhotst < rhomin )
	    continue;
	  
	  rhoran = (rhomax - rhomin)*ran1(&seed);
	  if( rhoran > rhotst )
	    continue;
	  phi = 2.0*M_PI*ran1(&seed);
	  x = R*cos(phi);
	  y = R*sin(phi);
	  psi = pot_(&R,&z);
	  if (psi < psic)
	    continue;
       	  vmax2 = 2*(psi-psic);
      	  vmax = sqrt(vmax2);
	  f0 = 0.0; frand = 1.0; /* dummy starters */
	  j = 0;
	  fmax = dfhalo_(&psi) - fcut_halo;
	  while( frand > f0 ) {
	    v2 = 1.1*vmax2;
	    while( v2 > vmax2 ) {
	      vR = 2*vmax*(ran1(&seed) - 0.5);
	      vp = 2*vmax*(ran1(&seed) - 0.5);
	      vz = 2*vmax*(ran1(&seed) - 0.5);
	      v2 = vR*vR + vp*vp + vz*vz;
	      E = psi - 0.5*v2;
	    }
	    f0 = dfhalo_(&E)-fcut_halo;
	    frand = fmax*ran1(&seed);
		if( j>10000 ) {
			fprintf(stderr,"j>10000 - E=%g psi=%g frand=%g f0=%g\n",
				E,psi,frand,f0);
			goto restart;
		}
	    j++;
	  }
	  
	  if( ran1(&seed) < stream )
	    vp = fabs(vp);
	  else
	    vp = -fabs(vp);
	  
	  cph = x/R; sph = y/R;
	  vx = vR*cph - vp*sph;
	  vy = vR*sph + vp*cph;
	  
	  r[i].x = (float) x;
	  r[i].y = (float) y;
	  r[i].z = (float) z;
	  r[i].vx = (float)vx;
	  r[i].vy = (float)vy;
	  r[i].vz = (float)vz;
	  i++;ran_seed(SKIP*(lseed + i));
	  if( i % 1000 == 0 ) fprintf(stderr,".");
	}
	fprintf(stderr,"\n");
	
	if( icofm ) {
	  xcm = ycm =zcm = vxcm =vycm =vzcm = 0;
	  for(i=0; i<nobj; i++) {
	    xcm += r[i].x;
	    ycm += r[i].y;
	    zcm += r[i].z;
	    vxcm += r[i].vx;
	    vycm += r[i].vy;
	    vzcm += r[i].vz;
	  }
	  xcm /= nobj; ycm /=nobj; zcm /= nobj;
	  vxcm /= nobj; vycm /=nobj; vzcm /= nobj;
	  for(i=0; i<nobj; i++) {
	    r[i].x -= xcm;
	    r[i].y -= ycm;
	    r[i].z -= zcm;
	    r[i].vx -= vxcm;
	    r[i].vy -= vycm;
	    r[i].vz -= vzcm;
	  }
	}
	t = 0.0;
#ifdef ASCII
	fprintf(stdout,"%d %f\n",nobj,t);
	for(i=0; i<nobj; i++) {
	  fprintf(stdout,"% 15.7e % 15.7e % 15.7e % 15.7e % 15.7e % 15.7e % 15.7e\n", mass, r[i].x, r[i].y, r[i].z, r[i].vx, r[i].vy, r[i].vz);
	}
#else
	for(i=0; i<nobj; i++) {
	  fwrite(&mass,sizeof(float),1,stdout);
      fwrite(r+i,sizeof(phase),1,stdout);
    }
#endif
	exit(0);
}

