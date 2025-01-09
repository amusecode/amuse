#include "main.h"
#include "query.h"
#include "readdenspsi.h"
#include "readharmfile.h"
#include "readmassrad.h"

// random skip factor, must be larger than the number of random draws per particle
#define SKIP  100000

int main(int argc, char **argv)
{
	int i, j, k, nobj=10000;
  long long lseed=0;
	int seed= -123;
	int nrtrunc;
	float rtrunc, rtrunc1;
	float rhoran, massapp;
	float x, y, z, vx, vy, v, v2, R, vmax, vmax2;
	float phi, cph, sph, cth, sth, vR, vp, vz;
	float E, rad, rad2;
	float f0, frand, fmax, psi;
	float dfbulge_(), bulgedens_(), pot_();
	float ran1();
  void ran_seed();
	float dr, rhomax1;
	float t, mass;
	float u1, v1, u1max, v1max;
	float xcm, ycm, zcm, vxcm, vycm, vzcm;
	float zip = 0.0;
	float rhomax, rhomin, rhotst;
	float rhocur, rcyl;
	float fr, fz;
	float stream=0.5;
	float psid,psic,nnn;
	float dfmax;
	float b, c, d, int_const, sig2, aeff;
	float v02, vb2, coef1, coef2;
	int icofm=1;
	int jj;
	FILE *rhoout;
	char harmfile[80];
	float fcut_bulge, fcut_halo;
	float kinetic, potential;

	strcpy(harmfile,"dbh.dat");
	fquery("Enter streaming fraction",&stream);
	iquery("Enter the number of particles",&nobj);
	llquery("Enter negative integer seed",&lseed);
	iquery("Center particles (0=no, 1=yes)",&icofm);
	readmassrad(); /* reads in mass and truncation radius of components */
	readharmfile_(harmfile,&gparam,&cparam,&bparam);
	readdenspsibulge_();

	chalo = gparam.chalo; v0 = gparam.v0; a = gparam.a; 
	nnn = gparam.nnn; v0bulge=gparam.v0bulge;
	psi0 = gparam.psi0; haloconst = gparam.haloconst;

	psic = cparam.psic; psid = cparam.psid;

	fcut_bulge = dfbulge_(&psic);
	fcut_bulge = 0.;

	rtrunc = bulgeedge;
	mass = bulgemass/nobj;
	sig2 = psi0;
	u1max = rtrunc;
	v1max = 0.5*M_PI;

	r = (phase *) calloc(nobj,sizeof(phase));

/* `Find maximum of rho*r^2 */
	rhoout = fopen("rhobulge.dat","w");

	nrtrunc = 10000;
	dr = rtrunc/nrtrunc;
	rhomax1 = 0.0;
  	for(i=1; i<nrtrunc; i++) {
  	  float z, rcyl, rhocur; 
  	  rcyl = i*dr; 
  	  z = 0; 
  	  rhocur = bulgedens_(&rcyl,&z); 
  	  rhocur *= (rcyl*rcyl); 
	  fprintf(rhoout,"%g %g \n",log10(rcyl),log10(rhocur));
  	  if( rhocur > rhomax1  )
  	    rhomax1 = rhocur; 
  	} 
  	fprintf(stderr,"rhomax on R-axis is %g\n",rhomax1); 
  	j = 0; 
  	for(i=1; i<nrtrunc; i++) { 
  	  rcyl = 0; 
  	  z = i*dr; 
  	  rhocur = bulgedens_(&rcyl,&z);
  	  rhocur *= (z*z); 
  	  if( rhocur > rhomax1  )
  	    rhomax1 = rhocur; 
  	} 
  	fprintf(stderr,"rhomax1 = %g\n",rhomax1); 
  	rhomax = 1.2*rhomax1; 
  	rhomin = 1.0e-10*rhomax; 

	kinetic = potential = 0.;

  	fprintf(stderr,"Calculating bulge positions and velocities\n");
	if(lseed<0) lseed=-lseed;
  ran_seed(SKIP*lseed);
  for(i=0; i<nobj;) {
restart:
	  u1 = u1max*ran1(&seed);
	  v1 = 2.0*v1max*(ran1(&seed) - 0.5);
	  R = u1;
	  z = R*tan(v1);
	  
	  rhotst = bulgedens_(&R,&z);
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
	  fmax = dfbulge_(&psi) - fcut_bulge;
	  while( frand > f0 ) {
	    v2 = 1.1*vmax2;
	    while( v2 > vmax2 ) {
	      vR = 2*vmax*(ran1(&seed) - 0.5);
	      vp = 2*vmax*(ran1(&seed) - 0.5);
	      vz = 2*vmax*(ran1(&seed) - 0.5);
	      v2 = vR*vR + vp*vp + vz*vz;
	      E = psi - 0.5*v2;
	    }
	    f0 = dfbulge_(&E)-fcut_bulge;
	    frand = fmax*ran1(&seed);
		if( j > 10000 ) {
			fprintf(stderr,"j > 10000 E=%g psic=%g psi=%g frand=%g fmax=%g f0=%g\n",
				E,psic,psi,frand,fmax,f0);
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
	return 0;
}

