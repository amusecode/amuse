#include "main.h"
#ifndef NOMPI
#include <mpi.h>
#endif

float rad;
float Rc, zc;

main(argc,argv)
int argc;
char **argv;
{
	int i, j, k, nobj=10000;
	int seed= -123;
	int icofm=1;
	int nrtrunc;
	float q=0.9, rtrunc=5.0;
	float rhoran, massapp;
	float x, y, z, vx, vy, v, v2, R, vmax, vmax2;
	float phi, cph, sph, cth, sth, vR, vp, vz;
	float E, Lz, rad, rad2;
	float f0, frand, fmax, fmax1, vphimax1, psi, fnorm;
	float fr;
	float Fdisk();
	float ran1();
	float invu();
	float dr, rhomax1;
	float t, mass;
	float u1, v1, u1max, v1max;
	float xcm, ycm, zcm, vxcm, vycm, vzcm;
	float zip = 0.0, psi0;
	float rhomax=0.15, rhomin, rhotst;
	float rhoguess, zcon;
	float stream=0.5;
	float con, outdisk, drtrunc;
	float omega, kappa;
	float dv, v0;
	float broadcon=1.0;
	float vfacR, vfacz;
	float f1;
	float gr, gp, gz, g2;
	float diskdensf_(), sigr2_(), sigz2_(), pot();
	float FindMax(), FindMax1(), Fmax();
	float simpson(), massring(); 
	float force_();
	float rcirc_();
	float vcirc, rpot, fz;
	float vphimaxold;
	float rcircular, angmom;
	float vsigzz;
	FILE *rhoout;
	FILE *errfile;
	char harmfile[80];
	float kinetic, potential, rcyl, kineticz, potentialz;
	int myid, numprocs;
	FILE *infile, *rvfile;
	char rvname[80];

#ifndef NOMPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#else
	numprocs = 1;
	myid = 0;
#endif
	if( argc == 3 ) {
		infile = fopen(argv[1],"r");
		if( myid == 0 )
			rvfile = fopen(argv[2],"w");
	}
	else {
		fprintf(stderr,"usage: pgendisk in.disk rvfilename\n");
		exit(0);
	}

	errfile = fopen("errmsg.dat","w");
	strcpy(harmfile,"dbh.dat");
#ifdef never
	iquery("Enter the number of particles",&nobj);
	iquery("Enter negative integer seed",&seed);
	iquery("Center the simulation (0=no,1=yes)",&icofm);
	cquery("Enter harmonics file",harmfile);
#endif

	fscanf(infile,"%d\n",&nobj);
	fscanf(infile,"%d\n",&seed);
	fscanf(infile,"%d\n",&icofm);
/*
	fscanf(infile,"%s\n",harmfile);
*/
	seed = seed*(1+myid); /* different random seed for each process */

	r = (phase *) calloc(nobj,sizeof(phase));

	readdiskdf_(harmfile,&gparam);
	mdisk = gparam.mdisk; rdisk = gparam.rdisk; zdisk = gparam.zdisk;
	outdisk = gparam.outdisk; drtrunc = gparam.drtrunc;

	rd = 1.2*rdisk;
	zd = 1.2*zdisk;
	rtrunc = (outdisk + 2.0*drtrunc);
	diskmass = 4.0*M_PI*simpson(massring,0.0,rtrunc,128);
	mass   = diskmass/(nobj*numprocs);

	nrtrunc = 10000;
	dr = rtrunc/nrtrunc;
	rhomax = 0.0;
#ifdef DEBUG
	rhoout = fopen("rhotst.dat","w");
#endif
	for(i=0; i<nrtrunc; i++) {
	  R = i*dr;
	  z = 0.0;
	  rhoguess = exp(-R/rd);
	  rhotst = diskdensf_(&R,&z);
	  rhotst /= rhoguess;
	  if( rhotst > rhomax ) rhomax = rhotst;
	  //	  fprintf(rhoout,"%g %g\n",R, rhotst);
	}
	rhomin = 0.0;
	rhomax *=1.2;

	kinetic = potential = 0.;

	fprintf(stderr,"Calculating disk positions and velocities - seed=%d\n",seed);
	for(i=0, j=0, k=0; i<nobj;) {
	  fmax = -1.;
	  while (fmax <= 0.) {
	    R = 2.*rtrunc;
	    while (R > rtrunc) {
	      u1 = -ran1(&seed);
	      v1 = 2.0*(ran1(&seed) - 0.5);
	      R = rd*invu(u1);
	      z = zd*atanh(v1);
	    }
	    zcon = cosh(z/zd);
	    rhoguess = exp(-R/rd)/(zcon*zcon);
	    rhotst = diskdensf_(&R,&z);
	    rhotst /= rhoguess;

	    k++;
	    if( rhotst < rhomin )
	      continue;
	    
	    rhoran = (rhomax - rhomin)*ran1(&seed);
	    if( rhoran > rhotst )
	      continue;
	    phi = 2.0*M_PI*ran1(&seed);
	    x = R*cos(phi);
	    y = R*sin(phi);
	    omekap_(&R, &omega, &kappa);
	    vphimax = omega*R;	  
	    vsigR = sqrt(sigr2_(&R));
	    vsigp = kappa/(2.0*omega)*vsigR;
	    vsigz = sqrt(sigz2_(&R));
	    vphimaxold = vphimax;
	    fmax = 1.1*FindMax1(R,z,&vphimax);
#ifdef DEBUG
	    fprintf(rhoout,"%d %g %g %g %g\n",i,fmax,R,z,vphimax);
#endif
	  }
	  vphimax = vphimaxold;

	  f0 = 0.0; frand = 1.0; /* dummy starters */
	  while( frand > f0 ) {
	    g2 = 999.;
	    while( g2 > 1.0) {
	      gr = 8.0*(ran1(&seed) - 0.5);
	      gp = 16.0*(ran1(&seed) - 0.5);
	      gz =  8.0*(ran1(&seed) - 0.5);
	      g2 = gr*gr/16. + gp*gp/64. + gz*gz/16.;
	    }
	    vR = broadcon*vsigR*gr;
	    vp = vphimax + broadcon*vsigp*gp;
	    vz = broadcon*vsigz*gz;

	    f0 = Fdisk(vR, vp, vz, R, z);
	    frand = fmax*ran1(&seed);
	    if( f0 > fmax ) {
	      float vpmax;
	      fprintf(errfile,"f0 > fmax at R=%g z=%g\nvr=%g vp=%g, vz=%g vphimax=%g f0=%g, fmax=%g\n", 
		      R,z, 
		      vR*broadcon/vsigR, 
		      (vp - vphimax)*broadcon/vsigp, 
		      vz*broadcon/vsigz, vphimax, f0, fmax);
	      fflush(errfile);
	      }
	    j++;
	  }
	  vphimax = vp - broadcon*vsigp*gp;
	  vp = vphimax + broadcon*vsigp*gp;
	  cph = x/R; sph = y/R;
	  vx = vR*cph - vp*sph;
	  vy = vR*sph + vp*cph;

	  r[i].x = (float) x;
	  r[i].y = (float) y;
	  r[i].z = (float) z;
	  r[i].vx = (float)vx;
	  r[i].vy = (float)vy;
	  r[i].vz = (float)vz;
	  i++;
	  if( i % 1000 == 0 ) {
	    fprintf(stderr,".");
	    fflush(stderr);
	  }
	}
	fprintf(stderr,"\n");
	fprintf(stderr,"number of density trials %d\n",k);
	fprintf(stderr,"number of velocity trials %d\n",j);
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
	fprintf(stdout,"%d %15.7e\n",nobj,t);
	for(i=0; i<nobj; i++) {
	  fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
		  mass, r[i].x, r[i].y, r[i].z, r[i].vx, r[i].vy, r[i].vz);
		  }
#else
	{
		if( myid == 0 ) {
			for(i=0; i<nobj; i++) {
				fwrite(&mass,sizeof(float),1,rvfile);
				fwrite(r+i,sizeof(phase),1,rvfile);
			}

#ifndef NOMPI
			for(j=1; j<numprocs; j++) {
				MPI_Status receive_status;

				fprintf(stderr,"Receiving %d\n",j);
				MPI_Recv(r,nobj*sizeof(phase), 
					MPI_BYTE,j,999,MPI_COMM_WORLD,&receive_status);
				for(i=0; i<nobj; i++) {
					fwrite(&mass,sizeof(float),1,rvfile);
					fwrite(r+i,sizeof(phase),1,rvfile);
				}
			}
#endif
			fclose(rvfile);
		}
		else {
#ifndef NOMPI
			fprintf(stderr,"Sending from %d\n",myid);
			MPI_Send(r, nobj*sizeof(phase),MPI_BYTE,0,999,MPI_COMM_WORLD); 
#endif
		}
	}
#endif
#ifndef NOMPI
	MPI_Finalize();
#endif
	exit(0);
}

float Fdisk(vr, vp, vz, R, z)
float vr, vp, vz, R, z;
{
  float vr0, vp0, vz0, R0, z0;
  float diskdf5ez_();
  
  vr0 = vr; vp0 = vp; vz0 = vz; R0 = R; z0 = z;
  return diskdf5ez_(&vr0, &vp0, &vz0, &R0, &z0);
}

float invu(u)
     float u;
{
  /* invert the u function to find R */
  int i;
  double dr, rg, rnew, eps=1.0e-8;
  
  rg = 1.0; 
  dr = 1.0;
  eps = 1.0e-8;
  for(i=0; (i<20 && dr > eps); i++) {
    rnew = rg - (-(1.0 + rg)*exp(-rg) - u)
      /(rg*exp(-rg));
    dr = fabs(rnew - rg);
    rg = rnew;
  }
  if( i == 20 ) fprintf(stderr,"Warning R did not converge\n");
  return (float) rg;
}

/* A approximate estimate of the local maximum of the distribution function */

float Fmax(R,z,vR,vz,vlimit, vpmax)
     float R, z, vR, vz;
     float *vpmax, vlimit;
{
  int i;
  float v, dv, fmax, vmax, v0;
  float Fdisk();
  
  dv = 2.0*vlimit/100;
  v0 = *vpmax - 50*dv;
  fmax = Fdisk(vR, v0, vz, R, z);
  for(i=0; i<100; i++) {
    float ftry, vtry;
    
    vtry = v0 + i*dv;
    ftry = Fdisk(vR, vtry, vz, R, z);
    if( ftry > fmax ) {
      fmax = ftry;
      vmax = vtry;
    }
  }
  *vpmax = vmax;
  
  return fmax;
}

float FindMax1(R,z,vpmax)

     float R, z, *vpmax;
{
  int upflag, flag;
  float dv, vpm, vpmold, v0, v1;
  float f0, f1, ftmp, fmid;
  float Fdisk();
  float zero;
  
  zero = 0.0;
  dv = 0.1*vsigp;
  vpm = *vpmax;
  v0 = vpm - dv;
  v1 = vpm + dv;
  
  f0 = Fdisk(0.0,v0,0.0,R,z);
  fmid = Fdisk(0.0,vpm,0.0,R,z);
  f1 = Fdisk(0.0,v1,0.0,R,z);
  if( fmid >= f0 && fmid >= f1 )
    return fmid;
  
  if( f0 > f1 ) {
    ftmp = f0;
    f0 = f1;
    f1 = ftmp;
    v1 = v0;
    dv = -dv;
  }
  vpm = v1;
  
  flag = 1;
  while( flag ) {
    
    dv *= 2.0;
    vpmold = vpm;
    vpm += dv;
    f0 = f1;
    f1 = Fdisk(0.0,vpm,0.0,R,z);
    flag = (f1 > f0);
  }
  *vpmax = vpmold;
  return f0;
}

float FindMax(R,z,vpmax)
     float R, z, *vpmax;
{
  float vpm;
  float v0, v1;
  float eps=0.01;
  float fmax, golden(), fgold();
  
  Rc = R; zc = z;
  
  vpm = *vpmax;
  v0 = *vpmax - 3.0*vsigp;
  v1 = *vpmax + 3.0*vsigp;
  fmax = -golden(v0,vpm,v1, fgold, eps, vpmax);
  
  return fmax;
  
}

float fgold(vp)
     float vp;
{
  float R, z;
  float zero;
  zero = 0.0;
  R = Rc; z = zc;
    return( -Fdisk(0.0,vp,0.0,R,z) );
}

float massring(r)
     float r;
{
  float denz();
  float simpson();
  
  rad = r;
  return( rad*simpson(denz,0.0,5.0*zdisk,128) );
}

float denz(z)
     float z;
{
  float z0;
  float diskdensf_();
  
  z0 = z;
  return diskdensf_(&rad, &z0);
}

