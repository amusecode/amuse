#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float ra[1000], va[1000];
float rh[1000], vh[1000];
float rd[1000], vd[1000];
float rb[1000], vb[1000];
float A, B, C, D1, D2, D3, q2, v0, v02, psi0;
float rho1=100, sigb=0.3, sigb2, psicutb= -1.0, fbulgeconst;
float rscale=.15, vscale=.8164866, mscale=0.1;

main(argc,argv)
int argc;
char **argv;
{
	int i;
	float dr;
	char harmfile[80];

	if( argc != 3 ) {
		fprintf(stderr,"usage: vtst r_scale mass_scale > vr.p\n");
		exit(0);
	}
	rscale = atof(argv[1]);
	mscale = atof(argv[2]);
	vscale = sqrt(mscale/rscale);
#ifdef PGPLOT
	pgbegin(0,"vr.p",1,1);
	pgscf(2);
    pgsetc(1.5);
	pgenv(0.0,10.,0.0,1.0,0,0);
	pglabel("R","V","");
#endif
    plotv(ra,va,1,"dbh.dat");
    plotv(rb,vb,2,"b.dat");
    plotv(rh,vh,3,"h.dat");
	for(i=0; i<1000; i++) {
		vd[i] = sqrt(va[i]*va[i] - vb[i]*vb[i] - vh[i]*vh[i]);
		rd[i] = ra[i];
		fprintf(stdout,"%g %g %g %g %g\n",rd[i],vd[i],vb[i],vh[i],va[i]);
	}
#ifdef PGPLOT
	pgsls(4);
	pgline(1000,rd,vd);

	pgend();
#endif
}

plotv(rp,vp,line,file)
float *rp, *vp;
int line;
char *file;
{
	int i;
	float dr;
	float v02;
	float q;
	float rad, z, fr, fz, pot;

	readharmfile_(file,&A,&B,&C,&v0,&q,&psi0);
	v02 = v0*v0;

	dr = .05;
	for(i=0; i<1000; i++) {
		rad = i*dr;
		z = 0.0;
		force_(&rad,&z,&fr,&fz,&pot);
		rp[i] = rad;
		vp[i] = sqrt(-rp[i]*fr);
		rp[i] *= rscale;
		vp[i] *= vscale;
	}
#ifdef PGPLOT
	pgsls(line);
	pgline(1000,rp,vp);
#endif
}
