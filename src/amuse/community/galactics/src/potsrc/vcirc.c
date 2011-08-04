#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float ra[10000], va[10000];
float rh[10000], vh[10000];
float rd0[10000], vd[10000], vd0[10000];
float rb[10000], vb[10000];

main(argc,argv)
int argc;
char **argv;
{
	int i;
	float dr;
	char harmfile[80];
	FILE *vrfile;

    plotv(ra,va,1,"dbh.dat");
    plotv(rb,vb,2,"b.dat");
    plotv(rh,vh,3,"h.dat");
/*
	plotv(rd0,vd0,4,"d.dat");
*/
	for(i=0; i<5000; i++) {
		float vd2;
		vd2 = (va[i]*va[i] - vb[i]*vb[i] - vh[i]*vh[i]);
		if( vd2 > 0 ) {
			vd[i] = sqrt(va[i]*va[i] - vb[i]*vb[i] - vh[i]*vh[i]);
		}
		else {
			vd[i] = 0;
		}
		rd0[i] = ra[i];
		fprintf(stdout,"%g %g %g %g %g\n",rd0[i],vd[i],vb[i],vh[i],va[i]);
	}
	exit(0);
}

plotv(rp,vp,line,file)
float *rp, *vp;
int line;
char *file;
{
	int i;
	float dr;
	float q;
	float rad, z, fr, fz, pot;

	InitRigidPotential(file);
	
	dr = 0.1;
	for(i=0; i<5000; i++) {
		rad = i*dr;
		z = 0.0;
		force_(&rad,&z,&fr,&fz,&pot);
		if( fr > 0 ) fr = 0;
		rp[i] = rad;
		vp[i] = sqrt(-rp[i]*fr);
	}
#ifdef PGPLOT
	pgsls(line);
	pgline(1000,rp,vp);
#endif
}





