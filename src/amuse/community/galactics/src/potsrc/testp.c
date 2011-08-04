#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(argc,argv)
int argc;
char **argv;
{

	int i, nobj;
	float x, y, z, ax, ay, az, pot, phi;
	float t, m, x0, y0, z0, dx, dy, dz, r, r2, r3;
	float eps;
	float vx0, vy0, vz0;


	x = 10;
	y = 1; 
	z = 2;

	ax = ay = az = pot = 0;
/*
	fscanf(stdin,"%d %f\n",&nobj,&t);
*/
	while( fread(&m,sizeof(float),1,stdin) ) {
		fread(&x0,sizeof(float),1,stdin);
		fread(&y0,sizeof(float),1,stdin);
		fread(&z0,sizeof(float),1,stdin);
		fread(&vx0,sizeof(float),1,stdin);
		fread(&vy0,sizeof(float),1,stdin);
		fread(&vz0,sizeof(float),1,stdin);
		dx = x - x0;
		dy = y - y0;
		dz = z - z0;

		r2 = dx*dx + dy*dy + dz*dz;
		r = sqrt(r2);
		r3 = r2*r;
		phi = -m/r;
		ax += phi*dx/r2;
		ay += phi*dy/r2;
		az += phi*dz/r2;
		pot += phi;
	}

	fprintf(stdout,"x y z %g %g %g\n",x,y,z);
	fprintf(stdout,"ax ay az %g %g %g\n",ax,ay,az);
	fprintf(stdout,"pot %g\n",pot);

}
