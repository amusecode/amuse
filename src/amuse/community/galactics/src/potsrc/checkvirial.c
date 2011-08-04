#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(argc,argv)
int argc;
char **argv;
{

	float m, x, y, z, vx, vy, vz, ax, ay, az, pot;
	double ke, pe;

	if( argc != 2 ) {
		fprintf(stderr,"usage: checkvirial potential_filename < rvfile\n");
		exit(0);
	}

	InitRigidPotential(argv[1]);

	ke = 0.0;
	pe = 0.0;
	while( fread(&m,sizeof(float),1,stdin) ) {
		float x, y, z, vx, vy, vz;
		float R;

		fread(&x,sizeof(float),1,stdin);
		fread(&y,sizeof(float),1,stdin);
		fread(&z,sizeof(float),1,stdin);
		fread(&vx,sizeof(float),1,stdin);
		fread(&vy,sizeof(float),1,stdin);
		fread(&vz,sizeof(float),1,stdin);
		getforce(x,y,z,&ax,&ay,&az,&pot);
		
		ke += (double) (0.5*m*(vx*vx + vy*vy + vz*vz));
		pe += (double) (0.5*m*pot);
	}
	fprintf(stderr,"T %g W %g -2T/W %g\n",ke,pe,-2*ke/pe);
	exit(0);
}

