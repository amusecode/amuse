#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(argc,argv)
int argc;
char **argv;
{

	float x, y, z, ax, ay, az, pot;

	if( argc != 2 ) {
		fprintf(stderr,"usage: testforce potential_filename\n");
		exit(0);
	}

	InitRigidPotential(argv[1]);

	x = 154;
	y = -137.58; 
	z = -2.78;
	getforce(x,y,z,&ax,&ay,&az,&pot);
	fprintf(stderr,"%g %g %g\n",x,y,z);
	fprintf(stderr,"%g %g %g %g\n",ax,ay,az,sqrt(2*pot));
	exit(0);

	for(x=0.01; x<400; x += 1 ) {
		getforce(x,y,z,&ax,&ay,&az,&pot);
		fprintf(stdout,"%g %g %g\n",x,sqrt(-x*ax),pot);
	}
exit(0);

	fprintf(stdout,"x y z %g %g %g\n",x,y,z);
	fprintf(stdout,"ax ay az %g %g %g\n",ax,ay,az);
	fprintf(stdout,"pot %g\n",pot);

}
