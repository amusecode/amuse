#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void InitRigidPotential(potname)
char *potname;
{

	readharmfile_(potname);
	fixpot_();

}

void getforce(x,y,z,ax,ay,az,pot)
float x, y, z, *ax, *ay, *az, *pot;
{
	float r0, x0, y0, z0, ax0, ay0, az0, ar0, pot0;

	x0 = x; y0 = y; z0 = z;

	r0 = sqrt(x0*x0 + y0*y0);
	force_(&r0,&z0,&ar0,&az0,&pot0);

	*pot = -pot0;	
	*ax = ar0*x0/r0;
	*ay = ar0*y0/r0;
	*az = az0;

}
