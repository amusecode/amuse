#include <stdio.h>
#include <stdlib.h>

extern float diskmass, diskedge;
extern float bulgemass, bulgeedge;
extern float halomass, haloedge;

readmassrad()
{
	FILE *mrfile;

	if( !(mrfile = fopen("mr.dat","r")) ) {
		fprintf(stderr,"Cannot open file mr.dat - contains component masses\n");
		exit(1);
	}

	fscanf(mrfile,"%f %f\n",&diskmass,&diskedge);
	fscanf(mrfile,"%f %f\n",&bulgemass,&bulgeedge);
	fscanf(mrfile,"%f %f\n",&halomass,&haloedge);
	fclose(mrfile);
}
