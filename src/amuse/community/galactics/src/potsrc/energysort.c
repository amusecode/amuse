#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	float m, x, y, z, vx, vy, vz;
} phase;

int *indx;
double *be;

main(argc,argv)
int argc;
char **argv;
{

	int i, j, nobj;
	int cmp();
	float m, x, y, z, vx, vy, vz, ax, ay, az, pot;
	double ke, pe;
	FILE *rvfile, *rvfile1;
	phase *r, *r0;

	if( argc != 4 ) {
		fprintf(stderr,"usage: energysort potential_filename init_rvfile final_rvfile > sorted_rvfile\n");
		exit(0);
	}

	InitRigidPotential(argv[1]);
	rvfile = fopen(argv[2],"r");
	fseek(rvfile,0,2); nobj = ftell(rvfile)/28; fseek(rvfile,0,0);
	rvfile1 = fopen(argv[3],"r");
	fseek(rvfile1,0,2); nobj = ftell(rvfile1)/28; fseek(rvfile1,0,0);

	r = (phase *) calloc(nobj,sizeof(phase));
	r0 = (phase *) calloc(nobj,sizeof(phase));
	be = (double *) calloc(nobj,sizeof(double));
	indx = (int *) calloc(nobj,sizeof(int));
	fread(r,sizeof(phase),nobj,rvfile);
	fclose(rvfile);

	ke = 0.0;
	pe = 0.0;
	for(i=0; i<nobj; i++) {
		x = r[i].x; y = r[i].y; z = r[i].z;
		vx = r[i].vx; vy = r[i].vy; vz = r[i].vz;

		getforce(x,y,z,&ax,&ay,&az,&pot);
		
		ke = (double) 0.5*(vx*vx + vy*vy + vz*vz);
		pe = (double) pot;
		be[i] = ke + pe;
		indx[i] = i;
		
	}
	qsort(indx,nobj,sizeof(int),cmp);

	fread(r,sizeof(phase),nobj,rvfile1);
	fclose(rvfile1);

	for(i=0; i<nobj; i++) {
		j = indx[i];
		r0[i] = r[j];
	}
	fwrite(r0,sizeof(phase),nobj,stdout);

	exit(0);
}

int cmp(a,b)
int *a, *b;
{
	if( be[*a] > be[*b]) return 1;
	if( be[*a] < be[*b]) return -1;
	return 0;
}
