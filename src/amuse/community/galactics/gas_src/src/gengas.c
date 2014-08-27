#include "main.h"

void generate_gas_(int *nobj, float* vdisp, long long * seed, 
 float *m,float *x,float *y,float *z,float *vx,float *vy,float *vz, float *cs);

main(int argc,char** argv)
{
	int i, j, k, nobj=10000;
  long long lseed=0;
  float vdisp,t;
  float *m,*x,*y,*z,*vx,*vy,*vz, *cs;
  int err=0;
  
  iquery("Enter the number of particles",&nobj);
  fquery("Enter velocity dispersion (can be zero)",&vdisp);
  llquery("Enter negative integer seed",&lseed);
  if(lseed<0) lseed*=-1;
  m=(float*) malloc(sizeof(float)*nobj);
  x=(float*) malloc(sizeof(float)*nobj);
  y=(float*) malloc(sizeof(float)*nobj);
  z=(float*) malloc(sizeof(float)*nobj);
  vx=(float*) malloc(sizeof(float)*nobj);
  vy=(float*) malloc(sizeof(float)*nobj);
  vz=(float*) malloc(sizeof(float)*nobj);
  cs=(float*) malloc(sizeof(float)*nobj);
  if(!m || !x || !y || !z || !vx || !vy || !vz || !cs)
  {
    fprintf(stderr,"memory allocation error");
    return 1;
  }
  generate_gas_(&nobj,&vdisp,&lseed,m,x,y,z,vx,vy,vz,cs);
  t=0.0;
#ifdef ASCII
	fprintf(stdout,"%d %15.7e\n",nobj,t);
	for(i=0; i<nobj; i++) {
	  fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
		  m[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], cs[i]);
		  }
#else
	for(i=0; i<nobj; i++) {
	  fwrite(m+i,sizeof(float),1,stdout);
	  fwrite(x+i,sizeof(float),1,stdout);
	  fwrite(y+i,sizeof(float),1,stdout);
	  fwrite(z+i,sizeof(float),1,stdout);
	  fwrite(vx+i,sizeof(float),1,stdout);
	  fwrite(vy+i,sizeof(float),1,stdout);
	  fwrite(vz+i,sizeof(float),1,stdout);
	  fwrite(cs+i,sizeof(float),1,stdout);
	}
#endif
  return 0;
}
