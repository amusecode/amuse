
/* calculates potential of thin axisymmetric disk */
/*  with exponential or user supplied density     */

/* gcc pot_schijf.cc -I/software/local/include -lgslcblas -L/software/local/lib -lgsl -lm
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>


#define 	SIG0	0.135
#define		RS	1.
#define 	NTHETA	100
#define		NR	100
#define		RMAX	15*RS


#define GSL_FN_EVAL(F,x)         (*((F)->function))(x,(F)->params)

class thinDisk {
 public:
  thinDisk(int nth=NTHETA,int nr=NR,double rmax=RMAX,double (*sig)(double)=&(sigma0));
  ~thinDisk(){ freemem();}
  double potentialC(double r, double z);
  double potentialP(double r, double th);
  static double integrantWrapper(double ra, void*p);
  
 private:
  int ntheta,nrad;
  double dtheta,drad,rmax;
  double (*sigma)(double r);
  double **potarray;
  static double sigma0(double x);    
  double integrant(double ra, void * p);
  double diskpotential( double r, double z,double rmax);
  void ToPolar(double &r, double &z);
  void freemem();
  void fillarray();
  void getmem();
};

struct int_params { thinDisk *Object; double r;double z;};

 thinDisk::thinDisk(int nth,int nr,double rmx,double (*sig)(double)):
  ntheta(nth),nrad(nr),dtheta(M_PI_2/double(nth)),drad(rmx/double(nr)),rmax(rmx),sigma(sig) 
{
 int i,j,k;
 getmem();
 fillarray();  
}

double thinDisk::integrantWrapper(double ra, void*p)
{
struct int_params * params
          = (struct int_params *)p;
  thinDisk* self=(thinDisk*) params->Object;
  return self->integrant(ra,params);
   
}


double thinDisk::sigma0(double x)
{
 double sigma;
 sigma=SIG0*exp(-x/RS);
 return sigma;
}


double thinDisk::potentialC(double r, double z)
{
 ToPolar(r,z);
 return potentialP(r,z);
}


void thinDisk::ToPolar(double &r, double &z)
{
 double x;
 x=sqrt(r*r+z*z);
 r=x; if(r==0){z=0.;return;};
 z=2*atan(z/(x+r));
}

double thinDisk::potentialP(double r, double th)
{
 int i0,j0,i1,j1;
 double t,u,result;

 j0=int(r/drad);i0=int(th/dtheta);
 j1=j0+1;i1=i0+1;
 if(i0<0) {i0=0;i1=0;};
 if(j0<0) {j0=0;j1=0;};
 if(j0>=nrad) {j0=nrad;j1=nrad;};
 if(i0>=ntheta) {i0=ntheta;i1=ntheta;};

 t=r/drad-j0;
 u=th/dtheta-i0;

 result=(1-t)*(1-u)*potarray[i0][j0]+(1-t)*u*potarray[i1][j0]
               +t*(1-u)*potarray[i0][j1]+t*u*potarray[i1][j1];
 return result;
}



void thinDisk::fillarray(){
 int i,j,k;
 double theta,r;
 printf(" Initializing thin disk potential");
 for(i=0;i<(ntheta+1);i++) {
  for(j=0;j<(nrad+1);j++){
  theta=i*M_PI_2/double(ntheta);
  r=j*rmax/double(nrad);
  potarray[i][j]=diskpotential( r*cos(theta),r*sin(theta),rmax);
  };
  if(i%(ntheta/10)==0){printf(".");fflush(stdout);}
 };
 printf("done\n");
}


void thinDisk::getmem()
{
 int i;
 potarray=new double*[(ntheta+1)];
 for(i=0;i<ntheta+1;i++){
  potarray[i]=new double[nrad+1];
 }
}


void thinDisk::freemem()
{
 int i;
 for(i=0;i<ntheta+1;i++){
 delete[] potarray[i];
 }
 delete[] potarray;
}


double thinDisk::integrant(double ra, void * p)
{
 struct int_params * params
          = (struct int_params *)p;

 double r = params->r;
 double z = params->z;
 double temp;
 double integrant;
 if(r==ra && z==0) {
  return sigma(ra);
 }
 temp=1/sqrt((ra+r)*(ra+r)+z*z);
 return sigma(ra)*ra*temp*
      gsl_sf_ellint_Kcomp(sqrt(4*r*ra)*temp,GSL_PREC_DOUBLE);
}

double thinDisk::diskpotential( double r, double z,double rmax)
{
double result, error;
gsl_function F;
struct int_params params ={ this,r, z};

gsl_integration_workspace * w = gsl_integration_workspace_alloc(100);

F.function=&integrantWrapper;
F.params=&params;

if(z!=0) {
gsl_integration_qag (&F, 0, rmax, 0, 1e-5, 100,GSL_INTEG_GAUSS15,w, &result, &error);
}
else
{
double interval[3]={0,r,rmax};
int ni=3;
gsl_integration_qagp (&F, &interval[0], ni, 0, 1e-5, 100,w, &result, &error);
}

gsl_integration_workspace_free(w);

return -4*result;
}
 
extern "C" float flatdiskpot_(float *s,float *z)
{
static thinDisk disk;
return disk.potentialC(*s,*z);
}

 
