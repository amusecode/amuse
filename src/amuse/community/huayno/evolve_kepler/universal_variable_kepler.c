#include <stdio.h>
#include <math.h>

#define TOLERANCE  1.e-15
#define ORDER  5
#define MAXITER 30

#define SIGN(x)  ((x)>=0? 1:-1)

DOUBLE stumpff_C(DOUBLE z)
{
  if(z>0) return (1-DBL(cos)(DBL(sqrt)(z)))/z;  
  if(z<0) return -(DBL(cosh)(DBL(sqrt)(-z))-1)/z;  
  return 1/2.;
}

DOUBLE stumpff_S(DOUBLE z)
{
  double sz;
  if(z>0)
  {
    sz=DBL(sqrt)(z);  
    return (sz-DBL(sin)(sz))/(sz*z);  
  }
  if(z<0)
  {
    sz=DBL(sqrt)(-z);  
    return (sz-DBL(sinh)(sz))/(sz*z);  
  }
  return 1/6.;
}

DOUBLE stumpff_C_prime(DOUBLE z)
{
  return (stumpff_C(z)-3*stumpff_S(z))/2./z; 
}

DOUBLE stumpff_S_prime(DOUBLE z)
{
  return (1-stumpff_S(z)-2*stumpff_C(z))/2./z; 
}

DOUBLE lagrange_f(DOUBLE xi, DOUBLE r0, DOUBLE vr0, DOUBLE smu, DOUBLE alpha)
{
  return 1.-xi*xi/r0*stumpff_C(alpha*xi*xi);    
}

DOUBLE lagrange_dfdxi(DOUBLE xi, DOUBLE r0, DOUBLE vr0, DOUBLE smu, DOUBLE alpha)
{
  DOUBLE z=alpha*xi*xi;
  return xi/r0*(z*stumpff_S(z)-1);    
}

DOUBLE lagrange_g(DOUBLE xi, DOUBLE r0, DOUBLE vr0, DOUBLE smu, DOUBLE alpha)
{
  DOUBLE z=alpha*xi*xi, xi_smu=xi/smu;
  return r0*vr0*xi_smu*xi_smu*stumpff_C(z)-r0*xi_smu*z*stumpff_S(z)+r0*xi_smu;
}

DOUBLE lagrange_dgdxi(DOUBLE xi, DOUBLE r0, DOUBLE vr0, DOUBLE smu, DOUBLE alpha)
{
  DOUBLE z=alpha*xi*xi,_smu=1./smu;
  return r0*vr0*_smu*_smu*xi*(1-z*stumpff_S(z))-z*r0*_smu*stumpff_C(z)+r0*_smu;
}

DOUBLE universal_kepler(DOUBLE xi,DOUBLE r0,DOUBLE vr0,DOUBLE smu,DOUBLE alpha)
{
  DOUBLE z=alpha*xi*xi;
  return (r0*vr0*xi*xi*stumpff_C(z)/smu+  
           (1-alpha*r0)*xi*xi*xi*stumpff_S(z)+r0*xi);
}

DOUBLE universal_kepler_dxi(DOUBLE xi,DOUBLE r0,DOUBLE vr0,DOUBLE smu,DOUBLE alpha)
{
  DOUBLE z=alpha*xi*xi;
  return (r0*vr0*xi*(1-z*stumpff_S(z))/smu +  
           (1-alpha*r0)*xi*xi*stumpff_C(z)+r0);
}

DOUBLE universal_kepler_dxidxi(DOUBLE xi,DOUBLE r0,DOUBLE vr0,DOUBLE smu,DOUBLE alpha)
{
  return -alpha*universal_kepler(xi,r0,vr0,smu,alpha)+r0*vr0/smu+xi;
}

typedef DOUBLE (ftype)(DOUBLE x, DOUBLE* arg);

int laguerre(DOUBLE x0, DOUBLE *x, DOUBLE *arg, 
                  ftype *f, ftype *fprime, ftype *fprimeprime)
{ 
  int i=0;
  DOUBLE fv,dfv,ddfv,delta;
  *x=x0;
  while( i< MAXITER)
  {
    fv=(*f)(*x,arg);
    dfv=(*fprime)(*x,arg);
    ddfv=(*fprimeprime)(*x,arg);
//    printf("a: %d %f %f %f %f\n", i,(float) DBL(fabs)(*x),(float)fv,(float)dfv,(float)ddfv);
    if(dfv==0 || ddfv==0) return -2;
    delta=-ORDER*fv/(dfv+
      SIGN(dfv)*DBL(sqrt)(DBL(fabs)((ORDER-1)*(ORDER-1)*dfv*dfv-ORDER*(ORDER-1)*fv*ddfv)));
//    delta=-fv/dfv;
    (*x)+=delta;
//    printf("%d %14.12Lg | %14.12Lg, %14.12Lg\n", i,DBL(fabs)(delta),*x,TOLERANCE*DBL(fabs)(*x));
    if(DBL(fabs)(delta)<TOLERANCE*DBL(fabs)(*x))
    {
      return 0;
    }
    i+=1;
  }
  return -1;
} 

DOUBLE f(DOUBLE xi, DOUBLE *arg)
{
  return universal_kepler(xi,arg[1],arg[2],arg[3],arg[4])-arg[0];
}

DOUBLE fprime(DOUBLE xi, DOUBLE *arg)
{
  return universal_kepler_dxi(xi,arg[1],arg[2],arg[3],arg[4]);
}

DOUBLE fprimeprime(DOUBLE xi, DOUBLE *arg)
{
  return universal_kepler_dxidxi(xi,arg[1],arg[2],arg[3],arg[4]);
}

int universal_variable_kepler_solver(DOUBLE dt,DOUBLE mu,DOUBLE pos0[3], 
             DOUBLE vel0[3],DOUBLE pos[3], DOUBLE vel[3])
{
  DOUBLE smu=DBL(sqrt)(mu);
  DOUBLE r0=DBL(sqrt)(pos0[0]*pos0[0]+pos0[1]*pos0[1]+pos0[2]*pos0[2]);
  DOUBLE v0=DBL(sqrt)(vel0[0]*vel0[0]+vel0[1]*vel0[1]+vel0[2]*vel0[2]);
  DOUBLE vr0=(pos0[0]*vel0[0]+pos0[1]*vel0[1]+pos0[2]*vel0[2])/r0;
  DOUBLE alpha=2./r0-v0*v0/mu;
  DOUBLE xi0,dxi0,arg[5],xi;
  int err;
//  printf("in: %f %f %f %f\n", (float) dt,(float)mu,(float)pos0[0],(float)vel0[0]);

  if(alpha > 0)
  {
    xi0=smu*alpha*dt;
  } else
  {
    xi0=SIGN(dt)/DBL(sqrt)(-alpha)*log(1-2*mu*dt*alpha/((vr0*r0)+  
          SIGN(dt)*smu/DBL(sqrt)(-alpha)*(1-r0*alpha)));
// this last formula is 4.5.11 in bate et al., fundamentals of astrodynamics 
// with +1 in the logarithm
    dxi0=smu/r0*dt;
    if(DBL(fabs)(alpha*dxi0*dxi0)<1) xi0=dxi0;
  }
  
  arg[0]=smu*dt;
  arg[1]=r0;
  arg[2]=vr0;
  arg[3]=smu;
  arg[4]=alpha;

  err=laguerre(xi0, &xi, arg, &f,&fprime,&fprimeprime);
  if(err !=0) {
   printf("%f arg: %f %f %f %f %f\n", (float) xi0,(float) arg[0], (float) arg[1], 
  (float) arg[2]  , (float) arg[3], (float) arg[4]);

  return err;
  }
  {
    int i;
    DOUBLE r;
    DOUBLE lf=lagrange_f(xi,r0,vr0,smu,alpha);
    DOUBLE lg=lagrange_g(xi,r0,vr0,smu,alpha);
    DOUBLE ldf=lagrange_dfdxi(xi,r0,vr0,smu,alpha);
    DOUBLE ldg=lagrange_dgdxi(xi,r0,vr0,smu,alpha);
    for(i=0;i<3;i++) pos[i]=pos0[i]*lf+vel0[i]*lg;
    r=DBL(sqrt)(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    for(i=0;i<3;i++) vel[i]=pos0[i]*smu/r*ldf+vel0[i]*smu/r*ldg;  
  }
//  printf("out: %f %f %f %f\n", (float) dt,(float)mu,(float)pos[0],(float)vel[0]);

  return 0;
}
