#include <stdio.h>
#include <tgmath.h>
#include "evolve.h"

#define TOLERANCE  (sizeof(DOUBLE)<8? 1.e-6:1.e-15)
#define ORDER  4
#define MAXITER 60

DOUBLE stumpff_C(DOUBLE z)
{
  if(z>0.04) return (1-cos(sqrt(z)))/z;  
  if(z<-0.04) return -(cosh(sqrt(-z))-1)/z;  
  return 1.L/2-z/24+z*z/720-z*z*z/40320+z*z*z*z/3628800-z*z*z*z*z/479001600;
}

DOUBLE stumpff_S(DOUBLE z)
{
  DOUBLE sz;
  if(z>0.04)
  {
    sz=sqrt(z);  
    return (sz-sin(sz))/(sz*z);  
  }
  if(z<-0.04)
  {
    sz=sqrt(-z);  
    return (sz-sinh(sz))/(sz*z);  
  }
  return 1.L/6-z/120+z*z/5040-z*z*z/362880+z*z*z*z/39916800-z*z*z*z*z/6227020800;
}

DOUBLE stumpff_C_prime(DOUBLE z)
{
  return (stumpff_C(z)-3*stumpff_S(z))/2/z; 
}

DOUBLE stumpff_S_prime(DOUBLE z)
{
  return (1-stumpff_S(z)-2*stumpff_C(z))/2/z; 
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

int laguerre(DOUBLE x0, DOUBLE *x, DOUBLE *arg,DOUBLE xtol,DOUBLE ytol, 
                  ftype *f, ftype *fprime, ftype *fprimeprime)
{ 
  int i=0;
  DOUBLE fv,dfv,ddfv,delta;//,g,h,hg;
  *x=x0;
  while( i< MAXITER)
  {
    fv=(*f)(*x,arg);
    if(fabs(fv) <= ytol) return 0;
    dfv=(*fprime)(*x,arg);
    ddfv=(*fprimeprime)(*x,arg);
    if(dfv==0) return -2;
    if(i<=MAXITER/4)
    {
    delta=-ORDER*fv/(dfv+
      SIGN(dfv)*sqrt(fabs((ORDER-1)*(ORDER-1)*dfv*dfv-ORDER*(ORDER-1)*fv*ddfv)));
    } else{
    delta=-fv/dfv;
    }
    (*x)+=delta;
//    if(i>MAXITER/4) printf("%d %14.12Lg | %14.12Lg, %14.12Lg\n", i,delta,*x,xtol);
//    if(i>MAXITER/4) printf("   %14.12Lg %14.12Lg, %14.12Lg\n", fv,dfv,ddfv);
    if(fabs(delta)<xtol)
    {
      return 0;
    }
    i+=1;
  }
  return -1;
} 

#define FAC 3

/* brackets root */
int bracketroot(DOUBLE *x1, DOUBLE *x2, DOUBLE *arg,
                  ftype *f)
{
  DOUBLE f1,f2;
  int i; //s
  if((*x2)==(*x1))return -2;
  f1=(*f)(*x1,arg);
  f2=(*f)(*x2,arg);
  for(i=0;i<MAXITER;i++)
  {
    if(f1*f2 < 0) return 0;  
    if(fabs(f1)<fabs(f2))
    {
      (*x1)+=FAC*(*x1-*x2);
      f1=(*f)(*x1,arg);    
    } else
    {
      (*x2)+=FAC*(*x2-*x1);
      f2=(*f)(*x2,arg);    
    }
  }
  return -1;
}

int saferoot(DOUBLE *x, DOUBLE x1, DOUBLE x2, DOUBLE *arg,
                  DOUBLE xtol,DOUBLE ytol, ftype *f, ftype *fprime)
{
  int i;
  DOUBLE fv,df,f1,f2,xl,xh,dx,dxold,temp;
  f1=(*f)(x1,arg);
  f2=(*f)(x2,arg);
  if(f1==0){(*x)=x1;return 0;}
  if(f2==0){(*x)=x2;return 0;}
  if(f1*f2>0) return -1;
  if(f1<0) {xl=x1;xh=x2;}
  else {xl=x2;xh=x1;}
  (*x)=(xl+xh)/2;
  dxold=fabs(xl-xh);
  dx=dxold;
  fv=(*f)(*x,arg);
  df=(*fprime)(*x,arg);
  for(i=0;i<MAXITER;i++)
  {
    if( (((*x)-xh)*df-fv)*(((*x)-xl)*df-fv) > 0   || 
        fabs(2*fv) > fabs(dxold*df))
    {
      dxold=dx;
      dx=(xh-xl)/2;
      (*x)=xl+dx;
      if(xl==*x || xh==*x) return 0;
    } else
    {
      dxold=dx;
      dx=fv/df;
      temp=(*x);
      (*x)-=dx;
      if(temp==(*x)) return 0;
    }
    if(fabs(dx)<xtol) return 0;
    fv=(*f)(*x,arg);
    if(fabs(fv)<=ytol) return 0;    
    df=(*fprime)(*x,arg);
    if(fv<0) xl=*x;
    else xh=*x;
  } 
  return -2;
}

int findroot(DOUBLE x0, DOUBLE *x, DOUBLE *arg,DOUBLE xtol,DOUBLE ytol, 
                  ftype *f, ftype *fprime, ftype *fprimeprime)
{
 int err;
 DOUBLE x1,x2;
 err=laguerre(x0, x, arg,xtol,ytol, f, fprime, fprimeprime);
 if(err==0) return err;
 if(x0==0) x0=1.;
 x1=0.5*x0;x2=x0;
 err=bracketroot(&x1,&x2,arg,f);
 if(err!=0) return err;
 return saferoot(x,x1,x2,arg,xtol,ytol,f,fprime);
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
  DOUBLE smu=sqrt(mu);
  DOUBLE r0=sqrt(pos0[0]*pos0[0]+pos0[1]*pos0[1]+pos0[2]*pos0[2]);
  DOUBLE v0=sqrt(vel0[0]*vel0[0]+vel0[1]*vel0[1]+vel0[2]*vel0[2]);
  DOUBLE vr0=(pos0[0]*vel0[0]+pos0[1]*vel0[1]+pos0[2]*vel0[2])/r0;
  DOUBLE alpha=2./r0-v0*v0/mu;
  DOUBLE xi0,arg[5],xi,xtol,ytol,dxi0;
  int err;

  if(alpha > 0)
  {
    xi0=smu*dt*alpha;
  } else
  {
    xi0=SIGN(dt)/sqrt(-alpha)*log(1-2*mu*dt*alpha/((vr0*r0)+  
          SIGN(dt)*smu/sqrt(-alpha)*(1-r0*alpha)));
// this last formula is 4.5.11 in bate et al., fundamentals of astrodynamics 
// with +1 in the logarithm
    dxi0=smu/r0*dt;
    if(fabs(alpha*dxi0*dxi0)<1) xi0=dxi0;
  }

//  xi0=smu*dt/r0;
  
  arg[0]=smu*dt;
  arg[1]=r0;
  arg[2]=vr0;
  arg[3]=smu;
  arg[4]=alpha;
  ytol=fabs(TOLERANCE*smu*dt);
  xtol=fabs(TOLERANCE*smu*dt/r0);

  err=findroot(xi0, &xi, arg, xtol,ytol,&f,&fprime,&fprimeprime);
  if(err !=0 || SIGN(xi)!=SIGN(dt)) {
   printf("xtol,ytol: %g %g\n",xtol,ytol);
   printf("err: %d %g %g\n",err,smu*dt/r0,alpha); 
   printf("%20.16g %20.16g arg: %20.16g %20.16g %20.16g %20.16g %20.16g\n", 
            xi, xi0, arg[0], arg[1], arg[2], arg[3], arg[4]);
   printf("%20.16g %20.16g %20.16g %20.16g %20.16g %20.16g\n", 
            pos0[0], pos0[1], pos0[2], vel0[0], vel0[1], vel0[2]);

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
    r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    for(i=0;i<3;i++) vel[i]=pos0[i]*smu/r*ldf+vel0[i]*smu/r*ldg;  
  }
//  printf("out: %f %f %f %f\n", (float) dt,(float)mu,(float)pos[0],(float)vel[0]);

  return 0;
}
