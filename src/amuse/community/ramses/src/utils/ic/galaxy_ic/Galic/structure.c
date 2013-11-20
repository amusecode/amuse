#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"




#define  N  4000   /* number of mass bins for halo */



/* a number of tables */

static double *rt,*mt,*m2r,*r2m;
static double *rfinal;  
static double *mfinal;
static double *rhofinal,*rho2r;  

static int III;   /* a dummy variable */



static double r200_halo,a_halo,rho0_halo,rhot_halo,gam1_halo,p1_halo;






structure()
{
  double disk_angmomentum(void);
  double additional_mass_in_halo_cutoff(void);
  double fc(double);
  double gc(double);
  int i;
  double jhalo,jdisk,jd;
  double hnew,dh;
 




  R200=V200/(10*H0);
  M200=pow(V200,3)/(10*G*H0);
  RS=R200/CC;
  M_DISK=MD*M200;
  M_BULGE=MB*M200;
  M_GAS=M_DISK*GasFraction;


  M_TOTAL=M200 + additional_mass_in_halo_cutoff();   

  rho0_halo=M200/( 4*PI*(log(1+CC)-CC/(1+CC)) * RS*RS*RS);  

  M_HALO=M_TOTAL-M_DISK-M_BULGE;
  
  jhalo=LAMBDA*sqrt(G)*pow(M200,1.5)*sqrt(2*R200/fc(CC));
  jdisk=JD*jhalo;


  printf("fc=%g\n",fc(CC));

  halo_spinfactor=1.5*LAMBDA*sqrt(2*CC/fc(CC))*
                  pow(log(1+CC)-CC/(1+CC),1.5)/gc(CC);


  printf("halo_spinfactor %g\n",halo_spinfactor);
  printf("R200: %g\n",R200);
  printf("M200: %g\n",M200);
  printf("Total mass will be: %g\n",M_TOTAL);




  setup_massprofile();  /* sets up tabulated cumulative mass profile */


  H=sqrt(2.0)/2.0*LAMBDA/fc(CC) * R200;  /* first guess for disk scale length */

  Z0=DiskHeight*H;   /* sets disk thickness */
  A=BulgeSize*H;     /* sets bulge size */


  do
    {
      solve_mass_shells();  /* computes new bin-radii for mass shells */
      
      jd=disk_angmomentum(); /* computes disk momentum */
      
      hnew=jdisk/jd*H;

      dh=hnew-H;

      if(fabs(dh)>0.5*H)
	{
	  dh= 0.5*H *dh/fabs(dh);
	}
      else
	dh=dh*0.1;

      H=H+dh;

      printf("Jd/J=%g   hnew: %g  \n",jd/jhalo,H);


      Z0=DiskHeight*H;   /* sets disk thickness */
      A=BulgeSize*H;     /* sets bulge size */
    }
  while(fabs(dh)/H > 1e-5);





  prepare_cumlative_profile();  /* prepare cumulative profile of dark mass */
}




prepare_cumlative_profile()
{
  double mass_cumulative_disk(double);
  double mass_cumulative_bulge(double);
  int i;

  

  for(i=2,mfinal[1]=0;i<=N;i++)
    {
      mfinal[i]= rt[i]/rfinal[i]* mt[i] - mass_cumulative_disk(rfinal[i]) - mass_cumulative_bulge(rfinal[i]);
    }

  for(i=1;i<=(N-1);i++)
    {
      rhofinal[i]= (mfinal[i+1] - mfinal[i])/
	( 4.0*PI/3*( pow(rfinal[i+1],3) - pow(rfinal[i],3) ));


    }
  rhofinal[N]=0;

  spline(rfinal,mfinal,N,1e40,1e40,m2r);
  spline(mfinal,rfinal,N,1e40,1e40,r2m);
  spline(rfinal,rhofinal,N,1e40,1e40,rho2r);


  LL= rfinal[N];
}




double halo_mass(double r)
{
  double x;
  
  if(r>rfinal[N])
    x=mfinal[N];
  else
    splint(rfinal,mfinal,m2r,N,r,&x);

  return x;
}


double halo_q_to_r(double q)
{
  double m,x;

  m=mfinal[N]*q;

  splint(mfinal,rfinal,r2m,N,m,&x);

  return x;
}


double halo_rho(double r)
{
  double x;
  
  if(r>rfinal[N])
    x=0;
  else
    splint(rfinal,rhofinal,rho2r,N,r,&x);

  return x;
}









double disk_angmomentum(void)
{
  double jdisk_int(double);
  double dmin(double a,double b);

  return M_DISK*qromb(jdisk_int,0,dmin(30*H,R200));
}

double dmin(double a,double b)
{
  if(a<b)
    return a;
  else
    return b;
}



double jdisk_int(double x)
{
  double vc2,Sigma0,vc,y;
  double mass_cumulative_disk(double);
  double mass_total_dark(double radius);

  if(x>1e-20)
    {
      vc2= G * (mass_total_dark(x) + mass_cumulative_bulge(x)) /x;
    }
  else
    vc2=0;


  if(vc2<0)
    {
      vc2=0;
      printf("wwww\n");
      exit(0);
    }

  Sigma0=(M_DISK)/(2*PI*H*H);
  y=x/(2*H);

  if(y>1e-4) 
    vc2+= x * 2*PI*G*Sigma0*y*(bessi0(y)*bessk0(y)-bessi1(y)*bessk1(y));
    	
  vc=sqrt(vc2);

  return pow(x/H,2)*vc*exp(-x/H);
}






solve_mass_shells()
{
  int i;
  double zriddr(double (*func)(double), double x1, double x2, double xacc);
  double masszero(double);


  for(i=2;i<=(N-1);i++)
    {
      III=i;
      rfinal[i]=zriddr(masszero,0,rt[N],1e-6*rt[2]);
    }


  for(i=2,mfinal[1]=0;i<=N;i++)
    {
      mfinal[i]= rt[i]/rfinal[i]* mt[i] - mass_cumulative_disk(rfinal[i]) - mass_cumulative_bulge(rfinal[i]);
    }

  spline(rfinal,mfinal,N,1e40,1e40,m2r);
  spline(mfinal,rfinal,N,1e40,1e40,r2m);

}


double masszero(double rf)
{
  double mi,ri;
  double mass_cumulative_disk(double r);
  double mass_cumulative_bulge(double r);
  
  mi=mt[III];
  ri=rt[III];

  return mi*ri-rf*( (1- M200/mt[N]*(MD+MB))*mi + mass_cumulative_disk(rf) + mass_cumulative_bulge(rf));
}












double mass_total_dark(double radius)
{
  double x;

  if(radius>rfinal[N])
    x=mfinal[N];
  else
    splint(rfinal,mfinal,m2r,N,radius,&x);

  return x;
}




setup_massprofile()
{
  int i;
  double q,s,qq,f,f_,ds,r;


  rt=dvector(1,N);
  mt=dvector(1,N);
  mfinal=dvector(1,N);
  r2m=dvector(1,N);
  m2r=dvector(1,N);
  rfinal=dvector(1,N);
  rhofinal=dvector(1,N);
  rho2r=dvector(1,N);  



  for(i=2, mt[1]=rt[1]=0 ;i<=N;i++)
    {
      mt[i]=(i-1.01)*(M_TOTAL/(N-1));

      q=mt[i]/M_TOTAL;

      if(q< (M200/M_TOTAL))
	{
	  s=1.0;
	  qq=q*M_TOTAL/(4*PI*rho0_halo*RS*RS*RS);
	  do
	    {
	      f=log(1+s)-s/(1+s)-qq;
	      f_=s/(1+s)/(1+s);
	      ds=-f/f_;
	      if(fabs(ds)/s > 0.1)
		ds=s*0.1*ds/fabs(ds);
	      s+=ds;
	    }
	  while(fabs(ds/s)>1e-8);
	}
      else
	{
	  s=R200/RS;
	  do
	    {
	      f=gam1_halo*gammp(3+a_halo,s)-p1_halo-q*M_TOTAL+M200;
	      f_=gam1_halo*exp((a_halo+2)*log(s)-s-gammln(a_halo+3));
	      ds=-f/f_;
	      s+=ds;
	    }
	  while(fabs(ds/s)>1e-8);
	}
	
      r=s*RS;
      rt[i]=r;
    }

  for(i=1;i<=N;i++)
    rfinal[i]=rt[i];

  spline(rfinal,mt,N,1e40,1e40,m2r);
  spline(mt,rfinal,N,1e40,1e40,r2m);
}








double additional_mass_in_halo_cutoff(void)
{
  a_halo=CC-(1+3*CC)/(1+CC);
  rho0_halo=M200/( 4*PI*(log(1+CC)-CC/(1+CC)) * RS*RS*RS);
  rhot_halo=rho0_halo/(CC*(1+CC)*(1+CC));
  gam1_halo=4*PI*rhot_halo*R200*R200*R200*exp(CC+gammln(3+a_halo)-(3+a_halo)*log(CC));
  p1_halo=gam1_halo*gammp(3+a_halo,CC);
  
  return gam1_halo-p1_halo;
}





double fc(double c)
{
  return c*(0.5-0.5/pow(1+c,2)-log(1+c)/(1+c))/pow(log(1+c)-c/(1+c),2);
}

double gc(double c)
{
  double gc_int(double);
  return qromb(gc_int,0,c);
}

double gc_int(double x)
{
  return pow(log(1+x)-x/(1+x),0.5)*pow(x,1.5)/pow(1+x,2);
}





















write_cumulative_mass()
{
  FILE *fd;
  int i;

  fd=fopen("cummass.txt","w");

  for(i=1;i<=N;i++)
    {
      fprintf(fd,"%g %g %g %g %g \n", rt[i]/R200, rfinal[i]/R200, mt[i]/M200, 
	      mass_cumulative_disk(rfinal[i])/M200 , mt[i]/M200  );
    }
  fclose(fd);
}



