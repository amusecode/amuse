#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"




static double R,z;

compute_velocity_dispersions_disk()
{
  int i,j;
  double z,R;
  double rho;

  if(N_DISK==0 && N_GAS==0) return;

  printf("disk velocity dispersion field...\n"); fflush(stdout);  

  for(i=0;i<=RSIZE;i++)
    {
      printf("disk A, %d\n",i);
      
      for(j=0;j<=ZSIZE;j++)
	{

	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z[i][j] * comp_rho_disk(list_R[i],list_z[j]);
	}

      spline(xl,yl,ZSIZE+1, 1e40, 1e40, D2yl);

      for(j=ZSIZE - 1, VelDispRz_disk[i][ZSIZE]=0  ;j>=0;j--)
	{
	  VelDispRz_disk[i][j] =  VelDispRz_disk[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_disk[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
    }

  for(i=0;i<=RSIZE;i++)
    {
      printf("disk B, %d\n",i);

      for(j=0;j<=ZSIZE;j++)
	{
	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z_dR[i][j] * comp_rho_disk(list_RplusdR[i],list_z[j]);

	}

      spline(xl,yl,ZSIZE+1, 1e40,1e40,D2yl);
      
      for(j=ZSIZE - 1, VelDispRz_dR_disk[i][ZSIZE]=0  ;j>=0;j--)
	{
	  VelDispRz_dR_disk[i][j] =  VelDispRz_dR_disk[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_dR_disk[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
      
    }

  
  for(i=0;i<=RSIZE;i++)
    {
      for(j=0;j<=ZSIZE;j++)
	{

	  R=list_R[i];
	  z=list_z[j];
	  
	  rho = comp_rho_disk(R,z);

	  if(rho>0)
	    {
	      if(i>0)
		VelDispPhi_disk[i][j]=R/rho * (VelDispRz_dR_disk[i][j]-VelDispRz_disk[i][j])/(list_RplusdR[i]-list_R[i]);
	      else
		VelDispPhi_disk[i][j]=0;

	      VelDispRz_disk[i][j]/=rho;
	    }
	  else
	    VelDispRz_disk[i][j]=VelDispPhi_disk[i][j]=0;

	  VelVc2_disk[i][j]=R*Dphi_R[i][j];
  
	  VelDispPhi_disk[i][j]+=VelVc2_disk[i][j]+VelDispRz_disk[i][j];
	  

	  if(VelDispPhi_disk[i][j]>VelDispRz_disk[i][j]/epi_gamma2[i])
	    {
	      VelStreamPhi_disk[i][j]=sqrt(VelDispPhi_disk[i][j]-VelDispRz_disk[i][j]/epi_gamma2[i]);
	      VelDispPhi_disk[i][j] = VelDispRz_disk[i][j]/epi_gamma2[i];
	    }
	  else
	    {
	      VelStreamPhi_disk[i][j]=0;
	      VelDispPhi_disk[i][j] = VelDispRz_disk[i][j]/epi_gamma2[i];
	    }



	  /*** alternativ ******/

	  /*	  if(VelVc2_disk[i][j]>0)
	    VelStreamPhi_disk[i][j]=sqrt(VelVc2_disk[i][j]);
	  else
	    VelStreamPhi_disk[i][j]=0;

	    VelDispPhi_disk[i][j] = VelDispRz_disk[i][j]/epi_gamma2[i]; */

	  /************/


	  if(VelDispRz_disk[i][j]<0)
	    VelDispRz_disk[i][j]=0;
	  
	  if(VelDispPhi_disk[i][j]<0)
	    VelDispPhi_disk[i][j]=0;
	}
    }

  printf("done.\n"); fflush(stdout);  
  
  /*  test(); */

}







double intz_di(double k);
double intz_di_abs(double);

double intR_di(double k);
double intR_di_abs(double k);


double comp_Dphi_z_disk(double RR,double zz)
{
  double comp_Dphi_z_disk_sph(double RR,double zz);
  double comp_Dphi_z_disk_exact(double RR,double zz);

  if(sqrt(RR*RR+zz*zz)>10*H)
    return comp_Dphi_z_disk_sph(RR,zz);
  else
    return comp_Dphi_z_disk_exact(RR,zz);
 
}

double comp_Dphi_z_disk_sph(double RR,double zz)
{
  double m;
  double r;


  r=sqrt(RR*RR+zz*zz);

  m=M_DISK*( 1- (1+r/H)*exp(-r/H) );

  return G*zz/(r*r*r)*m;
}



double comp_Dphi_z_disk_exact(double RR,double zz)
{
  int i;
  double dphiz,dphiz2,e2;
  double Sigma0;
  double in1,in2,in3,bb;
  double deltaz,zpos;


  if(N_DISK==0 && N_GAS==0) return 0;

  if(fabs(zz)<4*Z0)
    {
      deltaz=(6.0*Z0)/NSHEETS;

      dphiz=0;
      
      for(i=0;i<NSHEETS;i++)
	{
	  zpos=-3.0*Z0 + (i+0.5)*Z0*6.0/NSHEETS;
	  
	  R=RR;
	  z=zz-zpos;
	  
	  Sigma0=(M_DISK)/(2*PI*H*H) * deltaz/(2*Z0) * pow(2/(exp(zpos/Z0)+exp(-zpos/Z0)),2);

	  in1=qromb(intz_di,0,2/H);  
	  
	  bb=2;
	  do
	    {
	      in2=qromb(intz_di,bb/H,(bb+2)/H);  
	      in3=qromb(intz_di_abs,bb/H,(bb+2)/H);
	      in1+=in2;
	      bb+=2;
	    }
	  while(fabs(in3/in1)>1e-2);
	  
	  dphiz += 2*PI*G*Sigma0*H*H*( in1 );
	}
      return dphiz;
    }
  else
    {
      R=RR;
      z=zz;
      Sigma0=(M_DISK)/(2*PI*H*H);

      in1=qromb(intz_di,0,2/H);  
      
      bb=2;
      do
	{
	  in2=qromb(intz_di,bb/H,(bb+2)/H);  
	  in3=qromb(intz_di_abs,bb/H,(bb+2)/H);
	  in1+=in2;
	  bb+=2;
	}
      while(fabs(in3/in1)>1e-2);
	  
      dphiz = 2*PI*G*Sigma0*H*H*( in1 );

      return dphiz;

    }
}







double comp_Dphi_R_disk(double RR,double zz)
{
  double comp_Dphi_R_disk_sph(double RR,double zz);
  double comp_Dphi_R_disk_exact(double RR,double zz);

  if(RR>0)
    {
      if(sqrt(RR*RR+zz*zz)>10*H)
	return comp_Dphi_R_disk_sph(RR,zz);
      else
	return comp_Dphi_R_disk_exact(RR,zz);
    }
  else
    return 0;
}





double comp_Dphi_R_disk_sph(double RR,double zz)
{
  double m;
  double r;


  r=sqrt(RR*RR+zz*zz);

  m=M_DISK*( 1- (1+r/H)*exp(-r/H) );

  return G*RR/(r*r*r)*m;
}



double comp_Dphi_R_disk_exact(double RR,double zz)
{
  int i;
  double dphiR,e2;
  double Sigma0;
  double in1,in2,in3,bb;
  double deltaz,zpos;


  if(N_DISK==0 && N_GAS==0) return 0;

  if(fabs(zz)<4*Z0)
    {
      deltaz=(6.0*Z0)/NSHEETS;

      dphiR=0;
      
      for(i=0;i<NSHEETS;i++)
	{
	  zpos=-3.0*Z0 + (i+0.5)*Z0*6.0/NSHEETS;
	  
	  R=RR;
	  z=zz-zpos;
	  
	  Sigma0=(M_DISK)/(2*PI*H*H) * deltaz/(2*Z0) * pow(2/(exp(zpos/Z0)+exp(-zpos/Z0)),2);

	  in1=qromb(intR_di,0,2/H);  
	  
	  bb=2;
	  do
	    {
	      in2=qromb(intR_di,bb/H,(bb+2)/H);  
	      in3=qromb(intR_di_abs,bb/H,(bb+2)/H);
	      in1+=in2;
	      bb+=2;
	    }
	  while(fabs(in3/in1)>1e-2);
	  
	  dphiR += 2*PI*G*Sigma0*H*H*( in1 );
	}
      return dphiR;
    }
  else
    {
      R=RR;
      z=zz;
      Sigma0=(M_DISK)/(2*PI*H*H);

      in1=qromb(intR_di,0,2/H);  
      
      bb=2;
      do
	{
	  in2=qromb(intR_di,bb/H,(bb+2)/H);  
	  in3=qromb(intR_di_abs,bb/H,(bb+2)/H);
	  in1+=in2;
	  bb+=2;
	}
      while(fabs(in3/in1)>1e-2);
	  
      dphiR = 2*PI*G*Sigma0*H*H*( in1 );

      return dphiR;

    }
}






double comp_Dphi_R_disk_razorthin(double RR,double zz)
{
  double Sigma0,y;
  double dphidR;

  if(RR>0)
    {

      Sigma0=(M_DISK)/(2*PI*H*H);
      y=RR/(2*H);

      if(y>1e-4) 
	dphidR = 2*PI*G*Sigma0*y*(bessi0(y)*bessk0(y)-bessi1(y)*bessk1(y));
      else
	dphidR =0;

      return dphidR;
    }
  else
    return 0;
}









double comp_rho_disk(double R,double z)
{
  double x;

  x=(M_DISK)/(4*PI*H*H*Z0)*exp(-R/H)*pow(2/(exp(z/Z0)+exp(-z/Z0)),2);
  
  if(fabs(z)>6*Z0)
    x=0;
  
  return x;
}



double intz_di(double k)
{
  if(z>0)
    return ( bessj0(k*R)*k*exp(-z*k)/pow(1+k*k*H*H,1.5));
  else
    return (-bessj0(k*R)*k*exp(z*k)/pow(1+k*k*H*H,1.5));
}


double intz_di_abs(double k)
{
  if(z>0)
    return fabs( bessj0(k*R)*k*exp(-z*k)/pow(1+k*k*H*H,1.5));
  else
    return fabs(-bessj0(k*R)*k*exp(z*k)/pow(1+k*k*H*H,1.5));
}






double intR_di(double k)
{
  if(z>=0)
    return bessj1(k*R)*k*exp(-z*k)/pow(1+k*k*H*H,1.5);
  else
    return bessj1(k*R)*k*exp( z*k)/pow(1+k*k*H*H,1.5);
}

double intR_di_abs(double k)
{
  if(z>=0)
    return fabs(bessj1(k*R)*k*exp(-z*k)/pow(1+k*k*H*H,1.5));
  else
    return fabs(bessj1(k*R)*k*exp( z*k)/pow(1+k*k*H*H,1.5));
}





double mass_cumulative_disk(double R)
{
  return M_DISK*(1-(1+R/H)*exp(-R/H));
}


test()
{
  FILE *fd;
  int i,j;
  
  fd=fopen("out.dat","w");
  printf("VelDispRz_disk[0][0]= %14.8g \n", VelDispRz_disk[0][0]);
  printf("VelDispRz_disk[0][1]= %14.8g \n", VelDispRz_disk[0][1]);
  printf("VelDispRz_disk[0][2]= %14.8g \n", VelDispRz_disk[0][2]);
  printf("VelDispRz_disk[1][0]= %14.8g \n", VelDispRz_disk[1][0]);
  printf("VelDispRz_disk[1][1]= %14.8g \n", VelDispRz_disk[1][1]);
  printf("VelDispRz_disk[1][2]= %14.8g \n", VelDispRz_disk[1][2]);
  printf("VelDispRz_disk[2][0]= %14.8g \n", VelDispRz_disk[2][0]);
  printf("VelDispRz_disk[2][1]= %14.8g \n", VelDispRz_disk[2][1]);
  printf("VelDispRz_disk[2][2]= %14.8g \n", VelDispRz_disk[2][2]);



  printf("RSIZE= %d   ZSIZE=%d\n", RSIZE, ZSIZE);

  fwrite(&VelDispRz_disk[0][0], (RSIZE+1)*(ZSIZE+1),  sizeof(double), fd);
  


  fclose(fd);
}
