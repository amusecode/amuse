#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"




compute_velocity_dispersions_halo()
{
  int i,j;
  double z,R;
  double rho;
 

  printf("halo velocity dispersion field...\n"); fflush(stdout);  

  for(i=0;i<=RSIZE;i++)
    {
      printf("halo A, %d\n",i);
      
      for(j=0;j<=ZSIZE;j++)
	{
	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z[i][j] * comp_rho_halo(list_R[i],list_z[j]);
	}

      spline(xl,yl,ZSIZE+1,1e40,1e40,D2yl);


      for(j=ZSIZE - 1, VelDispRz_halo[i][ZSIZE]=0  ;j>=0;j--)
	{

	  VelDispRz_halo[i][j] =  VelDispRz_halo[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_halo[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
    }

  for(i=0;i<=RSIZE;i++)
    {
      printf("halo B, %d\n",i);

      for(j=0;j<=ZSIZE;j++)
	{
	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z_dR[i][j] * comp_rho_halo(list_RplusdR[i],list_z[j]);

	}

      spline(xl,yl,ZSIZE+1,1e40,1e40,D2yl);
      
      for(j=ZSIZE - 1, VelDispRz_dR_halo[i][ZSIZE]=0  ;j>=0;j--)
	{
	  VelDispRz_dR_halo[i][j] =  VelDispRz_dR_halo[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_dR_halo[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
      
    }

  
  for(i=0;i<=RSIZE;i++)
    {
      for(j=0;j<=ZSIZE;j++)
	{
	  R=list_R[i];
	  z=list_z[j];
	  
	  rho = comp_rho_halo(R,z);

	  if(rho>0)
	    {
	      if(i>0)
		VelDispPhi_halo[i][j]=R/rho * (VelDispRz_dR_halo[i][j]-VelDispRz_halo[i][j])/(list_RplusdR[i]-list_R[i]);
	      else
		VelDispPhi_halo[i][j]=0;

	      VelDispRz_halo[i][j]/=rho;
	    }
	  else
	    VelDispRz_halo[i][j]=VelDispPhi_halo[i][j]=0;

	  
	  VelVc2_halo[i][j]=R*Dphi_R[i][j];
	  
	  VelDispPhi_halo[i][j]+=VelVc2_halo[i][j]+VelDispRz_halo[i][j];

	  VelStreamPhi_halo[i][j]=halo_spinfactor * sqrt(VelVc2_halo[i][j]);
  
	  VelDispPhi_halo[i][j]-=VelStreamPhi_halo[i][j]*VelStreamPhi_halo[i][j];


	  if(VelDispRz_halo[i][j]<0)
	    VelDispRz_halo[i][j]=0;

	  if(VelDispPhi_halo[i][j]<0)
	    VelDispPhi_halo[i][j]=0;
	}
    }

  
  printf("done.\n"); fflush(stdout);  
}


double comp_Dphi_z_halo(double R,double z)
{
  double r,M_r;
  double halo_mass(double r);  

  r=sqrt(R*R+z*z);
  
  M_r=halo_mass(r);

  if(r>0)
    return G*z/(r*r*r)*M_r;
  else
    return 0;
}


double comp_Dphi_R_halo(double R,double z)
{
  double r,M_r;
  double halo_mass(double r);  


  r=sqrt(R*R+z*z);

  M_r=halo_mass(r);

  if(r>0)
    return G*R/(r*r*r)*M_r;
  else
    return 0;
}


double comp_rho_halo(double R,double z)
{
  double r,rr;
  double halo_rho(double);
 

  r=sqrt(R*R+z*z);

  return halo_rho(r);
}

