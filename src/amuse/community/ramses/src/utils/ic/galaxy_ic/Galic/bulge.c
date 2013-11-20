#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"


static double R,z;


compute_velocity_dispersions_bulge()
{
  int i,j;
  double z,R;
  double rho;

  if(N_BULGE==0) return;


  printf("bulge velocity dispersion field...\n"); fflush(stdout);  

  for(i=0;i<=RSIZE;i++)
    {
      printf("bulge A, %d\n",i);
      
      for(j=0;j<=ZSIZE;j++)
	{
	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z[i][j] * comp_rho_bulge(list_R[i],list_z[j]);
	}

      spline(xl,yl,ZSIZE+1,1e40,1e40,D2yl);

      for(j=ZSIZE - 1, VelDispRz_bulge[i][ZSIZE]=0  ;j>=0;j--)
	{
	  VelDispRz_bulge[i][j] =  VelDispRz_bulge[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_bulge[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
      
    }


  for(i=0;i<=RSIZE;i++)
    {
      printf("bulge B, %d\n",i);

      for(j=0;j<=ZSIZE;j++)
	{
	  xl[j+1]=list_z[j];
	  yl[j+1]=Dphi_z_dR[i][j] * comp_rho_bulge(list_RplusdR[i],list_z[j]);

	}

      spline(xl,yl,ZSIZE+1,1e40,1e40,D2yl);
      
      for(j=ZSIZE - 1, VelDispRz_dR_bulge[i][ZSIZE]=0  ;j>=0;j--)
	{
	  VelDispRz_dR_bulge[i][j] =  VelDispRz_dR_bulge[i][j+1];
	  if(fabs(yl[j+2])>1e-100 && fabs(yl[j+1])>1e-100)
	    VelDispRz_dR_bulge[i][j]+=
	      qromb(splint_xl_yl_D2yl,list_z[j],list_z[j+1]);
	}
      
    }

  
  for(i=0;i<=RSIZE;i++)
    {
      for(j=0;j<=ZSIZE;j++)
	{

	  R=list_R[i];
	  z=list_z[j];
	  
	  rho = comp_rho_bulge(R,z);

	  if(rho>0)
	    {
	      if(i>0)
		VelDispPhi_bulge[i][j]=R/rho * (VelDispRz_dR_bulge[i][j]-VelDispRz_bulge[i][j])/(list_RplusdR[i]-list_R[i]);
	      else
		VelDispPhi_bulge[i][j]=0;

	      VelDispRz_bulge[i][j]/=rho;
	    }
	  else
	    VelDispRz_bulge[i][j]=VelDispPhi_bulge[i][j]=0;

	  
	  VelVc2_bulge[i][j]=R*Dphi_R[i][j];
	  
	  VelDispPhi_bulge[i][j]+=VelVc2_bulge[i][j]+VelDispRz_bulge[i][j];

	  VelStreamPhi_bulge[i][j]= 0 ;
  
	  VelDispPhi_bulge[i][j]-=VelStreamPhi_bulge[i][j]*VelStreamPhi_bulge[i][j];


	  if(VelDispRz_bulge[i][j]<0)
	    VelDispRz_bulge[i][j]=0;
	  
	  if(VelDispPhi_bulge[i][j]<0)
	    VelDispPhi_bulge[i][j]=0;
	}
    }

  
  printf("done.\n"); fflush(stdout);  

}






double comp_Dphi_z_bulge(double RR,double zz)
{
  double m;
  double r;

  r=sqrt(RR*RR+zz*zz);

  m=M_BULGE * r*r / pow(A+r, 2);

  if(r>0)
    return G*zz/(r*r*r)*m;
  else
    return 0;
}



double comp_Dphi_R_bulge(double RR,double zz)
{
  double m;
  double r;

  r=sqrt(RR*RR+zz*zz);

  m=M_BULGE * r*r / pow(A+r, 2);


  if(r>0)
    return G*RR/(r*r*r)*m;
  else
    return 0;
}





double comp_rho_bulge(double R,double z)
{
  double r;

  r=sqrt(R*R+z*z);

  if(r<1e-6*A)
    r=1e-6*A;

  return  M_BULGE/(2*PI) * A/r/pow(A+r, 3);
}



double mass_cumulative_bulge(double R)
{
  return M_BULGE*pow ( R/(A+R) , 2 );
}


