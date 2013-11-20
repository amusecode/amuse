#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"

 




compute_force_field()
{
  int i,j;
  double k2,dphi_R_dr;


  printf("Start computing force field.\n");

  for(i=1;i<=RSIZE;i++)
    {
      printf("force field %d(%d)\n",i,RSIZE);

      for(j=0;j<=ZSIZE;j++)
	{
	  if(j==0)
	    {
	      Dphi_z[i][j] = 0;
	      Dphi_z_dR[i][j] =0;
	    }
	  else
	    {
	      Dphi_z[i][j]    = comp_Dphi_z( list_R[i], list_z[j] );
	      Dphi_z_dR[i][j] = comp_Dphi_z( list_RplusdR[i] , list_z[j] );
	    }
	  if(i==0)
	    Dphi_R[i][j] = 0;
	  else
	    Dphi_R[i][j] = comp_Dphi_R( list_R[i], list_z[j]);
	}

    }


  
  for(j=0;j<=ZSIZE;j++)
    {
      Dphi_z[0][j]    = Dphi_z[1][j];
      Dphi_z_dR[0][j] = Dphi_z_dR[1][j];
    }


  for(i=1,epi_gamma2[0]=1;i<=RSIZE;i++)
    {
      dphi_R_dr = comp_Dphi_R( list_RplusdR[i], 0);
      
      k2=3/list_R[i] *  Dphi_R[i][0] + (dphi_R_dr - Dphi_R[i][0]) /(list_RplusdR[i]- list_R[i]);
      
      epi_gamma2[i]=4/list_R[i]*Dphi_R[i][0]/k2;
      epi_kappa2[i]=k2;
    }

  epi_kappa2[0]=epi_kappa2[1];

  printf("Force field finished.\n");
}






double comp_Dphi_z(double R,double z)
{
  return comp_Dphi_z_halo(R,z)
        +comp_Dphi_z_bulge(R,z) 
        +comp_Dphi_z_disk(R,z) ; 
}

double comp_Dphi_R(double R,double z)
{
  return comp_Dphi_R_halo(R,z)
        +comp_Dphi_R_bulge(R,z) 
        +comp_Dphi_R_disk(R,z);
}


