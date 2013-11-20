
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"





init_units()
{

  UnitMass_in_g= 1.0e9 * SOLAR_MASS;                                      /* 1.0e9 M_sun    */
  UnitLength_in_cm= CM_PER_MPC / 1000;                                    /* 1 kpc          */
  UnitTime_in_s= sqrt(pow(UnitLength_in_cm,3)/(GRAVITY * UnitMass_in_g)); /* ~= 14.922 Myr  */
  UnitVelocity_in_cm_per_s= UnitLength_in_cm / UnitTime_in_s;             /* ~= 65.6 km/sec */

  //UnitMass_in_g= 2.3262e5 * SOLAR_MASS;    /* ... */
  //UnitLength_in_cm= CM_PER_MPC/1000;   /* 1 kpc */
  //UnitVelocity_in_cm_per_s=100000;        /* 1 km/sec */

  //UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  H0=HUBBLE * 100*1e5 / CM_PER_MPC / UnitVelocity_in_cm_per_s * UnitLength_in_cm;
  /* V200 : units conversion (km/s => user units) */
  V200=V200 * 1.0E5 / UnitVelocity_in_cm_per_s;
  G=GRAVITY/pow(UnitLength_in_cm,3)*UnitMass_in_g*pow(UnitTime_in_s,2);

  UnitTime_in_Megayears=UnitTime_in_s/SEC_PER_MEGAYEAR;

  printf("G = %g, H0 = %g, unit_time = %g Myr \n",G,H0,UnitTime_in_Megayears);
}







init()
{
  int i;
  double baselen;


  if(N_HALO>0)
    {
      xp_halo=vector(1,N_HALO);
      yp_halo=vector(1,N_HALO);
      zp_halo=vector(1,N_HALO);
      mp_halo=vector(1,N_HALO);
      vxp_halo=vector(1,N_HALO);
      vyp_halo=vector(1,N_HALO);
      vzp_halo=vector(1,N_HALO);
      vmax2_halo=vector(1,N_HALO);
    }

  if(N_DISK>0)
    {
      xp_disk=vector(1,N_DISK);
      yp_disk=vector(1,N_DISK);
      zp_disk=vector(1,N_DISK);
      mp_disk=vector(1,N_DISK);
      vxp_disk=vector(1,N_DISK);
      vyp_disk=vector(1,N_DISK);
      vzp_disk=vector(1,N_DISK);
      vmax2_disk=vector(1,N_DISK);      
    }

  if(N_BULGE>0)
    {
      xp_bulge=vector(1,N_BULGE);
      yp_bulge=vector(1,N_BULGE);
      zp_bulge=vector(1,N_BULGE);
      mp_bulge=vector(1,N_BULGE);
      vxp_bulge=vector(1,N_BULGE);
      vyp_bulge=vector(1,N_BULGE);
      vzp_bulge=vector(1,N_BULGE);
      vmax2_bulge=vector(1,N_BULGE);      
    }

  if(N_GAS>0)
    {
      xp_gas=vector(1,N_GAS);
      yp_gas=vector(1,N_GAS);
      zp_gas=vector(1,N_GAS);
      mp_gas=vector(1,N_GAS);
      vxp_gas=vector(1,N_GAS);
      vyp_gas=vector(1,N_GAS);
      vzp_gas=vector(1,N_GAS);
      u_gas=vector(1,N_GAS);
      vmax2_gas=vector(1,N_GAS);    
    }



  Dphi_z    = matrix(0,RSIZE,0,ZSIZE);
  Dphi_R    = matrix(0,RSIZE,0,ZSIZE);
  Dphi_z_dR = matrix(0,RSIZE,0,ZSIZE);

  VelDispRz_halo=matrix(0,RSIZE,0,ZSIZE);
  VelDispRz_dR_halo=matrix(0,RSIZE,0,ZSIZE);
  VelDispPhi_halo=matrix(0,RSIZE,0,ZSIZE);
  VelVc2_halo=matrix(0,RSIZE,0,ZSIZE);
  VelStreamPhi_halo=matrix(0,RSIZE,0,ZSIZE);

  VelDispRz_disk=matrix(0,RSIZE,0,ZSIZE);
  VelDispRz_dR_disk=matrix(0,RSIZE,0,ZSIZE);
  VelDispPhi_disk=matrix(0,RSIZE,0,ZSIZE);
  VelVc2_disk=matrix(0,RSIZE,0,ZSIZE);
  VelStreamPhi_disk=matrix(0,RSIZE,0,ZSIZE);


  VelDispRz_dR_bulge=matrix(0,RSIZE,0,ZSIZE);
  VelDispRz_bulge=matrix(0,RSIZE,0,ZSIZE);
  VelDispPhi_bulge=matrix(0,RSIZE,0,ZSIZE);
  VelVc2_bulge=matrix(0,RSIZE,0,ZSIZE);
  VelStreamPhi_bulge=matrix(0,RSIZE,0,ZSIZE);








  list_R=vector(0,RSIZE);
  list_RplusdR=vector(0,RSIZE);
  list_z=vector(0,ZSIZE);


  baselen= LL / ( (pow(FR,RSIZE)-1)/(FR-1));

  for(i=1,list_R[0]=0;i<=RSIZE;i++)
    {
      list_R[i]=list_R[i-1] + baselen;
      baselen*=FR;
    }

  for(i=1,list_RplusdR[0]=0;i<=RSIZE;i++)
    {
      list_RplusdR[i]= list_R[i] +  (list_R[i] - list_R[i-1])*dRfac;
    }





  baselen= LL / ( (pow(FZ,ZSIZE)-1)/(FZ-1));

  for(i=1,list_z[0]=0;i<=ZSIZE;i++)
    {
      list_z[i]=list_z[i-1] + baselen;
      baselen*=FZ;
    }



  epi_gamma2=dvector(0,RSIZE);   /* epicycle gamma^2  */ 
  epi_kappa2=dvector(0,RSIZE); 


  xl=vector(1,ZSIZE+1);
  yl=vector(1,ZSIZE+1);
  D2yl=vector(1,ZSIZE+1);
}
