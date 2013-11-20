#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"


static double R,z;



int set_gas_velocities()
{
  int   i;
  double q,R,phi,theta;
  long  dum;
  int   iz,ir;
  double ur,uz;
  double vdisp_rz,vdisp_phi,vstream_phi;
  double vr,vphi;
  double vx,vy,vz;



  if(N_GAS==0) return;

  dum=drand48()*1e8;

  printf("set gas velocities..."); fflush(stdout);  
  
  for(i=1;i<=N_GAS;i++)
    {
      R=sqrt(xp_gas[i]*xp_gas[i] + yp_gas[i]*yp_gas[i]);
      z=zp_gas[i];


      ir=(int)( log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR));
      ur=( log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR)) - ir;

      iz=(int)( log(fabs(z)/LL*(pow(FZ,ZSIZE)-1)+1)/log(FZ));
      uz=( log(fabs(z)/LL*(pow(FZ,ZSIZE)-1)+1)/log(FZ)) - iz;


   
      vdisp_rz= VelDispRz_disk[ir][iz]*(1-ur)*(1-uz)
               +VelDispRz_disk[ir+1][iz]*(ur)*(1-uz)
   	       +VelDispRz_disk[ir][iz+1]*(1-ur)*(uz) 
               +VelDispRz_disk[ir+1][iz+1]*(ur)*(uz);

      vdisp_phi=VelDispPhi_disk[ir][iz]*(1-ur)*(1-uz)
	       +VelDispPhi_disk[ir+1][iz]*(ur)*(1-uz)
	       +VelDispPhi_disk[ir][iz+1]*(1-ur)*(uz) 
               +VelDispPhi_disk[ir+1][iz+1]*(ur)*(uz);

      vstream_phi=VelStreamPhi_disk[ir][iz]*(1-ur)*(1-uz)
	       +VelStreamPhi_disk[ir+1][iz]*(ur)*(1-uz)
	       +VelStreamPhi_disk[ir][iz+1]*(1-ur)*(uz) 
	       +VelStreamPhi_disk[ir+1][iz+1]*(ur)*(uz);


      if(vdisp_rz<0)
	{
	  printf("in gas: vdisp_rz:%g   %g %g %d %d \n",vdisp_rz,ur,uz,ir,iz);
	  vdisp_rz=-vdisp_rz;
	}

      if(vdisp_phi<0)
	{
	  printf("in gas: vdisp_phi:%g  %g %g %d %d\n",vdisp_phi,ur,uz,ir,iz);
	  
	  vdisp_phi=0;
	}


      vr=0;
      vz=0;

      vphi=vstream_phi;

      vx=vr*xp_gas[i]/R - vphi*yp_gas[i]/R;
      vy=vr*yp_gas[i]/R + vphi*xp_gas[i]/R;

      vxp_gas[i]=vx;
      vyp_gas[i]=vy;
      vzp_gas[i]=vz;

      u_gas[i]=vdisp_rz/(GAMMA-1);

      
      if((vx*vx+vy*vy+vz*vz)>0.95*vmax2_gas[i])
	{
	  /* printf("%d Gas velocity rejected\n",i); */
	  i--;
	}
    }

  printf("done.\n"); fflush(stdout);
}





int set_gas_positions()
{
  int   i,countr,countz;
  double q,R,f,f_,Rold,phi,theta;

  int n_disk,n_HI;

  if(N_GAS==0) return;

  srand48(22277);
  
  if(HI_GasMassFraction>0)
    {
      n_disk= (1-HI_GasMassFraction)*N_GAS;
      n_HI= N_GAS - n_disk;
    }
  else
    {
      n_disk= N_GAS;
      n_HI= 0;
    }


  for(i=1,countr=countz=0;i<=n_disk;)
    {
      q=drand48();

      zp_gas[i]=Z0/2*log(q/(1-q));

      q=drand48();
      
      R=1.0;
      do
	{
	  f=(1+R)*exp(-R)+q-1;
	  f_=-R*exp(-R);
	  
	  Rold=R;
	  R=R-f/f_;
	}
      while(fabs(R-Rold)/R> 1e-6);

      R*=H;
      
      phi=drand48()*PI*2;
	  
      xp_gas[i]=R*cos(phi);
      yp_gas[i]=R*sin(phi);
      

      if(R>LL || fabs(zp_gas[i])>LL)
	countr++;
      else 
	i++;
    }


  for(i=1+n_disk;i<=N_GAS;)
    {
      q=drand48();

      zp_gas[i]=Z0/2*log(q/(1-q));

      q=drand48();
      
      R= H * sqrt(q) * HI_GasDiskScaleLength;
      
      phi=drand48()*PI*2;
	  
      xp_gas[i]=R*cos(phi);
      yp_gas[i]=R*sin(phi);
      

      if(R>LL || fabs(zp_gas[i])>LL)
	countr++;
      else 
	i++;
    }


  for(i=1;i<=N_GAS;i++)
    mp_gas[i]=M_GAS/N_GAS;
}




