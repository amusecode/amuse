#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



static double R,z;



double set_halo_velocities(void)
{
  int   i;
  double q,R,phi,theta;
  long  dum;
  int   iz,ir;
  double ur,uz;
  double vdisp_rz,vdisp_phi,vstream_phi;
  double vr,vphi;
  double vx,vy,vz;


  if(N_HALO==0) return 0;

  dum=drand48()*1e8;

  printf("set halo velocities..."); fflush(stdout);
  
  for(i=1;i<=N_HALO;i++)
    {
      R=sqrt(xp_halo[i]*xp_halo[i] + yp_halo[i]*yp_halo[i]);
      z=zp_halo[i];


      ir=(int)( log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR));
      ur=( log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR)) - ir;

      iz=(int)( log(fabs(z)/LL*(pow(FZ,ZSIZE)-1)+1)/log(FZ));
      uz=( log(fabs(z)/LL*(pow(FZ,ZSIZE)-1)+1)/log(FZ)) - iz;

   
      vdisp_rz= VelDispRz_halo[ir][iz]*(1-ur)*(1-uz)
               +VelDispRz_halo[ir+1][iz]*(ur)*(1-uz)
   	       +VelDispRz_halo[ir][iz+1]*(1-ur)*(uz) 
               +VelDispRz_halo[ir+1][iz+1]*(ur)*(uz);

      vdisp_phi=VelDispPhi_halo[ir][iz]*(1-ur)*(1-uz)
	       +VelDispPhi_halo[ir+1][iz]*(ur)*(1-uz)
	       +VelDispPhi_halo[ir][iz+1]*(1-ur)*(uz) 
               +VelDispPhi_halo[ir+1][iz+1]*(ur)*(uz);

      vstream_phi=VelStreamPhi_halo[ir][iz]*(1-ur)*(1-uz)
	       +VelStreamPhi_halo[ir+1][iz]*(ur)*(1-uz)
	       +VelStreamPhi_halo[ir][iz+1]*(1-ur)*(uz) 
               +VelStreamPhi_halo[ir+1][iz+1]*(ur)*(uz);
      
      if(vdisp_rz<0)
	{
	  printf("in halo: vdisp_rz:%g   %g %g %d %d \n",vdisp_rz,ur,uz,ir,iz);
	  vdisp_rz=-vdisp_rz;
	}
      if(vdisp_phi<0)
	{
	  printf("in halo: vdisp_phi:%g  %g %g %d %d\n",vdisp_phi,ur,uz,ir,iz);
	  
	  vdisp_phi=-vdisp_phi;
	}

      vr=gasdev(&dum)*sqrt(vdisp_rz);
      vz=gasdev(&dum)*sqrt(vdisp_rz);

      vphi=vstream_phi + gasdev(&dum)*sqrt(vdisp_phi);


      vx=vr*xp_halo[i]/R - vphi*yp_halo[i]/R;
      vy=vr*yp_halo[i]/R + vphi*xp_halo[i]/R;
      
      vxp_halo[i]=vx;
      vyp_halo[i]=vy;
      vzp_halo[i]=vz;
      
      if((vx*vx+vy*vy+vz*vz)>0.95*vmax2_halo[i])
	{
	  /* printf("%d Halo velocity rejected\n",i);*/
	  i--;
	}
    }
  
  printf("done.\n"); fflush(stdout);
  
  return 0;
}






double set_halo_positions(void)
{
  int    i,countr,countz;
  double s,ds,q,qq,R,dR,phi,theta,f,f_;
  double halo_q_to_r(double q);

  if(N_HALO==0) return 0;
  
  srand48(22);
  
  for(i=1,countr=countz=0;i<=N_HALO;)
    {
      q=drand48();
      
      R=halo_q_to_r(q);

      if(R>LL)
	{
	  continue;
	}
      
      phi=drand48()*PI*2;
      theta=acos(drand48()*2-1);
      
      xp_halo[i]=R*sin(theta)*cos(phi);
      yp_halo[i]=R*sin(theta)*sin(phi);
      zp_halo[i]=R*cos(theta);

      i++;
    }

  for(i=1;i<=N_HALO;i++)
    mp_halo[i]=M_HALO/N_HALO;

  return 0;
}















