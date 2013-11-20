
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"





double epicyclic_kappa2(double R)
{
  double dR,dphi,dphi_;

  if(R>0)
    {
      dR=R*0.02;

      dphi=comp_Dphi_R_halo(R,0)+comp_Dphi_R_disk(R,0);

      dphi_=comp_Dphi_R_halo(R+dR,0)+comp_Dphi_R_disk(R+dR,0);

      return 3*dphi/R + (dphi_-dphi)/dR;
      
    }
  else return 0;
}



int plot_toomre_stability(FILE *fd)
{
  int dummy;
  double *Q;
  int i;
  double Sigma0;
  int count;


  for(i=0,count=0;i<=RSIZE;i++)
    if(list_R[i]<=6*H)
      count++;
  


  Sigma0=(M_DISK)/(2*PI*H*H);
  Q=dvector(0,RSIZE);


  for(i=0;i<count;i++)
    {
      Q[i]=Qstabilizefactor*sqrt(VelDispRz_disk[i][0])*sqrt(epi_kappa2[i])/(3.36*G*Sigma0*exp(-list_R[i]/H));
    }


  fprintf(fd,"\n%d\n",count);

  for(i=0;i<count;i++)
    {
      fprintf(fd,"%g ",list_R[i]);
      fprintf(fd,"%g\n",Q[i]);
    }
}







int plot_circular_speeds(FILE *fd)
{
  int i;
  double R;
  double RMAX;
  double Sigma0;
  double vcd,vch,vcb,vc2;
  double vdisp;
  double vrot;
  int ir;
  double ur;

#define POINTS 1000



  RMAX=R200; 

  Sigma0=(M_DISK)/(2*PI*H*H);
  /*
    RMAX=10*H;
  */


  //fprintf(fd,"%d\n",POINTS);
  fprintf(fd,"0.0 0.0\n");

  for(i=1;i<=POINTS;i++)
    {
      R=(RMAX/POINTS)*i;
      ir=(int)(log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR));
      ur=( log(R/LL*(pow(FR,RSIZE)-1)+1)/log(FR)) - ir;
      vdisp=sqrt(VelDispRz_disk[ir][0]*(1-ur)+VelDispRz_disk[ir+1][0]*(ur));
      vrot=VelStreamPhi_disk[ir][0]*(1-ur)+VelStreamPhi_disk[ir+1][0]*(ur);      
      vcd=sqrt(R*comp_Dphi_R_disk_razorthin(R,0));
      vch=sqrt(R*comp_Dphi_R_halo(R,0));
      vcb=sqrt(R*comp_Dphi_R_bulge(R,0));
      vc2=sqrt(vcd*vcd+vch*vch+vcb*vcb);

      R = R * UnitLength_in_cm * 1.0E3 / CM_PER_MPC;
      vc2 = vc2 * UnitVelocity_in_cm_per_s / 1.0E5;
	  
      fprintf(fd,"%f ",R);
      fprintf(fd,"%f\n",vc2);
      //fprintf(fd,"%f ",vch);
      //fprintf(fd,"%f ",vcd);
      //fprintf(fd,"%f ",vcb);
      //fprintf(fd,"%f ",vrot);
      //fprintf(fd,"%f\n",vdisp);

    }

#undef POINTS 
}
