#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



int main(int argc,char *argv[])
{
  int i;
  char filedir[256];
  char fname[256];
  FILE *fd;
  double convert_unit;
  
  /*******************************************/

  CC=     10.;       /* halo concentration   (NFW DM Halo only)   */
					 /* halo concentration      */
  V200=   150.;      /* circular velocity v_200 (in km/sec) (NFW DM Halo only) */
					 /* circular velocity v_200 (in km/sec) (NFW DM Halo + baryons) */
  LAMBDA= 0.04;      /* spin parameter          */
  MD=     0.04;      /* disk mass fraction      */
  MB=     0.004;     /* bulge mass fraction     */
  JD= MD;            /* disk spin fraction      */

  GasFraction= 0.2;  /* relative content of gas in the disk*/ 
  DiskHeight=  0.1;  /* thickness of disk in units of radial scale length */
  BulgeSize=   0.1;  /* bulge scale length in units of disk scale length  */

  N_HALO= 100000;    /* desired number of particles in dark halo */
  N_DISK= 100000;    /* desired number of collisionless particles in disk */
  N_GAS=  100000;    /* number of gas particles in disk */ 
  N_BULGE= 10000;    /* number of bulge particles */ 

  HI_GasMassFraction=    0.1;     /* in terms of the total gas mass */
  HI_GasDiskScaleLength= 1.0;    /* in terms of scale length of the disk */ 

  Qstabilizefactor=1.3;

  /**********************************************************/


  if(argc!=2)
    {
      fprintf(stderr,"\n\nwrong argument(s).  Specify an output directory.\n\n");
      exit(0);
    }
  strcpy(filedir, argv[1]);


  init_units();        /* set system of units */
  structure();         /* determine structure of halo, disk, and bulge */
  init();              /* allocate arrays */


  set_halo_positions();
  set_disk_positions();
  set_bulge_positions();
  set_gas_positions();

  
  compute_force_field();

  compute_velocity_dispersions_disk();
  compute_velocity_dispersions_halo();  
  compute_velocity_dispersions_bulge();  
 
  compute_local_escape_speed();

  set_halo_velocities();
  set_disk_velocities();    
  set_gas_velocities();
  set_bulge_velocities();


  
  save_particles(filedir);


  strcpy(fname,filedir);
  strcat(fname,"Vcirc.dat");
  if(fd=fopen(fname,"w"))
    {
      printf("writing circular velocity curve\n");
      plot_circular_speeds(fd);
      fclose(fd);
    }
  else
    {
      fprintf(stderr,"Can't open file '%s'.\n",fname);
      exit(0);
    }

  /* Parameter file */
  strcpy(fname,filedir);
  strcat(fname,"params.txt");
  if(fd=fopen(fname,"w"))
    {
      printf("writing ic generation parameters\n");// + Toomre's Q\n");
      write_header(fd);
      //plot_toomre_stability(fd);
      fclose(fd);
    }
  else
    {
      fprintf(stderr,"Can't open file '%s'.\n",fname);
      exit(0);
    }
  printf("done.\n");
  
  convert_unit = H * UnitLength_in_cm * 1.0E3 / CM_PER_MPC;
  printf("Disk scale length (kpc) : %g\n",convert_unit);
  convert_unit = R200 * UnitLength_in_cm * 1.0E3 / CM_PER_MPC;
  printf("R200 (kpc) : %g\n",convert_unit);
}



int write_header(FILE *fd)
{
  double convert_unit;
  
  fprintf(fd,"CC=%g\n",CC);
  convert_unit= V200 * UnitVelocity_in_cm_per_s / 1.0E5;
  fprintf(fd,"V200(km/s)=%g\n",convert_unit);
  fprintf(fd,"LAMBDA=%g\n",LAMBDA);
  fprintf(fd,"MD=%g\n",MD);
  fprintf(fd,"JD=%g\n",JD);
  fprintf(fd,"MB=%g\n",MB);
  fprintf(fd,"DiskHeight=%g\n",DiskHeight);
  fprintf(fd,"BulgeSize=%g\n",BulgeSize);
  convert_unit = R200 * UnitLength_in_cm * 1.0E3 / CM_PER_MPC;
  fprintf(fd,"\nR200=%g kpc\n",convert_unit);
  convert_unit = H * UnitLength_in_cm * 1.0E3 / CM_PER_MPC;
  fprintf(fd,"H=%g kpc\n\n",convert_unit);
}
