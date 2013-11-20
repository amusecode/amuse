#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"
#include "prototypes.h"
#include "globvars.h"


#ifdef T3E
  typedef short int int4byte;   /* Note: int has 8 Bytes on the T3E ! */
#else
  typedef int int4byte;
#endif





struct io_header_1
{
  int4byte npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int4byte flag_sfr;
  int4byte flag_feedback;
  int4byte npartTotal[6];
  int4byte flag_cooling;
  int4byte num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;







void save_particles(char *fdir)
{
  /* ************************************************************************************************************
   * Save particle information into the RAMSES (v>=3.08, merger patch) ic_part format :
   * x (in kpc)  |  y (in kpc)  |  z (in kpc)  |  vx (in km/s)  | vy (in km/s)  |  vz (in km/s)  |  m (in GMsun)
   *    -----         -----         -----           -----           -----            -----           ------
   **************************************************************************************************************
  */
  FILE *fd_ic_part, *fd_info;
  char fname[256];
  int i,d,ntot;
  float xyz[3];
  double t, vfact_km_s;
  double mtot,mtot_sum,lb;
  int4byte blklen;
#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, fd);

  strcpy(fname,fdir);
  strcat(fname,"ic_part");
  
  if(!(fd_ic_part=fopen(fname,"w")))
    {
      printf("error opening file '%s'\n",fname);
      exit(0);
    }
  printf("saving initial conditions to file '%s'\n",fname);

  strcpy(fname,fdir);
  strcat(fname,"info.txt");
  
  if(!(fd_info=fopen(fname,"w")))
    {
      printf("error opening file '%s'\n",fname);
      exit(0);
    }
  printf("saving ic information to file '%s'\n\n",fname);
  fprintf(fd_info,"h                = 1.0\n");
  lb = R200 * 4.0;
  fprintf(fd_info,"Lbox             = %f\n",lb);
  mtot_sum = 0.0;
  vfact_km_s = UnitVelocity_in_cm_per_s / 1.0E5;


  /* Dark matter particles */
  printf("Number of particles in the halo %d\n",N_HALO);
  fprintf(fd_info," -Ndark_matter   = %d\n",N_HALO);
  mtot = 0.0;
  for(i=1;i<=N_HALO;i++)
    {
      fprintf(fd_ic_part," %g",xp_halo[i]);
      fprintf(fd_ic_part," %g",yp_halo[i]);
      fprintf(fd_ic_part," %g",zp_halo[i]);
      fprintf(fd_ic_part," %g",vxp_halo[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vyp_halo[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vzp_halo[i]*vfact_km_s);
      fprintf(fd_ic_part," %g\n",mp_halo[i]);
	  mtot +=mp_halo[i];
	}
  fprintf(fd_info," -Mass_dm        = %f\n",mtot);
  mtot_sum += mtot;
 
  /* Bulge star particles */
  printf("Number of particles in the bulge %d\n",N_BULGE);
  fprintf(fd_info," -Nstars_bulge   = %d\n",N_BULGE);
  mtot = 0.0;
  for(i=1;i<=N_BULGE;i++)
    {
      fprintf(fd_ic_part," %g",xp_bulge[i]);
      fprintf(fd_ic_part," %g",yp_bulge[i]);
      fprintf(fd_ic_part," %g",zp_bulge[i]);
      fprintf(fd_ic_part," %g",vxp_bulge[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vyp_bulge[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vzp_bulge[i]*vfact_km_s);
      fprintf(fd_ic_part," %g\n",mp_bulge[i]);
	  mtot +=mp_bulge[i];
    }
  fprintf(fd_info," -Mass_starsbulge= %f\n",mtot);
  mtot_sum += mtot;
  
  
  /* Disk star particles */
  printf("Number of particles in the disk %d\n",N_DISK);
  fprintf(fd_info," -Nstars_disk    = %d\n",N_DISK);
  mtot = 0.0;
  for(i=1;i<=N_DISK;i++)
    {
      fprintf(fd_ic_part," %g",xp_disk[i]);
      fprintf(fd_ic_part," %g",yp_disk[i]);
      fprintf(fd_ic_part," %g",zp_disk[i]);
      fprintf(fd_ic_part," %g",vxp_disk[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vyp_disk[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vzp_disk[i]*vfact_km_s);
      fprintf(fd_ic_part," %g\n",mp_disk[i]);
	  mtot +=mp_disk[i];
    }      
  fprintf(fd_info," -Mass_starsdisk = %f\n",mtot);
  mtot_sum += mtot;
   
  
  /* Gas particles */
/*   printf("Number of particles in the gas %d (not written)\n",N_GAS); */
  printf("Number of particles in the gas %d\n",N_GAS);
  fprintf(fd_info," -Ngaz           = %d\n", N_GAS);
  mtot = 0.0;
  for(i=1;i<=N_GAS;i++)
    {
      fprintf(fd_ic_part," %g",xp_gas[i]);
      fprintf(fd_ic_part," %g",yp_gas[i]);
      fprintf(fd_ic_part," %g",zp_gas[i]);
      fprintf(fd_ic_part," %g",vxp_gas[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vyp_gas[i]*vfact_km_s);
      fprintf(fd_ic_part," %g",vzp_gas[i]*vfact_km_s);
      fprintf(fd_ic_part," %g\n",mp_gas[i]);
	  mtot += mp_gas[i];
    }
  fprintf(fd_info," -Mass_gaz       = %f\n", mtot);
  printf("Total mass of gaseous disk :  %f\n",mtot);
  mtot_sum += mtot;

 
  /* New formed stars */
  fprintf(fd_info," -Nstar_form     = 0\n");
  fprintf(fd_info," -Mass_sf        = 0.0\n");
  
  /* Black holes */
  fprintf(fd_info," -Nblack_holes   = 0\n");
  fprintf(fd_info," -Mass_bh        = 0.0\n");

  
  ntot = N_BULGE + N_DISK + N_HALO + N_GAS;
  /* Total number of written particles */
  fprintf(fd_info," -Ntot           = %d\n",ntot);
  fprintf(fd_info," -Mass_tot       = %f\n",mtot_sum);
  
  
  fclose(fd_ic_part);
  fclose(fd_info);
}


