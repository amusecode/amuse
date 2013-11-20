


/***********  INPUT PARAMETERS *********/

 double  CC;      /* halo concentration */
 double  V200;    /* circular velocity v_200 */
 double  LAMBDA;    /* spin parameter  */
 double  MD;        /* disk mass fraction */
 double  JD;        /* disk spin fraction */
 double  MB;        /* bulge mass fraction */
 double  GasFraction;  
 double  DiskHeight; 
 double  BulgeSize;
 double  HI_GasMassFraction;    /* in terms of the total gas mass */
 double  HI_GasDiskScaleLength;  /* in terms of scale length of the disk */ 

 double  Qstabilizefactor;

 int     N_HALO;    /* desired number of particles in halo */
 int     N_DISK;    /* desired number of collsionless particles in disk */
 int     N_GAS;     /* number of gas particles in stellar disk */ 
 int     N_BULGE;   /* number of gas particles in stellar disk */ 

/*********************************************/







 double  M200;     /* virial mass */
 double  M_TOTAL;  /* total mass */

 double  RS;       /* scale radius for halo */
 double  R200;     /* virial radius */
 double  H;        /* disk scale length */
 double  Z0;       /* disk thickness */
 double  A;        /* bulge scale radius */


 double  M_HALO;   /* total dark mass */
 double  M_DISK;   /* mass of stellar disk (collisionless part) */
 double  M_GAS;    /* gas mass in disk */
 double  M_BULGE;  /* mass of bulge */

 double  halo_spinfactor;  /* computed streamin of dark matter */






 double G;            /* gravitational constant */
 double H0;           /* Hubble constant */
 double UnitTime_in_s;
 double UnitMass_in_g;
 double UnitLength_in_cm;
 double UnitVelocity_in_cm_per_s;
 double UnitTime_in_Megayears;



/* particle data */

 double    *vmax2_halo,*vmax2_disk,*vmax2_bulge,*vmax2_gas;

 double    *xp_halo,*yp_halo,*zp_halo,*mp_halo;
 double    *xp_disk,*yp_disk,*zp_disk,*mp_disk;
 double    *xp_bulge,*yp_bulge,*zp_bulge,*mp_bulge;
 double    *xp_gas,*yp_gas,*zp_gas,*mp_gas,*u_gas;

 double    *vxp_halo,*vyp_halo,*vzp_halo;
 double    *vxp_disk,*vyp_disk,*vzp_disk;
 double    *vxp_bulge,*vyp_bulge,*vzp_bulge;
 double    *vxp_gas,*vyp_gas,*vzp_gas;











double  LL,HR; /* LL = extension of fields in R and z.
		  HR = extension of high resolution region in z */
double  dR;    /* delta R */





 double **Dphi_z,**Dphi_R,**Dphi_z_dR;  /* derivatives of total potential */
 double *epi_gamma2,*epi_kappa2;   /* epicycle gamma^2  */ 



/* halo velocity fields */

 double **VelDispRz_halo;
 double **VelDispPhi_halo;
 double **VelVc2_halo;
 double **VelStreamPhi_halo;
 double **VelDispRz_dR_halo;


/* bulge velocity fields */

 double **VelDispRz_bulge;
 double **VelDispPhi_bulge;
 double **VelVc2_bulge;
 double **VelStreamPhi_bulge;
 double **VelDispRz_dR_bulge;


/* disk velocity fields */

 double **VelDispRz_disk;
 double **VelDispPhi_disk;
 double **VelVc2_disk;
 double **VelStreamPhi_disk;
 double **VelDispRz_dR_disk;




/* auxiliary field */

 double *xl,*yl,*D2yl;
 double *list_z,*list_R,*list_RplusdR;






