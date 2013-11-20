
#define  GAMMA (5.0/3)
#define  PI  3.1415926


#define  GRAVITY     6.67428e-8
#define  SOLAR_MASS  1.9889e33
#define  CM_PER_MPC  3.085677e24
#define  SEC_PER_MEGAYEAR   3.15576e13

#define  HUBBLE      0.73     /* Hubble constant in 100km/sec/Mpc */  




/*** integration parameters ***/

#define FR 1.5
#define FZ 1.5

#define dRfac 0.1     /* delta R factor*/

#define  RSIZE 50  /* 50 */   /* size of force field array */
#define  ZSIZE 50  /* 50 */

#define NSHEETS 100  /* 100 */





/***********  INPUT PARAMETERS *********/

extern double  CC;      /* halo concentration */
extern double  V200;    /* circular velocity v_200 */
extern double  LAMBDA;    /* spin parameter  */
extern double  MD;        /* disk mass fraction */
extern double  JD;        /* disk spin fraction */
extern double  MB;        /* bulge mass fraction */
extern double  GasFraction;  
extern double  DiskHeight; 
extern double  BulgeSize;
extern double  HI_GasMassFraction;    /* in terms of the total gas mass */
extern double  HI_GasDiskScaleLength;  /* in terms of scale length of the disk */ 

extern double  Qstabilizefactor;

extern int     N_HALO;    /* desired number of particles in halo */
extern int     N_DISK;    /* desired number of collsionless particles in disk */
extern int     N_GAS;     /* number of gas particles in stellar disk */ 
extern int     N_BULGE;   /* number of gas particles in stellar disk */ 

/*********************************************/







extern double  M200;     /* virial mass */
extern double  M_TOTAL;  /* total mass */

extern double  RS;       /* scale radius for halo */
extern double  R200;     /* virial radius */
extern double  H;        /* disk scale length */
extern double  Z0;       /* disk thickness */
extern double  A;        /* bulge scale radius */


extern double  M_HALO;   /* total dark mass */
extern double  M_DISK;   /* mass of stellar disk (collisionless part) */
extern double  M_GAS;    /* gas mass in disk */
extern double  M_BULGE;  /* mass of bulge */

extern double  halo_spinfactor;  /* computed streamin of dark matter */






extern double G;            /* gravitational constant */
extern double H0;           /* Hubble constant */
extern double UnitTime_in_s;
extern double UnitMass_in_g;
extern double UnitLength_in_cm;
extern double UnitVelocity_in_cm_per_s;
extern double UnitTime_in_Megayears;



/* particle data */

extern double    *vmax2_halo,*vmax2_disk,*vmax2_bulge,*vmax2_gas;

extern double    *xp_halo,*yp_halo,*zp_halo,*mp_halo;
extern double    *xp_disk,*yp_disk,*zp_disk,*mp_disk;
extern double    *xp_bulge,*yp_bulge,*zp_bulge,*mp_bulge;
extern double    *xp_gas,*yp_gas,*zp_gas,*mp_gas,*u_gas;

extern double    *vxp_halo,*vyp_halo,*vzp_halo;
extern double    *vxp_disk,*vyp_disk,*vzp_disk;
extern double    *vxp_bulge,*vyp_bulge,*vzp_bulge;
extern double    *vxp_gas,*vyp_gas,*vzp_gas;











double  LL;       /* LL = extension of fields in R and z. */






extern double **Dphi_z,**Dphi_R,**Dphi_z_dR;  /* derivatives of total potential */

extern double *epi_gamma2,*epi_kappa2;  /* epicycle gamma^2  */ 



/* halo velocity fields */

extern double **VelDispRz_halo;
extern double **VelDispPhi_halo;
extern double **VelVc2_halo;
extern double **VelStreamPhi_halo;
extern double **VelDispRz_dR_halo;


/* bulge velocity fields */

extern double **VelDispRz_bulge;
extern double **VelDispPhi_bulge;
extern double **VelVc2_bulge;
extern double **VelStreamPhi_bulge;
extern double **VelDispRz_dR_bulge;


/* disk velocity fields */

extern double **VelDispRz_disk;
extern double **VelDispPhi_disk;
extern double **VelVc2_disk;
extern double **VelStreamPhi_disk;
extern double **VelDispRz_dR_disk;




/* auxiliary field */

extern double *xl,*yl,*D2yl;
extern double *list_z,*list_R,*list_RplusdR;
