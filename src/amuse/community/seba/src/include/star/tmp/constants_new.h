
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  constants.h: header file for constants.
 *           
 *.............................................................................
 *    version 1:  May 1993   Simon F. Portegies Zwart
 *    version 2:  Jan 1998   Simon F. Portegies Zwart
 *.............................................................................
 *     This file includes:
 *  1) definition of physical constants, assumed numbers and pre set
 *     numbers.
 *
 *.............................................................................
 */
#ifndef     _CONSTANTS
#   define  _CONSTANTS

#include "star_support.h"
#include "double_support.h"
#include "star_support.h"

#define TRUE 1
#define FALSE 0

#define SEBA_VERSION 2.0

enum kick_distribution = {Delta_function,
			  Gaussian,
			  Maxwellian, 
			  Paczinski};

// Mathematical constants.
class {
  public:
  struct {
d    const real one_third = 0.333333333333333;
d    const real two_third = 0.666666666666666;
  } math;
  struct  {
    const real G = 6.67e-8;                      // [cgs]
    const real C = 2.9979e+10;                   // [cm/s]

d    const real Myear = 3.15e+13;                 // [s]
    const real days   = 86400;                   // [s]
    const real km_s   = 1.0e+5;                  // [cm/s]
  } phys;
  struct {
    const real mass   = 1.989e+33;               // [gr]
    const real radius = 6.96e+10;                // [cm]
    const real lumi   = 3.862e+33;               // [erg/s]
  }  sun;
  struct  {
    const real parsec  = 3.0857e+18; 	         // [cm]
    const real AU      = 1.496e+13; 	         // [cm]
    const real Rsun2AU = 4.65e-3; 		 // [Rsun]
    const real Rsun2pc = 2.26e-8; 		 // [Rsun]
  } astro;
  struct  {
    const real nucleair_efficiency = 0.0007;  // nucleair energy production eff. 
    const real black_hole   = 0.1;      // 0.1 Mass fraction accretes to core.
    const real neutron_star = 0.05;     // 0.01 Mass fraction accretes to core.
    const real white_dwarf  = 0.01;     // 0.01 Mass fraction accretes to core.
  } accretion;
  struct  {
    const real chandrasekhar_mass  =  1.44;      // [Msun].
    const real M_NS                =  1.35;      // [Msun].
    const real R_NS                =  1.5e-5;    // [Rsun].
    const real M_NSmax             =  2.0;	 // [Msun].
    const real M_NSmin             =  8.0;       // [Msun].
    const real M_HE_NS 	           =  2.2;       // [Msun].
    const real M_BHmin             =  40.0;      // [Msun].
    const real M_HE_BH	           =  15.0;      // [Msun].
    const real M_MSmax	           = 100.0;      // [Msun].
    const real M_MSmin	           =   0.075;    // [Msun]. (Pop III: 0.095)
    const real M_HEmin      	   =   2.20;     // [Msun]. (3: v/d Heuvel' 90)
    const real M_HE_WDmin          =   0.33;     // [Msun].
c    const real LOW_MASS_STAR       =   1.5;      // [Msun].
    const real M_MSW               =   0.6;      // Magnetic Stellar Wind AML
c    const real MEDIUM_MASS_STAR    =  15.0;      // [Msun].
  } mass;
  struct  {
    const bool exp_HE_lifetime = FALSE;	// exponential he star lifetime fits.
    const bool Hyper_critical_accretion = FALSE;  
    bool GBROWN  = FALSE // Special neutron star in common envelope treatment.
                  // According to Gerry Brown, Neutron star collapse to 
                  // black holes during the spiral-in.
    const bool AARSETH = FALSE; //This parameter is defined in star.C
    const real mean_kick  = 450.0;			// km/s
    const real sigma_kick = 90.0;			// km/s
    const distribution distr_kick = Paczinski;
  } evolution;
  struct  {
    const real  COROTATION_E_LIMIT = 0.001;  // minimum eccentricity for corotation.
    const real  TNF_MASS_FR        = 0.05;   // 0.05Envelope fraction causes TNF.(CV)
    const real  Be_limit = 0.10;      // Accrete fraction Msun for Be star.
    const real  Ba_limit = 0.01;      // Accrete fraction of mass for Ba.
    const real  Bs_limit = 0.05;      // fraction of TO mass for Bs detection.
    const real  ALPHA     = 1.0;     // 0<alpha<1: masstransfer efficiency.
    const real  LAMBDA    = 0.5;     // envelope binding constant (0.5).
    const real  GAMMA_mb  = 2.5;     // (or 4) magnetic braking index.
    const real  OVERSHOOT = 0.125;   // 0.125 convective overshooting parameter.
    const real  TEN_PERCENT = 0.1;
    const real  STABLE_ZETA = 15;
    const real  J_SPI_FR    = 0.5;	//Pols 93 advises 1/3
    const real  TIDAL_RANGE = 4.0;	// tidal synchronization distance in [Rstar].
    const real  SPI_TIME	     = 1.e-7;   // Spiral in duration time.
    const real  POST_AGB_TIME      = 0.005;   // time to lose giant envelope.
  } model;
c  struct  {
c    const real MIN_DT_FACTOR       = 0.01;
c    const real timestep_factor     = 0.01;
c                           // WARNING Tricky number: enlargment causes errors!!
c                           // (0.01) as does making the number too small!!!!
c    const real  ABSOLUTE_DT_MIN    = 1.e-7;  
c    const real  ABSOLUTE_MDOT_MIN  = 1.e-10;
c    const real  Maximal_timestep   = 2.5;     // Safety parameter mostly 2.5.
c    const int  MAXIMUM_RECURSIVE_CALLS = 1000;
c  } safety;
} static cnsts;

#define N_FATE 8

enum binary_history  {ms_ms = 0, bs_ms, bs_bs, he_ms, heN_ms, 
                      wr_ms, he_he, rscvn, 
		      wuma, wdxb, 
                      lmxb, mmxb, hmxb, spi, spi_2, 
                      wd_ms, wd_he, wd_wd, wd_ns, wdXns, 
                      pr_ms, pr_he, pr_st, prXst, pr_pr, prXpr,
                      ns_ms, ns_he, ns_st, nsXst, ns_ns, nsXns, st_st, 
                      no_of_type};


enum Be_binaries    {be_ms=0, be_he, be_wd, be_ns, no_of_be_binary};

#endif		// _CONSTANTS


