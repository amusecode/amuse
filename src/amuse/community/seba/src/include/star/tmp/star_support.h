
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  star_support.h: derived class for element evolution systems.
 *          functions as derived class for the real elements.
 *.............................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of stellar types and spectral classes
 *
 *.............................................................................
 */
#ifndef  _STAR_SUPPORT
#  define  _STAR_SUPPORT

#include "stdinc.h"

enum stellar_type {Star_Cluster = -4, 
		   Static_Star = -3, SPZDCH_Star = -2, NAS = -1, 
		   Proto_Star = 0, Planet, Brown_Dwarf,
		   Main_Sequence, Hyper_Giant, Hertzsprung_Gap,
		   Sub_Giant, Horizontal_Branch, Super_Giant=8,
		   Carbon_Star, Helium_Star, Helium_Giant, 
                   Carbon_Dwarf, Helium_Dwarf, Oxygen_Dwarf,
		   Thorn_Zytkow=15,
		   Xray_Pulsar, Radio_Pulsar, Neutron_Star, Black_Hole,
		   Disintegrated, Double, no_of_stellar_type};

enum stellar_type_summary {ZAMS=0, Early_Giant, Late_Giant,
			   Helium_Remnant, White_Dwarf, 
			   Neutron_Remnant, Inert_Remnant,
			   Unspecified, Undefined, no_of_star_type_summ};

enum spectral_class {O5=0, O6, O7, O8, O9, O95, 
                     B0=6, B05, B1, B2, B3, B5, B6, B7, B8, B9, B95,
                     A0=17, A1, A2, A3, A4, A5, A7,
                     F0=24, F2, F3, F5, F6, F7, F8,
                     G0=31, G1, G2, G5, K0, K5, M0, M5, M8,
                     he=40, wd, ns, bh, bd, di,
                     binary_star, no_spectral_class};

enum star_type_spec {NAC=0, Emission, Blue_Straggler, Barium, 
                     Rl_filling, Runaway, 
                     Merger, Accreting, Dsntgr, no_of_spec_type};

enum luminosity_class {I=0, II, III, IV, V,
		       no_luminosity_class};
				 
enum mass_transfer_type {Unknown=0, Nuclear, AML_driven, Thermal, Dynamic};

enum supernova_type {NAT=0, SN_Ia, SN_Ib, SN_Ic, 
		     SN_II, SN_IIL, SN_IIP, SN_IV,
		     no_of_supernova_type};

//              polytropic indices.
//enum polytrope {convective=1, composite, radiative};
// index          3             2          1.5
				 
const char * type_string(stellar_type);
const char * type_string(stellar_type_summary);
const char * type_string(spectral_class);
const char * type_string(luminosity_class);
const char * type_string(star_type_spec);
const char * type_string(mass_transfer_type);
const char * type_string(supernova_type sn_type);
const char * type_short_string(star_type_spec);
const char * type_short_string(stellar_type);
const char * type_short_string(stellar_type_summary);
const char * type_short_string(mass_transfer_type);
const char * type_short_string(spectral_class);

stellar_type extract_stellar_type_string(const char*);
star_type_spec extract_stellar_spec_summary_string(const char*);

stellar_type_summary extract_stellar_type_summary_string(const char*);
stellar_type_summary summarize_stellar_type(stellar_type);

// Conversion routines for SeBa to BSE
int convert_SeBa_stellar_type_to_BSE(stellar_type tpe);
stellar_type convert_BSE_to_SeBa_stellar_type(int tpe);

bool remmant(stellar_type);
bool post_supernova_star(stellar_type);
supernova_type type_of_supernova(stellar_type progenitor);

void combine_ubvri(real Up, real Bp, real Vp, real Rp, real Ip,
                   real Us, real Bs, real Vs, real Rs, real Is,
		   real &U, real &B, real &V, real &R, real &I);

//Because of the dyn*, the following functions are in starbase.h 
//void get_ubvri_star(dyn *bi, stellar_type& stype,
//		    real& U, real& B, real& V, real& R, real& I);
//  
//void get_Lubvri_star(dyn *bi, stellar_type& stype,
//		     real& Lu, real& Lb, real& Lv, real& Lr, real& Li);

void ltm_to_ubvri(const real logl,
		  const real logt,
		  const real mass,
		  real& U,
		  real& B,
		  real& V,
		  real& R,
		  real& I);

void ltm_to_ubvri(const real,
		  const real,
		  const real,
		  const real,
		  const real,
		  const real,
		  real& U,
		  real& B,
		  real& V,
		  real& R,
		  real& I);

spectral_class get_spectral_class(const real temperature);
luminosity_class get_luminosity_class(const real temperature, const real lum);
real lum_class_limit(const real, luminosity_class);

#endif          // _STAR_SUPPORT




