/*
 *    star_to_dyn.h: transformation function from stellar evolution
 *		     package to dynamical package.
 *
 *.....................................................................
 *    version 1:  Aug 1996   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) Transformation routines from stellar evolution to n-body dynamcs.
 *
 *....................................................................
 */

#ifndef     _SSTAR_TO_DYN
#   define  _SSTAR_TO_DYN

#include  "starbase.h"
#include  "dyn.h"

#include "stdinc.h"
#include "starlab_constants.h"

/*-----------------------------------------------------------------------------
 *  sstar_to_dyn  --  
 *-----------------------------------------------------------------------------
 */

void  addstar(dyn*,
	      real t_rel=0,
	      stellar_type type=Main_Sequence,
	      bool verbose = false);
void sstar_stats(dyn*, bool, vec, bool);

vec conv_v_star_to_dyn(vec&, real, real);
vec anomalous_velocity(dyn*);
real get_effective_radius(dyn*);
real get_total_mass(dyn*);
stellar_type get_element_type(dyn*);
bool has_sstar(dyn*);

//void change_dynamical_mass_of_the_stars(dyn*);
bool stellar_evolution(dyn* b);
real sudden_mass_loss(dyn* b);

// Function for output.
star_state make_star_state(dyn*);


void extract_story_chapter(char* type, real&, real&,
                           real&, real&, real&,
                           real&, real&,
			   real&, real&, story& s);
void extract_line_text(char*, real&, real&, real&,
                       real&, real&, real&, real&,
		       real&, real&, story&);

bool merge_with_primary(star* primary, star *secondary);

void print_sstar_time_scales(dyn*);


real compute_projected_luminosity_radii(dyn * b, int axis,
					bool Llagr = true,
					int nzones = 10,
					bool nonlin=false,
					boolfn bf=NULL);

enum radius_fraction {rcore_rvirial, r10_r90_light, proj_r10_r90_light};
void print_fitted_king_model(const real fractional_radius,
			     const radius_fraction
			           which_fraction = proj_r10_r90_light);

void put_ubvri(dyn*);

// reduction software
void new_camera_position(vec&, vec, vec, int, int);
void rdc_and_wrt_movie(dyn*, bool, real, real, char);
void print_povray_header(dyn*, vec, int, real, real, int, int, real, 
			 int, int, bool);
void print_mpegplayer_param_file(ostream&,
				 int first_frame =1,
				 int last_frame = 100,
				 int horizontal = 120,
				 int vertical = 90,
				 int GOP_size = 10);


//The following functions are in ltm_to_ubvri
void get_ubvri_star(dyn *bi, stellar_type& stype,
		    real& U, real& B, real& V, real& R, real& I);
  
void get_Lubvri_star(dyn *bi, stellar_type& stype,
		     real& Lu, real& Lb, real& Lv, real& Lr, real& Li);

				 
#endif          // _SSTAR_TO_DYN
