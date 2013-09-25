/*
 * thorne_zytkow.h: derived class for Thorne-Zytkow
 *
 *.....................................................................
 *     version 1:  Jan 1994   Simon F. Portegies Zwart
 *     version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class thorne_zytkow 
 *
 *....................................................................
 */

#ifndef    _THORNE_ZYTKOW 
#   define _THORNE_ZYTKOW

#include "single_star.h"
#include "super_giant.h"
#include "neutron_star.h"
#include "black_hole.h"

class main_sequence;

/*-----------------------------------------------------------------------------
 *  thorne_zytkow  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class thorne_zytkow : public single_star {
      public :

         thorne_zytkow(main_sequence &);
         thorne_zytkow(node* n) : single_star(n) {}
         ~thorne_zytkow() {}

	 stellar_type get_element_type() {return Thorn_Zytkow;}
	 bool giant_star()             {return true;}

	 void instantaneous_element();
	 void evolve_element(const real);
         void accrete_from_envelope(const real);
         void update();

//           Mass transfer utilities.
	 real add_mass_to_accretor(real, bool, const real = -1.);
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
        real mdot_limit(const real);
  //        void adjust_accretor_age(const real, const bool=true);
	void stellar_wind(const real);
        void adjust_next_update_age();

	real gyration_radius_sq();

	friend neutron_star::neutron_star(thorne_zytkow &);
	friend black_hole::black_hole(thorne_zytkow &);
        friend void super_giant::adjust_accretor_age(const real, 
						const bool rejuvenate);
};
#endif 		// _THORNE_ZYTKOW

//         void adjust_initial_star();
//	 real stellar_radius(const real m, const real t);
