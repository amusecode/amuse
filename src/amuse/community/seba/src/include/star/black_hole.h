/*
 * black_hole.h: derived class for "evolution" of black hole
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class black_hole
 *
 *....................................................................
 */

#ifndef    _BLACK_HOLE 
#   define _BLACK_HOLE

#include "single_star.h"

class        hyper_giant;
class        super_giant;
class        thorne_zytkow;
class        helium_giant;
class        neutron_star;
/*-----------------------------------------------------------------------------
 *  black_hole  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class black_hole : public single_star {
      private:

	real suddenly_lost_mass;

	real black_hole_mass();

      public :
         black_hole(hyper_giant &);
         black_hole(super_giant &);
         black_hole(thorne_zytkow &);
         black_hole(helium_giant &);
         black_hole(neutron_star &);
         black_hole(node* n) : single_star(n) {suddenly_lost_mass=0;}

         ~black_hole() {}
    
        real get_evolve_timestep();
        stellar_type get_element_type() {return Black_Hole;}
        bool remnant() {return true;}
        bool hydrogen_envelope_star() {return false;}
	
		// Member function definition.
         void evolve_element(const real);
         void update();

        void instantaneous_element();
        star* reduce_mass(const real);
        real accretion_limit(const real, const real);
        real mdot_limit(const real);
        star* subtrac_mass_from_donor(const real, real&);
        real add_mass_to_accretor(real, bool, const real = -1.);
        void accrete_from_envelope(const real);
        star* merge_elements(star*);

        bool super_nova();
        void direct_hit();
        real aic_binding_energy();
	
	real gyration_radius_sq();
	real angular_momentum();
	
        real sudden_mass_loss();

	real get_radius();
	real get_effective_radius() {return get_radius();}

	// Friend functions.
	     // poor guy.
   };
#endif 		// _BLACK_HOLE
