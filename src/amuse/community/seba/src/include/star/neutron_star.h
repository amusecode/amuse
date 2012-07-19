/*
 * neutron_star.h: derived class for neutron stars.
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class neutron_star
 *
 *....................................................................
 */

#ifndef    _NEUTRON_STAR 
#   define _NEUTRON_STAR

#include "single_star.h"
#include "black_hole.h"

class hyper_giant;
class super_giant;
class thorne_zytkow;
class helium_giant;
class white_dwarf;

/*-----------------------------------------------------------------------------
 *  neutron_star  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class neutron_star : public single_star {
      private:

         real suddenly_lost_mass;

      public :

         neutron_star(hyper_giant &);
         neutron_star(super_giant &);
         neutron_star(thorne_zytkow &);
         neutron_star(helium_giant &);
         neutron_star(white_dwarf &);
         neutron_star(node* n) : single_star(n) {suddenly_lost_mass=0;}

         ~neutron_star() {}
	 
        stellar_type get_element_type();
        bool remnant() {return true;}
        bool hydrogen_envelope_star() {return false;}

	void instantaneous_element();
	void evolve_element(const real);

	void update();
	bool super_nova();
	void direct_hit();
	real aic_binding_energy();

        star* reduce_mass(const real);
        real accretion_limit(const real, const real);
        star* subtrac_mass_from_donor(const real, real&);
        real add_mass_to_accretor(const real);
        real add_mass_to_accretor(real, const real);
        void accrete_from_envelope(const real);

        star* merge_elements(star*);

        real neutron_star_mass(stellar_type);
	real neutron_star_radius();
	real accretion_luminosity(const real, const real);
	real magnetic_moment();
	real magnetic_field_decay(const real, const real);
      real magnetic_field_decay(const real, const real, const real);
      real magnetic_field_strength(const real, const real);
      real pulsar_spin_up(const real, const real);
      real pulsar_spin_down(const real);
      bool propeller(const real, const real);
      bool dead_pulsar();
      
      real fastness_parameter(const real);
      real dimension_less_accretion_torque(const real);
      real pulsar_propeller_torque(const real, const real);
      
      real moment_of_inertia();
      real period_derivative();
      real spindown_luminosity();
      real gyration_radius_sq();
      
      real sudden_mass_loss();
      
      friend black_hole::black_hole(neutron_star &);

   };
#endif 		// _NEUTRON_STAR
