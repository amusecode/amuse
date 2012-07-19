/*
 * white_dwarf.h: derived class for white dwarfs
 *
 *.....................................................................
 *     version 1:  Jan 1994   Simon F. Portegies Zwart
 *     version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class white_dwarf
 *
 *....................................................................
 */

#ifndef    _WHITE_DWARF
#   define _WHITE_DWARF

#include "single_star.h"
#include "neutron_star.h"

class hertzsprung_gap;
class super_giant;
class sub_giant;
class helium_star;
class helium_giant;

/*-----------------------------------------------------------------------------
 *  white_dwarf  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class white_dwarf : public single_star {
      public :

         white_dwarf(hertzsprung_gap &);
         white_dwarf(sub_giant &);
         white_dwarf(super_giant &);
         white_dwarf(helium_star &);
         white_dwarf(helium_giant &);
         white_dwarf(node* n) : single_star(n) {}
	 
         ~white_dwarf() {}

	 stellar_type get_element_type();
	 bool remnant() {return true;}
	 bool hydrogen_envelope_star() {return false;}

	 void instantaneous_element();
	 void evolve_element(const real);

         void update();
         void thermo_nucleair_flash(const real);
         void nova(const real);
         void common_envelope(const real);

//              Mass transfer utilities.
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
        real accretion_limit(const real, const real); 
        real minimum_steady_burning(const real); 
        real maximum_steady_burning(const real); 
        real add_mass_to_accretor(const real);
        real add_mass_to_accretor(real, const real);
        void adjust_accretor_age(const real, const bool=true);
        void adjust_next_update_age() {/* do nothing */}
        void accrete_from_envelope(const real);

 
        star* merge_elements(star*);

//		Mass transfer stability
        real zeta_thermal();
        real gyration_radius_sq();

        friend neutron_star::neutron_star(white_dwarf&);
      };

#endif 		// _WHITE_DWARF

