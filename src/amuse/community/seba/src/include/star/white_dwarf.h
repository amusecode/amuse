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
    protected:
    
        stellar_type white_dwarf_type;
    
      public :

         white_dwarf(hertzsprung_gap &, stellar_type);
         white_dwarf(sub_giant &, stellar_type);
         white_dwarf(super_giant &, stellar_type);
         white_dwarf(helium_star &, stellar_type);
         white_dwarf(helium_giant &, stellar_type);
         white_dwarf(node* n) : single_star(n) {}
	 
         ~white_dwarf() {}

    real get_evolve_timestep();
	 //stellar_type get_element_type();
    stellar_type get_element_type() {return white_dwarf_type;}

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
        real add_mass_to_accretor(const real, bool);
        real add_mass_to_accretor(real, const real, bool);
        void adjust_accretor_age(const real, const bool=true);
        void adjust_next_update_age() {/* do nothing */}
        void accrete_from_envelope(const real);

 
        star* merge_elements(star*);

//		Mass transfer stability
        real zeta_thermal();
        real zeta_adiabatic();
        real gyration_radius_sq();

        friend neutron_star::neutron_star(white_dwarf&);
      };

#endif 		// _WHITE_DWARF

