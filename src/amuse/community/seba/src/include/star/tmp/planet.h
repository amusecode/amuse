/*
 * planet.h: derived class for evolution of planet
 *     is non helium burning extreme low mass star.
 *
 *.....................................................................
 *     version 1:  Jan 1994   Simon F. Portegies Zwart
 *     version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class planet
 *
 *....................................................................
 */

#ifndef    _PLANET
#   define _PLANET

//#include <stdlib.h>
#include "stdinc.h"
#include "starlab_constants.h"
#include "star.h"
//#include "ioserror.h"

		// Known class declarations.
class base_element;

/*-----------------------------------------------------------------------------
 *  planet  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class planet : public star {
      protected :
//         static ioserror error;
      public :
		// Constructors.
         planet() : star() {}
         planet(super_giant &);
         planet(main_sequence &);
         planet(base_element * b) : star(b) {}
		// Destructors.
         ~planet() {}
//              Member acces functions.
        stellar_type get_element_type() {return Planet;}
        bool remnant() {return FALSE;}

		// Member function definition.
        void adjust_initial_star() {
             relative_age = max(current_time, 0);}
        void  initialize_element();
        void instantaneous_element(const real);
        void evolve_element(const real);
        void update();
        real helium_core_mass();
        void thermo_nucleair_flash();
        void dwarf_nova();
        void common_envelope();
//              Mass transfer utilities.
        real mass_transfer_timescale();
        void reduce_donor_mass(const real);
        real subtrac_mass_from_donor(const real);
        real accretion_limit(const real, const real); 
        real add_mass_to_accretor(const real);
        real add_mass_to_accretor(real, const real);
        void adjust_accretor_age(const real);
        void adjust_next_update_age() {/* do nothing */}
        void accrete_from_envelope(const real);
//		Merge routines.
        void merge_elements(base_element*);
//		Mass transfer stability
        real zeta_thermal();
        void err_dump();


		// Debugging utilities.
        	// Friend functions.
//	In a future implementation AIC could occur.
//	Then the white dwarf makes friends.
        friend neutron_star::neutron_star(white_dwarf&);
      };
//base_element *wd = (base_element*)new white_dwarf(exemplar());
#endif 		// _PLANET
