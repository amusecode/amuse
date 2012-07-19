/*
 * brown_dwarf.h: derived class for evolution for non-hydrogen core
 *                burning stars.
 *
 *.....................................................................
 *     version 1:  Jan 1994   Simon F. Portegies Zwart
 *     version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class brown_dwarf
 *
 *....................................................................
 */

#ifndef    _BROWN_DWARF
#   define _BROWN_DWARF

#include "single_star.h"
//#include "proto_star.h"

class    main_sequence;
class    proto_star;

/*-----------------------------------------------------------------------------
 *  brown_dwarf  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class brown_dwarf : public single_star {
      private:

          real brown_dwarf_core_mass();

      public :

         brown_dwarf(node* n) : single_star(n) {}
         brown_dwarf(main_sequence &);
         brown_dwarf(proto_star &);
         ~brown_dwarf() {}

        stellar_type get_element_type();
        bool remnant() {return false;}

        void instantaneous_element();
        void evolve_element(const real);
        void update();

        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
        real accretion_limit(const real, const real); 
        real add_mass_to_accretor(const real);
        real add_mass_to_accretor(real, const real);
        void adjust_next_update_age() {/* do nothing */}
        void accrete_from_envelope(const real);

        star* merge_elements(star*);
        real zeta_thermal();
	real gyration_radius_sq();

      };
#endif 		// _BROWN_DWARF
