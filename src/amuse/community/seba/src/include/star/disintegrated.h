/*
 * disintegrated.h: derived class for evolution of nothing
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class disintegrated
 *
 *....................................................................
 */

#ifndef    _DISINTEGRATED 
#   define _DISINTEGRATED

#include "stdinc.h"
#include "single_star.h"
#include "black_hole.h"

class          super_giant;
class          helium_giant;
class          white_dwarf;
/*-----------------------------------------------------------------------------
 *  disintegrated  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class disintegrated : public single_star {
      protected :
         real suddenly_lost_mass;
      public :
         disintegrated(super_giant &);
         disintegrated(helium_giant &);
         disintegrated(white_dwarf &);
         disintegrated(node* n) : single_star(n) {suddenly_lost_mass=0;}
         ~disintegrated() {}

         stellar_type get_element_type() {return Disintegrated;}
	 bool remnant() {return false;}
	 bool hydrogen_envelope_star() {return false;}
         int no_of_elements() {return 0;}

         void instantaneous_element();
	 void evolve_element(const real);
         void update();
         bool super_nova();
         void adjust_next_update_age() {}

        star* reduce_mass(const real);
        real mass_transfer_timescale(mass_transfer_type &type);
        star* subtrac_mass_from_donor(const real, real&);
        real add_mass_to_accretor(real, bool, const real = -1.);

        star* merge_elements(star*);

        real temperature();
        real magnitude();
        real bolometric_correction();
        void stellar_wind(const real);
        real gyration_radius_sq();

        real sudden_mass_loss();

   };
#endif 		// _DISINTEGRATED        
