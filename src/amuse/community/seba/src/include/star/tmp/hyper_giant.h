/*
 * hyper_giant.h: derived class for evolution of Hyper_giants.
 *                was wolf_rayet.h
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class hyper_giant
 *
 *....................................................................
 */

#ifndef   _HYPER_GIANT
#   define _HYPER_GIANT

#include "single_star.h"
#include "helium_star.h"
#include "neutron_star.h"
#include "black_hole.h"

class main_sequence;

/*-----------------------------------------------------------------------------
 *  hyper_giant  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class hyper_giant : public single_star {
      private:

          real hyper_giant_core_mass();

      public :

         hyper_giant(main_sequence &);
         hyper_giant(node* n) : single_star(n) {}
         ~hyper_giant() {}

	 
         stellar_type get_element_type() {return Hyper_Giant;}
         bool giant_star()             {return true;}
         bool hydrogen_envelope_star() {return false;}

         void instantaneous_element();
         void evolve_element(const real);
         void create_remnant();
	 
//            Mass transfer utilities.
	 star* subtrac_mass_from_donor(const real, real&);
         real accretion_limit(const real, const real);

         void adjust_accretor_age(const real, const bool);
         void adjust_next_update_age();
         star* reduce_mass(const real);

//              Stability rourines.
         real zeta_thermal();
         real zeta_adiabatic();
	 real gyration_radius_sq();

	 friend helium_star::helium_star(hyper_giant &);
	 friend neutron_star::neutron_star(hyper_giant &);
	 friend black_hole::black_hole(hyper_giant &);
    };
#endif 		// _HYPER_GIANT


//         real zeta_adiabatic();
//         void adjust_initial_star();
//         real stellar_radius(const real m, const real t) 
//              {return 5;}
