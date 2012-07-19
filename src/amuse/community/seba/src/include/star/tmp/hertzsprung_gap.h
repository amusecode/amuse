/*
 * hertzsprung_gap.h: derived class for Hydrogen shell burning stars.
 *
 *.....................................................................
 *     version 1:  Jan 1994   Simon F. Portegies Zwart
 *     version 2:  
 *...................................................................
 *     This file includes:
 *  1) definition of class hertzsprung_gap
 *
 *....................................................................
 */

#ifndef    _HERTZSPRUNG_GAP
#   define _HERTZSPRUNG_GAP

#include "single_star.h"
#include "sub_giant.h"
#include "horizontal_branch.h"
#include "helium_star.h"

class main_sequence;

/*-----------------------------------------------------------------------------
 *  hertzsprung_gap  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class hertzsprung_gap : public single_star 
    {
      public :
         hertzsprung_gap(main_sequence &);
         hertzsprung_gap(node* n) : single_star(n) {}
         ~hertzsprung_gap() {}
	 
        stellar_type get_element_type() {return Hertzsprung_Gap;}
        bool giant_star()             {return true;}
	
        bool magnetic() {
	  return (low_mass_star() && get_total_mass()>=
		  cnsts.parameters(magnetic_mass_limit))?true:false;
	}

        void instantaneous_element();
        void evolve_element(const real);
        real gyration_radius_sq();

//           Mass transfer utilities.
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
         void adjust_accretor_age(const real, const bool);
        void adjust_next_update_age();
	void update_wind_constant();

//           Mass transfer Stability.
        real zeta_adiabatic();
        real zeta_thermal();

//           Core mass determination
	real TAMS_helium_core_mass();
	void evolve_core_mass(const real dt);

//	     Spectral type feature detection.
        void detect_spectral_features();

//              Spiral in and common envelope.

//           Friend constructors.
	 friend sub_giant::sub_giant(hertzsprung_gap &);
	 friend horizontal_branch::horizontal_branch(hertzsprung_gap &);
	 friend helium_star::helium_star(hertzsprung_gap &);
	 friend white_dwarf::white_dwarf(hertzsprung_gap &);
    };

#endif 		// _HERTZSPRUNG_GAP

//real stellar_radius(const real, const real);
//        void adjust_initial_star();
