/*
 * helium_star.h: derived class for evolution of naked helium core
 *                burning stars.
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class helium_star
 *
 *....................................................................
 */

#ifndef    _HELIUM_STAR 
#   define _HELIUM_STAR

#include "single_star.h"
#include "white_dwarf.h"
#include "helium_giant.h"

class         main_sequence;
class         hyper_giant;
class         hertzsprung_gap;
class         sub_giant;
class         horizontal_branch;

/*-----------------------------------------------------------------------------
 *  helium_star  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class helium_star : public single_star {
      private :

        //real final_core_mass;
        //real CO_core_mass();// the core of a He MS is not defined yet
      
      public :

         //helium_star(main_sequence &);
         //helium_star(hyper_giant &);
         helium_star(hertzsprung_gap &);
         helium_star(sub_giant &);
         helium_star(horizontal_branch &);
         helium_star(node* n) : single_star(n) {}

         ~helium_star() {}

         //bool star_with_COcore() {return true;}   

	 stellar_type get_element_type();

	 void adjust_next_update_age(); 

	 void instantaneous_element();
	 void evolve_element(const real);
    real nucleair_evolution_timescale();

	 void update();
	 //void stellar_wind(const real);
	 void update_wind_constant();
	 void create_remnant();
	 bool hydrogen_envelope_star() {return false;}
	 real temperature();
	 
//		Mass transfer utilities.
	 real accretion_limit(const real, const real);
     star* subtrac_mass_from_donor(const real, real&);
	 star* reduce_mass(const real);
     real add_mass_to_accretor(const real, bool);
	 real add_mass_to_accretor(real, const real, bool);
     void adjust_accretor_age(const real, const bool=true);
     void adjust_age_after_wind_mass_loss(const real mdot,
                                          const bool rejuvenate);
    
    
//              Stability rourines.
        real zeta_adiabatic();
        real zeta_thermal();
        real gyration_radius_sq();

        // not private because of super_giant::initial_CO_core_mass()
        real final_CO_core_mass(const real initial_mass);
	
	friend helium_giant::helium_giant(helium_star &);
	friend white_dwarf::white_dwarf(helium_star &, stellar_type);
   };
#endif 		// _HELIUM_STAR
