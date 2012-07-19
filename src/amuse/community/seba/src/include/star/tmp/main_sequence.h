/*
 * main_sequence.h: derived class for evolution of stars in the 
 *     core-hydrogen burning stage.
 *
 *.....................................................................
 *     version 1:  Jan 1994   Simon F. Portegies Zwart
 *     version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class main_sequence
 *
 *....................................................................
 */

#ifndef    _MAIN_SEQUENCE 
#   define _MAIN_SEQUENCE

#include "single_star.h"
//#include "proto_star.h"
#include "brown_dwarf.h"
#include "hyper_giant.h"
#include "hertzsprung_gap.h"
#include "thorne_zytkow.h"
#include "helium_star.h"

		// Known class declarations.
class proto_star;
/*-----------------------------------------------------------------------------
 *  main_sequence  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class main_sequence : public single_star {
      private:
  
          real main_sequence_core_mass();
	  real main_sequence_core_radius();
	  void adjust_donor_age(const real mdot);
  
      public :

         main_sequence(node* n) : single_star(n) {}
         main_sequence(proto_star& p);
         ~main_sequence() {}

        stellar_type get_element_type() {return Main_Sequence;}
        bool remnant() {return false;}
        bool magnetic() {
	  return (low_mass_star() && get_total_mass()>=
		  cnsts.parameters(magnetic_mass_limit))?true:false;
	}

	void adjust_next_update_age();
	real nucleair_evolution_timescale();
	void instantaneous_element();
        void evolve_element(const real);
        real bolometric_correction();
        void detect_spectral_features();
	void update_wind_constant();
	void update();
        real final_core_mass();

	
//		Mass transfer utilities.
        star* subtrac_mass_from_donor(const real, real&);
        star* reduce_mass(const real);
        void adjust_accretor_age(const real, const bool);
	

//		Mass transfer stability
        real zeta_adiabatic();
        real zeta_thermal();
	real gyration_radius_sq();

//		Spiral in and common envelope.
	star* merge_elements(star*);

	
//              Friend functions.
        friend brown_dwarf::brown_dwarf(main_sequence &);
	friend hertzsprung_gap::hertzsprung_gap(main_sequence &);
	friend hyper_giant::hyper_giant(main_sequence &);
	friend thorne_zytkow::thorne_zytkow(main_sequence &);
	friend helium_star::helium_star(main_sequence &);
      };

#endif 		// _MAIN_SEQUENCE


        //real stellar_radius(const real, const real);
	//real add_mass_to_accretor(const real);
	//real add_mass_to_accretor(real, const real);
