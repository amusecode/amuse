/*
 * super_giant.h: derived class for evolution of double shell burning stars.
 *
 *..........................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *..........................................................................
 *     This file includes:
 *  1) definition of class super_giant
 *
 *..........................................................................
 */

#ifndef    _SUPER_GIANT 
#   define _SUPER_GIANT

#include "single_star.h"
#include "helium_star.h"
#include "white_dwarf.h"
#include "neutron_star.h"
//#include "disintegrated.h"
#include "black_hole.h"

		// Known class declarations.
class horizontal_branch;

/*-----------------------------------------------------------------------------
 *  super_giant  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class super_giant : public single_star {
      protected:
	real second_dredge_up_time;
	real initial_CO_core_mass(const real initial_mass);

      public :

	super_giant(horizontal_branch & h);
        super_giant(node* n) : single_star(n) {
	  second_dredge_up_time = -1;}
        ~super_giant() {}
	
        stellar_type get_element_type() {return Super_Giant;}
        bool giant_star()             {return true;}
        bool star_with_COcore() {return true;}   

	void instantaneous_element();
	void evolve_element(const real);
	void create_remnant();

//           Mass transfer stability.
        real zeta_thermal();
	real gyration_radius_sq();
	
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
         void adjust_accretor_age(const real, const bool);
        void adjust_next_update_age();
	void update_wind_constant();

//              Friend constructors
	 friend helium_giant::helium_giant(super_giant &);
         friend white_dwarf::white_dwarf(super_giant &);
         friend neutron_star::neutron_star(super_giant &);
         friend black_hole::black_hole(super_giant &);
         friend disintegrated::disintegrated(super_giant &);

      };
#endif 		// _SUPER_GIANT


//	real stellar_radius(const real, const real);
//	void adjust_initial_star();
//	real helium_core_mass();
//	void evolve_without_wind(const real); 
//        real zeta_adiabatic();
