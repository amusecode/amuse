/*
 * horizontal_branch.h: derived class for core-helium burning stars.
 *
 *.....................................................................
 *    version 1:    Jan 1994   Simon F. Portegies Zwart
 *    version 1.1:  Jan 1998   Simon F. Portegies Zwart
 *...................................................................
 *     This file includes:
 *  1) definition of class horizontal_branch
 *
 *....................................................................
 */

#ifndef    _HORIZONTAL_BRANCH 
#   define _HORIZONTAL_BRANCH

#include "single_star.h"
#include "super_giant.h"
#include "helium_star.h"

class hertzsprung_gap;
class sub_giant;
class single_star;

/*-----------------------------------------------------------------------------
 *  horizontal_branch  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class horizontal_branch : public single_star {
      public :

         horizontal_branch(sub_giant&);
         horizontal_branch(hertzsprung_gap&);
         horizontal_branch(node* n) : single_star(n) {}
         ~horizontal_branch() {}
	 
         stellar_type get_element_type() {return Horizontal_Branch;}
	 bool giant_star()             {return true;}
         bool star_with_COcore() {return true;}   

	 void instantaneous_element();
         void evolve_element(const real);
  void evolve_core_mass(const real dt) {cerr<<"void evolve core mass"<<endl;} 

//            Mass transfer utilities.
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
         void adjust_accretor_age(const real, const bool);
        void adjust_next_update_age();
	void update_wind_constant();

//	     Mass transfer stability
        real zeta_adiabatic();
        real zeta_thermal();
        real gyration_radius_sq();

//              Friend constructors
	 friend super_giant::super_giant(horizontal_branch &);
	 friend helium_star::helium_star(horizontal_branch &);


      void evolve_core_mass(const real time,
			    const real mass,
			    const real z);
      void evolve_core_mass();
      void update();
    
    //Small envelope behaviour
    real small_envelope_core_luminosity(const real time, 
                             const real mass, const real m_core, const real z);
    real small_envelope_core_luminosity();
    real small_envelope_core_radius(const real time, 
                            const real mass, const real m_core, const real z);
    real small_envelope_core_radius();
    
    real helium_core_radius(const real time, 
                                    const real mass, const real m_core, const real z);
    real helium_core_radius();
    
    };
#endif 		// _HORIZONTAL_BRANCH

//        void adjust_initial_star();
//         real stellar_radius(const real, const real);

