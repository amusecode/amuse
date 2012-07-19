/*
 * sub_giant.h: derived class for evolution of 
 *              convective hydrogen shell burning stars.
 *
 *.....................................................................
 *    version 1:    Jan 1994   Simon F. Portegies Zwart
 *    version 1.1:  Jan 1998   Simon F. Portegies Zwart
 *...................................................................
 *     This file includes:
 *  1) definition of class sub_giant
 *
 *....................................................................
 */

#ifndef    _SUB_GIANT
#   define _SUB_GIANT

#include "single_star.h"
#include "helium_star.h"
#include "horizontal_branch.h"

		// Known class declarations.
class           hertzsprung_gap;

/*-----------------------------------------------------------------------------
 * sub_giant  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class sub_giant : public single_star 
    {
      public :

	sub_giant(hertzsprung_gap &);
        sub_giant(node* n) : single_star(n) {}
        ~sub_giant() {}

        stellar_type get_element_type() {return Sub_Giant;}
        bool giant_star()             {return true;}
        bool magnetic() {
	  return (low_mass_star() && get_total_mass()>=
		  cnsts.parameters(magnetic_mass_limit))?true:false;
	}

        void instantaneous_element();
        void evolve_element(const real);
        real get_evolve_timestep();
        
//            Mass transfer utilities.
        real  add_mass_to_accretor(const real, bool );
        real  add_mass_to_accretor(real, const real, bool);
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);

         void adjust_accretor_age(const real, const bool);
        void adjust_next_update_age();
        void detect_spectral_features();
	void update_wind_constant();

//           Mass transfer Stability.
        real zeta_adiabatic();
        real zeta_thermal();

	real gyration_radius_sq();
         
//              Friend constructors
	 friend horizontal_branch::horizontal_branch(sub_giant &);
	 friend helium_star::helium_star(sub_giant &);
	 friend white_dwarf::white_dwarf(sub_giant &, stellar_type);

      void evolve_core_mass(const real time,
			    const real mass,
			    const real z);
      void evolve_core_mass();
      real sub_giant_core_mass(const real time,
                            const real mass,
                            const real z);      
        //Small envelope behaviour
        real small_envelope_core_luminosity(const real mass, const real m_core, const real z);
        real small_envelope_core_luminosity();
        real small_envelope_core_radius(real mass, const real m_core, const real z);
        real small_envelope_core_radius();
        real helium_core_radius(const real mass, const real m_core, const real z);
        real helium_core_radius();
        
    };
#endif 		// _SUB_GIANT
//        void adjust_initial_star();
//        real stellar_radius(const real, const real);

