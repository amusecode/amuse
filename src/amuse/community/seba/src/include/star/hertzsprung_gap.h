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
        void adjust_accretor_age(const real mdot, const bool rejuvenate);
        void adjust_age_after_mass_loss(const real mdot, const bool rejuvenate);
        void adjust_next_update_age();
        void update_wind_constant();
        real  add_mass_to_accretor(const real, bool);
        real  add_mass_to_accretor(real, const real, bool);

//           Mass transfer Stability.
        real zeta_adiabatic();
        real zeta_thermal();

//	     Spectral type feature detection.
        void detect_spectral_features();

//              Spiral in and common envelope.

//           Friend constructors.
	 friend sub_giant::sub_giant(hertzsprung_gap &);
	 friend horizontal_branch::horizontal_branch(hertzsprung_gap &);
	 friend helium_star::helium_star(hertzsprung_gap &);
	 friend white_dwarf::white_dwarf(hertzsprung_gap &, stellar_type);



      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //Metalicity dependent functions:
      //      real hertzsprung_gap_time(const real mass, const real z);
      //      real hertzsprung_gap_time();
      real terminal_hertzsprung_gap_luminosity(const real mass, 
					       const real z);
      real terminal_hertzsprung_gap_luminosity();
      real terminal_hertzsprung_gap_radius(const real mass, 
					   const real mass_tot, const real z);
      real terminal_hertzsprung_gap_radius();
      real hertzsprung_gap_luminosity(const real time,
				      const real mass, 
				      const real z);
      real hertzsprung_gap_luminosity(const real time);
      real hertzsprung_gap_luminosity();
      real hertzsprung_gap_radius(const real time,
				  const real mass,
                  const real mass_tot, 
				  const real z);
      real hertzsprung_gap_radius(const real time);
      real hertzsprung_gap_radius();
        
      real hertzsprung_gap_time(const real mass, const real z);
      real hertzsprung_gap_time();


      void evolve_core_mass(const real time,
			    const real mass,
			    const real z, const real m_core_old);
      void evolve_core_mass();
      real hertzsprung_gap_core_mass(const real time, 
				     const real mass,
				     const real z, const real m_core_old);

    //Small envelope behaviour
        real small_envelope_core_luminosity(const real mass, const real m_core, const real z);
        real small_envelope_core_luminosity();
        real small_envelope_core_radius(const real mass, const real m_core, const real z);
        real small_envelope_core_radius();
        real helium_core_radius(const real mass, const real m_core, const real z);
        real helium_core_radius();

    
    };

#endif 		// _HERTZSPRUNG_GAP

//        void adjust_initial_star();
