/*
 * helium_giant.h: derived class for evolution of helium giants
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class helium_giant
 *
 *....................................................................
 */

#ifndef    _HELIUM_GIANT
#   define _HELIUM_GIANT

#include "single_star.h"
#include "disintegrated.h"
#include "white_dwarf.h"
#include "neutron_star.h"
#include "black_hole.h"


class         helium_star;
class         super_giant;

/*-----------------------------------------------------------------------------
 *  helium_giant  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class helium_giant : public single_star {
     private:
  
         real CO_core_mass();
         real relative_helium_mass;
    
      public :
         helium_giant(super_giant &);
         helium_giant(helium_star &);
         ~helium_giant() {}

        stellar_type get_element_type();
        bool giant_star();
        bool remnant();
        bool hydrogen_envelope_star() {return false;}
        bool star_with_COcore() {return true;}   

        void evolve_element(const real);
        void instantaneous_element();
    
        real helium_giant_luminosity_core_mass_relation(
                                                    const real time, const real mass, const real z);
        real helium_giant_age_core_mass_relation(const real m_core, const real mass);
        real helium_giant_end_time(const real mass, const real mass_tot); 
        real small_envelope_mu(const real lum, const real mass_tot, const real m_core);
        void small_envelope_perturbation();
    
        real maximum_helium_giant_core_mass(const real mass);
        void evolve_core_mass(const real time,
                      const real mass, const real mass_tot);
        void evolve_core_mass();
        real helium_giant_core_mass(const real time,
                         const real mass);
         
        void stellar_wind(const real);
        real nucleair_evolution_timescale(); 
    
        void update();
        void update_wind_constant();
        real temperature();

        void create_remnant(const real mass, const real mass_tot, const real mc_core);
        real  add_mass_to_accretor(real, bool, const real = -1. );
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
        real accretion_limit(const real, const real);

//        void adjust_accretor_age(const real, const bool=true);
        void adjust_next_update_age();
        void update_relative_helium_mass(const real new_relative_helium_mass);

        real zeta_adiabatic();
        real zeta_thermal();
        real gyration_radius_sq();

    //Small envelope behaviour
    real small_envelope_core_luminosity();
    real small_envelope_core_radius(const real m_core);
    real small_envelope_core_radius();
    real co_core_radius(const real m_core);
    real co_core_radius();
    
    
	friend disintegrated::disintegrated(helium_giant &);
	friend white_dwarf::white_dwarf(helium_giant &, stellar_type);
	friend neutron_star::neutron_star(helium_giant &);
	friend black_hole::black_hole(helium_giant &);
   };
#endif 		// _HELIUM_GIANT
