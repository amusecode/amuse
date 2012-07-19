/*
 * SPZDCH_star.h: derived class for evolution of sero radius stars 
 *                which lose mass proportional to wind_constant.
 *                This definition is used for the ski-jump problem.
 *
 *.....................................................................
 *    version 1:  Jan 1999   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class SPZDCH_star
 *
 *....................................................................
 */

#ifndef   _SPZDCH_STAR 
#   define _SPZDCH_STAR 

#include "single_star.h"

/*-----------------------------------------------------------------------------
 *  SPZDCH_star  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class SPZDCH_star : public single_star {
      private:

      public :

         SPZDCH_star(node* n) : single_star(n) {}
         ~SPZDCH_star() {}

	 
         stellar_type get_element_type() {return SPZDCH_Star;}
         bool hydrogen_envelope_star() {return false;}

         void instantaneous_element();
         void evolve_element(const real);
         void update();
         void update_wind_constant(const real);
         void stellar_wind(const real dt);

//            Mass transfer utilities.
        real mdot_limit(const real);
	 star* subtrac_mass_from_donor(const real, real&);
         real accretion_limit(const real, const real);

         void adjust_accretor_age(const real, const bool);
         void adjust_next_update_age();
         star* reduce_mass(const real);

//              Stability rourines.
         real zeta_thermal();
	 real gyration_radius_sq();

    };
#endif 		// _SPZDCH_STAR 


