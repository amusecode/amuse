/*
 * star_cluster.h: derived class for evolution of sero radius stars 
 *                which lose mass proportional to wind_constant.
 *                This definition is used for the ski-jump problem.
 *
 *.....................................................................
 *    version 1:  April 2005   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class star_cluster
 *
 *....................................................................
 */

#ifndef   _STAR_CLUSTER
#   define _STAR_CLUSTER

#include "single_star.h"

#define Rtidal_over_Rvir_KingW09 8.38684

/*-----------------------------------------------------------------------------
 *  star_cluster  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class star_cluster : public single_star {
      private:

        real black_hole_mass;
	real trlx_initial;

      protected:
	real rvir;
	real nstar;
	real x_imf;
	real minimum_mass;

      public :

         star_cluster(node* n) : single_star(n) { }
         ~star_cluster() {}

         stellar_type get_element_type() {return Star_Cluster;}
         bool hydrogen_envelope_star() {return false;}

         void instantaneous_element();
         void evolve_element(const real);
         void update();
         void update_wind_constant(const real);
         void stellar_wind(const real dt);

//            Mass transfer utilities.
	 star* subtrac_mass_from_donor(const real, real&);
         real accretion_limit(const real, const real);

         void adjust_accretor_age(const real, const bool);
         void adjust_next_update_age();
         star* reduce_mass(const real);

//              Stability rourines.
         real zeta_thermal();
	 real gyration_radius_sq();

//              Specific star cluster functions
	 void initialize_star_cluster();
	 real mean_stellar_mass(const real t, const real mmin);
	 real mean_stellar_mass();
	 real mass_loss_by_evolution(const real dt);
	 real turn_off_mass(const real t);
	 real tidal_radius();
	 real relaxation_time(real n, real m, real r);

	 void print(ostream& s);

    };
#endif 		// _STAR_CLUSTER


