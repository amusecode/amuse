/*
 * proto_star.h: derived class for evolution of proto stars.
 *
 *..........................................................................
 *    version 1:  May 1999   Simon F. Portegies Zwart
 *    version 2:
 *..........................................................................
 *     This file includes:
 *  1) definition of class proto_star
 *
 *..........................................................................
 */

#ifndef    _PROTO_STAR
#   define _PROTO_STAR

#include "single_star.h"
#include "brown_dwarf.h"
#include "main_sequence.h"

		// Known class declarations.

class brown_dwarf;
class main_sequence;
/*-----------------------------------------------------------------------------
 *  proto_star  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class proto_star : public single_star {
      private:
   
         vec angular_momentum;

      public :

        proto_star(node* n) : single_star(n) {}
        ~proto_star() {}
	
        stellar_type get_element_type() {return Proto_Star;}
        bool giant_star()             {return true;}

	void instantaneous_element();
	void evolve_element(const real);
        void create_zero_age_object();
        void create_binary_from_proto_star();

//           Mass transfer stability.
        real zeta_thermal();
	real gyration_radius_sq();
	
        star* reduce_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
        void adjust_accretor_age(const real, const bool);
        void adjust_next_update_age();
        void update_wind_constant();

        void stellar_wind(const real);
        real helium_core_mass();
       
	vec get_angular_momentum()         {return angular_momentum;}
	void   set_angular_momentum(vec v) {angular_momentum = v;}

//              Friend constructors
//	 friend planet::planet(proto_star &);
	 friend brown_dwarf::brown_dwarf(proto_star &);
	 friend main_sequence::main_sequence(proto_star &);

      };
#endif 		// _PROTO_STAR



