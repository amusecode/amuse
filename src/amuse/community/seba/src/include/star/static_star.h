/*
 * static_star.h: derived class for evolution of zero radius stars 
 *                which lose mass proportional to wind_constant.
 *                This definition is used for the ski-jump problem.
 *
 *.....................................................................
 *    version 1:  Jan 1999   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class static_star
 *
 *....................................................................
 */

#ifndef   _static_star 
#   define _static_star 

#include "single_star.h"

/*-----------------------------------------------------------------------------
 *  static_star  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class static_star : public single_star {
      private:

      public :

         static_star(node* n) : single_star(n) {}
         ~static_star() {}

	 
         stellar_type get_element_type() {return Static_Star;}
         bool hydrogen_envelope_star() {return false;}

         void instantaneous_element();
         void evolve_element(const real);
         void update();
      
         real gyration_radius_sq() {return 0;}

    };
#endif 		// _static_star 


