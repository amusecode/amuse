/*
 * giant.h: derived class for evolution of stars
 *    in the core helium burning stage.
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class giant
 *
 *....................................................................
 */

#ifndef    _GIANT 
#   define _GIANT

#include "single_star.h"
#include "horizontal_branch.h"
#include "super_giant.h"
#include "helium_star.h"

		// Known class declarations.
class sub_giant; 

/*-----------------------------------------------------------------------------
 *  giant  --  a derived class for element evolution.
 *-----------------------------------------------------------------------------
 */
class giant : public single_star {
      protected :
//         static ioserror error;
      public :
		// Constructors.
//         giant() : single_star() {}
         giant(sub_giant &);
         giant(node* n) : single_star(n) {}
		// Destructors.
         ~giant() {}
//              Member acces functions.
        stellar_type get_element_type() {return Giant;}

		// Member function definition.
        void adjust_initial_star();
         void evolve_element(const real);
         real stellar_radius(const real, const real);
//              Mass transfer utilities.
        star* reduce_donor_mass(const real);
        star* subtrac_mass_from_donor(const real, real&);
        void adjust_accretor_age(const real);
        void adjust_next_update_age();

//		Spiral in and common envelope.

		// Debugging utilities.
        	// Friend functions.
	 friend horizontal_branch::horizontal_branch(giant &);
	 friend super_giant::super_giant(giant &);
	 friend helium_star::helium_star(giant &);
      };
#endif 		// _GIANT
