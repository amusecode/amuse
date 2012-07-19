
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\
 //                                                       //            _\|/_
//=======================================================//              /|\

/*
 *  star_state.h: derived class for element evolution systems.
 *          functions as derived class for the real elements.
 *.............................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class star_support
 *
 *.............................................................................
 */
#ifndef     _STAR_STATE
#   define  _STAR_STATE

#include "stdinc.h"
#include "star_support.h"
#include "constants.h"

class starbase;	        // Must be known for print_star
class star; 		// Must be known for star_state
class single_star;
class node;
single_star * new_single_star(stellar_type type = Main_Sequence, 
			      int id=0, real z = cnsts.parameters(solar_metalicity),
			      real t_cur = 0, real t_rel=0,
                              real m_rel=1, real m_tot=1, 
			      real m_core=0, real co_core=0,
			      real p_rot=0, real b_fld=0,
                              node* n=0);

/*-----------------------------------------------------------------------------
 * star_hist -- the standard class for stellar evolution, with core and envelope
 *-----------------------------------------------------------------------------
 */
struct star_hist {
    public:

        stellar_type star_type;

        real  metalicity;
        real  current_time;
        real  last_update_age;
        real  next_update_age;
        real  relative_age;
        real  relative_mass;
        real  envelope_mass;
        real  core_mass;
	real  COcore_mass;
        real  radius;
        real  effective_radius;
        real  magnetic_field;
        real  rotation_period;
        real  birth_mass;

        void put_star_hist();
        void clean() {
           star_type=NAS;
	   metalicity = 0;
           current_time=last_update_age=relative_age=relative_mass
           = next_update_age=envelope_mass=core_mass=COcore_mass
           = radius=effective_radius=0;
	   magnetic_field=rotation_period=birth_mass=0;
        }
    };

/*-----------------------------------------------------------------------------
 * star_state -- state of star.
 *-----------------------------------------------------------------------------
 */

//		Was once called star_appeal.
struct star_state {
   public:

      int identity;

      stellar_type   type;
      spectral_class class_tpe;
      luminosity_class lclass;
      int class_spec[no_of_spec_type]; //keep star_type_spec's

      real metalicity;
      real mass;
      real radius;
      real velocity;
      real mdot; 
      
      star_state();
      star_state(star*);

      void put_star_state(ostream & s = cerr);

      void make_star_state(star*);
      void init_star_state(star*);
      bool special();

      void clean() {
         type=NAS;
         class_tpe=no_spectral_class;
         lclass=no_luminosity_class;
         for (int i=(int)NAC; i<(int)no_of_spec_type; i++)
             class_spec[(star_type_spec)i]=0;
	 metalicity=0;
         radius=velocity=0;
         mass=mdot=0;
      }
   };

void put_state(star_state, ostream & s=cerr);
void put_short_state(star_state, ostream & s=cerr);
char* type_dominant_state(star_state);
void print_star(starbase*, ostream & s=cerr);
void pretty_print_star(starbase*, int depth_level, ostream & s=cerr);
void pretty_print_star(starbase*, ostream & s=cerr);
bool remnant(stellar_type);


#endif 		// _STAR_STATE
