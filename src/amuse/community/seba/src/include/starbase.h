
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file starbase.h  Underlying class for star systems, on same level as node.
//
//  version 1:  Apr 1993   Piet Hut & Steve McMillan
//  version 2:
//
//  This file includes:
//  1) definition of class starbase

#ifndef  STARLAB_STARBASE_H
#  define  STARLAB_STARBASE_H

#include  "starlab_vector.h"

		    	// node.h is commented out deliberately,
			// because node.h includes starbase.h.

			// #include  "node.h"

#include  "story.h"
#include  "star/star_support.h"
#include  "star/double_support.h"

//#include  "star/seba.h"


class node;		// not included but defined as being some class
class star;

/// \a starbase:  The underlying class for stars.

class  starbase
{
    protected:

        node  * the_node;
	story * star_story;

        						
        static real  m_conv_star_to_dyn;   ///< Mass conversion factor.
        static real  r_conv_star_to_dyn;   ///< Length conversion factor.
        static real  t_conv_star_to_dyn;   ///< Time conversion factor.

        static bool  use_hdyn;		   ///< Using hdyn (kira) data.

        // static seba_counters* sbc;

    public:

        starbase(node* n=0);
        starbase(starbase& sb); 

	virtual ~starbase() {delete star_story;}

        //seba_counters* get_seba_counters() {return sbc;} 
        //void set_seba_counters(seba_counters *sb) {sbc=sb;} 
	    
        node  *get_node()		       {return the_node;}
	story *get_star_story()		       {return star_story;}

        void  set_node(node* n)		       {the_node = n;}
	void  set_star_story(story * ss)       {star_story = ss;}

	virtual ostream & print_star_story(ostream&,
					   int short_output = 0);
	virtual istream & scan_star_story(istream&, int level = 0);

	bool get_use_hdyn();
	void set_use_hdyn(bool);

	virtual void dump(ostream&, bool);
	
	//              Stellar functions
	//              REMOVE for standard STAREV
        virtual stellar_type get_element_type();
        virtual real get_total_mass();
        virtual real get_effective_radius();
        virtual real get_current_time();
        virtual real get_relative_age();
	
        virtual real temperature();
        virtual real get_luminosity();

        virtual vec get_anomal_velocity();
        virtual void set_anomal_velocity(const vec);
        virtual void evolve_element(const real);
        virtual real get_evolve_timestep();
        virtual star* merge_elements(star*);
        
        virtual real sudden_mass_loss();

        virtual binary_type get_bin_type();
        virtual real get_semi();
        virtual void set_semi(real);
        virtual real get_eccentricity();
        virtual void set_eccentricity(real);
        
// AMUSE
        
        virtual real get_time_offset();
        virtual void set_time_offset(real value);

//	Scaling:

        void set_stellar_evolution_scaling(real, real, real);
        bool get_stellar_evolution_scaling();
        void print_stellar_evolution_scaling(ostream&);

        real conv_m_star_to_dyn(real);
        real conv_r_star_to_dyn(real);
        real conv_t_star_to_dyn(real);

        real conv_m_dyn_to_star(real);
        real conv_r_dyn_to_star(real);
        real conv_t_dyn_to_star(real);

};

typedef  starbase *(*sbpfp)();

inline  starbase * new_starbase()    {return  new starbase;}

// Initialization functions for mass function, 
//                              mass ratio distribution,
//                              eccentricity distribution,
//                              semi-major axis distribution.

// Mass function.
enum mass_function {Unknown_MF=-1, 
		    Equal_Mass, mf_Power_Law, Miller_Scalo, Scalo, Kroupa,
		    GdeMarchi, KTG91, TwoComponent};

real get_random_stellar_mass(real m_lower, real m_upper, 
			     mass_function mf, real exponent);
char* type_string(mass_function mf);
mass_function extract_mass_function_type_string(char* type_string);
real general_power_law(real lowerl, real upperl, real exponent);

// Mass ratio distribution
enum mass_ratio_distribution {Unknown_qf=-1,
			      Equal_q, Flat_q, qf_Power_Law, Hogeveen};
char* type_string(mass_ratio_distribution qf);
mass_ratio_distribution 
    extract_mass_ratio_distribution_type_string(char* type_string);
real get_random_mass_ratio(real q_lower, real q_upper, 
			   mass_ratio_distribution qf, 
			   real exponent);
			      
enum sma_distribution {Unknown_smaf=-1, 
		       Equal_sma, sma_Power_Law, Duquennoy_Mayor, Eggleton};
char* type_string(sma_distribution smaf);
sma_distribution 
    extract_semimajor_axis_distribution_type_string(char* type_string);
real get_random_semimajor_axis(real a_lower, real a_upper, 
			       sma_distribution smaf, 
			       real exponent, real m_prim, real m_sec);

enum ecc_distribution {Unknown_eccf, 
		       Equal_ecc, ecc_Power_Law, Thermal_Distribution};
char* type_string(ecc_distribution eccf);
ecc_distribution 
    extract_eccentricity_distribution_type_string(char* type_string);
real get_random_eccentricity(real e_lower, real e_upper, 
			     ecc_distribution eccf, 
			       real exponent);

void mkrandom_binary(real m_min,  real m_max,
		     mass_function mf,  real m_exp,
		     real q_min,  real q_max,
		     mass_ratio_distribution qf,  real q_exp,
		     real a_min,  real a_max,
		     sma_distribution af,  real a_exp,
		     real e_min,  real e_max,
		     ecc_distribution ef,  real e_exp,
		     real &m_prim, real &m_sec, real &semi,
		     real &ecc, real z);

void print_initial_binary_distributions(real m_min,  real m_max,
					mass_function mf,  real m_exp,
					real q_min,  real q_max,
					mass_ratio_distribution qf,  
					real q_exp,
					real a_min,  real a_max,
					sma_distribution af,  real a_exp,
					real e_min,  real e_max,
					ecc_distribution ef,  real e_exp);


#endif
 
