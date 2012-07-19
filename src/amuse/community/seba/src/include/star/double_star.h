
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  double_star.h: derived class for synchronous element evolution.
 *           
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class double_star
 *
 *....................................................................
 */
#ifndef  _DOUBLE_STAR 
#  define  _DOUBLE_STAR

#include "double_support.h"
#include "main_sequence.h"
#include "single_star.h"

#define MAX_CIRC_ECC 1.e-6	// Maximum eccentricity that will be forced to
				// zero in kepler::initialize_from_pos_and_vel.
				// Consistent with value previously hard-coded
				// in kepler.C

				// Added by Steve 6/98 -- this may not be the
				// right place to define it...

class node;
class dyn;
class hdyn;
/*-----------------------------------------------------------------------------
 *  double_star  --  a derived class for base of element evolution.
 *-----------------------------------------------------------------------------
 */
class double_star : public star
    {
    protected:
        int recursive_counter;
       
        real  semi;
        real  eccentricity;

        binary_type   bin_type;

        int identity;
        real  binary_age;
        real  minimal_timestep;

        real  velocity;

        int  donor_identity;
        stellar_type donor_type;
        real donor_timescale;
	mass_transfer_type current_mass_transfer_type;

	bool first_contact;

        double_hist previous;
        double_init initial;

    private:

    public:
	
        double_star(node*);

        ~double_star() {}

        virtual  istream& scan_star_story(istream&);
        virtual  ostream& print_star_story(ostream&, int short_output = 0);

        int get_identity()		{return identity;}
        void set_identity(const int i)		{identity=i;}
        binary_type get_bin_type()	  	{return bin_type;}
	binary_type obtain_binary_type();
        real get_current_time()		{return binary_age;}
        real get_eccentricity()		{return eccentricity;}
        real get_semi()			{return semi;}
        real get_radius()               {return semi*(1-pow(eccentricity, 2));}
        real get_velocity()		{return velocity;}
	real get_effective_radius()     {
             real r_eff = semi;
             return r_eff;
        }
	mass_transfer_type get_current_mass_transfer_type()
	  {return current_mass_transfer_type;}

        void set_current_mass_transfer_type(mass_transfer_type type)	{current_mass_transfer_type = type;}

        real get_evolve_timestep();
 
        void set_bin_type(binary_type type)	{bin_type = type;}
        void set_eccentricity(real e)	{eccentricity = e;}
        void set_velocity(const real v)		{velocity = v;}
        void set_semi(real a)		{semi = a;}
        void set_current_time(const real t){binary_age = t;}
        void set_relative_age(const real t){binary_age = t;}

	void set_first_contact(bool c) {first_contact=c;}

        real mass_ratio();
        real get_total_mass() {
             return get_primary()->get_total_mass() 
    		  + get_secondary()->get_total_mass();
        }

	double_init * get_initial_conditions() {return &initial;}

        int no_of_elements() 	{
            return get_primary()->no_of_elements() 
	         + get_secondary()->no_of_elements();}

//		Initialization.
        void initialize(binary_type, real, real, real, int);

//		Timestep determination.
        void determine_minimal_timestep(); 
        real determine_dt(const real, const real);
	real internal_time_step(const real evolve_timestep);

//		Binary evolution.
        real angular_momentum();
        void circularize();
        void force_circularization();
        void recursive_binary_evolution(real, const real); 

        void instantaneous_element();
 	void evolve_element(const real);
 	void evolve_the_binary(const real);
	void try_zero_timestep();         
						
        void perform_wind_loss(const real);
        void angular_momentum_loss(const real);
        void gravrad(const real);
        real de_dt_gwr(const real);
	real gwr_angular_momentum_loss(const real m_prim,
				       const real m_sec,
				       const real sma);
	
	real mb_angular_momentum_loss();
        void magnetic_stellar_wind(const real);
        void calculate_velocities();
        void semi_detached(star*, star*, real); 
        void contact(star*, star*, real); 
        void contact_binary(real); 
        //void common_envelope();
        //void binary_in_contact(const real); 
        //bool ready_for_mass_transfer(star*);

//		Mass transfer utilities.
        void perform_mass_transfer(const real, star*, star*); 
	void angular_momentum_envelope_ejection(star*, star*);
        void dynamic_mass_transfer();
        void tidal_instability();
        real mdot_according_to_roche_radius_change(star*, star*, const real);
        void adjust_binary_after_wind_loss(star*, 
                                           const real, const real);

//		Stability of system.
        bool  stable(star* st=NULL); 
        real zeta(star*, star*, const real, bool no_aml = true);
	real gyration_radius_sq();
	
//		Stellar merging.
        void spiral_in(star*, star*); 
        void merge_elements(star*, star*);
	void double_spiral_in();
	
//		Super Nova utilities.
        void post_sn_companion_velocity(star*, const real);

        real get_period();
	real roche_radius(star *);
	real roche_radius(const real, const real, const real);
	real circularization_radius(const real m1, const real m2);

//		Mass transfer utilities.
        real get_donor_timescale() {return donor_timescale;}
        void set_donor_timescale(star*, bool first = false);
	
//		History functions.
        void refresh_memory();
        void recall_memory();
//        void update_binary_parameters();

//     		I/O functions.

	void dump(ostream &, bool brief = true);
	void dump(char*, bool);
        void put_element();
        void print_status();
        void print_roche();
        void put_state();
        void put_hrd(ostream &);

//		Energy
        real potential_energy();
        real kinetic_energy();
        real total_energy();

        void enhance_cluster_profile(cluster_profile&);

//		Triple ad multiple routines
        bool remnant() {return FALSE;}
        real wind_velocity();
        real get_core_mass() {return get_total_mass();}
        star* reduce_mass(const real);

        star* subtrac_mass_from_donor(const real, real&);
        void stellar_wind(const real);
        real temperature();
        real bolometric_correction();
        real mass_ratio_mdot_limit(real);
        real kelvin_helmholds_timescale();
        real nucleair_evolution_timescale();
        real dynamic_timescale();
        real accretion_limit(const real, const real);
        real mass_transfer_timescale(mass_transfer_type &type);
        real zeta_adiabatic();
        real zeta_thermal();
	//        real add_mass_to_accretor(const real, bool);
        real add_mass_to_accretor(real, bool, const real = -1.);
        real accrete_from_stellar_wind(const real, const real);
        void adjust_triple_after_wind_loss(star*,
                                           const real, const real);

//		Unusual specific stellar functions.
        star_type_spec get_spec_type(star_type_spec i)     {return NAC;}
        void set_spec_type(star_type_spec s, bool on=true) {}
        stellar_type get_element_type() {return Double;}
        void set_rotational_velocity(const real);
        void set_previous_radius(const real);
        real get_core_radius() {return get_radius();}

        bool low_mass_star();
        bool medium_mass_star();
        bool high_mass_star();

        void evolve_core_mass(const real) {}

        real orbital_timescale();

	void detect_spectral_features() {}
    };

double_star * new_double_star(node*,
			      real, real, real age=0, int id=0,
			      binary_type type = Detached);

void adddouble(node * b, real time, binary_type=Detached, 
               bool random_initialization=false, 
               real a_min=1, real a_max=1.e+6,
               real e_min=0, real e_max=1) ;

bool has_dstar(node*);

void add_secondary(node*, real);
void mksecondary(node*, real, real);
void mkrandom_binary(const real lm_min, const real lm_max, const real m_exp,
		     const real la_min, const real la_max,
		     const bool e_flag,
		     real &m_prim, real &m_sec, real &semi,
		                                real &ecc);

real random_exponential_mass(const real m_min,
			     const real m_max,
			     const real m_alpha);

#define  put_double_star  put_node

#endif 		// _DOUBLE_STAR
