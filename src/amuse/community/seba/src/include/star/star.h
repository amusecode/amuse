/*
 *    star.h: base class for base class ::star
 *
 *.....................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of class base_element
 *
 *....................................................................
 */

#ifndef  _STAR
#  define  _STAR

#include  "starbase.h"
#include  "node.h"

#include "stdinc.h"
//#include "constants.h"
#include "cluster_support.h"

class single_star;

class seba_counters;	// (defined later)

/*-----------------------------------------------------------------------------
 *  star  --  a base class for all element evolution.
 *-----------------------------------------------------------------------------
 */
class  star : public starbase
    {
    protected:

        // ---------------------- Global variables! ----------------------
	// --------------- (initialized in util/hdyn_io.C) ---------------

	// Counters and diagnostics:

        static seba_counters	*sc;

    public:
   
        star(node* n);
        star(star& st);

        virtual ~star() {}

        virtual  istream& scan_star_story(istream&);
        virtual  ostream& print_star_story(ostream&,
					   int short_output = 0);

        bool  is_binary();
        bool  is_binary_component();
        bool  is_star_in_binary();
        star* get_binary();
        star* get_companion();
        star* get_companion(star*);
        star* get_primary();
        star* get_secondary();
        star* get_initial_primary();
        star* get_initial_secondary();

        real  get_m_conv_star_to_dyn()      {return m_conv_star_to_dyn;}
        real  get_r_conv_star_to_dyn()      {return r_conv_star_to_dyn;}
        real  get_t_conv_star_to_dyn()      {return t_conv_star_to_dyn;}

        void  set_m_conv_star_to_dyn(const real mf)
            {m_conv_star_to_dyn = mf;}
        void  set_r_conv_star_to_dyn(const real rf)
            {r_conv_star_to_dyn = rf;}
        void  set_t_conv_star_to_dyn(const real tf)
            {t_conv_star_to_dyn = tf;}


        inline seba_counters *get_seba_counters()  {return sc;}
        void set_seba_counters(seba_counters *s){sc = s;}

        virtual void  evolve_element(const real) {}
        virtual void instantaneous_element() {}

        virtual int no_of_elements() {return 0;}
        virtual int   get_identity(){return 0;}
    	virtual stellar_type get_element_type(){return NAS;}
    	virtual binary_type get_bin_type(){return Unknown_Binary_Type;}
        virtual void  set_identity(const int){}

    	virtual real  get_current_time(){return 0;}
        virtual real  get_relative_age(){return 0;}
    	virtual real  get_previous_current_time(){return 0;}
    	virtual real  get_previous_total_mass(){return 0;} //ss.h and s.h
        virtual void  set_relative_age(const real){}
    	virtual void  set_current_time(const real){}

        virtual star_type_spec get_spec_type(star_type_spec) {return NAC;}
        virtual void set_spec_type(star_type_spec, bool on = true) {}
	virtual void set_first_contact(bool) {}

    	virtual real  get_total_mass(){return 0;}
    	virtual real  get_envelope_mass(){return 0;}
	virtual real  get_radius(){return 0;}
	virtual real  get_effective_radius(){return 0;}
        virtual real  get_luminosity(){return 0;}
        virtual void  set_luminosity(real){}
	virtual real  get_magnetic_field(){return 0;}
	virtual void  set_magnetic_field(real){}
	virtual real  get_rotation_period(){return 0;}
	virtual void  set_rotation_period(real){}

        virtual mass_transfer_type get_current_mass_transfer_type() {
	  return Unknown;} 
	
        virtual void  set_velocity(const real){}
        virtual void  set_anomal_velocity(const vec){}
	virtual void  set_effective_radius(const real){}

        virtual void  set_previous_radius(const real){}
        virtual void  set_rotational_velocity(const real) {} 

        virtual real  get_evolve_timestep(){return 0;}

        virtual void  initialize(int, real, real, real, real){}
        virtual void  initialize(double_init&, const int){}

        virtual void  clean(){}

        virtual bool  remnant(){return false;}
        virtual bool  magnetic(){return false;}
        virtual bool  hydrogen_envelope_star() {return true;}
        virtual bool  giant_star()             {return false;}
	

        virtual void refresh_memory(){}
        virtual void recall_memory(){}

        virtual void put_element(){}
        virtual void read_element(){}
 
        virtual real sudden_mass_loss() {return 0;}

        virtual void adjust_binary_after_wind_loss(star *, 
				const real, const real) {}
        virtual real get_semi(){return 0;}
        virtual real get_eccentricity(){return 0;}
        virtual void set_semi(real){}
        virtual void set_eccentricity(real){}
        virtual real get_donor_timescale(){return 0;}
        virtual void calculate_velocities(){}
        virtual real get_period() {return 0;}
        virtual void set_bin_type(binary_type){}
        virtual void post_sn_companion_velocity(star*, const real){}

//		star virtual functions
        virtual real get_relative_mass(){return 0;}
        virtual real get_core_mass(){return 0;}
        virtual real get_COcore_mass(){return 0;}
        virtual real get_previous_radius(){return 0;}
        virtual real get_core_radius(){return 0;}
        virtual void set_core_mass(const real){}
        virtual void set_COcore_mass(const real){}
        virtual void set_envelope_mass(const real){}
        virtual void set_relative_mass(const real){}
        virtual real zeta(star*, star*){return 0;}
        virtual real zeta_adiabatic(){return 0;}
        virtual real zeta_thermal(){return 0;}
        virtual real angular_momentum() {return 0;}
        virtual void stellar_wind(const real){}
	virtual void update_wind_constant(){}

        virtual real wind_velocity(){return 0;}
        virtual real get_velocity(){return 0;}
        virtual vec get_anomal_velocity(){return 0;}
        virtual void set_star_id(const int){}
	
        virtual bool low_mass_star() {return false;}
        virtual bool medium_mass_star() {return false;}
        virtual bool high_mass_star() {return false;}


//		Virtual definitions for class star
        virtual real  nucleair_evolution_timescale(){return 0;}
        virtual real  nucleair_evolution_time(){return 0;}
        virtual real  nucleair_evolution_time(const real) {
                      return 0;}
        virtual real  kelvin_helmholds_timescale(){return 0;}
        virtual real  dynamic_timescale(){return 0;}
        virtual real  main_sequence_time(){return 0;}
        virtual real  hertzsprung_gap_time(const real mass, 
					   const real z){return 0;}
        virtual real  hertzsprung_gap_time()            {return 0;}
        virtual real  helium_giant_time(const real){return 0;}
        virtual real  base_giant_branch_time(const real){return 0;}
        virtual real  base_giant_time(const real){return 0;}
        virtual void  post_constructor() {}
	

//              Luminosities
        virtual real  base_main_sequence_luminosity(){return 0;}
        virtual real  base_giant_branch_luminosity(){return 0;}
        virtual real  giant_luminosity(){return 0;}
        virtual real  base_agb_luminosity(const real){return 0;}
        virtual real  agb_luminosity(){return 0;}
        virtual real  helium_giant_luminosity(){return 0;}
        virtual real  maximum_luminosity(){return 0;}

//              Appearence.
        virtual real stellar_radius(const real, const real){return 0;}
        virtual void update() {}
        virtual real bolometric_correction() {return 0;}
        virtual real temperature() {return 1;}
        virtual real magnitude() {return 0;}
	
//		binary evolution.
        virtual real mdot_according_to_roche_radius_change(
                        star*,star*) {return 0;}
        virtual real mass_transfer_timescale(mass_transfer_type&)
	                                     {return 0;}
        virtual real accretion_limit(const real, const real) {return 0;}
        virtual star* reduce_mass(const real) {return NULL;}
        virtual star* subtrac_mass_from_donor(const real, real&){return NULL;}

	virtual void add_mass_to_core(const real) {} 

        virtual real add_mass_to_accretor(real, const real, bool) {return 0;}
        virtual real add_mass_to_accretor(const real, bool) {return 0;}
        virtual void adjust_accretor_age(const real,
					 const bool rejuvenate=true){}
        virtual void adjust_next_update_age(){}
        virtual real roche_radius(star*){return 0;}
        virtual void spiral_in(star*, star*) {}

//		Spiral in and common envelope.
        virtual star* merge_elements(star*){return 0;}
        virtual void merge_elements(star*, star*){}
        virtual real accrete_from_stellar_wind(const real, const real) 
                     {return 0;}

//		Energy
        virtual real potential_energy() {return 0;}
        virtual real kinetic_energy() {return 0;}
        virtual real total_energy() {return 0;}
        virtual real tf2_energy_diss(const real) {return 0;}
        virtual real tf3_energy_diss(const real) {return 0;}

//              output
        virtual void dump(ostream &, bool brief = true){}
        virtual void dump(char*, bool brief = true){}
        virtual void print_status(){}
        virtual void print_roche(){}
        virtual void put_state(){}
        virtual void put_hrd(ostream&){}
    
//		Tripel evolution
        virtual void adjust_triple_after_wind_loss(star *,
                                       const real, const real) {}

//              STARLAB routines.
        virtual real get_total_mass(node *) {return 0;}
        virtual void set_total_mass(node *, const real) {}

	//IO functions for starlab
	virtual void first_roche_lobe_contact_story(stellar_type) {}

	// Pure virtual functions
        virtual void detect_spectral_features()=0;
        virtual void evolve_core_mass(const real time,
				      const real mass,
				      const real z) {};
        virtual void evolve_core_mass() {};

        virtual real gyration_radius_sq()=0;

    };

#include  "seba.h"

// Shorthand for conversion from node pointer to hydro pointer:

#define N_S_PTR ((star *)n->get_starbase())

// note: automatic conversion from star to dyn scaling

inline  real get_m_conv_star_to_dyn(node * n)
    {return N_S_PTR->get_m_conv_star_to_dyn();}
inline  real get_r_conv_star_to_dyn(node * n)
    {return N_S_PTR->get_r_conv_star_to_dyn();}
inline  real get_t_conv_star_to_dyn(node * n)
    {return N_S_PTR->get_t_conv_star_to_dyn();}

// note: automatic conversion from dyn to hydro scaling

#endif 		// _STAR
