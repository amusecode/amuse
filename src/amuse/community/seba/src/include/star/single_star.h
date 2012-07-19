
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

/*
 *  single_star.h: derived class for single star base element.
 *
 *.............................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class estar
 *
 *.............................................................................
 */
#ifndef     _SINGLE_STAR
#   define  _SINGLE_STAR

#include "star.h"
#include "star_state.h"

/*-----------------------------------------------------------------------------
 *  star  --  the standard class for stellar evolution, with core and envelope
 *-----------------------------------------------------------------------------
 */
class  single_star : public star
    {
    protected:
       
        int   identity;
        stellar_type star_type; 
        star_type_spec spec_type[no_of_spec_type];

        real metalicity;

        real  current_time;
        real  relative_age;
        real  last_update_age;
        real  next_update_age;

	real  relative_mass;
        real  envelope_mass;         
        real  core_mass;         
        real  COcore_mass;

        real  radius;
        real  core_radius;
	real  effective_radius;

        real  luminosity;

        real  velocity;
        vec anomal_velocity;
	
        real birth_mass;
	real magnetic_field;            //[log G]
	real rotation_period;           //[s]
	
        real wind_constant;
        real accreted_mass;

        star_hist previous;

    public:
   
	single_star(node*);
        single_star(single_star&);
        ~single_star()	{
	  get_seba_counters()->del_sstar++;
	}
           
        int no_of_elements() {return 1;}
        star_type_spec get_spec_type(star_type_spec i) {
	  return spec_type[i];}
	
        int   get_identity()		   	{return identity;}
       
        real  get_metalicity()                  {return metalicity;}
        void  set_metalicity(const real z)      {metalicity = z;}

        real  get_previous_current_time() {return previous.current_time;}
        real  get_previous_total_mass()   {return previous.envelope_mass
                                                + previous.core_mass;}
        real  get_current_time() 		{return current_time;}
        real  get_relative_age()             {return relative_age;}
	real  get_effective_radius()         {return effective_radius;}
    	real  get_last_update_age()        {return last_update_age;}
    	real  get_next_update_age()        {return next_update_age;}
    	real  get_relative_mass()          {return relative_mass;}
    	real  get_envelope_mass()          {return envelope_mass;}
    	real  get_core_mass()              {return core_mass;}
    	real  get_COcore_mass()              {return COcore_mass;}
    	real  get_core_radius()              {return core_radius;}
	real  get_radius()                {return radius;}
        real  get_luminosity()		   {return luminosity;}
        real  get_velocity()	   	{return velocity;}
        vec  get_anomal_velocity()  	{return anomal_velocity;}
	real  get_magnetic_field()      {return magnetic_field;}
	real  get_rotation_period()     {return rotation_period;}
        real  get_total_mass()		{return envelope_mass + core_mass;}

        void  set_current_time(real t) 		{current_time=t;}
        void  set_relative_age(real t)          {relative_age=t;}

        void  set_luminosity(real l)		{luminosity = l;}
	void  set_magnetic_field(real b){magnetic_field=b;}
	void  set_rotation_period(real p){rotation_period=p;}
        void  set_identity(int i)		{identity=i;}
        void  set_envelope_mass(const real m) {envelope_mass = m;}
        void  set_core_mass(const real m) {core_mass = m;}
        void  set_COcore_mass(const real m) {COcore_mass = m;}
        void  set_spec_type(star_type_spec s, bool on=true);
	void  set_effective_radius(const real r)  {effective_radius=r;}
	void  set_last_update_age(const real t) {last_update_age = t;}
	void  set_next_update_age(const real t) {next_update_age = t;}
        void  set_previous_radius(const real r)	{previous.radius = r;}
        void  set_velocity(const real v)	{velocity = v;}
        void  set_anomal_velocity(const vec v)	{anomal_velocity = v;}

	real  magnitude();              
    	real  temperature();            

//        void clean();

      // Replaced see below
      //        bool low_mass_star();
      //        bool medium_mass_star();
      //        bool high_mass_star();
	
        void initialize(int id, real z, real t_cur,
			real t_rel, real m_rel, real m_tot,
			real m_core, real co_core);

	
//	      Time scales
        real  nucleair_evolution_timescale();
        real  kelvin_helmholds_timescale();
        real  dynamic_timescale();
        real  nucleair_evolution_time();
        real  nucleair_evolution_time(const real, const real, const real);
        real  get_evolve_timestep();
        real mass_transfer_timescale(mass_transfer_type &type);

//		Luminosities
        real  bolometric_correction();
        bool  remnant() {return false;}   // default: not a remnant.
        bool  magnetic() {return false;}  // default: no magnetic stellar wind.
        bool  hydrogen_envelope_star() {return true;}
        bool  giant_star()             {return false;}
        bool  star_with_COcore()       {return false;}   
     
        // real helium_core_radius();
 
        void refresh_memory();
        void recall_memory();

//	     stellar wind routines.
        void update_wind_constant();
        void stellar_wind(const real);
        real wind_velocity();
        real accrete_from_stellar_wind(const real, const real);
	
//		Stability rourines.
        real zeta_adiabatic();
        real zeta_thermal();
        real angular_momentum();

//           radius change for donors/accretors
	void adjust_donor_radius(const real);
        void adjust_accretor_radius(const real, const real);

//		Mass transfer routines.
	void add_mass_to_core(const real);

        real  add_mass_to_accretor(const real, bool );
        real  add_mass_to_accretor(real, const real, bool);
//        star* reduce_mass(const real);
        real  rejuvenation_fraction(const real);
	void  update_relative_mass(const real);
        void  lose_envelope_decent();
        star*  merge_elements(star*);
	
        real mass_ratio_mdot_limit(real);
        real accretion_limit(const real, const real);
        real expansionA(const real);
        real expansionB(const real);

        void update();
        void detect_spectral_features();

//           Input/Output
        void read_element();
        void put_element();
        void dump(ostream&, bool brief = true);
        void dump(char*, bool brief = true);
        void print_status();
        void print_roche();
        void put_state();
        void put_hrd(ostream &);


        real tf2_energy_diss(const real);
        real tf3_energy_diss(const real);

	real potential_energy();
	real kinetic_energy();
	real total_energy();
	void post_constructor();

        void star_transformation_story(stellar_type);
        void post_supernova_story();
        void first_roche_lobe_contact_story(stellar_type);

	virtual  istream& scan_star_story(istream&);
        virtual  ostream& print_star_story(ostream&,
					   int short_output = 0);

    real linear_function_inversion(real (single_star::*fptr)(real, const real),
                        const real x_guess, const real y_value, const real z = 0, 
                                   const real xmin = cnsts.parameters(minimum_main_sequence), 
                                   const real xmax = cnsts.parameters(maximum_main_sequence));
      real update_core_and_envelope_mass(const real m_core);
      real update_core_and_envelope_mass_TPAGB(const real m_core);
      real update_COcore_mass(const real mco_core);
//      real get_relative_mass_from_core_mass(char * input_filename, const real m_c,
//                                    const real m_r, const real z);
      real get_relative_mass_from_core_mass(const char * input_filename, const real m_c,
                                                    const real m_r, const real z);
        real get_relative_mass_from_table(const char * input_filename, const real m_c, const real m_r, const real z);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // The new metalicity dependencies come here:

      real get_zeta(const real z);
      real helium_flash_mass(const real z);
      real helium_ignition_mass(const real z);
      real helium_ignition_luminosity(const real mass, 
				      const real z);
      real helium_ignition_radius(const real mass, const real mass_tot, const real z);
      real helium_ignition_core_mass(const real mass, 
				     const real z); // Eq.67

      bool low_mass_star(const real mass, 
			 const real z);
      bool low_mass_star(); 
      bool medium_mass_star();
      bool intermediate_mass_star(const real mass,
				  const real z);
      bool high_mass_star(const real mass,
			  const real z);
      bool high_mass_star();
      real base_giant_branch_time(const real mass,
				  const real z);
      //      real base_giant_branch_time(const real mass);
      real base_giant_branch_time();
      real base_giant_branch_luminosity(const real mass,
					const real z);
      real base_giant_branch_luminosity(const real mass);
      real base_giant_branch_luminosity();
      real convective_envelope_mass(const real z);
      real FGB_mass(const real z);
      real get_hydrogen_fraction(const real z);


      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // These functions are part of main_sequence.C
      real main_sequence_hook_mass(const real z);
      //      friend main_sequence::main_sequence_hook_mass(const real z);

      //  real helium_flash_mass(const real z);
      real main_sequence_time(const real mass,
			      const real z);
      real main_sequence_time();
      real main_sequence_hook_time(const real mass, 
				   const real z);
      real stars_without_main_sequence_hook(const real mass,
					    const real z); 
      real stars_with_main_sequence_hook(const real mass,
					 const real z); 
      real terminal_main_sequence_luminosity(const real mass,
					     const real z);  
      real terminal_main_sequence_radius(const real mass,
					 const real z);  
      //real base_main_sequence_radius(const real mass);
      real base_main_sequence_radius(const real mass, const real z);

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Hertzsprung-gap
        real initial_hertzsprung_gap_core_mass(const real mass, const real z);

        
        real terminal_hertzsprung_gap_core_mass(const real mass, 
                                                const real z);

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //Functions for core helium burning
      real giant_branch_radius(const real l_hb,
			       const real mass_tot, 
			       const real z);
        real minimum_blue_loop_radius(const real mass, const real mass_tot, const real z);
        real base_horizontal_branch_luminosity(const real mass, const real z);
        real minimum_horizontal_branch_luminosity(const real mass, const real z);
        
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Core helium burning without envelope (helium star)
      real helium_main_sequence_time_for_solar_metalicity(const real mass);
      real helium_star_luminosity_for_solar_metalicity(const real mass);
      real helium_star_radius_for_solar_metalicity(const real mass_tot);
      real helium_main_sequence_luminosity(const real time, const real mass);
      real helium_main_sequence_radius(const real time, const real mass, const real mass_tot);
      real terminal_helium_main_sequence_luminosity(const real mass);
      real helium_giant_luminosity_from_core_mass(
                                const real m_core, const real mass, const real z);
      real helium_giant_radius(const real lum, const real mass, 
                                 const real mass_tot, const real z);
      real helium_giant_x_mass(const real mass);
      real helium_giant_x_luminosity(const real mass);
      real helium_giant_B_factor();
      real helium_giant_D_factor(const real mass);
      real helium_giant_p_parameter();
      real helium_giant_q_parameter();
      real helium_giant_initial_core_mass(const real mass, const real z = 0);
    
     
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Small envelope behaviour
      
      void small_envelope_perturbation();
      real small_envelope_mu(const real lum, const real mass_tot, 
                                const real m_core);
      real perturb_luminosity(const real lum, const real lum_c, 
                                const real mass_tot, const real m_core, const real mu);
      real perturb_radius(const real rad, const real rad_c, const real mass_tot, 
                                const real m_core, const real mu);
      real s_helper(const real mass_tot, const real mu);
      real r_helper(const real rad, const real rad_c, const real mass_tot, 
                                const real mu);

        virtual real small_envelope_core_luminosity(){};
        virtual real small_envelope_core_radius(){};
        virtual real helium_core_radius(){};
        
      //white dwarf      
      real white_dwarf_radius(real mass, real time);
  
      // for virtual declaration in 
      //      void evolve_core_mass(const real) {}

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Shell hydrogen burning with envelope (sub giant)
      real hydrogen_rate_constant(const real mass); 
      real sub_giant_Ah_estimator(const real mass); 
      real sub_giant_B_factor(const real mass); 
      real sub_giant_D_factor(const real mass); 
      real sub_giant_D_factor(const real mass, 
			      const real z);
      real sub_giant_p_parameter(const real mass, 
				 const real z);
      real sub_giant_q_parameter(const real mass, 
				 const real z);
//      real sub_giant_core_mass(const real time,
//			       const real mass, 
//      			       const real z);
      real determine_core_mass(const real time,
			       const real mass, 
			       const real z, 
			       const real A,  
			       const real t_b,  
			       const real l_b);

      real determine_age(const real m_core, const real mass,
                         const real z, const real A, 
                         const real t_b, const real l_b);
            
      real FGB_x_luminosity(const real mass, const real z);
      real FGB_x_mass(const real mass, const real z);
      real FGB_core_mass_luminosity_relation(const real lum,
					     const real mass,
                                             const real z); 
      //real FGB_core_mass_luminosity_relation(const real time, 
	//				     const real mass, 
	//				     const real z);
      real FGB_luminosity_core_mass_relation(const real time, 
					     const real mass, 
					     const real z);
      real FGB_luminosity_core_mass_relation(const real m_core, 
					     const real mass);

//      real sub_giant_luminosity(const real time,
//				const real mass, 
//				const real z);
//      real specific_time_x_factor(const real mass,
//				  const real z);

      real specific_time_boundary(const real mass,
                                    const real A,
                                    const real t_b,
                                    const real l_b, 
                                    const real D,
                                    const real p,
                                    const real l_x);
        
  //      real specific_time_boundary(const real mass,
  //				  const real z);
      real specific_time_limit(const real A,
			       const real t_min,
			       const real DB,
			       const real l_min,
			       const real pq);

      //      real specific_upper_time_limit(const real mass,
      //				     const real z,
      //           real (single_star::*fptr_time)(const real, const real),
      //           real (single_star::*fptr_luminosity)(const real, const real));

      real helium_ignition_time(const real mass, 
				const real z);
      real helium_ignition_time();

      real base_giant_branch_core_mass(const real mass, 
      				       const real z);
      real base_horizontal_branch_radius(const real mass, 
					 const real mass_tot, const real z);
      real core_helium_burning_timescale(const real mass, 
					 const real z); //Eq.57
      real core_helium_burning_timescale();

      real base_AGB_luminosity(const real mass, const real z); // Eq.56
      real base_AGB_time(const real mass, const real z);
      real base_AGB_time();

      real AGB_A_He_estimator(); //Eq.68
      real AGB_radius(const real lum, const real mass, const real mass_tot, const real z);

      real base_AGB_core_mass(const real mass, const real z); // Eq.66
      real base_AGB_relative_mass(const real m_core, const real z); // Eq.66 inverted
      real TAGB_time(const real mass, const real mass_tot, const real z);
      real dredge_up_core_mass(const real mass, const real z);
      real TPAGB_AH_He_estimator();
      real maximum_AGB_core_mass(const real mass,
                               const real z); //Eq.75
      real AGB_luminosity(const real CO_core_mass, const real mass, 
                        const real z); //Eq.37
      real dredge_up_time(const real mass, const real z);// Eq.70
      real dredge_up_luminosity(const real mass, const real z);
  
    };

// Shorthand for conversion from node pointer to star pointer:

#define N_SS_PTR ((single_star *)n->get_starbase())

// note: automatic conversion from star to dyn scaling

node* mkstar(int, real, stellar_type type=Main_Sequence);

void extract_story_chapter(stellar_type&, real& z, 
			   real& t_cur, real& t_ral,
                           real&, real&, real&, 
                           real&, real&,
			   real&, real&, real&,
			   story&);
void extract_line_text(stellar_type&, real& z, real&, real&, real&,
                       real&, real&, real&, real&, real&,
		       real&, real&, story&);

bool merge_elements(single_star* primary, single_star *secondary);

void  addstar(node*,
	      real t_rel=0,
	      stellar_type type=Main_Sequence,
	      real z = cnsts.parameters(Zsun),
              int id = 1,
	      bool verbose = false);


#endif		// _SINGLE_STAR








