
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

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
        ~single_star()	{}
           
        int no_of_elements() {return 1;}
        star_type_spec get_spec_type(star_type_spec i) {
	  return spec_type[i];}
	
        int   get_identity()		   	{return identity;}
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

        bool low_mass_star();
        bool medium_mass_star();
        bool high_mass_star();
	
        void initialize(int, real, real, real, real, real, real);
	
//	      Time scales
        real  nucleair_evolution_timescale();
        real  kelvin_helmholds_timescale();
        real  dynamic_timescale();
        real  main_sequence_time();
        real  main_sequence_time(const real);
        real  hertzsprung_gap_time(const real);
        real  hertzsprung_gap_time(const real, const real);
        real  helium_giant_time(const real);
        real  helium_giant_time(const real, const real);
        real  nucleair_evolution_time();
        real  nucleair_evolution_time(const real);
        real  base_giant_branch_time(const real);
        real  base_giant_branch_time(const real, const real);
        real  base_giant_time(const real);
        real  base_giant_time(const real, const real);
        real  helium_time();
        real  get_evolve_timestep();

        real mass_transfer_timescale(mass_transfer_type &type);

//		Luminosities
        real  base_main_sequence_luminosity(const real);
        real  base_main_sequence_luminosity();
        real  base_giant_branch_luminosity();
        real  base_giant_branch_luminosity(const real);
        real  giant_luminosity();
        real  giant_luminosity(const real);
        real  base_agb_luminosity(const real);
        real  base_agb_luminosity(const real, const real);
        real  agb_luminosity();
        real  agb_luminosity(const real);
        real  helium_giant_luminosity();
        real  helium_giant_luminosity(const real);
        real  maximum_luminosity();
        real  maximum_luminosity(const real);

        real  bolometric_correction();
        bool  remnant() {return false;}   // default: not a remnant.
        bool  magnetic() {return false;}  // default: no magnetic stellar wind.
        bool  hydrogen_envelope_star() {return true;}
        bool  giant_star()             {return false;}
        bool  star_with_COcore()       {return false;}   
     
        real helium_core_radius();
        void evolve_core_mass(const real) {}
        real final_core_mass();

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

        real  add_mass_to_accretor(const real);
        real  add_mass_to_accretor(real, const real);
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
        void dump(const char*, bool brief = true);
        void print_status();
        void print_roche();
        void put_state();
        void put_hrd(ostream &);

	real get_dlogR_dT();

        real tf2_energy_diss(const real);
        real tf3_energy_diss(const real);

	real potential_energy();
	real kinetic_energy();
	real total_energy();
	void post_constructor();


        void star_transformation_story(stellar_type);
	void merge_two_stars_story(stellar_type);
        void post_supernova_story();
        void first_roche_lobe_contact_story(stellar_type);

	virtual  istream& scan_star_story(istream&);
        virtual  ostream& print_star_story(ostream&,
					   int short_output = 0);
	
    };

// Shorthand for conversion from node pointer to star pointer:

#define N_SS_PTR ((single_star *)n->get_starbase())

// note: automatic conversion from star to dyn scaling

node* mkstar(int, real, stellar_type type=Main_Sequence);

void extract_story_chapter(stellar_type&, real&, real&,
                           real&, real&, real&, 
                           real&, real&,
			   real&, real&, real&,
			   story&);
void extract_line_text(stellar_type&, real&, real&, real&,
                       real&, real&, real&, real&, real&,
		       real&, real&, story&);

bool merge_elements(single_star* primary, single_star *secondary);

void  addstar(node*,
	      real t_rel=0,
	      stellar_type type=Main_Sequence,
	      bool verbose = false);

void  addstar1(node * b, real t_rel, stellar_type type,
	       real mf=1, real rf=1, real tf=1);


//real get_stellar_radius(const real mass, const real time, stellar_type& type=Main_Sequence);
#endif		// _SINGLE_STAR








