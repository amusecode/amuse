
       //=======================================================//   _\|/_
      //  __  _____           ___                    ___       //     /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //         _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //           /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


/// @file dyn.h  Base class for N-body systems.

//  version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
//  version 2:
//
//  This file includes:
//  1) definition of class dyn

#ifndef  STARLAB_DYN_H
#  define  STARLAB_DYN_H

#include <vector>
#include  "util_math.h"
#include  "node.h"
#include  "kepler.h"

#define BIN_INDENT	21	// for print_pert() and print_binary_params()

static const real rnull = 0.0;

//#include  "dyn_kepler.h"   NOTE: this file is included at the end instead,
//                           because it needs to know about the dyn class

/// \a dyn: The simplest class of dynamical particles.

class  dyn : public node
{
    protected:

        // Global variables:

    	static xreal system_time;	///< Extended-precision time.
    	static real real_system_time;	///< Real copy of system_time.

    	static bool use_sstar;  	///< Activate single star evolution.

	/// Flag to neglect all internal interactions.

	static bool ignore_internal;

	//------------------------------------------------------------

	/// External field specification.

	// (Moved from hdyn to dyn by Steve, 7/01):

	static unsigned int external_field;
					// 0 ==> none
					// 1 ==> tidal (bit 0)
					// 2 ==> Plummer (bit 1)
					// 3 ==> Plummer+tidal (bits 0+1)
					// 4 ==> power-law (bit 2)
					// etc.

	/// Tidal field specification.

	/// 0 => no tidal field <br>
	/// 1 => point-mass field <br>
	/// 2 => isothermal halo field <br>
	/// 3 => disk field (Oort constants)

	// (Changed by Steve, 6/99)
	// Note that alpha and omega are in general independent.

	static int  tidal_type;

        static real omega;		///< System (circular) angular speed.
        static real omega_sq;		// omega squared (probably not needed)

	static real alpha1;		///< Tidal x-acceleration is -alpha1*x.
	static real alpha3;		///< Tidal x-acceleration is alpha3*z.

	static vec tidal_center;	///< (Fixed) center of tidal force.
					//				(7/01)

	// Confining Plummer model parameters (Steve, 7/01):

	static real p_mass;		///< Total Plummer-sphere mass.
	static real p_scale_sq;		///< Plummer scale radius squared.

	static vec p_center;		///< (Fixed) center of Plummer forces.
	static bool p_friction;		///< Flag to include dynamical friction.
					// (currently for Plummer only)

	// Confining power-law model parameters (Steve, 7/01).
	// New implementation (Steve, 2/04) means that Plummer is
	// no longer a special case (with exponent = 0).

	static real pl_coeff;		///< Overall power-law coefficient.
	static real pl_scale;		///< Power-law scale/cutoff radius.
	static real pl_exponent;	///< Power-law exponent.

	static vec pl_center;		///< (Fixed) center of power-law forces.

	static FILE *ifp;  ///< If set, use file instead of C++ stream input.
	static FILE *ofp;  ///< If set, use file instead of C++ stream ouput.

	static bool col_output;  	///< If true, output is in col format.

	//------------------------------------------------------------

	vec  pos;	         	///< Position (3-D Cartesian vector).
	vec  vel;	         	///< Velocity: (d/dt) pos.
	vec  acc;	         	///< Acceleration: (d/dt) vel.
	kepler * kep;		 	///< Pointer to a kepler orbit object.

    public:

	static FILE* set_ifp(FILE* p)		{return ifp = p;}
	static FILE* get_ifp()			{return ifp;}
	void clear_ifp()			{ifp = NULL;}
	static FILE* set_ofp(FILE* p)		{return ofp = p;}
	static FILE* get_ofp()			{return ofp;}
	void clear_ofp()			{ofp = NULL;}

	/// Col output is a simpler data format suitable for flat trees.

	static bool set_col_output(bool b = true)	{return col_output = b;}
	static bool get_col_output()			{return col_output;}

	/// Reset dynamical quantities.

	inline void dyn_init() {
	    pos = vel = acc = 0.0;
	    kep = NULL;
	}

        dyn(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase,
	    bool use_stories = true)
	   : node(the_hbpfp, the_sbpfp, use_stories)	{dyn_init();}
 
	/// Delete kepler structure and pointer.

	inline void rmkepler() {
	    if (kep) {
		dyn* s = get_binary_sister();

		if (s && s->kep == kep)
		    s->kep = NULL;

		// changed (unsigned int) to unsigned long
		if ((unsigned long) kep != 1 && (unsigned long) kep != 2)
		    delete kep;		// in case kep is used as a flag

		kep = NULL;
	    }
	}

	virtual ~dyn() {

	    // Note added by Steve 7/7/98:
	    //
	    // Some care is required when deleting the kepler structure,
	    // as binary components in kira share a kepler.  In that case,
	    // when the kepler of the first component is deleted, the
	    // kepler pointer of the other component must also be set
	    // NULL.  Otherwise, an attempt will be made to delete a
	    // nonexistent structure.  Do this via another function so
	    // that the same procedure can be used from other classes.

	    rmkepler();
	}

	inline xreal get_system_time()	const	{return system_time;}
	inline real get_real_system_time() const {return real_system_time;}
	void set_system_time(xreal t)		{system_time = t;
					         real_system_time = t;}

	void set_ignore_internal(bool i = true)	{ignore_internal = i;}
	bool get_ignore_internal()	const	{return ignore_internal;}

	//-----------------------------------------------------------------
	// External field:
  
	inline unsigned int get_external_field() const {return external_field;}
	inline void set_external_field(unsigned int e) {external_field = e;}

	/// Tidal field (external field #1, bit 0).

	inline void set_tidal_field(int t)	{tidal_type = t;
					         if (t > 0)
						     SETBIT(external_field, 0);
					         else
					     	     CLRBIT(external_field, 0);}
	inline int get_tidal_field()	const	{return tidal_type;}
	inline void unset_tidal_field(int t)	{set_tidal_field(0);}

        void set_omega(real o)			{omega = o;
						 omega_sq = o*o;}
        inline real get_omega()		const	{return omega;}
        inline real get_omega_sq()	const	{return omega_sq;}

	void set_alpha(real a1, real a3)	{alpha1 = a1;
						 alpha3 = a3;}
	inline real get_alpha1()	const	{return alpha1;}
	inline real get_alpha3()	const	{return alpha3;}
	inline void set_tidal_center(vec c)	{tidal_center = c;}
	inline vec get_tidal_center()	const	{return tidal_center;}

	/// Set tidal field flag and alpha parameters.

	void set_tidal_field(real alpha1, real alpha3,
			     int type = 0, real omega = 0, vec c = 0.0) {
	    if (type > 0) set_tidal_field(type);
	    set_alpha(alpha1, alpha3);
	    set_tidal_center(c);
	}

	/// Plummer field (external field #2, bit 1).

	inline void set_plummer()		{SETBIT(external_field, 1);}
	inline bool get_plummer()	const	{return (GETBIT(external_field,
								1) > 0);}
	inline void unset_plummer()		{CLRBIT(external_field, 1);}

	inline real get_p_mass()	const	{return p_mass;}
	inline void set_p_mass(real m)		{p_mass = m;}
	inline real get_p_scale_sq()	const	{return p_scale_sq;}
	inline void set_p_scale_sq(real r2)	{p_scale_sq = r2;}
	inline void set_p_center(vec c) 	{p_center = c;}
	inline vec  get_p_center()	const	{return p_center;}
	inline bool get_p_friction()	const	{return p_friction;}
	inline void set_p_friction(bool f)	{p_friction = f;}

	/// Set Plummer field flag and parameters.

	inline void set_plummer(real m, real r2,
				vec c = 0.0, bool f = false) {  // all in one
	    set_plummer();
	    p_mass = m;
	    p_scale_sq = r2;
	    p_center = c;
	    p_friction = f;
	}

	/// Power-law field (external field #3, bit 2).

	inline void set_pl()			{SETBIT(external_field, 2);}
	inline bool get_pl()		const	{return (GETBIT(external_field,
								2) > 0);}
	inline void unset_pl()			{CLRBIT(external_field, 2);}

	inline real get_pl_coeff()	const	{return pl_coeff;}
	inline void set_pl_coeff(real c)	{pl_coeff = c;}
	inline real get_pl_scale()	const	{return pl_scale;}
	inline void set_pl_scale(real r)	{pl_scale = r;}
	inline real get_pl_exponent()	const	{return pl_exponent;}
	inline void set_pl_exponent(real e)	{pl_exponent = e;}
	inline void set_pl_center(vec c) 	{pl_center = c;}
	inline vec get_pl_center()	const	{return pl_center;}

	/// Set power-law field flag and parameters.

	inline void set_pl(real A, real a,
			   real e, vec c = 0.0) {	  // all in one

	    set_pl();

	    pl_coeff = A;
	    pl_scale = a;
	    pl_exponent = e;
	    pl_center = c;
	}

	/// Generic accessor for (single) external field (see dyn_external.C):

	real get_external_scale_sq(int bit = -1);

	/// Generic accessor for (single) external field (see dyn_external.C):

	vec get_external_center(int bit = -1);

	//-----------------------------------------------------------------

        void  set_pos(const vec& new_pos)      {pos = new_pos;}
        void  set_vel(const vec& new_vel)      {vel = new_vel;}
        void  set_acc(const vec& new_acc)      {acc = new_acc;}

	void  clear_pos()		       {pos = 0.0;}    ///< Set pos = 0
	void  clear_vel()                      {vel = 0.0;}    ///< Set vel = 0
	void  clear_acc()                      {acc = 0.0;}    ///< Set acc = 0

	inline void  inc_pos(const vec& d_pos) {pos += d_pos;} ///< pos += del
	inline void  inc_vel(const vec& d_vel) {vel += d_vel;} ///< vel += del
	inline void  inc_acc(const vec& d_acc) {acc += d_acc;} ///< acc += del

	/// Scale pos by the specified factor.

	inline void  scale_pos(const real scale_factor)  {pos *= scale_factor;}

	/// Scale vel by the specified factor.

	inline void  scale_vel(const real scale_factor)  {vel *= scale_factor;}

	/// Scale acc by the specified factor.

	inline void  scale_acc(const real scale_factor)  {acc *= scale_factor;}

	inline vec  get_pos()		const	{return pos;}
	inline vec  get_vel()		const	{return vel;}
	inline vec  get_acc()		const	{return acc;}

	inline dyn * get_parent() const
	    {return (dyn*) node::get_parent();}
	inline dyn * get_oldest_daughter() const
	    {return (dyn*)node::get_oldest_daughter();}
	inline dyn * get_younger_sister() const
	    {return (dyn*) node::get_younger_sister();}
	inline dyn * get_elder_sister() const
	    {return (dyn*) node::get_elder_sister();}

	/// Set or find the root node pointer.

        inline dyn * get_root() const
            {return (dyn*) node::get_root();}

	/// Return the top-level node of this node.

        inline dyn * get_top_level_node() const
            {return static_cast<dyn*>(node::get_top_level_node());}

	/// Find the binary sister of this node.

        inline dyn * get_binary_sister()
            {return (dyn*) node::get_binary_sister();}

	/// Compute the gravitational acceleration on this due to the system.

	void  calculate_acceleration(dyn *, real);

	inline kepler * get_kepler()	const	    {return kep;}
	void  set_kepler(kepler * new_kep)	    {kep = new_kep;}

	virtual void null_pointers();
	virtual void print_static(ostream &s = cerr);

        virtual istream& scan_dyn_story(istream&);
	virtual bool check_and_correct_node(bool verbose = true);

	virtual ostream& print_dyn_story(ostream &s,
					 bool print_xreal = true,
					 int short_output = 0);

        void to_com();		///< Transform to center-of-mass coordinates.
	void set_com(vec r = 0.0, vec v = 0.0);

	/// Place the root node at the origin, keeping the absolute pos
	/// and vel of all top-level nodes unchanged.

	void offset_com();

	/// Place the root node at the center of mass, keeping the
	/// absolute pos and vel of all top-level nodes unchanged.

	void reset_com();
        int flatten_node();	///< Remove all binary substructure.

	// Declaration of member function defined in dyn_tt.C.
	// There is no radius member datum in the dyn class, and 
	// get_radius will return 0 unless an R_eff entry is found
	// in the particle's dyn story.  For efficiency, the story
	// is only checked if explicitly requested.

        /// Optionally get R_eff from story; otherwise, (dyn version) return 0.

	virtual real get_radius(bool check_story = false);

        /// Optionally set R_eff in story; otherwise, (dyn version) do nothing.

	virtual void set_radius(real r, bool check_story = false);

        /// Scale R_eff in leaf dyn story, if present.

	virtual void scale_radius(real rfac);

	// Clump radius is just half the maximum diameter of a clump,
	// and has nothing directly to do with stellar radii in the
	// dyn case.

	real get_clump_radius();

	// Scale story R_eff if it exists; in _dyn_ scales radius too.

	// virtual void scale_radius(real fac);

        bool get_use_sstar()			{return use_sstar;}
	void set_use_sstar(bool u)		{use_sstar = u;}

	// Placeholders (~null in dyn, realized in hdyn)

	virtual bool nn_stats(real energy_cutoff, real kT,
			      vec center, bool verbose,
			      bool long_binary_output = true,
			      int which = 0);
	virtual real print_pert(bool long_binary_output = true,
				int indent = BIN_INDENT);
};

typedef dyn * dynptr;  // to enable dynamic array declarations such as
                       //    dyn** dyn_ptr_array = new dynptr[n];
                       // (note that the following expression is illegal:
                       //    dyn** dyn_ptr_array = new (dyn *)[n];)

typedef vec (dyn::*dyn_VMF_ptr)(void) const;	// vec member function pointer
typedef void (dyn::*dyn_MF_ptr)(const vec &);	// member function pointer

inline  node * new_dyn(hbpfp the_hbpfp,
		       sbpfp the_sbpfp,
		       bool use_stories)
    {return (node *) new dyn(the_hbpfp, the_sbpfp, use_stories);}

inline  dyn * mkdyn(int n, hbpfp the_hbpfp = new_hydrobase,
	                   sbpfp the_sbpfp = new_starbase)
    {return  (dyn *) mk_flat_tree(n, new_dyn, the_hbpfp, the_sbpfp);}

/// Read col format data from a stream.

dyn *get_col(istream& s = cin,
	     npfp the_npfp = new_dyn,
	     hbpfp the_hbpfp = new_hydrobase,
	     sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true);

/// Write col format data to a stream.

void put_col(dyn*, ostream& s = cout, bool put_time = true);

/// Read a standard format dyn snapshot from a stream.

dyn *get_dyn(istream & s = cin,
	     hbpfp the_hbpfp = new_hydrobase,
	     sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true);

/// Write a standard format dyn snapshot to a stream.

inline void put_dyn(dyn *b,			// note: not put_node, now!!
		    ostream &s = cout,
		    bool print_xreal = true,
		    int short_output = 0) {
    return dyn::get_col_output() ? put_col(b, s) :
	put_node(b, s, print_xreal, short_output);
}

/// Fast read, switching between file and pipe input.

dyn* fget_dyn(FILE* fp = dyn::get_ifp()?:stdin);

/// Return a pointer to the common ancestor of two nodes.

inline dyn * common_ancestor(dyn * bi, dyn * bj)
    {return (dyn*) common_ancestor((node*)bi, (node*)bj);}

/// Compute a vector quantity relative to the root node.

vec something_relative_to_root(dyn * b, dyn_VMF_ptr mf);

void dbg_message(char*, dyn*);
void dbg_message(char*, dyn*, dyn*);

// From dyn/init (used in sdyn3/evolve/bound3.C):

/// Create a uniform sphere of dyn particles.

void  makesphere(dyn * root, int n, real R = 1, int u_flag = 0);

// From dyn/util:

/// Compute the potential of a node in an N-body system.

real pot_on_general_node(dyn * bj, dyn * bi, real eps2, bool cm);

/// Calculate the potential, kinetic, and total energies of an N-body system.

void calculate_energies(dyn * root, real eps2,
			real & epot, real & ekin, real & etot,
			bool cm = false);

/// Comute and print the energies.

void print_recalculated_energies(dyn *, int mode = 0,
				 real eps2 = 0, real e_corr = 0);

/// Initialize the energy diagnostics and print the energies.

void initialize_print_energies(dyn *, real eps2 = 0);

/// Compute the densities associated with all nodes under b.

void compute_density(dyn* b,
		     int k = 12,
		     bool use_mass = false,
		     dyn** list = NULL,
		     int n_list = 0);

/// Combine two binary nodes into one.

void merge_low_level_nodes(dyn* b, real frac = 1, int option = 1);

// Special-case functions (superceded):

vector<real>& get_radial_densities(dyn *, vec, vector<real>&,
				   bool (*)(dyn*) = 0);

vector<real>& get_radial_numdensities(dyn *, vec, vector<real>&,
				      bool (*)(dyn*) = 0);

/// Compute the velocity dispersion in radial bins.

int get_radial_vdisp(dyn *b, vec cpos, vec cvel,
		     int n_zones, real r[], real v2[]);

/// Compute the density in radial bins.

int get_density_profile(dyn *b, vec cpos,
			int n_zones, real r[], real rho[],
			bool (*select)(dyn*) = NULL);
/// Compute the distribution of quantity q in radial bins.

int get_profile(dyn *b, vec cpos,
		int n_zones, real r[], real q[],
		real (*Q)(dyn*));

// Overloaded functions:

void compute_com(dyn*, vec&, vec&);	///< Calculate the center of mass.
void compute_com(dyn*);			///< Calculate the center of mass.

/// Calculate the modified center of mass.

void compute_mcom(dyn*, vec&, vec&, real f = 0.9, int n_iter = 2);

/// Calculate the modified center of mass.

void compute_mcom(dyn*, real f = 0.9, int n_iter = 2);

/// Compute the "standard" center (density, modified, normal center of mass.).

int  get_std_center(dyn*, vec&, vec&, bool verbose = false);

/// Compute the "standard" center (density, modified, normal center of mass.).

int  get_std_center(dyn*, bool verbose = false);

/// Compute the center defined by the particle with maximum density.

void compute_max_cod(dyn*, vec&, vec&);

/// Compute the center defined by the particle with maximum density.

void compute_max_cod(dyn*);

/// Compute the center, weighted by density.

void compute_mean_cod(dyn*, vec&, vec&);

/// Compute the center, weighted by density.

void compute_mean_cod(dyn*);

// From lagrad.C:

/// Compute Lagrangian radii.

void compute_mass_radii(dyn*);

/// Compute Lagrangian radii for 10-percentiles.

void compute_mass_radii_percentiles(dyn*);

/// Compute Lagrangian radii for quartiles.

void compute_mass_radii_quartiles(dyn*);

/// Set limiting masses for Lagrangian radii.

void set_lagr_cutoff_mass(dyn *b, real frac_low, real frac_high = 1.0);

/// Set limiting masses for Lagrangian radii.

void reset_lagr_cutoff_mass(dyn *b, real frac_low, real frac_high = 1.0);

/// Return lower limiting mass for Lagrangian radii.

real get_lagr_cutoff_mass();

/// Return lower limiting mass for Lagrangian radii.

real get_lagr_cutoff_mass_low();

/// Return upper limiting mass for Lagrangian radii.

real get_lagr_cutoff_mass_high();

/// Compute and print Lagrangian radii, write to root dyn story.

real print_lagrangian_radii(dyn* b,
			    int which_lagr = 2,
			    bool verbose = true,
			    int which_star = 0,
			    bool print = true);
/// Compute, don't print.

real compute_lagrangian_radii(dyn* b,
			      int which_lagr = 2,
			      bool verbose = true,
			      int which_star = 0);
/// Compute, don't print.

real compute_lagrangian_radii(dyn* b,
			      real *arr, int narr,
			      bool verbose = true,
			      int which_star = 0);
/// Compute, don't print.

real compute_lagrangian_radii(dyn* b,
			      real r,
			      bool verbose = true,
			      int which_star = 0);

/// Compute arbitrary Lagrangian radii.

typedef bool boolfn(dyn*);
bool compute_general_mass_radii(dyn*, int,
				bool nonlin = false,
				boolfn *bf = NULL,
				bool verbose = true);

// From sys_stats.C:

/// Search for bound pairs and print their properties.

void search_for_binaries(dyn* b, real energy_cutoff = 0, real kT = 0,
			 vec center = 0, bool verbose = true,
			 bool long_binary_output = true);

/// Parse the input for sys_stats and hsys_stats.

bool parse_sys_stats_main(int argc, char *argv[],
			  int  &which_lagr,
			  bool &binaries,
			  bool &short_output,
			  bool &B_flag,
			  bool &calc_e,
			  bool &n_sq,
			  bool &out,
			  int  &verbose,
			  char *cvs_id, char *source);
void check_addstar(dyn* b);

/// General statistics and diagnostics on an N-body system.

void sys_stats(dyn* root,
	       real energy_cutoff = 1,
	       int  verbose = 2,
	       bool binaries = true,
	       bool long_binary_output = false,
	       int which_lagr = 2,
	       bool print_time = false,
	       bool compute_energy = false,
	       bool allow_n_sq_ops = false,
	       void (*compute_energies)(dyn*, real, real&, real&, real&, bool)
			= calculate_energies,
	       void (*dstar_params)(dyn*) = NULL,
	       bool (*print_dstar_stats)(dyn*, bool, vec, bool) = NULL);

/// Determine the self-consistent mass of a cluster in a tidal field.

void refine_cluster_mass(dyn *b, int verbose = 0);

/// Determine the self-consistent mass of a cluster in a tidal field.

void refine_cluster_mass2(dyn *b, int verbose = 0);

// From dyn_stats.C:

/// Print the properties of a binary defined by a kepler.

real print_binary_params(kepler* k, real m1, real kT,
			 real dist_from_cod,
			 bool verbose = true,
			 bool long_output = true,
			 int init_indent = 0,
			 int indent = BIN_INDENT);

/// Return the total energy of the binary defined by two nodes.

real get_total_energy(dyn* bi, dyn* bj);

/// Return the semi-major axis of the binary defined by two nodes.

real get_semi_major_axis(dyn* bi, dyn* bj);

/// Return the period of the binary defined by two nodes.

real get_period(dyn* bi, dyn* bj);

/// Return the total energy and period of the binary defined by two nodes.

void get_total_energy_and_period(dyn* bi, dyn* bj, real& E, real& P);

/// Create a kepler structure from a pair of nodes.

void initialize_kepler_from_dyn_pair(kepler& k, dyn* bi, dyn* bj,
				     bool minimal = false);

/// Print the properties of a binary defined by a pair of nodes.

void print_binary_from_dyn_pair(dyn* bi, dyn* bj,
				real kT = 0,
				vec center = vec(0,0,0),
				bool verbose = true,
				bool long_output = true);

/// Recursively print the properties of binary node b.

real print_structure_recursive(dyn* b,
			       void (*dstar_params)(dyn*),
			       int& n_unp, real& e_unp,
			       real kT = 0.0,
			       vec center = vec(0,0,0),
			       bool verbose = true,
			       bool long_output = true,
			       int indent = 0);

/// Recursively print the properties of binary node b, simpler calling sequence.

real print_structure_recursive(dyn* b,
			       real kT = 0.0,
			       vec center = vec(0,0,0),
			       bool verbose = true,
			       bool long_output = true,
			       int indent = 2);

/// Compute core radius, mass, etc.

void compute_core_parameters(dyn* b, int k,
                             bool allow_n_sq_ops,
                             vec& center,
                             real& rcore, int& ncore, real& mcore);


// From plot_stars.C:

/// Pseudoplot a patch of an N-body system.

void plot_stars(dyn * bi,
		int n = 5,
		int k = 3);

// From scale.C:

/// Compute the total mass of a node.

real get_mass(dyn *b);

/// Scale the mass of a node (and descendents) to the specified value.

void scale_mass(dyn *b, real mscale);

/// Scale the positions of a node and descendents by the specified value.

void scale_pos(dyn *b, real rscale, vec com_pos = 0.0);

/// Scale the velocities of a node and descendents by the specified value.

void scale_vel(dyn *b, real vscale, vec com_vel = 0.0);

/// Compute the total kinetic energy of all top-level nodes under b.

real get_top_level_kinetic_energy(dyn *b);

/// Compute the total kinetic energy of all leaves under b.

real get_kinetic_energy(dyn *b);

/// Compute the total energies of all top-level nodes under b.

void get_top_level_energies(dyn *b, real eps2, real &potential, real &kinetic);

/// Scale velocities to produce the desired virial ratio.

void scale_virial(dyn *b, real q, real potential, real& kinetic,
		  vec com_vel = 0.0);

/// Scale the total energy of a system to the desired value.

real scale_energy(dyn *b, real e, real& energy,
		  vec com_pos = 0.0,
		  vec com_vel = 0.0);

/// Parse the command line in scale and hscale.

bool parse_scale_main(int argc, char *argv[],
		      real& eps, bool& c_flag,
		      bool& e_flag, real& e,
		      bool& m_flag, real& m,
		      bool& q_flag, real& q,
		      bool& r_flag, real& r,
		      bool& debug,
		      char *cvs_id, char *source);

/// Scale masses, radii, and velocities according to various criteria.

void scale(dyn *b, real eps,
	   bool c_flag,
	   bool e_flag, real e,
	   bool m_flag, real m,
	   bool q_flag, real q,
	   bool r_flag, real r,
	   bool debug = false,
	   void (*top_level_energies)(dyn*, real, real&, real&)
			= get_top_level_energies);

// From dyn_external.C:

/// Compute the acceleration and jerk due to an external field.

void get_external_acc(dyn * b,
		      vec pos,
		      vec vel,
		      real& pot,
		      vec& acc,
		      vec& jerk,
		      bool pot_only = false);

/// Compute the potential due to an external field.

real get_external_pot(dyn * b,
		      void (*pot_func)(dyn *, real) = NULL);

/// Compute the circular velocity in a (spherical) external potential.

real vcirc(dyn *b, vec r);

// Return the potential due to an external tidal field.

real get_tidal_pot(dyn *b);

// Return the potential due to an external Plummer field.

real get_plummer_pot(dyn *b);

// Return the potential due to an external power-law field.

real get_power_law_pot(dyn *b);

// Return the virial term due to an external field.

real get_external_virial(dyn * b);

// Print the properties of an external tidal field.

void print_external(unsigned int ext, bool shortbits = false);

/// Set half-mass radius used internally in dyn_external.

void set_new_rhalf(bool s = true);

void set_friction_beta(real b);	///< Set dynamical friction scaling parameter.
void set_friction_mass(real m);	///< Set mass determining dynamical friction.
void set_friction_vel(vec v);	///< Set vel determining dynamical friction.
void set_friction_acc(dyn *b, real r);	///< Dynamical friction acceleration.

// From add_plummer.C and add_power_law.C:

/// Determine physical scales associated with an N-body system.

bool get_physical_scales(dyn *b, real& mass, real& length, real& time);

/// Add parameters for an external Plummer field to an N-body system.

void add_plummer(dyn *b,
		 real coeff, real scale,
		 vec center = 0.0,
		 bool n_flag = false,
		 bool verbose = false,
		 bool fric_flag = false);

/// Turn dynamical friction on/off.

void toggle_plummer_friction(dyn *b);

/// Add parameters for an external power-law field to an N-body system.

void add_power_law(dyn *b,
		   real coeff, real exponent, real scale,
		   vec center = 0.0,
		   bool n_flag = false,
		   bool verbose = false,
		   bool G_flag = false);

// From dyn_story.C:

/// Read the initial mass from the root dyn story.

real get_initial_mass(dyn* b,
		      bool verbose = false,
		      bool mass_set = false,
		      real input_mass = 0);

/// Read the initial virial radius from the root dyn story.

real get_initial_virial_radius(dyn* b,
			       bool verbose = false,
			       bool r_virial_set = false,
			       real input_r_virial = 0);

/// Read the initial jacobi radius from the root dyn story.

real get_initial_jacobi_radius(dyn* b,
			       real r_virial,
			       bool verbose = false,
			       bool r_jacobi_set = false,
			       real input_r_jacobi = 0);

/// Read and set the tidal parameters from the root dyn story.

void set_tidal_params(dyn* b,
		      bool verbose,
		      real initial_r_jacobi,
		      real initial_mass,
		      int& tidal_field_type);

/// Check the tidal parameters.

void test_tidal_params(dyn* b,
		       bool verbose,
		       real initial_r_jacobi,
		       real initial_r_virial,
		       real initial_mass);

/// Check and set the tidal parameters.

int check_set_tidal(dyn *b, bool verbose = false);

// Check for and set Plummer parameters from input snapshot.

void check_set_plummer(dyn *b, bool verbose = false);

// Check for and set power-law parameters from input snapshot.

void check_set_power_law(dyn *b, bool verbose = false);

// Check for and set external field parameters from input snapshot.

void check_set_external(dyn *b, bool verbose = false, int fric_int = -1);

// Check for and set ignore_internal flag from input snapshot.

void check_set_ignore_internal(dyn *b, bool verbose = false);

// Optionally modify external parameters.

void update_external(dyn *b);

//----------------------------------------------------------------------
//
// Standard wrappers (for starcluster).
// See dyn/util/wrappers.C for functions not specified inline.

// "System" quantities (b is root node):

int  bound_number(dyn *b);	///< Number of bound stars.
real bound_mass(dyn *b);	///< Mass of bound stars.

/// Total angular momentum of the system.

vec  total_angular_momentum(dyn *b, vec x = 0, vec v = 0);
real total_energy(dyn *b);	///< Total energy of the system.
real core_radius(dyn *b);	///< Cluster core radius.
real core_mass(dyn *b);		///< Cluster core mass.
int  core_number(dyn *b);	///< Cluster core number.
real virial_radius(dyn *b);	///< Cluster virial radius.
real tidal_radius(dyn *b);	///< Cluster tidal radius.
real core_density(dyn *b);	///< Cluster core density.

// Quantities defined (mainly) for individual nodes:

inline real	mass(dyn *b)
		{
		    if (b->get_oldest_daughter()) {
			if (b->get_mass() > 0)
			    return b->get_mass();	// assume OK
			else
			    return get_mass(b);	// recompute
		    } else
			return b->get_mass();
		}

inline vec	pos(dyn *b, vec x = 0)	{return b->get_pos()-x;}
inline vec	vel(dyn *b, vec v = 0)	{return b->get_vel()-v;}

/// |pos|
inline real	distance(dyn *b, vec x = 0)
					{return abs(b->get_pos()-x);}
/// pos^2
inline real	distance_sq(dyn *b, vec x = 0)
					{return square(b->get_pos()-x);}
/// |vel|
inline real	speed(dyn *b, vec v = 0)
					{return abs(b->get_vel()-v);}
/// vel^2
inline real	speed_sq(dyn *b, vec v = 0)
					{return square(b->get_vel()-v);}

/// Squared radial velocity.

inline real	v_rad_sq(dyn *b, vec x = 0, vec v = 0)
					{
					    vec dx = b->get_pos()-x;
					    vec dv = b->get_vel()-v;
					    real dr2 = square(dx);
					    if (dr2 > 0)
						return square(dx*dv)/dr2;
					    else
						return 0;
					}
/// Radial velocity.

inline real	v_rad(dyn *b, vec x = 0, vec v = 0)
					{return sqrt(v_rad_sq(b, x, v));}

/// Squared transverse velocity.

inline real	v_trans_sq(dyn *b, vec x = 0, vec v = 0)
					{
					    vec dx = b->get_pos()-x;
					    vec dv = b->get_vel()-v;
					    real dr2 = square(dx);
					    if (dr2 > 0)
						return dv*dv
							 - square(dx*dv)/dr2;
					    else
						return 0;
					}

/// Transverse velocity.

inline real	v_trans(dyn *b, vec x = 0, vec v = 0)
					{return sqrt(v_trans_sq(b, x, v));}

/// Total energy.

inline real	energy(dyn *b, vec v = 0)
					{
					    if (find_qmatch(b->get_dyn_story(),
							    "pot")) {
						if (!b->get_parent())	// root
						    return total_energy(b);
						else {
						    real pot =
						      getrq(b->get_dyn_story(),
							    "pot");
						    vec dv = b->get_vel()-v;
						    return b->get_mass()
							     *(pot+0.5*dv*dv);
						}
					    } else
						return 0;
					}

/// Total angular momentum.

inline vec	angular_momentum(dyn *b, vec x = 0, vec v = 0)
					{
					    if (b->get_parent()) {
						vec dx = b->get_pos()-x;
						vec dv = b->get_vel()-v;
						return b->get_mass()*dx^dv;
					    } else
						return
						    total_angular_momentum(b);
					}

//----------------------------------------------------------------------

#include  "dyn_kepler.h"

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/dyn.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
