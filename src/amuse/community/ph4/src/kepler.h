
// kepler.h: Class for managing the internal dynamics of a 2-body system.
//
//  version 1:  Nov 1993   Piet Hut & Steve McMillan, for Starlab
//  version 2:	Feb 2011   Steve McMillan, imported to AMUSE

#ifndef  KEPLER_H
#  define  KEPLER_H

#include  "vec.h"
class hdyn;

class  kepler
{
    protected:

	// Constants of the orbit:

	real  total_mass;                    ///< sum of the two masses

	real  energy;		             ///< per unit reduced mass
	real  angular_momentum;	             ///< per unit reduced mass

	real  mean_motion;		     // keep these separate to
	real  period;			     // handle special cases

	vec  normal_unit_vector;	     ///< norm = (long ^ trans)
	vec  longitudinal_unit_vector;       ///< along the major axis
					     //   pointing from the focus
					     //   to the pericenter
	vec  transverse_unit_vector;	     ///< along the minor axis
					     //   pointing from the focus to
					     //   a true anomaly of 90 degrees

	real  semi_major_axis;
	real  eccentricity;
	real  periastron;

	/// Time of the most recent periastron in the elliptical case.

	real  time_of_periastron_passage;    

	// Dynamical variables:

        real  time;                          ///< current time

	vec  rel_pos;                        ///< r_2 - r_1
	vec  rel_vel;                        ///< v_2 - v_1

	real  separation;                    ///< separation between components
	                                     //   NOTE: ALWAYS update separation
					     //   when updating rel_pos.

	real  true_anomaly;                  // helpful when starting with a
	                                     //   prescribed separation.
	real  mean_anomaly;                  // helpful when starting with a
                                             //   prescribed time.
	// Predicted quantities:

        real pred_time;                      ///< predicted value
	vec  pred_rel_pos;                   ///< predicted value
	vec  pred_rel_vel;                   ///< predicted value
	real pred_separation;                ///< predicted value
	real pred_true_anomaly;              ///< predicted value
	real pred_mean_anomaly;              ///< predicted value


    	real circular_binary_limit;	     // Use to allow parts of code to
					     // make nearly circular binaries
					     // exactly circular, but don't
					     // force it on other functions.

	// Local functions:

	void fast_to_pred_rel_pos_vel(real r_cos_true_an,
				      real r_sin_true_an);
	void fast_pred_mean_anomaly_to_pos_and_vel();

        void to_pred_rel_pos_vel(real, real);
        void to_pred_rel_pos_vel_linear(real);
	void pred_true_to_mean_anomaly();
	void mean_anomaly_to_periastron_passage();
        void update_time_of_periastron();
	void set_real_from_pred();
        void pred_mean_anomaly_to_pos_and_vel();

    public:

	// Variables to manage perturbations due to partial
	// unperturbed motion in an external (tidal) field.

	real tidal_potential;		     // tidal potential at the start
					     // of unperturbed motion

	real delta_angular_momentum;	     // tidal change in angular
					     // momentum change during
					     // unperturbed motion; only the
					     // component perpendicular to
					     // the orbit plane is considered

        kepler();
	virtual ~kepler() {}

        void set_time(const real new_time)    	  {time = new_time;}

	/// Offset all times by the specified amount.

        void offset_time(const real dt) {
	    time += dt;
	    time_of_periastron_passage += dt;
	}
        void set_total_mass(const real new_mass)  {total_mass = new_mass;}
        void set_rel_pos(const vec& new_pos)	  {rel_pos = new_pos;
						   separation = abs(rel_pos);}
        void set_rel_vel(const vec& new_vel)	  {rel_vel = new_vel;}

	void set_energy(const real E)		  {energy = E;}
	void set_semi_major_axis(const real a)	  {semi_major_axis = a;}
	void set_angular_momentum(const real h)   {angular_momentum = h;}
	void set_eccentricity(const real e)	  {eccentricity = e;}
	void set_periastron(const real q)	  {periastron = q;}

	void set_mean_anomaly(const real M)	  {mean_anomaly = M;}

	/// Set all vectors defining the orbit.

        void set_orientation(vec& l, vec& t, vec& n) {
	    longitudinal_unit_vector = l;
	    transverse_unit_vector = t;
	    normal_unit_vector = n;
	}

	/// Set orientation vectors along x, y, z axes.

	void  align_with_axes(int);

	real get_circular_binary_limit()	{return circular_binary_limit;}
	void set_circular_binary_limit(real e)	{circular_binary_limit = e;}

        void initialize_from_pos_and_vel(bool minimal = false,
					 bool verbose = true);
        void initialize_from_shape_and_phase();
	void initialize_from_integrals_and_separation(bool receding = true);

        real get_time()    	  		{return time;}
        real get_pred_time()    		{return pred_time;}
        real get_total_mass()			{return total_mass;}
        vec  get_rel_pos()  			{return rel_pos;}
        vec  get_rel_vel()  			{return rel_vel;}

        real get_energy()			{return energy;}
	real get_semi_major_axis()		{return semi_major_axis;}
        real get_angular_momentum()		{return angular_momentum;}
	real get_eccentricity()		  	{return eccentricity;}
	real get_periastron()			{return periastron;}
	real get_apastron()			{return semi_major_axis
						       * (1 + eccentricity);}
	real get_period()			{return period;}

	real get_mean_motion()			{return mean_motion;}

	vec  get_normal_unit_vector()		{return normal_unit_vector;}
	vec  get_longitudinal_unit_vector()
		{return longitudinal_unit_vector;}
	vec  get_transverse_unit_vector()
		{return transverse_unit_vector;}

	real get_separation()			{return separation;}
	real get_true_anomaly()			{return true_anomaly;}
	real get_mean_anomaly()			{return mean_anomaly;}
	real get_time_of_periastron_passage()
		{return time_of_periastron_passage;}

	// The following functions transform the relative position and
	// velocity, as well as the other internal variables (`radius'
	// is the separation between the two stars).

        real advance_to_radius(real);		///< Transform forward in time.
        real return_to_radius(real);		///< Transform backward in time.

        real advance_to_periastron();
        real return_to_periastron();
        real advance_to_apastron();
        real return_to_apastron();

        real transform_to_time(real);

	// The following functions transform only the predicted relative
	// position and velocity, and predicted mean and true anomaly, but
	// leave all other internal variables unchanged.

        real pred_transform_to_radius(real, int);
        real pred_advance_to_radius(real);	///< Transform forward in time.
        real pred_return_to_radius(real);	///< Transform backward in time.

        real pred_advance_to_periastron();
        real pred_return_to_periastron();
        real pred_advance_to_apastron();
        real pred_return_to_apastron();

        real pred_transform_to_time(real);

	void print_dyn(ostream & s = cout);
	void print_elements(ostream & s = cout);
	void print_all(ostream & s = cout);
};

// Functions in kepler.cc:

void set_kepler_tolerance(int);  //	0: quit on any trig error
                                 //	1: quit on large trig error
                                 //	2: force to periastron or apastron on
                                 //        trig error

void set_kepler_print_trig_warning(bool);

void make_standard_kepler(kepler &k, real t, real mass,
			  real energy, real eccentricity,
			  real q, real mean_anomaly,
			  int align_axis = 0);

kepler* hdyn_to_kepler(hdyn * com, real t = 0);
void new_kepler(hdyn * com, real t = 0);

void set_random_orientation(kepler &k, int planar = 0);
void print_orbital_elements(hdyn *bi, hdyn *bj, bool verbose = true);

//  The following functions allow the user to access the "fast" solver
//  for the elliptical Kepler's equation:

void set_kepler_fast_flag(bool flag = true);
void clear_kepler_fast_flag();
bool get_kepler_fast_flag();

void set_kepler_fast_tol(real tol = 1.e-6);
void reset_kepler_fast_tol();
real get_kepler_fast_tol();

void set_kepler_print_trig_warning(bool p);

#endif
