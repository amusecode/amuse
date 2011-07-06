#include "src/stdinc.h"
#include "src/kepler.h"

static kepler *k;

int initialize_code()
{
    k = NULL;
    return 0;
}

int commit_parameters()
{
    return 0;
}

int recommit_parameters()
{
    return 0;
}

int cleanup_code()
{
    return 0;
}

int initialize_from_dyn(double mass,
			double x, double y, double z,
			double vx, double vy, double vz,
			double time)
{
    if (k) delete k;
    k = new kepler;
    if (!k) return -1;
    k->set_time(time);
    k->set_total_mass(mass);
    k->set_rel_pos(vec(x,y,z));
    k->set_rel_vel(vec(vx,vy,vz));
    k->initialize_from_pos_and_vel();
    return 0;
}

int initialize_from_elements(double mass, double semi, double ecc,
			     double time)
{
    if (k) delete k;
    k = new kepler;
    if (!k) return -1;

    // Standard orbit will be in the x-y plane, with long axis along
    // x.  Later, we will add the option of randomizing the
    // orientation.  Code based on make_standard_kepler().  Mean
    // anomaly is 0 (i.e. periapsis) at the specified time.  Note that
    // ecc = 1 is a special case that isn't addressed here (more info
    // is needed).

    k->set_time(time);
    k->set_total_mass(mass);
    k->set_semi_major_axis(semi);
    k->set_eccentricity(ecc);
    k->set_mean_anomaly(0);
    k->align_with_axes(1);
    k->initialize_from_shape_and_phase();	// expects a, e [, q [, E]]
    return 0;
}

int transform_to_time(double time)
{
    k->transform_to_time(time);
    return 0;
}

// All advance/return functions should check the allowability of the
// command and return the appropriate status.  The kepler functions
// will produce warning messages and do something, but won't return an
// error indicator.

int advance_to_radius(double radius)
{
    k->advance_to_radius(radius);
    return 0;
}

int return_to_radius(double radius)
{
    k->return_to_radius(radius);
    return 0;
}

int advance_to_periastron()
{
    k->advance_to_periastron();
    return 0;
}

int advance_to_apastron()
{
    k->advance_to_apastron();
    return 0;
}

int return_to_periastron()
{
    k->return_to_periastron();
    return 0;
}

int return_to_apastron()
{
    k->return_to_apastron();
    return 0;
}

int get_elements(double * semi, double * ecc)
{
    *semi = k->get_semi_major_axis();
    *ecc = k->get_eccentricity();
    return 0;
}

int get_separation(double * x, double * y, double * z)
{
    vec r = k->get_rel_pos();
    *x = r[0];
    *y = r[1];
    *z = r[2];
    return 0;
}
