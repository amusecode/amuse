#include "src/stdinc.h"
#include "src/kepler.h"

static kepler *k;

int initialize_code()
{
    if (k) delete k;
    k = new kepler;
    if (!k) return -1;
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
    k->set_time(time);
    k->set_total_mass(mass);
    k->set_rel_pos(vec(x,y,z));
    k->set_rel_vel(vec(vx,vy,vz));
    k->initialize_from_pos_and_vel();
    return 0;
}

int initialize_from_elements(double mass, double semi, double ecc,
			     double mean_anomaly, double time,
			     double periastron)
{
    // Standard orbit will be in the x-y plane, with long axis along
    // x.  However, we can set the orientation separately, before or
    // after calling this function.  The code here is based on that
    // found in Starlab::make_standard_kepler().  Default mean anomaly
    // is 0 (i.e. periapsis) at the specified time (default 0).  Note
    // that ecc = 1 is a special case that requires periastron also to
    // be set.

    k->set_time(time);
    k->set_total_mass(mass);
    k->set_semi_major_axis(semi);
    k->set_eccentricity(ecc);
    k->set_mean_anomaly(mean_anomaly);
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

int get_total_mass(double * mass)
{
    *mass = k->get_total_mass();
    return 0;
}

int get_time(double * time)
{
    *time = k->get_time();
    return 0;
}

int get_period(double * period)
{
    *period = k->get_period();
    return 0;
}

int get_elements(double * semi, double * ecc)
{
    *semi = k->get_semi_major_axis();
    *ecc = k->get_eccentricity();
    return 0;
}

int get_integrals(double * energy, double * angular_momentum)
{
    *energy = k->get_energy();
    *angular_momentum = k->get_angular_momentum();
    return 0;
}

int get_separation_vector(double * x, double * y, double * z)
{
    vec r = k->get_rel_pos();
    *x = r[0];
    *y = r[1];
    *z = r[2];
    return 0;
}

int get_separation(double * r)
{
    *r = k->get_separation();
    return 0;
}

int set_periastron(double p)
{
    k->set_periastron(p);
    return 0;
}

int get_periastron(double * p)
{
    *p = k->get_periastron();
    return 0;
}

int get_apastron(double * a)
{
    *a = k->get_apastron();
    return 0;
}

int get_velocity_vector(double * vx, double * vy, double * vz)
{
    vec v = k->get_rel_vel();
    *vx = v[0];
    *vy = v[1];
    *vz = v[2];
    return 0;
}

int get_angles(double * mean_anomaly, double * true_anomaly)
{
    *mean_anomaly = k->get_mean_anomaly();
    *true_anomaly = k->get_true_anomaly();
    return 0;
}

// Remains TBD if this is actually the best way to handle orbit orientation.

#define TOL 1.e-10

int set_longitudinal_unit_vector(double x, double y, double z)
{
    vec l = k->get_longitudinal_unit_vector();
    vec t = k->get_transverse_unit_vector();
    vec n = k->get_normal_unit_vector();
    real r = sqrt(x*x+y*y+z*z);
    if (r <= 0) return -1;
    l = vec(x,y,z)/r;
    if (fabs(l*n) < TOL*abs(l)*abs(n)) t = n^l;
    k->set_orientation(l, t, n);
    return 0;
}

int set_transverse_unit_vector(double x, double y, double z)
{
    vec l = k->get_longitudinal_unit_vector();
    vec t = k->get_transverse_unit_vector();
    vec n = k->get_normal_unit_vector();
    real r = sqrt(x*x+y*y+z*z);
    if (r <= 0) return -1;
    t = vec(x,y,z)/r;
    if (fabs(t*n) < TOL*abs(t)*abs(n)) t = n^t;
    k->set_orientation(l, t, n);
    return 0;
}

int set_normal_unit_vector(double x, double y, double z)
{
    vec l = k->get_longitudinal_unit_vector();
    vec t = k->get_transverse_unit_vector();
    vec n = k->get_normal_unit_vector();
    real r = sqrt(x*x+y*y+z*z);
    if (r <= 0) return -1;
    n = vec(x,y,z)/r;
    if (fabs(l*n) < TOL*abs(l)*abs(n)) t = n^l;
    k->set_orientation(l, t, n);
    return 0;
}

int get_longitudinal_unit_vector(double * x, double * y, double * z)
{
    vec l = k->get_longitudinal_unit_vector();
    *x = l[0];
    *y = l[1];
    *z = l[2];
    return 0;
}

int get_transverse_unit_vector(double * x, double * y, double * z)
{
    vec t = k->get_transverse_unit_vector();
    *x = t[0];
    *y = t[1];
    *z = t[2];
    return 0;
}


int get_normal_unit_vector(double * x, double * y, double * z)
{
    vec n = k->get_normal_unit_vector();
    *x = n[0];
    *y = n[1];
    *z = n[2];
    return 0;
}
