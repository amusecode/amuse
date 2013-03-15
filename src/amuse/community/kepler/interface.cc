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
			     double periastron, int random_orientation)
{
    // Standard orbit will be in the x-y plane, with long axis along
    // x.  However, we can set the orientation separately, before or
    // after calling this function.  The code here is based on that
    // found in Starlab::make_standard_kepler().  The mean anomaly on
    // entry is usually taken to be 0 (i.e. periapsis) at the
    // specified time (default 0), unless the orbit is linear.
    //
    // Note that ecc = 1 is a special case that requires periastron
    // also to be set.  In that case, the sign of semi determines
    // whether the orbit is bound or unbound.

    k->set_time(time);
    k->set_total_mass(mass);
    k->set_semi_major_axis(semi);
    k->set_eccentricity(ecc);
    k->set_periastron(periastron);
    k->set_mean_anomaly(mean_anomaly);

    if (random_orientation)
	set_random_orientation(*k, 0);

    k->initialize_from_shape_and_phase();

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
    //if (fabs(l*n) < TOL*abs(l)*abs(n)) t = n^l;	// no side effect
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
    //if (fabs(t*n) < TOL*abs(t)*abs(n)) t = n^t;	// no side effect
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
    //if (fabs(l*n) < TOL*abs(l)*abs(n)) t = n^l;	// no side effect
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

int print_all()
{
    k->print_all();
    return 0;
}

//--------------------------------------------------------------------
// From Amanda Benson:  Extra code to initialize a 3-body encounter,
// which is too slow in python relative to the cost of a typical
// 3-body scattering.

static void set_inner_orbit(real ecc)
{
    real mass = 1;	// standard assumption from
    real semi = 1;	// Hut & Bahcall 1983
    k->align_with_axes(1);
    real mean_an = 2*M_PI*randinter(0, 1);
    initialize_from_elements(mass, semi, ecc, mean_an, 0.0, 0.0, 0);
    cout << "inner "; PRC(semi); PRL(ecc);
}

static void set_outer_orbit(real m, real M,
			    real v_inf, real impact_parameter,
			    int planar)
{
    real mtotal = 1 + M;

    // Multiply the incoming velocity by the critical value (see Hut &
    // Bahcall 1983).

    v_inf *= sqrt((1-m)*m*mtotal/M);
    real energy3 = 0.5*v_inf*v_inf;
    cout << "m1, m2, m3 = " << 1-m << " " << m << " " << M << endl << flush;
    PRC(v_inf); PRC(energy3); PRL(impact_parameter);

    real semi, ang_mom3, ecc, periastron;
    if (energy3 > 0) {
	semi = -0.5*mtotal/energy3;
	ang_mom3 = impact_parameter*v_inf;
	ecc = sqrt(1+2*energy3*pow(ang_mom3/mtotal, 2));
	periastron = -semi*fmax(ecc-1, 0.0);	// not used
    } else {
	semi = 0;				// not used
	ecc = 1;
	periastron = impact_parameter;
    }

    // Orientation:

    vec longitudinal(1,0,0), transverse(0,1,0), normal(0,0,1);
    if (planar < 0) {
	transverse = -transverse;
	normal = -normal;
    } else if (planar == 0) {
	real costheta = randinter(-1, 1);
	real sintheta = sqrt(fmax(0., 1-pow(costheta, 2)));
	real phi = 2*M_PI*randinter(-1, 1);
	longitudinal = vec(sintheta*cos(phi), sintheta*sin(phi), costheta);
	vec temp(1,0,0);
	if (abs(longitudinal[0]) > 0.5) temp = vec(0,1,0);
	transverse = longitudinal^temp;
	transverse /= abs(transverse);
	normal = longitudinal^transverse;

	real psi = randinter(0,2*M_PI);
	real cospsi = cos(psi);
	real sinpsi = sin(psi);
	normal = cospsi*normal + sinpsi*transverse;
	transverse = normal^longitudinal;
    }
    k->set_orientation(longitudinal, transverse, normal);

    real time = 0;
    real mean_anomaly = 0;
    if (periastron == 0) mean_anomaly = -1.e-3;
    PRL(mean_anomaly);

    cout << "outer "; PRC(semi); PRL(ecc);
    initialize_from_elements(mtotal, semi, ecc, mean_anomaly,
			     time, periastron, 0);
    PRL(k->get_normal_unit_vector());
    PRL(k->get_periastron());
    // k->print_all();
}

int set_random(int seed)
{
    cout << "set_random: seed = " << seed << endl;
    srandinter(seed, 0);	// kepler built-in
    return 0;
}

// Binary scattering in standard (Hut & Bahcall 1983) units.  Binary
// mass = 1, semimajor axis = 1, t = 0 at pericenter if appropriate.

int make_binary_scattering(real m, real ecc,
			   real M, real v_inf, real impact_parameter,
			   real gamma, int planar,
			   double * time,
			   double * m1, double * m2, double * m3,
			   double * x1, double * x2, double * x3,
			   double * y1, double * y2, double * y3,
			   double * z1, double * z2, double * z3,
			   double * vx1, double * vx2, double * vx3,
			   double * vy1, double * vy2, double * vy3,
			   double * vz1, double * vz2, double * vz3)
{
    // Inner orbit (1,2).

    *m1 = 1 - m;
    *m2 = m;

    set_inner_orbit(ecc);
    vec rel_pos = k->get_rel_pos();
    vec rel_vel = k->get_rel_vel();
    vec pos1 = -m*rel_pos;
    vec pos2 = (1-m)*rel_pos;
    vec vel1 = -m*rel_vel;
    vec vel2 = (1-m)*rel_vel;

    // Outer orbit ((1,2),3).

    *m3 = M;

    set_outer_orbit(m, M, v_inf, impact_parameter, planar);
    // cout << "----------" << endl << flush;
    // k->print_all();
    k->return_to_radius(pow(gamma/M, -1./3));
    // cout << "----------" << endl << flush;
    // k->print_all();
    // cout << "----------" << endl << flush;
    *time = k->get_time();

    rel_pos = k->get_rel_pos();
    rel_vel = k->get_rel_vel();
    real f = M/(1+M);
    vec pos12 = -f*rel_pos, pos3 = (1-f)*rel_pos;
    vec vel12 = -f*rel_vel, vel3 = (1-f)*rel_vel;

    PRC(k->get_separation()); PRL(k->get_time());

    pos1 += pos12;
    pos2 += pos12;
    vel1 += vel12;
    vel2 += vel12;

    PRL(*time);

    // Remaining return arguments:

    *x1 = pos1[0];
    *x2 = pos2[0];
    *x3 = pos3[0];
    *y1 = pos1[1];
    *y2 = pos2[1];
    *y3 = pos3[1];
    *z1 = pos1[2];
    *z2 = pos2[2];
    *z3 = pos3[2];
    *vx1 = vel1[0];
    *vx2 = vel2[0];
    *vx3 = vel3[0];
    *vy1 = vel1[1];
    *vy2 = vel2[1];
    *vy3 = vel3[1];
    *vz1 = vel1[2];
    *vz2 = vel2[2];
    *vz3 = vel3[2];

    return 0;
}
