
//=============================================================================
//
//       N-body integration module with shared but variable time step
//       (the same for all particles but its size changing in time),
//       using the Hermite integration scheme.
//                        
//       ref.: Hut, P., Makino, J. & McMillan, S., 1995,
//             Astrophysical Journal Letters 443, L93-L96.
//                
//  All module functions are included in this file.
//_____________________________________________________________________________
//
//    version 1:  Jan 2002   Piet Hut, Jun Makino
//    version 2:  Mar 2006   Jun Makino, Steve McMillan
//    version 3:  Jan 2007   Steve McMillan (MODEST-7b)
//    version 3a: Aug 2007   Steve McMillan (MODEST-7a)
//=============================================================================

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>

#include "muse_dynamics.h"

#include <vector>
#include <algorithm>
#include <string>
#include <map>

using namespace std;
typedef double  real;

#include "vec3.h"		// borrowed from starlab

#define PR(x)  *sout << #x << " = " << x << " " << flush
#define PRC(x) *sout << #x << " = " << x << ",  " << flush
#define PRL(x) *sout << #x << " = " << x << endl << flush

// (Global!) static data:

const int NDIM = 3;		// number of spatial dimensions

// N-body data:

real  t = 0;
real t_evolve = t;		// Time requested by evolve.  Control returns
				// when t <= t_evolve and t + dt > t_evolve,
				// and this is assumed when the state of the
				// system is computed by extrapolation.

vector<int>  ident;
vector<real> mass, radius;
vector<vec>  pos, vel, acc, jerk;

// Control parameters:

const real DT_PARAM = 0.03;
const real DT_DIA = 1;

real dt_param = DT_PARAM;	// control parameter to determine time step size
real dt_dia = DT_DIA;		// time interval between diagnostic output

bool x_flag = false;		// set true for serious debugging only!

int nsteps = 0;			// number of integration time steps completed
real einit = 0;			// initial total energy of the system
bool init_flag = false;
real t_dia = 0;
real eps2 = 0;

bool flag_collision = true;
bool reeval = false;
bool test_mode = false;
ostream* sout = &cout;

// Accessors for use by the C++ main driver only:

void set_t(real tt)		{t = tt;}
void set_dt_param(real dt)	{dt_param = dt;}
void set_dt_dia(real dt)	{dt_dia = dt;}
void set_eps(real eps) 		{eps2 = eps*eps;}



// (Most of) the original code (some arguments replaced by global data):

//-----------------------------------------------------------------------------
//  write_diagnostics  --  writes diagnostics on the output stream cout:
//                         current time; number of integration steps so far;
//                         kinetic, potential, and total energy; absolute and
//                         relative energy errors since the start of the run.
//                         If x_flag (x for eXtra data) is true, all internal
//                         data are dumped for each particle (mass, position,
//                         velocity, acceleration, and jerk).
//
//  Note: the kinetic energy is calculated here, while the potential
//  energy is calculated in the function get_acc_jerk_pot_coll().
//-----------------------------------------------------------------------------

void write_diagnostics(real epot, ostream& s = cout)
{
    int n = ident.size();
    real ekin = 0;			// kinetic energy
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            ekin += 0.5 * mass[i] * vel[i][k] * vel[i][k];

    real etot = ekin + epot;		// total energy

    if (!init_flag) {			// on first pass, set
        einit = etot;			// the initial energy
	init_flag = true;
	return;				// suppress initial output
    }

    s << "    internal diagnostics at time t = " << t
	 << " after " << nsteps << " steps"
	 << endl
         << "        E_kin = " << ekin
         << "  E_pot = " << epot
         << "  E_tot = " << etot << endl;
    s << "        "
         << "absolute energy error  E_tot - E_init = "
         << etot - einit << endl;
    s << "        "
         << "relative energy error  (E_tot - E_init) / E_init = "
         << (etot - einit) / einit << endl;

    if (x_flag) {
	s << "        system dump, n = " << n << endl;
        for (int i = 0; i < n ; i++){
            s << "        data for particle " << ident[i]
		 << ": " << endl;
            s << "            "; s << mass[i] << endl;
            s << "            "; s << radius[i] << endl;
	    s << "           "; 
            for (int k = 0; k < NDIM; k++)
                s << ' ' << pos[i][k];
            s << endl;
	    s << "           "; 
            for (int k = 0; k < NDIM; k++)
                s << ' ' << vel[i][k];
            s << endl;
	    s << "           "; 
            for (int k = 0; k < NDIM; k++)
                s << ' ' << acc[i][k];
            s << endl;
	    s << "           "; 
            for (int k = 0; k < NDIM; k++)
                s << ' ' << jerk[i][k];
            s << endl;
        }
    }
}
    
//-----------------------------------------------------------------------------
//  predict_step  --  take the first approximation of one Hermite integration
//                    step, advancing the positions and velocities through a
//                    Taylor series development up to the order of the jerks.
//		      Note that all pos and vel are overwritten.
//-----------------------------------------------------------------------------

void predict_step(real dt)
{
    if (dt <= 0) return;

    int n = ident.size();
    for (int i = 0; i < n ; i++) {
        for (int k = 0; k < NDIM ; k++){
            pos[i][k] += vel[i][k]*dt + acc[i][k]*dt*dt/2
                                      + jerk[i][k]*dt*dt*dt/6;
            vel[i][k] += acc[i][k]*dt + jerk[i][k]*dt*dt/2;
        }
    }
}

//-----------------------------------------------------------------------------
//  get_acc_jerk_pot_coll  --  calculate accelerations and jerks, and as
//                             side effect also calculates potential energy
//                             and the time scale coll_time for significant
//                             changes in local configurations to occur.
//                                                  __                     __
//                                                 |          -->  -->       |
//               M                           M     |           r  . v        |
//   -->          j    -->       -->          j    | -->        ji   ji -->  |
//    a   ==  --------  r    ;    j   ==  -------- |  v   - 3 ---------  r   |
//     ji     |-->  |3   ji        ji     |-->  |3 |   ji      |-->  |2   ji |
//            | r   |                     | r   |  |           | r   |       |
//            |  ji |                     |  ji |  |__         |  ji |     __|

//  Note: it would be cleaner to calculate potential energy and
//	  collision time in a separate function.  However, the current
//	  function is by far the most time consuming part of the whole
//	  program, with a double loop over all particles that is executed
//	  every time step.  Splitting off some of the work to another
//	  function would significantly increase the total computer time
//	  (by an amount close to a factor two).
//
//  With softening, the r^3 in the denominator becomes (r^2 + e^2)^{1.5}.
//
//  We determine the values of all four quantities of interest by
//  walking through the system in a double {i,j} loop.  The first
//  three, acceleration, jerk, and potential energy, are calculated by
//  adding successive terms; the last, the estimate for the collision
//  time, is found by determining the minimum value over all particle
//  pairs and over the two choices of collision time,
//  position/velocity and sqrt(position/acceleration), where position
//  and velocity indicate their relative values between the two
//  particles, while acceleration indicates their pairwise
//  acceleration.  At the start, the first three quantities are set to
//  zero, to prepare for accumulation, while the last one is set to a
//  very large number, to prepare for minimization.
//
//  The integration loops only over half of the pairs, with j > i,
//  since the contributions to the acceleration and jerk of particle j
//  on particle i are the same as those of particle i on particle j,
//  apart from a minus sign and a different mass factor.
//-----------------------------------------------------------------------------

static int id_coll_primary = -1, id_coll_secondary = -1;

void get_acc_jerk_pot_coll(real & epot, real & coll_time)
{
    int n = ident.size();
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            acc[i][k] = jerk[i][k] = 0;

    epot = 0;
    const real VERY_LARGE_NUMBER = 1e300;
    real coll_time_q = VERY_LARGE_NUMBER;      // collision time to 4th power
    real coll_est_q;                           // collision time scale estimate
                                               // to 4th power (quartic)
    id_coll_primary = id_coll_secondary = -1;

    for (int i = 0; i < n ; i++){
        for (int j = i+1; j < n ; j++){            // rji[] is the vector from
            real rji[NDIM];                        // particle i to particle j
            real vji[NDIM];                        // vji[] = d rji[] / d t
            for (int k = 0; k < NDIM ; k++){
                rji[k] = pos[j][k] - pos[i][k];
                vji[k] = vel[j][k] - vel[i][k];
            }
            real r2 = 0;                           // | rji |^2
            real v2 = 0;                           // | vji |^2
            real rv_r2 = 0;                        // ( rij . vij ) / | rji |^2
            for (int k = 0; k < NDIM ; k++){
                r2 += rji[k] * rji[k];
                v2 += vji[k] * vji[k];
                rv_r2 += rji[k] * vji[k];
            }
            rv_r2 /= r2;

	    if (id_coll_primary < 0) {
		real rsum = radius[i] + radius[j];
		if (r2 <= rsum*rsum) {
		    if (mass[i] >= mass[j]) {
			id_coll_primary = ident[i];
			id_coll_secondary = ident[j];
		    } else {
			id_coll_primary = ident[j];
			id_coll_secondary = ident[i];
		    }
		}
	    }

	    r2 += eps2;				   // | rji |^2 + eps^2
            real r = sqrt(r2);                     // | rji |
            real r3 = r * r2;                      // | rji |^3

	    // Add the {i,j} contribution to the total potential energy.

            epot -= mass[i] * mass[j] / r;

	    // Add the {j (i)} contribution to the {i (j)} acc and jerk.

            real da[NDIM];                         // main terms in pairwise
            real dj[NDIM];                         // acceleration and jerk
            for (int k = 0; k < NDIM ; k++){
                da[k] = rji[k] / r3;                          // see equations
                dj[k] = (vji[k] - 3 * rv_r2 * rji[k]) / r3;   // in the header
            }
            for (int k = 0; k < NDIM ; k++){
                acc[i][k]  += mass[j] * da[k];                // using symmetry
                acc[j][k]  -= mass[i] * da[k];                // find pairwise
                jerk[i][k] += mass[j] * dj[k];                // acceleration
                jerk[j][k] -= mass[i] * dj[k];                // and jerk
            }

	    // First collision time estimate is based on unaccelerated
	    // linear motion.

            coll_est_q = (r2*r2) / (v2*v2);
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q;

	    // Second collision time estimate is based on free fall.

            real da2 = 0;                                  // da2 becomes the 
            for (int k = 0; k < NDIM ; k++)                // square of the 
                da2 += da[k] * da[k];                      // pairwise accel-
            real mij = mass[i] + mass[j];                  // eration between
            da2 *= mij * mij;                              // particles i and j

            coll_est_q = r2/da2;
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q;
        }                                     
    }                                               // from q for quartic back
    coll_time = sqrt(sqrt(coll_time_q));            // to linear collision time
}                                             

//-----------------------------------------------------------------------------
//  correct_step  --  take one iteration to improve the new values of
//                    position and velocities, constructing a higher-order
//                    Taylor series from the terms up to jerk at the
//                    beginning and the end of the time step.  This symmetric
//		      formulation is not the original one from Makino and
//		      Aarseth; it comes from ACS (Hut & Makino).
//-----------------------------------------------------------------------------

void correct_step(const real old_pos[][NDIM], const real old_vel[][NDIM], 
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  real dt)
{
    int n = ident.size();
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
            vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2
                                      + (old_jerk[i][k] - jerk[i][k])*dt*dt/12;
            pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2
                                      + (old_acc[i][k] - acc[i][k])*dt*dt/12;
        }
}

//-----------------------------------------------------------------------------
//  evolve_step  --  take one integration step for an N-body system, using
//                   the Hermite algorithm.
//-----------------------------------------------------------------------------

void evolve_step(real dt, real & epot, real & coll_time)
{
    int n = ident.size();
    real (* old_pos)[NDIM] = new real[n][NDIM];
    real (* old_vel)[NDIM] = new real[n][NDIM];
    real (* old_acc)[NDIM] = new real[n][NDIM];
    real (* old_jerk)[NDIM] = new real[n][NDIM];

    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
	    old_pos[i][k] = pos[i][k];
	    old_vel[i][k] = vel[i][k];
	    old_acc[i][k] = acc[i][k];
	    old_jerk[i][k] = jerk[i][k];
        }

    predict_step(dt);
    get_acc_jerk_pot_coll(epot, coll_time);
    correct_step(old_pos, old_vel, old_acc, old_jerk, dt);
    if (reeval) get_acc_jerk_pot_coll(epot, coll_time);
    t += dt;
    t_evolve = t;

    delete[] old_pos;
    delete[] old_vel;
    delete[] old_acc;
    delete[] old_jerk;
}

real calculate_step(real coll_time)
{
    // Determine the new system time step from coll_time.

    real step = dt_param;
    if (!test_mode) {

	step *= coll_time;
#if 0
	// Round down to the next power of 2.

	real step2 = 1;
	while (step2 > step) step2 /= 2;
	step = step2;
#endif
    }

    return step;
}

void compute_nn()
{
    // For debugging only.

    real rijmin2 = 1.e30;
    int imin = -1, jmin = -1;
    int n = ident.size();

    for (int i = 0; i < n-1; i++) {
	for (int j = i+1; j < n; j++) {
	    real rij2 = 0;
	    for (int k = 0; k < NDIM; k++){
		real dx = pos[i][k]-pos[j][k];
		rij2+=dx*dx;
	    }
	    if (rij2 <= rijmin2) {
	        rijmin2 = rij2;
		imin = i;
		jmin = j;
	    }
	}
    }

    real rsum = radius[imin]+radius[jmin];
    if (rsum <= 0) rsum = 1;
    real rijmin = sqrt(rijmin2);
    *sout << "closest: " << ident[imin] << " " << ident[jmin]
	  << " " << rijmin << " " << rijmin/rsum << endl << flush;
}

// New/modified integrator code:

//-----------------------------------------------------------------------------
//  evolve_system -- integrate an N-body system to time t_end.
//                   Diagnostics are sent to cout every time interval dt_dia.
//
//  Note: the integration time step, shared by all particles at any given time,
//        is variable.  Before each integration step we use coll_time (short
//        for collision time, an estimate of the time scale for any significant
//        change in configuration to happen), multiplying it by dt_param (the
//        accuracy parameter governing the size of dt in units of coll_time),
//        to obtain the new time step size.
//
//  Before moving any particles, we start with an initial diagnostics output
//  and snapshot output if desired.  In order to write the diagnostics, we
//  first have to calculate the potential energy, with get_acc_jerk_pot_coll().
//  That function also calculates accelerations, jerks, and an estimate for the
//  collision time scale, all of which are needed before we can enter the main
//  integration loop below.
//
//  The main loop takes as many integration time steps as needed to
//  reach the next output time, does the output required, then
//  continues taking integration steps and invoking output in this way
//  until the final time is reached, which triggers a `break' out of
//  the infinite loop set up with `while(true)'.
//
//-----------------------------------------------------------------------------

int evolve_system(real t_end, int sync)
{
    real epot;			// potential energy of the n-body system
    real coll_time;		// collision (close encounter) time scale

    // May be overkill to compute acc and jerk at start and end of
    // this routine, as usually the stars won't have changed on
    // return.  This way, however, we can guarantee that the particles
    // can be extrapolated when the interface calls for positions and
    // velocities.

    get_acc_jerk_pot_coll(epot, coll_time);
    real dt = calculate_step(coll_time);

    if (!init_flag) {
	write_diagnostics(epot, *sout);
	t_dia = t + dt_dia;	// next time for diagnostics output
    }

    // Don't flag a collision if no step is to be taken (presumably
    // just handled?).

    if (t + dt > t_end) return -1;

    while (true) {

	while (t < t_dia && t+dt <= t_end){
	    dt = calculate_step(coll_time);
            evolve_step(dt, epot, coll_time);	// sets t, t_evolve
	    if (test_mode) {
		real E = get_kinetic_energy() + epot;
		if (!init_flag) {
		    einit = E;
		    init_flag = true;
		}
		cout << t << " " << pos[0][0] << " " << pos[0][1]
		     << " " << E - einit << endl;
	    }
            nsteps++;

	    // compute_nn();

	    if (flag_collision && id_coll_primary >= 0) break;
        }

        if (t >= t_dia){
            write_diagnostics(epot, *sout);
            t_dia += dt_dia;
        }

	if (flag_collision && id_coll_primary >= 0) break;
        if (t+dt > t_end) break;

    }

    if (!flag_collision || id_coll_primary  < 0) {
	if (sync && t < t_end) {
	    evolve_step(t_end-t, epot, coll_time);
            nsteps++;
	}
	get_acc_jerk_pot_coll(epot, coll_time);
	dt = calculate_step(coll_time);
	t_evolve = t_end;
    }

    // Note: On exit, under all circumstances, the system is
    // synchronized at time t, with t <= t_evolve < t + dt.  If a
    // collision has been detected, we return with t_evolve = t;
    // otherwise, we set t_evolve = t_end.  If sync is set, we will
    // have t = t_end.

    return id_coll_primary;
}

int get_n_steps()
{
    return nsteps;
}



// Interface functions (headers, in order, in hermite.h):

int setup_module(bool in_reeval, bool in_test_mode)	// defaults = false
{
    reeval = in_reeval;
    test_mode = in_test_mode;
    if (test_mode) sout = &cerr;
    return 1;
}

int cleanup_module()
{
    return 1;
}

int initialize_particles(double t0)
{
    real epot, coll_time, dt;

    t = t0;
    t_evolve = t0;
    get_acc_jerk_pot_coll(epot, coll_time);
    dt = calculate_step(coll_time);

    return 0;
}

int reinitialize_particles()
{
    real epot, coll_time, dt;

    get_acc_jerk_pot_coll(epot, coll_time);
    dt = calculate_step(coll_time);

    return 0;
}


int add_particle(int id, double _mass, double _radius, double x, double y, double z, double vx, double vy, double vz)
    // add d to the dynamical system
{
    unsigned int i = 0; //find(ident.begin(), ident.end(), id) - ident.begin();
    if (0 && i < ident.size()) {

        // Particle already exists.  Do nothing.

        *sout << "add_particle: " << id
              << " already exists.  Use set_particle." << endl << flush;
    } else {

	    // Always add to the end of the list.

        ident.push_back(id);		// generally want to specify id
        mass.push_back(_mass);
        radius.push_back(_radius);
        pos.push_back(vec(x, y, z));
        vel.push_back(vec(vx, vy, vz));
        acc.push_back(vec(0,0,0));
        jerk.push_back(vec(0,0,0));
    }

    return ident.size();
}

int set_particle(int id, double _mass, double _radius, double x, double y, double z, double vx, double vy, double vz)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();

    if (i < ident.size()) {

        // Particle already exists.  Change it.

        mass[i] = _mass;
        radius[i] = _radius;
        pos[i] = vec(x, y, z);
        vel[i] = vec(vx, vy, vz);
        acc[i] = vec(0.0);
        jerk[i] = vec(0.9);

        return 0;

    } else {

        *sout << "set_particle: " << id
              << " doesn't exist.  Use add_particle." << endl << flush;
        return -1;
    }
}

int set_mass(int id, double m)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	mass[i] = m;
	return 0;
    } else 
	return -1;
}

int set_radius(int id, double r)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	radius[i] = r;
	return 0;
    } else 
	return -1;
}

int set_pos(int id, double x[])
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	pos[i] = vec(x[0], x[1], x[2]);
	return 0;
    } else 
	return -1;
}

int set_vel(int id, double v[])
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	vel[i] = vec(v[0], v[1], v[2]);
	return 0;
    } else 
	return -1;
}

int remove_particle(int id)		// remove id from the dynamical system
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	ident.erase(ident.begin()+i);
	mass.erase(mass.begin()+i);
	radius.erase(radius.begin()+i);
	pos.erase(pos.begin()+i);
	vel.erase(vel.begin()+i);
	acc.erase(acc.begin()+i);
	jerk.erase(jerk.begin()+i);
    }
    return ident.size();		// deletion never leaves empty space
}

int get_number()
{
    return ident.size();
}

double get_dynamical_time_scale()
{
    real mtot = 0, ekin = 0, epot, coll_time;
    for (unsigned int i = 0; i < ident.size(); i++) {
	mtot += mass[i];
	real dekin = 0;
	for (int k = 0; k < NDIM; k++) dekin += pow(vel[i][k],2);
	ekin += 0.5*mass[i]*dekin;
    }
    get_acc_jerk_pot_coll(epot, coll_time);

    real tdyn = (-0.5*mtot*mtot/epot) / sqrt(2*ekin/mtot);
    return tdyn;
}

real get_time()
{
    return t;
}

double get_time_step()
{
    real epot, coll_time;
    get_acc_jerk_pot_coll(epot, coll_time);
    return calculate_step(coll_time);
}

int initialize_time_step() {return 0;}
int finalize_time_step() {return 0;}

int evolve(double t_end, int sync)	// default sync = 0
{
  return evolve_system(t_end, sync);
}

int get_colliding_primary()
{
    id_coll_primary = id_coll_secondary = -1;

    int n = ident.size();
    for (int i = 0; i < n-1; i++) {	// inefficient; better to flag
					// collisions during the force loop
	for (int j = i+1; j < n; j++) {
	    real r2 = 0;
	    for(int k = 0; k < NDIM; k++){
		real dx = pos[i][k]-pos[j][k];
		r2+=dx*dx;
	    }
	    real rsum = radius[i] + radius[j];
	    if (r2 <= rsum*rsum) {
		if (mass[i] >= mass[j]) {
		    id_coll_primary = ident[i];
		    id_coll_secondary = ident[j];
		} else {
		    id_coll_primary = ident[j];
		    id_coll_secondary = ident[i];
		}
		return id_coll_primary;
	    }
	}
    }
    return -1;
}

int get_colliding_secondary(int id1)
{
    if (id1 >= 0 && id1 == id_coll_primary && id_coll_secondary >= 0)
	return id_coll_secondary;
    else
	return -1;
}

void get_state(int id, int * id_out,  double * _mass, double * _radius, double * x, double * y, double * z, double * vx, double * vy, double * vz) 
{
    *id_out = -1;

    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {

        real del = t_evolve - t;

        vec position = pos[i] + vel[i]*del + acc[i]*del*del/2
                                        + jerk[i]*del*del*del/6;
        vec velocity = vel[i] + acc[i]*del + jerk[i]*del*del/2;

        *id_out = id;
        *_mass = mass[i];
        *_radius = radius[i];
        *x = position[0];
        *y = position[1];
        *z = position[2];
        *vx = velocity[0];
        *vy = velocity[1];
        *vz = velocity[2];
    }
}

double get_mass(int id)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) 
	return mass[i];
    else
	return -1;
}

double get_radius(int id)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    //cerr << "hermite0: get_radius: "; PRC(id); PRL(i); cerr << flush;
    if (i < ident.size()) 
	return radius[i];
    else
	return -1;
}

double get_kinetic_energy()
{
    int n = ident.size();
    real ekin = 0;
    real dt = t_evolve-t;
    for (int i = 0; i < n ; i++)
	for (int k = 0; k < NDIM ; k++) {
	    real v = vel[i][k] + acc[i][k]*dt + jerk[i][k]*dt*dt/2;
	    ekin += mass[i] * v * v;
	}
    return 0.5*ekin;
}

double get_potential_energy()
{
    int n = ident.size();
    real (* save_pos)[NDIM] = new real[n][NDIM];
    real (* save_vel)[NDIM] = new real[n][NDIM];

    real dt = t_evolve-t;

    if (dt > 0) {
	for (int i = 0; i < n ; i++)
	    for (int k = 0; k < NDIM ; k++){
	        save_pos[i][k] = pos[i][k];
		save_vel[i][k] = vel[i][k];
	    }
	predict_step(t_evolve-t);
    }
    
    real epot, coll_time;
    get_acc_jerk_pot_coll(epot, coll_time);

    if (dt > 0) {
	for (int i = 0; i < n ; i++)
	    for (int k = 0; k < NDIM ; k++){
		pos[i][k] = save_pos[i][k];
		vel[i][k] = save_vel[i][k];
	    }
    }

    return epot;
}

int get_escaper()		{return -1;}	// not implemented yet
