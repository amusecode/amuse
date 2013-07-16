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

#ifndef NOMPI
#include <mpi.h>
#endif

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>

#include "worker_code.h"

#include <vector>
#include <algorithm>
#include <string>
#include <map>


// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>
#include <time.h>


using namespace std;
typedef double  real;

#include "vec3.h"                // borrowed from starlab

#define PR(x)  *sout << #x << " = " << x << " " << flush
#define PRC(x) *sout << #x << " = " << x << ",  " << flush
#define PRL(x) *sout << #x << " = " << x << endl << flush



// (Global!) static data:

const int NDIM = 3;                // number of spatial dimensions

/*------------- MPI data ---------------------*/
int mpi_rank = 0;
int mpi_size = 1;
/*------------- MPI data ---------------------*/

// N-body data:


real  t = 0;
real t_evolve = t;                // Time requested by evolve.  Control returns
                                // when t <= t_evolve and t + dt > t_evolve,
                                // and this is assumed when the state of the
                                // system is computed by extrapolation.
real t_wanted = 0;
static double begin_time = 0;
static double end_time_accuracy_factor = 0.5;

vector<int>  ident;
vector<real> mass, radius, potential;
vector<vec>  pos, vel, acc, jerk;
vector<vec>  acc_reduced, jerk_reduced;
vector<real> potential_reduced;

// Control parameters:

const real DT_PARAM = 0.03;
const real DT_DIA = 1;

real dt_param = DT_PARAM;        // control parameter to determine time step size
real dt_dia = DT_DIA;                // time interval between diagnostic output

bool x_flag = false;                // set true for serious debugging only!
bool is_time_reversed_allowed = false;

int nsteps = 0;                        // number of integration time steps completed
real einit = 0;                        // initial total energy of the system
bool init_flag = false;
real t_dia = 0;
real eps2 = 0;

bool flag_collision = true;
bool reeval = false;
bool test_mode = false;
ostream* sout = &cout;

// Accessors for use by the C++ main driver only:

//void set_t(real tt)                {t = tt;}
//void set_dt_param(real dt)        {dt_param = dt;}
//void set_dt_dia(real dt)        {dt_dia = dt;}
//void set_eps(real eps)                 {eps2 = eps*eps;}

int get_eps2(double *_epsilon_squared){
    *_epsilon_squared = eps2;
    return 0;
}
int set_eps2(double _epsilon_squared){
    eps2 = _epsilon_squared;
    return 0;
}
int get_dt_dia(double *_dt_dia){
    *_dt_dia = dt_dia;
    return 0;
}
int set_dt_dia(double _dt_dia){
    dt_dia = _dt_dia;
    return 0;
}

int get_is_time_reversed_allowed(int *value){
    *value = is_time_reversed_allowed ? 1 : 0;
    return 0;
}
int set_is_time_reversed_allowed(int value){
    is_time_reversed_allowed = value == 1;
    return 0;
}


int get_dt_param(double *_dt_param)
{
  *_dt_param = dt_param;
  return 0;
}
int set_dt_param(double _dt_param)
{
  dt_param = _dt_param;
  return 0;
}
int get_time(double *_t)
{
  *_t = t;
  return 0;
}


int set_begin_time(double input) {
    begin_time = input;
    return 0;
}

int get_begin_time(double * output) {
    *output = begin_time;
    return 0;
}


// (Most of) the original code (some arguments replaced by global data):

//-----------------------------------------------------------------------------
//  write_diagnostics  --  writes diagnostics on the output stream cout:
//                         current time; total mass; number of
//                         integration steps so far; kinetic,
//                         potential, and total energy; absolute and
//                         relative energy errors since the start of
//                         the run.  If x_flag (x for eXtra data) is
//                         true, all internal data are dumped for each
//                         particle (mass, position, velocity,
//                         acceleration, and jerk).
//
//  Note: the kinetic energy is calculated here, while the potential
//  energy is calculated in the function get_acc_jerk_pot_coll().
//-----------------------------------------------------------------------------

void write_diagnostics(real epot, ostream& s = cout)
{
    int n = ident.size();
    real total_mass = 0;
    real ekin = 0;                        // kinetic energy
    for (int i = 0; i < n ; i++) {
        total_mass += mass[i];
        for (int k = 0; k < NDIM ; k++)
            ekin += 0.5 * mass[i] * vel[i][k] * vel[i][k];
    }

    real etot = ekin + epot;                // total energy

    if (!init_flag)                     // on first pass, set
      {
        einit = etot;                   // the initial energy
        init_flag = true;
        return;                         // suppress initial output
      }

    s << "    internal diagnostics at time t = " << t
      << " after " << nsteps << " steps"
      << endl
      << "        total mass = " << total_mass
      << "  initial energy E_init = " << einit << endl
      << "        E_kin = " << ekin
      << "  E_pot = " << epot
      << "  E_tot = " << etot << endl;
    s << "        "
      << "absolute energy error  E_tot - E_init = "
      << etot - einit << endl;
    s << "        "
      << "relative energy error  (E_tot - E_init) / E_init = "
      << (etot - einit) / einit << endl;

    if (x_flag)
      {
        s << "        system dump, n = " << n << endl;
        for (int i = 0; i < n ; i++)
          {
            s << "        data for particle " << ident[i]
	      << ": " << endl;
            s << "            "; s << mass[i] << endl;
            s << "            "; s << radius[i] << endl;
            s << "           ";
            for (int k = 0; k < NDIM; k++)
              {
                s << ' ' << pos[i][k];
              }
            s << endl;
            s << "           ";
            for (int k = 0; k < NDIM; k++)
              {
                s << ' ' << vel[i][k];
              }
            s << endl;
            s << "           ";
            for (int k = 0; k < NDIM; k++)
              {
                s << ' ' << acc[i][k];
              }
            s << endl;
            s << "           ";
            for (int k = 0; k < NDIM; k++)
              {
                s << ' ' << jerk[i][k];
              }
            s << endl;
          }
      }
}

//-----------------------------------------------------------------------------
//  predict_step  --  take the first approximation of one Hermite integration
//                    step, advancing the positions and velocities through a
//                    Taylor series development up to the order of the jerks.
//                    Note that all pos and vel are overwritten.
//-----------------------------------------------------------------------------

void predict_step(real dt)
{
    if (!is_time_reversed_allowed && dt <= 0) return;

    int n = ident.size();
    for (int i = 0; i < n ; i++)
      {
        for (int k = 0; k < NDIM ; k++)
          {
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
//          collision time in a separate function.  However, the current
//          function is by far the most time consuming part of the whole
//          program, with a double loop over all particles that is executed
//          every time step.  Splitting off some of the work to another
//          function would significantly increase the total computer time
//          (by an amount close to a factor two).
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

inline int mpi_distribute_data(int n) {

#ifndef NOMPI
    MPI_Bcast(&n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    
    if(mpi_rank) {
        vel.resize(n+1);
        pos.resize(n+1);
        acc.resize(n+1);
        jerk.resize(n+1);
        potential.resize(n+1, 0.0);
        ident.resize(n+1);
        mass.resize(n+1, 0.0);
        radius.resize(n+1, 0.0);
    }
    
    
    MPI_Bcast(&vel[0], n * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&pos[0], n * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mass[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&radius[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ident[0], n, MPI_INTEGER, 0, MPI_COMM_WORLD);
    
    mpi_distribute_stopping_conditions();
#endif

    return n;
}

inline void mpi_collect_data(int n, real *epot, real *coll_time_q_out, real coll_time_q_in) {

#ifndef NOMPI
    real summed = 0.0;
    MPI_Reduce(&coll_time_q_in, &summed, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    *coll_time_q_out = summed;
    
    summed = 0.0;
    MPI_Reduce(epot, &summed, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    *epot = summed;
    
    if(!mpi_rank) {
        acc_reduced.resize(n+1);
        jerk_reduced.resize(n+1);
        potential_reduced.resize(n+1);
    }
    
    MPI_Reduce(&acc[0], &acc_reduced[0], n * 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&jerk[0], &jerk_reduced[0], n * 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potential[0], &potential_reduced[0], n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(!mpi_rank) {
        acc = acc_reduced;
        jerk = jerk_reduced;
        potential = potential_reduced;
    }
    
    mpi_collect_stopping_conditions();
#endif
}

void get_acc_jerk_pot_coll(real *epot, real *coll_time)
{
    int n = 0;
    int error;
    int is_collision_detection_enabled;

    if(mpi_rank == 0){
        n = ident.size();
    }
    
    n = mpi_distribute_data(n);
    
    for (int i = 0; i < n ; i++)
      {
        for (int k = 0; k < NDIM ; k++)
          {
            acc[i][k] = jerk[i][k] = 0;
          }
        potential[i]= 0;
      }

    *epot = 0;
    const real VERY_LARGE_NUMBER = 1e300;
    real coll_time_q = VERY_LARGE_NUMBER;      // collision time to 4th power
    real coll_est_q;                           // collision time scale estimate
                                               // to 4th power (quartic)
    id_coll_primary = id_coll_secondary = -1;

    reset_stopping_conditions();
    //int npart = n / mpi_size;
    //int nleft = n % mpi_size;
    //if(mpi_rank == mpi_size - 1) {iend = n;}
    //cerr << istart <<", "<<iend<<" , "<<nleft<<endl;
    
    // divide the work over every process, division of work should be somewhat equal
    int istart = 0 + mpi_rank;
    int iend = n ; //(mpi_rank + 1) * npart;
    
    error = is_stopping_condition_enabled(
                COLLISION_DETECTION, 
                &is_collision_detection_enabled
    );
    for (int i = istart; i < iend ; i+= mpi_size)
      {
        for (int j = i+1; j < n ; j++)             // rji[] is the vector from
          {
            real rji[NDIM];                        // particle i to particle j
            real vji[NDIM];                        // vji[] = d rji[] / d t
            for (int k = 0; k < NDIM ; k++)
              {
                rji[k] = pos[j][k] - pos[i][k];
                vji[k] = vel[j][k] - vel[i][k];
              }
            real r2 = 0;                           // | rji |^2
            real v2 = 0;                           // | vji |^2
            real rv_r2 = 0;                        // ( rij . vij ) / | rji |^2
            for (int k = 0; k < NDIM ; k++)
              {
                r2 += rji[k] * rji[k];
                v2 += vji[k] * vji[k];
                rv_r2 += rji[k] * vji[k];
              }
            rv_r2 /= r2;
 
            if(is_collision_detection_enabled) {  
              real rsum = radius[i] + radius[j];
              if (r2 <= rsum*rsum) {
                int stopping_index  = next_index_for_stopping_condition();
                if(stopping_index < 0)
                {
                    
                }
                else
                {
                    set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
                    set_stopping_condition_particle_index(stopping_index, 0, ident[i]);
                    set_stopping_condition_particle_index(stopping_index, 1, ident[j]);
                }
              }
            }
        
            
            r2 += eps2;                            // | rji |^2 + eps^2
            real r = sqrt(r2);                     // | rji |
            real r3 = r * r2;                      // | rji |^3

            // Add the {i,j} contribution to the total potential energy.

            *epot -= mass[i] * mass[j] / r;
            potential[i] = -mass[j] / r;
            potential[j] = -mass[i] / r;

            // Add the {j (i)} contribution to the {i (j)} acc and jerk.

            real da[NDIM];                         // main terms in pairwise
            real dj[NDIM];                         // acceleration and jerk
            for (int k = 0; k < NDIM ; k++)
              {
                da[k] = rji[k] / r3;                          // see equations
                dj[k] = (vji[k] - 3 * rv_r2 * rji[k]) / r3;   // in the header
              }
            for (int k = 0; k < NDIM ; k++)
              {
                acc[i][k]  += mass[j] * da[k];                // using symmetry
                acc[j][k]  -= mass[i] * da[k];                // find pairwise
                jerk[i][k] += mass[j] * dj[k];                // acceleration
                jerk[j][k] -= mass[i] * dj[k];                // and jerk
              }

            // First collision time estimate is based on unaccelerated
            // linear motion.

            coll_est_q = (r2*r2) / (v2*v2);
            if (coll_time_q > coll_est_q)
              {
                coll_time_q = coll_est_q;
              }
            // Second collision time estimate is based on free fall.

            real da2 = 0;                                  // da2 becomes the
            for (int k = 0; k < NDIM ; k++)                // square of the
              {
                da2 += da[k] * da[k];                      // pairwise accel-
              }
            real mij = mass[i] + mass[j];                  // eration between
            da2 *= mij * mij;                              // particles i and j

            coll_est_q = r2/da2;
            if (coll_time_q > coll_est_q)
              {
                coll_time_q = coll_est_q;
              }
          }
      }
    mpi_collect_data(n, epot, &coll_time_q, coll_time_q);
    
                                                 // from q for quartic back
     *coll_time = sqrt(sqrt(coll_time_q));            // to linear collision time
}

//-----------------------------------------------------------------------------
//  correct_step  --  take one iteration to improve the new values of
//                    position and velocities, constructing a higher-order
//                    Taylor series from the terms up to jerk at the
//                    beginning and the end of the time step.  This symmetric
//                      formulation is not the original one from Makino and
//                      Aarseth; it comes from ACS (Hut & Makino).
//-----------------------------------------------------------------------------

void correct_step(const real old_pos[][NDIM], const real old_vel[][NDIM],
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  real dt)
{
    int n = ident.size();
    for (int i = 0; i < n ; i++)
      {
        for (int k = 0; k < NDIM ; k++)
          {
            vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2
                                      + (old_jerk[i][k] - jerk[i][k])*dt*dt/12;
            pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2
                                      + (old_acc[i][k] - acc[i][k])*dt*dt/12;
          }
      }
}

//-----------------------------------------------------------------------------
//  evolve_step  --  take one integration step for an N-body system, using
//                   the Hermite algorithm.
//-----------------------------------------------------------------------------

void evolve_step(real dt, real *epot, real *coll_time)
{
    int n = ident.size();
    real (* old_pos)[NDIM] = new real[n][NDIM];
    real (* old_vel)[NDIM] = new real[n][NDIM];
    real (* old_acc)[NDIM] = new real[n][NDIM];
    real (* old_jerk)[NDIM] = new real[n][NDIM];

    for (int i = 0; i < n ; i++)
    {
        for (int k = 0; k < NDIM ; k++)
        {
            old_pos[i][k] = pos[i][k];
            old_vel[i][k] = vel[i][k];
            old_acc[i][k] = acc[i][k];
            old_jerk[i][k] = jerk[i][k];
        }
    }

    predict_step(dt);
    get_acc_jerk_pot_coll(epot, coll_time);
    
/*****
 * put back in previous positions if a collision was detected
 * for experiment, need to remove
    if(set_conditions & enabled_conditions) {
        for (int i = 0; i < n ; i++)
        {
            for (int k = 0; k < NDIM ; k++)
            {
                pos[i][k] = old_pos[i][k];
                vel[i][k] = old_vel[i][k];
                acc[i][k] = old_acc[i][k];
                jerk[i][k] = old_jerk[i][k];
            }
        }
        delete[] old_pos;
        delete[] old_vel;
        delete[] old_acc;
        delete[] old_jerk;
        return;
    }
 *
 ****/
    correct_step(old_pos, old_vel, old_acc, old_jerk, dt);
    
    if (reeval)
    {
        get_acc_jerk_pot_coll(epot, coll_time);
    }
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
    if (!test_mode)
      {

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

    for (int i = 0; i < n-1; i++)
      {
        for (int j = i+1; j < n; j++)
          {
            real rij2 = 0;
            for (int k = 0; k < NDIM; k++)
              {
                real dx = pos[i][k]-pos[j][k];
                rij2+=dx*dx;
              }
            if (rij2 <= rijmin2)
              {
                rijmin2 = rij2;
                imin = i;
                jmin = j;
              }
          }
      }

    real rsum = radius[imin]+radius[jmin];

    if (rsum <= 0)
      {
        rsum = 1;
      }
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

int evolve_system(real t_end)
{
    real epot;                     // potential energy of the n-body system
    real coll_time;                // collision (close encounter) time scale
    int nest_err = 0;              // error of subprocedure
    int must_run = 1;

    int is_number_of_steps_detection_enabled;
    int is_out_of_box_detection_enabled;
    int number_of_steps_innerloop = 0;
    int max_number_of_steps;
    double box_size;
    double sqr_distance_wrt_origin;
    double total, center_of_mass[NDIM];
    int error;
    int n = ident.size();//no particles
    int n_particles_out_of_box = 0;
    int i, k; 
    // May be overkill to compute acc and jerk at start and end of
    // this routine, as usually the stars won't have changed on
    // return.  This way, however, we can guarantee that the particles
    // can be extrapolated when the interface calls for positions and
    // velocities.
    
#ifndef NOMPI
    MPI_Bcast(&must_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif

    get_acc_jerk_pot_coll(&epot, &coll_time);
    
    bool is_time_reversed = false;
    
    if(t_end < t && is_time_reversed_allowed) {
        is_time_reversed = true;
    }
    
    real dt = calculate_step(coll_time);
    
    if(is_time_reversed) {
        dt *= -1;
    }
    
    std::cout<<"t"<<t<<", DT:"<<dt<<", coll_time:"<<coll_time<<std::endl;
    t_wanted = t_end;
    if (!init_flag)
      {
        write_diagnostics(epot, *sout);
        t_dia = t + dt_dia;        // next time for diagnostics output
      }
      
    
    if(end_time_accuracy_factor == 0.0) {
        if(!is_time_reversed && t < t_end && t + dt > t_end) {
            dt = t_end - t;
        } 
        else if(is_time_reversed && t_end < t  &&  t_end > t + dt ) {
            dt = t_end - t;
        }
    }
   
    std::cout<<"t0"<<t<<", DT:"<<dt<<", coll_time:"<<coll_time<<std::endl;
    
    if (
        (!is_time_reversed && t + dt > t_end + (end_time_accuracy_factor * dt))
        ||
        
        (is_time_reversed && t_end - (end_time_accuracy_factor * dt) >  t + dt )
    )
    {
        must_run = 0;
#ifndef NOMPI
        MPI_Bcast(&must_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
        return -1;
    
    }
    
    
    // AMUSE STOPPING CONDITIONS
    time_t starttime, currenttime;
    time(&starttime);

    int timeout_detection;
    //
    error = is_stopping_condition_enabled(TIMEOUT_DETECTION, 
                    &timeout_detection);
    error = is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION, 
                    &is_number_of_steps_detection_enabled);
    error = is_stopping_condition_enabled(OUT_OF_BOX_DETECTION,
                    &is_out_of_box_detection_enabled);
    get_stopping_condition_number_of_steps_parameter(&max_number_of_steps);    
    get_stopping_condition_out_of_box_parameter(&box_size);    
    // AMUSE STOPPING CONDITIONS
    
    while (true) {
        while (
            (!is_time_reversed && t < t_dia && t + dt <= t_end + (end_time_accuracy_factor*dt))
            ||
            (is_time_reversed && t_end - (end_time_accuracy_factor*dt) < t + dt )
        ) {
            
            #ifndef NOMPI
            MPI_Bcast(&must_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
            #endif

            evolve_step(dt, &epot, &coll_time);        // sets t, t_evolve

            dt = calculate_step(coll_time);
            
            if(is_time_reversed) {
                dt *= -1;
            }
            
            std::cout<<"t1"<<t<<", DT:"<<dt<<", coll_time:"<<coll_time<<std::endl;
            if(
                (!is_time_reversed && end_time_accuracy_factor == 0.0 && t < t_end && t + dt > t_end)
                ||
                (is_time_reversed && end_time_accuracy_factor == 0.0 && t_end < t &&  t_end > t + dt )
            ) {
                dt = t_end - t;
            }
            
            std::cout<<"t2"<<t<<", DT:"<<dt<<", coll_time:"<<coll_time<<std::endl;
            if (test_mode) {
                real E = 0.0;
                nest_err = get_kinetic_energy(&E);
                E += epot;
                if (!init_flag) {
                    einit = E;
                    init_flag = true;
                }
                cout << t << " " << pos[0][0] << " " << pos[0][1]
                 << " " << E - einit << endl;
            }
            nsteps++;

            // compute_nn();

            // AMUSE STOPPING CONDITIONS
            if(timeout_detection) {
                time(&currenttime);
                cerr << currenttime << " : " << starttime << " : " << timeout_parameter << " : " << (currenttime - starttime) << endl;
                if((currenttime - starttime) > timeout_parameter) {
                    int stopping_index  = next_index_for_stopping_condition();
                    set_stopping_condition_info(stopping_index, TIMEOUT_DETECTION);
                }
            }
            
            if(is_number_of_steps_detection_enabled) {
                number_of_steps_innerloop++;
                if (number_of_steps_innerloop > max_number_of_steps) {
                  int stopping_index  = next_index_for_stopping_condition();
                  set_stopping_condition_info(stopping_index, 
                              NUMBER_OF_STEPS_DETECTION);
                }
            }
            if (is_out_of_box_detection_enabled) {
                for (k = 0; k < NDIM; k++) {
                    center_of_mass[k] = 0.0;
                }
                total = 0.0;
                for (i = 0; i < n; i++) {
                    for (k = 0; k < NDIM; k++) {
                        center_of_mass[k] += mass[i] * pos[i][k] * pos[i][k];
                        total += mass[i];
                    }
                }
                for (k = 0; k < NDIM; k++) {
                    center_of_mass[k] /= total;
                }
            }
            if (is_out_of_box_detection_enabled) {
                for (i = 0; i < n; i++) {
                    sqr_distance_wrt_origin = 0.0;
                    for (k = 0; k < NDIM; k++) {
                        sqr_distance_wrt_origin += (pos[i][k] - center_of_mass[k])*(pos[i][k] - center_of_mass[k]);
                    }
                    if (sqr_distance_wrt_origin > box_size*box_size) {
                        int stopping_index = next_index_for_stopping_condition();
                        set_stopping_condition_info(stopping_index, OUT_OF_BOX_DETECTION);
                        if (n_particles_out_of_box < 10) {
                            set_stopping_condition_particle_index(stopping_index,
                                                  n_particles_out_of_box,
                                                  ident[i]);
                            n_particles_out_of_box++;
                        }
                        else {
                            printf("Run out of storable out of box events\n");
                        }
                    }
                }
            }
            
            if(set_conditions & enabled_conditions) {
                break;
            }
        }

        if (t >= t_dia) {
          write_diagnostics(epot, *sout);
          t_dia += dt_dia;
        }

            //is_stopping_condition_set???
        if (set_conditions & enabled_conditions) {
          break;
        }
        
        if (
                (!is_time_reversed && t + dt >= t_end + (end_time_accuracy_factor * dt))
                ||
                (is_time_reversed && t_end - (end_time_accuracy_factor * dt) >= t + dt )
        ) {
          break;
        }
    }
    
    if (!(set_conditions & enabled_conditions))
      {

#ifndef NOMPI
        MPI_Bcast(&must_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
        get_acc_jerk_pot_coll(&epot, &coll_time);
        dt = calculate_step(coll_time);
        if(is_time_reversed) {
            dt *= -1;
        }
        t_evolve = t;
      } else {
        t_evolve = t;
      }

    // Note: On exit, under all circumstances, the system is
    // synchronized at time t, with t <= t_evolve < t + dt.  If a
    // collision has been detected, we return with t_evolve = t;
    // otherwise, we set t_evolve = t_end. 
    
    
    must_run = 0;
#ifndef NOMPI
    MPI_Bcast(&must_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif

    return nest_err;
}

int get_n_steps()
{
  return nsteps;
}


static int max_identifier = 0;

int cleanup_code()
{
    if(mpi_rank) {
        return 0;
    }
    
    vel.clear();
    pos.clear();
    mass.clear();
    acc.clear();
    jerk.clear();
    potential.clear();
    radius.clear();
    max_identifier = 0;
    return 0;
}



int new_particle(int *id, double _mass,
                 double x, double y, double z,
                 double vx, double vy, double vz, double _radius)
// add d to the dynamical system
// cello, proj 1
// make test, not done yet
{
    
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    
  int new_element;

  // Always add to the end of the list.
  // Allways at the end anyway

  new_element = max_identifier++;

  ident.push_back(new_element);                // generally want to specify id
  mass.push_back(_mass);
  radius.push_back(_radius);
  pos.push_back(vec(x, y, z));
  vel.push_back(vec(vx, vy, vz));
  acc.push_back(vec(0,0,0));
  jerk.push_back(vec(0,0,0));
  potential.push_back(0);

  *id = new_element;

  //return ident.size();
  return 0;
}

int delete_particle(int id)
{
    
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();

  if (i < ident.size())
    {
      ident.erase(ident.begin()+i);
      mass.erase(mass.begin()+i);
      radius.erase(radius.begin()+i);
      pos.erase(pos.begin()+i);
      vel.erase(vel.begin()+i);
      acc.erase(acc.begin()+i);
      jerk.erase(jerk.begin()+i);
      potential.erase(potential.begin()+i);
      return 0;
    }
  else
    {
      return -1;
    }
}


/*
int remove_particle(int id)                // remove id from the dynamical system
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
    return ident.size();                // deletion never leaves empty space
}
*/

int get_state(int id, double *_mass, double *x, double *y, double *z,
              double *vx, double *vy, double *vz, double *_radius)
// cello, proj1 changed return type void-->int, only OK no errors yet?
// todo: replace find fction.
{
  if(mpi_rank) {return 0;}
  //*id_out = -1;

  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  
  if (i < (int) ident.size())
    {
      vec position = pos[i];
      vec velocity = vel[i];

      //*id_out = id;
      *_mass = mass[i];
      *_radius = radius[i];
      *x = position[0];
      *y = position[1];
      *z = position[2];
      *vx = velocity[0];
      *vy = velocity[1];
      *vz = velocity[2];
      return 0;
    }
  else
    {
      return -1;
    }
}

int set_state(int id, double _mass, double x, double y, double z,
              double vx, double vy, double vz, double _radius)
//cello, proj1,
{
  if(mpi_rank) {return 0;}
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  
  if (i < ident.size())
    {
      mass[i] = _mass;
      radius[i] = _radius;
      pos[i] = vec(x,y,z);
      vel[i] = vec(vx, vy, vz);
      return 0;
    }
  else
    {
      return -1;
    }

}

int get_mass(int id, double *_mass)
//cello, proj1,
{
    if(mpi_rank) {return 0;}
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size())
      {
        *_mass = mass[i];
        return 0;
      }
    else
      {
        return -1;
      }
}

int set_mass(int id, double _mass)
{
    if(mpi_rank) {return 0;}
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size())
      {
        mass[i] = _mass;
        return 0;
      } else
      return -1;
}

int get_radius(int id, double *_radius)
{
    if(mpi_rank) {return 0;}
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    //cerr << "hermite0: get_radius: "; PRC(id); PRL(i); cerr << flush;
    if (i < ident.size())
      {
        *_radius = radius[i];
        return 0;
      }
    else
      {
        return -1;
      }
}

int set_radius(int id, double _radius)
{
    if(mpi_rank) {return 0;}
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size())
    {
      radius[i] = _radius;;
      return 0;
    }
  else
    {
      return -1;
    }
}

int get_position(int id, double *x, double *y, double *z)
{
    if(mpi_rank) {return 0;}
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size())
    {
      *x = pos[i][0];
      *y = pos[i][1];
      *z = pos[i][2];
      return 0;
    }
  else
    {
      return -2;
    }
}

int set_position(int id, double x, double y, double z)
{
    if(mpi_rank) {return 0;}
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size())
    {
      pos[i] = vec(x, y, z);
      return 0;
    }
  else
    {
      return -1;
    }
}

int get_velocity(int id, double *vx, double *vy, double *vz)
{
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size())
    {
      *vx = vel[i][0];
      *vy = vel[i][1];
      *vz = vel[i][2];
      return 0;
    }
  else
    {
      return -2;
    }
}

int set_velocity(int id, double vx, double vy, double vz)
{
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();

  if (i < ident.size())
    {
      vel[i] = vec(vx,vy,vz);
      return 0;
    }
  else
    {
      return -1;
    }
  return 0;
}

int get_acceleration(int id, double *ax, double *ay, double *az)
{
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size())
    {
      *ax = acc[i][0];
      *ay = acc[i][1];
      *az = acc[i][2];
      return 0;
    }
  else
    {
      return -2;
    }
}

int set_acceleration(int id, double ax, double ay, double az)
{
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size())
    {
      acc[i][0] = ax;
      acc[i][1] = ay;
      acc[i][2] = az;
      return 0;
    }
  else
    {
      return -2;
    }

  return -3;
}

int evolve_not_on_root() {

#ifndef NOMPI
    int must_run = 1;
    int mpi_error = MPI_Bcast(&must_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    if(mpi_error) {
        return -1;
    }
    while(must_run) {
        real epot, coll_time;
        get_acc_jerk_pot_coll(&epot, &coll_time);
        int mpi_error = MPI_Bcast(&must_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
        if(mpi_error) {
            return -1;
        }
    }
#endif
    return 0;
}

int evolve_model(double t_end)
{
    if(mpi_rank)     {
        evolve_not_on_root();
    } else {
        evolve_system(t_end);
    }
    return 0;
}


int initialize_code()
{
    begin_time = 0.0;
    is_time_reversed_allowed = false;
    
#ifndef NOMPI
    int error = 0;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(error) {
        cerr << error << endl;
        return -1;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if(error) {
        cerr << error << endl;
        return -1;
    }
#else
    mpi_rank = 0;
    mpi_size = 1;
#endif
    //cerr <<"mpi rank: "<<mpi_rank<<", mpi size: "<<mpi_size<<endl;
    
    // AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(COLLISION_DETECTION);
    set_support_for_condition(TIMEOUT_DETECTION);
    set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    set_support_for_condition(OUT_OF_BOX_DETECTION);
    // -----------------------

#ifndef NOMPI
    mpi_setup_stopping_conditions();
#endif

    return 0;
}


int get_kinetic_energy(double *kinetic_energy)
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    int n = ident.size();
    real ekin = 0;
    real dt = t_evolve-t;
    for (int i = 0; i < n ; i++)
    {
        for (int k = 0; k < NDIM ; k++)
        {
            real v = vel[i][k] + acc[i][k]*dt + jerk[i][k]*dt*dt/2;
            ekin += mass[i] * v * v;
        }
    }
    *kinetic_energy = 0.5*ekin;
    return 0;
}

int get_potential_energy(double *potential_energy)
{
    int n = 0;
    real (* save_pos)[NDIM] = 0;
    real (* save_vel)[NDIM] = 0;

    real dt = t_evolve-t;
    if(mpi_rank == 0) {
        n = ident.size();
        save_pos = new real[n][NDIM];
        save_vel = new real[n][NDIM];
        dt = t_evolve-t;

        if (dt > 0)
        {
          for (int i = 0; i < n ; i++)
            {
              for (int k = 0; k < NDIM ; k++)
                {
                  save_pos[i][k] = pos[i][k];
                  save_vel[i][k] = vel[i][k];
                }
            }
          predict_step(t_evolve-t);
        }
    }
    real epot, coll_time;
    get_acc_jerk_pot_coll(&epot, &coll_time);
    if(mpi_rank == 0) {
        if (dt > 0)
        {
            for (int i = 0; i < n ; i++)
              {
                for (int k = 0; k < NDIM ; k++)
                  {
                    pos[i][k] = save_pos[i][k];
                    vel[i][k] = save_vel[i][k];
                  }
              }
        }
        delete[] save_pos;
        delete[] save_vel;
    }
    *potential_energy = epot;
    return 0;
}

int get_potential(int id,  double *value)
{
	if(mpi_rank) {return 0;}
	unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
	if (i < ident.size())
	  {
		*value = potential[i];
		return 0;
	  }
	else
	  {
		*value = 0.0;
		return -1;
	  }
}

int get_potential_at_point(double eps, double x, double y, double z, double *phi)
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    int n = ident.size();
    double rx,ry,rz,r;
    *phi = 0.0;

    for (int i = 0; i< n; i++)
    {
        rx = pos[i][0]-x;
        ry = pos[i][1]-y;
        rz = pos[i][2]-z;
        r = sqrt(rx*rx+ry*ry+rz*rz + eps2);
        *phi -= mass[i]/r;
    }

    return 0;
}

int get_gravity_at_point(
    double eps, 
    double x,double y, double z,
    double *ax, double *ay, double *az
    )
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    int n = ident.size();
    double rx, ry, rz, r3, r2, r, F;

    *ax = 0;
    *ay = 0;
    *az = 0;

    for (int i = 0; i<n; i++)
    {
        rx = pos[i][0]-x;
        ry = pos[i][1]-y;
        rz = pos[i][2]-z;
        r2 = (rx*rx+ry*ry+rz*rz + eps2);
        r = sqrt(r2);
        r3 = r2*r;
        F = mass[i]/r3;
        *ax += F * rx;
        *ay += F * ry;
        *az += F * rz;
    }

    return 0;
}



int get_total_mass(double *_mass)
// cello, proj1
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    int n = ident.size();
    *_mass = 0.0;

    for (int i = 0; i< n; i++)
    {
        *_mass += mass[i];
    }

    return 0;
}

int get_center_of_mass_position(double *x, double *y, double *z)
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    int n = ident.size();
    double M = 0;

    *x=0; *y=0; *z=0;

    get_total_mass(&M);

    for (int i = 0; i<n; i++)
    {
        *x += mass[i]*pos[i][0];
        *y += mass[i]*pos[i][1];
        *z += mass[i]*pos[i][2];
    }

    *x /= M;
    *y /= M;
    *z /= M;

    return 0;
}

int get_center_of_mass_velocity(double *vx, double *vy,double *vz)
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    int n = ident.size();
    double M = 0;

    *vx=0; *vy=0; *vz=0;

    get_total_mass(&M);

    for (int i = 0; i<n; i++)
    {
        *vx += mass[i]*vel[i][0];
        *vy += mass[i]*vel[i][1];
        *vz += mass[i]*vel[i][2];
    }

    *vx /= M;
    *vy /= M;
    *vz /= M;

    return 0;
}


int get_total_radius(double *radius)
{
  // not implemented yet
  return -1;
}

int get_number_of_particles(int *no_parts)
//cello proj1
//used to be get_number assuming this is OK
{
  *no_parts = ident.size();
  return 0;
}

int get_index_of_first_particle(int * index_of_first_particle)
{
  *index_of_first_particle = ident.front();
  return 0;
}

int get_index_of_next_particle(int id, int *index_of_the_next_particle)
{

  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();

  if (i < ident.size()-1)
    {
      *index_of_the_next_particle = ident.at(i+1);
      return 0;
    }
  else
    {
      if (i == (ident.size()-1))
        {
          return 1;
        }
      else
        {
          return -1;
        }
    }

}

int set_particle(int id, double _mass, double _radius, double x, double y, double z, double vx, double vy, double vz)
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();

    if (i < ident.size())
      {
        // Particle already exists.  Change it.
        mass[i] = _mass;
        radius[i] = _radius;
        pos[i] = vec(x, y, z);
        vel[i] = vec(vx, vy, vz);
        acc[i] = vec(0.0);
        jerk[i] = vec(0.9);
        potential[i] = 0;
        return 0;
      }
    else
      {
        *sout << "set_particle: " << id
              << " doesn't exist.  Use add_particle." << endl << flush;
        return -1;
      }
}

int get_dynamical_time_scale(double *ts)
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    real mtot = 0, ekin = 0, epot, coll_time;
    for (unsigned int i = 0; i < ident.size(); i++)
      {
        mtot += mass[i];
        real dekin = 0;
        for (int k = 0; k < NDIM; k++) dekin += pow(vel[i][k],2);
        ekin += 0.5*mass[i]*dekin;
      }

    get_acc_jerk_pot_coll(&epot, &coll_time);

    real tdyn = (-0.5*mtot*mtot/epot) / sqrt(2*ekin/mtot);
    return tdyn;
}

int get_time_step(double *time_step)
{
    if(mpi_rank)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    real epot, coll_time;
    get_acc_jerk_pot_coll(&epot, &coll_time);
    *time_step = calculate_step(coll_time);
    return 0;
}

int initialize_time_step() {return -2;}

int finalize_time_step() {return -2;}

int get_escaper(){return -2;}   // not implemented yet



int set_reeval(int value) {
    reeval = value;
    return 0;
}

int get_reeval(int * value) {
    *value = reeval;
    return 0;
}

int set_testmode(int value) {
    test_mode = value;
    return 0;
}

int get_testmode(int * value) {
    *value = test_mode;
    return 0;
}

int recommit_particles(){
    real epot, coll_time;

    get_acc_jerk_pot_coll(&epot, &coll_time);

    return 0;
}

int recommit_parameters(){
    real epot, coll_time;

    get_acc_jerk_pot_coll(&epot, &coll_time);

    return 0;
}


int commit_particles()
{
    real epot, coll_time;
    get_acc_jerk_pot_coll(&epot, &coll_time);
    return 0;
}

int commit_parameters(){
    
    t = begin_time;
    t_evolve = t;
    t_wanted = t;
    t_dia = begin_time;
    
    if(test_mode) {
        sout = &cerr;
    }
    return 0;
}

int synchronize_model() {
    return 0;
}

int set_end_time_accuracy_factor(double value)
{
    if(value < -1.0 || value > 1.0) {
        return -1;
    }
    end_time_accuracy_factor = value;
    return 0;
}
int get_end_time_accuracy_factor(double * value)
{
    *value = end_time_accuracy_factor;
    return 0;
}
