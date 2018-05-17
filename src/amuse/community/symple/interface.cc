
// Symplectic N-body integration modules of orders 1, 2, 3, 4, 5, 6,
// 8, and 10, with fixed, shared time step.  Setting eta > 0 enables
// variable shared steps (but the integrators then are no longer
// symplectic). Individual particle continuous mass loss is also
// included.
//
//					Steve McMillan, 11/2017
//
// Added mass loss (dmdt).				 5/2018
// Expanded to multiple integrators.			 5/2018
//
// Notes:
//
// 1. Symplectic integrators can still lead to precession in an
// elliptical 2-body orbit (no conserved quantity to constrain
// it). Only the 1 and 2 integrators show any measurable precession.
//
// 2. Symplectic integrators should show no long-term growth in the
// energy of a periodic orbit (e.g. a simple two-body system of
// moderate eccentricity). However, the argument leading to this
// conclusion assumes that the period is commensurate with the period.
// Non commensurate steps lead to long-period modulations in the
// energy error, but no long-term drift. The amplitude off the
// modulation can be reduced by making the step an integral fraction
// of the period, and reduced more by making it an integral fraction
// of the actual (numerical) period, and apparently better still by
// making that integer a power of 2 (presumably reducing the rounding
// error in the solution near pericenter).  Modulation amplitude drops
// off with decreasing time step.
//
// Two-body "standard" problem results (m1 = m2 = 0.5, a = 0.57, e =
// 0.75, period P = 2.71):
//
// 1. Implicit Euler
//
// 	- no significant long-term energy drift (but large orbital
//           variations)
// 	- small amplitude modulation for longer time steps
// 	- step for ~1.e-6 max error in standard test = P/64M (est.)
//
// 2. Second-order predictor-corrector
//
// 	- max error/orbit has a long-term, apparently periodic modulation
// 	- not significantly affected by commensurate/non-commensurate step
// 	- may be driven by step not commensurate with *numerical* period
// 	- modulation form and period depends on precise choice of step
// 	- making the step commensurate with the numerical period substantially
// 	  reduces the modulation amplitude
// 	- no net drift of energy on long time scales
// 	- step for ~1.e-6 max error in standard test = P/32k (2.e-6)
// 	- no significant drift/modulation at this step
//
// 3. Third order scheme
//
// 	- step for ~1.e-6 max error in standard test = P/1024 (1.2e-6)
// 	- no significant drift/modulation at this step
// 	- power of 2 is *very* important -- substantially larger
// 	  drift if not used
//
// 4. Fourth order scheme
//
// 	- step for ~1.e-6 max error in standard test = P/1024 (1.3e-6)
// 	- no significant drift/modulation at this step
// 	- power of 2 not important for drift, but do see small
// 	  drift if not used
//
// 5. Fifth order scheme
//
// 	- step for ~1.e-6 max error in standard test = P/256 (2.5e-6)
// 	- no significant drift/modulation at this step
// 	- power of 2 is *very* important -- much larger drift/
// 	  modulation if not used (P/300 much worse than P/256)
//
// 6. Sixth order scheme
//
// 	- step for ~1.e-6 max error in standard test = P/256 (8.e-7)
// 	- no significant drift/modulation at this step
// 	- power of 2 is *very* important -- much larger drift
// 	  if not used (P/300 has substantial drift)
//
// 7. Eighth order scheme
//
// 	- step for << 1.e-6 max error in standard test = P/256 (1.5e-7)
// 	- no significant drift/modulation at this step
// 	- P/128 error >> 1.e-6 and noticeable drift
// 	- power of 2 not very important for drift; see error ~ 1.e-6
// 	  and only small drift with P/200 drift if not used (drift
// 	  per orbit is small fraction of the peak error per orbit)
//
// 8. Tenth order scheme
//
// 	- step for << 1.e-6 max error in standard test = P/128 (3.e-9)
// 	- no significant drift/modulation at this step
// 	- P/64 has substantial drift and modulation (and error
// 	  >> 1.e-6)
// 	- power of 2 is *very* important -- much larger drift
// 	  if not used (P/100 has substantial drift)
//
//				function      evals/period
//	integrator	order	evals/step	(const. dE)
//		 1	    1		 1		-- 
//		 2	    2		 1	       32k
//		 3	    3		 3		3k
//		 4	    4		 3		3k
//		 5	    5		 7	      1792
//		 6	    6		 7	      1792
//		 8	    8		15	      3840
//		10	   10		33	      4224
//
// Someone should parallelize the acceleration calculation...

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

#include "vec3.h"               // borrowed from starlab

#define PR(x)  *sout << #x << " = " << x << " " << flush
#define PRC(x) *sout << #x << " = " << x << ",  " << flush
#define PRL(x) *sout << #x << " = " << x << endl << flush

// (Global!) static data:

const int NDIM = 3;             // number of spatial dimensions
static int which_int;
static void (*integrate)(real, real*, real*);

// N-body data (structure copied from hermite0):

static real  t = 0;
static double begin_time = 0;

static vector<int>  ident;
static vector<real> mass, dmdt, radius, potential;
static vector<vec>  pos, vel, acc;

static int nsteps = 0;		// number of time steps completed
static real einit = 0;          // initial total energy of the system
static bool init_flag = false;
static real eps2 = 0;

static const real ETA = 0.05;
static real eta = ETA;		// allow variable shared time steps if > 0
static real timestep = 0;	// fixed time step off by default

static bool flag_collision = true;
static bool reeval = false;
static ostream* sout = &cout;

static real potential_energy = 0.0;
static int max_identifier = 0;

// Collisions:

static int id_coll_primary = -1, id_coll_secondary = -1;



//-----------------------------------------------------------------------
//
// Integrator functions:
//
//-----------------------------------------------------------------------

real calculate_step(real coll_time)
{
    // Determine the new system time step from coll_time.

    real step = timestep;
    if (eta > 0)
	step = eta*coll_time;

    return step;
}

// Collision detection code is copied from hermite0, but the details
// are not implemented or tested.

void get_acc_pot_coll(real *epot, real *coll_time,
		      bool get_coll=false)
{
    int n = ident.size();
    
    *coll_time = 0.0;
    *epot = 0;

    const real VERY_LARGE_NUMBER = 1e300;
    real coll_time_q = VERY_LARGE_NUMBER;      // collision time to 4th power
    real coll_est_q;                           // collision time scale estimate
                                               // to 4th power (quartic)
    id_coll_primary = id_coll_secondary = -1;

    reset_stopping_conditions();
    int is_collision_detection_enabled;
    int error = is_stopping_condition_enabled(COLLISION_DETECTION, 
					      &is_collision_detection_enabled);

    for (int i = 0; i < n ; i++) {
        for (int k = 0; k < NDIM; k++) acc[i][k] = 0;
        potential[i]= 0;
    }
 
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            real rji[NDIM];
            real vji[NDIM];
            for (int k = 0; k < NDIM ; k++) {
                rji[k] = pos[j][k] - pos[i][k];
                vji[k] = vel[j][k] - vel[i][k];
	    }
            real r2 = 0;
            real v2 = 0;
            for (int k = 0; k < NDIM ; k++) {
                r2 += rji[k] * rji[k];
                if (get_coll) v2 += vji[k] * vji[k];
	    }
 
            if (is_collision_detection_enabled) {	// not working/tested
		real rsum = radius[i] + radius[j];
		if (r2 <= rsum*rsum) {
		    int stopping_index  = next_index_for_stopping_condition();
		    if(stopping_index >= 0) {
			set_stopping_condition_info(stopping_index,
						    COLLISION_DETECTION);
			set_stopping_condition_particle_index(stopping_index,
							      0, ident[i]);
			set_stopping_condition_particle_index(stopping_index,
							      1, ident[j]);
		    }
		}
            }
            
            r2 += eps2;
            real r = sqrt(r2);
            real r3 = r * r2;

            // Add the {i,j} contribution to the total potential energy.
	    
            *epot -= mass[i] * mass[j] / r;
            potential[i] = -mass[j] / r;
            potential[j] = -mass[i] / r;
	    
            // Add the {j (i)} contribution to the {i (j)} acc.
	    
            real da[NDIM];
            for (int k = 0; k < NDIM ; k++) {
                da[k] = rji[k] / r3;
	    }
            for (int k = 0; k < NDIM ; k++) {
                acc[i][k]  += mass[j] * da[k];
                acc[j][k]  -= mass[i] * da[k];
	    }

	    if (get_coll) {

		// First collision time estimate is based on unaccelerated
		// linear motion.
	    
		coll_est_q = (r2*r2) / (v2*v2);
		if (coll_time_q > coll_est_q)
		    coll_time_q = coll_est_q;
	    
		// Second collision time estimate is based on free-fall.
	    
		real da2 = 0;
		for (int k = 0; k < NDIM ; k++)
		    da2 += da[k] * da[k];
	    
		real mij = mass[i] + mass[j];
		da2 *= mij * mij;
	    
		coll_est_q = r2/da2;
		if (coll_time_q > coll_est_q)
		    coll_time_q = coll_est_q;
	    }
	}
    }

    if (get_coll)
	*coll_time = pow(coll_time_q, 0.25);
}

inline void step_symp(int n, int kk, real *c, real *d, real dt,
		      real *epot, real *coll_time)
{
    // Generic code to take a symplectic step (see 3, 4, 5, 6, 8
    // below).  Note that, in *all* cases, sum(c) = sum(d) = 1.  The
    // check is below, but the result is important in determining how
    // to distribute mass loss.

#if 0
    static int count = 0;
    if (count == 0) {
	real ctot = 0, dtot = 0;
	for (int k = 0; k < kk; k++) {
	    ctot += c[k];
	    dtot += d[k];
	}
	PRC(ctot); PRL(dtot);
	count = 1;
    }
#endif

    for (int j = 0; j < kk; j++) {
	if (d[j] != 0) {
	    real ddt = d[j]*dt;
	    get_acc_pot_coll(epot, coll_time, (eta > 0 && j == kk-1));
	    for (int i = 0; i < n ; i++)
		for (int k = 0; k < NDIM ; k++)
		    vel[i][k] += ddt*acc[i][k];
	}
	real cdt = c[j]*dt;
	for (int i = 0; i < n ; i++) {
	    mass[i] += cdt*dmdt[i];		// mass follows pos
	    for (int k = 0; k < NDIM ; k++)
		pos[i][k] += cdt*vel[i][k];
	}
    }
    t += dt;
}

void evolve_step1(real dt, real *epot, real *coll_time)
{
    // Implicit Euler scheme.

    static const int k1 = 1;
    static real c1[k1] = {1};
    static real d1[k1] = {1};

    step_symp(ident.size(), k1, c1, d1, dt, epot, coll_time);
}

void evolve_step2(real dt, real *epot, real *coll_time)
{
    // Second-order predictor-corrector scheme.  Assume acc is already
    // set on entry.
    
    int n = ident.size();
    real (*old_acc)[NDIM] = new real[n][NDIM];

    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            old_acc[i][k] = acc[i][k];
    
    // Predict.
    
    for (int i = 0; i < n ; i++) {
	mass[i] += dmdt[i]*dt;
        for (int k = 0; k < NDIM ; k++)
            pos[i][k] += (vel[i][k] + acc[i][k]*dt/2) * dt;
    }

    get_acc_pot_coll(epot, coll_time, (eta > 0));

    // Correct.  Note that old_acc used mass at the start of the step,
    // acc uses mass at the end.
    
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            vel[i][k] += (old_acc[i][k] + acc[i][k]) * dt/2;

    t += dt;
    delete[] old_acc;
}

void evolve_step3(real dt, real *epot, real *coll_time)
{
    // Third-order symplectic scheme (Ruth 1983).

    static const int k3 = 3;
    static real c3[k3] = {2./3, -2./3, 1.0};
    static real d3[k3] = {7./24, 3./4, -1./24};

    step_symp(ident.size(), k3, c3, d3, dt, epot, coll_time);
}

void evolve_step4(real dt, real *epot, real *coll_time)
{
    // Fourth-order symplectic scheme (Yoshida 1990).

    static const int k4 = 4;
    static const real p = pow(2., 1./3);
    static const real q = 1 - p;
    static const real r = 4 - 2*p;
    static real c4[k4] = {1/r, q/r, q/r, 1/r};
    static real d4[k4] = {0., 2/r, -2*p/r, 2/r};

    step_symp(ident.size(), k4, c4, d4, dt, epot, coll_time);
}

void evolve_step5(real dt, real *epot, real *coll_time)
{
    // Fifth-order symplectic scheme (Tselios & Simos 2013).
    // Optimized for orbital applications.

    static const int k5 = 7;
    static real c5[k5]
	= {-0.402020995028838599420412333241250172914690575978880873429,
	    0.345821780864741783378055242038676806930765132085822482512,
	    0.400962967485371350147918025877657753577504227492190779513,
	    0.980926531879316517259793318227431991923428491844523669724,
	   -1.362064898669775624786044007840908597402026042205084284026,
	    0.923805029000837468447500070054064432491178527428114178991,
	    0.112569584468347104973189684884327785393840239333314075493};
    static real d5[k5]
	= { 0.47501834514453949720351208570106713494289203770372938037,
	    0.021856594741098449005512783774683495267598355789295971623,
	   -0.334948298035883491345320878224434762455516821029015086331,
	    0.51263817465269673604202785657395553607442158325539698102,
	   -0.011978701020553903586622444048386301410473649207894475166,
	   -0.032120004263046859169923904393901683486678946201463277409,
	    0.36953388878114957185081450061701658106775743968995046842};

    step_symp(ident.size(), k5, c5, d5, dt, epot, coll_time);
}

void initialize_yosh(int mm, real *w, int kk, real *c, real *d)
{
    // Yoshida initialization of high-order integration coefficients.
    // Array w is of length mm; arrays c and d are of length kk.
    
    real w0 = 1;
    for (int m = 0; m < mm; m++) w0 -= 2*w[m];

    c[0] =c[kk-1] = w[mm-1]/2;
    for (int m = 1; m < mm; m++)
	c[m] =  c[kk-m-1] = (w[mm-m]+w[mm-m-1])/2;
    c[mm] = c[mm+1] = (w[0]+w0)/2;

    // First element of dd is undefined by Yoshida -- 0 seems to work.
    
    for (int m = 1; m <= mm; m++)
	d[m] = d[kk-m] = w[mm-m];
    d[mm+1] = w0;
}

void evolve_step6(real dt, real *epot, real *coll_time)
{
    // Sixth-order symplectic scheme (Yoshida 1990).  Yoshida gives
    // three solutions; the first (A) is most (x10) accurate.

    static const int m6 = 3;
    static real w6[m6]
	= {-0.117767998417887e1, 0.235573213359357e0, 0.784513610477560e0};
    // w6 = {-0.213228522200144e1, 0.426068187079180e-2, 0.143984816797678e1};
    // w6 = {0.152886228424922e-2, -0.214403531630539e1, 0.144778256239930e1};
    
    static const int k6 = 2*m6+2;
    static real c6[k6] = {0};
    static real d6[k6] = {};
    if (c6[0] == 0) initialize_yosh(m6, w6, k6, c6, d6);

    step_symp(ident.size(), k6, c6, d6, dt, epot, coll_time);
}

void evolve_step8(real dt, real *epot, real *coll_time)
{
    // Eighth-order symplectic (Yoshida 1990).  As listed, they are
    // approximately in order of increasing accuracy.  Choice 1
    // (Yoshida solution A) is substantially less accurate than the
    // others, and is actually less accurate than the 6th-order scheme
    // for reasonable step choices.  Choices 4 and 5 (Yoshida
    // solutions D and E) are roughly 1000 times more accurate than
    // choice 1 for my standard test problem (Duffing, dx = 0.01).
    // Choice 4 is 10 times more accurate for longer steps (dx = 0.1 -
    // 0.03)

    static const int m8 = 7;

    // w8 = {-0.161582374150097e1, -0.244699182370524e1, -0.716989419708120e-2,
    //        0.244002732616735e1,  0.157739928123617e0,  0.182020630970714e1,
    //        0.10424262086999e1}

    // w8 = {-0.169248587770116e-2, 0.289195744315849e1,  0.378039588360192e-2,
    //       -0.289688250328827e1,  0.289105148970595e1, -0.233864815101035e1,
    //        0.148819229202922e1}

    // w8 = { 0.311790812418427e0, -0.155946803821447e1, -0.167896928259640e1,
    //        0.166335809963315e1, -0.106458714789183e1,  0.136934946416871e1,
    //        0.629030650210433e0}

    static real w8[m8]
	= { 0.102799849391985e0, -0.196061023297549e1,  0.193813913762276e1,
           -0.158240635368243e0, -0.144485223686048e1,  0.253693336566229e0,
            0.914844246229740e0};

    // w8 = { 0.227738840094906e-1,  0.252778927322839e1, -0.719180053552772e-1,
    //        0.536018921307285e-2, -0.204809795887393e1,  0.107990467703699e0,
    //        0.130300165760014e1}

    static const int k8 = 2*m8+2;
    static real c8[k8] = {0};
    static real d8[k8] = {};
    if (c8[0] == 0) initialize_yosh(m8, w8, k8, c8, d8);

    step_symp(ident.size(), k8, c8, d8, dt, epot, coll_time);
}

void evolve_step10(real dt, real *epot, real *coll_time)
{
    // Tenth-order symplectic (Tsitouras 1998).

    static const int m10 = 16;
    static real w10[m10] = {
	 0.02690013604768968151437422144441685467297755661,
	 0.939801567135683337900037741418097674632505563,
	-0.00803583920385358749646880826318191806393661063,
	-0.866485197373761372803661767454208401679117010,
	 0.1023112911193598731078563285067131541328142449,
	-0.1970772151393080101376018465105491958660525085,
	 0.617877713318069357335731125307019691019646679,
	 0.1907272896000121001605903836891198270441436012,
	 0.2072605028852482559382954630002620777969060377,
	-0.395006197760920667393122535979679328161187572,
	-0.582423447311644594710573905438945739845940956,
	 0.742673314357319863476853599632017373530365297,
	 0.1643375495204672910151440244080443210570501579,
	-0.615116639060545182658778437156157368647972997,
	 0.2017504140367640350582861633379013481712172488,
	 0.45238717224346720617588658607423353932336395045};

    static const int k10 = 2*m10+2;
    static real c10[k10] = {0};
    static real d10[k10] = {};
    if (c10[0] == 0) initialize_yosh(m10, w10, k10, c10, d10);

    step_symp(ident.size(), k10, c10, d10, dt, epot, coll_time);
}

int evolve_system(real t_end)
{
    real epot;                     // potential energy of the n-body system
    real coll_time;                // collision (close encounter) time scale
    int n = ident.size();
    int i, k;

    get_acc_pot_coll(&epot, &coll_time);

    while (t < t_end) {
	real dt = calculate_step(coll_time);
	if (t+dt > t_end) dt = t_end - t;
	integrate(dt, &epot, &coll_time);
    }

    return 0;
}

int get_n_steps()
{
    return nsteps;
}

int cleanup_code()
{
    reset_stopping_conditions();

    ident.clear();
    vel.clear();
    pos.clear();
    mass.clear();
    acc.clear();
    potential.clear();
    radius.clear();
    
    begin_time = 0.0;
    max_identifier = 0;
    nsteps = 0;
    eta = ETA;

    einit = 0;
    init_flag = false;
    eps2 = 0;
   
    id_coll_primary = -1;
    id_coll_secondary = -1;

    flag_collision = true;
    reeval = false;

    potential_energy = 0.0;
    t = 0;

    return 0;
}


//-----------------------------------------------------------------------
//
// AMUSE interface functions:
//
//-----------------------------------------------------------------------

int get_integrator(int *_i)
{
    *_i = which_int;
    return 0;
}

void set_integration_scheme() {
    if (which_int == 1)
	integrate = &evolve_step1;
    else if (which_int == 2)
	integrate = &evolve_step2;
    else if (which_int == 3)
	integrate = &evolve_step3;
    else if (which_int == 4)
	integrate = &evolve_step4;
    else if (which_int == 5)
	integrate = &evolve_step5;
    else if (which_int == 6)
	integrate = &evolve_step6;
    else if (which_int == 8)
	integrate = &evolve_step8;
    else if (which_int == 10)
	integrate = &evolve_step10;
    else {
	cout << "warning: unknown integrator; using evolve_step2"
	     << endl << flush;
	integrate = &evolve_step2;
    }
}

int set_integrator(int _i)
{
    which_int = _i;
    set_integration_scheme();
    return 0;
}

int get_eps2(double *_epsilon_squared)
{
    *_epsilon_squared = eps2;
    return 0;
}

int set_eps2(double _epsilon_squared)
{
    eps2 = _epsilon_squared;
    return 0;
}

int get_eta(double *_eta)
{
  *_eta = eta;
  return 0;
}

int set_eta(double _eta)
{
  eta = _eta;
  timestep = 0;
  return 0;
}

int get_timestep(double *_timestep)
{
  *_timestep = timestep;
  return 0;
}

int set_timestep(double _timestep)
{
  timestep = _timestep;
  eta = 0;
  return 0;
}

int get_time(double *_t)
{
  *_t = t;
  return 0;
}

int set_begin_time(double input)
{
    begin_time = input;
    return 0;
}

int get_begin_time(double * output)
{
    *output = begin_time;
    return 0;
}

int new_particle(int *id, double _mass,
                 double x, double y, double z,
                 double vx, double vy, double vz, double _radius)
{
    int new_element;

    // Always add to the end of the list.

    new_element = max_identifier++;
    ident.push_back(new_element);	// generally want to specify id

    mass.push_back(_mass);
    dmdt.push_back(0.0);		// default = 0; set with set_dmdt
    radius.push_back(_radius);
    pos.push_back(vec(x, y, z));
    vel.push_back(vec(vx, vy, vz));
    acc.push_back(vec(0,0,0));
    potential.push_back(0);

    *id = new_element;

    return 0;
}

int delete_particle(int id)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();

    if (i < ident.size()) {
	ident.erase(ident.begin()+i);
	mass.erase(mass.begin()+i);
	dmdt.erase(mass.begin()+i);
	radius.erase(radius.begin()+i);
	pos.erase(pos.begin()+i);
	vel.erase(vel.begin()+i);
	acc.erase(acc.begin()+i);
	potential.erase(potential.begin()+i);
	return 0;
    } else
	return -1;
}

int get_state(int id, double *_mass, double *x, double *y, double *z,
              double *vx, double *vy, double *vz, double *_radius)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  
    if (i < (int) ident.size()) {
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
    } else
	return -1;
}

int set_state(int id, double _mass, double x, double y, double z,
              double vx, double vy, double vz, double _radius)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  
    if (i < ident.size()) {
	mass[i] = _mass;
	radius[i] = _radius;
	pos[i] = vec(x,y,z);
	vel[i] = vec(vx, vy, vz);
	return 0;
    } else
	return -1;
}

int get_mass(int id, double *_mass)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
        *_mass = mass[i];
        return 0;
    } else
        return -1;
}

int set_mass(int id, double _mass)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
        mass[i] = _mass;
        return 0;
    } else
        return -1;
}

int get_dmdt(int id, double *_dmdt)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
        *_dmdt = dmdt[i];
        return 0;
    } else
        return -1;
}

int set_dmdt(int id, double _dmdt)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
        dmdt[i] = _dmdt;
        return 0;
    } else
	return -1;
}

int get_radius(int id, double *_radius)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
        *_radius = radius[i];
        return 0;
    } else
        return -1;
}

int set_radius(int id, double _radius)
{
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size()) {
      radius[i] = _radius;
      return 0;
  } else
      return -1;
}

int get_position(int id, double *x, double *y, double *z)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	*x = pos[i][0];
	*y = pos[i][1];
	*z = pos[i][2];
	return 0;
    } else
	return -2;
}

int set_position(int id, double x, double y, double z)
{
  unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
  if (i < ident.size()) {
      pos[i] = vec(x, y, z);
      return 0;
  } else
      return -1;
}

int get_velocity(int id, double *vx, double *vy, double *vz)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	*vx = vel[i][0];
	*vy = vel[i][1];
	*vz = vel[i][2];
	return 0;
    } else
	return -2;
}

int set_velocity(int id, double vx, double vy, double vz)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	vel[i] = vec(vx,vy,vz);
	return 0;
    } else
	return -1;
}

int get_acceleration(int id, double *ax, double *ay, double *az)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	*ax = acc[i][0];
	*ay = acc[i][1];
	*az = acc[i][2];
	return 0;
    } else
	return -2;
}

int set_acceleration(int id, double ax, double ay, double az)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	acc[i][0] = ax;
	acc[i][1] = ay;
	acc[i][2] = az;
	return 0;
    } else
	return -2;
}

int evolve_model(double t_end)
{
    evolve_system(t_end);
    return 0;
}

int initialize_code()
{
    which_int = 2;
    set_integration_scheme();
    begin_time = 0.0;

    initialize_stopping_conditions();
        
    // AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(COLLISION_DETECTION);
    set_support_for_condition(TIMEOUT_DETECTION);
    set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    set_support_for_condition(OUT_OF_BOX_DETECTION);
    // -----------------------

    return 0;
}

int get_kinetic_energy(double *kinetic_energy)
{
    int n = ident.size();
    real ekin = 0;
    for (int i = 0; i < n ; i++) {
        for (int k = 0; k < NDIM ; k++) {
            real v = vel[i][k];
            ekin += mass[i] * v * v;
        }
    }
    *kinetic_energy = 0.5*ekin;
    return 0;
}

int get_potential_energy(double *value)
{
    int n = ident.size();
    real epot, coll_time;
    get_acc_pot_coll(&potential_energy, &coll_time);
    *value = potential_energy;
    
    return 0;
}

int get_potential(int id, double *value)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();
    if (i < ident.size()) {
	*value = potential[i];
	return 0;
    } else {
	*value = 0.0;
	return -1;
    }
}

int get_potential_at_point(double eps,
			   double x, double y, double z,
			   double *phi)
{
    int n = ident.size();
    double rx,ry,rz,r;
    *phi = 0.0;

    for (int i = 0; i< n; i++) {
        rx = pos[i][0]-x;
        ry = pos[i][1]-y;
        rz = pos[i][2]-z;
        r = sqrt(rx*rx+ry*ry+rz*rz + eps2);
        *phi -= mass[i]/r;
    }

    return 0;
}

int get_gravity_at_point(double eps, 
			 double x, double y, double z,
			 double *ax, double *ay, double *az)
{
    int n = ident.size();
    double rx, ry, rz, r3, r2, r, F;

    *ax = 0;
    *ay = 0;
    *az = 0;

    for (int i = 0; i<n; i++) {
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
{
    int n = ident.size();
    *_mass = 0.0;

    for (int i = 0; i< n; i++) {
        *_mass += mass[i];
    }

    return 0;
}

int get_total_radius(double *_radius)
{
    int n = ident.size();
    *_radius = 0.0;

    return 0;
}

int get_center_of_mass_position(double *x, double *y, double *z)
{
    int n = ident.size();
    double M = 0;

    *x=0; *y=0; *z=0;

    get_total_mass(&M);

    for (int i = 0; i<n; i++) {
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
    int n = ident.size();
    double M = 0;

    *vx=0; *vy=0; *vz=0;

    get_total_mass(&M);

    for (int i = 0; i<n; i++) {
        *vx += mass[i]*vel[i][0];
        *vy += mass[i]*vel[i][1];
        *vz += mass[i]*vel[i][2];
    }

    *vx /= M;
    *vy /= M;
    *vz /= M;

    return 0;
}

int get_number_of_particles(int *no_parts)
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

    if (i < ident.size()-1) {
	*index_of_the_next_particle = ident.at(i+1);
	return 0;
    } else {
	if (i == (ident.size()-1)) {
	    return 1;
	} else {
	    return -1;
	}
    }
}

int set_particle(int id, double _mass, double _radius,
		 double x, double y, double z,
		 double vx, double vy, double vz)
{
    unsigned int i = find(ident.begin(), ident.end(), id) - ident.begin();

    if (i < ident.size()) {

        // Particle already exists.  Change it.
	
        mass[i] = _mass;
	dmdt[i] = 0.0;
        radius[i] = _radius;
        pos[i] = vec(x, y, z);
        vel[i] = vec(vx, vy, vz);
        acc[i] = vec(0.0);
        potential[i] = 0;
        return 0;
    } else {
        *sout << "set_particle: " << id
              << " doesn't exist.  Use add_particle." << endl << flush;
        return -1;
    }
}

int get_dynamical_time_scale(double *ts)
{
    real mtot = 0, ekin = 0, epot, coll_time;
    for (unsigned int i = 0; i < ident.size(); i++) {
        mtot += mass[i];
        real dekin = 0;
        for (int k = 0; k < NDIM; k++) dekin += pow(vel[i][k],2);
        ekin += 0.5*mass[i]*dekin;
    }

    get_acc_pot_coll(&potential_energy, &coll_time);

    real tdyn = (-0.5*mtot*mtot/epot) / sqrt(2*ekin/mtot);
    return tdyn;
}

int get_time_step(double *time_step)
{
    real epot, coll_time;
    get_acc_pot_coll(&potential_energy, &coll_time);
    *time_step = calculate_step(coll_time);
    return 0;
}

int initialize_time_step() {return -2;}

int finalize_time_step() {return -2;}

int set_reeval(int value)
{
    reeval = value;
    return 0;
}

int get_reeval(int * value)
{
    *value = reeval;
    return 0;
}

int recommit_particles()
{
    real epot, coll_time;
    get_acc_pot_coll(&potential_energy, &coll_time);

    return 0;
}

int recommit_parameters()
{
    real epot, coll_time;
    get_acc_pot_coll(&potential_energy, &coll_time);

    return 0;
}

int commit_particles()
{
    real epot, coll_time;
    get_acc_pot_coll(&potential_energy, &coll_time);

    return 0;
}

int commit_parameters()
{
    t = begin_time;
    return 0;
}

int synchronize_model()
{
    return 0;
}
