#ifndef  NBODY_PARTICLE_H
#  define   NBODY_PARTICLE_H
/*-----------------------------------------------------------------------------
 *  nbody-particle : basic class for simple nbody implementation
 *  J. Makino 1998/11/29
 *-----------------------------------------------------------------------------
 */
#include "stdinc.h"
#include "vec.h"

#ifndef ONED
#define THREED
#endif
#define REAL_GRAVITY

class nbody_particle
    {
    private:
        vec pos;
	vec vel;
	vec acc_gravity;
	real phi_gravity;
	real mass;
	real radius;			// added by SLWM, 6/07
	int index;

    public:
	nbody_particle(){
	    pos = 0.0;
	    vel = 0.0;
	    acc_gravity = 0.0;
	    phi_gravity = mass = radius = 0.0;
	    index = 0;
	}
        void  set_pos(const vec& new_pos)      {pos = new_pos;}
        void  set_vel(const vec& new_vel)      {vel = new_vel;}
        void  set_acc_gravity(const vec& new_acc)      {acc_gravity = new_acc;}
        void  set_phi_gravity(real new_phi)      {phi_gravity = new_phi;}

	void  clear_pos()                         {pos = 0.0;}
	void  clear_vel()                         {vel = 0.0;}
	void  clear_acc_phi_gravity()      {acc_gravity = 0.0;phi_gravity = 0.0;}
	void  correct_phi_self_gravity(real epsinv)      {phi_gravity += mass*epsinv;}

	void  inc_pos(const vec& d_pos)        {pos += d_pos; }
	void  inc_vel(const vec& d_vel)        {vel += d_vel;}
	void  update_vel(real dt)        {vel = dt*acc_gravity;}
	void  update_pos(real dt)        {pos = (pos+dt*vel).readjust();}

	void  scale_pos(const real scale_factor)  {pos *= scale_factor; }
	void  scale_vel(const real scale_factor)  {vel *= scale_factor; }

	vec  get_pos()                         {return pos;}
	vec  get_vel()                         {return vel;}
	real  get_phi_gravity()                         {return phi_gravity;}
	vec  get_acc_gravity()                         {return acc_gravity;}

	real get_mass()			{return mass;}
	void set_mass(real m)		{mass = m;}


	// added by SLWM, 6/07:

	real get_radius()		{return radius;}
	void set_radius(real r)		{radius = r;}


	void set_index(int i){index = i;}
	int get_index(){return index;}
	void predict(real dt){
	    real dt2 = dt*dt*0.5;
	    pos = (pos + dt*vel + dt2*acc_gravity).readjust();
	    vel +=  (dt*0.5)*acc_gravity;
	}
	void correct(real dt){
	    vel +=  (dt*0.5)*acc_gravity;
	}

	void read(istream & );
	void write(ostream & );
	void dump();

	void plot(real parm);

	real kinetic_energy();
	real energy();
	real get_ke(){ return 0.5*mass*vel*vel;}

	void friend accumulate_mutual_gravity(nbody_particle & p1,
					      nbody_particle & p2,
					      real eps2);
	void calculate_gravity_using_tree(real eps2, real theta2);
	void correct_gravity(real);

};

typedef vec (nbody_particle::*nbody_VMF_ptr)(void);     // vec member function pointer
typedef void (nbody_particle::*nbody_MF_ptr)(const vec &);     // member function pointer

typedef void (nbody_particle::*nbody_VF_ptr)(void);     // void member function pointer
typedef void (nbody_particle::*nbody_RF_ptr)(real);     // void member function
						    // pointer with real arg
typedef void (nbody_particle::*nbody_RRF_ptr)(real,real);     // void member function
						    // pointer with two real args
class nbody_system
    {
    private:


	int nsize;
	nbody_particle * pb;

    public:
	int n;
	real time;
	real timestep;
	real eps2_for_gravity;
	int use_self_gravity;
	vec pos;
	vec vel;
	real   mass;
	real plot_xmax;
	real theta_for_tree;
	int ncrit_for_tree;
	    

	nbody_system(){
	    n = 0;
	    nsize = 0;
	    time = 0;
	    timestep = 0;
	    pb = NULL;
	    }
 
	void calculate_uncorrected_gravity();
	void calculate_uncorrected_gravity_using_grape4();
	void calculate_uncorrected_gravity_direct();

	void set_nsize(int n) {nsize = n;}

	void read(istream & );
	void write(ostream & );
	void atos(istream & );
	void stoa(ostream & );
	void dump();
	void  setup_tree();
	void generate_cube(int nx);
	void apply_vf(nbody_VF_ptr);
	void apply_vf(nbody_RF_ptr, real);
	void apply_vf(nbody_RRF_ptr, real, real);
	void plot(real time);
	real kinetic_energy();
	real energy();
	void show_energy();
	void calculate_gravity();
	void create_uniform_sphere(int nbody, real power_index, real r0);

	void evolve( real dt, real tend);
	void evolve_onestep( real dt);
	void integrate( real dt);

	void calculate_cmterms();
	void make_collision(vec relpos, vec relv);
	nbody_particle * get_particle_pointer(){return pb;}
	void friend copy_nbody_particles(nbody_system * source,
				       nbody_system * desitination);
	void correct_gravity();

	// routine added by SPZ (Nov 2006)
	void set_particle_pointer(nbody_particle * np){pb = np;}
};

#endif


