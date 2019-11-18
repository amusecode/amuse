#include "interface.h"

// A stub of this file is machine generated, but the content is
// hand-coded.  SAVE A COPY to avoid accidental overwriting!

#include "src/smallN_interface.h"

static hdyn *b = NULL;		// root node for all smallN data
static hdyn *b_copy = NULL;

static double begin_time = 0;
static real smalln_dtlog = _INFINITY_;
static int smalln_verbose = 0;
static string smalln_outfile("");

class UpdatedParticle {

  // AMUSE bookkeeping.  A particle on the UpdatedParticle list has
  // been added or removed internally.  Use this class to communicate
  // this information back to the python level.
    
  public:

    int index_of_particle;
    int status;			// 1 ==> deleted, 2 ==> added
    
    UpdatedParticle():index_of_particle(-1),status(0) {}
    UpdatedParticle(int i, int s):index_of_particle(i), status(s) {}
    UpdatedParticle(const UpdatedParticle& src)
	:index_of_particle(src.index_of_particle), status(src.status) {}
};

#include <vector>
vector<UpdatedParticle> UpdatedParticles;

// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>

// Setup and parameters.

int cleanup_code()
{
    // Clean up at the end of the calculation.

    if (b) {
        rmtree(b);		// deletes b
        b = NULL;
    }
    if (b_copy) {
        rmtree(b_copy);		// deletes b_copy
        b_copy = NULL;
    }
    begin_time = 0.0;
    real smalln_dtlog = _INFINITY_;
    smalln_verbose = 0;
    smalln_outfile = string("");
    UpdatedParticles.clear();
    return 0;
}

int initialize_code()
{
    initialize_stopping_conditions();    
    // Begin the initialization by creating a basic hdyn data structure.

    cleanup_code();
    
    b = new hdyn;
    b_copy = NULL;
    
    b->set_eta(0.14);
    b->set_gamma(1e-6);
    b->set_system_time(0.0);
    b->set_allow_full_unperturbed(1);
    b->set_cm_index(100000);
    
    // AMUSE stopping conditions support.

    set_support_for_condition(COLLISION_DETECTION);
    set_support_for_condition(INTERACTION_OVER_DETECTION);
    mpi_setup_stopping_conditions();

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

int set_eps2(double softening_parameter_sq)		// not used
{
    return 0;
}

int get_eps2(double * softening_parameter_sq)
{
    *softening_parameter_sq = 0;
    return 0;
}

int set_eta(double timestep_parameter)
{
    b->set_eta(timestep_parameter);
    return 0;
}

int get_eta(double * timestep_parameter)
{
    *timestep_parameter = b->get_eta();
    return 0;
}

int set_gamma(double unperturbed_parameter)
{
    b->set_gamma(unperturbed_parameter);
    return 0;
}

int get_gamma(double * unperturbed_parameter)
{
    *unperturbed_parameter = b->get_gamma();
    return 0;
}

int set_allow_full_unperturbed(int allow_full_unperturbed)
{
    b->set_allow_full_unperturbed(allow_full_unperturbed);
    return 0;
}

int get_allow_full_unperturbed(int * allow_full_unperturbed)
{
    *allow_full_unperturbed = b->get_allow_full_unperturbed();
    return 0;
}

int set_cm_index(int cm_index)
{
    b->set_cm_index(cm_index);
    // cout << "interface set_cm_index " << cm_index << endl << flush;
    return 0;
}

int get_cm_index(int * cm_index)
{
    *cm_index = b->get_cm_index();
    return 0;
}

int set_time(double sys_time)
{
    b->set_system_time(sys_time);
    return 0;
}

int get_time(double * sys_time)
{
    *sys_time = b->get_system_time();
    return 0;
}

int set_outfile(char * file)
{
    smalln_outfile = string(file);
    // cout << "file = " << file << endl;
    // cout << "smalln_outfile = " << smalln_outfile << endl;
    return 0;
}

int get_outfile(char ** file)
{
    *file = (char*)smalln_outfile.c_str();
    return 0;
}

int commit_parameters()
{
    // Perform any needed setup after initial code parameters have been set.

    if (b->get_system_time() == 0) {
        b->set_system_time(begin_time);
    }
    return 0;
}

int recommit_parameters()
{
    // Perform any needed changes after code parameters have been reset.

    return 0;
}

int commit_particles()
{
    // Complete the initialization, after all particles have been loaded.

    return 0;
}

int recommit_particles()
{
    // Reinitialize/reset the system after particles have been added
    // or removed.

    return 0;
}

// Setters and getters for individual particles.

int new_particle(int * index_of_the_particle,
		 double mass, 
		 double x, double y, double z,
		 double vx, double vy, double vz,
		 double radius, int index_to_set)
{
    // Add a particle to the system.  Let the module set the id, or
    // force the index to index_to_set if >= 0 and allowed.

    *index_of_the_particle = add_particle(b, mass, radius,
					  vec(x,y,z), vec(vx,vy,vz),
					  index_to_set);
    //cout << "new_particle: " << index_to_set
    //	 << " " << *index_of_the_particle
    //	 << endl << flush;

    return 0;
}

int delete_particle(int index_of_the_particle)
{
    // Remove a (top-level) particle from the system.

    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    else return remove_particle(bb);
}

int get_index_of_next_particle(int index_of_the_particle, 
			       int * index_of_the_next_particle)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb)
	return -1;
    else if (!bb->get_younger_sister())
	return 1;
    else
	 *index_of_the_next_particle = bb->get_younger_sister()->get_index();
    return 0;
}

int set_state(int index_of_the_particle,
	      double mass, 
	      double x, double y, double z,
	      double vx, double vy, double vz,
	      double radius)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    bb->set_mass(mass);
    bb->set_radius(radius);
    bb->set_pos(vec(x,y,z));
    bb->set_vel(vec(vx,vy,vz));
    return 0;
}

int get_state(int index_of_the_particle,
	      double * mass, 
	      double * x, double * y, double * z,
	      double * vx, double * vy, double * vz,
	      double * radius)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    *mass = bb->get_mass();
    *radius = bb->get_radius();
    vec pos = bb->get_pos(), vel = bb->get_vel();
    hdyn *p = bb->get_parent();
    while (p && p->get_parent()) {
	pos += p->get_pos();		// return absolute values
	vel += p->get_vel();
	p = p->get_parent();
    }
    *x = pos[0];
    *y = pos[1];
    *z = pos[2];
    *vx = vel[0];
    *vy = vel[1];
    *vz = vel[2];
    return 0;
}

int set_mass(int index_of_the_particle, double mass)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    bb->set_mass(mass);
    return 0;
}

int get_mass(int index_of_the_particle, double * mass)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    *mass = bb->get_mass();
    return 0;
}

int set_radius(int index_of_the_particle, double radius)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    bb->set_radius(radius);
    return 0;
}

int get_radius(int index_of_the_particle, double * radius)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb){*radius = 0.0; return -1;}
    *radius = bb->get_radius();
    return 0;
}

int set_position(int index_of_the_particle,
		 double x, double y, double z)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    bb->set_pos(vec(x,y,z));
    return 0;
}

int get_position(int index_of_the_particle,
		 double * x, double * y, double * z)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    vec pos = bb->get_pos();
    hdyn *p = bb->get_parent();
    while (p && p->get_parent()) {
	pos += p->get_pos();		// return absolute value
	p = p->get_parent();
    }
    *x = pos[0];
    *y = pos[1];
    *z = pos[2];
    return 0;
}

int set_velocity(int index_of_the_particle,
		 double vx, double vy, double vz)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    bb->set_vel(vec(vx,vy,vz));
    return 0;
}

int get_velocity(int index_of_the_particle,
		 double * vx, double * vy, double * vz)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    vec vel = bb->get_vel();
    hdyn *p = bb->get_parent();
    while (p && p->get_parent()) {
	vel += p->get_vel();		// return absolute value
	p = p->get_parent();
    }
    *vx = vel[0];
    *vy = vel[1];
    *vz = vel[2];
    return 0;
}

int set_acceleration(int index_of_the_particle,
		     double ax, double ay, double az)
{
    return -1;
}

int get_acceleration(int index_of_the_particle,
		     double * ax, double * ay, double * az)
{
    return -1;
}

int get_potential(int index_of_the_particle, double * pot)
{
    return -1;
}


// System-wide operations.

static real break_scale_sq = _INFINITY_;
static real structure_check_interval = _INFINITY_;

int set_break_scale(real r) {
    break_scale_sq = r*r;
    if (break_scale_sq <= 0) break_scale_sq = _INFINITY_;
    return 0;
}

int set_structure_check_interval(real dt) {
    structure_check_interval = dt;
    return 0;
}

void myspaces(int n) {for (int i = 0; i < n; i++) cout << " ";}

void myprint(hdyn *b, int level = 0)
{
    myspaces(4*level);
    cout << b->get_index() << "  " << "mass = " << b->get_mass();
    if (b->get_kepler()) cout << "  kepler";
    cout << endl;
    myspaces(4*level);
    cout << "    pos = " << b->get_pos() << endl;
    for_all_daughters(hdyn, b, bb)
	myprint(bb, level+1);
}

int evolve_model(double to_time)
{
    // On return, system_time will be greater than or equal to the
    // specified time.  All particles j will have time[j] <=
    // system_time < time[j] + timestep[j].

    reset_stopping_conditions();    
    UpdatedParticles.clear();

    // Default calling sequence to smallN_evolve forces it to
    // integrate to the specified time, even if it thinks the
    // interaction is over.

    // int status = 
    smallN_evolve(b, to_time, break_scale_sq, structure_check_interval,
		  smalln_dtlog, smalln_verbose, smalln_outfile);
    //myprint(b);

    // Note: return status is independent of the termination reason...

    return 0;
}

int is_over(int * over, real rlimit, int verbose)
{
    real rlimit2 = rlimit*rlimit;
    if (rlimit2 <= 0) rlimit2 = _INFINITY_;
    *over = check_structure(b, rlimit2, verbose);

    return 0;
}

int synchronize_model()
{
    // Synchronize all particles at the current system time.
    // No action needed.

    return 0;
}

int get_time_step(double * time_step)
{
    *time_step = calculate_top_level_acc_and_jerk(b);
    return 0;
}

int get_index_of_first_particle(int * index_of_the_particle)
{
    *index_of_the_particle = -1;
    hdyn *od = b->get_oldest_daughter();
    if (!od) return -1;
    *index_of_the_particle = od->get_index();
    return 0;
}

int get_indices_of_colliding_particles(int * index_of_particle1, 
				       int * index_of_particle2)
{
    // TODO.

    *index_of_particle1 = -1;
    *index_of_particle2 = -1;
    return -1;
}

int get_number_of_particles(int * number_of_particles)
{
    *number_of_particles = 0;
    for_all_daughters(hdyn, b, bb) (*number_of_particles)++;
    return 0;
}

int get_total_mass(double * mass)
{
    *mass = 0;
    for_all_daughters(hdyn, b, bb) (*mass) += bb->get_mass();
    return 0;
}

int get_potential_energy(double * potential_energy)
{
    real kin, pot;
    get_energies(b, kin, pot);
    *potential_energy = pot;
    return 0;
}

int get_kinetic_energy(double * kinetic_energy)
{
    real kin, pot;
    get_energies(b, kin, pot);
    *kinetic_energy = kin;
    return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z)
{
    real mtot = 0;
    vec cmpos(0,0,0);
    for_all_daughters(hdyn, b, bb) {
	mtot += bb->get_mass();
	cmpos += bb->get_mass()*bb->get_pos();
    }
    *x = cmpos[0]/mtot;
    *y = cmpos[1]/mtot;
    *z = cmpos[2]/mtot;
    return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz)
{
    real mtot = 0;
    vec cmvel(0,0,0);
    for_all_daughters(hdyn, b, bb) {
	mtot += bb->get_mass();
	cmvel += bb->get_mass()*bb->get_vel();
    }
    *vx = cmvel[0]/mtot;
    *vy = cmvel[1]/mtot;
    *vz = cmvel[2]/mtot;
    return 0;
}

int get_total_radius(double * radius)
{
    vec cmpos;
    get_center_of_mass_position(&cmpos[0], &cmpos[1], &cmpos[2]);
    real r2max = 0;
    for_all_daughters(hdyn, b, bb) {
	real r2 = square(bb->get_pos()-cmpos);
	if (r2 > r2max) r2max = r2;
    }
    *radius = sqrt(r2max);
    return 0;
}

// Optional (not implemented; needed for bridge):

int get_potential_at_point(double eps,
			   double x, double y, double z, 
			   double * phi)
{
    return -1;
}

int get_gravity_at_point(double eps, double x, double y, double z, 
			 double * forcex, double * forcey, double * forcez)
{
    return -1;
}



int update_particle_tree(int over)	// default = 0
{
    b_copy = b;			// save the smallN tree until restored
    b = get_tree(b_copy);	// b now is the structured tree
				//	 	(why use b? - see below)

    if (over > 1) {

        // Over = 2 means a quasistable system has been detected.
        // Over = 3 means that the interaction is over (because the
        // system exceeded the size limit).  In either case,
        // get_tree2() may not have returned usable hierarchical
        // structure (e.g. an outer body escapes during an inner
        // triple resonance.  Impose usable structure by relaxing the
        // perturbation criterion on the components.  By effectively
        // ignoring all perturbations we ensure that top-level nodes
        // are unbound, leaving all bound substructure in a clump
        // which will be treatewd as a single object.

	rmtree(b);
	b = get_tree(b_copy, 1.e6);	// 1.e6 means ignore perturbations

	// Goals: (1) substructure, (2) top-level nodes are unbound.
	// Should check both here, but not clear what to do if they
	// aren't met...

        int nmul = 0;
	for_all_daughters(hdyn, b, bi)
	  if (bi->get_oldest_daughter()) nmul++;

        if (nmul == 0) {

	    // No substructure.  Flag for now.

	    cout << "*** SmallN: nmul = 0" << endl;
	}

	int nbound = 0;
	for_all_daughters(hdyn, b, bi) {
	    for_all_daughters(hdyn, b, bj) {
	        if (bj > bi) {
	            real m = bi->get_mass()+bj->get_mass();
		    vec dx = bi->get_pos()-bj->get_pos();
		    vec dv = bi->get_vel()-bj->get_vel();
		    if (0.5*dv*dv - m/sqrt(dx*dx) < 0) nbound++;
		}
	    }
	}

        if (nbound > 0) {

	    // Top level not unbound.  Flag for now.

	    cout << "*** SmallN: "; PRL(nbound);
	}
    }

    UpdatedParticles.clear();

    // Create indices for the CM nodes and add them to the updated list.

    if (b->get_cm_index() <= 0) {
	int max_index = -1;
	for_all_leaves(hdyn, b, bb)
	    if (bb->get_index() > max_index) max_index = bb->get_index();
	b->set_cm_index(pow(10, (int)log10((real)max_index+1) + 2.) + 1);
    }

    for_all_nodes(hdyn, b, bb)
	if (bb->get_parent() && bb->get_oldest_daughter()) {
	    int newindex = b->get_cm_index();
	    b->set_cm_index(newindex+1);
	    bb->set_index(newindex);
	    // cout << "push_back " << bb->get_index() << endl << flush;
	    UpdatedParticles.push_back(UpdatedParticle(newindex, 1));
	}

    // ARJEN: I'm pretty sure that the memory leak occurs because each
    // invocation to get_tree() above creates a new data structure,
    // but it is never deleted.  If I could do it locally I would
    // clean up as follows, but this causes an error, since
    // synchronize_to() in the driving script needs the structured
    // data.  It looks like restore_particle_tree() below was intended
    // to do just that, but it is not accessible from AMUSE and never
    // invoked.  Maybe it can be called through your state model?
    // -- Steve
    //
    // rmtree(b);
    // b = b_copy;

    return 0;
}

int restore_particle_tree()	// *** never used ***
{
    // Clean up temporary data structures.

    if (!b_copy) return -1;
    rmtree(b);
    b = b_copy;			// b is now the smallN tree again
    b_copy = NULL;
    UpdatedParticles.clear();
    return 0;
}

int get_number_of_particles_added(int * n_added)
{
    *n_added = (int)UpdatedParticles.size();
    return 0;
}

int get_id_of_added_particle(int index, int * index_of_particle)
{
    if (index < 0 || index > (int) UpdatedParticles.size())
        return -1;
    *index_of_particle = UpdatedParticles[index].index_of_particle;
    return 0;
}

int get_children_of_particle(int index_of_the_particle,
			     int * child1, int * child2)
{
    hdyn *bb = particle_with_index(b, index_of_the_particle);
    if (!bb) return -1;
    hdyn *od = bb->get_oldest_daughter();
    if (od && od->get_younger_sister()) {
	*child1 = od->get_index();
	*child2 = od->get_younger_sister()->get_index();
    } else {
	*child1 = *child2 = -1;
    }
    return 0;
}
