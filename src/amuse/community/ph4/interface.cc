#include "interface.h"

// A stub of this file is machine generated, but the content is
// hand-coded.  SAVE A COPY (here interface.cc.1) to avoid accidental
// overwriting!

#include "src/stdinc.h"
#include "src/jdata.h"
#include "src/idata.h"
#include "src/scheduler.h"

// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>

static jdata *jd = NULL;
static idata *id = NULL;
static scheduler *s = NULL;
static double begin_time = 0;

// Setup and parameters.

int initialize_code()
{
    // Begin the initialization by creating the jdata data structure.

    //PRL(1);
    jd = new jdata;

#ifndef NOMPI
    jd->setup_mpi(MPI::COMM_WORLD);
#endif
    //PRL(2);
    jd->setup_gpu();
    //PRL(3);
    if (jd->mpi_rank == 0) {
	cout << "initialize_code: ";
	PRC(jd->mpi_size); PRL(jd->have_gpu);
    }
    
    //PRL(4);
    begin_time = 0.0;
    jd->system_time = 0;		// ? TBD

    // AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(COLLISION_DETECTION);
    //PRL(5);
    mpi_setup_stopping_conditions();
    
    //PRL(6);
    jd->set_manage_encounters(4);	// 4 ==> enable AMUSE suport
    //PRL(7);

    return 0;
}

int set_eps2(double epsilon_squared)
{
    jd->eps2 = epsilon_squared;
    return 0;
}

int get_eps2(double * epsilon_squared)
{
    *epsilon_squared = jd->eps2;
    return 0;
}

int set_eta(double timestep_parameter)
{
    jd->eta = timestep_parameter;
    return 0;
}

int get_eta(double * timestep_parameter)
{
    *timestep_parameter = jd->eta;
    return 0;
}

int set_gpu(int gpu)
{
    jd->use_gpu = jd->have_gpu && (gpu == 1);
    return 0;
}

int get_gpu(int * gpu)
{
    *gpu = jd->use_gpu;
    return 0;
}

int set_gpu_id(int gpu_id)
{
    jd->gpu_id = gpu_id;
    return 0;
}

int get_gpu_id(int * gpu_id)
{
    *gpu_id = jd->gpu_id;
    return 0;
}

int set_manage_encounters(int m)
{
    jd->set_manage_encounters(m);
    return 0;
}

int get_manage_encounters(int * m)
{
    *m = jd->manage_encounters;
    return 0;
}

int set_time(double sys_time)
{
    jd->system_time = sys_time - begin_time;
    return 0;
}

int get_time(double * sys_time)
{
    *sys_time = jd->system_time + begin_time;
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


int commit_parameters()
{
    // Perform any needed setup after initial code parameters have been set.

    // Consistency check:

    if (jd->use_gpu && !jd->have_gpu) jd->use_gpu = false;

    if(jd->system_time == 0) {
        jd->system_time = 0;
    }
    
    if (jd->mpi_rank == 0) {
	cout << "commit_parameters: ";
	PRC(jd->have_gpu); PRL(jd->use_gpu);
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

    jd->initialize_arrays();
    id = new idata(jd);	  // set up idata data structures (sets acc and jerk)
    jd->set_initial_timestep();		// set timesteps (needs acc and jerk)
    s = new scheduler(jd);
#if 10
    cout << "commit_particles:";
    for (int j = 0; j < jd->nj; j++) cout << " " << jd->id[j];
    cout << endl << flush;
#endif
    return 0;
}

int recommit_particles()
{
    // Reinitialize/reset the system after particles have been added
    // or removed.  The system should be synchronized at some reasonable
    // system_time, so we just need to recompute forces and update the
    // GPU and scheduler.  Note that we don't resize the jdata or
    // idata arrays.  To resize idata, just delete and create a new
    // one.  Resizing jdata is more complicated -- defer for now.

    //cout << "recommit_particles" << endl << flush;
    if (!jd->use_gpu)
	jd->predict_all(jd->system_time, true);	// set pred quantities
    else
	jd->initialize_gpu(true);		// reload the GPU
    id->setup();				// compute acc and jerk
    jd->set_initial_timestep();			// set timesteps if not set
    s->initialize();				// reconstruct the scheduler
    //s->print(true);
    return 0;
}

int recompute_timesteps()
{
    // Same as recommit_particles(), except that the name isn't
    // reserved for the state model and we always recompute the time
    // steps.

    //cout << "recompute_timesteps" << endl << flush;
    if (!jd->use_gpu)
	jd->predict_all(jd->system_time, true);	// set pred quantities
    else
	jd->initialize_gpu(true);		// reload the GPU
    id->setup();				// compute acc and jerk
    jd->force_initial_timestep();
    s->initialize();				// reconstruct the scheduler
    return 0;
}

int cleanup_code()
{
    // Clean up at the end of the calculation.

    if (jd != NULL)
        jd->cleanup();
    if (id != NULL)
        id->cleanup();
    if (s != NULL)
        s->cleanup();
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

    *index_of_the_particle = jd->add_particle(mass, radius,
					      vec(x,y,z), vec(vx,vy,vz),
					      index_to_set);
#if 0
    int p = cout.precision(15);
    cout << "new_particle " << *index_of_the_particle << " " << mass
	    << " " << x << " " << y << " " << z   << endl;
    cout.precision(p);
#endif
    return 0;
}

int delete_particle(int index_of_the_particle)
{
    // Remove a particle from the system.

    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    jd->remove_particle(j);
    //cout << "delete_particle " << index_of_the_particle << endl;
    return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
			       int * index_of_the_next_particle)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    else if (j >= jd->nj-1) return 1;
    else *index_of_the_next_particle = jd->id[j+1];
    return 0;
}

int set_state(int index_of_the_particle,
	      double mass, 
	      double x, double y, double z,
	      double vx, double vy, double vz, double radius)
{
    //cout << "set_state" << endl << flush;
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    jd->mass[j] = mass;
    jd->radius[j] = radius;
    jd->pos[j][0] = x;
    jd->pos[j][1] = y;
    jd->pos[j][2] = z;
    jd->vel[j][0] = vx;
    jd->vel[j][1] = vy;
    jd->vel[j][2] = vz;
    return 0;
}

int get_state(int index_of_the_particle,
	      double * mass, 
	      double * x, double * y, double * z,
	      double * vx, double * vy, double * vz, double * radius)
{
    //cout << "get_state" << endl << flush;
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    *mass = jd->mass[j];
    *radius = jd->radius[j];
    *x = jd->pos[j][0];
    *y = jd->pos[j][1];
    *z = jd->pos[j][2];
    *vx = jd->vel[j][0];
    *vy = jd->vel[j][1];
    *vz = jd->vel[j][2];
    return 0;
}

int set_mass(int index_of_the_particle, double mass)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    jd->mass[j] = mass;
    return 0;
}

int get_mass(int index_of_the_particle, double * mass)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    *mass = jd->mass[j];
    return 0;
}

int set_radius(int index_of_the_particle, double radius)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    jd->radius[j] = radius;
    return 0;
}

int get_radius(int index_of_the_particle, double * radius)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    *radius = jd->radius[j];
    return 0;
}

int set_position(int index_of_the_particle,
		 double x, double y, double z)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    //cout << "set_position " << j<< endl << flush;
    if (j < 0) return -1;
    jd->pos[j][0] = x;
    jd->pos[j][1] = y;
    jd->pos[j][2] = z;
    return 0;
}

int get_position(int index_of_the_particle,
		 double * x, double * y, double * z)
{
    //cout << "get_position " << index_of_the_particle << endl << flush;
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    *x = jd->pos[j][0];
    *y = jd->pos[j][1];
    *z = jd->pos[j][2];
    return 0;
}

int set_velocity(int index_of_the_particle,
		 double vx, double vy, double vz)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    jd->vel[j][0] = vx;
    jd->vel[j][1] = vy;
    jd->vel[j][2] = vz;
    return 0;
}

int get_velocity(int index_of_the_particle,
		 double * vx, double * vy, double * vz)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    *vx = jd->vel[j][0];
    *vy = jd->vel[j][1];
    *vz = jd->vel[j][2];
    return 0;
}

int set_acceleration(int index_of_the_particle,
		     double ax, double ay, double az)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    jd->acc[j][0] = ax;
    jd->acc[j][1] = ay;
    jd->acc[j][2] = az;
    return 0;
}

int get_acceleration(int index_of_the_particle,
		     double * ax, double * ay, double * az)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    *ax = jd->acc[j][0];
    *ay = jd->acc[j][1];
    *az = jd->acc[j][2];
    return 0;
}

int get_potential(int index_of_the_particle, double * pot)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) return -1;
    *pot = jd->pot[j];
    return 0;
}

int get_particle_timestep(int index_of_the_particle, double * timestep)
{
    int j = jd->get_inverse_id(index_of_the_particle);
    if (j < 0) {
        return -1;
    }
    *timestep = jd->timestep[j];
    return 0;
}


// System-wide operations.

int evolve_model(double to_time)
{
    // On return, system_time will be greater than or equal to the
    // specified time.  All particles j will have time[j] <=
    // system_time < time[j] + timestep[j].  If synchronization is
    // needed, do it with synchronize_model().  The function breaks
    // out of the jd->advance() loop if an encounter is detected.

    if (jd->mpi_rank == 0) {
	//cout << "in evolve_model: "; PRC(to_time); PRL(jd->nj);
	for (int j = 0; j < jd->nj; j++)
	    if (jd->id[j] <= 0) {PRC(j); PRC(jd->mass[j]); PRL(jd->id[j]);}
    }

    reset_stopping_conditions();    
    jd->UpdatedParticles.clear();
    while (jd->system_time < (to_time - begin_time))
        if (jd->advance_and_check_encounter()) break;

    return 0;
}

int synchronize_model()
{
    // Synchronize all particles at the current system time.  The
    // default is not to reinitialize the scheduler, as this will be
    // handled later, in recommit_particles().

    //cout << "synchronize_model" << endl << flush;
    jd->UpdatedParticles.clear();
    jd->synchronize_all();
    return 0;
}

int get_time_step(double * time_step)
{
    *time_step = 0;			// not relevant here
    return -1;
}

int get_index_of_first_particle(int * index_of_the_particle)
{
    *index_of_the_particle = jd->id[0];
    return 0;
}

int get_indices_of_colliding_particles(int * index_of_particle1, 
				       int * index_of_particle2)
{
    *index_of_particle1 = jd->coll1;
    *index_of_particle2 = jd->coll2;
    return 0;
}

int get_number_of_particles(int * number_of_particles)
{
    *number_of_particles = jd->nj;
    return 0;
}

int get_total_mass(double * mass)
{
    *mass = jd->get_total_mass();
    return 0;
}

int get_potential_energy(double * potential_energy)
{
    *potential_energy = jd->get_pot();
    return 0;
}

int get_kinetic_energy(double * kinetic_energy)
{
    *kinetic_energy = jd->get_kin();
    return 0;
}

int get_binary_energy(double * binary_energy)	// new interface function
{						// see interface.py
    *binary_energy = jd->Emerge;		// not quite the same thing...
    return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z)
{
    // (Could also use jdata::get_com.)

    real mtot = 0;
    vec cmx(0,0,0);
    for (int j = 0; j < jd->nj; j++) {
	mtot += jd->mass[j];
	for (int k = 0; k < 3; k++) cmx[k] += jd->mass[j]*jd->pos[j][k];
    }
    *x = cmx[0]/mtot;
    *y = cmx[1]/mtot;
    *z = cmx[2]/mtot;
    return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz)
{
    // (Could also use jdata::get_com.)

    real mtot = 0;
    vec cmv(0,0,0);
    for (int j = 0; j < jd->nj; j++) {
	mtot += jd->mass[j];
	for (int k = 0; k < 3; k++) cmv[k] += jd->mass[j]*jd->vel[j][k];
    }
    *vx = cmv[0]/mtot;
    *vy = cmv[1]/mtot;
    *vz = cmv[2]/mtot;
    return 0;
}

int get_total_radius(double * radius)
{
    real r2max = 0;
    vec pos;
    get_center_of_mass_position(&pos[0], &pos[1], &pos[2]);
    for (int j = 0; j < jd->nj; j++)  {
	real r2 = 0;
	for (int k = 0; k < 3; k++) r2 += pow(jd->pos[j][k]-pos[k], 2);
	if (r2 > r2max) r2max = r2;
    }
    *radius = sqrt(r2max);
    return 0;
}

int get_potential_at_point(double * eps,
			   double * x, double * y, double * z, 
			   double * phi, int n)
{
    idata * tmp_idata = new idata();
    tmp_idata -> set_ni(n);
    tmp_idata -> ni = n; // set_ni only allocates buffers, does not set ni!
    tmp_idata -> jdat = jd;
    for(int i = 0; i < n; i++)
    {
        tmp_idata -> iid[i] = -1;
        tmp_idata -> ilist[i] = 0;
        tmp_idata -> inn[i] = 0;
        tmp_idata -> pnn[i] = 0;
        
        tmp_idata -> imass[i] = 1;
        tmp_idata -> iradius[i] = eps[i];
        tmp_idata -> itime[i] = 0.0;
        tmp_idata -> itimestep[i] = 0.0;
        tmp_idata -> ipot[i] = 0.0;
        tmp_idata -> ppot[i] = 0.0;
        tmp_idata -> idnn[i] = 0.0;
        tmp_idata -> ipos[i][0] = x[i];
        tmp_idata -> ipos[i][1] = y[i];
        tmp_idata -> ipos[i][2] = z[i];
        tmp_idata -> ivel[i][0] = 0.0;
        tmp_idata -> ivel[i][1] = 0.0;
        tmp_idata -> ivel[i][2] = 0.0;
    }
    tmp_idata -> get_acc_and_jerk();
    if (jd->mpi_rank == 0) {
        for(int i = 0; i < n; i++)
        {
            phi[i] = tmp_idata->ipot[i];
        }
    }
    delete tmp_idata;
    return 0;
}

int get_gravity_at_point(double * eps, double * x, double * y, double * z, 
			 double * forcex, double * forcey, double * forcez, int n)
{
    idata * tmp_idata = new idata();
    tmp_idata -> set_ni(n);
    tmp_idata -> ni = n; // set_ni only allocates buffers, does not set ni!
    tmp_idata -> jdat = jd;
    for(int i = 0; i < n; i++)
    {
        tmp_idata -> iid[i] = -1;
        tmp_idata -> ilist[i] = 0;
        tmp_idata -> inn[i] = 0;
        tmp_idata -> pnn[i] = 0;
        
        tmp_idata -> imass[i] = 1;
        tmp_idata -> iradius[i] = eps[i];
        tmp_idata -> itime[i] = 0.0;
        tmp_idata -> itimestep[i] = 0.0;
        tmp_idata -> ipot[i] = 0.0;
        tmp_idata -> ppot[i] = 0.0;
        tmp_idata -> idnn[i] = 0.0;
        tmp_idata -> ipos[i][0] = x[i];
        tmp_idata -> ipos[i][1] = y[i];
        tmp_idata -> ipos[i][2] = z[i];
        tmp_idata -> ivel[i][0] = 0.0;
        tmp_idata -> ivel[i][1] = 0.0;
        tmp_idata -> ivel[i][2] = 0.0;
    }
    tmp_idata -> get_acc_and_jerk();
    if (jd->mpi_rank == 0) {
        for(int i = 0; i < n; i++)
        {
            forcex[i] = tmp_idata->iacc[i][0];
            forcey[i] = tmp_idata->iacc[i][1];
            forcez[i] = tmp_idata->iacc[i][2];
        }
    }
    delete tmp_idata;
    return 0;
}

//----------------------------------------------------------------------
//
// From Arjen:

int get_number_of_particles_updated(int * value)
{
    *value = jd->UpdatedParticles.size();
    return 0;
}

int get_id_of_updated_particle(int index, int * index_of_particle, int * status)
{
    if (index < 0 || index > (int) jd->UpdatedParticles.size())
        return -1;
    *index_of_particle = jd->UpdatedParticles[index].index_of_particle;
    *status = jd->UpdatedParticles[index].status;
    return 0;
}
