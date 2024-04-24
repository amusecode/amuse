#include "src/Hermite_GRX/src/Model.h"
#include "src/Hermite_GRX/src/Perturbation.h"
#include "src/Hermite_GRX/src/Integrator.h"

#include "src/Hermite_GRX/src/HermiteIntegrator.h"
#include "src/Hermite_GRX/src/RegularizedHermiteIntegrator.h"

#include "src/Hermite_GRX/src/Pairwise1PNPerturbation.h"
#include "src/Hermite_GRX/src/EIHPerturbation.h"

#include <iostream>
#include <map>
#include <string>

using namespace std;

Model *model = nullptr;
Integrator *integrator = nullptr;
Perturbation *perturbation = nullptr;

map<string, Integrator *> integrators;
map<string, Perturbation *> perturbations;

string current_integrator = "Hermite";
string current_perturbation = "None";

Real time_step_parameter = 0.03;
Real speed_of_light = 1.0;
size_t num_threads = 1;
Real epsilon2 = 0;

int id_counter = 0;

int initialize_code()
{
	if (model)
		return 0;

	model = new Model();

	integrators["Hermite"] = new HermiteIntegrator(false);
	integrators["SymmetrizedHermite"] = new HermiteIntegrator(true);
	integrators["RegularizedHermite"] = new RegularizedHermiteIntegrator(false);
	integrators["SymmetrizedRegularizedHermite"] = new RegularizedHermiteIntegrator(true);

	perturbations["None"] = new Perturbation();
	perturbations["1PN_Pairwise"] = new Pairwise1PNPerturbation();
	perturbations["1PN_EIH"] = new EIH1PNPerturbation();
	perturbations["2.5PN_EIH"] = new EIH2dot5PNPerturbation();
	

	return 0;
}

int commit_parameters()
{
	try
	{
		integrator = integrators.at(current_integrator);
		perturbation = perturbations.at(current_perturbation);

		integrator->SetPerturbation(perturbation);
		integrator->SetNumThreads(num_threads);
		integrator->SetTimeStepParameter(time_step_parameter);
		perturbation->SetSpeedOfLight(speed_of_light);
	}
	catch (...)
	{
		integrator = nullptr;
		perturbation = nullptr;
		return -1;
	}

	return 0;
}

int recommit_parameters()
{
	return commit_parameters();
}

int commit_particles()
{
	return 0;
}

int recommit_particles()
{
	return commit_particles();
}

int get_eps2(double *eps2)
{
	*eps2 = epsilon2;
	return 0;
}

int set_eps2(double eps2)
{
	epsilon2 = eps2;
	return 0;
}

int get_dt_param(double *dt_param)
{
	*dt_param = time_step_parameter;
	return 0;
}

int set_dt_param(double dt_param)
{
	time_step_parameter = dt_param;
	return 0;
}

int get_time(double *t)
{
	*t = model->GetTime();
	return 0;
}

int set_time(double t)
{
	model->SetTime(t);
	return 0;
}

int get_light_speed(double *c)
{
	*c = speed_of_light;
	return 0;
}

int set_light_speed(double c)
{
	speed_of_light = c;
	return 0;
}

int get_integrator(char **i)
{
	// Casting away constness
	*i = (char *)current_integrator.c_str();
	return 0;
}

int set_integrator(char *i)
{
	current_integrator = string(i);
	return 0;
}

int get_perturbation(char **i)
{
	// Casting away constness
	*i = (char *)current_perturbation.c_str();
	return 0;
}

int set_perturbation(char *i)
{
	current_perturbation = string(i);
	return 0;
}

int get_num_threads(int *thread_number)
{
	*thread_number = (int)num_threads;
	return 0;
}

int set_num_threads(int thread_number)
{
	if (thread_number <= 0)
		return -1;

	num_threads = thread_number;
	return 0;
}

int cleanup_code()
{
	if (model)
		delete model;
	model = nullptr;

	for (auto &p : perturbations)
		delete p.second;
	perturbations.clear();
	perturbation = nullptr;

	for (auto &i : integrators)
		delete i.second;
	integrators.clear();
	integrator = nullptr;

	return 0;
}

int new_large_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
	*id = id_counter++;
	Particle p = Particle(*id, mass, Vec(x, y, z), Vec(vx, vy, vz), radius);
	model->AddLargeParticle(p);
	return 0;
}

int new_small_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
	*id = id_counter++;
	Particle p = Particle(*id, mass, Vec(x, y, z), Vec(vx, vy, vz), radius);
	model->AddSmallParticle(p);
	return 0;
}

int new_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
	return new_small_particle(id, mass, x, y, z, vx, vy, vz, radius);
}

int delete_particle(int id)
{
	model->RemoveParticle(id);
	return 0;
}

int get_state(int id, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *radius)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		*mass = p->mass;
		*x = p->pos[0];
		*y = p->pos[1];
		*z = p->pos[2];
		*vx = p->vel[0];
		*vy = p->vel[1];
		*vz = p->vel[2];
		*radius = p->radius;
		return 0;
	}
	else
	{
		return -1;
	}
}

int set_state(int id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		p->mass = mass;
		p->pos = Vec(x, y, z);
		p->vel = Vec(vx, vy, vz);
		p->radius = radius;
		return 0;
	}
	else
	{
		return -1;
	}
}

int get_mass(int id, double *mass)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		*mass = p->mass;
		return 0;
	}
	else
	{
		return -1;
	}
}

int set_mass(int id, double mass)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		p->mass = mass;
		return 0;
	}
	else
	{
		return -1;
	}
}

int get_radius(int id, double *radius)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		*radius = p->radius;
		return 0;
	}
	else
	{
		return -1;
	}
}

int set_radius(int id, double radius)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		p->radius = radius;
		return 0;
	}
	else
	{
		return -1;
	}
}

int get_position(int id, double *x, double *y, double *z)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		*x = p->pos[0];
		*y = p->pos[1];
		*z = p->pos[2];
		return 0;
	}
	else
	{
		return -1;
	}
}

int set_position(int id, double x, double y, double z)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		p->pos = Vec(x, y, z);
		return 0;
	}
	else
	{
		return -1;
	}
}

int get_velocity(int id, double *vx, double *vy, double *vz)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		*vx = p->vel[0];
		*vy = p->vel[1];
		*vz = p->vel[2];
		return 0;
	}
	else
	{
		return -1;
	}
}

int set_velocity(int id, double vx, double vy, double vz)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		p->vel = Vec(vx, vy, vz);
		return 0;
	}
	else
	{
		return -1;
	}
}

int get_acceleration(int id, double *ax, double *ay, double *az)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		Vec acc = p->acc_newton + p->acc_pert;
		*ax = acc[0];
		*ay = acc[1];
		*az = acc[2];
		return 0;
	}
	else
	{
		return -1;
	}
}

int set_acceleration(int id, double ax, double ay, double az)
{
	return -2;
}

int get_jerk(int id, double *sx, double *sy, double *sz)
{
	Particle *p = model->GetParticle(id);
	if (p)
	{
		Vec jerk = p->jerk_newton + p->jerk_pert;
		*sx = jerk[0];
		*sy = jerk[1];
		*sz = jerk[2];
		return 0;
	}
	else
	{
		return -1;
	}
}

int set_jerk(int id, double jx, double jy, double jz)
{
	return -2;
}

int evolve_model(double t)
{
	if (num_threads > model->GetNumParticles())
		return -1;

	integrator->Evolve(model, t);
	return 0;
}

int get_kinetic_energy(double *ekin)
{
	*ekin = model->GetNewtonianKineticEnergy();
	return 0;
}

int get_potential_energy(double *epot)
{
	*epot = model->GetNewtonianPotentialEnergy();
	return 0;
}

int get_total_energy_with(char *type, double *etot)
{
	try
	{
		Perturbation *pert = perturbations.at(string(type));
		pert->SetSpeedOfLight(speed_of_light);

		*etot = model->GetNewtonianEnergy();
		*etot += pert->GetEnergy(model);
	}
	catch (...)
	{
		return -1;
	}

	return 0;
}

int get_total_linear_momentum_with(char *type, double *px, double *py, double *pz)
{
	try
	{
		Perturbation *pert = perturbations.at(string(type));
		pert->SetSpeedOfLight(speed_of_light);

		Vec p = model->GetNewtonianLinearMomentum();
		p += pert->GetLinearMomentum(model);

		*px = p[0];
		*py = p[1];
		*pz = p[2];
	}
	catch (...)
	{
		return -1;
	}

	return 0;
}

int get_time_step(double *dt)
{
	return -2;
}

int get_potential(int id, double *val)
{
	return -2;
}

int get_potential_at_point(double eps, double x, double y, double z, double *phi)
{
	return -2;
}

int get_gravity_at_point(double eps, double x, double y, double z, double *ax, double *ay, double *az)
{
	return -2;
}

int get_total_mass(double *mass)
{
	ParticleSetView &all = model->GetAllParticles();
	Real m = 0;

	for (Particle &p : all)
		m += p.mass;

	*mass = m;
	return 0;
}

int get_center_of_mass_position(double *x, double *y, double *z)
{
	return -2;
}

int get_center_of_mass_velocity(double *vx, double *vy, double *vz)
{
	return -2;
}

int get_total_radius(double *radius)
{
	return -2;
}

int get_number_of_particles(int *num_parts)
{
	*num_parts = model->GetNumLargeParticles() + model->GetNumSmallParticles();
	return 0;
}

int get_index_of_first_particle(int *index)
{
	return -2;
}

int get_index_of_next_particle(int id, int *index)
{
	return -2;
}

int synchronize_model()
{
	return 0;
}

int get_begin_time(double *t)
{
	*t = model->GetTime();
	return 0;
}

int set_begin_time(double t)
{
	model->SetTime(t);
	return 0;
}
