#pragma once
#include <stopcond.h>

#ifdef __cplusplus
extern "C" {
#endif

//  GravitationalDynamicsCode

int initialize_code();

int cleanup_code();

int commit_parameters();

int recommit_parameters();

int new_particle(int* index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double radius);

int delete_particle(int index_of_the_particle);

int get_state(int index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, double * radius);

int set_state(int index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double radius);

int get_mass(int index_of_the_particle, double * mass);

int set_mass(int index_of_the_particle, double mass);

int get_radius(int index_of_the_particle, double * radius);

int set_radius(int index_of_the_particle, double radius);

int get_position(int index_of_the_particle, double * x, double * y, double * z);

int set_position(int index_of_the_particle, double x, double y, double z);

int get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz);

int set_velocity(int index_of_the_particle, double vx, double vy, double vz);

int get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az);

int set_acceleration(int index_of_the_particle, double ax, double ay, double az);

int get_potential(int index_of_the_particle, double * potential);

int evolve_model(double time);

int commit_particles();

int synchronize_model();

int recommit_particles();

int reconstruct_particle_list();

int get_eps2(double * epsilon_squared);

int set_eps2(double epsilon_squared);

int get_kinetic_energy(double * kinetic_energy);

int get_potential_energy(double * potential_energy);

int get_time(double * time);

int get_begin_time(double * time);

int set_begin_time(double time);

int get_time_step(double * time_step);

int get_total_mass(double * mass);

int get_center_of_mass_position(double * x, double * y, double * z);

int get_center_of_mass_velocity(double * vx, double * vy, double * vz);

int get_total_radius(double * radius);

int get_number_of_particles(int * number_of_particles);

int get_index_of_first_particle(int * index_of_the_particle);

int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle);

// GravityFieldCode

int get_gravity_at_point(double * eps, double * x, double * y, double * z, double * ax, double * ay, double * az, int npoints);

int get_potential_at_point(double * eps, double * x, double * y, double * z, double * phi, int npoints);

//int get_eta(double * eta);

//int set_eta(double eta);

#ifdef __cplusplus
}
#endif

