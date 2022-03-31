#include "worker_code.h"

int get_mass(int index_of_the_particle, double * mass){
  return 0;
}

int commit_particles(){
  return 0;
}

int get_time(double * time){
  return 0;
}

int set_mass(int index_of_the_particle, double mass){
  return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
  return 0;
}

int get_total_radius(double * radius){
  return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz, double radius){
  return 0;
}

int get_total_mass(double * mass){
  return 0;
}

int evolve_model(double time){
  return 0;
}

int set_eps2(double epsilon_squared){
  return 0;
}

int get_begin_time(double * time){
  return 0;
}

int get_eps2(double * epsilon_squared){
  return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
  int * index_of_the_next_particle){
  return 0;
}

int delete_particle(int index_of_the_particle){
  return 0;
}

int get_potential(int index_of_the_particle, double * potential){
  return 0;
}

int synchronize_model(){
  return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y, 
  double z, double vx, double vy, double vz, double radius){
  return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz, 
  double * radius){
  return 0;
}

int get_time_step(double * time_step){
  return 0;
}

int recommit_particles(){
  return 0;
}

int get_kinetic_energy(double * kinetic_energy){
  return 0;
}

int get_number_of_particles(int * number_of_particles){
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
  return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  return 0;
}

int get_radius(int index_of_the_particle, double * radius){
  return 0;
}

int set_begin_time(double time){
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  return 0;
}

int cleanup_code(){
  return 0;
}

int recommit_parameters(){
  return 0;
}

int initialize_code(){
  return 0;
}

int get_potential_energy(double * potential_energy){
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
  double * vz){
  return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
  double * z){
  return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
  double * az){
  return 0;
}

int commit_parameters(){
  return 0;
}

int set_parameters(char * param_file){
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
  double vz){
  return 0;
}

