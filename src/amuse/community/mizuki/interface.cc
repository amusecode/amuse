#include "worker_code.h"
#include "stopcond.h"
#include "src/mizuki/user_defined.hpp"
#include "src/mizuki/mizuki.hpp"

static Mizuki* mizuki = new Mizuki;

int initialize_code(){
    mizuki->initialize();
    return 0;
}

int commit_particles(){
    mizuki->commit_particles();
    return 0;
}

int new_sph_particle(int * index_of_the_particle, double mass, double x,
    double y, double z, double vx, double vy, double vz, double eng, double h){
    *index_of_the_particle = mizuki->add_sph_particle(
            mass, x, y, z, vx, vy, vz, eng, h*2
            );
    return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, 
    double y, double z, double vx, double vy, double vz, double radius){
    return 0;
}

int get_mass(int index_of_the_particle, double * mass){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *mass = p->mass;
    return 0;
}

int set_mass(int index_of_the_particle, double mass){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    p->mass = mass;
    return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
    double * z){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *x = p->pos.x;
    *y = p->pos.y;
    *z = p->pos.z;
    return 0;
}

int set_position(int index_of_the_particle, double x, double y, 
    double z){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    p->pos.x = x;
    p->pos.y = y;
    p->pos.z = z;
    return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
    double * vz){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *vx = p->vel.x;
    *vy = p->vel.y;
    *vz = p->vel.z;
    return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
    double vz){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    p->vel.x = vx;
    p->vel.y = vy;
    p->vel.z = vz;
    return 0;
}

int get_internal_energy(int index_of_the_particle, double * eng){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *eng = p->eng;
    return 0;
}

int set_internal_energy(int index_of_the_particle, double eng){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    p->eng = eng;
    return 0;
}

int get_smoothing_length(int index_of_the_particle, double * smth){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *smth = p->smth/2;
    return 0;
}

int set_smoothing_length(int index_of_the_particle, double smth){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    p->smth = smth*2;
    return 0;
}

int get_radius(int index_of_the_particle, double * radius){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *radius = p->smth/2;
    return 0;
}

int set_radius(int index_of_the_particle, double radius){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    p->smth = radius*2;
    return 0;
}

int get_density(int index_of_the_particle, double * dens){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *dens = p->dens;
    return 0;
}

int get_pressure(int index_of_the_particle, double * pres){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *pres = p->pres;
    return 0;
}

int get_state_sph(int index_of_the_particle, double * mass, double * x, 
    double * y, double * z, double * vx, double * vy, double * vz, 
    double * eng, double * smth){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *mass = p->mass;
    *x = p->pos.x;
    *y = p->pos.y;
    *z = p->pos.z;
    *vx = p->vel.x;
    *vy = p->vel.y;
    *vz = p->vel.z;
    *eng = p->eng;
    *smth = p->smth/2;
    return 0;
}

int set_state_sph(int index_of_the_particle, double mass, double x, 
    double y, double z, double vx, double vy, double vz, 
    double eng, double smth){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    p->mass = mass;
    p->pos.x = x;
    p->pos.y = y;
    p->pos.z = z;
    p->vel.x = vx;
    p->vel.y = vy;
    p->vel.z = vz;
    p->eng = eng;
    p->smth = smth*2;
    return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
    double * y, double * z, double * vx, double * vy, double * vz, 
    double * radius){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *mass = p->mass;
    *x = p->pos.x;
    *y = p->pos.y;
    *z = p->pos.z;
    *vx = p->vel.x;
    *vy = p->vel.y;
    *vz = p->vel.z;
    *radius = 0;
    return 0;
}

int set_stopping_condition_timeout_parameter(double value){
    return 0;
}

int get_time(double * time){
    *time = mizuki->time;
    return 0;
}

int get_stopping_condition_maximum_density_parameter(double * value){
    return 0;
}

int set_stopping_condition_out_of_box_use_center_of_mass_parameter(
    _Bool value){
    return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
    return 0;
}

int get_stopping_condition_out_of_box_use_center_of_mass_parameter(
    _Bool * value){
    return 0;
}

int enable_stopping_condition(int type){
    return 0;
}

int get_total_radius(double * radius){
    return 0;
}

int get_stopping_condition_maximum_internal_energy_parameter(
    double * value){
    return 0;
}

int get_potential_at_point(double eps, double x, double y, double z, 
    double * phi, int npoints){
    return 0;
}

int get_number_of_stopping_conditions_set(int * result){
    return 0;
}

int is_stopping_condition_set(int type, int * result){
    return 0;
}

int get_total_mass(double * mass){
    return 0;
}

int evolve_model(double time){
    mizuki->evolve_model(time);
    return 0;
}

int set_stopping_condition_out_of_box_parameter(double value){
    return 0;
}

int set_eps2(double epsilon_squared){
    return 0;
}

int set_stopping_condition_number_of_steps_parameter(int value){
    return 0;
}

int get_stopping_condition_timeout_parameter(double * value){
    return 0;
}

int get_begin_time(double * time){
    return 0;
}

int get_eps2(double * epsilon_squared){
    *epsilon_squared = mizuki->epsilon_gravity * mizuki->epsilon_gravity;
    return 0;
}

int get_stopping_condition_minimum_internal_energy_parameter(
    double * value){
    return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
    int * index_of_the_next_particle){
    return 0;
}

int delete_particle(int index_of_the_particle){
    return 0;
}

int is_stopping_condition_enabled(int type, int * result){
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

int get_stopping_condition_minimum_density_parameter(double * value){
    return 0;
}

int get_time_step(double * time_step){
    *time_step = mizuki->dt_max;
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

int get_stopping_condition_number_of_steps_parameter(int * value){
    return 0;
}

int disable_stopping_condition(int type){
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

int set_stopping_condition_minimum_internal_energy_parameter(
    double value){
    return 0;
}

int set_begin_time(double time){
    return 0;
}

int set_stopping_condition_minimum_density_parameter(double value){
    return 0;
}

int has_stopping_condition(int type, int * result){
    return 0;
}

int cleanup_code(){
    return 0;
}

int set_stopping_condition_maximum_density_parameter(double value){
    return 0;
}

int recommit_parameters(){
    return 0;
}

int get_potential_energy(double * potential_energy){
    return 0;
}

int get_gravity_at_point(double eps, double x, double y, double z, 
    double * ax, double * ay, double * az, int npoints){
    return 0;
}

int get_stopping_condition_out_of_box_parameter(double * value){
    return 0;
}

int set_stopping_condition_maximum_internal_energy_parameter(
    double value){
    return 0;
}

int get_stopping_condition_info(int index, int * type, 
    int * number_of_particles){
    return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
    double * az){
    FP_sph* p = &(mizuki->psys_sph[index_of_the_particle]);
    *ax = p->acc_hydro.x + p->acc_grav.x;
    *ay = p->acc_hydro.y + p->acc_grav.y;
    *az = p->acc_hydro.z + p->acc_grav.z;
    return 0;
}

int commit_parameters(){
    return 0;
}

int get_stopping_condition_particle_index(int index, 
    int index_of_the_column, int * index_of_particle){
    return 0;
}
