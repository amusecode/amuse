import numpy
cimport numpy
cdef extern from "worker_code.h":
    pass
    
cdef extern from "stopcond.h":
    pass
cdef extern from "mpi.h":
    pass
cdef extern from "amuse_mpi.h":
    pass
    
cimport mpi4py.MPI 

cdef extern from "worker_code.h":
  
  int c_set_comm_world "set_comm_world" (mpi4py.MPI.MPI_Comm world);
  
  int c_get_mass "get_mass" (int, double *);
  
  
  int c_commit_particles "commit_particles" ();
  
  
  int c_set_stopping_condition_timeout_parameter "set_stopping_condition_timeout_parameter" (double);
  
  
  int c_get_time "get_time" (double *);
  
  
  int c_set_mass "set_mass" (int, double);
  
  
  int c_get_stopping_condition_maximum_density_parameter "get_stopping_condition_maximum_density_parameter" (double *);
  
  
  int c_set_stopping_condition_out_of_box_use_center_of_mass_parameter "set_stopping_condition_out_of_box_use_center_of_mass_parameter" (int);
  
  
  int c_get_index_of_first_particle "get_index_of_first_particle" (int *);
  
  
  int c_get_stopping_condition_out_of_box_use_center_of_mass_parameter "get_stopping_condition_out_of_box_use_center_of_mass_parameter" (int *);
  
  
  int c_get_dt_dia "get_dt_dia" (double *);
  
  
  int c_enable_stopping_condition "enable_stopping_condition" (int);
  
  
  int c_set_end_time_accuracy_factor "set_end_time_accuracy_factor" (double);
  
  
  int c_get_total_radius "get_total_radius" (double *);
  
  
  int c_get_stopping_condition_maximum_internal_energy_parameter "get_stopping_condition_maximum_internal_energy_parameter" (double *);
  
  
  int c_get_potential_at_point "get_potential_at_point" (double, double, double, double, double *);
  
  
  int c_get_number_of_stopping_conditions_set "get_number_of_stopping_conditions_set" (int *);
  
  
  int c_is_stopping_condition_set "is_stopping_condition_set" (int, int *);
  
  
  int c_new_particle "new_particle" (int *, double, double, double, double, double, double, double, double);
  
  
  int c_get_end_time_accuracy_factor "get_end_time_accuracy_factor" (double *);
  
  
  int c_get_total_mass "get_total_mass" (double *);
  
  
  int c_evolve_model "evolve_model" (double);
  
  
  int c_set_stopping_condition_out_of_box_parameter "set_stopping_condition_out_of_box_parameter" (double);
  
  
  int c_set_eps2 "set_eps2" (double);
  
  
  int c_set_stopping_condition_number_of_steps_parameter "set_stopping_condition_number_of_steps_parameter" (int);
  
  
  int c_get_stopping_condition_timeout_parameter "get_stopping_condition_timeout_parameter" (double *);
  
  
  int c_get_begin_time "get_begin_time" (double *);
  
  
  int c_get_eps2 "get_eps2" (double *);
  
  
  int c_set_dt_dia "set_dt_dia" (double);
  
  
  int c_get_stopping_condition_minimum_internal_energy_parameter "get_stopping_condition_minimum_internal_energy_parameter" (double *);
  
  
  int c_get_index_of_next_particle "get_index_of_next_particle" (int, int *);
  
  
  int c_delete_particle "delete_particle" (int);
  
  
  int c_is_stopping_condition_enabled "is_stopping_condition_enabled" (int, int *);
  
  
  int c_get_is_time_reversed_allowed "get_is_time_reversed_allowed" (int *);
  
  
  int c_set_is_time_reversed_allowed "set_is_time_reversed_allowed" (int);
  
  
  int c_get_potential "get_potential" (int, double *);
  
  
  int c_synchronize_model "synchronize_model" ();
  
  
  int c_set_state "set_state" (int, double, double, double, double, double, double, double, double);
  
  
  int c_get_stopping_condition_minimum_density_parameter "get_stopping_condition_minimum_density_parameter" (double *);
  
  
  int c_get_state "get_state" (int, double *, double *, double *, double *, double *, double *, double *, double *);
  
  
  int c_get_time_step "get_time_step" (double *);
  
  
  int c_recommit_particles "recommit_particles" ();
  
  
  int c_get_kinetic_energy "get_kinetic_energy" (double *);
  
  
  int c_get_number_of_particles "get_number_of_particles" (int *);
  
  
  int c_get_stopping_condition_number_of_steps_parameter "get_stopping_condition_number_of_steps_parameter" (int *);
  
  
  int c_disable_stopping_condition "disable_stopping_condition" (int);
  
  
  int c_set_acceleration "set_acceleration" (int, double, double, double);
  
  
  int c_get_center_of_mass_position "get_center_of_mass_position" (double *, double *, double *);
  
  
  int c_get_center_of_mass_velocity "get_center_of_mass_velocity" (double *, double *, double *);
  
  
  int c_get_radius "get_radius" (int, double *);
  
  
  int c_set_stopping_condition_minimum_internal_energy_parameter "set_stopping_condition_minimum_internal_energy_parameter" (double);
  
  
  int c_set_begin_time "set_begin_time" (double);
  
  
  int c_set_stopping_condition_minimum_density_parameter "set_stopping_condition_minimum_density_parameter" (double);
  
  
  int c_set_radius "set_radius" (int, double);
  
  
  int c_has_stopping_condition "has_stopping_condition" (int, int *);
  
  
  int c_cleanup_code "cleanup_code" ();
  
  
  int c_set_stopping_condition_maximum_density_parameter "set_stopping_condition_maximum_density_parameter" (double);
  
  
  int c_recommit_parameters "recommit_parameters" ();
  
  
  int c_initialize_code "initialize_code" ();
  
  
  int c_get_potential_energy "get_potential_energy" (double *);
  
  
  int c_get_gravity_at_point "get_gravity_at_point" (double, double, double, double, double *, double *, double *);
  
  
  int c_get_velocity "get_velocity" (int, double *, double *, double *);
  
  
  int c_get_stopping_condition_out_of_box_parameter "get_stopping_condition_out_of_box_parameter" (double *);
  
  
  int c_set_dt_param "set_dt_param" (double);
  
  
  int c_get_position "get_position" (int, double *, double *, double *);
  
  
  int c_set_stopping_condition_maximum_internal_energy_parameter "set_stopping_condition_maximum_internal_energy_parameter" (double);
  
  
  int c_get_dt_param "get_dt_param" (double *);
  
  
  int c_set_position "set_position" (int, double, double, double);
  
  
  int c_get_stopping_condition_info "get_stopping_condition_info" (int, int *, int *);
  
  
  int c_get_acceleration "get_acceleration" (int, double *, double *, double *);
  
  
  int c_commit_parameters "commit_parameters" ();
  
  
  int c_get_stopping_condition_particle_index "get_stopping_condition_particle_index" (int, int, int *);
  
  
  int c_set_velocity "set_velocity" (int, double, double, double);
  

def set_comm_world(mpi4py.MPI.Comm comm not None):
    return c_set_comm_world(comm.ob_mpi)
    
def get_mass(index_of_the_particle, mass):
  
  cdef double output_mass
  cdef int __result__ = c_get_mass(index_of_the_particle, &output_mass);
  mass.value = output_mass
  return __result__


def commit_particles():
  
  cdef int __result__ = c_commit_particles();
  return __result__


def set_stopping_condition_timeout_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_timeout_parameter(value);
  return __result__


def get_time(time):
  
  cdef double output_time
  cdef int __result__ = c_get_time(&output_time);
  time.value = output_time
  return __result__


def set_mass(index_of_the_particle, mass):
  
  cdef int __result__ = c_set_mass(index_of_the_particle, mass);
  return __result__


def get_stopping_condition_maximum_density_parameter(value):
  
  cdef double output_value
  cdef int __result__ = c_get_stopping_condition_maximum_density_parameter(&output_value);
  value.value = output_value
  return __result__


def set_stopping_condition_out_of_box_use_center_of_mass_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_out_of_box_use_center_of_mass_parameter(value);
  return __result__


def get_index_of_first_particle(index_of_the_particle):
  
  cdef int output_index_of_the_particle
  cdef int __result__ = c_get_index_of_first_particle(&output_index_of_the_particle);
  index_of_the_particle.value = output_index_of_the_particle
  return __result__


def get_stopping_condition_out_of_box_use_center_of_mass_parameter(value):
  
  cdef int output_value
  cdef int __result__ = c_get_stopping_condition_out_of_box_use_center_of_mass_parameter(&output_value);
  value.value = output_value
  return __result__


def get_dt_dia(dt_dia):
  
  cdef double output_dt_dia
  cdef int __result__ = c_get_dt_dia(&output_dt_dia);
  dt_dia.value = output_dt_dia
  return __result__


def enable_stopping_condition(type):
  
  cdef int __result__ = c_enable_stopping_condition(type);
  return __result__


def set_end_time_accuracy_factor(value):
  
  cdef int __result__ = c_set_end_time_accuracy_factor(value);
  return __result__


def get_total_radius(radius):
  
  cdef double output_radius
  cdef int __result__ = c_get_total_radius(&output_radius);
  radius.value = output_radius
  return __result__


def get_stopping_condition_maximum_internal_energy_parameter(value):
  
  cdef double output_value
  cdef int __result__ = c_get_stopping_condition_maximum_internal_energy_parameter(&output_value);
  value.value = output_value
  return __result__


def get_potential_at_point(eps, x, y, z, phi):
  
  cdef double output_phi
  cdef int __result__ = c_get_potential_at_point(eps, x, y, z, &output_phi);
  phi.value = output_phi
  return __result__


def get_number_of_stopping_conditions_set(result):
  
  cdef int output_result
  cdef int __result__ = c_get_number_of_stopping_conditions_set(&output_result);
  result.value = output_result
  return __result__


def is_stopping_condition_set(type, result):
  
  cdef int output_result
  cdef int __result__ = c_is_stopping_condition_set(type, &output_result);
  result.value = output_result
  return __result__


def new_particle(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius):
  
  cdef int output_index_of_the_particle
  cdef int __result__ = c_new_particle(&output_index_of_the_particle, mass, x, y, z, vx, vy, vz, radius);
  index_of_the_particle.value = output_index_of_the_particle
  return __result__


def get_end_time_accuracy_factor(value):
  
  cdef double output_value
  cdef int __result__ = c_get_end_time_accuracy_factor(&output_value);
  value.value = output_value
  return __result__


def get_total_mass(mass):
  
  cdef double output_mass
  cdef int __result__ = c_get_total_mass(&output_mass);
  mass.value = output_mass
  return __result__


def evolve_model(time):
  
  cdef int __result__ = c_evolve_model(time);
  return __result__


def set_stopping_condition_out_of_box_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_out_of_box_parameter(value);
  return __result__


def set_eps2(epsilon_squared):
  
  cdef int __result__ = c_set_eps2(epsilon_squared);
  return __result__


def set_stopping_condition_number_of_steps_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_number_of_steps_parameter(value);
  return __result__


def get_stopping_condition_timeout_parameter(value):
  
  cdef double output_value
  cdef int __result__ = c_get_stopping_condition_timeout_parameter(&output_value);
  value.value = output_value
  return __result__


def get_begin_time(time):
  
  cdef double output_time
  cdef int __result__ = c_get_begin_time(&output_time);
  time.value = output_time
  return __result__


def get_eps2(epsilon_squared):
  
  cdef double output_epsilon_squared
  cdef int __result__ = c_get_eps2(&output_epsilon_squared);
  epsilon_squared.value = output_epsilon_squared
  return __result__


def set_dt_dia(dt_dia):
  
  cdef int __result__ = c_set_dt_dia(dt_dia);
  return __result__


def get_stopping_condition_minimum_internal_energy_parameter(value):
  
  cdef double output_value
  cdef int __result__ = c_get_stopping_condition_minimum_internal_energy_parameter(&output_value);
  value.value = output_value
  return __result__


def get_index_of_next_particle(index_of_the_particle, index_of_the_next_particle):
  
  cdef int output_index_of_the_next_particle
  cdef int __result__ = c_get_index_of_next_particle(index_of_the_particle, &output_index_of_the_next_particle);
  index_of_the_next_particle.value = output_index_of_the_next_particle
  return __result__


def delete_particle(index_of_the_particle):
  
  cdef int __result__ = c_delete_particle(index_of_the_particle);
  return __result__


def is_stopping_condition_enabled(type, result):
  
  cdef int output_result
  cdef int __result__ = c_is_stopping_condition_enabled(type, &output_result);
  result.value = output_result
  return __result__


def get_is_time_reversed_allowed(value):
  
  cdef int output_value
  cdef int __result__ = c_get_is_time_reversed_allowed(&output_value);
  value.value = output_value
  return __result__


def set_is_time_reversed_allowed(value):
  
  cdef int __result__ = c_set_is_time_reversed_allowed(value);
  return __result__


def get_potential(index_of_the_particle, potential):
  
  cdef double output_potential
  cdef int __result__ = c_get_potential(index_of_the_particle, &output_potential);
  potential.value = output_potential
  return __result__


def synchronize_model():
  
  cdef int __result__ = c_synchronize_model();
  return __result__


def set_state(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius):
  
  cdef int __result__ = c_set_state(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius);
  return __result__


def get_stopping_condition_minimum_density_parameter(value):
  
  cdef double output_value
  cdef int __result__ = c_get_stopping_condition_minimum_density_parameter(&output_value);
  value.value = output_value
  return __result__


def get_state(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius):
  
  cdef double output_mass
  cdef double output_x
  cdef double output_y
  cdef double output_z
  cdef double output_vx
  cdef double output_vy
  cdef double output_vz
  cdef double output_radius
  cdef int __result__ = c_get_state(index_of_the_particle, &output_mass, &output_x, &output_y, &output_z, &output_vx, &output_vy, &output_vz, &output_radius);
  mass.value = output_mass
  x.value = output_x
  y.value = output_y
  z.value = output_z
  vx.value = output_vx
  vy.value = output_vy
  vz.value = output_vz
  radius.value = output_radius
  return __result__


def get_time_step(time_step):
  
  cdef double output_time_step
  cdef int __result__ = c_get_time_step(&output_time_step);
  time_step.value = output_time_step
  return __result__


def recommit_particles():
  
  cdef int __result__ = c_recommit_particles();
  return __result__


def get_kinetic_energy(kinetic_energy):
  
  cdef double output_kinetic_energy
  cdef int __result__ = c_get_kinetic_energy(&output_kinetic_energy);
  kinetic_energy.value = output_kinetic_energy
  return __result__


def get_number_of_particles(number_of_particles):
  
  cdef int output_number_of_particles
  cdef int __result__ = c_get_number_of_particles(&output_number_of_particles);
  number_of_particles.value = output_number_of_particles
  return __result__


def get_stopping_condition_number_of_steps_parameter(value):
  
  cdef int output_value
  cdef int __result__ = c_get_stopping_condition_number_of_steps_parameter(&output_value);
  value.value = output_value
  return __result__


def disable_stopping_condition(type):
  
  cdef int __result__ = c_disable_stopping_condition(type);
  return __result__


def set_acceleration(index_of_the_particle, ax, ay, az):
  
  cdef int __result__ = c_set_acceleration(index_of_the_particle, ax, ay, az);
  return __result__


def get_center_of_mass_position(x, y, z):
  
  cdef double output_x
  cdef double output_y
  cdef double output_z
  cdef int __result__ = c_get_center_of_mass_position(&output_x, &output_y, &output_z);
  x.value = output_x
  y.value = output_y
  z.value = output_z
  return __result__


def get_center_of_mass_velocity(vx, vy, vz):
  
  cdef double output_vx
  cdef double output_vy
  cdef double output_vz
  cdef int __result__ = c_get_center_of_mass_velocity(&output_vx, &output_vy, &output_vz);
  vx.value = output_vx
  vy.value = output_vy
  vz.value = output_vz
  return __result__


def get_radius(index_of_the_particle, radius):
  
  cdef double output_radius
  cdef int __result__ = c_get_radius(index_of_the_particle, &output_radius);
  radius.value = output_radius
  return __result__


def set_stopping_condition_minimum_internal_energy_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_minimum_internal_energy_parameter(value);
  return __result__


def set_begin_time(time):
  
  cdef int __result__ = c_set_begin_time(time);
  return __result__


def set_stopping_condition_minimum_density_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_minimum_density_parameter(value);
  return __result__


def set_radius(index_of_the_particle, radius):
  
  cdef int __result__ = c_set_radius(index_of_the_particle, radius);
  return __result__


def has_stopping_condition(type, result):
  
  cdef int output_result
  cdef int __result__ = c_has_stopping_condition(type, &output_result);
  result.value = output_result
  return __result__


def cleanup_code():
  
  cdef int __result__ = c_cleanup_code();
  return __result__


def set_stopping_condition_maximum_density_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_maximum_density_parameter(value);
  return __result__


def recommit_parameters():
  
  cdef int __result__ = c_recommit_parameters();
  return __result__


def initialize_code():
  
  cdef int __result__ = c_initialize_code();
  return __result__


def get_potential_energy(potential_energy):
  
  cdef double output_potential_energy
  cdef int __result__ = c_get_potential_energy(&output_potential_energy);
  potential_energy.value = output_potential_energy
  return __result__


def get_gravity_at_point(eps, x, y, z, ax, ay, az):
  
  cdef double output_ax
  cdef double output_ay
  cdef double output_az
  cdef int __result__ = c_get_gravity_at_point(eps, x, y, z, &output_ax, &output_ay, &output_az);
  ax.value = output_ax
  ay.value = output_ay
  az.value = output_az
  return __result__


def get_velocity(index_of_the_particle, vx, vy, vz):
  
  cdef double output_vx
  cdef double output_vy
  cdef double output_vz
  cdef int __result__ = c_get_velocity(index_of_the_particle, &output_vx, &output_vy, &output_vz);
  vx.value = output_vx
  vy.value = output_vy
  vz.value = output_vz
  return __result__


def get_stopping_condition_out_of_box_parameter(value):
  
  cdef double output_value
  cdef int __result__ = c_get_stopping_condition_out_of_box_parameter(&output_value);
  value.value = output_value
  return __result__


def set_dt_param(dt_dia):
  
  cdef int __result__ = c_set_dt_param(dt_dia);
  return __result__


def get_position(index_of_the_particle, x, y, z):
  
  cdef double output_x
  cdef double output_y
  cdef double output_z
  cdef int __result__ = c_get_position(index_of_the_particle, &output_x, &output_y, &output_z);
  x.value = output_x
  y.value = output_y
  z.value = output_z
  return __result__


def set_stopping_condition_maximum_internal_energy_parameter(value):
  
  cdef int __result__ = c_set_stopping_condition_maximum_internal_energy_parameter(value);
  return __result__


def get_dt_param(dt_dia):
  
  cdef double output_dt_dia
  cdef int __result__ = c_get_dt_param(&output_dt_dia);
  dt_dia.value = output_dt_dia
  return __result__


def set_position(index_of_the_particle, x, y, z):
  
  cdef int __result__ = c_set_position(index_of_the_particle, x, y, z);
  return __result__


def get_stopping_condition_info(index, type, number_of_particles):
  
  cdef int output_type
  cdef int output_number_of_particles
  cdef int __result__ = c_get_stopping_condition_info(index, &output_type, &output_number_of_particles);
  type.value = output_type
  number_of_particles.value = output_number_of_particles
  return __result__


def get_acceleration(index_of_the_particle, ax, ay, az):
  
  cdef double output_ax
  cdef double output_ay
  cdef double output_az
  cdef int __result__ = c_get_acceleration(index_of_the_particle, &output_ax, &output_ay, &output_az);
  ax.value = output_ax
  ay.value = output_ay
  az.value = output_az
  return __result__


def commit_parameters():
  
  cdef int __result__ = c_commit_parameters();
  return __result__


def get_stopping_condition_particle_index(index, index_of_the_column, index_of_particle):
  
  cdef int output_index_of_particle
  cdef int __result__ = c_get_stopping_condition_particle_index(index, index_of_the_column, &output_index_of_particle);
  index_of_particle.value = output_index_of_particle
  return __result__


def set_velocity(index_of_the_particle, vx, vy, vz):
  
  cdef int __result__ = c_set_velocity(index_of_the_particle, vx, vy, vz);
  return __result__

