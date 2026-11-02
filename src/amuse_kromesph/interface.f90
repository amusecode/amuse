function initialize_code() result(ret)
  use chem_mod, only : chem_initialize
  implicit none
  integer :: ret
  ret = chem_initialize()
end function initialize_code

function cleanup_code() result(ret)
  use chem_mod, only : chem_end
  implicit none
  integer :: ret
  ret = chem_end()
end function cleanup_code

function commit_parameters() result(ret)
  use chem_mod, only : chem_commit_parameters
  implicit none
  integer :: ret
  ret = chem_commit_parameters()
end function commit_parameters

function recommit_parameters() result(ret)
  use chem_mod, only : chem_commit_parameters
  implicit none
  integer :: ret
  ret = chem_commit_parameters()
end function recommit_parameters

function commit_particles() result(ret)
  use chem_mod, only : chem_commit_particles
  implicit none
  integer :: ret
  ret = chem_commit_particles()
end function commit_particles

function recommit_particles() result(ret)
  use chem_mod, only : chem_commit_particles
  implicit none
  integer :: ret
  ret = chem_commit_particles()
end function recommit_particles

function evolve_model(time) result(ret)
  use chem_mod, only : evolve_chem_model
  implicit none
  double precision, intent(in) :: time
  integer :: ret
  ret = evolve_chem_model(time)
end function evolve_model

function new_particle(index_of_the_particle, rho, u, gamma, mu, ionrate) result(ret)
  use chem_mod, only : add_particle
  implicit none
  integer, intent(out) :: index_of_the_particle
  double precision, intent(in) :: rho, u, gamma, mu, ionrate
  integer :: ret
  ret = add_particle(index_of_the_particle, rho, u, gamma, mu, ionrate)
end function new_particle

function delete_particle(index_of_the_particle) result(ret)
  use chem_mod, only : remove_particle
  implicit none
  integer, intent(in) :: index_of_the_particle
  integer :: ret
  ret = remove_particle(index_of_the_particle)
end function delete_particle

function get_state(index_of_the_particle, rho, u, gamma, mu, ionrate) result(ret)
  use chem_mod, only : get_particle_state
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(out) :: rho, u, gamma, mu, ionrate
  integer :: ret
  ret = get_particle_state(index_of_the_particle, rho, u, gamma, mu, ionrate)
end function get_state

function set_state(index_of_the_particle, rho, u, gamma, mu, ionrate) result(ret)
  use chem_mod, only : set_particle_state
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(in) :: rho, u, gamma, mu, ionrate
  integer :: ret
  ret = set_particle_state(index_of_the_particle, rho, u, gamma, mu, ionrate)
end function set_state

function get_density(index_of_the_particle, rho) result(ret)
  use chem_mod, only : get_particle_rho
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(out) :: rho
  integer :: ret
  ret = get_particle_rho(index_of_the_particle, rho)
end function get_density

function set_density(index_of_the_particle, rho) result(ret)
  use chem_mod, only : set_particle_rho
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(in) :: rho
  integer :: ret
  ret = set_particle_rho(index_of_the_particle, rho)
end function set_density

function get_internal_energy(index_of_the_particle, u) result(ret)
  use chem_mod, only : get_particle_internal_energy
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(out) :: u
  integer :: ret
  ret = get_particle_internal_energy(index_of_the_particle, u)
end function get_internal_energy

function set_internal_energy(index_of_the_particle, u) result(ret)
  use chem_mod, only : set_particle_internal_energy
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(in) :: u
  integer :: ret
  ret = set_particle_internal_energy(index_of_the_particle, u)
end function set_internal_energy

function get_gamma(index_of_the_particle, gamma) result(ret)
  use chem_mod, only : get_particle_gamma
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(out) :: gamma
  integer :: ret
  ret = get_particle_gamma(index_of_the_particle, gamma)
end function get_gamma

function set_gamma(index_of_the_particle, gamma) result(ret)
  use chem_mod, only : set_particle_gamma
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(in) :: gamma
  integer :: ret
  ret = set_particle_gamma(index_of_the_particle, gamma)
end function set_gamma

function get_mu(index_of_the_particle, mu) result(ret)
  use chem_mod, only : get_particle_mu
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(out) :: mu
  integer :: ret
  ret = get_particle_mu(index_of_the_particle, mu)
end function get_mu

function set_mu(index_of_the_particle, mu) result(ret)
  use chem_mod, only : set_particle_mu
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(in) :: mu
  integer :: ret
  ret = set_particle_mu(index_of_the_particle, mu)
end function set_mu

function get_ionrate(index_of_the_particle, ionrate) result(ret)
  use chem_mod, only : get_particle_ionrate
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(out) :: ionrate
  integer :: ret
  ret = get_particle_ionrate(index_of_the_particle, ionrate)
end function get_ionrate

function set_ionrate(index_of_the_particle, ionrate) result(ret)
  use chem_mod, only : set_particle_ionrate
  implicit none
  integer, intent(in) :: index_of_the_particle
  double precision, intent(in) :: ionrate
  integer :: ret
  ret = set_particle_ionrate(index_of_the_particle, ionrate)
end function set_ionrate

function get_abundance(index_of_the_particle, species_index, abundance) result(ret)
  use chem_mod, only : get_particle_abundance
  implicit none
  integer, intent(in) :: index_of_the_particle, species_index
  double precision, intent(out) :: abundance
  integer :: ret
  ret = get_particle_abundance(index_of_the_particle, species_index, abundance)
end function get_abundance

function set_abundance(index_of_the_particle, species_index, abundance) result(ret)
  use chem_mod, only : set_particle_abundance
  implicit none
  integer, intent(in) :: index_of_the_particle, species_index
  double precision, intent(in) :: abundance
  integer :: ret
  ret = set_particle_abundance(index_of_the_particle, species_index, abundance)
end function set_abundance

function set_abundances(index_of_the_particle, abundances, N) result(ret)
  use chem_mod, only : set_particle_abundances
  implicit none
  integer, intent(in) :: N
  integer, intent(in) :: index_of_the_particle(N)
  double precision, intent(in) :: abundances(N)
  integer :: ret
  ret = set_particle_abundances(index_of_the_particle, abundances, N)
end function set_abundances

function get_firstlast_species_index(first, last) result(ret)
  use chem_mod, only : krome_nmols
  implicit none
  integer, intent(out) :: first, last
  integer :: ret
  first = 1
  last = krome_nmols ! this is the last species defined in amuse_helpers.f90
  ret = 0
end function get_firstlast_species_index

function get_species_index(name, species_index) result(ret)
  use chem_mod, only : krome_get_index
  implicit none
  character(len=16), intent(in) :: name
  integer, intent(out) :: species_index
  integer :: ret
  species_index = krome_get_index(name)
  ret = 0
end function get_species_index

function get_species_name(species_index, name) result(ret)
  use chem_mod, only : krome_get_names, krome_nmols
  implicit none
  integer, intent(in) :: species_index
  character(len=16), intent(out) :: name
  integer :: ret
  character(len=16) :: names(krome_nmols)
  if (species_index .LT. 1 .OR. species_index .GT. krome_nmols) then
    ret = -1
    return
  endif
  names = krome_get_names()
  name = names(species_index)
  ret = 0
end function get_species_name

function get_time(time) result(ret)
  use chem_mod, only : chem_model_time
  implicit none
  double precision, intent(out) :: time
  integer :: ret
  ret = chem_model_time(time)
end function get_time

function get_number_of_particles(number_of_particles) result(ret)
  use chem_mod, only : nparticle
  implicit none
  integer, intent(out) :: number_of_particles
  integer :: ret
  number_of_particles = nparticle
  ret = 0
end function get_number_of_particles

function set_amu_in_g(amu_in_g) result(ret)
  use krome_constants_amuse, only : amu_2_g
  implicit none
  double precision, intent(in) :: amu_in_g
  integer :: ret
  amu_2_g = amu_in_g
  ret = 0
end function set_amu_in_g
