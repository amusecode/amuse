import os
import os.path
import shutil
import math
import subprocess
import numpy
import pickle

from amuse.units import units, constants, nbody_system
from amuse.units.core import enumeration_unit
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.datamodel import Particles, Particle, ParticlesSuperset
from amuse.io import write_set_to_file, read_set_from_file
from amuse.support.exceptions import AmuseException

from amuse.ext.star_to_sph import convert_stellar_model_to_SPH, StellarModelInSPH
from amuse.couple.bridge import Bridge, CalculateFieldForParticles
from amuse.community.evtwin.interface import EVtwin
from amuse.community.mesa.interface import MESA
from amuse.community.gadget2.interface import Gadget2
from amuse.community.twobody.twobody import TwoBody
from amuse.community.huayno.interface import Huayno

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel, plot,loglog,semilogx,semilogy, sph_particles_plot
from amuse.plot import pynbody_column_density_plot, HAS_PYNBODY


def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print("Created new directory for output:", new_directory)
    os.mkdir(os.path.join(new_directory, "plots"))
    os.mkdir(os.path.join(new_directory, "snapshots"))
    shutil.copy(__file__, new_directory)
    if not os.path.exists(os.path.join(current_directory, "giant_models")):
        os.mkdir(os.path.join(current_directory, "giant_models"))
    os.chdir(new_directory)

def get_relative_velocity(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 + ecc)/(1.0 - ecc)) / semimajor_axis).sqrt()

def set_up_initial_conditions():
    stars = set_up_inner_binary()
    giant = set_up_outer_star(stars.total_mass())
    view_on_giant = stars.add_particle(giant)
    stars.move_to_center()
    return stars, view_on_giant

def set_up_inner_binary():
    semimajor_axis = 0.133256133158 | units.AU
    eccentricity = 0
    masses = [3.2, 3.1] | units.MSun
    orbital_period = (4 * numpy.pi**2 * semimajor_axis**3 / 
        (constants.G * masses.sum())).sqrt().as_quantity_in(units.day)
    
    print("   Initializing inner binary")
    print("   Orbital period inner binary:", orbital_period)
    stars =  Particles(2)
    stars.mass = masses
    stars.position = [0.0, 0.0, 0.0] | units.AU
    stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
    stars[0].x = semimajor_axis
    stars[0].vy = get_relative_velocity(stars.total_mass(), semimajor_axis, eccentricity)
    stars.move_to_center()
    return stars

def set_up_outer_star(inner_binary_mass):
    semimajor_axis = 1.22726535008 | units.AU
    eccentricity = 0.15
    inclination = math.radians(9.0)
    
    print("   Initializing outer star")
    giant = Particle()
    giant.mass = 5.5 | units.MSun
    giant.position = semimajor_axis * ([math.cos(inclination), 0, math.sin(inclination)] | units.none)
    giant.velocity = get_relative_velocity(giant.mass + inner_binary_mass, 
        semimajor_axis, eccentricity) * ([0, 1, 0] | units.none)
    return giant

def estimate_roche_radius(triple, view_on_giant):
    # 'mass ratio' of giant to inner binary
    q = (view_on_giant.mass / (triple-view_on_giant).total_mass())
    # Assume ~ circular orbit:
    a = (view_on_giant.position - (triple-view_on_giant).center_of_mass()).length()
    q13 = q**(1./3.)
    q23 = q13**2
    return (a*(0.49*q23/(0.6*q23+math.log(1+q13)))).as_quantity_in(units.RSun)

def evolve_stars(triple, view_on_giant, stellar_evolution_code, radius_factor):
    stop_radius = radius_factor * estimate_roche_radius(triple, view_on_giant)
    stellar_evolution = stellar_evolution_code(redirection='file', redirect_file='stellar_evolution_code_out.log')
    se_giant = stellar_evolution.particles.add_particle(view_on_giant)
    while (se_giant.radius < stop_radius):
        stellar_evolution.evolve_model(keep_synchronous = False)
    
    se_binary = stellar_evolution.particles.add_particles(triple - view_on_giant)
    for particle in se_binary:
        particle.evolve_for(se_giant.age)
    return stellar_evolution.particles, stellar_evolution

def convert_giant_to_sph(view_on_se_giant, number_of_sph_particles):
    giant_in_sph = convert_stellar_model_to_SPH(
        view_on_se_giant,
        number_of_sph_particles,
        with_core_particle = True,
        target_core_mass  = 1.4 | units.MSun,
        do_store_composition = False
    )
    return giant_in_sph

def relax_in_isolation(giant_in_sph, sph_code, output_base_name):
    total_mass = giant_in_sph.gas_particles.total_mass() + giant_in_sph.core_particle.mass
    total_radius = max(giant_in_sph.gas_particles.position.lengths_squared()).sqrt()
    dynamical_timescale = (total_radius**3 / (2 * constants.G * total_mass)).sqrt().as_quantity_in(units.day)
    t_end = (10.0 * dynamical_timescale).as_quantity_in(units.day)
    n_steps = 100
    hydro_code_options = dict(number_of_workers=2, redirection='file', redirect_file='hydrodynamics_code_relax_out.log')
    
    unit_converter = ConvertBetweenGenericAndSiUnits(total_radius, total_mass, t_end)
    hydrodynamics = sph_code(unit_converter, **hydro_code_options)
    hydrodynamics.parameters.epsilon_squared = giant_in_sph.core_radius**2
    hydrodynamics.parameters.max_size_timestep = t_end
    hydrodynamics.parameters.time_max = 1.1 * t_end
    hydrodynamics.parameters.time_limit_cpu = 7.0 | units.day
    hydrodynamics.gas_particles.add_particles(giant_in_sph.gas_particles)
    hydrodynamics.dm_particles.add_particle(giant_in_sph.core_particle)
    
    potential_energies = hydrodynamics.potential_energy.as_vector_with_length(1).as_quantity_in(units.erg)
    kinetic_energies = hydrodynamics.kinetic_energy.as_vector_with_length(1).as_quantity_in(units.erg)
    thermal_energies = hydrodynamics.thermal_energy.as_vector_with_length(1).as_quantity_in(units.erg)
    
    print("   Relaxing for", t_end, "(10 * dynamical timescale)")
    times = (t_end * list(range(1, n_steps+1)) / n_steps).as_quantity_in(units.day)
    for i_step, time in enumerate(times):
        hydrodynamics.evolve_model(time)
        print("   Relaxed for:", time)
        potential_energies.append(hydrodynamics.potential_energy)
        kinetic_energies.append(hydrodynamics.kinetic_energy)
        thermal_energies.append(hydrodynamics.thermal_energy)
    
    hydrodynamics.gas_particles.copy_values_of_attributes_to(
        ['mass', 'x','y','z', 'vx','vy','vz', 'u', 'h_smooth'], 
        giant_in_sph.gas_particles)
    giant_in_sph.core_particle.position = hydrodynamics.dm_particles[0].position
    giant_in_sph.core_particle.velocity = hydrodynamics.dm_particles[0].velocity
    hydrodynamics.stop()
    
    snapshotfile = output_base_name + "_gas.amuse"
    write_set_to_file(giant_in_sph.gas_particles, snapshotfile, format='amuse')
    shutil.copy(snapshotfile, os.path.join("..", "giant_models"))
    
    snapshotfile = output_base_name + "_core.amuse"
    # temporarily store core_radius on the core particle
    giant_in_sph.core_particle.radius = giant_in_sph.core_radius
    write_set_to_file(giant_in_sph.core_particle.as_set(), snapshotfile, format='amuse')
    giant_in_sph.core_particle.radius = 0.0 | units.m
    shutil.copy(snapshotfile, os.path.join("..", "giant_models"))
    
    energy_evolution_plot(times, kinetic_energies, potential_energies, thermal_energies, 
        figname = output_base_name + "_energy_evolution.png")

def load_giant_model(file_base_name):
    snapshotfile = os.path.join("..", "giant_models", file_base_name + "_gas.amuse")
    sph_particles = read_set_from_file(snapshotfile, format='amuse')
    snapshotfile = os.path.join("..", "giant_models", file_base_name + "_core.amuse")
    gd_particles = read_set_from_file(snapshotfile, format='amuse')
    core_radius = gd_particles[0].radius
    gd_particles[0].radius =  0.0 | units.m
    giant_model = StellarModelInSPH(
        gas_particles = sph_particles, 
        core_particle = gd_particles[0], 
        core_radius = core_radius)
    return giant_model

def prepare_binary_system(dynamics_code, binary_particles):
    unit_converter = nbody_system.nbody_to_si(
        binary_particles.total_mass(),
        (binary_particles[0].position - binary_particles[1].position).length())
    system = dynamics_code(unit_converter, redirection="none")
    system.particles.add_particles(binary_particles)
    return system

def prepare_giant_system(sph_code, giant_model, view_on_giant, time_unit, n_steps):
    hydro_code_options = dict(number_of_workers=2, 
        redirection='file', redirect_file='hydrodynamics_code_out.log')
    unit_converter = ConvertBetweenGenericAndSiUnits(
        1.0 | units.RSun, 
        giant_model.gas_particles.total_mass() + giant_model.core_particle.mass, 
        time_unit)
    system = sph_code(unit_converter, **hydro_code_options)
    system.parameters.epsilon_squared = giant_model.core_radius**2
    system.parameters.max_size_timestep = time_unit / n_steps
    system.parameters.time_max = 1.1 * time_unit
    system.parameters.time_limit_cpu = 7.0 | units.day

    giant_model.gas_particles.position += view_on_giant.position
    giant_model.gas_particles.velocity += view_on_giant.velocity
    
    giant_model.core_particle.position += view_on_giant.position
    giant_model.core_particle.velocity += view_on_giant.velocity
    
    system.gas_particles.add_particles(giant_model.gas_particles)
    system.dm_particles.add_particle(giant_model.core_particle)
    return system

def calculate_orbital_elements(m1, m2, pos1, pos2, vel1, vel2, m3, pos3, vel3):
    print("   Calculating semimajor axis and eccentricity evolution of the giant's orbit")

    m12 = m1+m2
    rel_position = (m1 * pos1 + m2 * pos2)/m12 - pos3
    rel_velocity = (m1 * vel1 + m2 * vel2)/m12 - vel3
    mtot = m12 + m3
    separation = rel_position.lengths()
    speed_squared = rel_velocity.lengths_squared()
    
    # Now calculate the important quantities:
    semimajor_axis = (constants.G * mtot * separation / 
        (2 * constants.G * mtot - separation * speed_squared)).as_quantity_in(units.AU)
    eccentricity = numpy.sqrt(1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=1) / 
        (constants.G * mtot * semimajor_axis))
    return semimajor_axis, eccentricity

def evolve_coupled_system(binary_system, giant_system, t_end, n_steps, 
        do_energy_evolution_plot, previous_data=None):
    directsum = CalculateFieldForParticles(particles=giant_system.particles, gravity_constant=constants.G)
    directsum.smoothing_length_squared = giant_system.parameters.gas_epsilon**2
    coupled_system = Bridge(timestep=(t_end / (2 * n_steps)), verbose=False, use_threading=True)
    coupled_system.add_system(binary_system, (directsum,), False)
    coupled_system.add_system(giant_system, (binary_system,), False)
    
    times = (t_end * list(range(1, n_steps+1)) / n_steps).as_quantity_in(units.day)
    
    if previous_data:
        with open(previous_data, 'rb') as file:
            (all_times, potential_energies, kinetic_energies, thermal_energies, 
                giant_center_of_mass, ms1_position, ms2_position, 
                giant_center_of_mass_velocity, ms1_velocity, ms2_velocity) = pickle.load(file)
        all_times.extend(times + all_times[-1])
    else:
        all_times = times
        if do_energy_evolution_plot:
            potential_energies = coupled_system.particles.potential_energy().as_vector_with_length(1).as_quantity_in(units.erg)
            kinetic_energies = coupled_system.particles.kinetic_energy().as_vector_with_length(1).as_quantity_in(units.erg)
            thermal_energies = coupled_system.gas_particles.thermal_energy().as_vector_with_length(1).as_quantity_in(units.erg)
        else:
            potential_energies = kinetic_energies = thermal_energies = None
        
        giant_center_of_mass = [] | units.RSun
        ms1_position = [] | units.RSun
        ms2_position = [] | units.RSun
        giant_center_of_mass_velocity = [] | units.km / units.s
        ms1_velocity = [] | units.km / units.s
        ms2_velocity = [] | units.km / units.s
    i_offset = len(giant_center_of_mass)
    
    giant_total_mass = giant_system.particles.total_mass()
    ms1_mass = binary_system.particles[0].mass
    ms2_mass = binary_system.particles[1].mass
    
    print("   Evolving for", t_end)
    for i_step, time in enumerate(times):
        coupled_system.evolve_model(time)
        print("   Evolved to:", time, end=' ')
        
        if do_energy_evolution_plot:
            potential_energies.append(coupled_system.particles.potential_energy())
            kinetic_energies.append(coupled_system.particles.kinetic_energy())
            thermal_energies.append(coupled_system.gas_particles.thermal_energy())
        
        giant_center_of_mass.append(giant_system.particles.center_of_mass())
        ms1_position.append(binary_system.particles[0].position)
        ms2_position.append(binary_system.particles[1].position)
        giant_center_of_mass_velocity.append(giant_system.particles.center_of_mass_velocity())
        ms1_velocity.append(binary_system.particles[0].velocity)
        ms2_velocity.append(binary_system.particles[1].velocity)

        a_giant, e_giant = calculate_orbital_elements(ms1_mass, ms2_mass,
                                                      ms1_position, ms2_position,
                                                      ms1_velocity, ms2_velocity,
                                                      giant_total_mass,
                                                      giant_center_of_mass,
                                                      giant_center_of_mass_velocity)
        print("Outer Orbit:", time.in_(units.day),  a_giant[-1].in_(units.AU), e_giant[-1], ms1_mass.in_(units.MSun), ms2_mass.in_(units.MSun), giant_total_mass.in_(units.MSun))
        
        if i_step % 10 == 9:
            snapshotfile = os.path.join("snapshots", "hydro_triple_{0:=04}_gas.amuse".format(i_step + i_offset))
            write_set_to_file(giant_system.gas_particles, snapshotfile, format='amuse')
            snapshotfile = os.path.join("snapshots", "hydro_triple_{0:=04}_core.amuse".format(i_step + i_offset))
            write_set_to_file(giant_system.dm_particles, snapshotfile, format='amuse')
            snapshotfile = os.path.join("snapshots", "hydro_triple_{0:=04}_binary.amuse".format(i_step + i_offset))
            write_set_to_file(binary_system.particles, snapshotfile, format='amuse')
            
            datafile = os.path.join("snapshots", "hydro_triple_{0:=04}_info.amuse".format(i_step + i_offset))
            with open(datafile, 'wb') as outfile:
                pickle.dump((all_times[:len(giant_center_of_mass)], potential_energies, 
                    kinetic_energies, thermal_energies, giant_center_of_mass, 
                    ms1_position, ms2_position, giant_center_of_mass_velocity,
                    ms1_velocity, ms2_velocity), outfile)
       
        figname1 = os.path.join("plots", "hydro_triple_small{0:=04}.png".format(i_step + i_offset))
        figname2 = os.path.join("plots", "hydro_triple_large{0:=04}.png".format(i_step + i_offset))
        print("  -   Hydroplots are saved to: ", figname1, "and", figname2)
        for plot_range, plot_name in [(8|units.AU, figname1), (40|units.AU, figname2)]:
            if HAS_PYNBODY:
                pynbody_column_density_plot(coupled_system.gas_particles, width=plot_range, vmin=26, vmax=32)
                scatter(coupled_system.dm_particles.x, coupled_system.dm_particles.y, c="w")
            else:
                pyplot.figure(figsize = [16, 16])
                sph_particles_plot(coupled_system.gas_particles, gd_particles=coupled_system.dm_particles, 
                    view=plot_range*[-0.5, 0.5, -0.5, 0.5])
            pyplot.savefig(plot_name)
            pyplot.close()
        
    
    coupled_system.stop()
    
    make_movie()
    
    if do_energy_evolution_plot:
        energy_evolution_plot(all_times[:len(kinetic_energies)-1], kinetic_energies, 
            potential_energies, thermal_energies)
    
    print("   Calculating semimajor axis and eccentricity evolution for inner binary")
    # Some temporary variables to calculate semimajor_axis and eccentricity evolution
    total_mass = ms1_mass + ms2_mass
    rel_position = ms1_position - ms2_position
    rel_velocity = ms1_velocity - ms2_velocity
    separation_in = rel_position.lengths()
    speed_squared_in = rel_velocity.lengths_squared()
    
    # Now calculate the important quantities:
    semimajor_axis_binary = (constants.G * total_mass * separation_in / 
        (2 * constants.G * total_mass - separation_in * speed_squared_in)).as_quantity_in(units.AU)
    eccentricity_binary = numpy.sqrt(1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=1) / 
        (constants.G * total_mass * semimajor_axis_binary))

    print("   Calculating semimajor axis and eccentricity evolution of the giant's orbit")
    # Some temporary variables to calculate semimajor_axis and eccentricity evolution
    rel_position = ((ms1_mass * ms1_position + ms2_mass * ms2_position)/total_mass - 
        giant_center_of_mass)
    rel_velocity = ((ms1_mass * ms1_velocity + ms2_mass * ms2_velocity)/total_mass - 
        giant_center_of_mass_velocity)
    total_mass += giant_total_mass
    separation = rel_position.lengths()
    speed_squared = rel_velocity.lengths_squared()
    
    # Now calculate the important quantities:
    semimajor_axis_giant = (constants.G * total_mass * separation / 
        (2 * constants.G * total_mass - separation * speed_squared)).as_quantity_in(units.AU)
    eccentricity_giant = numpy.sqrt(1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=1) / 
        (constants.G * total_mass * semimajor_axis_giant))
    
    orbit_parameters_plot(semimajor_axis_binary, semimajor_axis_giant, all_times[:len(semimajor_axis_binary)])
    orbit_ecc_plot(eccentricity_binary, eccentricity_giant, all_times[:len(eccentricity_binary)])
    
    orbit_parameters_plot(separation_in.as_quantity_in(units.AU), 
        separation.as_quantity_in(units.AU), 
        all_times[:len(eccentricity_binary)], 
        par_symbol="r", par_name="separation")
    orbit_parameters_plot(speed_squared_in.as_quantity_in(units.km**2 / units.s**2), 
        speed_squared.as_quantity_in(units.km**2 / units.s**2), 
        all_times[:len(eccentricity_binary)], 
        par_symbol="v^2", par_name="speed_squared")

def make_movie():
    print("   Creating movie from snapshots")
    try:
        subprocess.call(['mencoder', "mf://hydro_triple_small*.png", '-ovc', 'lavc', 
            '-o', '../hydro_triple_small.avi', '-msglevel', 'all=1'], cwd="./plots")
        subprocess.call(['mencoder', "mf://hydro_triple_large*.png", '-ovc', 'lavc', 
            '-o', '../hydro_triple_large.avi', '-msglevel', 'all=1'], cwd="./plots")
    except Exception as exc:
        print("   Failed to create movie, error was:", str(exc))

def continue_evolution(sph_code, dynamics_code, t_end, n_steps, 
        relaxed_giant_output_base_name, do_energy_evolution_plot):
    print("Loading snapshots...", end=' ')
    files = os.listdir("snapshots")
    files.sort()
    files = files[-4:]
    print(files)
    binary = read_set_from_file(os.path.join("snapshots", files[0]), format='amuse')
    gd_particles = read_set_from_file(os.path.join("snapshots", files[1]), format='amuse')
    sph_particles = read_set_from_file(os.path.join("snapshots", files[2]), format='amuse')
    
    snapshotfile = os.path.join("..", "giant_models", relaxed_giant_output_base_name + "_core.amuse")
    core_particle = read_set_from_file(snapshotfile, format='amuse')
    
    giant_model = StellarModelInSPH(
        gas_particles = sph_particles, 
        core_particle = gd_particles[0], 
        core_radius = core_particle.radius)
    
    view_on_giant = Particle()
    view_on_giant.position = [0]*3 | units.m
    view_on_giant.velocity = [0]*3 | units.m / units.s
    
    print("\nSetting up {0} to simulate inner binary system".format(dynamics_code.__name__))
    binary_system = prepare_binary_system(dynamics_code, binary)
    
    print("\nSetting up {0} to simulate giant in SPH".format(sph_code.__name__))
    giant_system = prepare_giant_system(sph_code, giant_model, view_on_giant, t_end, n_steps)
    
    print("\nEvolving with bridge between", sph_code.__name__, "and", dynamics_code.__name__)
    evolve_coupled_system(binary_system, giant_system, t_end, n_steps, do_energy_evolution_plot, 
        previous_data = os.path.join("snapshots", files[3]))
    print("Done")

def energy_evolution_plot(time, kinetic, potential, thermal, figname = "energy_evolution.png"):
    time.prepend(0.0 | units.day)
    pyplot.figure(figsize = (5, 5))
    plot(time, kinetic, label='K')
    plot(time, potential, label='U')
    plot(time, thermal, label='Q')
    plot(time, kinetic + potential + thermal, label='E')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(prop={'size':"x-small"}, loc=4)
    pyplot.savefig(figname)
    pyplot.close()

def orbit_ecc_plot(eccentricity_in,eccentricity_out,time):
    figure = pyplot.figure(figsize = (10, 6), dpi = 100)
    subplot = figure.add_subplot(2, 1, 1)
    plot(time,eccentricity_in)
    xlabel('t')
    ylabel('e$_\mathrm{binary}$')
    
    subplot = figure.add_subplot(2, 1, 2)
    plot(time,eccentricity_out )  
    xlabel('t')
    ylabel('e$_\mathrm{giant}$')
    pyplot.minorticks_on()
    pyplot.savefig("eccentricity_evolution.png")
    pyplot.close()


def orbit_parameters_plot(semi_major_in,semi_major_out, time, par_symbol="a", par_name="semimajor_axis"):
    figure = pyplot.figure(figsize = (10, 6), dpi = 100)
    subplot = figure.add_subplot(2, 1, 1)
    plot(time,semi_major_in )
    xlabel('t')
    ylabel('$'+par_symbol+'_\mathrm{binary}$')
    
    subplot = figure.add_subplot(2, 1, 2)
    plot(time,semi_major_out )
    xlabel('t')
    ylabel('$'+par_symbol+'_\mathrm{giant}$')
    pyplot.minorticks_on()
    pyplot.savefig(par_name+"_evolution.png")
    pyplot.close()


if __name__ == "__main__":
    stellar_evolution_code = MESA
    sph_code = Gadget2
    #dynamics_code = TwoBody
    dynamics_code = Huayno

    number_of_sph_particles = 50000
    # Stop stellar evolution when giant's radius is (radius_factor * Roche lobe radius)
    radius_factor = 1.0
    relaxed_giant_output_base_name = "relaxed_giant_" + str(number_of_sph_particles) + "_" + str(radius_factor)
    t_end = 300.0 | units.day
    n_steps = 3000
    
    do_energy_evolution_plot = False
    
    if os.path.exists("snapshots"):
        print("Found snapshots folder, continuing evolution of previous run")
        continue_evolution(sph_code, dynamics_code, t_end, n_steps, 
            relaxed_giant_output_base_name, do_energy_evolution_plot)
        exit(0)
    
    new_working_directory()
    
    print("Initializing triple")
    triple, view_on_giant = set_up_initial_conditions()
    print("\nInitialization done:\n", triple)
    
    print("\nEvolving with", stellar_evolution_code.__name__)
    se_stars, se_code_instance = evolve_stars(triple, view_on_giant, stellar_evolution_code, radius_factor)
    triple.radius = se_stars.radius
    print("\nStellar evolution done:\n", se_stars)
    
    if os.path.exists(os.path.join("..", "giant_models", relaxed_giant_output_base_name + "_gas.amuse")):
        print("\nLoading SPH model for giant from:", end=' ') 
        print(os.path.join("..", "giant_models", relaxed_giant_output_base_name + "_gas.amuse"))
        giant_model = load_giant_model(relaxed_giant_output_base_name)
        se_code_instance.stop()
    else:
        print("\nConverting giant to", number_of_sph_particles, "SPH particles")
        view_on_se_giant = view_on_giant.as_set().get_intersecting_subset_in(se_stars)[0]
        giant_model = convert_giant_to_sph(view_on_se_giant, number_of_sph_particles)
        se_code_instance.stop()
        print("Relaxing giant with", sph_code.__name__)
        relax_in_isolation(giant_model, sph_code, relaxed_giant_output_base_name)
    
    print("\nSetting up {0} to simulate inner binary system".format(dynamics_code.__name__))
    binary_system = prepare_binary_system(dynamics_code, triple - view_on_giant)
    
    print("\nSetting up {0} to simulate giant in SPH".format(sph_code.__name__))
    giant_system = prepare_giant_system(sph_code, giant_model, view_on_giant, t_end, n_steps)
    
    print("\nEvolving with bridge between", sph_code.__name__, "and", dynamics_code.__name__)
    evolve_coupled_system(binary_system, giant_system, t_end, n_steps, do_energy_evolution_plot)
    print("Done")
