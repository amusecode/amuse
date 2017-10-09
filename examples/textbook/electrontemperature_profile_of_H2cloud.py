"""
Evolves the steady state solution of a star irridiating a H2 region.
"""
from __future__ import print_function

import numpy
import os
from matplotlib import pyplot
from optparse import OptionParser

from amuse.lab import *

from amuse.community.mocassin.interface import Mocassin, mocassin_rydberg_unit

from amuse.units import units
from amuse.units import nbody_system

from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particle
from amuse.datamodel import Grid
from amuse.io import write_set_to_file

from distinct_colours import get_distinct

def make_grid(number_of_grid_cells, length, constant_hydrogen_density, inner_radius, outer_radius):
    grid = Grid.create([number_of_grid_cells] * 3, length.as_vector_with_length(3))
    
    grid.radius = grid.position.lengths()
    grid.hydrogen_density = constant_hydrogen_density
    grid.hydrogen_density[grid.radius <= inner_radius] = 0 | units.cm ** -3
    grid.hydrogen_density[grid.radius >= outer_radius] = 0 | units.cm ** -3
    return grid
    
def setup_abundancies(code):
    table = code.abundancies_table()
    for atom in table.keys():
        table[atom] = 0.0
    table['H'] = 1.0
    table['He'] = 0.1
    table['C'] = 2.2e-4
    table['N'] = 4.e-5
    table['O'] = 3.3e-4
    table['Ne'] = 5.e-5
    table['S'] = 9.e-6

def plot_temperature_line(radii, electron_temperatures, ci):
  
    colors = get_distinct(4)
    s = 100
    if ci==1:
        s = 50
    
    pyplot.scatter(radii.value_in(units.parsec),
                   electron_temperatures.value_in(units.K),
                   c=colors[ci], lw=0, s=s)
    pyplot.xlabel('R [pc]')
    pyplot.ylabel('T [K]')
    pyplot.xlim(0.8, 3.2)
#    pyplot.ylim(6000,9000)
#    pyplot.show()
#    pyplot.savefig("electrontemperature_profile_of_H2cloud")
    
def evolve_star(mass, age, position):
    star=Particles(1)
    star.position = position
    star.mass = mass
    stellar = SeBa()
    stellar.particles.add_particles(star)
    stellar.evolve_model(age)
    star.luminosity = stellar.particles.luminosity
    star.temperature = stellar.particles.temperature
    stellar.stop()
    return star

def setup_grid(radiative_transfer, outer_radius, Ngrid):
    radiative_transfer.parameters.length_x = outer_radius
    radiative_transfer.parameters.length_y = outer_radius
    radiative_transfer.parameters.length_z = outer_radius
    radiative_transfer.parameters.mesh_size = [Ngrid, Ngrid, Ngrid]

def initiate_radiative_transfer_code(outer_radius, Ngrid):
    radiative_transfer = Mocassin(number_of_workers = 4) #, debugger = "xterm")
    
    radiative_transfer.set_input_directory(radiative_transfer.get_default_input_directory())
    radiative_transfer.set_mocassin_output_directory(radiative_transfer.output_directory  + os.sep)
    radiative_transfer.initialize_code()
    radiative_transfer.set_symmetricXYZ(True)

    setup_grid(radiative_transfer, outer_radius, Ngrid)
    setup_abundancies(radiative_transfer)
    
    radiative_transfer.parameters.initial_nebular_temperature = 6000.0 | units.K
    radiative_transfer.parameters.high_limit_of_the_frequency_mesh = 15 | mocassin_rydberg_unit
    radiative_transfer.parameters.low_limit_of_the_frequency_mesh = 1.001e-5 | mocassin_rydberg_unit
    
    radiative_transfer.parameters.total_number_of_photons = 10000000
    radiative_transfer.parameters.total_number_of_points_in_frequency_mesh = 600
    radiative_transfer.parameters.convergence_limit = 0.09
    radiative_transfer.parameters.number_of_ionisation_stages = 6
    radiative_transfer.commit_parameters()

    return radiative_transfer

def main(number_of_grid_cells = 15, min_convergence = 20):

    cloud_center = [0.0, 0.0, 0.0] | units.AU
    star = evolve_star(120|units.MSun, 3.3|units.Myr, cloud_center)

    outer_radius = 3.0 | units.parsec    
    grid=make_grid(number_of_grid_cells = number_of_grid_cells,
                   length = outer_radius,
                   constant_hydrogen_density = 100 | units.cm**-3,
                   inner_radius = 1.0 | units.parsec,
                   outer_radius = outer_radius)
    radiative_transfer = initiate_radiative_transfer_code(outer_radius, number_of_grid_cells)
    
    radiative_transfer.grid.hydrogen_density = grid.hydrogen_density
    radiative_transfer.commit_grid()
    radiative_transfer.particles.add_particle(star)
    radiative_transfer.commit_particles()
    
    max_number_of_photons = radiative_transfer.parameters.total_number_of_photons * 100
    percentage_converged = previous_percentage_converged = 0.0


    grid.electron_temperature = radiative_transfer.grid.electron_temperature
    radius = grid.radius.flatten()
    electron_temperature = grid.electron_temperature.flatten()
    selection = electron_temperature > 0 | units.K
    pyplot.rcParams.update({'font.size': 30})
    figure = pyplot.figure(figsize=(16, 12))
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)
    plot_temperature_line(radius[selection], electron_temperature[selection], 0)

    step = 0
    while percentage_converged < min_convergence:

        radiative_transfer.step()
        percentage_converged = radiative_transfer.get_percentage_converged()
        print("percentage converged :", percentage_converged, ", step :", step, ", photons:", radiative_transfer.parameters.total_number_of_photons)
        
        if previous_percentage_converged > 5 and percentage_converged < 95:
            convergence_increase = (percentage_converged-previous_percentage_converged)/previous_percentage_converged
            if convergence_increase < 0.2 and radiative_transfer.parameters.total_number_of_photons < max_number_of_photons:
                radiative_transfer.parameters.total_number_of_photons *= 2
        step += 1
        previous_percentage_converged = percentage_converged    

    grid.electron_temperature = radiative_transfer.grid.electron_temperature
    radius = grid.radius.flatten()
    electron_temperature = grid.electron_temperature.flatten()
    selection = electron_temperature > 0 | units.K
    plot_temperature_line(radius[selection], electron_temperature[selection], 1)
    pyplot.savefig("electrontemperature_profile_of_H2cloud")
    #pyplot.show()
    write_set_to_file(grid, 'h2region.h5', 'amuse')
    

def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n", "--gridcells", 
        default = 13,
        dest="number_of_grid_cells",
        help="number of cells in each direction",
        type="int"
    )
    result.add_option(
        "-c", "--min-convergence", 
        default = 60,
        dest="min_convergence",
        help="stop the iteratation when the solution is converged to the given percentage (in whole numbers between 10 and 100)",
        type="int"
    )
    
    return result
    
if __name__ == "__plot__":
    main(13, 30)
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
    
