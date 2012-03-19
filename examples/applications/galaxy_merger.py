"""
Simple Galaxy Merger

example simulating a merger of two Milky Way-type galaxies, collisionless particles (i.e. no gas)
"""

import os
import os.path
import numpy

from amuse.units import nbody_system, units, generic_unit_converter, constants
from amuse.datamodel import Particles
from amuse.io import read_set_from_file, write_set_to_file
from amuse.plot import plot, native_plot, pynbody_column_density_plot, HAS_PYNBODY
from amuse.community.gadget2.interface import Gadget2
from amuse.ext.galactics_model import new_galactics_model

def make_plots(all_particles, disk_only, i=0):
    for j, particles in enumerate([all_particles, disk_only]):
        if HAS_PYNBODY:
            temp_particles = particles.copy_to_memory()
            temp_particles.u = 1 | units.ms**2
            temp = Gadget2()
            temp.gas_particles.add_particles(temp_particles)
            temp.commit_particles()
            pynbody_column_density_plot(temp.gas_particles, width=300|units.kpc, units = 'm_p cm^-2', vmin=17, vmax=23)
            temp.stop()
        else:
            native_plot.gca().set_aspect("equal")
            native_plot.xlim(-150, 150)
            native_plot.ylim(-150, 150)
            plot(particles.x.as_quantity_in(units.kpc), particles.y.as_quantity_in(units.kpc), 'r.')
    
        native_plot.savefig(os.path.join("plots", "plot_galaxy_merger_{0}_{1:=04}.png".format("disk" if j else "total", i)))
        native_plot.close()

def make_galaxies():
    # File 'disk_galactICs.amuse' contains 30k particles, (10k disk, no bulge, 20k halo)
    if os.path.exists('disk_galactICs.amuse'):
        galaxy1 = read_set_from_file('disk_galactICs.amuse', 'amuse')
    else:
        halo_number_of_particles = 20000
        converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(constants.G, 1.0e12 | units.MSun, 50.0 | units.kpc)
        galaxy1 = new_galactics_model(halo_number_of_particles, disk_number_of_particles=10000, 
            generate_bulge_flag=False, do_scale=True, unit_system_converter=converter)
        write_set_to_file(galaxy1, 'disk_galactICs.amuse', 'amuse')
    
    galaxy2 = Particles(len(galaxy1))
    galaxy2.mass = galaxy1.mass
    galaxy2.position = galaxy1.position
    galaxy2.velocity = galaxy1.velocity
    
    galaxy1.rotate(0.0, numpy.pi/4, 0.0)
    galaxy2.rotate(numpy.pi/6, 0.0, 0.0)
    galaxy1.position += [100.0, 0, 0] | units.kpc
    galaxy2.position -= [100.0, 0, 0] | units.kpc
    galaxy1.velocity += [0.0, 50.0, 0] | units.km/units.s
    galaxy2.velocity -= [0.0, 50.0, 0] | units.km/units.s
    return galaxy1, galaxy2

def simulate_merger(galaxy1, galaxy2):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    dynamics = Gadget2(converter, number_of_workers=2)
    dynamics.parameters.epsilon_squared = 0.0000001 | nbody_system.length**2
    set1 = dynamics.particles.add_particles(galaxy1)
    set2 = dynamics.particles.add_particles(galaxy2)
    
    galaxies_without_halo = set1[:10000] + set2[:10000]
    make_plots(dynamics.particles, galaxies_without_halo)
    
    for i in range(1, 101):
        dynamics.evolve_model(i * (25 | units.Myr))
        print dynamics.model_time.as_quantity_in(units.Myr)
        make_plots(dynamics.particles, galaxies_without_halo, i)
    
    dynamics.stop()

if __name__ == '__main__':
    galaxy1, galaxy2 = make_galaxies()
    if not os.path.exists("plots"):
        os.mkdir("plots")
    simulate_merger(galaxy1, galaxy2)
