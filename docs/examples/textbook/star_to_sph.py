import os
import os.path
import shutil
import numpy

from amuse.lab import *
"""
from amuse.units import units, constants
from amuse.datamodel import Particles, Particle
from amuse.io import write_set_to_file
"""
from amuse.community.mesa.interface import MESA as stellar_evolution_code
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.relax_sph import relax

def evolve_star(mass, age):
    stellar = MESA()
    stellar.parameters.metallicity = 0.02
    star = stellar.particles.add_particle(Particle(mass=mass))
    stellar.evolve_model(age)
    return stellar

def evolve_star_and_convert_to_sph(stellar, omega, Nsph):
    new_age = stellar.model_time
    star = stellar.particles[0]
    if star.core_mass>zero:
        with_core_particle = True
        target_core_mass  = coremass
    else:
        with_core_particle = False
        target_core_mass  = 0|units.MSun

    m_sph = (star.mass-target_core_mass)/float(Nsph)
    print "N_sph=", Nsph, m_sph.in_(units.MSun)
    print "Target core mass:", target_core_mass
    model = convert_stellar_model_to_SPH(
        star, 
        Nsph,
        with_core_particle = with_core_particle,
        target_core_mass  = target_core_mass,
        do_store_composition = False,
        base_grid_options=dict(type="fcc")
    )
    print "Final star:", star
    age = stellar.model_time
    M_scale = star.mass
    R_scale = star.radius
    stellar.stop()

    core = model.core_particle
    print "core", core
    if core is None:
        print "Make zero mass core"
        core = Particle(mass=0|units.MSun, radius=1|units.RSun)
        core.position = (0,0,0) | units.AU
        core.velocity = (0,0,0) | units.kms
    sph_star = model.gas_particles
    print "Add Spin to star:", omega
    sph_star.add_spin(omega)

    relaxed_sph_star, core_particles = relax_sph_realization(sph_star, core)
    return relaxed_sph_star, core_particles, age

def relax_sph_realization(sph_star, core):
    
    dynamical_timescale = sph_star.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1|units.RSun)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.parameters.time_max = 3 * dynamical_timescale
    hydro.parameters.max_size_timestep = dynamical_timescale / 100
    hydro.parameters.time_limit_cpu = 1.0 | units.Gyr

    t_end_in_t_dyn = 2.5 # Relax for this many dynamical timescales
    t_end = t_end_in_t_dyn * sph_star.dynamical_timescale(mass_fraction=0.9)
    n_steps = 250
    velocity_damp_factor = 1.0 - (2.0*numpy.pi*t_end_in_t_dyn)/n_steps # Critical damping

    hydro.gas_particles.add_particles(sph_star)
    to_framework = hydro.gas_particles.new_channel_to(sph_star)

    if core.mass>0|units.MSun:
        hydro.dm_particles.add_particles(core)
    for i_step, time in enumerate(t_end * numpy.linspace(1.0/n_steps, 1.0, n_steps)):
        hydro.evolve_model(time)
        to_framework.copy()
        hydro.gas_particles.velocity = velocity_damp_factor * hydro.gas_particles.velocity
    
    return hydro.gas_particles, hydro.dm_particles

    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", 
                      dest="Nsph", type="int", 
                      default = 1000,
                      help="number of SPH particles[%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="mass", type="float", 
                      default = 0.6 | units.MSun,
                      help="stellar mass [%default]")
    result.add_option("-o", unit=units.day**-1,
                      dest="omega", type="float", 
                      default = 24 | units.day**-1,
                      help="stellar rotation [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="age", type="float", 
                      default = 10 | units.Myr,
                      help="stellar age [%default]")
    return result

if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()
    stellar = evolve_star(o.mass, o.age)
    star, core, age = evolve_star_and_convert_to_sph(stellar, o.omega, o.Nsph)
    print "age=", age.in_(units.Myr)
    print star
    print star.mass.sum().in_(units.MSun), core.mass
    filename = 'Hydro_M%2.2dMSun.h5'%o.mass.value_in(units.MSun)
    print filename

    write_set_to_file(star, filename, format='hdf5', append_to_file=False)
