try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
import numpy

from amuse.support.data.core import Particle, Grid
from amuse.support.units import units, constants
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.community.evtwin.interface import EVtwin
from amuse.community.fi.interface import Fi
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH


def create_particles():
    star =  Particle()
    star.mass = 3.0 | units.MSun
    
    stellar_evolution = EVtwin()
    stellar_evolution.initialize_code()
    stellar_evolution.commit_parameters() 
    se_star = stellar_evolution.particles.add_particle(star)
    stellar_evolution.commit_particles()
    
    print "Evolving", star.mass, "star with", stellar_evolution.__class__.__name__, "up to", 100 | units.Myr
    stellar_evolution.evolve_model(100 | units.Myr)
    
    print "Creating SPH particles from the (1D) stellar evolution model"
    sph_particles = convert_stellar_model_to_SPH(
        se_star, 
        1000
    ).gas_particles
    stellar_evolution.stop()
    return sph_particles


def hydro_plot(sph_particles):
    print "Creating figure of the 'hydro' star (RGB = log density, log speed, log internal energy)"
    shape = (300, 300, 1)
    size = shape[0] * shape[1]
    grid = Grid.create(shape, (4*max(sph_particles.x), 4*max(sph_particles.x), 0 | units.m))
    grid.x -= 2*max(sph_particles.x)
    grid.y -= 2*max(sph_particles.x)
    speed = ([0.0] * size) | units.m / units.s
    
    unit_converter = ConvertBetweenGenericAndSiUnits(constants.G, max(sph_particles.x), sph_particles.mass.sum())
    hydro_code = Fi(unit_converter)
    hydro_code.parameters.n_smooth = 96 | units.none
    hydro_code.gas_particles.add_particles(sph_particles)
    rho, rhovx, rhovy, rhovz, rhoe = hydro_code.get_hydro_state_at_point(grid.x.reshape(size), 
        grid.y.reshape(size), grid.z.reshape(size), speed, speed, speed)
    
    min_rho = min(hydro_code.gas_particles.rho)
    max_rho = max(hydro_code.gas_particles.rho)
    hydro_code.stop()
    min_v =  10.0 | units.km / units.s
    max_v = 1000.0 | units.km / units.s
    min_E = min(sph_particles.u)*3
    max_E = max(sph_particles.u)
    
    v_sqr = (rhovx**2 + rhovy**2 + rhovz**2) / rho**2
    E = rhoe / rho
    log_v = numpy.log((v_sqr / min_v**2).value_in(units.none)) / numpy.log((max_v**2 / min_v**2).value_in(units.none))
    log_rho = numpy.log((rho / min_rho).value_in(units.none)) / numpy.log((max_rho / min_rho).value_in(units.none))
    log_E = numpy.log((E / min_E).value_in(units.none)) / numpy.log((max_E / min_E).value_in(units.none))
    
    red   = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_rho)).reshape(shape)
    green = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_v)).reshape(shape)
    blue  = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_E)).reshape(shape)
    alpha = numpy.minimum(numpy.ones_like(log_v), numpy.maximum(numpy.zeros_like(log_v), 
        numpy.log((rho / (10*min_rho)).value_in(units.none)))).reshape(shape)
    
    pyplot.figure(figsize = (6, 6), dpi = 50)
    im = pyplot.figimage(numpy.concatenate((red, green, blue, alpha), axis = 2), origin='lower')
    pyplot.show()


if __name__ in ["__main__", "__plot__"]:
    sph_particles = create_particles()
    if HAS_MATPLOTLIB:
        hydro_plot(sph_particles)
    else:
        print sph_particles
