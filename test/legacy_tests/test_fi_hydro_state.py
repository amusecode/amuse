import numpy

from amuse.legacy.fi.interface import Fi
from amuse.ext.spherical_model import new_uniform_spherical_particle_distribution
from amuse.support.units import nbody_system as nbody
from amuse.support.units import units
    
def test15():
        print "Testing Fi get_hydro_state_at_point II: uniform sphere"
        number_sph_particles = 10000
        convert_nbody = nbody.nbody_to_si(1.0 | units.kpc, 1.0e9 | units.MSun)
        gas = new_uniform_spherical_particle_distribution(number_sph_particles, 1.0 | units.kpc, 1.0e9 | units.MSun, seed = 1234)
        gas.velocity = [0.0, 0.0, 0.0] | units.m / units.s
        gas.h_smooth = 0.01 | nbody.length
        gas.u = 0.05 | nbody.specific_energy
        instance = Fi(convert_nbody)
        instance.parameters.n_smooth     =   64 | units.none
        instance.parameters.n_smooth_tol = 0.2 | units.none
        instance.gas_particles.add_particles(gas)
#        print instance.gas_particles.h_smooth.as_quantity_in(units.kpc)
        coords = 0.0 | units.kpc
        speeds = [0.0 | units.m / units.s]*3
        hydro_state = instance.get_hydro_state_at_point(coords, coords, coords, *speeds)
        print hydro_state
        expected = [ 3.5469e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     8.6426e-10 | units.kg * units.m**-1 * units.s**-2]
        
        coords2 = 0.1 | units.kpc
        hydro_state = instance.get_hydro_state_at_point(coords2, coords2, coords2, *speeds)
        print hydro_state
        expected = [ 4.1456e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     9.9316e-10 | units.kg * units.m**-1 * units.s**-2]
        
        print ((1.0e10 | units.MSun) / (4.0/3.0 * numpy.pi * (1.0 | units.kpc)**3)).as_quantity_in(units.kg/units.m**3)
        instance.stop()
    
if __name__ == '__main__':
    test15()
