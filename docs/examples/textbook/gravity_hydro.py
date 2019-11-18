"""
   example code for bridging a gravity solver with a hydrodynamics solver
"""
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary

###BOOKLISTSTART1###
class BaseCode:
    def __init__(self, code, particles, eps=0|units.RSun):

        self.local_particles = particles
        m = self.local_particles.mass.sum()
        l = self.local_particles.position.length()
        self.converter = nbody_system.nbody_to_si(m, l)
        self.code = code(self.converter)
        self.code.parameters.epsilon_squared = eps**2

    def evolve_model(self, time):
        self.code.evolve_model(time)
    def copy_to_framework(self):
        self.channel_to_framework.copy()
    def get_gravity_at_point(self, r, x, y, z):
        return self.code.get_gravity_at_point(r, x, y, z)
    def get_potential_at_point(self, r, x, y, z):
        return self.code.get_potential_at_point(r, x, y, z)
    def get_timestep(self):
        return self.code.parameters.timestep
    @property
    def model_time(self):            
        return self.code.model_time
    @property
    def particles(self):
        return self.code.particles
    @property
    def total_energy(self):
        return self.code.kinetic_energy + self.code.potential_energy
    @property
    def stop(self):
        return self.code.stop
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
class Gravity(BaseCode):
    def __init__(self, code, particles, eps=0|units.RSun):
        BaseCode.__init__(self, code, particles, eps)
        self.code.particles.add_particles(self.local_particles)
        self.channel_to_framework \
            = self.code.particles.new_channel_to(self.local_particles)
        self.channel_from_framework \
            = self.local_particles.new_channel_to(self.code.particles)
        self.initial_total_energy = self.total_energy
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
class Hydro(BaseCode):
    def __init__(self, code, particles, eps=0|units.RSun,
                 dt=None, Rbound=None):
        BaseCode.__init__(self, code, particles, eps)
        self.channel_to_framework \
            = self.code.gas_particles.new_channel_to(self.local_particles)
        self.channel_from_framework \
            = self.local_particles.new_channel_to(self.code.gas_particles)
        self.code.gas_particles.add_particles(particles)
        m = self.local_particles.mass.sum()
        l = self.code.gas_particles.position.length()
        if Rbound is None:
            Rbound = 10*l
        self.code.parameters.periodic_box_size = Rbound
        if dt is None:
            dt = 0.01*numpy.sqrt(l**3/(constants.G*m))
        self.code.parameters.timestep = dt/8.
        self.initial_total_energy = self.total_energy
    @property
    def total_energy(self):
        return self.code.kinetic_energy \
            + self.code.potential_energy \
            + self.code.thermal_energy
###BOOKLISTSTOP3###
        
###BOOKLISTSTART4###
def gravity_hydro_bridge(Mprim, Msec, a, ecc, t_end, n_steps,
                         Rgas, Mgas, Ngas):

    stars = new_binary_from_orbital_elements(Mprim, Msec, a, ecc,
                                             G=constants.G)
    eps = 1 | units.RSun
    gravity = Gravity(ph4, stars, eps)

    converter = nbody_system.nbody_to_si(1.0|units.MSun, Rgas)
    ism = new_plummer_gas_model(Ngas, convert_nbody=converter)
    ism.move_to_center()
    ism = ism.select(lambda r: r.length()<2*a,["position"])
    hydro = Hydro(Fi, ism, eps)
    model_time = 0 | units.Myr
    filename = "gravhydro.hdf5"
    write_set_to_file(stars.savepoint(model_time), filename, 'amuse')
    write_set_to_file(ism, filename, 'amuse')

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,))
    gravhydro.timestep = 2*hydro.get_timestep()

    while model_time < t_end:
        orbit = orbital_elements_from_binary(stars, G=constants.G)
        dE_gravity = gravity.initial_total_energy/gravity.total_energy
        dE_hydro = hydro.initial_total_energy/hydro.total_energy
        print "Time:", model_time.in_(units.yr), \
              "ae=", orbit[2].in_(units.AU), orbit[3], \
              "dE=", dE_gravity, dE_hydro
        
        model_time += 10*gravhydro.timestep
        gravhydro.evolve_model(model_time)
        gravity.copy_to_framework()
        hydro.copy_to_framework()
        write_set_to_file(stars.savepoint(model_time), filename, 'amuse')
        write_set_to_file(ism, filename, 'amuse')
        print "P=", model_time.in_(units.yr), gravity.particles.x.in_(units.au)
    gravity.stop()
    hydro.stop()
###BOOKLISTSTOP4###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 1000,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-N", dest="Ngas", type="int", default = 1024,
                      help="number of gas particles [%default]")
    result.add_option("--Mprim", unit=units.MSun,
                      dest="Mprim", type="float", default = 3.2|units.MSun,
                      help="Mass of the primary star [%default]")
    result.add_option("--Msec", unit=units.MSun,
                      dest="Msec", type="float", default = 3.1|units.MSun,
                      help="Mass of the secondary star [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mgas", type="float", default = 1|units.MSun,
                      help="Mass of the gas [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rgas", type="float", default = 10|units.AU,
                      help="Size of the gas distribution [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="a", type="float", default = 1|units.AU,
                      help="initial orbital separation [%default]")
    result.add_option("-e", dest="ecc", type="float", default = 0.6,
                      help="initial orbital eccentricity [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 10000|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    numpy.random.seed(123)
    gravity_hydro_bridge(**o.__dict__)
