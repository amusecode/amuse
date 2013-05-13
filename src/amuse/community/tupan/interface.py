import numpy as np

from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode
from amuse.rfi.core import PythonCodeInterface

try:
    from tupan.integrator import Integrator
    from tupan.particles.allparticles import System
    from tupan.particles.star import Stars
    MODULES_MISSING = False
except ImportError:
    MODULES_MISSING = True

"""
MyCodeImplementation is what needs to be adapted to each specific
community code, MyCodeInterface and MyCode do not need to be changed
for standard dynamics codes (except for changing the name).
"""



class TupanImplementation(object):

    def __init__(self):
        self.eta = 0.03125
        self.current_time = 0.0
        self.eps2 = 0.0
        self.time_begin = 0.0
        self.integrator_method = "sia.dkd21hcc"
        self.particles = []
        self.particles_initialized = False


    def initialize_code(self):
        return 0

    def cleanup_code(self):
        return 0

    def commit_parameters(self):
        return 0



    def commit_particles(self):
        num = len(self.particles)
        particles = System(nstars=num)
        for (i, p) in enumerate(self.particles):
            particles.stars.id[i] = i
            particles.stars.mass[i] = p.mass
            particles.stars.radius[i] = p.radius   # XXX: 'radius' is not yet used in Tupan.
            particles.stars.eps2[i] = self.eps2/2
            particles.stars.rx[i] = p.rx
            particles.stars.ry[i] = p.ry
            particles.stars.rz[i] = p.rz
            particles.stars.vx[i] = p.vx
            particles.stars.vy[i] = p.vy
            particles.stars.vz[i] = p.vz
        self.integrator = Integrator(self.eta, self.time_begin, particles, method=self.integrator_method)
        return 0

    def synchronize_model(self):
        return 0

    def new_particle(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        star = Stars(1)
        star.mass = mass
        star.radius = radius
        star.rx = x
        star.ry = y
        star.rz = z
        star.vx = vx
        star.vy = vy
        star.vz = vz
        self.particles.append(star)
        index_of_the_particle.value = len(self.particles)-1
        return 0


#    def new_black_hole(self, index_of_the_particle, mass, x, y, z, vx, vy, vz, spin):
#        self.particle_buffer.append(
#            {
#            'mass': mass,
#            'radius' : radius,
#            'x' : x,
#            'y' : y,
#            'z' : z,
#            'vx' : vx,
#            'vy' : vy,
#            'vz' : vz,
#            }
#        )
#        index_of_the_particle.value = len(self.particle_buffer)-1
#        return 0


    def set_state(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            p.mass[i] = mass
            p.radius[i] = radius
            p.rx[i] = x
            p.ry[i] = y
            p.rz[i] = z
            p.vx[i] = vx
            p.vy[i] = vy
            p.vz[i] = vz
            return 0
        except Exception as exc:
            print str(exc)
            return -1

    def set_mass(self, index_of_the_particle, mass):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            p.mass[i] = mass
            return 0
        except Exception as exc:
            print str(exc)
            return -1

    def set_radius(self, index_of_the_particle, radius):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            p.radius[i] = radius
            return 0
        except Exception as exc:
            print str(exc)
            return -1

    def set_position(self, index_of_the_particle, x, y, z):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            p.rx[i] = x
            p.ry[i] = y
            p.rz[i] = z
            return 0
        except:
            return -1

    def set_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            p.vx[i] = vx
            p.vy[i] = vy
            p.vz[i] = vz
            return 0
        except:
            return -1


    def get_state(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            mass.value = p.mass[i]
            radius.value = p.radius[i]
            x.value, y.value, z.value = p.rx[i], p.ry[i], p.rz[i]
            vx.value, vy.value, vz.value = p.vx[i], p.vy[i], p.vz[i]
            return 0
        except:
            return -1

    def get_mass(self, index_of_the_particle, mass):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            mass.value = p.mass[i]
            return 0
        except:
            return -1

    def get_radius(self, index_of_the_particle, radius):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            radius.value = p.radius[i]
            return 0
        except:
            return -1

    def get_position(self, index_of_the_particle, x, y, z):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            x.value, y.value, z.value = p.rx[i], p.ry[i], p.rz[i]
            return 0
        except:
            return -1

    def get_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            p = particles.stars
            vx.value, vy.value, vz.value = p.vx[i], p.vy[i], p.vz[i]
            return 0
        except:
            return -1


    def get_kinetic_energy(self, kinetic_energy):
        particles = self.integrator.particles
        ke = particles.kinetic_energy
        kinetic_energy.value = ke
        return 0

    def get_potential_energy(self, potential_energy):
        particles = self.integrator.particles
        pe = particles.potential_energy
        potential_energy.value = pe
        return 0


    def get_total_mass(self, total_mass):
        particles = self.integrator.particles
        mtot = particles.total_mass
        total_mass.value = mtot
        return 0

    def get_center_of_mass_position(self, x, y, z):
        particles = self.integrator.particles
        rcom = particles.rcom
        x.value, y.value, z.value = rcom
        return 0

    def get_center_of_mass_velocity(self, vx, vy, vz):
        particles = self.integrator.particles
        vcom = particles.vcom
        vx.value, vy.value, vz.value = vcom
        return 0


    def get_gravity_at_point(self, eps, x, y, z, ax, ay, az, length):
        ax.value = 0.0
        ay.value = 0.0
        az.value = 0.0
        return -2 # Not implemented

    def get_potential_at_point(self, eps, x, y, z, phi, length):
        phi.value = 0.0
        return -2 # Not implemented


    def evolve_model(self, t_end):
        while (abs(self.integrator.time) < t_end):
            self.integrator.evolve_step(t_end)
        self.integrator.finalize(t_end)
        self.current_time = self.integrator.time
        return 0


    def set_eta(self, eta):
        self.eta = eta
        return 0

    def get_eta(self, eta):
        eta.value = self.eta
        return 0


    def set_time(self, time):
        self.current_time = time
        return 0

    def get_time(self, time):
        time.value = self.current_time
        return 0

    def set_time(self, time):
        self.current_time = time
        return 0


    def set_eps2(self, epsilon_squared):
        self.eps2 = epsilon_squared
        return 0

    def get_eps2(self, epsilon_squared):
        epsilon_squared.value = self.eps2
        return 0


    def set_begin_time(self, time_begin):
        self.time_begin = time_begin
        return 0

    def get_begin_time(self, time_begin):
        time_begin.value = self.time_begin
        return 0


    def set_integrator_method(self, integrator_method):
        self.integrator_method = integrator_method
        return 0

    def get_integrator_method(self, integrator_method):
        integrator_method.value = self.integrator_method
        return 0



class TupanInterface(PythonCodeInterface, GravitationalDynamicsInterface, SinglePointGravityFieldInterface):

    def __init__(self, **options):
        PythonCodeInterface.__init__(self, TupanImplementation, 'tupan_worker', **options)


    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle", default = 0)
        function.result_type = 'int32'
        return function


#    @legacy_function
#    def new_black_hole():
#        function = LegacyFunctionSpecification()
#        function.can_handle_array = True
#        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
#        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
#        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
#        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
#        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
#        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
#        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
#        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
#        function.addParameter('spin', dtype='float64', direction=function.IN, description = "The spin of the particle", default = 0)
#        function.result_type = 'int32'
#        return function


    @legacy_function
    def set_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_integrator_method():
        function = LegacyFunctionSpecification()
        function.addParameter('integrator_method', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_integrator_method():
        function = LegacyFunctionSpecification()
        function.addParameter('integrator_method', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function



class Tupan(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody = None, **options):
        nbody_interface = TupanInterface(**options)

        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
            **options
        )

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        GravityFieldCode.define_state(self, object)


    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eta",
            "set_eta",
            "timestep_parameter",
            "timestep parameter",
            default_value = 0.01
        )

        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )

        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            "get_eta",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            "set_eta",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_time",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )

        object.add_method(
            "set_begin_time",
            (nbody_system.time,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_begin_time",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )

        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )

        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )

