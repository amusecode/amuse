import numpy as np

from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics
from amuse.rfi.core import PythonCodeInterface

try:
    from pynbody.integrator import Integrator
    from pynbody.particles import Particles
    from pynbody.particles.body import Body
    MODULES_MISSING = False
except ImportError:
    MODULES_MISSING = True

"""
MyCodeImplementation is what needs to be adapted to each specific
community code, MyCodeInterface and MyCode do not need to be changed
for standard dynamics codes (except for changing the name).
"""



class PyNbodyImplementation(object):

    def __init__(self):
        self.eta = 0.01
        self.current_time = 0.0
        self.eps2 = 0.0
        self.time_begin = 0.0
        self.integrator_method = "hts"
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
        particles = Particles({"body": num})
        for (i, p) in enumerate(self.particles):
            particles["body"].id[i] = i
            particles["body"].mass[i] = p.mass
            particles["body"].radius[i] = p.radius   # XXX: 'radius' is not used in PyNbody yet.
            particles["body"].eps2[i] = self.eps2/2
            particles["body"].pos[i] = p.pos
            particles["body"].vel[i] = p.vel
        self.integrator = Integrator(self.eta, self.time_begin, particles, method=self.integrator_method)
        return 0

    def synchronize_model(self):
        return 0

    def new_particle(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        body = Body(1)
        body.mass = mass
        body.radius = radius
        body.pos = [x, y, z]
        body.vel = [vx, vy, vz]
        self.particles.append(body)
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
            particle = particles['body']
            particle.mass[i] = mass
            particle.radius[i] = radius
            particle.pos[i] = [x, y, z]
            particle.vel[i] = [vx, vy, vz]
            return 0
        except Exception as exc:
            print str(exc)
            return -1

    def set_mass(self, index_of_the_particle, mass):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            particle.mass[i] = mass
            return 0
        except Exception as exc:
            print str(exc)
            return -1

    def set_radius(self, index_of_the_particle, radius):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            particle.radius[i] = radius
            return 0
        except Exception as exc:
            print str(exc)
            return -1

    def set_position(self, index_of_the_particle, x, y, z):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            particle.pos[i] = [x, y, z]
            return 0
        except:
            return -1

    def set_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            particle.vel[i] = [vx, vy, vz]
            return 0
        except:
            return -1


    def get_state(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            mass.value = particle.mass[i]
            radius.value = particle.radius[i]
            x.value, y.value, z.value = particle.pos[i]
            vx.value, vy.value, vz.value = particle.vel[i]
            return 0
        except:
            return -1

    def get_mass(self, index_of_the_particle, mass):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            mass.value = particle.mass[i]
            return 0
        except:
            return -1

    def get_radius(self, index_of_the_particle, radius):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            radius.value = particle.radius[i]
            return 0
        except:
            return -1

    def get_position(self, index_of_the_particle, x, y, z):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            x.value, y.value, z.value = particle.pos[i]
            return 0
        except:
            return -1

    def get_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            i = index_of_the_particle
            particles = self.integrator.particles
            particle = particles['body']
            vx.value, vy.value, vz.value = particle.vel[i]
            return 0
        except:
            return -1


    def get_kinetic_energy(self, kinetic_energy):
        particles = self.integrator.particles
        ke = particles.get_total_kinetic_energy()
        kinetic_energy.value = ke
        return 0

    def get_potential_energy(self, potential_energy):
        particles = self.integrator.particles
        particles.update_phi(particles)
        pe = particles.get_total_potential_energy()
        potential_energy.value = pe
        return 0


    def get_total_mass(self, total_mass):
        particles = self.integrator.particles
        mtot = particles.get_total_mass()
        total_mass.value = mtot
        return 0

    def get_center_of_mass_position(self, x, y, z):
        particles = self.integrator.particles
        rcom = particles.get_center_of_mass_position()
        x.value, y.value, z.value = rcom
        return 0

    def get_center_of_mass_velocity(self, vx, vy, vz):
        particles = self.integrator.particles
        vcom = particles.get_center_of_mass_velocity()
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
        while (self.integrator.current_time < t_end):
            self.integrator.step(t_end)
        self.current_time = self.integrator.current_time
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


    def set_eps2(self, epsilon_squared):
        self.eps2 = epsilon_squared
        return 0

    def get_eps2(self, epsilon_squared):
        epsilon_squared.value = self.eps2
        return 0


    def set_time_begin(self, time_begin):
        self.time_begin = time_begin
        return 0

    def get_time_begin(self, time_begin):
        time_begin.value = self.time_begin
        return 0


    def set_integrator_method(self, integrator_method):
        self.integrator_method = integrator_method
        return 0

    def get_integrator_method(self, integrator_method):
        integrator_method.value = self.integrator_method
        return 0



class PyNbodyInterface(PythonCodeInterface, GravitationalDynamicsInterface):

    def __init__(self, **options):
        PythonCodeInterface.__init__(self, PyNbodyImplementation, **options)


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
    def set_time_begin():
        function = LegacyFunctionSpecification()
        function.addParameter('time_begin', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time_begin():
        function = LegacyFunctionSpecification()
        function.addParameter('time_begin', dtype='float64', direction=function.OUT)
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



class PyNbody(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **options):
        nbody_interface = PyNbodyInterface(**options)

        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
            **options
        )


    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eta",
            "set_eta",
            "timestep_parameter",
            "timestep parameter",
            default_value = 0.01 | units.none
        )

        object.add_method_parameter(
            "get_time",
            "set_time",
            "time",
            "current simulation time",
            default_value = 0.0 | nbody_system.time
        )

        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )


    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            "get_eta",
            (),
            (units.none, object.ERROR_CODE,)
        )

        object.add_method(
            "set_eta",
            (units.none,),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_time",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )

        object.add_method(
            "set_time",
            (nbody_system.time,),
            (object.ERROR_CODE,)
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

