from __future__ import print_function

from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode
from amuse.rfi.core import PythonCodeInterface

try:
    from tupan.integrator import Integrator
    from tupan.particles.allparticles import ParticleSystem
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
        self.integrator_method = "sia21h.dkd"
        self.pn_order = 0
        self.clight = None
        self.particles = []
        self.particles_initialized = False

    def initialize_code(self):
        return 0

    def cleanup_code(self):
        return 0

    def commit_parameters(self):
        if not self.integrator_method in Integrator.PROVIDED_METHODS:
            msg = "Unknown integrator: {0}. Provided methods are: {1}."
            print(msg.format(self.integrator_method,
                             Integrator.PROVIDED_METHODS))
            return -1
        if self.pn_order > 0 and self.clight is None:
            print("'clight' is None. Please set the speed of light "
                  "parameter 'clight' when using 'pn_order' > 0.")
            return -1
        return 0

    def commit_particles(self):
        num = len(self.particles)
        ps = ParticleSystem(nstars=num)
        for (i, p) in enumerate(self.particles):
            ps.id[i] = i
            ps.mass[i] = p.mass
            ps.radius[i] = p.radius   # XXX: 'radius' is not yet used in Tupan.
            ps.eps2[i] = self.eps2/2
            ps.rx[i] = p.rx
            ps.ry[i] = p.ry
            ps.rz[i] = p.rz
            ps.vx[i] = p.vx
            ps.vy[i] = p.vy
            ps.vz[i] = p.vz
        self.integrator = Integrator(self.eta,
                                     self.time_begin,
                                     ps,
                                     method=self.integrator_method,
                                     pn_order=self.pn_order,
                                     clight=self.clight)
        return 0

    def synchronize_model(self):
        return 0

    def new_particle(self, index_of_the_particle,
                     mass, radius, x, y, z, vx, vy, vz):
        ps = ParticleSystem(nstars=1)
        ps.mass[0] = mass
        ps.radius[0] = radius
        ps.rx[0] = x
        ps.ry[0] = y
        ps.rz[0] = z
        ps.vx[0] = vx
        ps.vy[0] = vy
        ps.vz[0] = vz
        self.particles.append(ps)
        index_of_the_particle.value = len(self.particles)-1
        return 0

    def set_state(self, index_of_the_particle,
                  mass, radius, x, y, z, vx, vy, vz):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            ps.mass[i] = mass
            ps.radius[i] = radius
            ps.rx[i] = x
            ps.ry[i] = y
            ps.rz[i] = z
            ps.vx[i] = vx
            ps.vy[i] = vy
            ps.vz[i] = vz
            return 0
        except Exception as exc:
            print(str(exc))
            return -1

    def set_mass(self, index_of_the_particle, mass):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            ps.mass[i] = mass
            return 0
        except Exception as exc:
            print(str(exc))
            return -1

    def set_radius(self, index_of_the_particle, radius):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            ps.radius[i] = radius
            return 0
        except Exception as exc:
            print(str(exc))
            return -1

    def set_position(self, index_of_the_particle, x, y, z):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            ps.rx[i] = x
            ps.ry[i] = y
            ps.rz[i] = z
            return 0
        except:
            return -1

    def set_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            ps.vx[i] = vx
            ps.vy[i] = vy
            ps.vz[i] = vz
            return 0
        except:
            return -1

    def get_state(self, index_of_the_particle,
                  mass, radius, x, y, z, vx, vy, vz):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            mass.value = ps.mass[i]
            radius.value = ps.radius[i]
            x.value, y.value, z.value = ps.rx[i], ps.ry[i], ps.rz[i]
            vx.value, vy.value, vz.value = ps.vx[i], ps.vy[i], ps.vz[i]
            return 0
        except:
            return -1

    def get_mass(self, index_of_the_particle, mass):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            mass.value = ps.mass[i]
            return 0
        except:
            return -1

    def get_radius(self, index_of_the_particle, radius):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            radius.value = ps.radius[i]
            return 0
        except:
            return -1

    def get_position(self, index_of_the_particle, x, y, z):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            x.value, y.value, z.value = ps.rx[i], ps.ry[i], ps.rz[i]
            return 0
        except:
            return -1

    def get_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            i = index_of_the_particle
            ps = self.integrator.particle_system
            vx.value, vy.value, vz.value = ps.vx[i], ps.vy[i], ps.vz[i]
            return 0
        except:
            return -1

    def get_kinetic_energy(self, kinetic_energy):
        ps = self.integrator.particle_system
        ke = ps.kinetic_energy
        kinetic_energy.value = ke
        return 0

    def get_potential_energy(self, potential_energy):
        ps = self.integrator.particle_system
        pe = ps.potential_energy
        potential_energy.value = pe
        return 0

    def get_total_mass(self, total_mass):
        ps = self.integrator.particle_system
        mtot = ps.total_mass
        total_mass.value = mtot
        return 0

    def get_center_of_mass_position(self, x, y, z):
        ps = self.integrator.particle_system
        com_r = ps.com_r
        x.value, y.value, z.value = com_r
        return 0

    def get_center_of_mass_velocity(self, vx, vy, vz):
        ps = self.integrator.particle_system
        com_v = ps.com_v
        vx.value, vy.value, vz.value = com_v
        return 0

    def get_gravity_at_point(self, eps, x, y, z, ax, ay, az, length):
        ax.value = 0.0
        ay.value = 0.0
        az.value = 0.0
        return -2  # Not implemented

    def get_potential_at_point(self, eps, x, y, z, phi, length):
        phi.value = 0.0
        return -2  # Not implemented

    def evolve_model(self, t_end):
        while (abs(self.integrator.time) < abs(t_end)):
            self.integrator.evolve_step(t_end)
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

    def set_pn_order(self, pn_order):
        self.pn_order = pn_order
        return 0

    def get_pn_order(self, pn_order):
        pn_order.value = self.pn_order
        return 0

    def set_clight(self, clight):
        self.clight = clight
        return 0

    def get_clight(self, clight):
        clight.value = self.clight
        return 0


class TupanInterface(PythonCodeInterface,
                     GravitationalDynamicsInterface,
                     SinglePointGravityFieldInterface):

    def __init__(self, **options):
        PythonCodeInterface.__init__(
            self,
            TupanImplementation,
            'tupan_worker',
            **options)

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle',
            dtype='int32',
            direction=function.OUT)
        function.addParameter(
            'mass',
            dtype='float64',
            direction=function.IN,
            description="The mass of the particle")
        function.addParameter(
            'x',
            dtype='float64',
            direction=function.IN,
            description="The initial position vector of the particle")
        function.addParameter(
            'y',
            dtype='float64',
            direction=function.IN,
            description="The initial position vector of the particle")
        function.addParameter(
            'z',
            dtype='float64',
            direction=function.IN,
            description="The initial position vector of the particle")
        function.addParameter(
            'vx',
            dtype='float64',
            direction=function.IN,
            description="The initial velocity vector of the particle")
        function.addParameter(
            'vy',
            dtype='float64',
            direction=function.IN,
            description="The initial velocity vector of the particle")
        function.addParameter(
            'vz',
            dtype='float64',
            direction=function.IN,
            description="The initial velocity vector of the particle")
        function.addParameter(
            'radius',
            dtype='float64',
            direction=function.IN,
            description="The radius of the particle",
            default=0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eta():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'eta',
            dtype='float64',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'eta',
            dtype='float64',
            direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'time',
            dtype='float64',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'time',
            dtype='float64',
            direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'epsilon_squared',
            dtype='float64',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'epsilon_squared',
            dtype='float64',
            direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_integrator_method():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'integrator_method',
            dtype='string',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_integrator_method():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'integrator_method',
            dtype='string',
            direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_pn_order():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'pn_order',
            dtype='int32',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_pn_order():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'pn_order',
            dtype='int32',
            direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_clight():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'clight',
            dtype='float64',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_clight():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'clight',
            dtype='float64',
            direction=function.OUT)
        function.result_type = 'int32'
        return function


class Tupan(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody=None, **options):
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
            default_value=0.01
        )

        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value=0.0 | nbody_system.length * nbody_system.length
        )

        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value=0.0 | nbody_system.time
        )

        object.add_method_parameter(
            "get_integrator_method",
            "set_integrator_method",
            "integrator_method",
            "The method to use to integrate the evolution of the system",
            default_value="sia21h.dkd"
        )

        object.add_method_parameter(
            "get_pn_order",
            "set_pn_order",
            "pn_order",
            "Order of the Post-Newtonian corrections \
                (choices: [0, 2, 4, 5, 6, 7])",
            default_value=0
        )

        object.add_method_parameter(
            "get_clight",
            "set_clight",
            "clight",
            "Speed of light to use in post-Newtonian corrections",
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

        object.add_method(
            "set_clight",
            (nbody_system.speed,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_clight",
            (),
            (nbody_system.speed, object.ERROR_CODE,)
        )


### end of file ###
