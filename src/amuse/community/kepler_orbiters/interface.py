from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics, GravityFieldCode
from amuse.support.options import option
from amuse.units import units
from amuse.datamodel import Particles
import os.path

class KeplerInterface(CodeInterface,
                       GravitationalDynamicsInterface):
    """
    kepler_orbiters is a code to solve Kepler's equations for number of
    two-body problems that share one of the bodies. It is Suited for problems
    with a dominant central mass (e.g., planetary systems, debris disks).
    
    The orbits are initialized from mass, position, and velocity.
    The shared body is represented by particle 'central_particle' and 
    defined as:
      instance.central_particle.add_particles(central_particle)
    and the other particles are 'orbiters':
      instance.orbiters.add_particles(orbiters)
    The postition and velocity vectors of the orbiters are integrated, while 
    these of the central particle stays fixed.
    
    The Kepler's equations are solved using universal variables with two 
    independent solvers provided and set by parameter 'method':
    * solver from the SAKURA code (Goncalves Ferrari et al., 2014), default,
      instance.parameters.method = 0
    * solver from the HUAYNO code (Pelupessy et al., 2012)
      instance.parameters.method = 1
    
    See the example kepler_orbiters_planets_around_sun.py.
    
    The relevant references are:
    .. [#] Goncalves Ferrari, Boekholt, Portegies Zwart; 2014 MNRAS, 440, 719 \
               (SAKURA, method 0 -- default) [2014MNRAS.440..719G]
    .. [#] Pelupessy, Janes, Portegies Zwart; 2012, New Astronomy, 17, 711 \
               (HUAYNO, method 1) [2012NewA...17..711P]
    """

    include_headers = ['interface.h']

    def __init__(self, **options):
        CodeInterface.__init__(self,
                               name_of_the_worker = "keplerorbiters_worker",
                               **options)

    @legacy_function
    def get_method():
        function = LegacyFunctionSpecification()
        function.addParameter('method', dtype='i', direction=function.OUT, default=0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_method():
        function = LegacyFunctionSpecification()
        function.addParameter('method', dtype='i', direction=function.IN, default=0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_central_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN, default=0)
        function.addParameter('mass', dtype='float64', direction=function.OUT,
                              unit = nbody_system.mass)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.IN,
                              unit = nbody_system.mass)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_central_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN, default=0)
        function.addParameter('radius', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        function.addParameter('radius', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_central_pos():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN, default=0)
        function.addParameter('x', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.addParameter('y', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.addParameter('z', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_pos():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('y', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('z', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_central_vel():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN, default=0)
        function.addParameter('vx', dtype='float64', direction=function.OUT,
                              unit = nbody_system.speed)
        function.addParameter('vy', dtype='float64', direction=function.OUT,
                              unit = nbody_system.speed)
        function.addParameter('vz', dtype='float64', direction=function.OUT,
                              unit = nbody_system.speed)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_vel():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vy', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vz', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_central_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.OUT)
        function.addParameter('mass', dtype='float64', direction=function.IN,
                              unit = nbody_system.mass)
        function.addParameter('x', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('y', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('z', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('vx', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vy', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vz', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('radius', dtype='float64', direction=function.IN,
                              unit = nbody_system.length,default=0.)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_central_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        function.result_type = 'int32'
        return function


    def get_gravity_at_point(self,radius,x,y,z):
        mass,err=self.get_central_mass()
        xc,yc,zc,err=self.get_central_pos()
        dr2=((x-xc)**2+(y-yc)**2+(z-zc)**2+radius**2)
        dr=dr2**0.5
        ax=-mass*(x-xc)/(dr2*dr)
        ay=-mass*(y-yc)/(dr2*dr)
        az=-mass*(z-zc)/(dr2*dr)
        return ax,ay,az

    def get_potential_at_point(self,radius,x,y,z):
        mass,err=self.get_central_mass()
        xc,yc,zc,err=self.get_central_pos()
        dr2=((x-xc)**2+(y-yc)**2+(z-zc)**2+radius**2)
        dr=dr2**0.5
        phi=-mass/dr
        return phi

class Kepler(GravitationalDynamics, GravityFieldCode):

    def __init__(self, unit_converter = None,  **options):
        self.unit_converter = unit_converter

        CommonCode.__init__(self,
                               KeplerInterface(**options),
                               **options)

    def define_converter(self, handler):
        if not self.unit_converter is None:
            handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())


    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_gravity_at_point',
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                nbody_system.acceleration,
                nbody_system.acceleration,
                nbody_system.acceleration,
            )
        )
        handler.add_method(
            'get_potential_at_point',
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                nbody_system.potential,
            )
        )


    def define_particle_sets(self, handler):
        handler.define_super_set('particles', ['central_particle','orbiters'],
            index_to_default_set = 1)

        handler.define_set('central_particle', 'index_of_the_particle')
        handler.set_new('central_particle', 'new_central_particle')
        handler.set_delete('central_particle', 'delete_central_particle')
        handler.add_setter('central_particle', 'set_central_mass')
        handler.add_getter('central_particle', 'get_central_mass')
        handler.add_setter('central_particle', 'set_central_radius')
        handler.add_getter('central_particle', 'get_central_radius')
        handler.add_setter('central_particle', 'set_central_pos')
        handler.add_getter('central_particle', 'get_central_pos')
        handler.add_setter('central_particle', 'set_central_vel')
        handler.add_getter('central_particle', 'get_central_vel')

        handler.define_set('orbiters', 'index_of_the_particle')
        handler.set_new('orbiters', 'new_particle')
        handler.set_delete('orbiters', 'delete_particle')
        handler.add_setter('orbiters', 'set_state')
        handler.add_getter('orbiters', 'get_state')
        handler.add_setter('orbiters', 'set_mass')
        handler.add_getter('orbiters', 'get_mass')
        handler.add_setter('orbiters', 'set_position')
        handler.add_getter('orbiters', 'get_position')
        handler.add_setter('orbiters', 'set_velocity')
        handler.add_getter('orbiters', 'get_velocity')


    def get_gravity_at_point(self,radius,x,y,z):
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z
        return self.overridden().get_gravity_at_point(radius,xx,yy,zz)

    def get_potential_at_point(self,radius,x,y,z):
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z
        return self.overridden().get_potential_at_point(radius,xx,yy,zz)


Kepler_orbiters = Kepler
