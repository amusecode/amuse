from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics

class MercuryInterface(CodeInterface):
    def __init__(self, **args):
        CodeInterface.__init__(self, name_of_the_worker = 'mercury_worker',**args)

    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function   
    def cleanup_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def commit_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function          
    def delete_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function    

    @legacy_function    
    def new_orbiter():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','density','x','y','z','vx','vy','vz','sx','sy','sz','celimit']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_orbiter_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz','sx','sy','sz','celimit']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_orbiter_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz','sx','sy','sz','celimit']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_density():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('density', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_density():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('density', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_mass():
        function = LegacyFunctionSpecification()   
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_mass():
        function = LegacyFunctionSpecification()   
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_radius():
        function = LegacyFunctionSpecification()   
        function.addParameter('radius', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_radius():
        function = LegacyFunctionSpecification()   
        function.addParameter('radius', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_oblateness():
        function = LegacyFunctionSpecification()   
        for x in ['j2','j4','j6']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_oblateness():
        function = LegacyFunctionSpecification()   
        for x in ['j2','j4','j6']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_spin():
        function = LegacyFunctionSpecification()   
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_spin():
        function = LegacyFunctionSpecification()   
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_number_of_orbiters():
        function = LegacyFunctionSpecification()   
        function.addParameter('norbiters', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_energy_deviation():
        function = LegacyFunctionSpecification()   
        function.addParameter('delta_e', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('ek', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_potential_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('ep', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_total_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('e_tot', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_total_angular_momentum():
        function = LegacyFunctionSpecification()   
        function.addParameter('am_tot', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

class MercuryDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_doc+"\n\n"+instance.parameters.__doc__

class Mercury(GravitationalDynamics):

    __doc__ = MercuryDoc()

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = MercuryInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )
        
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_time",
            None,
            "time",
            "current simulation time", 
            nbody_system.time, 
            0.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_central_mass",
            "set_central_mass",
            "central_mass",
            "mass of the central object", 
            nbody_system.mass, 
            1.0 | nbody_system.mass
        )
        object.add_method_parameter(
            "get_central_radius",
            "set_central_radius",
            "central_radius",
            "radius of the central object", 
            nbody_system.length, 
            1.0 | nbody_system.length
        )

        self.stopping_conditions.define_parameters(object)

    def define_particle_sets(self, object):
        print object
        object.define_super_set('particles', ['orbiters'],
                                 index_to_default_set = 0)
        object.define_set('orbiters', 'index_of_the_particle')
        object.set_new('orbiters', 'new_orbiter')
        object.set_delete('orbiters', 'delete_particle')
        object.add_setter('orbiters', 'set_orbiter_state')
        object.add_getter('orbiters', 'get_orbiter_state')
        object.add_setter('orbiters', 'set_mass')
        object.add_getter('orbiters', 'get_mass')
        object.add_setter('orbiters', 'set_density')
        object.add_getter('orbiters', 'get_density')
        object.add_setter('orbiters', 'set_position')
        object.add_getter('orbiters', 'get_position')
        object.add_setter('orbiters', 'set_velocity')
        object.add_getter('orbiters', 'get_velocity')

        #GravitationalDynamics.define_particle_sets(self, object)
        #self.stopping_conditions.define_particle_set(object, 'particles')

    def define_methods(self, object):
        #GravitationalDynamics.define_methods(self, object)
        object.add_method(
            'new_orbiter',
            (
                units.g,
                units.g/units.cm**3,
                units.cm,
                units.cm,
                units.cm,
                units.cm/units.s,
                units.cm/units.s,
                units.cm/units.s,
                units.s**-1,
                units.s**-1,
                units.s**-1,
                units.none,
            ),
            (
                units.none,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_mass",
            (
                object.NO_UNIT,
                generic_unit_system.mass,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_mass",
            (
                object.NO_UNIT,
            ),
            (
                generic_unit_system.mass,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_position",
            (
                object.NO_UNIT,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_position",
            (
                object.NO_UNIT,
            ),
            (
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_velocity",
            (
                object.INDEX,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_velocity",
            (
                object.INDEX,
            ),
            (
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_spin",
            (
                object.INDEX,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_spin",
            (
                object.INDEX,
            ),
            (
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            'set_density',
            (
                object.NO_UNIT,
                units.g/units.cm**3
            ),
            (
                object.ERROR_CODE
            )
        )    
        object.add_method(
            'get_density',
            (
                object.NO_UNIT,
            ),
            (
                units.g/units.cm**3,
                object.ERROR_CODE
            )
        )    

        object.add_method(
            'get_central_oblateness',
            (
                object.NO_UNIT,        
            ),
            (
                units.none,
                units.none,
                units.none,
                object.ERROR_CODE
            )
        )
        object.add_method(
            'set_central_oblateness',
            (
                units.none,
                units.none,
                units.none,
            ),
            (
                object.ERROR_CODE,
            )
        )
        
        self.stopping_conditions.define_methods(object)



"""
"""
