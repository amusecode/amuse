from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

class MercuryInterface(CodeInterface, LiteratureRefs):
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
    def commit_parameters():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function
    def recommit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
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
        for x in ['mass','density','x','y','z','vx','vy','vz','sx','sy','sz','celimit']:
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
    def get_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the position from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current x position of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current y position of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the position for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new x position of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new y position of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the position from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current x position of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current y position of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the position for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new x position of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new y position of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
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
    def get_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
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
    def get_celimit():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('celimit', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_celimit():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('celimit', dtype='d', direction=function.IN)
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
    def get_spin():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['sx','sy','sz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_spin():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['sx','sy','sz']:
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
        object.define_set('particles', 'id')
        object.set_new('particles', 'new_orbiter')
        object.set_delete('particles', 'delete_particle')
        object.add_setter('particles', 'set_mass')
        object.add_getter('particles', 'get_mass')
        object.add_setter('particles', 'set_density')
        object.add_getter('particles', 'get_density')
        object.add_setter('particles', 'set_position')
        object.add_getter('particles', 'get_position')
        object.add_setter('particles', 'set_velocity')
        object.add_getter('particles', 'get_velocity')
        object.add_setter('particles', 'set_spin')
        object.add_getter('particles', 'get_spin')
        object.add_setter('particles', 'set_celimit')
        object.add_getter('particles', 'get_celimit')

        #GravitationalDynamics.define_particle_sets(self, object)
        #self.stopping_conditions.define_particle_set(object, 'particles')

    def define_methods(self, object):
        #GravitationalDynamics.define_methods(self, object)
        object.add_method(
            'new_orbiter',
            (
                units.MSun,
                units.g/units.cm**3,
                units.AU,
                units.AU,
                units.AU,
                units.AU/units.day,
                units.cm/units.day,
                units.cm/units.day,
                units.day,
                units.day,
                units.day,
                units.none
            ),
            (
                object.INDEX, 
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_mass",
            (
                object.INDEX,
                units.MSun,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_mass",
            (
                object.INDEX,
            ),
            (
                units.MSun,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_celimit",
            (
                object.INDEX,
                units.none,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_celimit",
            (
                object.INDEX,
            ),
            (
                units.none,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "set_position",
            (
                object.INDEX,
                units.AU,
                units.AU,
                units.AU,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_position",
            (
                object.INDEX,
            ),
            (
                units.AU,
                units.AU,
                units.AU,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_velocity",
            (
                object.INDEX,
                units.AUd,
                units.AUd,
                units.AUd,
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
                units.AUd,
                units.AUd,
                units.AUd,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_spin",
            (
                object.INDEX,
                units.day,
                units.day,
                units.day,
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
                units.day,
                units.day,
                units.day,
                object.ERROR_CODE
            )
        )
        object.add_method(
            'set_density',
            (
                object.INDEX,
                units.g/units.cm**3
            ),
            (
                object.ERROR_CODE
            )
        )    
        object.add_method(
            'get_density',
            (
                object.INDEX,
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
