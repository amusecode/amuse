from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from numpy import pi
from amuse.units import constants

class MikkolaInterface(CodeInterface,
                       GravitationalDynamicsInterface,
                       StoppingConditionInterface):
    
    use_modules = ['Mikkola', 'StoppingConditions']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="mikkola_worker", **keyword_arguments)
    
    @legacy_function
    def set_time_step():
        """
        Set the model timestep.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "The current model timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time step was retrieved
        -1 - ERROR
            The code does not have support for querying the time
        """
        return function
        
    @legacy_function
    def get_lightspeed():
        """
        Get the lightspeed value, the lightspeed scales the units (like G=1)
        and limits the valid velocity terms
        """
        function = LegacyFunctionSpecification()
        function.addParameter('lightspeed', dtype='float64',
                              direction=function.OUT,
            description = "value for the lightspeed")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_lightspeed():
        """
        Set the lightspeed value, the lightspeed scales the units (like G=1)
        and limits the valid velocity terms
        """
        function = LegacyFunctionSpecification()
        function.addParameter('lightspeed', dtype='float64',
                              direction=function.IN,
            description = "value for the lightspeed")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
        
    @legacy_function
    def get_tolerance():
        """
        Retrieve the accurancy parameter for the evolve
        """
        function = LegacyFunctionSpecification()
        function.addParameter('tolerance', dtype='float64',
                              direction=function.OUT,
            description = "tolerance for the evolve")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_tolerance():
        """
        Set the accurancy parameter for the evolve
        """
        function = LegacyFunctionSpecification()
        function.addParameter('tolerance', dtype='float64',
                              direction=function.IN,
            description = "tolerance for the evolve")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_maximum_number_of_particles():
        """
        Retrieve the maximum number of particles that can be evolved with this code
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'maximum_number_of_particles', 
            dtype='int32',
            direction=function.OUT
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_maximum_number_of_particles():
        """
        Change the maximum number of particles that can be evolved with this code
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'maximum_number_of_particles', 
            dtype='int32',
            direction=function.IN
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_evolve_to_exact_time():
        """
        Evolve to model to the exact time given in the evolve_model call (can be slower)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'evolve_to_exact_time', 
            dtype='bool',
            direction=function.OUT
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_evolve_to_exact_time():
        """
        Evolve to model to the exact time given in the evolve_model call (can be slower)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'evolve_to_exact_time', 
            dtype='bool',
            direction=function.IN
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_radiated_gravitational_energy():
        """
        Retrieve the current radiated gravitational energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('radiated_gravitational_energy', dtype='float64', direction=function.OUT,
            description = "The energy radiated by gravitational waves")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the energy was set
        -1 - ERROR
            Energy could not be provided
        """
        return function
    
    
    @legacy_function
    def get_total_energy():
        """
        Retrieve the current total energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('energy', dtype='float64', direction=function.OUT,
            description = "The energy radiated by gravitational waves")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the energy was set
        -1 - ERROR
            Energy could not be provided
        """
        return function
        
    @legacy_function
    def get_number_of_particles_added():
        """
        Return the number of particles added during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_id_of_added_particle():
        """
        Return the id of the new particle in the code
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_add', dtype='int32',
                              direction=function.IN, 
                 description = 'index in the added particles list (0-n)')
        function.addParameter('index_of_particle', dtype='int32',
                              direction=function.OUT)
        function.can_handle_array = True
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_children_of_particle():
        """
        Return the number of particles added or deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32',
                              direction=function.IN, 
                 description = 'index of the parent particle',
                 unit = INDEX)
        function.addParameter('child1', dtype='int32', direction=function.OUT,
                description = 'index of the first child particle, -1 if none',
                unit = LINK('particles') )
        function.addParameter('child2', dtype='int32', direction=function.OUT,
                unit = LINK('particles'))
        function.can_handle_array = True
        function.result_type = 'int32'
        return function
        
class Mikkola(GravitationalDynamics):

    def __init__(self, convert_nbody=None, **options):
      
        self.stopping_conditions = StoppingConditions(self)
        GravitationalDynamics.__init__(
            self, 
            MikkolaInterface(**options),
            convert_nbody,
            **options
        )
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        
        if not self.unit_converter is None:
            value=self.unit_converter.to_nbody(constants.c)
            self.parameters._original.lightspeed = value

        return result
        
    def define_parameters(self, handler):
        #~ GravitationalDynamics.define_parameters(self, handler)
        
        handler.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "initial timestep for iteration", 
            default_value = 1.0 | nbody_system.time
        )
        handler.add_method_parameter(
            "get_lightspeed", 
            "set_lightspeed",
            "lightspeed", 
            "lightspeed used in the code", 
            default_value = 1.0 | nbody_system.length / nbody_system.time
        )
        handler.add_method_parameter(
            "get_tolerance", 
            "set_tolerance",
            "tolerance", 
            "tolerance used in the code", 
            default_value = 1e-13
        )
        handler.add_method_parameter(
            "get_maximum_number_of_particles", 
            "set_maximum_number_of_particles",
            "maximum_number_of_particles", 
            "the code will evolve this number of particles, please be sure to account for mergers", 
            default_value = 100
        )
        handler.add_boolean_parameter(
            "get_evolve_to_exact_time", 
            "set_evolve_to_exact_time",
            "evolve_to_exact_time", 
            "the code will evolve the model to the exact time given in evolve_model", 
            True
        )
        
    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        
        handler.add_method(
            "set_lightspeed",
            ( nbody_system.length / nbody_system.time, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_lightspeed",
            ( ),
            (nbody_system.length / nbody_system.time, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_time_step",
            (nbody_system.time, ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_tolerance",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_tolerance",
            ( ),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_radiated_gravitational_energy",
            (),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_total_energy",
            (),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, handler.ERROR_CODE,)
        )

    def define_properties(self, handler):
        
        GravitationalDynamics.define_properties(self, handler)
        handler.add_property("get_radiated_gravitational_energy")
    
    def update_particle_set(self):
        """
        update the particle set after changes in the code
        
        this implementation needs to move to the
        amuse.datamodel.incode_storage module, as
        it uses a lot of internal methods and info!
        
        """  
        
        number_of_added_particles = self.get_number_of_particles_added()
        if number_of_added_particles == 0:
            return
        
        indices_in_update_list = list(range(number_of_added_particles))
        indices_to_add = self.get_id_of_added_particle(indices_in_update_list)
        
        incode_storage = self.particles._private.attribute_storage
        
        if len(indices_to_add) > 0:
            incode_storage._add_indices(indices_to_add)


    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)
        
