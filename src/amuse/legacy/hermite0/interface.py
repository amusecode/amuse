from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs

class HermiteInterface(LegacyInterface, LiteratureRefs, GravitationalDynamics):
    """ 
        .. [#] Hut, P., Makino, J. & McMillan, S., *Astrophysical Journal Letters* , **443**, L93-L96 (1995)
    """		
    include_headers = ['hermite_code.h', 'parameters.h', 'local.h']

    t = legacy_global(name='t',id=20,dtype='d')
    dt_param = legacy_global(name='dt_param',id=21,dtype='d')
    dt_dia = legacy_global(name='dt_dia',id=22,dtype='d')
    eps2 = legacy_global(name='eps2',id=23,dtype='d')
    flag_collision = legacy_global(name='flag_collision',id=24,dtype='i')

    def __init__(self, convert_nbody = None):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code")       
        LiteratureRefs.__init__(self)

    @legacy_function
    def setup_module():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        return function

    @legacy_function
    def cleanup_module():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        return function

    @legacy_function  
    def reinitialize_particles():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_particle', dtype='int32', direction=function.IN,
            description = "Index of particle to delete")
        function.result_type = 'int32'
        function.resutl_doc = """
        0 - OK
           The particle was deleted
        """
        return function
        
class HermiteBinding(object):
    parameter_definitions = [
        parameters.ModuleAttributeParameterDefinition(
            "eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.0 | nbody_system.length * nbody_system.length
        )
    ]
    
    attribute_definitions = [
        attributes.ScalarAttributeDefinition_Next(
            "set_mass",
            None,
            "mass",
            "mass",
            "mass of a star",
             nbody_system.mass,
             1 | nbody_system.mass
        ),
        attributes.ScalarAttributeDefinition_Next(
            None,
            None,
            "radius",
            "radius",
            "radius of a star, used for collision detection",
             nbody_system.length,
             1 | nbody_system.length
        ),
        attributes.VectorAttributeDefinition_Next(
            None,
            None,
            ["x","y","z"],
            "position",
            "cartesian coordinates of a star",
             nbody_system.length,
             [0.0, 0.0, 0.0] | nbody_system.length
        ),
        attributes.VectorAttributeDefinition_Next(
            None,
            None,
            ["vx","vy","vz"],
            "velocity",
            "velocity of a star",
            nbody_system.speed,
            [0.0, 0.0, 0.0] | nbody_system.speed
        ),
    ]

    def __init__(self, convert_nbody = None):
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
        self.convert_nbody = convert_nbody
        
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        
    
    def setup_particles(self, particles):
        keyword_arguments = {}
        for attribute_definition in self.attribute_definitions:
            values = particles.get_values_of_attribute(attribute_definition.name)
            attribute_definition.set_keyword_arguments(self, values, keyword_arguments)
        
        indices, errors = self.new_particle(**keyword_arguments)
        
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not setup a particle")
        
        module_id = id(self)
        for particle, new_index in zip(particles, indices):
            particle._module_ids_to_index[module_id] = (self, new_index)
            

    def update_particles(self, particles):
        state = self.get_state(list(particles.ids_for_module_with_id(id(self))))
        time = self.convert_nbody.to_si( self.t | nbody_system.time)
        for attribute_definition in self.attribute_definitions:
            values = attribute_definition.get_keyword_results(self, state)
            particles.set_values_of_attribute(attribute_definition.name, time, values)

    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).value_in(nbody_system.time))
        return result
        
class Hermite(HermiteInterface, HermiteBinding):
    
    def __init__(self, convert_nbody = None):
        HermiteInterface.__init__(self)
        HermiteBinding.__init__(self, convert_nbody)
        

