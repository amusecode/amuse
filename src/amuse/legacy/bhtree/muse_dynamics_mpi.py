from amuse.legacy import *
from amuse.legacy.interface.gd import NBodyGravitationalDynamicsBinding

class BHTreeInterface(LegacyInterface):
    include_headers = ['muse_dynamics.h', 'parameters.h', 'local.h']
    
    timestep = legacy_global(name='timestep',id=21,dtype='d')
    eps2_for_gravity = legacy_global(name='eps2_for_gravity',id=22,dtype='d')
    theta_for_tree = legacy_global(name='theta_for_tree',id=23,dtype='d')
    
    use_self_gravity = legacy_global(name='use_self_gravity',id=24,dtype='i')
    ncrit_for_tree = legacy_global(name='ncrit_for_tree',id=25,dtype='i')
    
    dt_dia = legacy_global(name='dt_dia',id=246,dtype='d')
    
    
      
    def __init__(self, convert_nbody = None, **kwargs):
        LegacyInterface.__init__(self, **kwargs)
        

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
    def initialize_particles():
        function = LegacyFunctionSpecification() 
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
        
    @legacy_function    
    def add_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('id_out', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = None
        return function
        
    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function  
    def reinitialize_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def get_number():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function;
     
    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        return function;
    
    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function;
    
    @legacy_function      
    def get_kinetic_energy():
        function = LegacyFunctionSpecification() 
        function.result_type = 'd'
        return function

    @legacy_function      
    def get_potential_energy():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function
         
    
  

class BHTreeBinding(NBodyGravitationalDynamicsBinding):
    parameter_definitions = [
        parameters.ModuleAttributeParameterDefinition(
            "eps2_for_gravity",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.3 | nbody_system.length * nbody_system.length
        )
    ]
    
    attribute_definitions = [
        attributes.ScalarAttributeDefinition(
            "set_mass",
            None,
            "mass",
            "mass",
            "mass of a star",
             nbody_system.mass,
             1 | nbody_system.mass
        ),
        attributes.ScalarAttributeDefinition(
            None,
            None,
            "radius",
            "radius",
            "radius of a star, used for collision detection",
             nbody_system.length,
             1 | nbody_system.length
        ),
        attributes.VectorAttributeDefinition(
            None,
            None,
            ["x","y","z"],
            "position",
            "cartesian coordinates of a star",
             nbody_system.length,
             [0.0, 0.0, 0.0] | nbody_system.length
        ),
        attributes.VectorAttributeDefinition(
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
        NBodyGravitationalDynamicsBinding.__init__(self, convert_nbody)
        
    def current_model_time(self):
        return self.convert_nbody.to_si( self.get_time() | nbody_system.time)
    
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
            
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).value_in(nbody_system.time), 1)
        return result
        
        
    def add_particles(self, particles):
        keyword_arguments = {}
        for attribute_definition in self.attribute_definitions:
            values = particles.get_values_of_attribute(attribute_definition.name)
            attribute_definition.set_keyword_arguments(self, values, keyword_arguments)
        keyword_arguments['id'] = list(particles.ids)
        
        indices = self.add_particle(**keyword_arguments)
        
        
        module_id = id(self)
        for particle, new_index in zip(particles, indices):
            particle._module_ids_to_index[module_id] = (self, new_index)
        
    def get_energies(self):
        energy_unit = nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2
        kinetic_energy = self.get_kinetic_energy() | energy_unit
        potential_energy = self.get_potential_energy() | energy_unit
        return (self.convert_nbody.to_si(kinetic_energy), self.convert_nbody.to_si(potential_energy))
    
        
        
class BHTree(BHTreeInterface, BHTreeBinding):
    """ 
    """	
    
    def __init__(self, convert_nbody = None):
        BHTreeInterface.__init__(self)
        BHTreeBinding.__init__(self, convert_nbody)
