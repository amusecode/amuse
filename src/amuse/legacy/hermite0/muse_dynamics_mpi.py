from amuse.legacy import *

class Hermite(LegacyInterface):
    include_headers = ['muse_dynamics.h', 'parameters.h', 'local.h']
    
    t = legacy_global(name='t',id=20,dtype='d')
    dt_param = legacy_global(name='dt_param',id=21,dtype='d')
    dt_dia = legacy_global(name='dt_dia',id=22,dtype='d')
    eps2 = legacy_global(name='eps2',id=23,dtype='d')
    flag_collision = legacy_global(name='flag_collision',id=24,dtype='i')
            
            
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
    
    def __init__(self, convert_nbody = nbody_system.noconvert_nbody_to_si()):
        LegacyInterface.__init__(self)
        self.convert_nbody = convert_nbody
        self.parameters = parameters.Parameters(self.parameter_definitions, self)

    @legacy_function   
    def setup_module():
        function = RemoteFunction() 
        function.result_type = 'i'
        return function
        
    @legacy_function      
    def cleanup_module():
        function = RemoteFunction() 
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def initialize_particles():
        function = RemoteFunction() 
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
        
    @legacy_function    
    def add_particle():
        function = RemoteFunction() 
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_state():
        function = RemoteFunction()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('id_out', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = None
        return function
        
    @legacy_function    
    def evolve():
        function = RemoteFunction()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function  
    def reinitialize_particles():
        function = RemoteFunction() 
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def get_number():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function;
     
    @legacy_function
    def set_mass():
        function = RemoteFunction() 
        function.result_type = 'i'
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        return function;
        
    @legacy_function      
    def get_kinetic_energy():
        function = RemoteFunction() 
        function.result_type = 'd'
        return function

    @legacy_function      
    def get_potential_energy():
        function = RemoteFunction()  
        function.result_type = 'd'
        return function
         
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).value_in(nbody_system.time), 0.0)
        return result
        
    def set_eps2(self, eps2):
        self.eps2 = self.convert_nbody.to_nbody(eps2).number
    
    
    def update_particles(self, particles):
        state = self.get_state(list(particles.ids))
        time = self.convert_nbody.to_si( self.t | nbody_system.time)
        for attribute_definition in self.attribute_definitions:
            values = attribute_definition.get_keyword_results(self, state)
            particles.set_values_of_attribute(attribute_definition.name, time, values)
            
    def add_particles(self, particles):
        keyword_arguments = {}
        for attribute_definition in self.attribute_definitions:
            values = particles.get_values_of_attribute(attribute_definition.name)
            attribute_definition.set_keyword_arguments(self, values, keyword_arguments)
        keyword_arguments['id'] = list(particles.ids)
        self.add_particle(**keyword_arguments)
            
    def update_attributes(self, attributes):
        
        ids = []
        mass_of_particles = []
        for id, x in attributes:
            if x.name == 'mass':
                id.append(x.id)
                mass = self.convert_nbody.to_nbody(x.mass.value())
                mass_of_particles.append(mass.number)
        self.set_mass(ids, mass_of_particles)
        
    
    def get_energies(self):
        energy_unit = nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2
        kinetic_energy = self.get_kinetic_energy() | energy_unit
        potential_energy = self.get_potential_energy() | energy_unit
        return (self.convert_nbody.to_si(kinetic_energy), self.convert_nbody.to_si(potential_energy))
        
    
  
