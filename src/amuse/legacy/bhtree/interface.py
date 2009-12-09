# -*- coding: utf-8 -*-
from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.interface.gd import NBodyGravitationalDynamicsBinding
from amuse.legacy.support.lit import LiteratureRefs

class BHTreeInterface(LegacyInterface, LiteratureRefs, GravitationalDynamics):
    """
        .. [#] Barnes, J., Hut, P., A Hierarchical O(N log N) force-calculation algorithm, *Nature*, **4**, 324 (1986)   
    """
    include_headers = ['BHTree_code.h', 'parameters.h', 'local.h']
    
    timestep = legacy_global(name='timestep',id=21,dtype='d')
    eps2_for_gravity = legacy_global(name='eps2_for_gravity',id=22,dtype='d')
    theta_for_tree = legacy_global(name='theta_for_tree',id=23,dtype='d')
    
    use_self_gravity = legacy_global(name='use_self_gravity',id=24,dtype='i')
    ncrit_for_tree = legacy_global(name='ncrit_for_tree',id=25,dtype='i')
    
    dt_dia = legacy_global(name='dt_dia',id=246,dtype='d')

    def __init__(self, convert_nbody = None, **kwargs):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code", **kwargs)
        """
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
            
        self.convert_nbody = convert_nbody
        """
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
    def initialize_particles():
        function = LegacyFunctionSpecification() 
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
       
    @legacy_function  
    def reinitialize_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
        
    #def evolve_model(self, time_end):
    #    result = self.evolve(self.convert_nbody.to_nbody(time_end).value_in(nbody_system.time), 1)
    #    return result
            
    def add_particles(self, particles):
        keyword_arguments = {}
        for attribute_definition in self.attribute_definitions:
            values = particles.get_values_of_attribute(attribute_definition.name)
            attribute_definition.set_keyword_arguments(self, values, keyword_arguments)
        keyword_arguments['id'] = list(particles.ids)
        
        self.add_particle(**keyword_arguments)
            
    #def update_attributes(self, attributes):
    #    for id, x in attributes:
    #        if x.name == 'mass':
    #            self.set_mass(id, self.convert_nbody.to_nbody(x.value()).value_in(nbody_system.mass))
   
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
        attributes.AttributeDefinition(
            name = "mass",
            setup_parameters = ["mass"],
            setter = ("set_mass", ["mass"]),
            description = "mass of a star",
            unit = nbody_system.mass,
            default = 1 | nbody_system.mass          
        ),
        attributes.AttributeDefinition(
            name = "radius",
            setup_parameters = ["radius"],
            setter = ("set_radius", ["radius"]),
            description = "radius of a star",
            unit = nbody_system.length,
            default = 1 | nbody_system.length          
        ),
        attributes.AttributeDefinition(
            names = ["x","y","z"],
            setup_parameters = ["x","y","z"],
            setter = ("set_position", ["x","y","z"]),
            description = "coordinate of a star",
            unit = nbody_system.length,
            default = 0.0 | nbody_system.length          
        ),
        attributes.AttributeDefinition(
            names = ["vx","vy","vz"],
            setup_parameters = ["vx","vy","vz"],
            setter = ("set_velocity", ["vx","vy","vz"]),
            description = "coordinate of a star",
            unit = nbody_system.speed,
            default = 0.0 | nbody_system.speed          
        ),
    ]    

    def __init__(self, convert_nbody = None):
        NBodyGravitationalDynamicsBinding.__init__(self, convert_nbody)
        
    def current_model_time(self):
        return self.convert_nbody.to_si( self.get_time()['time'] | nbody_system.time)
    
    def new_particle(self, **keyword_arguments):
        x = keyword_arguments['x']
        keyword_arguments['id'] = numpy.arange(len(x))
       
        self.new_particle(**keyword_arguments)
        
        return keyword_arguments['id'],numpy.zeros(len(x))
       
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
