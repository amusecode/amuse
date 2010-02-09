
from amuse.support.data import parameters
from amuse.support.data import core
from amuse.support.data import values
from amuse.support.units import nbody_system


import inspect

class OldObjectsBindingMixin(object):
                
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def update_particles(self, particles):
        self.particles.copy_values_of_state_attributes_to(particles)
    
    def set_attribute(self, attribute_name, particles):
        particles.copy_values_of_attribute_to(attribute_name, self.particles)
        
    def update_attribute(self, attribute_name, particles):
        self.particles.copy_values_of_attribute_to(attribute_name, particles) 
        
        
        
class CodeInterface(OldObjectsBindingMixin):
    """Base class of next level interfaces to codes
    """
    
    parameter_definitions = []
    
    def __init__(self):
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        

class CodeInterfaceWithConvertedUnits(OldObjectsBindingMixin):
    class UnitsConvertionMethod(object):
        def __init__(self, real_method, converter):
            self.real_method = real_method
            self.converter = converter
            
        def __call__(self, *list_arguments, **keyword_arguments):
            converted_list_arguments = [self.from_source_to_target(x) for x in list_arguments]
            converted_keyword_arguments = {}
            for key, value in keyword_arguments:
                converted_keyword_arguments[key] = self.from_source_to_target(value)
            
            result = self.real_method(*converted_list_arguments, **converted_keyword_arguments)
            print result
            return self.from_target_to_source(result)
            
        def from_source_to_target(self, x):
            if isinstance(x, values.Quantity):
                return self.converter.from_source_to_target(x)
            else:
                return x
                
        def from_target_to_source(self, x):
            if isinstance(x, values.Quantity):
                return self.converter.from_target_to_source(x)
            elif hasattr(x, '__iter__'):
                return self.convert_and_iterate(x)
            else:
                return x
        
        def convert_and_iterate(self, iterable):
            for x in iterable:
                yield self.converter.from_target_to_source(x)
            
    def __init__(self, interface, converter):
        self.converter = converter
        self.interface = interface
        
    def __getattr__(self, name):
        attribute = getattr(self.interface, name)
        if inspect.ismethod(attribute):
            return self.UnitsConvertionMethod(attribute, self.converter)
        elif isinstance(attribute, core.AbstractParticleSet):
            result = core.ParticlesWithUnitsConverted(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif isinstance(attribute, values.Quantity):
            return self.converter.from_target_to_source(attribute)
        elif isinstance(attribute, parameters.Parameters):
            result = parameters.ParametersWithUnitsConverted(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif hasattr(attribute, '__iter__'):
            return self.convert_and_iterate(attribute)
        else:
            return attribute
            
    def convert_and_iterate(self, iterable):
        for x in iterable:
            yield self.converter.from_target_to_source(x)
            
    def __dir__(self):
        return dir(self.interface)
        
class CodeInterfaceWithNBodyUnitsConverted(CodeInterfaceWithConvertedUnits):
    
    def __init__(self, interface, convert_nbody):
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
        
        CodeInterfaceWithConvertedUnits.__init__(
            self,
            interface,
            convert_nbody.as_converter_from_si_to_nbody()
        )
        
        
        

