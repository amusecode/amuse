
from amuse.support.data import parameters


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
