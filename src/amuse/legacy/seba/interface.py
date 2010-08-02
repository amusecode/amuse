from amuse.legacy import *

class SebaInterface(LegacyInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        LegacyInterface.__init__(self, **options)
    
    @legacy_function
    def evolve_star():
        function = LegacyFunctionSpecification()  
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.addParameter('endtime', dtype='float64', direction=function.IN)
        function.addParameter('metal', dtype='float64', direction=function.IN)
        function.addParameter('resulttime', dtype='float64', direction=function.OUT)
        function.addParameter('end_mass', dtype='float64', direction=function.OUT)
        function.addParameter('end_radius', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
    
class Seba(CodeInterface):

    def __init__(self, **options):
        CodeInterface.__init__(self,  SebaInterface(**options), **options)
    
    def define_methods(self, object):
        object.add_method(
            "evolve_star",
            (units.MSun, units.Myr, units.none),
            (units.Myr, units.MSun, units.RSun, object.ERROR_CODE),
            public_name = 'evolve'
        )
    
    
    def define_particle_sets(self, object):
        object.define_inmemory_set('particles')
        
    def _evolve_particles(self, particles, end_time):
        attributes = (
            "age",
            "mass",
            "radius",
        )
        
        result = self.evolve(
            particles.mass,
            end_time.as_vector_with_length(len(particles)),
            particles.metallicity,
        )
        
        particles._set_values(particles._get_keys(), attributes, result)

    def evolve_model(self, end_time = None):
        if end_time is None:
            raise exceptions.LegacyException("Cannot determine end_time from the code yet, so end_time must be provided!")
            
        self._evolve_particles(self.particles, end_time)
