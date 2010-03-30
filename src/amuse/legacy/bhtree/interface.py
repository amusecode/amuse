from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs


class BHTreeInterface(LegacyInterface, LiteratureRefs, GravitationalDynamicsInterface):
    """
        .. [#] Barnes, J., Hut, P., A Hierarchical O(N log N) force-calculation algorithm, *Nature*, **4**, 324 (1986)   
    """
    include_headers = [ 'bhtree_code.h', 'parameters.h', 'worker_code.h', 'local.h']
    
    #timestep = legacy_global(name='timestep',id=21,dtype='d')
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

    def setup_module(self):
        self.initialize_code()
        self.commit_parameters()
        self.commit_particles()
       
    def cleanup_module(self):
        self.cleanup_code()
    
    @legacy_function  
    def reinitialize_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_time_step():
        """
        Update timestep.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function    
        
       
class BHTree(GravitationalDynamics):
    
    
    def __init__(self, convert_nbody = None):
        
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
       
        
        legacy_interface = BHTreeInterface()
        
        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
        )     
            
    def define_parameters(self, object):
        object.add_attribute_parameter(
            "eps2_for_gravity",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.3 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            nbody_system.time, 
            0.7 | nbody_system.time
        )
        object.add_attribute_parameter(
            "theta_for_tree",
            "openings_angle", 
            "openings angle for building the tree between 0 and 1", 
            units.none,
            0.5 | units.none
        )
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        
        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration, nbody_system.acceleration, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )
        
        object.add_method(
            'set_time_step',
            (nbody_system.time),
            (object.ERROR_CODE)
            )
