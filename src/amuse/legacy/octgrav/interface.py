from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs

class OctgravInterface(LegacyInterface, LiteratureRefs, GravitationalDynamicsInterface):
    """
        .. [#] Gaburov, Nitadori, Harfst, Portegies Zwart & Makino,"A gravitational tree code on graphics processing units:
               Implementation in CUDA", in preparetion; and main MUSE paper, arXiv/0807.1996
    """

    include_headers = ['octgrav_code.h', 'parameters.h', 'worker_code.h', 'local.h']

    def __init__(self, convert_nbody = None, **options):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code", **options)
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
    def get_theta_for_tree():
        """
        
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('openings_angle', dtype='float64', direction=function.OUT,
            description = "openings_angle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            xx
        -1 - ERROR
            xx
        """
        return function    

    @legacy_function
    def set_theta_for_tree():
        """
        
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('openings_angle', dtype='float64', direction=function.IN,
            description = "openings_angle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            xx
        -1 - ERROR
            xx
        """
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
        

class Octgrav(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **options):

        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()

        legacy_interface = OctgravInterface(**options)

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.01 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            nbody_system.time, 
            0.01 | nbody_system.time
        )
        object.add_method_parameter(
            "get_theta_for_tree",
            "set_theta_for_tree",
            "openings_angle",
            "openings angle for building the tree between 0 and 1", 
            units.none,
            0.8 | units.none
        )
