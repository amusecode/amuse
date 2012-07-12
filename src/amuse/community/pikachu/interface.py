from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

class PikachuInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, StoppingConditionInterface):
    """
    Pikachu - a.k.a. P^3 Tree
    Hybrid N-body module, combining a tree (Barnes & Hut) to approximate long-range 
    forces, with direct summation of the forces from neighbour particles.
    """
    include_headers = ['worker_code.h', 'stopcond.h', 'interface.h']
    
    MODE_NORMAL = 'normal'
    MODE_LARGE_N = 'large_n'
    
    def __init__(self, mode=MODE_NORMAL, **options):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)
    
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'pikachu_worker'
        elif mode == self.MODE_LARGE_N:
            return 'pikachu_worker_large_n'
        else:
            print "Warning: unknown mode: '{0}' - using default ('{1}').".format(mode, self.MODE_NORMAL)
            return 'pikachu_worker'
    
    @option(type="string", sections=('data',))
    def amuse_root_directory(self):
        """
        The root directory of AMUSE, used as default root for all data directories
        """
        return get_amuse_root_dir()
    
    @option(type="string", sections=('data',))
    def output_data_root_directory(self):
        """
        The root directory of the output data,
        read - write directory
        """
        return os.path.join(self.amuse_root_directory, 'data')
    
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, 'pikachu', 'output')
    
    @legacy_function
    def get_gravity_at_point():
        """
        Determine the gravitational force at a given point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for par in ['eps', 'x', 'y', 'z']:
            function.addParameter(par, dtype='float64', direction=function.IN)
        for par in ['ax', 'ay', 'az']:
            function.addParameter(par, dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_potential_at_point():
        """
        Determine the potential at a given point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for par in ['eps', 'x', 'y', 'z']:
            function.addParameter(par, dtype='float64', direction=function.IN)
        function.addParameter('phi', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    

class Pikachu(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = PikachuInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)
        self.stopping_conditions.define_parameters(object)
    
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        self.stopping_conditions.define_methods(object)
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object)
