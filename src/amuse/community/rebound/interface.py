from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class ReboundInterface(CodeInterface,
                       LiteratureReferencesMixIn,
                       GravitationalDynamicsInterface,
                       StoppingConditionInterface,
                       #SinglePointGravityFieldInterface
    ):
    """
    REBOUND - An open-source multi-purpose N-body code
    
    .. [#] Rein, H., Liu, S.F., *Astronomy and Astrophysics* , **Volume 537**, A128 (2012)
    
    """
    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="rebound_worker",
                                 **options)
        LiteratureReferencesMixIn.__init__(self)
        
    @legacy_function
    def _set_integrator():
        function = LegacyFunctionSpecification()  
        function.addParameter('integrator_name', dtype='i', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
        
    @legacy_function
    def _get_integrator():
        function = LegacyFunctionSpecification()  
        function.addParameter('integrator_name', dtype='i', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
    
    INTEGRATORS = {"ias15": 0, "whfast": 1, "sei": 2, "wh": 3, "leapfrog": 4, "hybrid": 5, "none": 6}
    def set_integrator(self, name):
        return self._set_integrator(self.INTEGRATORS[name])
    
    def get_integrator(self):
        value, error = self._get_integrator()
        for key, index in self.INTEGRATORS.iteritems():
            if value == index:
                return key, error
        return "none", -1
    
class Rebound(GravitationalDynamics, GravityFieldCode):


    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = ReboundInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        #GravityFieldCode.define_state(self, object)
        self.stopping_conditions.define_state(object)
        

    def define_parameters(self, object):
        self.stopping_conditions.define_parameters(object)

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        self.stopping_conditions.define_methods(object)
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        
        self.stopping_conditions.define_particle_set(object)
