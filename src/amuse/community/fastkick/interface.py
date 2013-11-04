from amuse.community import *
#~from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldCode, GravityFieldInterface
from amuse.community.interface.common import CommonCodeInterface, CommonCode

class FastKickInterface(CodeInterface, CommonCodeInterface, GravityFieldInterface):
    """
    """
    include_headers = ['worker_code.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="fastkick_worker", **options)

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def commit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def recommit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

class FastKickDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_doc+"\n\n"+instance.parameters.__doc__


class FastKick(CommonCode, GravityFieldCode):

    __doc__ = FastKickDoc()

    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        legacy_interface = FastKickInterface(**options)
        self.legacy_doc = legacy_interface.__doc__
        CommonCode.__init__(self, legacy_interface, **options)

    def define_methods(self, object):
        CommonCode.define_methods(self, object)
        object.add_method("new_particle", [nbody_system.mass] + [nbody_system.length]*3,
            (object.INDEX, object.ERROR_CODE))
        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )
    
    def define_state(self, object): 
        CommonCode.define_state(self, object)   
        object.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        object.add_transition('INITIALIZED', 'EDIT', 'commit_parameters')
        object.add_transition('RUN', 'CHANGE_PARAMETERS_RUN', 'before_set_parameter', False)
        object.add_transition('EDIT', 'CHANGE_PARAMETERS_EDIT', 'before_set_parameter', False)
        object.add_transition('UPDATE', 'CHANGE_PARAMETERS_UPDATE', 'before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN', 'RUN', 'recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT', 'EDIT', 'recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_UPDATE', 'UPDATE', 'recommit_parameters')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE', 'before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        object.add_method('UPDATE','before_get_parameter')
        
        object.add_method('EDIT', 'new_particle')
        object.add_method('EDIT', 'delete_particle')
        object.add_method('UPDATE', 'new_particle')
        object.add_method('UPDATE', 'delete_particle')
        object.add_transition('EDIT', 'RUN', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_particle', False)
        object.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        GravityFieldCode.define_state(self, object)
        
    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
            
    def commit_parameters(self):
        self.parameters.send_not_set_parameters_to_code()
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
        
    def cleanup_code(self):
        self.overridden().cleanup_code()
        handler = self.get_handler('PARTICLES')
        handler._cleanup_instances()

    def reset(self):
        parameters = self.parameters.copy()
        self.cleanup_code()
        self.initialize_code()
        self.parameters.reset_from_memento(parameters)
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
    
    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')
