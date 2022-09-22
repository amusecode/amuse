from amuse.community import *
#~from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldCode, GravityFieldInterface
from amuse.community.interface.common import CommonCodeInterface, CommonCode


class FastKickInterface(CodeInterface, CommonCodeInterface, GravityFieldInterface):
    """
    """
    include_headers = ['worker_code.h']

    MODE_CPU = 'cpu'
    MODE_GPU = 'gpu'

    def __init__(self, mode=MODE_CPU, **options):
        CodeInterface.__init__(self, name_of_the_worker=self.get_name_of_the_worker(mode), **options)

    def get_name_of_the_worker(self, mode):
        if mode == self.MODE_CPU:
            return "fastkick_worker"
        if mode == self.MODE_GPU:
            return "fastkick_worker_gpu"
        else:
            return "fastkick_worker"

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
    
    @legacy_function
    def get_potential_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy', dtype='float64', direction=function.OUT, unit=nbody_system.energy)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_mass():
        """
        Retrieve the mass of a particle. Mass is a scalar property of a particle,
        this function has one OUT argument.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_mass():
        """
        Update the mass of a particle. Mass is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        """
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

    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method("new_particle", [nbody_system.mass] + [nbody_system.length]*3,
            (handler.INDEX, handler.ERROR_CODE))
        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (handler.ERROR_CODE,)
        )
    
        handler.add_method(
            "set_mass",
            (
                handler.NO_UNIT,
                nbody_system.mass,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_mass",
            (
                handler.NO_UNIT,
            ),
            (
                nbody_system.mass,
                handler.ERROR_CODE
            )
        )
    def define_state(self, handler): 
        CommonCode.define_state(self, handler)   
        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        handler.add_transition('INITIALIZED', 'EDIT', 'commit_parameters')
        handler.add_transition('RUN', 'CHANGE_PARAMETERS_RUN', 'before_set_parameter', False)
        handler.add_transition('EDIT', 'CHANGE_PARAMETERS_EDIT', 'before_set_parameter', False)
        handler.add_transition('UPDATE', 'CHANGE_PARAMETERS_UPDATE', 'before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN', 'RUN', 'recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT', 'EDIT', 'recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_UPDATE', 'UPDATE', 'recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE', 'before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        handler.add_method('UPDATE','before_get_parameter')
        
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_method('UPDATE', 'new_particle')
        handler.add_method('UPDATE', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        
        GravityFieldCode.define_state(self, handler)
        handler.add_method('RUN', 'get_potential_energy')
        
    def define_converter(self, handler):
        if not self.unit_converter is None:
            handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())
            
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
    
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
    
    def define_particle_sets(self, handler):
        handler.define_set('particles', 'index_of_the_particle')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
    
    def define_properties(self, handler):
        handler.add_property("get_potential_energy")


Fastkick = FastKick
