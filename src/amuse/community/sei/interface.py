from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.units import nbody_system
from amuse.community.interface import common

class SeiInterface(CodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="sei_worker", **keyword_arguments)
    
    @legacy_function
    def initialization():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function

    @legacy_function
    def select_integrator():
        function = LegacyFunctionSpecification()  
        function.addParameter('i', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.IN)
        function.addParameter('vy', dtype='float64', direction=function.IN)
        function.addParameter('vz', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN, description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('t', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN, description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.IN)
        function.addParameter('vy', dtype='float64', direction=function.IN)
        function.addParameter('vz', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN, description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.addParameter('vx', dtype='float64', direction=function.OUT)
        function.addParameter('vy', dtype='float64', direction=function.OUT)
        function.addParameter('vz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve():
        function = LegacyFunctionSpecification()
        function.addParameter('t_end', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN, description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.result_type = 'int32'
        return function
        
class Sei(common.CommonCode):

    def __init__(self, unit_converter=None, **options):
        legacy_interface = SeiInterface(**options)
        self.unit_converter = unit_converter
        InCodeComponentImplementation.__init__(self, legacy_interface, **options)
        
    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')
        object.add_setter('particles', 'set_state')
        object.add_getter('particles', 'get_state')

    def define_methods(self, object):
        #GravitationalDynamics.define_methods(self, object)
        
        common.CommonCode.define_methods(self, object)
        object.add_method(
            'evolve',
            (nbody_system.time,),
            public_name = 'evolve_model'
        )

        object.add_method(
            "new_particle",
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )

        object.add_method(
            "delete_particle",
            (
                object.INDEX
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_state",
            (
                object.INDEX,
            ),
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_state",
            (
                object.INDEX,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                object.ERROR_CODE
            )
        )


    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
