from amuse.community.interface.chem import ChemicalModelingInterface
from amuse.community.interface.chem import ChemicalModeling
from amuse.community.interface.common import CommonCode
from amuse.community import *

class UclchemInterface(CodeInterface):
    def __init__(self, mode = 'cpu', **options):
        CodeInterface.__init__(
            self,
            name_of_the_worker='uclchem_worker',
            **options
        )
    
    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['number_density','temperature','ionrate']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['number_density','temperature','ionrate']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('aid', dtype='i', direction=function.IN)
        function.addParameter('abundance', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('aid', dtype='i', direction=function.IN)
        function.addParameter('abundance', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_firstlast_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('first', dtype='i', direction=function.OUT)
        function.addParameter('last', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT)
        for x in ['dens','temperature','ionrate']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function   
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function

class Uclchem(ChemicalModeling):
    def __init__(self, convert_nbody=None, **options):
        legacy_interface = UclchemInterface(**options)

        ChemicalModeling.__init__(self,legacy_interface)
    
    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method(
            "set_state",
            (
                handler.INDEX,
                units.cm**-3,
                units.K,
                units.s**-1,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state",
            (
                handler.INDEX,
            ),
            (
                units.cm**-3,
                units.K,
                units.s**-1,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_abundance",
            (
                handler.INDEX,
                handler.INDEX,
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_abundance",
            (
                handler.INDEX,
                handler.INDEX,
                handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_firstlast_abundance",
            (
            ),
            (
                handler.NO_UNIT,
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "new_particle",
            (
                units.cm**-3,
                units.K,
                units.s**-1,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "delete_particle",
            (
                handler.INDEX,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        
    
    def define_particle_sets(self, handler):
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_state')
        handler.add_getter('particles', 'get_state')
        handler.add_gridded_getter('particles', 'get_abundance','get_firstlast_abundance', names = ('abundances',))
        handler.add_gridded_setter('particles', 'set_abundance','get_firstlast_abundance', names = ('abundances',))
    
    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','PARAMETER_CHANGE_A','invoke_state_change2')
        handler.add_transition('EDIT','PARAMETER_CHANGE_B','invoke_state_change2')
        handler.add_transition('PARAMETER_CHANGE_A','RUN','recommit_parameters')
        handler.add_transition('PARAMETER_CHANGE_B','EDIT','recommit_parameters')
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_method('RUN', 'evolve_model')
        handler.add_method('RUN', 'get_state')
        handler.add_method('RUN', 'get_abundance')


