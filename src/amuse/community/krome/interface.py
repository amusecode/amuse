import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option

#(Grevesse & Sauval, 1998, Space Sci. Rev. 85, 161)
solar_abundances= dict(H=1.,
                       HE=.085,
                       C=3.31e-4,
                       N=8.3e-5,
                       O=6.76e-4,
                       Ne=1.2e-4,
                       SI=3.55e-5,
                       Fe=3.2e-5)
 
class KromeInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
    
    KROME - a package to embed chemistry in astrophysical simulations

    .. [#] Grassi, T.; Bovino, S.; Schleicher, D. R. G.; Prieto, J.; Seifried, D.; Simoncini, E.; Gianturco, F. A., MNRAS, 439, 3, p.2386-2419 [2014MNRAS.439.2386G]
    """

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'krome_worker'

    @legacy_function
    def commit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def commit_parameters():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    def recommit_parameters():
        return self.commit_parameters()

    @legacy_function
    def recommit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT)
        for x in ['number_density','temperature','ionrate']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

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
    def set_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('aid', dtype='i', direction=function.IN)
        function.addParameter('abundance', dtype='d', direction=function.IN)
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
    def get_firstlast_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('first', dtype='i', direction=function.OUT)
        function.addParameter('last', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_index_of_species():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('name', dtype='s', direction=function.IN)
        function.addParameter('index', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_name_of_species():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('name', dtype='s', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function   
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function

class Krome(CommonCode):

    def __init__(self,unit_converter = None, **options):
        
        if unit_converter is not None:
            raise Exception("krome uses predefined units")
           
        InCodeComponentImplementation.__init__(self, KromeInterface(**options))
        first,last=self.get_firstlast_abundance()
        self.species=dict()
        for i in range(first,last+1):
          self.species[self.get_name_of_species(i)]=i-1
#        self.set_data_directory(self.data_directory())
#        self.set_output_directory(self.output_directory())
                  
    def define_properties(self, handler):
        handler.add_property('get_time', public_name = "model_time")

    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method(
            'evolve_model', 
            (
                units.s,
            ), 
            (
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
            "set_state",
            (
                handler.NO_UNIT,
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
            "delete_particle",
            (
                handler.INDEX,
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
            "get_time",
            (
            ),
            (
                units.s,
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


