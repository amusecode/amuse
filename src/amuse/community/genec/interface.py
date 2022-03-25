import numpy
from amuse.datamodel import ParticlesWithFilteredAttributes
from amuse.community import CodeInterface
from amuse.community import LegacyFunctionSpecification
from amuse.community import legacy_function
from amuse.community import LiteratureReferencesMixIn
from amuse.community import CodeWithDataDirectories
from amuse.community import InCodeComponentImplementation
from amuse.community.interface.se import StellarEvolution
from amuse.community.interface.se import StellarEvolutionInterface
from amuse.community.interface.se import InternalStellarStructure
from amuse.community.interface.se import InternalStellarStructureInterface
from amuse.community.interface.stopping_conditions import (
    StoppingConditionInterface, StoppingConditions
)
from amuse.units import units


class GenecInterface(
    CodeInterface,
    LiteratureReferencesMixIn,
    StellarEvolutionInterface,
    InternalStellarStructureInterface,
    CodeWithDataDirectories,
    # StoppingConditionInterface,
):
    """
    GENEC is the Geneva Stellar Evolution Code

    References:
        .. [#] The Geneva Stellar Evolution Group
    """

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="genec_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)

    @legacy_function
    def finalize_stellar_model():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def write_genec_model():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_genec_path():
        """
        set path to Genec data directories
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'genec_path', dtype='string', direction=function.IN,
            description="Path to Genec",
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function

    # Parameters

    @legacy_function
    def get_par_ipoly():
        'get parameter ipoly'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ipoly', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_ipoly():
        'set parameter ipoly'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ipoly', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_nwseq():
        'get parameter nwseq'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nwseq', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_nwseq():
        'set parameter nwseq'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nwseq', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_modanf():
        'get parameter modanf'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'modanf', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_modanf():
        'set parameter modanf'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'modanf', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_nzmod():
        'get parameter nzmod'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nzmod', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_nzmod():
        'set parameter nzmod'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nzmod', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_irot():
        'get parameter irot'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'irot', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_irot():
        'set parameter irot'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'irot', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_isol():
        'get parameter isol'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'isol', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_isol():
        'set parameter isol'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'isol', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_imagn():
        'get parameter imagn'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'imagn', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_imagn():
        'set parameter imagn'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'imagn', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_ialflu():
        'get parameter ialflu'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ialflu', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_ialflu():
        'set parameter ialflu'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ialflu', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_ianiso():
        'get parameter ianiso'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ianiso', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_ianiso():
        'set parameter ianiso'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ianiso', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_ipop3():
        'get parameter ipop3'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ipop3', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_ipop3():
        'set parameter ipop3'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ipop3', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_ibasnet():
        'get parameter ibasnet'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ibasnet', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_ibasnet():
        'set parameter ibasnet'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ibasnet', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_phase():
        'get parameter phase'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'phase', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_phase():
        'set parameter phase'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'phase', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iopac():
        'get parameter iopac'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iopac', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iopac():
        'set parameter iopac'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iopac', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_ikappa():
        'get parameter ikappa'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ikappa', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_ikappa():
        'set parameter ikappa'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ikappa', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_idiff():
        'get parameter idiff'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idiff', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_idiff():
        'set parameter idiff'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idiff', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iadvec():
        'get parameter iadvec'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iadvec', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iadvec():
        'set parameter iadvec'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iadvec', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_istati():
        'get parameter istati'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'istati', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_istati():
        'set parameter istati'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'istati', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_icoeff():
        'get parameter icoeff'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'icoeff', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_icoeff():
        'set parameter icoeff'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'icoeff', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_igamma():
        'get parameter igamma'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'igamma', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_igamma():
        'set parameter igamma'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'igamma', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_idialo():
        'get parameter idialo'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idialo', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_idialo():
        'set parameter idialo'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idialo', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_idialu():
        'get parameter idialu'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idialu', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_idialu():
        'set parameter idialu'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idialu', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_imloss():
        'get parameter imloss'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'imloss', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_imloss():
        'set parameter imloss'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'imloss', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_ifitm():
        'get parameter ifitm'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ifitm', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_ifitm():
        'set parameter ifitm'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ifitm', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_nndr():
        'get parameter nndr'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nndr', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_nndr():
        'set parameter nndr'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nndr', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iledou():
        'get parameter iledou'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iledou', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iledou():
        'set parameter iledou'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iledou', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_idifcon():
        'get parameter idifcon'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idifcon', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_idifcon():
        'set parameter idifcon'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idifcon', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_my():
        'get parameter my'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'my', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_my():
        'set parameter my'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'my', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iover():
        'get parameter iover'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iover', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iover():
        'set parameter iover'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iover', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iunder():
        'get parameter iunder'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iunder', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iunder():
        'set parameter iunder'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iunder', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_nbchx():
        'get parameter nbchx'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nbchx', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_nbchx():
        'set parameter nbchx'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nbchx', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_nrband():
        'get parameter nrband'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nrband', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_nrband():
        'set parameter nrband'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'nrband', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_islow():
        'get parameter islow'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'islow', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_islow():
        'set parameter islow'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'islow', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_icncst():
        'get parameter icncst'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'icncst', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_icncst():
        'set parameter icncst'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'icncst', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_tauH_fit():
        'get parameter tauH_fit'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'tauH_fit', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_tauH_fit():
        'set parameter tauH_fit'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'tauH_fit', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iauto():
        'get parameter iauto'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iauto', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iauto():
        'set parameter iauto'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iauto', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iprn():
        'get parameter iprn'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iprn', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iprn():
        'set parameter iprn'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iprn', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_iout():
        'get parameter iout'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iout', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_iout():
        'set parameter iout'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'iout', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_itmin():
        'get parameter itmin'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'itmin', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_itmin():
        'set parameter itmin'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'itmin', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_idebug():
        'get parameter idebug'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idebug', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_idebug():
        'set parameter idebug'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'idebug', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_itests():
        'get parameter itests'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'itests', dtype='int32',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_itests():
        'set parameter itests'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'itests', dtype='int32',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_var_rates():
        'get parameter var_rates'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'var_rates', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_var_rates():
        'set parameter var_rates'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'var_rates', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_bintide():
        'get parameter bintide'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'bintide', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_bintide():
        'set parameter bintide'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'bintide', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_const_per():
        'get parameter const_per'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'const_per', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_const_per():
        'set parameter const_per'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'const_per', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_Add_Flux():
        'get parameter Add_Flux'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'Add_Flux', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_Add_Flux():
        'set parameter Add_Flux'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'Add_Flux', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_diff_only():
        'get parameter diff_only'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'diff_only', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_diff_only():
        'set parameter diff_only'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'diff_only', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_lowRSGMdot():
        'get parameter lowRSGMdot'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'lowRSGMdot', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_lowRSGMdot():
        'set parameter lowRSGMdot'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'lowRSGMdot', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_plot():
        'get parameter plot'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'plot', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_plot():
        'set parameter plot'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'plot', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_refresh():
        'get parameter refresh'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'refresh', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_refresh():
        'set parameter refresh'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'refresh', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_xyfiles():
        'get parameter xyfiles'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xyfiles', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_xyfiles():
        'set parameter xyfiles'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xyfiles', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_verbose():
        'get parameter verbose'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'verbose', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_verbose():
        'set parameter verbose'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'verbose', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_stop_deg():
        'get parameter stop_deg'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'stop_deg', dtype='bool',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_stop_deg():
        'set parameter stop_deg'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'stop_deg', dtype='bool',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_index_poly():
        'get parameter index_poly'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'index_poly', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_index_poly():
        'set parameter index_poly'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'index_poly', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function
    @legacy_function
    def get_par_binm2():
        'get parameter binm2'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'binm2', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_binm2():
        'set parameter binm2'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'binm2', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_periodini():
        'get parameter periodini'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'periodini', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_periodini():
        'set parameter periodini'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'periodini', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_zinit():
        'get parameter zinit'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'zinit', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_zinit():
        'set parameter zinit'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'zinit', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_zsol():
        'get parameter zsol'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'zsol', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_zsol():
        'set parameter zsol'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'zsol', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_z():
        'get parameter z'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'z', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_z():
        'set parameter z'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'z', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_fenerg():
        'get parameter fenerg'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fenerg', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_fenerg():
        'set parameter fenerg'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fenerg', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_richac():
        'get parameter richac'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'richac', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_richac():
        'set parameter richac'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'richac', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_frein():
        'get parameter frein'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'frein', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_frein():
        'set parameter frein'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'frein', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_K_Kawaler():
        'get parameter K_Kawaler'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'K_Kawaler', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_K_Kawaler():
        'set parameter K_Kawaler'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'K_Kawaler', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_Omega_saturation():
        'get parameter Omega_saturation'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'Omega_saturation', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_Omega_saturation():
        'set parameter Omega_saturation'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'Omega_saturation', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_rapcrilim():
        'get parameter rapcrilim'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'rapcrilim', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_rapcrilim():
        'set parameter rapcrilim'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'rapcrilim', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_vwant():
        'get parameter vwant'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'vwant', dtype='float64',
            direction=function.OUT,
            unit=units.kms
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_vwant():
        'set parameter vwant'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'vwant', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_xfom():
        'get parameter xfom'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xfom', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_xfom():
        'set parameter xfom'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xfom', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_omega():
        'get parameter omega'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'omega', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_omega():
        'set parameter omega'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'omega', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_xdial():
        'get parameter xdial'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xdial', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_xdial():
        'set parameter xdial'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xdial', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_B_initial():
        'get parameter B_initial'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'B_initial', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_B_initial():
        'set parameter B_initial'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'B_initial', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_add_diff():
        'get parameter add_diff'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'add_diff', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_add_diff():
        'set parameter add_diff'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'add_diff', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_fmlos():
        'get parameter fmlos'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fmlos', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_fmlos():
        'set parameter fmlos'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fmlos', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_fitm():
        'get parameter fitm'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fitm', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_fitm():
        'set parameter fitm'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fitm', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_fitmi():
        'get parameter fitmi'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fitmi', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_fitmi():
        'set parameter fitmi'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fitmi', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_fitmi_default():
        'get parameter fitmi_default'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fitmi_default', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_fitmi_default():
        'set parameter fitmi_default'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'fitmi_default', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_deltal():
        'get parameter deltal'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'deltal', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_deltal():
        'set parameter deltal'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'deltal', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_deltat():
        'get parameter deltat'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'deltat', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_deltat():
        'set parameter deltat'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'deltat', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_elph():
        'get parameter elph'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'elph', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_elph():
        'set parameter elph'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'elph', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dovhp():
        'get parameter dovhp'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dovhp', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dovhp():
        'set parameter dovhp'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dovhp', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dunder():
        'get parameter dunder'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dunder', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dunder():
        'set parameter dunder'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dunder', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_gkorm():
        'get parameter gkorm'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'gkorm', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_gkorm():
        'set parameter gkorm'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'gkorm', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_alph():
        'get parameter alph'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'alph', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_alph():
        'set parameter alph'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'alph', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_agdr():
        'get parameter agdr'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'agdr', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_agdr():
        'set parameter agdr'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'agdr', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_faktor():
        'get parameter faktor'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'faktor', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_faktor():
        'set parameter faktor'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'faktor', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dgrp():
        'get parameter dgrp'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgrp', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dgrp():
        'set parameter dgrp'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgrp', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dgrl():
        'get parameter dgrl'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgrl', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dgrl():
        'set parameter dgrl'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgrl', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dgry():
        'get parameter dgry'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgry', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dgry():
        'set parameter dgry'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgry', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dgrc():
        'get parameter dgrc'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgrc', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dgrc():
        'set parameter dgrc'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgrc', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dgro():
        'get parameter dgro'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgro', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dgro():
        'set parameter dgro'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgro', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_dgr20():
        'get parameter dgr20'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgr20', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_dgr20():
        'set parameter dgr20'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dgr20', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_xcn():
        'get parameter xcn'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xcn', dtype='float64',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_xcn():
        'set parameter xcn'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'xcn', dtype='float64',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    @legacy_function
    def get_par_starname():
        'get parameter starname'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'starname', dtype='string',
            direction=function.OUT,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            Got the value.
        -1 - ERROR
            Unable to get.
        '''
        return function

    @legacy_function
    def set_par_starname():
        'set parameter starname'
        function = LegacyFunctionSpecification()
        function.addParameter(
            'starname', dtype='string',
            direction=function.IN,
        )
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        '''
        return function

    # End parameters

    @legacy_function
    def commit_parameters():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_starname():
        """
        Set the star name (identical to AMUSE particle key?)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'index_of_the_star', dtype='int32',
            direction=function.IN,
            description="The star's key"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        """
        return function

    @legacy_function
    def new_particle():
        """
        Define a new star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.OUT,
            description=(
                "The new index for the star. This index can be used to refer "
                "to this star in other functions"
            )
        )
        function.addParameter(
            'mass', dtype='float64', direction=function.IN,
            description="The initial mass of the star")
        function.addParameter(
            'metallicity', dtype='float64', direction=function.IN,
            default=0.014,
            description="The initial metallicity of the star (default: 0.014)")
        function.addParameter(
            'starname', dtype='string', direction=function.IN,
            default='AmuseStar', description="The star's name")
        # function.addParameter(
        #     'age_tag', dtype='float64', direction=function.IN,
        #     description="Starting age of the star *to be specified exactly*")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
        -1 - ERROR
            New star could not be created.
        """
        return function

    @legacy_function
    def get_firstlast_zone():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('first', dtype='int32', direction=function.OUT)
        function.addParameter('last', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_luminosity_at_zone():
        """
        Retrieve the luminosity at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of"
        )
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of"
        )
        function.addParameter(
            'lum_i', dtype='float64', direction=function.OUT,
            description=(
                "The luminosity at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_mass_of_species():
        """
        Retrieve the mass number of the chemical abundance variable of the
        star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'species', dtype='int32', direction=function.IN,
            description="The species of the star to get the mass number of")
        function.addParameter(
            'species_mass', dtype='float64', direction=function.OUT,
            description=(
                "The mass number of the chemical abundance variable of "
                "the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function

    @legacy_function
    def get_mass_fraction_at_zone():
        """
        Retrieve the mass fraction at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'dq_i', dtype='float64', direction=function.OUT,
            description=(
                "The mass fraction at the specified zone/mesh-cell of "
                "the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_h_at_zone():
        """
        Retrieve the fractional abundance of h at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_he3_at_zone():
        """
        Retrieve the fractional abundance of he3 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_he_at_zone():
        """
        Retrieve the fractional abundance of he at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_c12_at_zone():
        """
        Retrieve the fractional abundance of c12 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_c13_at_zone():
        """
        Retrieve the fractional abundance of c13 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_n14_at_zone():
        """
        Retrieve the fractional abundance of n14 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_n15_at_zone():
        """
        Retrieve the fractional abundance of n15 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_o16_at_zone():
        """
        Retrieve the fractional abundance of o16 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_o17_at_zone():
        """
        Retrieve the fractional abundance of o17 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_o18_at_zone():
        """
        Retrieve the fractional abundance of o18 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_ne20_at_zone():
        """
        Retrieve the fractional abundance of ne20 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_ne22_at_zone():
        """
        Retrieve the fractional abundance of ne22 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_mg24_at_zone():
        """
        Retrieve the fractional abundance of mg24 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_mg25_at_zone():
        """
        Retrieve the fractional abundance of mg25 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_mg26_at_zone():
        """
        Retrieve the fractional abundance of mg26 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_c14_at_zone():
        """
        Retrieve the fractional abundance of c14 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_f18_at_zone():
        """
        Retrieve the fractional abundance of f18 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_f19_at_zone():
        """
        Retrieve the fractional abundance of f19 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_ne21_at_zone():
        """
        Retrieve the fractional abundance of ne21 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_na23_at_zone():
        """
        Retrieve the fractional abundance of na23 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_al26_at_zone():
        """
        Retrieve the fractional abundance of al26 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_al27_at_zone():
        """
        Retrieve the fractional abundance of al27 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_si28_at_zone():
        """
        Retrieve the fractional abundance of si28 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_neut_at_zone():
        """
        Retrieve the fractional abundance of neut at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_prot_at_zone():
        """
        Retrieve the fractional abundance of prot at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_bid_at_zone():
        """
        Retrieve the fractional abundance of bid at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_bid1_at_zone():
        """
        Retrieve the fractional abundance of bid1 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


class Genec(StellarEvolution, InternalStellarStructure):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(
            self,  GenecInterface(**options), **options)
        # self.stopping_conditions = StoppingConditions(self)
        self.model_time = 0.0 | units.yr

    def define_parameters(self, handler):
        handler.add_method_parameter(
            None,
            "set_genec_path",
            "path_to_data",
            "Path to the data directory",
            default_value=self.data_directory
        )

        handler.add_method_parameter(
            "get_par_ipoly",
            "set_par_ipoly",
            "ipoly",
            "GENEC parameter ipoly",
        )

        handler.add_method_parameter(
            "get_par_nwseq",
            "set_par_nwseq",
            "nwseq",
            "GENEC parameter nwseq",
        )

        handler.add_method_parameter(
            "get_par_modanf",
            "set_par_modanf",
            "modanf",
            "GENEC parameter modanf",
        )

        handler.add_method_parameter(
            "get_par_nzmod",
            "set_par_nzmod",
            "nzmod",
            "GENEC parameter nzmod",
        )

        handler.add_method_parameter(
            "get_par_irot",
            "set_par_irot",
            "irot",
            "GENEC parameter irot",
        )

        handler.add_method_parameter(
            "get_par_isol",
            "set_par_isol",
            "isol",
            "GENEC parameter isol",
        )

        handler.add_method_parameter(
            "get_par_imagn",
            "set_par_imagn",
            "imagn",
            "GENEC parameter imagn",
        )

        handler.add_method_parameter(
            "get_par_ialflu",
            "set_par_ialflu",
            "ialflu",
            "GENEC parameter ialflu",
        )

        handler.add_method_parameter(
            "get_par_ianiso",
            "set_par_ianiso",
            "ianiso",
            "GENEC parameter ianiso",
        )

        handler.add_method_parameter(
            "get_par_ipop3",
            "set_par_ipop3",
            "ipop3",
            "GENEC parameter ipop3",
        )

        handler.add_method_parameter(
            "get_par_ibasnet",
            "set_par_ibasnet",
            "ibasnet",
            "GENEC parameter ibasnet",
        )

        handler.add_method_parameter(
            "get_par_phase",
            "set_par_phase",
            "phase",
            "GENEC parameter phase",
        )

        handler.add_method_parameter(
            "get_par_iopac",
            "set_par_iopac",
            "iopac",
            "GENEC parameter iopac",
        )

        handler.add_method_parameter(
            "get_par_ikappa",
            "set_par_ikappa",
            "ikappa",
            "GENEC parameter ikappa",
        )

        handler.add_method_parameter(
            "get_par_idiff",
            "set_par_idiff",
            "idiff",
            "GENEC parameter idiff",
        )

        handler.add_method_parameter(
            "get_par_iadvec",
            "set_par_iadvec",
            "iadvec",
            "GENEC parameter iadvec",
        )

        handler.add_method_parameter(
            "get_par_istati",
            "set_par_istati",
            "istati",
            "GENEC parameter istati",
        )

        handler.add_method_parameter(
            "get_par_icoeff",
            "set_par_icoeff",
            "icoeff",
            "GENEC parameter icoeff",
        )

        handler.add_method_parameter(
            "get_par_igamma",
            "set_par_igamma",
            "igamma",
            "GENEC parameter igamma",
        )

        handler.add_method_parameter(
            "get_par_idialo",
            "set_par_idialo",
            "idialo",
            "GENEC parameter idialo",
        )

        handler.add_method_parameter(
            "get_par_idialu",
            "set_par_idialu",
            "idialu",
            "GENEC parameter idialu",
        )

        handler.add_method_parameter(
            "get_par_imloss",
            "set_par_imloss",
            "imloss",
            "GENEC parameter imloss",
        )

        handler.add_method_parameter(
            "get_par_ifitm",
            "set_par_ifitm",
            "ifitm",
            "GENEC parameter ifitm",
        )

        handler.add_method_parameter(
            "get_par_nndr",
            "set_par_nndr",
            "nndr",
            "GENEC parameter nndr",
        )

        handler.add_method_parameter(
            "get_par_iledou",
            "set_par_iledou",
            "iledou",
            "GENEC parameter iledou",
        )

        handler.add_method_parameter(
            "get_par_idifcon",
            "set_par_idifcon",
            "idifcon",
            "GENEC parameter idifcon",
        )

        handler.add_method_parameter(
            "get_par_my",
            "set_par_my",
            "my",
            "GENEC parameter my",
        )

        handler.add_method_parameter(
            "get_par_iover",
            "set_par_iover",
            "iover",
            "GENEC parameter iover",
        )

        handler.add_method_parameter(
            "get_par_iunder",
            "set_par_iunder",
            "iunder",
            "GENEC parameter iunder",
        )

        handler.add_method_parameter(
            "get_par_nbchx",
            "set_par_nbchx",
            "nbchx",
            "GENEC parameter nbchx",
        )

        handler.add_method_parameter(
            "get_par_nrband",
            "set_par_nrband",
            "nrband",
            "GENEC parameter nrband",
        )

        handler.add_method_parameter(
            "get_par_islow",
            "set_par_islow",
            "islow",
            "GENEC parameter islow",
        )

        handler.add_method_parameter(
            "get_par_icncst",
            "set_par_icncst",
            "icncst",
            "GENEC parameter icncst",
        )

        handler.add_method_parameter(
            "get_par_tauH_fit",
            "set_par_tauH_fit",
            "tauH_fit",
            "GENEC parameter tauH_fit",
        )

        handler.add_method_parameter(
            "get_par_iauto",
            "set_par_iauto",
            "iauto",
            "GENEC parameter iauto",
        )

        handler.add_method_parameter(
            "get_par_iprn",
            "set_par_iprn",
            "iprn",
            "GENEC parameter iprn",
        )

        handler.add_method_parameter(
            "get_par_iout",
            "set_par_iout",
            "iout",
            "GENEC parameter iout",
        )

        handler.add_method_parameter(
            "get_par_itmin",
            "set_par_itmin",
            "itmin",
            "GENEC parameter itmin",
        )

        handler.add_method_parameter(
            "get_par_idebug",
            "set_par_idebug",
            "idebug",
            "GENEC parameter idebug",
        )

        handler.add_method_parameter(
            "get_par_itests",
            "set_par_itests",
            "itests",
            "GENEC parameter itests",
        )

        handler.add_method_parameter(
            "get_par_var_rates",
            "set_par_var_rates",
            "var_rates",
            "GENEC parameter var_rates",
        )

        handler.add_method_parameter(
            "get_par_bintide",
            "set_par_bintide",
            "bintide",
            "GENEC parameter bintide",
        )

        handler.add_method_parameter(
            "get_par_const_per",
            "set_par_const_per",
            "const_per",
            "GENEC parameter const_per",
        )

        handler.add_method_parameter(
            "get_par_Add_Flux",
            "set_par_Add_Flux",
            "Add_Flux",
            "GENEC parameter Add_Flux",
        )

        handler.add_method_parameter(
            "get_par_diff_only",
            "set_par_diff_only",
            "diff_only",
            "GENEC parameter diff_only",
        )

        handler.add_method_parameter(
            "get_par_lowRSGMdot",
            "set_par_lowRSGMdot",
            "lowRSGMdot",
            "GENEC parameter lowRSGMdot",
        )

        handler.add_method_parameter(
            "get_par_plot",
            "set_par_plot",
            "plot",
            "GENEC parameter plot",
        )

        handler.add_method_parameter(
            "get_par_refresh",
            "set_par_refresh",
            "refresh",
            "GENEC parameter refresh",
        )

        handler.add_method_parameter(
            "get_par_xyfiles",
            "set_par_xyfiles",
            "xyfiles",
            "GENEC parameter xyfiles",
        )

        handler.add_method_parameter(
            "get_par_verbose",
            "set_par_verbose",
            "verbose",
            "GENEC parameter verbose",
        )

        handler.add_method_parameter(
            "get_par_stop_deg",
            "set_par_stop_deg",
            "stop_deg",
            "GENEC parameter stop_deg",
        )

        handler.add_method_parameter(
            "get_par_index_poly",
            "set_par_index_poly",
            "index_poly",
            "GENEC parameter index_poly",
        )

        handler.add_method_parameter(
            "get_par_binm2",
            "set_par_binm2",
            "binm2",
            "GENEC parameter binm2",
        )

        handler.add_method_parameter(
            "get_par_periodini",
            "set_par_periodini",
            "periodini",
            "GENEC parameter periodini",
        )

        handler.add_method_parameter(
            "get_par_zinit",
            "set_par_zinit",
            "zinit",
            "GENEC parameter zinit",
        )

        handler.add_method_parameter(
            "get_par_zsol",
            "set_par_zsol",
            "zsol",
            "GENEC parameter zsol",
        )

        handler.add_method_parameter(
            "get_par_z",
            "set_par_z",
            "z",
            "GENEC parameter z",
        )

        handler.add_method_parameter(
            "get_par_fenerg",
            "set_par_fenerg",
            "fenerg",
            "GENEC parameter fenerg",
        )

        handler.add_method_parameter(
            "get_par_richac",
            "set_par_richac",
            "richac",
            "GENEC parameter richac",
        )

        handler.add_method_parameter(
            "get_par_frein",
            "set_par_frein",
            "frein",
            "GENEC parameter frein",
        )

        handler.add_method_parameter(
            "get_par_K_Kawaler",
            "set_par_K_Kawaler",
            "K_Kawaler",
            "GENEC parameter K_Kawaler",
        )

        handler.add_method_parameter(
            "get_par_Omega_saturation",
            "set_par_Omega_saturation",
            "Omega_saturation",
            "GENEC parameter Omega_saturation",
        )

        handler.add_method_parameter(
            "get_par_rapcrilim",
            "set_par_rapcrilim",
            "rapcrilim",
            "GENEC parameter rapcrilim",
        )

        handler.add_method_parameter(
            "get_par_vwant",
            "set_par_vwant",
            "vwant",
            "GENEC parameter vwant",
        )

        handler.add_method_parameter(
            "get_par_xfom",
            "set_par_xfom",
            "xfom",
            "GENEC parameter xfom",
        )

        handler.add_method_parameter(
            "get_par_omega",
            "set_par_omega",
            "omega",
            "GENEC parameter omega",
        )

        handler.add_method_parameter(
            "get_par_xdial",
            "set_par_xdial",
            "xdial",
            "GENEC parameter xdial",
        )

        handler.add_method_parameter(
            "get_par_B_initial",
            "set_par_B_initial",
            "B_initial",
            "GENEC parameter B_initial",
        )

        handler.add_method_parameter(
            "get_par_add_diff",
            "set_par_add_diff",
            "add_diff",
            "GENEC parameter add_diff",
        )

        handler.add_method_parameter(
            "get_par_fmlos",
            "set_par_fmlos",
            "fmlos",
            "GENEC parameter fmlos",
        )

        handler.add_method_parameter(
            "get_par_fitm",
            "set_par_fitm",
            "fitm",
            "GENEC parameter fitm",
        )

        handler.add_method_parameter(
            "get_par_fitmi",
            "set_par_fitmi",
            "fitmi",
            "GENEC parameter fitmi",
        )

        handler.add_method_parameter(
            "get_par_fitmi_default",
            "set_par_fitmi_default",
            "fitmi_default",
            "GENEC parameter fitmi_default",
        )

        handler.add_method_parameter(
            "get_par_deltal",
            "set_par_deltal",
            "deltal",
            "GENEC parameter deltal",
        )

        handler.add_method_parameter(
            "get_par_deltat",
            "set_par_deltat",
            "deltat",
            "GENEC parameter deltat",
        )

        handler.add_method_parameter(
            "get_par_elph",
            "set_par_elph",
            "elph",
            "GENEC parameter elph",
        )

        handler.add_method_parameter(
            "get_par_dovhp",
            "set_par_dovhp",
            "dovhp",
            "GENEC parameter dovhp",
        )

        handler.add_method_parameter(
            "get_par_dunder",
            "set_par_dunder",
            "dunder",
            "GENEC parameter dunder",
        )

        handler.add_method_parameter(
            "get_par_gkorm",
            "set_par_gkorm",
            "gkorm",
            "GENEC parameter gkorm",
        )

        handler.add_method_parameter(
            "get_par_alph",
            "set_par_alph",
            "alph",
            "GENEC parameter alph",
        )

        handler.add_method_parameter(
            "get_par_agdr",
            "set_par_agdr",
            "agdr",
            "GENEC parameter agdr",
        )

        handler.add_method_parameter(
            "get_par_faktor",
            "set_par_faktor",
            "faktor",
            "GENEC parameter faktor",
        )

        handler.add_method_parameter(
            "get_par_dgrp",
            "set_par_dgrp",
            "dgrp",
            "GENEC parameter dgrp",
        )

        handler.add_method_parameter(
            "get_par_dgrl",
            "set_par_dgrl",
            "dgrl",
            "GENEC parameter dgrl",
        )

        handler.add_method_parameter(
            "get_par_dgry",
            "set_par_dgry",
            "dgry",
            "GENEC parameter dgry",
        )

        handler.add_method_parameter(
            "get_par_dgrc",
            "set_par_dgrc",
            "dgrc",
            "GENEC parameter dgrc",
        )

        handler.add_method_parameter(
            "get_par_dgro",
            "set_par_dgro",
            "dgro",
            "GENEC parameter dgro",
        )

        handler.add_method_parameter(
            "get_par_dgr20",
            "set_par_dgr20",
            "dgr20",
            "GENEC parameter dgr20",
        )

        handler.add_method_parameter(
            "get_par_xcn",
            "set_par_xcn",
            "xcn",
            "GENEC parameter xcn",
        )

        handler.add_method_parameter(
            "get_par_starname",
            "set_par_starname",
            "starname",
            "GENEC parameter starname",
        )

        handler.add_method_parameter(
            "get_min_timestep_stop_condition",
            "set_min_timestep_stop_condition",
            "min_timestep_stop_condition",
            "The minimum timestep stop condition",
            default_value=1e-4 | units.julianyr
        )

    def define_particle_sets(self, handler):

        # for particle_set in ['particles', 'fullparticles']:
        for set_name in ['fullparticles']:
            handler.define_set(set_name, 'index_of_the_star')
            InternalStellarStructure.define_particle_sets(
                self, handler, set_name=set_name
            )
            handler.set_new(set_name, 'new_particle')

            handler.add_getter(set_name, 'get_radius')
            handler.add_getter(set_name, 'get_mass')
            handler.add_getter(set_name, 'get_age')
            handler.add_getter(set_name, 'get_luminosity')
            handler.add_getter(set_name, 'get_temperature')
            handler.add_getter(set_name, 'get_time_step', names=('time_step',))
            handler.add_getter(
                set_name, 'get_number_of_zones', names=('n_zones',)
            )
            handler.add_getter(
                set_name, 'get_number_of_species', names=('n_species',)
            )

            # handler.add_method(set_name, 'get_number_of_zones')
            handler.add_method(set_name, 'get_number_of_zones')

            # handler.add_method(set_name, 'get_radius_profile')
            # handler.add_method(set_name, 'get_temperature_profile')
            handler.add_method(set_name, 'get_luminosity_profile')
            handler.add_method(set_name, 'get_mass_profile')
            handler.add_method(set_name, 'get_cumulative_mass_profile')

            handler.add_method(set_name, 'evolve_one_step')
            handler.add_method(set_name, 'evolve_for')
            handler.set_delete(set_name, 'delete_star')

            handler.add_gridded_getter(
                set_name,
                'get_radius_at_zone', 'get_firstlast_zone',
                names=('radius_profile',)
            )
            handler.add_gridded_setter(
                set_name,
                'set_radius_at_zone', 'get_firstlast_zone',
                names=('radius_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_temperature_at_zone', 'get_firstlast_zone',
                names=('temperature_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_density_at_zone', 'get_firstlast_zone',
                names=('density_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_luminosity_at_zone', 'get_firstlast_zone',
                names=('luminosity_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_pressure_at_zone', 'get_firstlast_zone',
                names=('pressure_profile',)
            )


            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_h_at_zone',
                'get_firstlast_zone',
                names=('abundance_h',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_he3_at_zone',
                'get_firstlast_zone',
                names=('abundance_he3',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_he_at_zone',
                'get_firstlast_zone',
                names=('abundance_he',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_c12_at_zone',
                'get_firstlast_zone',
                names=('abundance_c12',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_c13_at_zone',
                'get_firstlast_zone',
                names=('abundance_c13',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_n14_at_zone',
                'get_firstlast_zone',
                names=('abundance_n14',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_n15_at_zone',
                'get_firstlast_zone',
                names=('abundance_n15',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_o16_at_zone',
                'get_firstlast_zone',
                names=('abundance_o16',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_o17_at_zone',
                'get_firstlast_zone',
                names=('abundance_o17',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_o18_at_zone',
                'get_firstlast_zone',
                names=('abundance_o18',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_ne20_at_zone',
                'get_firstlast_zone',
                names=('abundance_ne20',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_ne22_at_zone',
                'get_firstlast_zone',
                names=('abundance_ne22',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_mg24_at_zone',
                'get_firstlast_zone',
                names=('abundance_mg24',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_mg25_at_zone',
                'get_firstlast_zone',
                names=('abundance_mg25',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_mg26_at_zone',
                'get_firstlast_zone',
                names=('abundance_mg26',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_c14_at_zone',
                'get_firstlast_zone',
                names=('abundance_c14',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_f18_at_zone',
                'get_firstlast_zone',
                names=('abundance_f18',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_f19_at_zone',
                'get_firstlast_zone',
                names=('abundance_f19',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_ne21_at_zone',
                'get_firstlast_zone',
                names=('abundance_ne21',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_na23_at_zone',
                'get_firstlast_zone',
                names=('abundance_na23',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_al26_at_zone',
                'get_firstlast_zone',
                names=('abundance_al26',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_al27_at_zone',
                'get_firstlast_zone',
                names=('abundance_al27',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_si28_at_zone',
                'get_firstlast_zone',
                names=('abundance_si28',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_neut_at_zone',
                'get_firstlast_zone',
                names=('abundance_neut',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_prot_at_zone',
                'get_firstlast_zone',
                names=('abundance_prot',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_bid_at_zone',
                'get_firstlast_zone',
                names=('abundance_bid',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_mass_fraction_of_bid1_at_zone',
                'get_firstlast_zone',
                names=('abundance_bid1',)
            )

            # for species in self.SPECIES_NAMES:
            #     getter = f'get_{species}_at_zone'
            #     handler.add_gridded_getter(
            #         particle_set,
            #         getter,
            #         'get_firstlast_zone',
            #     )
            # for species in species_name:
            #     handler.add_gridded_getter(
            #         particle_set,
            #         'get_species_at_zone', 'get_firstlast_zone',
            #     )
        # handler.add_gridded_getter(
        #     'particle'
        # )

    @property
    def particles(self):
        basic_attributes = ["age", "mass", "radius", "temperature", "luminosity",]
        return ParticlesWithFilteredAttributes(
            self.fullparticles,
            basic_attributes,
        )

    def define_state(self, handler):
        StellarEvolution.define_state(self, handler)
        # InternalStellarStructure.define_state(self, handler)

        # Only allow setting of starname in EDIT or UPDATE states
        # I.e. must do initialize_code and commit_parameters FIRST!

        # Initialized (initialize_code)
        handler.add_method

        # -> Edit (commit_parameters)
        handler.add_method('EDIT', 'set_starname')
        # handler.add_method('EDIT', 'new_particle')

        # -> Run (commit_particles)
        handler.add_transition('EDIT', 'RUN', 'commit_particles')

        handler.add_method('RUN', 'get_chemical_abundance_profiles')
        handler.add_method('RUN', 'get_mass_fraction_of_species_at_zone')
        handler.add_method('RUN', 'get_mu_at_zone')
        handler.add_method('RUN', 'get_number_of_species')
        handler.add_method('RUN', 'get_number_of_zones')
        handler.add_method('RUN', 'get_pressure_at_zone')
        handler.add_method('RUN', 'get_radius')
        handler.add_method('RUN', 'get_radius_at_zone')
        handler.add_method('RUN', 'get_temperature_at_zone')
        handler.add_method('RUN', 'get_density_at_zone')
        handler.add_method('RUN', 'get_luminosity_at_zone')
        handler.add_method('RUN', 'get_mass_fraction_at_zone')
        handler.add_method('RUN', 'get_time_step')

        # -> Update
        handler.add_transition('RUN', 'UPDATE', 'finalize_stellar_model')
        handler.add_method('UPDATE', 'set_chemical_abundance_profiles')
        handler.add_method('UPDATE', 'set_mass_fraction_of_species_at_zone')
        handler.add_method('UPDATE', 'set_number_of_zones')
        # handler.add_method('UPDATE', 'set_radius')
        handler.add_method('UPDATE', 'set_radius_at_zone')
        handler.add_method('UPDATE', 'set_temperature_at_zone')
        handler.add_method('UPDATE', 'set_density_at_zone')
        handler.add_method('UPDATE', 'set_luminosity_at_zone')
        handler.add_method('UPDATE', 'set_mass_fraction_at_zone')
        # handler.add_method('UPDATE', 'set_time_step')
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        # -> Run (recommit_particles)

        handler.add_method('UPDATE', 'set_starname')
        # handler.add_method('UPDATE', 'new_particle')

    def define_methods(self, handler):
        InternalStellarStructure.define_methods(self, handler)
        StellarEvolution.define_methods(self, handler)
        handler.add_method(
            "new_particle",
            (units.MSun, handler.NO_UNIT, handler.NO_UNIT),
            (handler.INDEX, handler.ERROR_CODE)
        )
        # handler.add_method(
        #     "get_radius",
        #     (handler.INDEX,),
        #     (units.RSun, handler.ERROR_CODE,)
        # )
        # handler.add_method(
        #     "get_number_of_zones",
        #     (handler.INDEX,),
        #     (handler.NO_UNIT, handler.ERROR_CODE,)
        # )
        # handler.add_method(
        #     "get_radius_at_zone",
        #     (handler.INDEX, handler.NO_UNIT,),
        #     (units.cm, handler.ERROR_CODE,)
        # )
        # handler.add_method(
        #     "set_radius_at_zone",
        #     (handler.INDEX, handler.NO_UNIT, units.cm),
        #     (handler.ERROR_CODE,)
        # )
        # handler.add_method(
        #     "get_temperature_at_zone",
        #     (handler.INDEX, handler.NO_UNIT,),
        #     (units.K, handler.ERROR_CODE,)
        # )
        # handler.add_method(
        #     "set_temperature_at_zone",
        #     (handler.INDEX, handler.NO_UNIT, units.K),
        #     (handler.ERROR_CODE,)
        # )
        handler.add_method(
            "get_luminosity_at_zone",
            (handler.INDEX, handler.NO_UNIT,),
            (units.erg/units.s, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_luminosity_at_zone",
            (handler.INDEX, handler.NO_UNIT,units.erg/units.s),
            (handler.ERROR_CODE,)
        )
    # def define_parameters(self, handler):

    def get_luminosity_profile(
            self,
            indices_of_the_stars,
            number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying luminosity profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_luminosity_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_mass_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying mass profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_mass_fraction_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_cumulative_mass_profile(
            self, indices_of_the_stars, number_of_zones=None):
        frac_profile = self.get_mass_profile(
            indices_of_the_stars, number_of_zones=number_of_zones)
        return frac_profile.cumsum()
