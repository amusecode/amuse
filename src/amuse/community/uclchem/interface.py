from amuse.community.interface.chem import ChemicalModelingInterface
from amuse.community.interface.chem import ChemicalModeling
from amuse.community.interface import common
from amuse.community import *

class UclchemInterface(CodeInterface):
    def __init__(self, mode = 'cpu', **options):
        CodeInterface.__init__(
            self,
            name_of_the_worker='uclchem_worker',
            **options
        )
    @remote_function(can_handle_array=True)
    def sim_cloud(outSpeciesIn='s', dictionary='s'):
        returns (abundance='d')

class Uclchem(ChemicalModeling):
    def __init__(self, convert_nbody=None, **options):
        legacy_interface = UclchemInterface(**options)

        ChemicalModeling.__init__(self,legacy_interface)
    
    def define_methods(self, handler):
        common.CommonCode.define_methods(self, handler)
        handler.add_method(
            'sim_cloud',
            (handler.NO_UNIT, handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE)
        )
