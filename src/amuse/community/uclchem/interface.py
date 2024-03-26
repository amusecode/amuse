from amuse.community.interface.chem import ChemicalModelingInterface
from amuse.community.interface.chem import ChemicalModeling
from amuse.community.interface import common
from amuse.community import *

class uclcheminterface(CodeInterface, ChemicalModelingInterface, StoppingConditionInterface):
    def __init__(self, mode = 'cpu', **options):
        CodeInterface.__init__(
            self,
            name_of_the_worker='uclchem_worker',
            **options
        )
    @remote_function
    def sim_cloud():
        returns ()

class uclchem(ChemicalModeling):
    def __init__(self, convert_nbody=None, **options):
        legacy_interface = uclcheminterface(**options)

        ChemicalModeling.__init__(legacy_interface)
    
    def define_methods(self, handler):
        common.CommonCode.define_methods(self, handler)
        handler.add_method(
            'sim_cloud',
            (),(handler.ERROR_CODE)
        )
