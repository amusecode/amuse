from amuse.community.interface.chem import ChemicalModelingInterface
from amuse.community.interface.chem import ChemicalModeling
from amuse.community import *

class uclcheminterface(CodeInterface, ChemicalModelingInterface, StoppingConditionInterface):
    def make_dict(*options):
        dictionary = dict(options)
        return dictionary