from amuse.legacy import *
from amuse.legacy.interface.common import CommonCodeInterface

class MpiAmrVacInterface(LegacyInterface, CommonCodeInterface):
    
    use_modules = ['mpiamrvac_interface']
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, name_of_the_worker="mpiamrvac_worker", **keyword_arguments)
    
    
        
    
class MpiAmrVac(CodeInterface):

    def __init__(self, **options):
        CodeInterface.__init__(self,  MpiAmrVacInterface(**options), **options)
    
