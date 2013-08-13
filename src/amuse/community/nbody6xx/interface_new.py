from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class Nbody6xxInterface(CodeInterface,
        GravitationalDynamicsInterface,
        GravitationalDynamics,
        GravityFieldInterface,
        GravityFieldCode):

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="nbody6xx_worker", **keyword_arguments)

    @legacy_function
    def main():
        function = LegacyFunctionSpecification()
        #function.addParameter('int_in', dtype='int32', direction=function.IN)
        #function.addParameter('int_out', dtype='int32', direction=function.OUT)
        #function.result_type = 'int32'
        function.can_handle_array = True
        return function


class Nbody6xx(InCodeComponentImplementation):

    def __init__(self):
        InCodeComponentImplementation.__init__(self,  nbody6xxInterface())

