from amuse.rfi.core import CodeInterface, remote_function
from amuse.units import units


class ForTestingInterface(CodeInterface):
    include_headers = ['c_interface.h']

    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    @remote_function
    def set0(a=0. | units.m):
        returns()

    @remote_function
    def get0():
        returns(a=0. | units.m)

    @remote_function
    def get_grid0_range():
        returns()
