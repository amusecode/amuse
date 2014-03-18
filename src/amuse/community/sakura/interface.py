from amuse.community.interface.gd import GravitationalDynamics
from amuse.rfi.core import PythonCodeInterface

from amuse.community.tupan.interface import TupanImplementation, TupanInterface, Tupan, MODULES_MISSING


class SakuraImplementation(TupanImplementation):

    def __init__(self):
        super(SakuraImplementation, self).__init__()
        self.integrator_method = "sakura"


class SakuraInterface(TupanInterface):

    def __init__(self, **options):
        PythonCodeInterface.__init__(
            self,
            SakuraImplementation,
            'sakura_worker',
            **options)


class Sakura(Tupan):

    def __init__(self, convert_nbody=None, **options):
        nbody_interface = SakuraInterface(**options)

        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        super(Sakura, self).define_parameters(object)

        object.add_method_parameter(
            "get_integrator_method",
            "set_integrator_method",
            "integrator_method",
            "The method to use to integrate the evolution \
             of the system (choices: [sakura, asakura]) \
             Note: 'asakura' is sakura with adaptive time-steps.",
            default_value="sakura"
        )





### end of file ###
