from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class Nbody6xxInterface(
        CodeInterface,
        GravitationalDynamicsInterface,
        GravityFieldInterface
    ):
    use_modules = ['AMUSE_INTERFACE']

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="nbody6xx_worker", **keyword_arguments)

    @legacy_function
    def main():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def set_kz():
        function = LegacyFunctionSpecification()
        function.addParameter('kz_option', dtype='int32', direction=function.IN)
        function.addParameter('val', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def get_kz():
        function = LegacyFunctionSpecification()
        function.addParameter('kz_option', dtype='int32', direction=function.IN)
        function.addParameter('val', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def set_eta():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_etai():
        """
        Set the current irregular force time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('etai', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_etai():
        """
        Get the current system irregular force time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('etai', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_etar():
        """
        Set the current regular force time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('etar', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_etar():
        """
        Get the current system regular force time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('etar', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rbar():
        """
        Set the scaling unit in parsec for one N-body unit of length.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64',
                              direction=function.IN, unit=units.parsec)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_rbar():
        """
        Get the scaling unit in parsec for one N-body unit of length.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64',
                              direction=function.OUT, unit=units.parsec)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_zmbar():
        """
        Set the scaling unit for average particle mass in solar masses.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64',
                              direction=function.IN,unit=units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_zmbar():
        """
        Get the scaling unit for average particle mass in solar masses.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64',
                              direction=function.OUT,unit=units.MSun)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_qe():
        """
        Set the relative energy error tolerance.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_qe():
        """
        Get the relative energy error tolerance.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function



class Nbody6xx(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody = None, **kargs):
        GravitationalDynamics.__init__(self,  Nbody6xxInterface(**kargs), convert_nbody, **kargs)
#        RBAR= virial radius in parsec
#        ZMBAR= mean particle mass in solar units
#        if convert_nbody is not None:
#          self.parameters.RBAR=convert_nbody.to_si(1 | nbody_system.length)
#          self.parameters.ZMBAR=convert_nbody.to_si(1 | nbody_system.mass)

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        handler.add_method('RUN', 'get_particle_timestep')
        GravityFieldCode.define_state(self, handler)

        handler.add_method('EDIT', 'set_state')
        handler.add_method('EDIT', 'set_velocity')
        handler.add_method('EDIT', 'set_mass')
        handler.add_method('EDIT', 'set_position')
        handler.add_method('CHANGED','before_get_parameter')

        handler.add_transition('RUN', 'CHANGED', 'set_state', False)
        handler.add_transition('RUN', 'CHANGED', 'set_velocity', False)
        handler.add_transition('RUN', 'CHANGED', 'set_mass', False)
        handler.add_transition('RUN', 'CHANGED', 'set_position', False)
        handler.add_transition('CHANGED', 'RUN', 'synchronize_model')
        handler.add_method('CHANGED', 'get_state')
        handler.add_method('CHANGED', 'get_mass')
        handler.add_method('CHANGED', 'get_position')
        handler.add_method('CHANGED', 'get_velocity')
        handler.add_method('CHANGED', 'get_particle_timestep')

    def define_parameters(self, handler):

        # Set/get parameters specific to the module, not part of the
        # standard interface.  Accessors used here must be defined
        # above and reflected in interface.cc.  Python access is
        # (e.g.)
        #
        #        ph4.parameters.timestep_parameter = xxx

        handler.add_method_parameter(
            "get_eta",                   # getter name in interface.cc
            "set_eta",                   # setter name in interface.cc
            "timestep_parameter",        # python parameter name
            "timestep parameter",        # description
            default_value = 0.05
        )

        handler.add_method_parameter(
            "get_eps2",                  # already defined in standard interface
            "set_eps2",                  # already defined in standard interface
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )


        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )

        handler.add_method_parameter(
            "get_qe",
            "set_qe",
            "QE",
            "Relative energy error tolerance",
            default_value = 2.0E-4
        )

        handler.add_method_parameter(
            "get_rbar",
            "set_rbar",
            "RBAR",
            "The scaling unit in parsec for one N-body unit of length",
            default_value = 1.0 | units.parsec
        )

        handler.add_method_parameter(
            "get_zmbar",
            "set_zmbar",
            "ZMBAR",
            "The scaling unit for average particle mass in solar masses",
            default_value = 0.7 | units.MSun
        )

        #self.stopping_conditions.define_parameters(handler)

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)

        # Turn interface functions into methods.

        handler.add_method(
            "set_eps2",
            (
                nbody_system.length * nbody_system.length
            ),
            (
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "get_eps2",
            (),
            (
                nbody_system.length * nbody_system.length,
                handler.ERROR_CODE
            )
        )

