from amuse.community import *
from amuse.community.interface.gd import *
from amuse.community.interface.stopping_conditions import *


class HermitegrxInterface(
    CodeInterface, GravitationalDynamicsInterface, StoppingConditionInterface
):
    include_headers = ["worker_code.h", "stopcond.h"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="hermite_grx_worker", **keyword_arguments
        )
        self.initialize_code()

    def reinitialize_particles(self):
        self.recommit_particles()

    @legacy_function
    def get_dt_param():
        """
        Get the timestep scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "dt_param",
            dtype="float64",
            direction=function.OUT,
            description="the timestep scaling factor",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def set_dt_param():
        """
        Set the timestep scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "dt_param",
            dtype="float64",
            direction=function.IN,
            description="the timestep scaling factor",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_time():
        """
        Get the current simulation time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "time",
            dtype="float64",
            direction=function.OUT,
            description="the current simulation time",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def set_time():
        """
        Set the current simulation time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "time",
            dtype="float64",
            direction=function.IN,
            description="the current simulation time",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_light_speed():
        """
        Get the speed of light.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "lightspeed",
            dtype="float64",
            direction=function.OUT,
            description="the current speed of light",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def set_light_speed():
        """
	    Set the speed of light.
	    """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "lightspeed",
            dtype="float64",
            direction=function.IN,
            description="the current speed of light",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_integrator():
        """
        Get the current integrator.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "integrator",
            dtype="string",
            direction=function.OUT,
            description="the current integrator",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def set_integrator():
        """
        Set the current force calculator.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "integrator",
            dtype="string",
            direction=function.IN,
            description="the current integrator",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_perturbation():
        """
        Get the current perturbation.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "interaction",
            dtype="string",
            direction=function.OUT,
            description="the current perturbation",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def set_perturbation():
        """
        Set the current perturbation.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "interaction",
            dtype="string",
            direction=function.IN,
            description="the current perturbation",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_num_threads():
        """
        Get the number of threads.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "num_threads",
            dtype="int32",
            direction=function.OUT,
            description="the number of threads",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def set_num_threads():
        """
        Set the number of threads.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "num_threads",
            dtype="int32",
            direction=function.IN,
            description="the number of threads",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_total_energy_with():
        """
        Get the total energy.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "type",
            dtype="string",
            direction=function.IN,
            description="the perturbation with which to calculate",
        )
        function.addParameter(
            "etot",
            dtype="float64",
            direction=function.OUT,
            description="the current total energy",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_total_linear_momentum_with():
        """
        Get the total linear momentum.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "type",
            dtype="string",
            direction=function.IN,
            description="the perturbation with which to calculate",
        )
        function.addParameter(
            "px",
            dtype="float64",
            direction=function.OUT,
            description="the momentum in x",
        )
        function.addParameter(
            "py",
            dtype="float64",
            direction=function.OUT,
            description="the momentum in y",
        )
        function.addParameter(
            "pz",
            dtype="float64",
            direction=function.OUT,
            description="the momentum in z",
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def new_large_particle():
        """
        Add a new black hole particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            "id", dtype="int32", direction=function.OUT, description="the identifier"
        )
        function.addParameter(
            "mass", dtype="float64", direction=function.IN, description="the mass"
        )
        function.addParameter(
            "x", dtype="float64", direction=function.IN, description="the x position"
        )
        function.addParameter(
            "y", dtype="float64", direction=function.IN, description="the y position"
        )
        function.addParameter(
            "z", dtype="float64", direction=function.IN, description="the z position"
        )
        function.addParameter(
            "vx", dtype="float64", direction=function.IN, description="the x velocity"
        )
        function.addParameter(
            "vy", dtype="float64", direction=function.IN, description="the y velocity"
        )
        function.addParameter(
            "vz", dtype="float64", direction=function.IN, description="the z velocity"
        )
        function.addParameter(
            "radius", dtype="float64", direction=function.IN, description="the radius"
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def new_small_particle():
        """
        Add a new star particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            "id", dtype="int32", direction=function.OUT, description="the identifier"
        )
        function.addParameter(
            "mass", dtype="float64", direction=function.IN, description="the mass"
        )
        function.addParameter(
            "x", dtype="float64", direction=function.IN, description="the x position"
        )
        function.addParameter(
            "y", dtype="float64", direction=function.IN, description="the y position"
        )
        function.addParameter(
            "z", dtype="float64", direction=function.IN, description="the z position"
        )
        function.addParameter(
            "vx", dtype="float64", direction=function.IN, description="the x velocity"
        )
        function.addParameter(
            "vy", dtype="float64", direction=function.IN, description="the y velocity"
        )
        function.addParameter(
            "vz", dtype="float64", direction=function.IN, description="the z velocity"
        )
        function.addParameter(
            "radius", dtype="float64", direction=function.IN, description="the radius"
        )
        function.result_type = "int32"
        return function

    @legacy_function
    def get_jerk():
        """
        Retrieve the jerk vector of a particle. Third time derivative
        of the position.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "index_of_the_particle",
            dtype="int32",
            direction=function.IN,
            description=(
                "Index of the particle to get the state from. This index must "
                "have been returned by an earlier call to :meth:`new_particle`"
            ),
        )
        function.addParameter(
            "jx",
            dtype="float64",
            direction=function.OUT,
            description="The current jerk vector of the particle",
        )
        function.addParameter(
            "jy",
            dtype="float64",
            direction=function.OUT,
            description="The current jerk vector of the particle",
        )
        function.addParameter(
            "jz",
            dtype="float64",
            direction=function.OUT,
            description="The current jerk vector of the particle",
        )
        function.result_type = "int32"
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def set_jerk():
        """
        Update the jerk of a particle.
        *Defined for symetry with the get_acceleration function.*
        *Should be removed if physaccily unsound*
        *Maybe moved to snapshot support functionality*
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            "index_of_the_particle",
            dtype="int32",
            direction=function.IN,
            description=(
                "Index of the particle for which the state is to be updated. "
                "This index must have been returned by an earlier call to "
                ":meth:`new_particle`"
            ),
        )
        function.addParameter(
            "ax",
            dtype="float64",
            direction=function.IN,
            description="The new jerk vector of the particle",
        )
        function.addParameter(
            "ay",
            dtype="float64",
            direction=function.IN,
            description="The new jerk vector of the particle",
        )
        function.addParameter(
            "az",
            dtype="float64",
            direction=function.IN,
            description="The new jerk vector of the particle",
        )
        function.result_type = "int32"
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function


class Hermitegrx(GravitationalDynamics):
    def __init__(self, convert_nbody=None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = HermitegrxInterface(**options)
        GravitationalDynamics.__init__(self, legacy_interface, convert_nbody, **options)

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value=0.0 | nbody_system.length ** 2,
        )
        object.add_method_parameter(
            "get_dt_param",
            "set_dt_param",
            "dt_param",
            "timestep scaling factor",
            default_value=0.03,
        )
        object.add_method_parameter(
            "get_light_speed",
            "set_light_speed",
            "light_speed",
            "the speed of light",
            default_value=1 | nbody_system.speed,
        )
        object.add_method_parameter(
            "get_integrator",
            "set_integrator",
            "integrator",
            "the integrator",
            default_value="Hermite",
        )
        object.add_method_parameter(
            "get_perturbation",
            "set_perturbation",
            "perturbation",
            "the perturbation",
            default_value="None",
        )
        object.add_method_parameter(
            "get_num_threads",
            "set_num_threads",
            "num_threads",
            "the number of threads",
            default_value=1,
        )

        self.stopping_conditions.define_parameters(object)

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        momentum = nbody_system.speed * nbody_system.mass

        object.add_method(
            "get_eps2", (), (nbody_system.length ** 2, object.ERROR_CODE,)
        )
        object.add_method("set_eps2", (nbody_system.length ** 2,), (object.ERROR_CODE,))
        object.add_method("get_dt_param", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_dt_param", (object.NO_UNIT,), (object.ERROR_CODE,))
        object.add_method("get_time", (), (nbody_system.time, object.ERROR_CODE,))
        object.add_method("set_time", (nbody_system.time,), (object.ERROR_CODE,))
        object.add_method(
            "get_light_speed", (), (nbody_system.speed, object.ERROR_CODE,)
        )
        object.add_method(
            "set_light_speed", (nbody_system.speed,), (object.ERROR_CODE,)
        )
        object.add_method(
            "get_total_energy_with",
            (object.NO_UNIT),
            (nbody_system.energy, object.ERROR_CODE,),
        )
        object.add_method(
            "get_total_linear_momentum_with",
            (object.NO_UNIT),
            (momentum, momentum, momentum, object.ERROR_CODE,),
        )
        object.add_method(
            "new_large_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
            ),
            (object.INDEX, object.ERROR_CODE),
        )
        object.add_method(
            "new_small_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
            ),
            (object.INDEX, object.ERROR_CODE),
        )
        object.add_method(
            "get_acceleration",
            (object.NO_UNIT,),
            (
                nbody_system.acceleration,
                nbody_system.acceleration,
                nbody_system.acceleration,
                object.ERROR_CODE,
            ),
        )
        object.add_method(
            "get_jerk",
            (object.NO_UNIT,),
            (
                nbody_system.length / nbody_system.time ** 3,
                nbody_system.length / nbody_system.time ** 3,
                nbody_system.length / nbody_system.time ** 3,
                object.ERROR_CODE,
            ),
        )

        self.stopping_conditions.define_methods(object)

    def define_particle_sets(self, object):
        object.define_super_set(
            "particles", ["large_particles", "small_particles"], index_to_default_set=1
        )

        object.define_set("large_particles", "index_of_the_particle")
        object.set_new("large_particles", "new_large_particle")
        object.set_delete("large_particles", "delete_particle")
        object.add_setter("large_particles", "set_state")
        object.add_getter("large_particles", "get_state")
        object.add_setter("large_particles", "set_mass")
        object.add_getter("large_particles", "get_mass", names=("mass",))
        object.add_setter("large_particles", "set_position")
        object.add_getter("large_particles", "get_position")
        object.add_setter("large_particles", "set_velocity")
        object.add_getter("large_particles", "get_velocity")
        object.add_setter("large_particles", "set_radius")
        object.add_getter("large_particles", "get_radius")
        object.add_getter("large_particles", "get_acceleration")
        object.add_getter("large_particles", "get_jerk")

        object.define_set("small_particles", "index_of_the_particle")
        object.set_new("small_particles", "new_small_particle")
        object.set_delete("small_particles", "delete_particle")
        object.add_setter("small_particles", "set_state")
        object.add_getter("small_particles", "get_state")
        object.add_setter("small_particles", "set_mass")
        object.add_getter("small_particles", "get_mass", names=("mass",))
        object.add_setter("small_particles", "set_position")
        object.add_getter("small_particles", "get_position")
        object.add_setter("small_particles", "set_velocity")
        object.add_getter("small_particles", "get_velocity")
        object.add_setter("small_particles", "set_radius")
        object.add_getter("small_particles", "get_radius")
        object.add_getter("small_particles", "get_acceleration")
        object.add_getter("small_particles", "get_jerk")

        self.stopping_conditions.define_particle_set(object)


# for backwards compatibility reasons
HermiteGRXInterface = HermitegrxInterface
HermiteGRX = Hermitegrx
