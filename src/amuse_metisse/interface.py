"""
Interface for metisse
"""

import os

from amuse.community import (
    CodeInterface,
    LiteratureReferencesMixIn,
    remote_function,
)
from amuse.community.interface import se

from amuse.units import units


# low level interface class
class MetisseInterface(
    CodeInterface,
    se.StellarEvolutionInterface,
    LiteratureReferencesMixIn,
):
    """
    Low level interface for METISSE

    Details in publication:
        .. [#] Agrawal, P. et al., 2020, https://doi.org/10.1093/mnras/staa2264
        .. [#] Agrawal, P. et al., 2023, https://doi.org/10.1093/mnras/stad2334
    """

    use_modules = ["metisseInterface"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="metisse_worker", **keyword_arguments
        )
        LiteratureReferencesMixIn.__init__(self)
        self.model_time = 0.0 | units.mega(units.julianyr)
        self.set_metallicity_dir(
            f"{os.path.dirname(__file__)}/data/sample_tracks_solarZ/Hydrogen"
        )
        self.set_metallicity_dir_he(
            f"{os.path.dirname(__file__)}/data/sample_tracks_solarZ/Helium"
        )

    # Remote functions - getters and setters
    # Note that we should maybe use SI units rather than derived (MSun etc), at
    # least while these are not certain to be the same in the code and in
    # AMUSE...
    @remote_function(can_handle_array=True)
    def get_initial_mass(index_of_the_star="i"):
        returns (mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_epoch(index_of_the_star="i"):
        returns (epoch="d" | units.mega(units.julianyr))

    @remote_function(can_handle_array=True)
    def get_core_mass(index_of_the_star="i"):
        returns (core_mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_core_radius(index_of_the_star="i"):
        returns (core_radius="d" | units.RSun)

    @remote_function(can_handle_array=True)
    def get_convective_envelope_mass(index_of_the_star="i"):
        returns (convective_envelope_mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_convective_envelope_radius(index_of_the_star="i"):
        returns (convective_envelope_radius="d" | units.RSun)

    @remote_function(can_handle_array=True)
    def get_CO_core_mass(index_of_the_star="i"):
        returns (CO_core_mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_main_sequence_lifetime(index_of_the_star="i"):
        returns (main_sequence_lifetime="d" | units.mega(units.julianyr))

    # getters and setters for tracks
    # metallicity_dir (string)
    # metallicity_dir_he (string)
    # z_accuracy_limit (float)
    # mass_accuracy_limit (float)

    @remote_function
    def get_metallicity_dir():
        returns (metallicity_dir="s")

    @remote_function
    def set_metallicity_dir(metallicity_dir="s"):
        returns ()

    @remote_function
    def get_metallicity_dir_he():
        returns (metallicity_dir_he="s")

    @remote_function
    def set_metallicity_dir_he(metallicity_dir_he="s"):
        returns ()

    @remote_function
    def get_z_accuracy_limit():
        returns (z_accuracy_limit="d")

    @remote_function
    def set_z_accuracy_limit(z_accuracy_limit="d"):
        returns ()

    @remote_function
    def get_mass_accuracy_limit():
        returns (mass_accuracy_limit="d")

    @remote_function
    def set_mass_accuracy_limit(mass_accuracy_limit="d"):
        returns ()

    # getters and setters for miscellaneous controls
    # verbose (bool)
    # construct_postagb_track (bool)

    @remote_function
    def get_verbose():
        returns (verbose="b")

    @remote_function
    def set_verbose(verbose="b"):
        returns ()

    @remote_function
    def get_construct_postagb_track():
        returns (construct_postagb_track="b")

    @remote_function
    def set_construct_postagb_track(construct_postagb_track="b"):
        returns ()

    # getters and setters for parameters
    # metallicity(real)
    # wd_mass_scheme (string, 256)
    # use_initial_final_mass_relation(bool)
    # bhns_mass_scheme (string, 256)
    # max_ns_mass (real)
    # allow_electron_capture (bool)

    @remote_function
    def get_metallicity():
        returns (metallicity="d")

    @remote_function
    def set_metallicity(metallicity="d"):
        returns ()

    @remote_function
    def get_wd_mass_scheme():
        returns (wd_mass_scheme="s")

    @remote_function
    def set_wd_mass_scheme(wd_mass_scheme="s"):
        returns ()

    @remote_function
    def get_use_initial_final_mass_relation():
        returns (use_initial_final_mass_relation="b")

    @remote_function
    def set_use_initial_final_mass_relation(use_initial_final_mass_relation="b"):
        returns ()

    @remote_function
    def get_bhns_mass_scheme():
        returns (bhns_mass_scheme="s")

    @remote_function
    def set_bhns_mass_scheme(bhns_mass_scheme="s"):
        returns ()

    @remote_function
    def get_max_ns_mass():
        returns (max_ns_mass="d")

    @remote_function
    def set_max_ns_mass(max_ns_mass="d"):
        returns ()

    @remote_function
    def get_allow_electron_capture():
        returns (allow_electron_capture="b")

    @remote_function
    def set_allow_electron_capture(allow_electron_capture="b"):
        returns ()

    @remote_function
    def get_time_step_pts_1():
        returns (fractional_time_step_1="d")

    @remote_function
    def set_time_step_pts_1(fractional_time_step_1="d"):
        returns ()

    @remote_function
    def get_time_step_pts_2():
        returns (fractional_time_step_2="d")

    @remote_function
    def set_time_step_pts_2(fractional_time_step_2="d"):
        returns ()

    @remote_function
    def get_time_step_pts_3():
        returns (fractional_time_step_3="d")

    @remote_function
    def set_time_step_pts_3(fractional_time_step_3="d"):
        returns ()


# high level interface class
class Metisse(se.StellarEvolution):
    """
    High level interface for METISSE
    """
    __interface__ = MetisseInterface

    def __init__(self, **options):
        # self.stopping_conditions = StoppingConditions(self)
        # self.stopping_conditions.supernova_detection = code.StoppingCondition('supernova_detection')
        se.StellarEvolution.__init__(self, MetisseInterface(**options), **options)
        self.model_time = 0.0 | units.mega(units.julianyr)

    # the definition of the state model of the code
    def define_state(self, handler):
        se.StellarEvolution.define_state(self, handler)
        handler.add_method('RUN', 'evolve_model')

    # the definition of any properties
    def define_properties(self, handler):
        # handler.add_property("name_of_the_getter", public_name="name_of_the_property")
        pass

    # the definition of the parameters
    def define_parameters(self, handler):
        # handler.add_method_parameter(
        #     "name_of_the_getter",
        #     "name_of_the_setter",
        #     "parameter_name",
        #     "description",
        #     default_value = <default value>
        # )

        # Track parameters
        handler.add_method_parameter(
            "get_metallicity_dir",
            "set_metallicity_dir",
            "metallicity_dir",
            "Location of the tracks",
            default_value=f"{os.path.dirname(__file__)}/data/Hydrogen",
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_metallicity_dir_he",
            "set_metallicity_dir_he",
            "metallicity_dir_he",
            "Location of the He tracks",
            default_value=f"{os.path.dirname(__file__)}/data/Helium",
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_z_accuracy_limit",
            "set_z_accuracy_limit",
            "z_accuracy_limit",
            "Metallicity accuracy limit",
            default_value=1.0e-2,
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_mass_accuracy_limit",
            "set_mass_accuracy_limit",
            "mass_accuracy_limit",
            "Mass accuracy limit",
            default_value=1.0e-4,
            must_set_before_get=True,
        )

        # handlers for parameters:
        # initial_metallicity
        # wd_mass_scheme
        # use_initial_final_mass_relation
        # bhns_mass_scheme
        # max_ns_mass
        # allow_electron_capture

        handler.add_method_parameter(
            "get_metallicity",
            "set_metallicity",
            "metallicity",
            "(initial) metallicity",
            default_value=-1.0,
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_wd_mass_scheme",
            "set_wd_mass_scheme",
            "wd_mass_scheme",
            (
                "White Dwarf (WD) luminosity calculation method:\n"
                "(1) \"Mestel\" - Shapiro S. L., Teukolsky S. A., 1983\n"
                "(2) \"Modified_mestel\" - Hurley J. R., Shara M. M., 2003"
            ),
            default_value="Modified_mestel",
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_use_initial_final_mass_relation",
            "set_use_initial_final_mass_relation",
            "use_initial_final_mass_relation",
            (
                "If True use the initial final mass relation for white dwarfs "
                "from Han, Z., Posialowski, P., Eggleton, P. P., 1995."
            ),
            default_value=False,
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_bhns_mass_scheme",
            "set_bhns_mass_scheme",
            "bhns_mass_scheme",
            (
                "Neutron Star/Black Hole (NS/BH) type and mass calculation method:\n"
                "(1) \"original_SSE\" - Hurley et al. 2000\n"
                "(2) \"Belczynski2002\" - Belczynski et al. 2002\n"
                "(3) \"Belczynski2008\" - Belczynski et al. 2008\n"
                "(4) \"Eldridge_Tout2004\" - Eldridge J. J., Tout C. A., 2004"
            ),
            default_value="Belczynski2008",
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_time_step_pts_1",
            "set_time_step_pts_1",
            "fractional_time_step_1",
            "Determine timestep for 95% of MS, and HeMS",
            default_value=0.05,
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_time_step_pts_2",
            "set_time_step_pts_2",
            "fractional_time_step_2",
            "Determine timestep for last 5% of MS, cHeBurn, HeHG, and HeGB",
            default_value=0.01,
            must_set_before_get=True,
        )

        handler.add_method_parameter(
            "get_time_step_pts_3",
            "set_time_step_pts_3",
            "fractional_time_step_3",
            "Determine timestep for HG, RGB, EAGB, and TPAGB",
            default_value=0.02,
            must_set_before_get=True,
        )

    def define_particle_sets(self, handler):
        handler.define_set("particles", "index_of_the_star")
        handler.set_new("particles", "new_particle")
        handler.set_delete("particles", "delete_star")

        handler.add_getter("particles", "get_mass", names=("mass",))
        handler.add_getter("particles", "get_radius", names=("radius",))
        handler.add_getter("particles", "get_age", names=("age",))
        handler.add_getter(
            "particles", "get_time_step", names=("time_step",)
        )
        handler.add_getter(
            "particles", "get_temperature", names=("temperature",)
        )
        handler.add_getter(
            "particles", "get_luminosity", names=("luminosity",)
        )
        handler.add_getter(
            "particles", "get_stellar_type", names=("stellar_type",)
        )
        # handler.add_getter("particles", "get_spin", names=("spin",))
        handler.add_getter("particles", "get_epoch", names=("epoch",))
        handler.add_getter(
            "particles",
            "get_main_sequence_lifetime",
            names=("main_sequence_lifetime",),
        )
        handler.add_getter(
            "particles", "get_core_mass", names=("core_mass",)
        )
        handler.add_getter(
            "particles", "get_CO_core_mass", names=("CO_core_mass",)
        )
        handler.add_getter(
            "particles", "get_core_radius", names=("core_radius",)
        )
        handler.add_getter(
            "particles",
            "get_convective_envelope_mass",
            names=("convective_envelope_mass",),
        )
        handler.add_getter(
            "particles",
            "get_convective_envelope_radius",
            names=("convective_envelope_radius",),
        )
        handler.add_getter(
            "particles", "get_initial_mass", names=("initial_mass",)
        )

        handler.add_method("particles", "evolve_one_step")
        handler.add_method("particles", "evolve_for")

    def evolve_model(self, end_time=None, keep_synchronous=True):
        if not keep_synchronous:
            for particle in self.particles:
                particle.evolve_one_step()
            return

        delta_time = (
            end_time-self.model_time
            if end_time
            else 0.99*min(self.particles.time_step)
        )
        for i, particle in enumerate(self.particles):
            delta_time = end_time - particle.age
            particle.evolve_for(delta_time)
        self.model_time += delta_time

    def define_methods(self, handler):
        se.StellarEvolution.define_methods(self, handler)

        # Metisse specific: uses Myr instead of yr
        handler.add_method(
            "evolve_for",
            (handler.INDEX, units.mega(units.yr)),
            (handler.ERROR_CODE,)
        )
