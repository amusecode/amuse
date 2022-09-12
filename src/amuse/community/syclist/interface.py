from amuse.units import units
from amuse.community import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
    LiteratureReferencesMixIn, CodeWithDataDirectories,
    InCodeComponentImplementation,
    NO_UNIT,
)
from amuse.community.interface.common import (
    CommonCodeInterface,
    CommonCode,
)
# from amuse.community.interface.se import (
#     StellarEvolutionInterface,
#     StellarEvolution
# )
from amuse.datamodel import Particles, ParticlesSubset


class SyclistInterface(
    CodeInterface,
    CommonCodeInterface,
    CodeWithDataDirectories,
    LiteratureReferencesMixIn,
):
    """
    You have used the SYCLIST population synthesis code.
        .. [#] [2014A&A...566A..21G] Georgy et al., 2014

    Isochrones are calculated on the basis of the following grids of stellar
    models:

    B-type stars grids:
    - M = 1.7 - 15 M⊙
    - Ω/Ωcrit = 0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95
    - Z = 0.014 (=Z⊙), 0.006, 0.002

        .. [#] Georgy et al. (2013), A&A 553 A24

    large grids:
    - M = 0.8 - 120 M⊙
    - V/Vcrit = 0.0, 0.4 (i.e. Ω/Ωcrit = 0.0, 0.568)
    - Z = 0.020, 0.014 (=Z⊙), 0.006, 0.002, 0.0004

    published in:
        .. [#] Yusof et al. (2022), MNRAS 511 2814 (Z=0.020)
        .. [#] Ekström et al. (2012), A&A 537 A146 (Z=0.014=Z⊙)
        .. [#] Eggenberger et al. (2021) A&A 652 A137 (Z=0.006)
        .. [#] Georgy et al. (2013) A&A 558 A103 (Z=0.002)
        .. [#] Groh et al. (2019) A&A 627 A24 (Z=0.0004)
        .. [#] Yusof et al. (2013) MNRAS 433 1114 (M≥150 M⊙)
    NB: the most massive models (150-500 M⊙) have been evolved without the
        insurance of angular momentum conservation.

    low-mass stars dense grids:
    - M = 0.5 - 3.5 M⊙
    - Z = 0.006, 0.010, 0.014 (=Z⊙), 0.020, 0.030, 0.040
    - no rotation

    published in:
        .. [#] Mowlavi et al. (2012), A&A 541 A41

    very low-mass stars grids:
    - M = 0.2 - 1.5 M⊙
    - [Fe/H] = -1.0, -0.5, -0.3, -0.15, 0, 0.15, 0.3
    - no rotation, slow-, medium-, fast-rotation
     
    published in:
        .. [#] Amard et al. (2019), A&A 631 A77
    """

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="syclist_worker",
            **keyword_arguments
        )
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def set_models_path():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'modelpath', dtype='string', direction=function.IN,
            description="Path from which the models are read",
        )
        function.result_type = 'int32'
        function.result_doc = ''
        return function

    @legacy_function
    def get_models_path():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'modelpath', dtype='string', direction=function.OUT,
            description="Path from which the models are read",
        )
        function.result_type = 'int32'
        function.result_doc = ''
        return function

    @legacy_function
    def evolve_star():
        function = LegacyFunctionSpecification()
        # function.name = 'evolve_star'
        # function.can_handle_array = True
        function.must_handle_array = True
        stellar_properties = {
            'age': [0., units.julianyr],
            'initial_mass': [0., units.MSun],
            'metallicity': [0.014, NO_UNIT],
            'omega': [0., NO_UNIT],
            'mass': [0., units.MSun],
            'luminosity': [0., units.LSun],
            'temperature': [0., units.K],
            'abundance_1H': [0., NO_UNIT],
            'abundance_4He': [0., NO_UNIT],
            'abundance_12C': [0., NO_UNIT],
            'abundance_13C': [0., NO_UNIT],
            'abundance_14N': [0., NO_UNIT],
            'abundance_16O': [0., NO_UNIT],
            'abundance_17O': [0., NO_UNIT],
            'abundance_18O': [0., NO_UNIT],
            'abundance_20Ne': [0., NO_UNIT],
            'abundance_22Ne': [0., NO_UNIT],
            'abundance_26Al': [0., NO_UNIT],
            'MccMt': [0., NO_UNIT],
            'temperature_wr': [0., units.K],
            'Md': [0., NO_UNIT],
            'rhoc': [0., NO_UNIT],
            'abundance_core_1H': [0., NO_UNIT],
            'abundance_core_4He': [0., NO_UNIT],
            'abundance_core_12C': [0., NO_UNIT],
            'abundance_core_13C': [0., NO_UNIT],
            'abundance_core_14N': [0., NO_UNIT],
            'abundance_core_16O': [0., NO_UNIT],
            'abundance_core_17O': [0., NO_UNIT],
            'abundance_core_20Ne': [0., NO_UNIT],
            'abundance_core_22Ne': [0., NO_UNIT],
            'abundance_core_26Al': [0., NO_UNIT],
            'omega_surface': [0., units.s**-1],
            'omega_core': [0., units.s**-1],
            'RpReq': [0., NO_UNIT],
            'MdMd0': [0., NO_UNIT],
            'v_crit1': [0., units.kms],
            'v_crit2': [0., units.kms],
            'v_equator': [0., units.kms],
            'OmOmcr': [0., NO_UNIT],
            'Gamma_Ed': [0., NO_UNIT],
            'MdotMech': [0., units.MSun/units.yr],
            'Ltot': [0., 10**53*units.g*units.cm**2/units.s],
        }
        for stellar_property in stellar_properties.items():
            function.addParameter(
                stellar_property[0],
                dtype='float64',
                direction=function.INOUT,
                unit=stellar_property[1][1],
            )
        function.addParameter(
            'n',
            dtype='int64',
            direction=function.LENGTH,
        )

        return function

    def recommit_particles(self):
        return 0

    def commit_particles(self):
        return 0


class SyclistParticles(
    Particles
):
    def __init__(self, code_interface, storage=None):
        Particles.__init__(self, storage=storage)
        self._private.code_interface = code_interface
        self.add_function_attribute(
            "evolve_for",
            self.particleset_evolve_for,
            self.evolve_for,
        )
    def evolve_one_step(self, particles, subset):
        self._private.code_interface._evolve_particles(
            subset.as_set(),
            subset.age + subset.time_step,
        )

    def add_particles_to_store(self, keys, attributes=[], values=[]):
        if len(keys) == 0:
            return

        all_attributes = []
        all_attributes.extend(attributes)
        all_values = []
        all_values.extend(values)

        mapping_from_attribute_to_default_value = {
            'age': 0. | units.julianyr,
            'metallicity': 0.014,
            'omega': 0. | units.none,
            'luminosity': 0. | units.LSun,
            'temperature': 0. | units.K,
            'abundance_1H': 0. | units.none,
            'abundance_4He': 0. | units.none,
            'abundance_12C': 0. | units.none,
            'abundance_13C': 0. | units.none,
            'abundance_14N': 0. | units.none,
            'abundance_16O': 0. | units.none,
            'abundance_17O': 0. | units.none,
            'abundance_18O': 0. | units.none,
            'abundance_20Ne': 0. | units.none,
            'abundance_22Ne': 0. | units.none,
            'abundance_26Al': 0. | units.none,
            'MccMt': 0. | units.none,
            'temperature_wr': 0. | units.K,
            'Md': 0. | units.none,
            'rhoc': 0. | units.none,
            'abundance_core_1H': 0. | units.none,
            'abundance_core_4He': 0. | units.none,
            'abundance_core_12C': 0. | units.none,
            'abundance_core_13C': 0. | units.none,
            'abundance_core_14N': 0. | units.none,
            'abundance_core_16O': 0. | units.none,
            'abundance_core_17O': 0. | units.none,
            'abundance_core_20Ne': 0. | units.none,
            'abundance_core_22Ne': 0. | units.none,
            'abundance_core_26Al': 0. | units.none,
            'omega_surface': 0. | units.s**-1,
            'omega_core': 0. | units.s**-1,
            'RpReq': 0. | units.none,
            'MdMd0': 0. | units.none,
            'v_crit1': 0. | units.kms,
            'v_crit2': 0. | units.kms,
            'v_equator': 0. | units.kms,
            'OmOmcr': 0. | units.none,
            'Gamma_Ed': 0. | units.none,
            'MdotMech': 0. | units.MSun/units.yr,
            'Ltot': 0. | 10**53*units.g*units.cm**2/units.s,
        }
        given_attributes = set(attributes)

        if "initial_mass" not in given_attributes:
            index_of_mass_attribute = attributes.index("mass")
            all_attributes.append("initial_mass")
            all_values.append(values[index_of_mass_attribute].in_(units.MSun))

        for attribute, default_value in (
                mapping_from_attribute_to_default_value.items()
        ):
            if attribute not in given_attributes:
                all_attributes.append(attribute)
                all_values.append(
                    default_value.as_vector_with_length(len(keys))
                )

        super(SyclistParticles, self).add_particles_to_store(
            keys, all_attributes, all_values
        )

        added_particles = ParticlesSubset(self, keys)
        # self._private.code_interface.evolve_star((
        #     added_particles, 0 | units.yr
        # )

    def particleset_evolve_one_step(self, particles):
        self._private.code_interface._evolve_particles(
            particles,
            particles.age + particles.time_step,
        )

    def evolve_for(self, particles, subset, delta_time):
        self._private.code_interface._evolve_particles(
            subset.as_set(),
            subset.age + delta_time,
        )
    
    def particleset_evolve_for(self, particles, delta_time):
        self._private.code_interface._evolve_particles(
            particles,
            particles.age + delta_time,
        )


class Syclist(
    CommonCode,
    # StellarEvolution,
    InCodeComponentImplementation
):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(
            self, SyclistInterface(**options), **options)
        self.model_time = 0.0 | units.yr

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_models_path",
            "set_models_path",
            "models_path",
            "Path where the models are stored",
            "./src/SYCLIST/inputs/"
        )
        # set defaults for metallicity, rotation etc

    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition(
            'INITIALIZED', 'RUN', 'commit_parameters'
        )
        handler.add_method('RUN', 'evolve_star')

        # Only allow changing of models path before they are read
        handler.add_method('INITIALIZED', 'set_models_path')

    def define_particle_sets(self, handler):
        handler.define_inmemory_set(
            'particles',
            SyclistParticles,
        )

    def evolve_model(
        self,
        end_time=None,
    ):
        if end_time is None:
            end_time = self.model_time
        self.particles.age = self.particles.age + end_time - self.model_time

        result = self.evolve_star(
            self.particles.age,  # 0
            self.particles.initial_mass,  # 1
            self.particles.metallicity,  # 2
            self.particles.omega,  # 3
            self.particles.mass,  # 4
            self.particles.luminosity,  # 5
            self.particles.temperature,  # 6
            self.particles.abundance_1H,  # 7
            self.particles.abundance_4He,  # 8
            self.particles.abundance_12C,  # 9
            self.particles.abundance_13C,  # 10
            self.particles.abundance_14N,  # 11
            self.particles.abundance_16O,  # 12
            self.particles.abundance_17O,  # 13
            self.particles.abundance_18O,  # 14
            self.particles.abundance_20Ne,  # 15
            self.particles.abundance_22Ne,  # 16
            self.particles.abundance_26Al,  # 17
            self.particles.MccMt,  # 18
            self.particles.temperature_wr,  # 19
            self.particles.Md,  # 20
            self.particles.rhoc,  # 21
            self.particles.abundance_core_1H,  # 22
            self.particles.abundance_core_4He,  # 23
            self.particles.abundance_core_12C,  # 24
            self.particles.abundance_core_13C,  # 25
            self.particles.abundance_core_14N,  # 26
            self.particles.abundance_core_16O,  # 27
            self.particles.abundance_core_17O,  # 28
            self.particles.abundance_core_20Ne,  # 29
            self.particles.abundance_core_22Ne,  # 30
            self.particles.abundance_core_26Al,  # 31
            self.particles.omega_surface,  # 32
            self.particles.omega_core,  # 33
            self.particles.RpReq,  # 34
            self.particles.MdMd0,  # 35
            self.particles.v_crit1,  # 36
            self.particles.v_crit2,  # 37
            self.particles.v_equator,  # 38
            self.particles.OmOmcr,  # 39
            self.particles.Gamma_Ed,  # 40
            self.particles.MdotMech,  # 41
            self.particles.Ltot,  # 42
            # ...
        )
        self.particles.age = result[0]
        self.particles.initial_mass = result[1]
        self.particles.metallicity = result[2]
        self.particles.omega = result[3]
        self.particles.mass = result[4]
        self.particles.luminosity = result[5]
        self.particles.temperature = result[6]
        self.particles.abundance_1H = result[7]
        self.particles.abundance_4He = result[8]
        self.particles.abundance_12C = result[9]
        self.particles.abundance_13C = result[10]
        self.particles.abundance_14N = result[11]
        self.particles.abundance_16O = result[12]
        self.particles.abundance_17O = result[13]
        self.particles.abundance_18O = result[14]
        self.particles.abundance_20Ne = result[15]
        self.particles.abundance_22Ne = result[16]
        self.particles.abundance_26Al = result[17]
        self.particles.MccMt = result[18]
        self.particles.temperature_wr = result[19]
        self.particles.Md = result[20]
        self.particles.rhoc = result[21]
        self.particles.abundance_core_1H = result[22]
        self.particles.abundance_core_4He = result[23]
        self.particles.abundance_core_12C = result[24]
        self.particles.abundance_core_13C = result[25]
        self.particles.abundance_core_14N = result[26]
        self.particles.abundance_core_16O = result[27]
        self.particles.abundance_core_17O = result[28]
        self.particles.abundance_core_20Ne = result[29]
        self.particles.abundance_core_22Ne = result[30]
        self.particles.abundance_core_26Al = result[31]
        self.particles.omega_surface = result[32]
        self.particles.omega_core = result[33]
        self.particles.RpReq = result[34]
        self.particles.MdMd0 = result[35]
        self.particles.v_crit1 = result[36]
        self.particles.v_crit2 = result[37]
        self.particles.v_equator = result[38]
        self.particles.OmOmcr = result[39]
        self.particles.Gamma_Ed = result[40]
        self.particles.MdotMech = result[41]
        self.particles.Ltot = result[42]
        self.model_time = end_time
