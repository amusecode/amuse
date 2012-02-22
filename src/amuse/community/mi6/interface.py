from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

class MI6Interface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, StoppingConditionInterface):
    """
    MI6 - Masaki's Integrator 6th order
    N-body module with mixed 4th and 6th order Hermite integration scheme, aimed 
    at simulating the Galactic center. It has different treatments for 
    supermassive black holes (SMBHs), intermediate mass black holes (IMBHs), and 
    normal "field stars". The post-newtonian effects are included up to 2.5.
    
    NOTE: there is always a SMBH fixed at the center, with a mass of unity (nbody 
    units). Masses of other particles are assumed to be much less than unity, 
    such that the SMBH remains at the origin. Therefore a test particle will 
    orbit around the SMBH with a velocity proportional to sqrt(GM/r), with GM=1, 
    since the mass of the test particle itself is neglected.
    Validity of the results is questionable when using particles with masses of 
    order unity or greater, or when the center-of-mass (-velocity) does not 
    coincide with the origin.
    """
    include_headers = ['worker_code.h', 'stopcond.h', 'interface.h']
    
    MODE_GPU = 'gpu'
    MODE_CPU = 'cpu'
    
    def __init__(self, mode = MODE_CPU, **options):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)
    
    def name_of_the_worker(self, mode):
        if mode == self.MODE_CPU:
            return 'mi6_worker'
        elif mode == self.MODE_GPU:
            return 'mi6_worker_gpu'
        else:
            print "Warning: unknown mode: '{0}' - using default ('{1}').".format(mode, self.MODE_CPU)
            return 'mi6_worker'

    @legacy_function
    def get_gravity_at_point():
        """
        Determine the gravitational force at a given point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for par in ['eps', 'x', 'y', 'z']:
            function.addParameter(par, dtype='float64', direction=function.IN)
        for par in ['ax', 'ay', 'az']:
            function.addParameter(par, dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_potential_at_point():
        """
        Determine the potential at a given point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for par in ['eps', 'x', 'y', 'z']:
            function.addParameter(par, dtype='float64', direction=function.IN)
        function.addParameter('phi', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eps2_fs_fs():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_fs', dtype='float64', direction=function.OUT,
            description = "The current value of the smooting parameter, squared, for star-star interactions.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eps2_fs_fs():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_fs', dtype='float64', direction=function.IN,
            description = "The new value of the smooting parameter, squared, for star-star interactions.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eps2_fs_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_bh', dtype='float64', direction=function.OUT,
            description = "The current value of the smooting parameter, squared, for star-blackhole interactions.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eps2_fs_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_bh', dtype='float64', direction=function.IN,
            description = "The new value of the smooting parameter, squared, for star-blackhole interactions.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eps2_bh_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_bh_bh', dtype='float64', direction=function.OUT,
            description = "The current value of the smooting parameter, squared, for blackhole-blackhole interactions.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eps2_bh_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_bh_bh', dtype='float64', direction=function.IN,
            description = "The new value of the smooting parameter, squared, for blackhole-blackhole interactions.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eta_s():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_start', dtype='float64', direction=function.OUT,
            description = "The current value of the initial timestep parameter.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eta_s():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_start', dtype='float64', direction=function.IN,
            description = "The new value of the initial timestep parameter.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eta_fs():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_field_star', dtype='float64', direction=function.OUT,
            description = "The current value of the timestep parameter for field stars.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eta_fs():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_field_star', dtype='float64', direction=function.IN,
            description = "The new value of the timestep parameter for field stars.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eta_smbh():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_supermassive_black_hole', dtype='float64', direction=function.OUT,
            description = "The current value of the timestep parameter for supermassive black holes.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eta_smbh():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_supermassive_black_hole', dtype='float64', direction=function.IN,
            description = "The new value of the timestep parameter for supermassive black holes.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eta_imbh():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_intermediate_mass_black_hole', dtype='float64', direction=function.OUT,
            description = "The current value of the timestep parameter for intermediate mass black holes.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eta_imbh():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_intermediate_mass_black_hole', dtype='float64', direction=function.IN,
            description = "The new value of the timestep parameter for intermediate mass black holes.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_max_relative_energy_error():
        function = LegacyFunctionSpecification()
        function.addParameter('max_relative_energy_error', dtype='float64', direction=function.OUT,
            description = "The current value of the maximum relative energy error per full step.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_max_relative_energy_error():
        function = LegacyFunctionSpecification()
        function.addParameter('max_relative_energy_error', dtype='float64', direction=function.IN,
            description = "The new value of the maximum relative energy error per full step.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_maximum_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('maximum_timestep', dtype='float64', direction=function.OUT,
            description = "The current value of the maximum timestep a particle may take.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_maximum_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('maximum_timestep', dtype='float64', direction=function.IN,
            description = "The new value of the maximum timestep a particle may take.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_smbh_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('smbh_mass', dtype='float64', direction=function.OUT,
            description = "The current value of the mass of the supermassive black hole at the center.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_smbh_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('smbh_mass', dtype='float64', direction=function.IN,
            description = "The new value of the mass of the supermassive black hole at the center.")
        function.result_type = 'int32'
        return function
    



class MI6(GravitationalDynamics):

##    __doc__ = MasakiDoc()

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = MI6Interface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)
        self.stopping_conditions.define_parameters(object)
        object.add_alias_parameter(
            "epsilon_squared", 
            "epsilon_squared_star_star", 
            "smoothing parameter for gravity calculations - star-star interactions only (alias for epsilon_squared_star_star)"
        )
        object.add_method_parameter(
            "get_eps2_fs_fs",
            "set_eps2_fs_fs", 
            "epsilon_squared_star_star", 
            "smoothing parameter for gravity calculations - star-star interactions only", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_eps2_fs_bh",
            "set_eps2_fs_bh", 
            "epsilon_squared_star_blackhole", 
            "smoothing parameter for gravity calculations - star-blackhole interactions only", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_eps2_bh_bh",
            "set_eps2_bh_bh", 
            "epsilon_squared_blackhole_blackhole", 
            "smoothing parameter for gravity calculations - blackhole-blackhole interactions only", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_eta_s",
            "set_eta_s", 
            "initial_timestep_parameter", 
            "initial timestep parameter (eta)", 
            default_value = 1.0e-4 | units.none
        )
        object.add_method_parameter(
            "get_eta_fs",
            "set_eta_fs", 
            "timestep_parameter", 
            "timestep parameter (eta) for field stars (alias for timestep_parameter_stars)", 
            default_value = 0.1 | units.none
        )
        object.add_method_parameter(
            "get_eta_fs",
            "set_eta_fs", 
            "timestep_parameter_stars", 
            "timestep parameter (eta) for field stars", 
            default_value = 0.1 | units.none
        )
        object.add_method_parameter(
            "get_eta_smbh",
            "set_eta_smbh", 
            "timestep_parameter_supermassive_black_holes", 
            "timestep parameter (eta) for supermassive black holes", 
            default_value = 0.4 | units.none
        )
        object.add_method_parameter(
            "get_eta_imbh",
            "set_eta_imbh", 
            "timestep_parameter_intermediate_mass_black_holes", 
            "timestep parameter (eta) for intermediate mass black holes", 
            default_value = 0.4 | units.none
        )
        object.add_method_parameter(
            "get_drink",
            None, 
            "drink", 
            "Order a drink at MI6", 
            default_value = ""
        )
        object.add_method_parameter(
            "get_max_relative_energy_error",
            "set_max_relative_energy_error", 
            "max_relative_energy_error", 
            "the maximum relative energy error per full step", 
            default_value = 5e-5 | units.none # or nbody_system.time**-1 ??? why /dt_max ???
        )
        object.add_method_parameter(
            "get_maximum_timestep",
            "set_maximum_timestep", 
            "maximum_timestep", 
            "the maximum timestep a particle may take", 
            default_value = 1.0/1024.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_smbh_mass",
            "set_smbh_mass", 
            "smbh_mass", 
            "the mass of the supermassive black hole at the center", 
            default_value = 1.0 | nbody_system.mass
        )
    
    def get_drink(self):
        return "Vodka martini. Shaken, not stirred."
    
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        self.stopping_conditions.define_methods(object)
        
        object.add_method("get_eps2", (), (nbody_system.length**2, object.ERROR_CODE))
        object.add_method("set_eps2", (nbody_system.length**2,), (object.ERROR_CODE,))
        object.add_method("get_eps2_fs_fs", (), (nbody_system.length**2, object.ERROR_CODE))
        object.add_method("set_eps2_fs_fs", (nbody_system.length**2,), (object.ERROR_CODE,))
        object.add_method("get_eps2_fs_bh", (), (nbody_system.length**2, object.ERROR_CODE))
        object.add_method("set_eps2_fs_bh", (nbody_system.length**2,), (object.ERROR_CODE,))
        object.add_method("get_eps2_bh_bh", (), (nbody_system.length**2, object.ERROR_CODE))
        object.add_method("set_eps2_bh_bh", (nbody_system.length**2,), (object.ERROR_CODE,))
        
        object.add_method("get_eta_s", (), (units.none, object.ERROR_CODE))
        object.add_method("set_eta_s", (units.none,), (object.ERROR_CODE,))
        object.add_method("get_eta_fs", (), (units.none, object.ERROR_CODE))
        object.add_method("set_eta_fs", (units.none,), (object.ERROR_CODE,))
        object.add_method("get_eta_smbh", (), (units.none, object.ERROR_CODE))
        object.add_method("set_eta_smbh", (units.none,), (object.ERROR_CODE,))
        object.add_method("get_eta_imbh", (), (units.none, object.ERROR_CODE))
        object.add_method("set_eta_imbh", (units.none,), (object.ERROR_CODE,))
        
        object.add_method("get_max_relative_energy_error", (), (units.none, object.ERROR_CODE))
        object.add_method("set_max_relative_energy_error", (units.none,), (object.ERROR_CODE,))
        object.add_method("get_maximum_timestep", (), (nbody_system.time, object.ERROR_CODE))
        object.add_method("set_maximum_timestep", (nbody_system.time,), (object.ERROR_CODE,))
        
        object.add_method("get_smbh_mass", (), (nbody_system.mass, object.ERROR_CODE))
        object.add_method("set_smbh_mass", (nbody_system.mass,), (object.ERROR_CODE,))
        
#~        object.add_method(
#~            'get_potential_at_point',
#~            (nbody_system.length, nbody_system.length,
#~             nbody_system.length, nbody_system.length),
#~            (nbody_system.potential, object.ERROR_CODE)
#~        )
#~
#~        object.add_method(
#~            'get_gravity_at_point',
#~            (nbody_system.length, nbody_system.length,
#~             nbody_system.length, nbody_system.length),
#~            (nbody_system.acceleration, nbody_system.acceleration,
#~             nbody_system.acceleration, object.ERROR_CODE)
#~        )
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object, 'particles')
