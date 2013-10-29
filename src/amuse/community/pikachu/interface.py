from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class PikachuInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, 
        StoppingConditionInterface, GravityFieldInterface, CodeWithDataDirectories):
    """
    Pikachu - a.k.a. P^3 Tree
    Hybrid N-body module, combining a tree (Barnes & Hut) to approximate long-range 
    forces, with direct summation of the forces from neighbour particles.
    """
    include_headers = ['worker_code.h', 'stopcond.h', 'interface.h']
    
    MODE_NORMAL = 'normal'
    MODE_LARGE_N = 'large_n'
    
    def __init__(self, mode=MODE_NORMAL, **options):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)
    
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'pikachu_worker'
        elif mode == self.MODE_LARGE_N:
            return 'pikachu_worker_large_n'
        else:
            print "Warning: unknown mode: '{0}' - using default ('{1}').".format(mode, self.MODE_NORMAL)
            return 'pikachu_worker'
    
    @option(type="string", sections=('data',))
    def default_kernel_directory(self):
        """
        The default directory containing the Sequoia kernels
        """
        return os.path.join(self.amuse_root_directory, 'src', 'amuse', 'community', 'pikachu')
    
    @legacy_function
    def get_kernel_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('kernel_directory', dtype='string', direction=function.OUT,
            description = "Name of the Sequoia kernel directory")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_kernel_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('kernel_directory', dtype='string', direction=function.IN,
            description = "Name of the Sequoia kernel directory")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eps2_fs_fs():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_fs', dtype='float64', direction=function.OUT,
            description = "The current value of the smooting parameter, squared, for star-star interactions.",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eps2_fs_fs():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_fs', dtype='float64', direction=function.IN,
            description = "The new value of the smooting parameter, squared, for star-star interactions.",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eps2_fs_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_bh', dtype='float64', direction=function.OUT,
            description = "The current value of the smooting parameter, squared, for star-blackhole interactions.",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eps2_fs_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_fs_bh', dtype='float64', direction=function.IN,
            description = "The new value of the smooting parameter, squared, for star-blackhole interactions.",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_eps2_bh_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_bh_bh', dtype='float64', direction=function.OUT,
            description = "The current value of the smooting parameter, squared, for blackhole-blackhole interactions.",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eps2_bh_bh():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared_bh_bh', dtype='float64', direction=function.IN,
            description = "The new value of the smooting parameter, squared, for blackhole-blackhole interactions.",
            unit = nbody_system.length * nbody_system.length)
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
            description = "The current value of the timestep parameter for black holes.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eta_smbh():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_supermassive_black_hole', dtype='float64', direction=function.IN,
            description = "The new value of the timestep parameter for black holes.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_theta_for_tree():
        function = LegacyFunctionSpecification()
        function.addParameter('opening_angle', dtype='float64', direction=function.OUT,
            description = "The opening angle, theta, for building the tree: between 0 and 1.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_theta_for_tree():
        function = LegacyFunctionSpecification()
        function.addParameter('opening_angle', dtype='float64', direction=function.IN,
            description = "The opening angle, theta, for building the tree: between 0 and 1.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_calculate_quadrupole_moments():
        function = LegacyFunctionSpecification()
        function.addParameter('calculate_quadrupole_moments', dtype='int32', direction=function.OUT,
            description = "Flag that specifies whether quadrupole moments are calculated for the tree")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_calculate_quadrupole_moments():
        function = LegacyFunctionSpecification()
        function.addParameter('calculate_quadrupole_moments', dtype='int32', direction=function.IN,
            description = "Flag that specifies whether quadrupole moments are calculated for the tree")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "The new value of the global timestep.",
            unit = nbody_system.time)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_search_factor():
        function = LegacyFunctionSpecification()
        function.addParameter('search_factor', dtype='float64', direction=function.OUT,
            description = "The search factor, if positive, determines rsearch = rcut_out + search_factor * velocity_dispersion * timestep")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_search_factor():
        function = LegacyFunctionSpecification()
        function.addParameter('opening_angle', dtype='float64', direction=function.IN,
            description = "The search factor, if positive, determines rsearch = rcut_out + search_factor * velocity_dispersion * timestep")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_vel_disp():
        function = LegacyFunctionSpecification()
        function.addParameter('vel_disp', dtype='float64', direction=function.OUT,
            description = "The velocity dispersion assumed when calculating rsearch",
            unit = nbody_system.speed)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_vel_disp():
        function = LegacyFunctionSpecification()
        function.addParameter('vel_disp', dtype='float64', direction=function.IN,
            description = "The velocity dispersion assumed when calculating rsearch",
            unit = nbody_system.speed)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_rcut_out_FS_FS():
        function = LegacyFunctionSpecification()
        function.addParameter('rcut_out_FS_FS', dtype='float64', direction=function.OUT,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_rcut_out_FS_FS():
        function = LegacyFunctionSpecification()
        function.addParameter('rcut_out_FS_FS', dtype='float64', direction=function.IN,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_rcut_out_FS_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rcut_out_FS_BH', dtype='float64', direction=function.OUT,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_rcut_out_FS_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rcut_out_FS_BH', dtype='float64', direction=function.IN,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_rcut_out_BH_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rcut_out_BH_BH', dtype='float64', direction=function.OUT,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_rcut_out_BH_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rcut_out_BH_BH', dtype='float64', direction=function.IN,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_rsearch_FS_FS():
        function = LegacyFunctionSpecification()
        function.addParameter('rsearch_FS_FS', dtype='float64', direction=function.OUT,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_rsearch_FS_FS():
        function = LegacyFunctionSpecification()
        function.addParameter('rsearch_FS_FS', dtype='float64', direction=function.IN,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_rsearch_FS_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rsearch_FS_BH', dtype='float64', direction=function.OUT,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_rsearch_FS_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rsearch_FS_BH', dtype='float64', direction=function.IN,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_rsearch_BH_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rsearch_BH_BH', dtype='float64', direction=function.OUT,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_rsearch_BH_BH():
        function = LegacyFunctionSpecification()
        function.addParameter('rsearch_BH_BH', dtype='float64', direction=function.IN,
            unit = nbody_system.length)
        function.result_type = 'int32'
        return function
    

class Pikachu(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = PikachuInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )
        
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        GravityFieldCode.define_state(self, object)


    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)
        self.stopping_conditions.define_parameters(object)
        object.add_method_parameter(
            "get_kernel_directory", 
            "set_kernel_directory",
            "kernel_directory", 
            "Name of the Sequoia kernel directory", 
            default_value = self.default_kernel_directory
        )
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
            default_value = 1.0e-8 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_eps2_fs_bh",
            "set_eps2_fs_bh", 
            "epsilon_squared_star_blackhole", 
            "smoothing parameter for gravity calculations - star-blackhole interactions only", 
            default_value = 1.0e-8 | nbody_system.length * nbody_system.length
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
            default_value = 0.005
        )
        object.add_alias_parameter(
            "timestep_parameter", 
            "timestep_parameter_stars", 
            "timestep parameter (eta) for field stars (alias for timestep_parameter_stars)"
        )
        object.add_method_parameter(
            "get_eta_fs",
            "set_eta_fs", 
            "timestep_parameter_stars", 
            "timestep parameter (eta) for field stars", 
            default_value = 0.025
        )
        object.add_method_parameter(
            "get_eta_smbh",
            "set_eta_smbh", 
            "timestep_parameter_black_holes", 
            "timestep parameter (eta) for black holes", 
            default_value = 0.025
        )
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "global timestep for iteration", 
            default_value = 1.0 / 2048.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_theta_for_tree",
            "set_theta_for_tree",
            "opening_angle", 
            "opening angle, theta, for building the tree: between 0 and 1", 
            default_value = 0.4
        )
        object.add_method_parameter(
            "get_search_factor",
            "set_search_factor",
            "search_factor", 
            "search factor, if positive, determines rsearch = rcut_out + search_factor * velocity_dispersion * timestep", 
            default_value = 3.0
        )
        object.add_method_parameter(
            "get_vel_disp",
            "set_vel_disp",
            "velocity_dispersion", 
            "velocity dispersion assumed when calculating rsearch", 
            default_value = 0.707106781 | nbody_system.speed
        )
        object.add_method_parameter(
            "get_rcut_out_FS_FS",
            "set_rcut_out_FS_FS",
            "rcut_out_star_star", 
            "cut-off radius beyond which direct force calculations smoothly transition into tree approximations", 
            default_value = 2.0e-3 | nbody_system.length
        )
        object.add_method_parameter(
            "get_rcut_out_FS_BH",
            "set_rcut_out_FS_BH",
            "rcut_out_star_blackhole", 
            "cut-off radius beyond which direct force calculations smoothly transition into tree approximations", 
            default_value = 2.0e-2 | nbody_system.length
        )
        object.add_method_parameter(
            "get_rcut_out_BH_BH",
            "set_rcut_out_BH_BH",
            "rcut_out_blackhole_blackhole", 
            "cut-off radius beyond which direct force calculations smoothly transition into tree approximations", 
            default_value = 1.0e5 | nbody_system.length
        )
        object.add_method_parameter(
            "get_rsearch_FS_FS",
            "set_rsearch_FS_FS",
            "rsearch_star_star", 
            "maximum radius for neighbour search, must be larger than rcut_out to "
            "provide a buffer for particles moving into rcut_out during a time step, "
            "only effective if search_factor <= 0", 
            default_value = 0.0 | nbody_system.length
        )
        object.add_method_parameter(
            "get_rsearch_FS_BH",
            "set_rsearch_FS_BH",
            "rsearch_star_blackhole", 
            "maximum radius for neighbour search, must be larger than rcut_out to "
            "provide a buffer for particles moving into rcut_out during a time step, "
            "only effective if search_factor <= 0", 
            default_value = 0.0 | nbody_system.length
        )
        object.add_method_parameter(
            "get_rsearch_BH_BH",
            "set_rsearch_BH_BH",
            "rsearch_blackhole_blackhole", 
            "maximum radius for neighbour search, must be larger than rcut_out to "
            "provide a buffer for particles moving into rcut_out during a time step, "
            "only effective if search_factor <= 0", 
            default_value = 0.0 | nbody_system.length
        )
        object.add_boolean_parameter(
            "get_calculate_quadrupole_moments",
            "set_calculate_quadrupole_moments",
            "calculate_quadrupole_moments",
            "Flag that specifies whether quadrupole moments are calculated for the tree",
            False
        )
    
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        self.stopping_conditions.define_methods(object)
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object)
