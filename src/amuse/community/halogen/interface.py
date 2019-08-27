import os.path
from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.options import option

class HalogenInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn,
        CodeWithDataDirectories):
    """
    This is a stripped-down version of Halogen, developed during the Modest 7a workshop 
    in Split, Croatia, in August 2007. This version can be used for generating 
    single-mass initial conditions.
    
    Halogen allows to generate spherical structures from the alpha-beta-gamma-model 
    family with an isotropic velocity tensor. Particles are sampled self-consistently 
    from the distribution function of the models.
    
    Relevant references:
        .. [#] Zemp M., Moore B., Stadel J., Carollo C.M. & Madau P. 2008, MNRAS, 386, 1543 
    """
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="halogen_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    @legacy_function
    def generate_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    new_particle = None
    
    def delete_particle(self, index_of_the_particle):
        return 0
    
    def invoke_state_change2(self):
        pass
    
    def invoke_state_change_updated(self):
        pass
    
    @legacy_function
    def get_number_of_particles_updated():
        """
        Return the number of particles added during the last generate_particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current x component of the position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current y component of the position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current z component of the position vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current x component of the velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current y component of the velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current z component of the velocity vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
    
    
    @legacy_function
    def get_model_alpha():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_model_alpha():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_model_beta():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_model_beta():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_model_gamma():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_model_gamma():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_total_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_total_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_scale_radius():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_scale_radius():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_cutoff_radius():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_cutoff_radius():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_target_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_target_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_black_hole_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_black_hole_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_do_exact_virial_radius_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_do_exact_virial_radius_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_outputgridr_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_outputgridr_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_outputgriddf_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_outputgriddf_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_write_output_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_write_output_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_output_path():
        function = LegacyFunctionSpecification()
        function.addParameter('output_directory', dtype='string', direction=function.OUT,
            description = "The path to the output directory.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_output_path():
        function = LegacyFunctionSpecification()
        function.addParameter('output_directory', dtype='string', direction=function.IN,
            description = "The path to the output directory.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_output_basename():
        function = LegacyFunctionSpecification()
        function.addParameter('output_basename', dtype='string', direction=function.OUT,
            description = "The basename of output files.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_output_basename():
        function = LegacyFunctionSpecification()
        function.addParameter('output_basename', dtype='string', direction=function.IN,
            description = "The basename of output files.")
        function.result_type = 'int32'
        return function
    



class Halogen(CommonCode):
    
    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        InCodeComponentImplementation.__init__(self, HalogenInterface(**options), **options)
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        self.parameters.set_defaults()
        self.parameters.output_directory = self.get_output_directory()
    
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_model_alpha",
            "set_model_alpha",
            "alpha",
            "alpha parameter in density profile (see amuse/community/halogen/src/doc for details)",
            default_value = -1.0
        )
        
        handler.add_method_parameter(
            "get_model_beta",
            "set_model_beta",
            "beta",
            "beta parameter in density profile (see amuse/community/halogen/src/doc for details)",
            default_value = -1.0
        )
        
        handler.add_method_parameter(
            "get_model_gamma",
            "set_model_gamma",
            "gamma",
            "gamma parameter in density profile (see amuse/community/halogen/src/doc for details)",
            default_value = -1.0
        )
        
        handler.add_method_parameter(
            "get_total_mass",
            "set_total_mass",
            "total_mass",
            "the total mass of the model",
            default_value = 1.0 | nbody_system.mass
        )
        
        handler.add_method_parameter(
            "get_scale_radius",
            "set_scale_radius",
            "scale_radius",
            "the scale radius of the density profile (see amuse/community/halogen/src/doc for details)",
            default_value = 1.0 | nbody_system.length
        )
        
        handler.add_method_parameter(
            "get_cutoff_radius",
            "set_cutoff_radius",
            "cutoff_radius",
            "the cutoff radius of the density profile (see amuse/community/halogen/src/doc for details)",
            default_value = -1.0 | nbody_system.length
        )
        
        handler.add_method_parameter(
            "get_target_number_of_particles",
            "set_target_number_of_particles",
            "number_of_particles",
            "the number of particles to be generated in the model",
            default_value = -1
        )
        
        handler.add_method_parameter(
            "get_black_hole_mass",
            "set_black_hole_mass",
            "black_hole_mass",
            "the mass of the central black hole",
            default_value = 0.0 | nbody_system.mass
        )
        
        handler.add_method_parameter(
            "get_random_seed",
            "set_random_seed",
            "random_seed",
            "the initial seed to be used by the random number generator",
            default_value = 42
        )
        
        handler.add_boolean_parameter(
            "get_do_exact_virial_radius_flag",
            "set_do_exact_virial_radius_flag",
            "do_exact_virial_radius_flag",
            "Flag specifying whether to calculate the virial radius exactly via N^2 sum - Warning: time consuming for large N!",
            False
        )
        
        handler.add_boolean_parameter(
            "get_outputgridr_flag",
            "set_outputgridr_flag",
            "outputgridr_flag",
            "Flag specifying whether to write a file outputting grid in r",
            False
        )
        
        handler.add_boolean_parameter(
            "get_outputgriddf_flag",
            "set_outputgriddf_flag",
            "outputgriddf_flag",
            "Flag specifying whether to write a file outputting grid for distribution function",
            False
        )
        
        handler.add_boolean_parameter(
            "get_write_output_flag",
            "set_write_output_flag",
            "write_output_flag",
            "Flag specifying whether to write the model and a log to file",
            False
        )
        
        handler.add_method_parameter(
            "get_output_path", 
            "set_output_path",
            "output_directory", 
            "The path to the output directory", 
            default_value = "./"
        )
        
        handler.add_method_parameter(
            "get_output_basename", 
            "set_output_basename",
            "output_basename", 
            "The basename of output files", 
            default_value = "halogen"
        )
    
    def define_errorcodes(self, handler):
        handler.add_errorcode(-1, 'Unspecified, other error.')
        handler.add_errorcode(-2, 'Missing or bad parameter for halo (see amuse/community/halogen/src/doc for details on required parameters).')
    
    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method("generate_particles", (), (handler.ERROR_CODE,))
        handler.add_method("get_number_of_particles_updated", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler.add_method("get_mass", (handler.INDEX,), 
            (nbody_system.mass, handler.ERROR_CODE)
        )
        handler.add_method("get_position", (handler.INDEX,), 
            (nbody_system.length, nbody_system.length, nbody_system.length, handler.ERROR_CODE)
        )
        handler.add_method("get_velocity", (handler.INDEX,), 
            (nbody_system.speed, nbody_system.speed, nbody_system.speed, handler.ERROR_CODE)
        )
        
        handler.add_method("get_model_alpha", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_model_alpha", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_model_beta", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_model_beta", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_model_gamma", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_model_gamma", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_total_mass", (), (nbody_system.mass, handler.ERROR_CODE,))
        handler.add_method("set_total_mass", (nbody_system.mass, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_scale_radius", (), (nbody_system.length, handler.ERROR_CODE,))
        handler.add_method("set_scale_radius", (nbody_system.length, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_cutoff_radius", (), (nbody_system.length, handler.ERROR_CODE,))
        handler.add_method("set_cutoff_radius", (nbody_system.length, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_target_number_of_particles", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_target_number_of_particles", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_black_hole_mass", (), (nbody_system.mass, handler.ERROR_CODE,))
        handler.add_method("set_black_hole_mass", (nbody_system.mass, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_random_seed", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_random_seed", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_output_path", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_output_path", (handler.NO_UNIT,), (handler.ERROR_CODE,))
        
        handler.add_method("get_output_basename", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_output_basename", (handler.NO_UNIT,), (handler.ERROR_CODE,))

    def define_converter(self, handler):
        if not self.unit_converter is None:
            handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())
    
    def define_particle_sets(self, handler):
        handler.define_set('particles', 'index_of_the_particle')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_getter('particles', 'get_mass')
        handler.add_getter('particles', 'get_position')
        handler.add_getter('particles', 'get_velocity')
    
    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        handler.add_method('UPDATE','before_get_parameter')
        
        handler.add_transition('EDIT', 'UPDATE', 'generate_particles', False)
        handler.add_transition('UPDATE', 'RUN', 'update_particle_set')
        handler.add_transition('RUN', 'EDIT', 'clear_particle_set')
        handler.add_method('RUN', 'invoke_state_change_updated')
        handler.add_method('EDIT', 'get_number_of_particles_updated')
        handler.add_method('UPDATE', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_mass')
        handler.add_method('RUN', 'get_position')
        handler.add_method('RUN', 'get_velocity')
    
    def generate_particles(self):
        result = self.overridden().generate_particles()
        self.invoke_state_change_updated()
    
    def update_particle_set(self):
        """
        update the particle set after changes in the code
        
        this implementation needs to move to the
        amuse.datamodel.incode_storage module, as
        it uses a lot of internal methods and info!
        """
        number_of_updated_particles = self.get_number_of_particles_updated()
        if number_of_updated_particles:
            self.particles._private.attribute_storage._add_indices(
                range(number_of_updated_particles)
            )
    
    def clear_particle_set(self):
        if len(self.particles):
            self.particles.remove_particles(self.particles)
    

