from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface
from amuse.community.interface.common import CommonCode
from amuse.units import units

class AsteriskInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
    Asterisk is a 3D visualization package for AMUSE simulations.
    
        .. [#] The Asterisk 3D visualization project is a collaboration between Sterrewacht Leiden and The Netherlands eScience Center.
    """

    classpath = 'src:src/lib/*:src/lib/jogl/*:src/dist/*:*:.'
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="asterisk_worker_java", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    @option(choices=['mpi','remote','ibis', 'sockets'], sections=("channel",))
    def channel_type(self):
        return 'sockets'
    
    @legacy_function
    def new_star_particle():
        """
        Define a new particle in the visualisation code. The particle is initialized with the provided
        radius, position and color. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.IN, 
                description = "The initial position vector of the particle")
        function.addParameter('radius', dtype='float64', unit=generic_unit_system.length, direction=function.IN, description = "The radius of the particle")
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The RGB color of the particle")
        function.addParameter("alpha", dtype='float64', direction=function.IN, description = "The opacity of the particle", default = 1.0)
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_gas_particle():
        """
        Define a new particle in the visualisation code. The particle is initialized with the provided
        radius, position and color. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.IN, 
                description = "The initial position vector of the particle")
        function.addParameter('radius', dtype='float64', unit=generic_unit_system.length, direction=function.IN, description = "The radius of the particle")
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The RGB color of the particle")
        function.addParameter("alpha", dtype='float64', direction=function.IN, description = "The opacity of the particle", default = 1.0)
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_sphere_particle():
        """
        Define a new particle in the visualisation code. The particle is initialized with the provided
        radius, position and color. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.IN, 
                description = "The initial position vector of the particle")
        function.addParameter('radius', dtype='float64', unit=generic_unit_system.length, direction=function.IN, description = "The radius of the particle")
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The RGB color of the particle")
        function.addParameter("alpha", dtype='float64', direction=function.IN, description = "The opacity of the particle", default = 1.0)
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_marker_particle():
        """
        Define a new particle in the visualisation code. The particle is initialized with the provided
        radius, position and color. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.IN, 
                description = "The initial position vector of the particle")
        function.addParameter('radius', dtype='float64', unit=generic_unit_system.length, direction=function.IN, description = "The radius of the particle")
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The RGB color of the particle")
        function.addParameter("alpha", dtype='float64', direction=function.IN, description = "The opacity of the particle", default = 1.0)
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        """
        Remove the definition of particle from the code. After calling this function the particle is
        no longer part of the model evolution. It is up to the code if the index will be reused.
        This function is optional.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        #function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to be removed. This index must have been returned by an earlier call to :meth:`new_particle`")

        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be removed
        -2 - ERROR
            not yet implemented
        """
        return function



    @legacy_function
    def get_radius():
        """
        Retrieve the radius of a particle. Radius is a scalar property of a particle,
        this function has one OUT argument.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        #function.can_handle_array = True
        function.must_handle_array = True
        return function


    @legacy_function
    def set_radius():
        """
        Set the radius of a particle. Radius is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('radius', dtype='float64', unit=generic_unit_system.length, direction=function.IN, description = "The new radius of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        #function.can_handle_array = True
        function.must_handle_array = True
        return function

    @legacy_function
    def commit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function

    @legacy_function
    def recommit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function

    
    @legacy_function
    def get_position():
        """
        Retrieve the position vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.OUT, 
                description = "The current position vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.IN, 
                description = "The new position vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def get_color():
        """
        Retrieve the RGB color vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.OUT, 
                description = "The current RGB color vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def set_color():
        """
        Update the RGB color of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The new RGB color vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def get_opacity():
        """
        Retrieve the alpha (opacity) of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter("alpha", dtype='float64', direction=function.OUT, 
            description = "The opacity of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def set_opacity():
        """
        Update the alpha (opacity) of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter("alpha", dtype='float64', direction=function.IN, 
            description = "The new opacity of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def get_use_star_shader_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_star_shader_flag", dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_use_star_shader_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_star_shader_flag", dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_use_octree_for_gas_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_octree_for_gas_flag", dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_use_octree_for_gas_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_octree_for_gas_flag", dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_rotation():
        function = LegacyFunctionSpecification()
        function.addParameter("x_rotation", dtype='float64', direction=function.IN)
        function.addParameter("y_rotation", dtype='float64', direction=function.IN)
        function.addParameter("z_rotation", dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_current_rotation():
        function = LegacyFunctionSpecification()
        function.addParameter("x_rotation", dtype='float64', direction=function.OUT)
        function.addParameter("y_rotation", dtype='float64', direction=function.OUT)
        function.addParameter("z_rotation", dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_translation():
        function = LegacyFunctionSpecification()
        function.addParameter("x_translation", dtype='float64', direction=function.IN)
        function.addParameter("y_translation", dtype='float64', direction=function.IN)
        function.addParameter("z_translation", dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_current_translation():
        function = LegacyFunctionSpecification()
        function.addParameter("x_translation", dtype='float64', direction=function.OUT)
        function.addParameter("y_translation", dtype='float64', direction=function.OUT)
        function.addParameter("z_translation", dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_camera_distance():
        function = LegacyFunctionSpecification()
        function.addParameter("camera_distance", dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_camera_distance():
        function = LegacyFunctionSpecification()
        function.addParameter("camera_distance", dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_scene_number():
        """
        Get the number of the currently displayed scene, which can be used to redisplay it later
        """
        function = LegacyFunctionSpecification()
        function.addParameter("scene_number", dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_scene():
        """
        Display the specified scene
        """
        function = LegacyFunctionSpecification()
        function.addParameter("scene_number", dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def screenshot():
        function = LegacyFunctionSpecification()  
        function.addParameter('screenshot_file_name', dtype='string', direction=function.IN,
            description = "Name of the file to write the (PNG) screenshot to")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def store_view():
        """
        Store and view the current model, corresponding to the given description.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='string', direction=function.IN,
            description = "The description of the scene.")
        function.result_type = 'int32'
        return function
    
class Asterisk(CommonCode):

    def __init__(self, unit_converter=None, **options):
        self.unit_converter = unit_converter
        CommonCode.__init__(self,  AsteriskInterface(**options), **options)
        
    def store_view(self, description=""):
        self.overridden().store_view(str(description))

    def define_properties(self, object):
        object.add_property("get_current_rotation")
        object.add_property("get_current_translation")
        object.add_property("get_current_rotation", public_name = "rotation")
        object.add_property("get_current_translation", public_name = "translation")
    
    def define_particle_sets(self, object):
        object.define_super_set('particles', ['star_particles', 'gas_particles', 'sphere_particles', 'marker_particles'], 
            index_to_default_set = 0)
        
        object.define_set('star_particles', 'index_of_the_particle')
        object.set_new('star_particles', 'new_star_particle')
        object.set_delete('star_particles', 'delete_particle')
        
        object.define_set('gas_particles', 'index_of_the_particle')
        object.set_new('gas_particles', 'new_gas_particle')
        object.set_delete('gas_particles', 'delete_particle')
        
        object.define_set('sphere_particles', 'index_of_the_particle')
        object.set_new('sphere_particles', 'new_sphere_particle')
        object.set_delete('sphere_particles', 'delete_particle')
        
        object.define_set('marker_particles', 'index_of_the_particle')
        object.set_new('marker_particles', 'new_marker_particle')
        object.set_delete('marker_particles', 'delete_particle')
        
        for particle_set_name in ['star_particles', 'gas_particles', 'sphere_particles', 'marker_particles']:
            object.add_getter(particle_set_name, 'get_position')
            object.add_setter(particle_set_name, 'set_position')
            object.add_getter(particle_set_name, 'get_color')
            object.add_setter(particle_set_name, 'set_color')
            object.add_getter(particle_set_name, 'get_opacity')
            object.add_setter(particle_set_name, 'set_opacity')
            object.add_getter(particle_set_name, 'get_radius')
            object.add_setter(particle_set_name, 'set_radius')

    def define_state(self, object): 
        CommonCode.define_state(self, object)   
        object.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        object.add_method('UPDATE','before_get_parameter')
        
        object.add_method('EDIT', 'new_star_particle')
        object.add_method('EDIT', 'new_gas_particle')
        object.add_method('EDIT', 'new_sphere_particle')
        object.add_method('EDIT', 'new_marker_particle')
        object.add_method('EDIT', 'delete_particle')
        object.add_method('UPDATE', 'new_star_particle')
        object.add_method('UPDATE', 'new_gas_particle')
        object.add_method('UPDATE', 'new_sphere_particle')
        object.add_method('UPDATE', 'new_marker_particle')
        object.add_method('UPDATE', 'delete_particle')
        object.add_transition('EDIT', 'RUN', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_star_particle', False)
        object.add_transition('RUN', 'UPDATE', 'new_gas_particle', False)
        object.add_transition('RUN', 'UPDATE', 'new_sphere_particle', False)
        object.add_transition('RUN', 'UPDATE', 'new_marker_particle', False)
        object.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        
        object.add_method('RUN', 'store_view')
        object.add_method('RUN', 'set_position')
        object.add_method('RUN', 'set_color')
        object.add_method('RUN', 'set_opacity')
        object.add_method('RUN', 'get_position')
        object.add_method('RUN', 'get_color')
        object.add_method('RUN', 'get_opacity')
        
    def define_parameters(self, object):
        object.add_boolean_parameter(
            "get_use_star_shader_flag",
            "set_use_star_shader_flag",
            "use_star_shader",
            "Use-star-shader flag. False means: plain spheres.",
            True
        )
        object.add_boolean_parameter(
            "get_use_octree_for_gas_flag",
            "set_use_octree_for_gas_flag",
            "use_octree_for_gas",
            "Use-octree-for-gas flag. True means: gas particles are divided over "
                "octree cells, and these cells will be visualized instead.",
            False
        )
        
        object.add_caching_parameter(
            "set_rotation", 
            "x_rotation",
            "rotation_x", 
            "Rotation of the scene about the x axis (degrees)",
            0.0,
        )
        object.add_caching_parameter(
            "set_rotation", 
            "y_rotation",
            "rotation_y", 
            "Rotation of the scene about the y axis (degrees)",
            0.0,
        )
        object.add_caching_parameter(
            "set_rotation", 
            "z_rotation",
            "rotation_z", 
            "Rotation of the scene about the z axis (degrees)",
            0.0,
        )
        object.add_vector_parameter(
            "rotation",
            "Rotation of the scene about the x, y, and z axes (degrees)",
            ("rotation_x", "rotation_y","rotation_z")
        )
        object.add_caching_parameter(
            "set_translation", 
            "x_translation",
            "translation_x", 
            "Translation of the scene, corresponding to a rotation of the scene w.r.t. the view point (degrees)",
            0.0,
        )
        object.add_caching_parameter(
            "set_translation", 
            "y_translation",
            "translation_y", 
            "Translation of the scene, corresponding to a rotation of the scene w.r.t. the view point (degrees)",
            0.0,
        )
        object.add_caching_parameter(
            "set_translation", 
            "z_translation",
            "translation_z", 
            "Translation of the scene, corresponding to a rotation of the scene w.r.t. the view point (degrees)",
            0.0,
        )
        object.add_vector_parameter(
            "translation",
            "Translation of the scene, corresponding to a rotation of the scene w.r.t. the view point (degrees)",
            ("translation_x", "translation_y","translation_z")
        )
        object.add_method_parameter(
            "get_camera_distance",
            "set_camera_distance", 
            "camera_distance", 
            "Distance from the view point to the scene", 
            default_value = 0.0 | nbody_system.length
        )
        object.add_method_parameter(
            "get_scene_number",
            "set_scene", 
            "scene", 
            "set: Set the scene to display; get: Get the current scene number, which can be used to display it later", 
            default_value = 0
        )
    
    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
            
