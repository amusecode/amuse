import os.path
from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

class BonsaiInterface2(CodeInterface, LiteratureReferencesMixIn, GravitationalDynamicsInterface, 
        StoppingConditionInterface, CodeWithDataDirectories):
    """
        .. [#] Bedorf J., Gaburov E., Fujii M. S., Nitadori K. Ishiyama T., Portegies Zwart S.,
        ...[#] "24.77 Pflops on a gravitational tree-code to simulate the Milky Way Galaxy
        ...[#] with 18600 GPUs", 2014, SC'14 proceedings, 54-65. https://doi.org/10.1109/SC.2014.10

        .. [#] Bedorf J., Gaburov E., Portegies Zwart S., "A sparse octree
        .. [#] gravitational N-body code that runs entirely on the GPU processor",
        .. [#] 2012, JCoPh, 231, 2825
    """

    include_headers = ['worker_code.h', 'stopcond.h']
    
    def name_of_worker(self,mode):
        return 'bonsai_worker'
    
    def __init__(self, mode=None, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
        self.set_src_directory(
            os.path.join(self.amuse_root_directory, 'src', 'amuse', 'community', 'bonsai2', 'src', '')
        )
    
    @legacy_function
    def set_src_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('src_directory', dtype='string', direction=function.IN,
            description = "The path to the Bonsai2 src directory.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "The current timestep for the system")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_mass():
        """
        Update the mass of a particle. Mass is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """        

        return function     
        
    @legacy_function
    def set_state():
        """
        Update the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, position and velocity) is updated.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The new radius of the particle")  
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function
  
        
    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        """
        return function       
        
    @legacy_function
    def set_velocity():
        """
        Set the velocity vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The current x component of the velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The current y component of the velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The current z component of the velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
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
    def set_acceleration():
        """
        Set the velocity vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('ax', dtype='float64', direction=function.IN, description = "The current x component of the velocity vector of the particle")
        function.addParameter('ay', dtype='float64', direction=function.IN, description = "The current y component of the velocity vector of the particle")
        function.addParameter('az', dtype='float64', direction=function.IN, description = "The current z component of the velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
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
    def get_theta_for_tree():
        """
        Get theta, the opening angle for building the tree: between 0 and 1.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('theta_for_tree', dtype='float64', direction=function.OUT,
            description = "theta, the opening angle for building the tree: between 0 and 1")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_theta_for_tree():
        """
        Set theta, the opening angle for building the tree: between 0 and 1.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('theta_for_tree', dtype='float64', direction=function.IN,
            description = "theta, the opening angle for building the tree: between 0 and 1")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function


#    @legacy_function
#    def get_eta():
#        """
#        Get eta, the block time-step parameter.
#        """
#        function = LegacyFunctionSpecification()
#        function.addParameter('eta', dtype='float64', direction=function.OUT,
#            description = "Eta, time-step parameter")
#        function.result_type = 'int32'
#        return function
#        
#    @legacy_function
#    def set_eta():
#        """
#        Set eta, the block time-step parameter.
#        """
#        function = LegacyFunctionSpecification()
#        function.addParameter('eta', dtype='float64', direction=function.IN,
#            description = "Eta, time-step parameter")
#        function.result_type = 'int32'
#        return function

class Bonsai2(GravitationalDynamics):
    
    def __init__(self, unit_converter = None, **options):
        self.stopping_conditions = StoppingConditions(self)
        
        
        legacy_interface = BonsaiInterface2(**options)
        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            unit_converter,
            **options
        )
    
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_eps2", 
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.05**2 | nbody_system.length**2
        )
        
        handler.add_method_parameter(
            "get_time_step", 
            "set_time_step", 
            "timestep", 
            "timestep for the system", 
            default_value = 1.0 / 64 | nbody_system.time
        )

#        handler.add_method_parameter(
#            "get_eta", 
#            "set_eta", 
#            "timestep_parameter", 
#            "timestep parameter for the block time-step", 
#            default_value = 0.1
#        )

        handler.add_method_parameter(
            "get_theta_for_tree",
            "set_theta_for_tree",
            "opening_angle", 
            "opening angle, theta, for building the tree: between 0 and 1", 
            default_value = 0.75
        )
 
        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )
        self.stopping_conditions.define_parameters(handler)
    
    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        handler.add_method(
            "set_time",
            (nbody_system.time,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_time_step",
            (nbody_system.time,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length**2, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_time",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_eps2",
            (nbody_system.length**2,),
            (handler.ERROR_CODE,)
        )
        self.stopping_conditions.define_methods(handler)
        
    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        
        self.stopping_conditions.define_particle_set(handler)
    
    def define_errorcodes(self, handler):
        handler.add_errorcode(-1, 'Unspecified, other error.')
        handler.add_errorcode(-2, 'Called function is not implemented.')
        handler.add_errorcode(-3, 'A particle with the given index was not found.')
        handler.add_errorcode(-4, 'The tree has become too deep, consider the removal of far away particles to prevent a too large box.')
    

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
          
        self.stopping_conditions.define_state(handler)         

