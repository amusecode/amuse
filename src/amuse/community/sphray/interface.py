import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option


class SPHRayInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """    
    
    SPHRAY is a smoothed particle hydrodynamics (SPH) ray tracer designed to
    solve the 3D, time-dependent, radiative transfer equation. SPHRAY relies on 
    a Monte Carlo (MC) ray-tracing scheme that does not interpolate the SPH 
    particles on to a grid but instead integrates directly through the SPH
    kernels.  A quick Axis Aligned Bounding Box (AABB) test taken from computer
    graphics applications allows for the acceleration of the ray-tracing
    component. 

    The relevant references are:
        .. [#] Altay, Gabriel; Croft, Rupert A. C.; Pelupessy, Inti, 2008, MNRAS 386
    """

    def __init__(self,  **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)
        
    
    @option(type="string", sections=('data',))
    def input_data_root_directory(self):
        """
        The root directory of the input data, read only directories
        """
        return os.path.join(get_amuse_root_dir(), 'data')
        
    @option(type="string", sections=('data',))
    def output_data_root_directory(self):
        """
        The root directory of the output data,
        read - write directory
        """
        return os.path.join(get_amuse_root_dir(), 'data')
        
    def name_of_the_worker(self):
        return 'sphray_worker'
        
    def data_directory(self):
        """
        Returns the root name of the directory for the 
        application data files.
        """
        return os.path.join(self.input_data_root_directory, 'sphray', 'input')

    def output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, 'sphray', 'output')

    def new_particle(self, mass, radius, x, y, z, vx, vy, vz):
        pass


    @legacy_function    
    def commit_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def recommit_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_gas_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT)
        for x in ['mass','hsml','x','y','z','rho','xe','u']:
            function.addParameter(x, dtype='float32', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_src_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT)
        for x in ['L','x','y','z','SpcType']:
            function.addParameter(x, dtype='float32', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','hsml','x','y','z','rho','xe','u']:
            function.addParameter(x, dtype='float32', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state_src():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        for x in ['L','x','y','z','SpcType']:
            function.addParameter(x, dtype='float32', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_state_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','hsml','x','y','z','rho','xe','u']:
            function.addParameter(x, dtype='float32', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_state_src():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        for x in ['L','x','y','z','SpcType']:
            function.addParameter(x, dtype='float32', direction=function.IN)
        function.result_type = 'i'
        return function


    @legacy_function    
    def remove_src_particle():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def remove_gas_particle():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_isothermal():
        """ set_isothermal([0,1]): isothermal if 1, Temperature evolution if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('isothermal_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_isothermal():
        """ get_isothermal(): isothermal if 1, Temperature evolution if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('isothermal_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;


    @legacy_function    
    def evolve_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_data_directory():
        """
        Update the path to the sphray database.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('data_directory', dtype='string', direction=function.IN,
            description = "Name of the data directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function

    @legacy_function
    def get_data_directory():
        """
        Retrieve the path to the database currently used.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('data_directory', dtype='string', direction=function.OUT,
            description = "Name of the data directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Value was retrieved
        -1 - ERROR
            Could not retrieve value
        """
        return function

    @legacy_function
    def set_output_directory():
        """
        Update the output path.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('output_directory', dtype='string', direction=function.IN,
            description = "Name of the output directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function

    @legacy_function
    def get_output_directory():
        """
        Retrieve the output path.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('output_directory', dtype='string', direction=function.OUT,
            description = "Name of the output directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Value was retrieved
        -1 - ERROR
            Could not retrieve value
        """
        return function

class SPHRay(CommonCode):
        
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, SPHRayInterface(**options))
        self.set_data_directory(self.data_directory())
        self.set_output_directory(self.output_directory())  


