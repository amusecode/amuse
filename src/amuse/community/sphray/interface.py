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
    def set_H_caseA():
        """ set_H_caseA([0,1]): use case A for H if 1, not if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('H_caseA_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_H_caseA():
        """ get_H_caseA(): use case A for H if 1, not if 0  """
        function = LegacyFunctionSpecification()  
        function.addParameter('H_caseA_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def set_He_caseA():
        """ set_He_caseA([0,1]): use case A for H if 1, not if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('He_caseA_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_He_caseA():
        """ get_He_caseA(): use case A for H if 1, not if 0  """
        function = LegacyFunctionSpecification()  
        function.addParameter('He_caseA_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_raynumber():
        """ set number of rays per evolve """
        function = LegacyFunctionSpecification()  
        function.addParameter('number_rays', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_raynumber():
        """ get number of rays per evolve  """
        function = LegacyFunctionSpecification()  
        function.addParameter('number_rays', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;
        
    @legacy_function   
    def set_iontempsolver():
        """ set solver to use (1,2 = euler, bdf) """
        function = LegacyFunctionSpecification()  
        function.addParameter('iontempsolver', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_iontempsolver():
        """ get solver to use (1,2 = euler, bdf)  """
        function = LegacyFunctionSpecification()  
        function.addParameter('iontempsolver', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_boundary():
        """ set boundary condition to use (-1,0,1 = reflective, vacuum, periodic) """
        function = LegacyFunctionSpecification()  
        function.addParameter('boundary', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_boundary():
        """ get boundary condition to use (-1,0,1 = reflective, vacuum, periodic)  """
        function = LegacyFunctionSpecification()  
        function.addParameter('boundary', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_boxsize():
        """ set box size in kpc  """
        function = LegacyFunctionSpecification()  
        function.addParameter('boxsize', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_boxsize():
        """ get box size in kpc   """
        function = LegacyFunctionSpecification()  
        function.addParameter('boxsize', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_globalHefraction():
        """ set He mass fraction (f_H+f_He=1 )  """
        function = LegacyFunctionSpecification()  
        function.addParameter('globalHefraction', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_globalHefraction():
        """ get He mass fraction (f_H+f_He=1 )   """
        function = LegacyFunctionSpecification()  
        function.addParameter('globalHefraction', dtype='d', direction=function.OUT)
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

    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")

    def define_parameters(self, object):
        object.add_boolean_parameter(
            "get_isothermal",
            "set_isothermal", 
            "isothermal_flag", 
            "whether to evolve temperature (0) or keep isothermal (1)", 
            default_value = False
        )

        object.add_method_parameter(
            "get_raynumber",
            "set_raynumber", 
            "number_of_rays", 
            "number of rays per evolve step", 
            default_value = 30000
        )

        object.add_method_parameter(
            "get_iontempsolver",
            "set_iontempsolver", 
            "ionization_temperature_solver", 
            "solver to use (1: euler, 2: bdf)", 
            default_value = 2
        )        

        object.add_method_parameter(
            "get_globalHefraction",
            "set_globalHefraction", 
            "global_helium_mass_Fraction", 
            "global helium mass fraction (f_H+f_He=1)", 
            default_value = 0.
        )

        object.add_method_parameter(
            "get_boxsize",
            "set_boxsize", 
            "box_size", 
            "simulation box size", 
            default_value = 13.2 | units.kpc
        )

        object.add_method_parameter(
            "get_boundary",
            "set_boundary", 
            "boundary_condition", 
            "simulation bundary condition (-1,0,1 = reflective, vacuum, periodic", 
            default_value = 0
        )

        object.add_boolean_parameter(
            "get_H_caseA",
            "set_H_caseA", 
            "hydrogen_case_A", 
            "flag for hydrogen case A recombination (1)", 
            default_value = True
        )

        object.add_boolean_parameter(
            "get_He_caseA",
            "set_He_caseA", 
            "helium_case_A", 
            "flag for  helium case A recombination (1)", 
            default_value = True
        )

