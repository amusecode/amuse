import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option


class SPHRayInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn, 
        CodeWithDataDirectories):
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
        CodeWithDataDirectories.__init__(self)
    
    def name_of_the_worker(self):
        return 'sphray_worker'
        
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
        for x in ['mass','h_smooth','x','y','z','rho','xion','u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='float64', direction=function.IN,default=0)
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_src_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT)
        for x in ['luminosity','x','y','z']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter("SpcType", dtype='float64', direction=function.IN,default=0)    
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','h_smooth','x','y','z','rho','xion','u']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state_src():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        for x in ['luminosity','x','y','z','SpcType']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_state_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','h_smooth','x','y','z','rho','xion','u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_pos_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_hsml_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['h_smooth']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_rho_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['rho']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_u_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_dudt_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['du_dt']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_vel_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_vel_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_dudt_gas():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['du_dt']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_state_src():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        for x in ['luminosity','x','y','z']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter("SpcType", dtype='float64', direction=function.IN,default=0.)    
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_luminosity_src():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        for x in ['luminosity']:
            function.addParameter(x, dtype='float64', direction=function.IN)
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
        return function
    @legacy_function   
    def get_isothermal():
        """ get_isothermal(): isothermal if 1, Temperature evolution if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('isothermal_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function   
    def set_H_caseA():
        """ set_H_caseA([0,1]): use case A for H if 1, not if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('H_caseA_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_H_caseA():
        """ get_H_caseA(): use case A for H if 1, not if 0  """
        function = LegacyFunctionSpecification()  
        function.addParameter('H_caseA_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function   
    def set_He_caseA():
        """ set_He_caseA([0,1]): use case A for H if 1, not if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('He_caseA_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_He_caseA():
        """ get_He_caseA(): use case A for H if 1, not if 0  """
        function = LegacyFunctionSpecification()  
        function.addParameter('He_caseA_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_raynumber():
        """ set number of rays per unit time """
        function = LegacyFunctionSpecification()  
        function.addParameter('number_rays', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_raynumber():
        """ get number of rays per unit time  """
        function = LegacyFunctionSpecification()  
        function.addParameter('number_rays', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_defaultspectype():
        """ set default specType (negative in units of rydbergs, or positive integer) """
        function = LegacyFunctionSpecification()  
        function.addParameter('defaultspectype', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_defaultspectype():
        """ get default specType (negative in units of rydbergs, or positive integer) """
        function = LegacyFunctionSpecification()  
        function.addParameter('defaultspectype', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

        
    @legacy_function   
    def set_iontempsolver():
        """ set solver to use (1,2 = euler, bdf) """
        function = LegacyFunctionSpecification()  
        function.addParameter('iontempsolver', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_iontempsolver():
        """ get solver to use (1,2 = euler, bdf)  """
        function = LegacyFunctionSpecification()  
        function.addParameter('iontempsolver', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_boundary():
        """ set boundary condition to use (-1,0,1 = reflective, vacuum, periodic) """
        function = LegacyFunctionSpecification()  
        function.addParameter('boundary', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_boundary():
        """ get boundary condition to use (-1,0,1 = reflective, vacuum, periodic)  """
        function = LegacyFunctionSpecification()  
        function.addParameter('boundary', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_boxsize():
        """ set box size in kpc  """
        function = LegacyFunctionSpecification()  
        function.addParameter('boxsize', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_boxsize():
        """ get box size in kpc   """
        function = LegacyFunctionSpecification()  
        function.addParameter('boxsize', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_globalHefraction():
        """ set He mass fraction (f_H+f_He=1 )  """
        function = LegacyFunctionSpecification()  
        function.addParameter('globalHefraction', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_globalHefraction():
        """ get He mass fraction (f_H+f_He=1 )   """
        function = LegacyFunctionSpecification()  
        function.addParameter('globalHefraction', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_momentum_kicks():
        """ set_momentume_kicks([0,1]): calc momentum kicks if 1, not if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('momentum_kicks', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_momentum_kicks():
        """ set_momentume_kicks([0,1]): calc momentum kicks if 1, not if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('momentum_kicks', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

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
    def set_sphray_data_directory():
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
    def get_sphray_data_directory():
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
    def set_spectra_file():
        """
        set the spectrafile to use(relative to data directory)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('spectra_file', dtype='string', direction=function.IN,
            description = "Name of the file")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function

    @legacy_function
    def get_spectra_file():
        """
        get the spectrafile to use(relative to data directory)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('spectra_file', dtype='string', direction=function.OUT,
            description = "Name of the file")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Value was retrieved
        -1 - ERROR
            Could not retrieve value
        """
        return function

    @legacy_function
    def set_sphray_output_directory():
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
    def get_sphray_output_directory():
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
        
    def __init__(self,unit_converter = None, **options):
        
        if unit_converter is not None:
            raise Exception("atm SPHRay uses predefined units converter")
           
        self.unit_converter = ConvertBetweenGenericAndSiUnits(
                1. | units.kpc,   
                1.e10 | units.MSun,
                1. | units.kms)

        InCodeComponentImplementation.__init__(self, SPHRayInterface(**options))
        self.set_sphray_data_directory(self.data_directory)
        self.set_sphray_output_directory(self.output_directory)
        
    def define_converter(self, handler):
        if self.unit_converter is None:
            raise Exception("something seriously wrong")
        
        handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())
          
    def define_properties(self, handler):
        handler.add_property('get_time', public_name = "model_time")

    def define_parameters(self, handler):
        handler.add_boolean_parameter(
            "get_isothermal",
            "set_isothermal", 
            "isothermal_flag", 
            "whether to evolve temperature (0) or keep isothermal (1)", 
            False
        )

        handler.add_boolean_parameter(
            "get_momentum_kicks",
            "set_momentum_kicks", 
            "momentum_kicks_flag", 
            "whether to use momentum kicks (1) not (0)", 
            False
        )

        handler.add_method_parameter(
            "get_raynumber",
            "set_raynumber", 
            "number_of_rays", 
            "number of rays per unit time", 
            default_value = 1000 | generic_unit_system.time**-1
        )

        handler.add_method_parameter(
            "get_iontempsolver",
            "set_iontempsolver", 
            "ionization_temperature_solver", 
            "solver to use (1: euler, 2: bdf)", 
            default_value = 2
        )        

        handler.add_method_parameter(
            "get_globalHefraction",
            "set_globalHefraction", 
            "global_helium_mass_Fraction", 
            "global helium mass fraction (f_H+f_He=1)", 
            default_value = 0.
        )

        handler.add_method_parameter(
            "get_boxsize",
            "set_boxsize", 
            "box_size", 
            "simulation box size", 
            default_value = 13.2 | units.kpc
        )

        handler.add_method_parameter(
            "get_defaultspectype",
            "set_defaultspectype", 
            "default_spectral_type", 
            "default src spectral type", 
            default_value = -1.
        )

        handler.add_method_parameter(
            "get_boundary",
            "set_boundary", 
            "boundary_condition", 
            "simulation bundary condition (-1,0,1 = reflective, vacuum, periodic", 
            default_value = 0
        )

        handler.add_boolean_parameter(
            "get_H_caseA",
            "set_H_caseA", 
            "hydrogen_case_A_flag", 
            "flag for hydrogen case A recombination (1)", 
            True
        )

        handler.add_boolean_parameter(
            "get_He_caseA",
            "set_He_caseA", 
            "helium_case_A_flag", 
            "flag for  helium case A recombination (1)", 
            True
        )

        handler.add_method_parameter(
            "get_spectra_file",
            "set_spectra_file", 
            "spectra_file", 
            "spectra file to use", 
            default_value = "./spectra/thermal30k.cdf"
        )

    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method('evolve_model', (generic_unit_system.time,), ( handler.ERROR_CODE, ))
        handler.add_method(
            "new_gas_particle",
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.density,
                handler.NO_UNIT,
                generic_unit_system.speed**2,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_state_gas",
            (
                handler.NO_UNIT,
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.density,
                handler.NO_UNIT,
                generic_unit_system.speed**2
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_vel_gas",
            (
                handler.NO_UNIT,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_pos_gas",
            (
                handler.NO_UNIT,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_rho_gas",
            (
                handler.NO_UNIT,
                generic_unit_system.density,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_hsml_gas",
            (
                handler.NO_UNIT,
                generic_unit_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_u_gas",
            (
                handler.NO_UNIT,
                generic_unit_system.speed**2
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_dudt_gas",
            (
                handler.NO_UNIT,
                generic_unit_system.length**2/generic_unit_system.time**3
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_dudt_gas",
            (
                handler.INDEX,
            ),
            (
                generic_unit_system.length**2/generic_unit_system.time**3,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state_gas",
            (
                handler.INDEX,
            ),
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.density,
                handler.NO_UNIT,
                generic_unit_system.speed**2,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_vel_gas",
            (
                handler.INDEX,
            ),
            (
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "new_src_particle",
            (
                1.e50 * units.s**-1,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                handler.NO_UNIT,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_state_src",
            (
                handler.NO_UNIT,
                1.e50 * units.s**-1,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_luminosity_src",
            (
                handler.NO_UNIT,
                1.e50 * units.s**-1,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state_src",
            (
                handler.INDEX,
            ),
            (
                1.e50 * units.s**-1,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "remove_src_particle",
            (
                handler.INDEX,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "remove_gas_particle",
            (
                handler.INDEX,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_boxsize",
            (
                generic_unit_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_defaultspectype",
            (
                handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_globalHefraction",
            (
               handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_boxsize",
            (
            ),
            (
                generic_unit_system.length,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_globalHefraction",
            (
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_defaultspectype",
            (
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_raynumber",
            (
               generic_unit_system.time**-1,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_raynumber",
            (
            ),
            (
                generic_unit_system.time**-1,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_iontempsolver",
            (
               handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_iontempsolver",
            (
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_boundary",
            (
               handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_boundary",
            (
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_time",
            (
            ),
            (
                generic_unit_system.time,
                handler.ERROR_CODE,
            )
        )



    def define_particle_sets(self, handler):
        handler.define_set('gas_particles', 'id')
        handler.set_new('gas_particles', 'new_gas_particle')
        handler.set_delete('gas_particles', 'remove_gas_particle')
        handler.add_setter('gas_particles', 'set_state_gas')
        handler.add_getter('gas_particles', 'get_state_gas')
        handler.add_setter('gas_particles', 'set_pos_gas')
        handler.add_setter('gas_particles', 'set_vel_gas')
        handler.add_getter('gas_particles', 'get_vel_gas')
        handler.add_setter('gas_particles', 'set_hsml_gas')
        handler.add_setter('gas_particles', 'set_rho_gas')
        handler.add_setter('gas_particles', 'set_u_gas')
        handler.add_setter('gas_particles', 'set_dudt_gas')
        handler.add_getter('gas_particles', 'get_dudt_gas')
        
                        
        handler.define_set('src_particles', 'id')
        handler.set_new('src_particles', 'new_src_particle')
        handler.set_delete('src_particles', 'remove_src_particle')
        handler.add_setter('src_particles', 'set_state_src')
        handler.add_setter('src_particles', 'set_luminosity_src')
        handler.add_getter('src_particles', 'get_state_src')
        
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
        handler.add_method('EVOLVED','before_get_parameter')
        
        handler.add_method('EDIT', 'new_gas_particle')
        handler.add_method('EDIT', 'remove_gas_particle')
        handler.add_method('EDIT', 'new_src_particle')
        handler.add_method('EDIT', 'remove_src_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_gas_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'remove_gas_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'new_src_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'remove_src_particle', False)
        
        handler.add_transition('RUN', 'UPDATE', 'set_pos_gas', False)
        handler.add_transition('RUN', 'UPDATE', 'set_rho_gas', False)
        handler.add_transition('RUN', 'UPDATE', 'set_hsml_gas', False)
        handler.add_transition('RUN', 'UPDATE', 'set_u_gas', False)
        handler.add_transition('RUN', 'UPDATE', 'set_dudt_gas', False)
        handler.add_transition('RUN', 'UPDATE', 'set_state_gas', False)
        handler.add_transition('RUN', 'UPDATE', 'set_state_src', False)
        handler.add_transition('RUN', 'UPDATE', 'set_luminosity_src', False)
        
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_method('RUN', 'evolve_model')
        handler.add_method('RUN', 'get_state_gas')
        handler.add_method('RUN', 'get_state_src')

        handler.add_method('INITIALIZED', 'set_momentum_kicks')
        handler.add_method('INITIALIZED', 'set_isothermal')
        handler.add_method('INITIALIZED', 'set_boxsize')
        handler.add_method('INITIALIZED', 'set_globalHefraction')
        handler.add_method('INITIALIZED', 'set_raynumber')
        handler.add_method('INITIALIZED', 'set_iontempsolver')
        handler.add_method('INITIALIZED', 'set_boundary')
        handler.add_method('INITIALIZED', 'set_H_caseA')
        handler.add_method('INITIALIZED', 'set_He_caseA')


Sphray = SPHRay
