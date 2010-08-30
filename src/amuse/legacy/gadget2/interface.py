import os
import numpy
from amuse.legacy.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics
from amuse.legacy import *

class Gadget2Interface(LegacyInterface, GravitationalDynamicsInterface, LiteratureRefs, StoppingConditionInterface):
    """
    GADGET-2 computes gravitational forces with a hierarchical tree 
    algorithm (optionally in combination with a particle-mesh 
    scheme for long-range gravitational forces, currently not 
    supported from the AMUSE interface) and represents fluids by 
    means of smoothed particle hydrodynamics (SPH). The code can 
    be used for studies of isolated systems, or for simulations 
    that include the cosmological expansion of space, both with 
    or without periodic boundary conditions. In all these types 
    of simulations, GADGET follows the evolution of a self-
    gravitating collisionless N-body system, and allows gas 
    dynamics to be optionally included. Both the force computation 
    and the time stepping of GADGET are fully adaptive, with a 
    dynamic range which is, in principle, unlimited. 
    
    The relevant references are:
        .. [#] Springel V., 2005, MNRAS, 364, 1105  (GADGET-2)
        .. [#] Springel V., Yoshida N., White S. D. M., 2001, New Astronomy, 6, 51  (GADGET-1)
    """
    include_headers = ['gadget_code.h', 'worker_code.h', 'stopcond.h']
    
    def __init__(self, **options):
        LegacyInterface.__init__(self, name_of_the_worker = 'worker_code', **options)
        LiteratureRefs.__init__(self)
    
    def get_data_directory(self):
        """
        Returns the root name of the directory for the Gadget2
        application data files.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'gadget2', 'input')
    
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'gadget2', 'output')
    
    @property
    def default_path_to_parameterfile(self):
        return os.path.join(self.get_data_directory(), 'defaults.param')
    
    @legacy_function
    def set_parameterfile_path():
        function = LegacyFunctionSpecification()
        function.addParameter('parameterfile_path', dtype='string', direction=function.IN,
            description = "The path to the default Gadget-2 Parameterfile.")
        function.result_type = 'int32'
        return function;
    
    @legacy_function
    def new_dm_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        for x in ['mass','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    def new_particle(self, mass, x, y, z, vx, vy, vz):
        return self.new_dm_particle(mass, x, y, z, vx, vy, vz)
    
    @legacy_function
    def new_sph_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_state():
        """
        Retrieve the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, position and velocity) is returned.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retrieved
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
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function
    @legacy_function    
    def get_state_sph():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_state_sph():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('u', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('u', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_thermal_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('thermal_energy', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    
    
# setting/ getting parameters
    @legacy_function
    def set_time_step():
        """ timestep (code units)"""
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "The current model timestep")
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_epsilon():
        """Get epsilon, a softening parameter for gravitational potentials with point particles."""
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.OUT,
            description = "epsilon, a softening parameter for gravitational potentials with point particles")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_epsilon():
        """Set epsilon, a softening parameter for gravitational potentials with point particles."""
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.IN,
            description = "epsilon, a softening parameter for gravitational potentials with point particles")
        function.result_type = 'int32'
        return function
    def set_epsilon_squared(self,epsilon_squared):
        return self.set_epsilon(epsilon_squared**0.5)
    def get_epsilon_squared(self):
        epsilon, err = self.get_epsilon()
        return epsilon**2, err
    @legacy_function
    def get_gadget_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('gadget_output_directory', dtype='string', direction=function.OUT,
            description = "The path to the Gadget-2 OutputDir.")
        function.result_type = 'int32'
        return function;
    @legacy_function
    def set_gadget_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('gadget_output_directory', dtype='string', direction=function.IN,
            description = "The path to the Gadget-2 OutputDir.")
        function.result_type = 'int32'
        return function;
    @legacy_function
    def get_unit_mass():
        """Get the code mass unit (in g/h, default: 1.989e43 g = 10^10 MSun)."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_mass_unit', dtype='float64', direction=function.OUT,
            description = "The code mass unit (in g/h, default: 1.989e43 g = 10^10 MSun).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_unit_mass():
        """Set the code mass unit (in g/h, default: 1.989e43 g = 10^10 MSun)."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_mass_unit', dtype='float64', direction=function.IN,
            description = "The code mass unit (in g/h, default: 1.989e43 g = 10^10 MSun).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_unit_length():
        """Get the code length unit (in cm/h, default: 3.085678e21 cm = 1 kpc)."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_length_unit', dtype='float64', direction=function.OUT,
            description = "The code length unit (in cm/h, default: 3.085678e21 cm = 1 kpc).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_unit_length():
        """Set the code length unit (in cm/h, default: 3.085678e21 cm = 1 kpc)."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_length_unit', dtype='float64', direction=function.IN,
            description = "The code length unit (in cm/h, default: 3.085678e21 cm = 1 kpc).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_unit_time():
        """Get the code time unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr). Implicitly changes velocity unit."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_time_unit', dtype='float64', direction=function.OUT,
            description = "The code time unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_unit_time():
        """Set the code time unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr)."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_time_unit', dtype='float64', direction=function.IN,
            description = "The code time unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_unit_velocity():
        """Get the code velocity unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr)."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_velocity_unit', dtype='float64', direction=function.OUT,
            description = "The code velocity unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_unit_velocity():
        """Set the code velocity unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr)."""
        function = LegacyFunctionSpecification()
        function.addParameter('code_velocity_unit', dtype='float64', direction=function.IN,
            description = "The code velocity unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr).")
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_nogravity():
        """ get_nogravity(): get no-gravity flag. True means: gravitational forces are switched of for all particles
        (read-only: makefile option NOGRAVITY)."""
        function = LegacyFunctionSpecification()
        function.addParameter('no_gravity_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_gdgop():
        """ set_gdgop([0,1]): use of gadget cell opening criterion if 1 """
        function = LegacyFunctionSpecification()
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function
    def get_gdgop():
        """ get_gdgop(): use of gadget cell opening criterion if 1 """
        function = LegacyFunctionSpecification()
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_epsgas():
        """ gas grav smoothing eps"""
        function = LegacyFunctionSpecification()
        function.addParameter('gas_epsilon', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function
    def get_epsgas():
        """ gas grav smoothing eps"""
        function = LegacyFunctionSpecification()
        function.addParameter('gas_epsilon', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def get_isotherm():
        """(default: False)
        True means: isothermal gas (read-only: makefile option ISOTHERM_EQS).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('isothermal_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def get_eps_is_h():
        """ get_eps_is_h(): gas particles grav. eps to SPH h if 1"""
        function = LegacyFunctionSpecification()
        function.addParameter('eps_is_h_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()
        function.addParameter('nsmooth', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function
    def get_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()
        function.addParameter('nsmooth', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """
        function = LegacyFunctionSpecification()
        function.addParameter('bh_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function
    def get_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """
        function = LegacyFunctionSpecification()
        function.addParameter('bh_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_gdgtol():
        """ Gadget cell openings criterion parameter  (unitless, .005) """
        function = LegacyFunctionSpecification()
        function.addParameter('gdgtol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function
    def get_gdgtol():
        """ Gadget cell openings criterion parameter  (unitless, .005) """
        function = LegacyFunctionSpecification()
        function.addParameter('gdgtol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_gamma():
        """ gas polytropic index (1.666667) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gamma', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('alpha', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('alpha', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('courant', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('courant', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
    @legacy_function
    def set_nsmtol():
        """ fractional tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()  
        function.addParameter('n_neighbour_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_nsmtol():
        """ fractional tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()  
        function.addParameter('n_neighbour_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
    

class Gadget2Doc(object):

    def __get__(self, instance, owner):
        return instance.legacy_interface.__doc__+"\n\n"+instance.parameters.__doc__

class Gadget2(GravitationalDynamics):
    
    __doc__ = Gadget2Doc()
    
    def __init__(self, unit_converter = None, **options):
        legacy_interface = Gadget2Interface(**options)
        if unit_converter is None:
            unit_converter = ConvertBetweenGenericAndSi(
                3.085678e21 | units.cm,   # 1.0 kpc
                1.989e43 | units.g,       # 1.0e10 solar masses
                1e5 | units.cm / units.s)# 1 km/sec

        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            unit_converter,
            **options
        )
    
    def initialize_code(self):
        self.set_parameterfile_path(self.default_path_to_parameterfile)
        result = self.overridden().initialize_code()
        self.set_gadget_output_directory(self.get_output_directory())
        self.set_unit_mass(self.unit_converter.to_si(generic_unit_system.mass).value_in(units.g))
        self.set_unit_length(self.unit_converter.to_si(generic_unit_system.length).value_in(units.cm))
        self.set_unit_velocity(self.unit_converter.to_si(
            generic_unit_system.length/generic_unit_system.time).value_in(units.cm/units.s))
        return result
    
    def define_properties(self, object):
        object.add_property("get_kinetic_energy", generic_unit_system.energy)
        object.add_property("get_potential_energy", generic_unit_system.energy)
        object.add_property("get_thermal_energy", generic_unit_system.energy)
        object.add_property("get_total_radius", generic_unit_system.length)
        object.add_property("get_center_of_mass_position", generic_unit_system.length)
        object.add_property("get_center_of_mass_velocity", generic_unit_system.speed)
        object.add_property("get_total_mass", generic_unit_system.mass)
        object.add_property('get_time', generic_unit_system.time, "model_time")
    
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        object.add_method('EDIT', 'new_dm_particle')
        object.add_transition('RUN', 'UPDATE', 'new_dm_particle', False)
        object.add_method('EDIT', 'new_sph_particle')
        object.add_transition('RUN', 'UPDATE', 'new_sph_particle', False)
        object.add_method('RUN', 'get_velocity')
        object.add_method('RUN', 'get_acceleration')
        object.add_method('RUN', 'get_internal_energy')
        
        object.add_method('RUN', 'get_kinetic_energy')
        object.add_method('RUN', 'get_potential_energy')
        object.add_method('RUN', 'get_thermal_energy')
        object.add_method('RUN', 'get_total_radius')
        object.add_method('RUN', 'get_center_of_mass_position')
        object.add_method('RUN', 'get_center_of_mass_velocity')
        object.add_method('RUN', 'get_total_mass')
        object.add_method('RUN', 'get_time')
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_epsilon_squared", 
            "set_epsilon_squared",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            generic_unit_system.length * generic_unit_system.length, 
            0.0001 | generic_unit_system.length * generic_unit_system.length
        )
        
        object.add_method_parameter(
            "get_time_step", 
            None,
            "timestep", 
            "timestep for system, Gadget2 calculates this by itself, based on particle acceleration.", 
            generic_unit_system.time, 
            1.0 | generic_unit_system.time
        ) 
        
        object.add_boolean_parameter(
            "get_nogravity",
            None,
            "no_gravity_flag",
            "No-gravity flag. True means: gravitational forces are switched of for all particles "
                "(read-only: makefile option NOGRAVITY).",
            False
        )
        
        object.add_boolean_parameter(
            "get_gdgop",
            "set_gdgop",
            "gadget_cell_opening_flag",
            "Gadget-cell-opening flag. True means: use of Gadget cell opening criterion; Barnes-Hut otherwise",
            True
        )
        
        object.add_boolean_parameter(
            "get_isotherm",
            None,
            "isothermal_flag",
            "Isothermal flag. True means: isothermal gas, u is interpreted as c_s^2 "
                "(read-only: makefile option ISOTHERM_EQS).",
            False
        )
        
        object.add_boolean_parameter(
            "get_eps_is_h",
            None,
            "eps_is_h_flag",
            "Eps-is-h flag. True means: set gas particles gravitational epsilon to h (SPH smoothing length) "
                "(read-only: makefile option ADAPTIVE_GRAVSOFT_FORGAS).",
            False
        )
        
        object.add_method_parameter(
            "get_nsmooth", 
            "set_nsmooth",
            "nsmooth", 
            "The target number of SPH neighbours.", 
            units.none,
            50 | units.none
        )
        
        object.add_method_parameter(
            "get_unit_mass", 
            None,
            "code_mass_unit", 
            "The code mass unit (in g/h, 1.989e43 g = 10^10 MSun standard).", 
            units.g,
            1.989e43 | units.g
        )
        
        object.add_method_parameter(
            "get_unit_length", 
            None,
            "code_length_unit", 
            "The code length unit (in cm/h, 3.085678e21 cm = 1 kpc standard).", 
            units.cm,
            3.085678e21 | units.cm
        )
        
        object.add_method_parameter(
            "get_unit_time", 
            None,
            "code_time_unit", 
            "The code time unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr).", 
            units.s,
            3.085678e16 | units.s
        )
        
        object.add_method_parameter(
            "get_unit_velocity", 
            None,
            "code_velocity_unit", 
            "The code velocity unit (in cm/s, default: 1e5 cm/s = 1 km/s).", 
            units.cm / units.s,
            1e5 | units.cm / units.s
        )
        
        
        object.add_method_parameter(
            "get_bh_tol", 
            "set_bh_tol",
            "opening_angle", 
            "Opening angle, theta, for building the tree: between 0 and 1 (unitless, 0.5).", 
            units.none,
            0.5 | units.none
        )
        
        object.add_method_parameter(
            "get_gdgtol", 
            "set_gdgtol",
            "gadget_cell_opening_constant", 
            "Gadget-cell-openings criterion parameter  (unitless, 0.005)", 
            units.none,
            0.005 | units.none
        )
        
        object.add_method_parameter(
            "get_epsgas", 
            "set_epsgas",
            "gas_epsilon", 
            "The gas gravitational smoothing epsilon.", 
            generic_unit_system.length,
            0.01 | generic_unit_system.length
        )
        
        object.add_method_parameter(
            "get_gamma", 
            None,
            "polytropic_index_gamma", 
            "gas polytropic index (1.6666667 or 1 for isothermal"
                "(read-only: makefile option ISOTHERM_EQS).", 
            units.none,
            (5.0/3) | units.none
        )
        
        object.add_method_parameter(
            "get_alpha", 
            "set_alpha",
            "artificial_viscosity_alpha", 
            "SPH artificial viscosity alpha parameter (0.5)", 
            units.none,
            0.5 | units.none
        )
        
        object.add_method_parameter(
            "get_courant", 
            "set_courant",
            "courant", 
            "SPH courant condition parameter (0.3). Note that we follow conventional smoothing length "
                "definitions, implying a factor 2 difference with Gadget's CourantFac parameter", 
            units.none,
            0.3 | units.none
        )
        
        object.add_method_parameter(
            "get_nsmtol", 
            "set_nsmtol",
            "n_neighbour_tol", 
            "fractional tolerance in number of SPH neighbours", 
            units.none,
            0.1 | units.none
        )
        
        object.add_method_parameter(
            "get_gadget_output_directory", 
            "set_gadget_output_directory",
            "gadget_output_directory", 
            "Name of the Gadget-2 OutputDir", 
            units.string,
            "" | units.string
        )
        
        self.stopping_conditions.define_parameters(object)        
        
    
    def define_particle_sets(self, object):
        object.define_super_set('particles', ['dm_particles','gas_particles'], 
            index_to_default_set = 0)
        
        object.define_set('dm_particles', 'index_of_the_particle')
        object.set_new('dm_particles', 'new_dm_particle')
        object.set_delete('dm_particles', 'delete_particle')
        object.add_setter('dm_particles', 'set_state')
        object.add_getter('dm_particles', 'get_state')
        object.add_setter('dm_particles', 'set_mass')
        object.add_getter('dm_particles', 'get_mass', names = ('mass',))
        object.add_setter('dm_particles', 'set_position')
        object.add_getter('dm_particles', 'get_position')
        object.add_setter('dm_particles', 'set_velocity')
        object.add_getter('dm_particles', 'get_velocity')
        object.add_getter('dm_particles', 'get_acceleration')
        
        object.define_set('gas_particles', 'index_of_the_particle')
        object.set_new('gas_particles', 'new_sph_particle')
        object.set_delete('gas_particles', 'delete_particle')
        object.add_setter('gas_particles', 'set_state_sph')
        object.add_getter('gas_particles', 'get_state_sph')
        object.add_setter('gas_particles', 'set_mass')
        object.add_getter('gas_particles', 'get_mass', names = ('mass',))
        object.add_setter('gas_particles', 'set_position')
        object.add_getter('gas_particles', 'get_position')
        object.add_setter('gas_particles', 'set_velocity')
        object.add_getter('gas_particles', 'get_velocity')
        object.add_getter('gas_particles', 'get_acceleration')
        object.add_setter('gas_particles', 'set_internal_energy')
        object.add_getter('gas_particles', 'get_internal_energy')

        #what are we going to do with this?...
        self.stopping_conditions.define_particle_set(object, 'dm_particles')
        self.stopping_conditions.define_particle_set(object, 'gas_particles')
    
    def define_methods(self, object):
        #GravitationalDynamics.define_methods(self, object)
        object.add_method(
            'evolve',
            (generic_unit_system.time,),
            public_name = 'evolve_model'
        )
        object.add_method(
            "new_particle",
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "delete_particle",
            (
                object.NO_UNIT,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_state",
            (
                object.NO_UNIT,
            ),
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_state",
            (
                object.NO_UNIT,
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_mass",
            (
                object.NO_UNIT,
                generic_unit_system.mass,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_mass",
            (
                object.NO_UNIT,
            ),
            (
                generic_unit_system.mass,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_radius",
            (
                object.NO_UNIT,
                generic_unit_system.length,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_radius",
            (
                object.NO_UNIT,
            ),
            (
                generic_unit_system.length,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_position",
            (
                object.NO_UNIT,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_position",
            (
                object.NO_UNIT,
            ),
            (
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_velocity",
            (
                object.INDEX,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_velocity",
            (
                object.INDEX,
            ),
            (
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "new_dm_particle",
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "new_sph_particle",
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.specific_energy,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_state_sph",
            (
                object.INDEX,
            ),
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.specific_energy,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_state_sph",
            (
                object.INDEX,
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.specific_energy,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "set_internal_energy",
            (
                object.INDEX,
                generic_unit_system.specific_energy,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_internal_energy",
            (
                object.INDEX,
            ),
            (
                generic_unit_system.specific_energy,
                object.ERROR_CODE
            )
        )
        
        object.add_method(
            'get_gravity_at_point',
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length),
            (generic_unit_system.acceleration, generic_unit_system.acceleration, generic_unit_system.acceleration, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_potential_at_point',
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length),
            (generic_unit_system.potential, object.ERROR_CODE)
        )
        
        self.stopping_conditions.define_methods(object)

