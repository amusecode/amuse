import os
import numpy
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics
from amuse.community import *
from amuse.support.options import option

class Gadget2Interface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, StoppingConditionInterface):
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
        .. [#] Durier F., Dalla Vecchia C., 2011, MNRAS (Time integration scheme fix)
    """
    include_headers = ['gadget_code.h', 'worker_code.h', 'stopcond.h']
    
    MODE_NORMAL = 'normal'
    MODE_PERIODIC_BOUNDARIES   = 'periodic'
    MODE_NOGRAVITY = 'nogravity'
    
    def __init__(self, mode = MODE_NORMAL,  **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)
        
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'gadget2_worker'
        elif mode == self.MODE_PERIODIC_BOUNDARIES:
            return 'gadget2_worker_periodic'
        elif mode == self.MODE_NOGRAVITY:
            return 'gadget2_worker_nogravity'
        else:
            return 'gadget2_worker'
    
            
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
        
    @option(type="string")
    def data_directory(self):
        """
        Returns the root name of the directory for the Gadget2
        application data files.
        """
        return os.path.join(self.input_data_root_directory, 'gadget2', 'input')
        
    
    @option(type="string")
    def output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, 'gadget2', 'output')        
    
    def get_data_directory(self):
        """
        Returns the root name of the directory for the Gadget2
        application data files.
        """
        return self.data_directory
    
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return self.output_directory
    
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
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_potential():
        """
        Retrieve the current potential of a particle. 
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the potential from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('Potential', dtype='float64', direction=function.OUT, description = "The current potential of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
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
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
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
    def get_mass():
        """
        Retrieve the mass of a particle. Mass is a scalar property of a particle,
        this function has one OUT argument.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
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
        -2 - ERROR
            code does not support updating of a particle
        """
        return function
    @legacy_function
    def get_position():
        """
        Retrieve the position vector of a particle. Position is a vector property,
        this function has 3 OUT arguments.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
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
    def get_velocity():
        """
        Retrieve the velocity vector of a particle. Position is a vector property,
        this function has 3 OUT arguments.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the velocity from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current x component of the position vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current y component of the position vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current z component of the position vector of the particle")
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
    def get_acceleration():
        """
        Retrieve the acceleration vector of a particle. Second time derivative of the position.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('ax', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('ay', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('az', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
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
    def get_state_sph():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_state_sph():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('u', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_internal_energy():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('u', dtype='float64', direction=function.IN)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_smoothing_length():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('h_smooth', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('rho', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_pressure():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('pressure', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_d_internal_energy_dt():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('du_dt', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_n_neighbours():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('num_neighbours', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_thermal_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('thermal_energy', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_epsilon_dm_part():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('radius', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_epsilon_gas_part():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('radius', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
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
        return function
        
    @legacy_function
    def set_gadget_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('gadget_output_directory', dtype='string', direction=function.IN,
            description = "The path to the Gadget-2 OutputDir.")
        function.result_type = 'int32'
        return function
    
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
        return function
        
    @legacy_function
    def set_gdgop():
        """ set_gdgop([0,1]): use of gadget cell opening criterion if 1 """
        function = LegacyFunctionSpecification()
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_gdgop():
        """ get_gdgop(): use of gadget cell opening criterion if 1 """
        function = LegacyFunctionSpecification()
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_epsgas():
        """ gas grav smoothing eps"""
        function = LegacyFunctionSpecification()
        function.addParameter('gas_epsilon', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_epsgas():
        """ gas grav smoothing eps"""
        function = LegacyFunctionSpecification()
        function.addParameter('gas_epsilon', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_isotherm():
        """(default: False)
        True means: isothermal gas (read-only: makefile option ISOTHERM_EQS).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('isothermal_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_eps_is_h():
        """ get_eps_is_h(): gas particles grav. eps to SPH h if 1"""
        function = LegacyFunctionSpecification()
        function.addParameter('eps_is_h_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()
        function.addParameter('nsmooth', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()
        function.addParameter('nsmooth', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """
        function = LegacyFunctionSpecification()
        function.addParameter('bh_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """
        function = LegacyFunctionSpecification()
        function.addParameter('bh_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_gdgtol():
        """ Gadget cell openings criterion parameter  (unitless, .005) """
        function = LegacyFunctionSpecification()
        function.addParameter('gdgtol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_gdgtol():
        """ Gadget cell openings criterion parameter  (unitless, .005) """
        function = LegacyFunctionSpecification()
        function.addParameter('gdgtol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_gamma():
        """ gas polytropic index (1.666667) """        
        function = LegacyFunctionSpecification()
        function.addParameter('gamma', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()
        function.addParameter('alpha', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()
        function.addParameter('alpha', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()
        function.addParameter('courant', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()
        function.addParameter('courant', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_nsmtol():
        """ fractional tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()
        function.addParameter('n_neighbour_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_nsmtol():
        """ fractional tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()
        function.addParameter('n_neighbour_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_energy_file():
        function = LegacyFunctionSpecification()
        function.addParameter('energy_file', dtype='string', direction=function.OUT,
            description = "The path to the Gadget-2 energy statistics output file.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_energy_file():
        function = LegacyFunctionSpecification()
        function.addParameter('energy_file', dtype='string', direction=function.IN,
            description = "The path to the Gadget-2 energy statistics output file.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_info_file():
        function = LegacyFunctionSpecification()
        function.addParameter('info_file', dtype='string', direction=function.OUT,
            description = "The path to the Gadget-2 info output file.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_info_file():
        function = LegacyFunctionSpecification()
        function.addParameter('info_file', dtype='string', direction=function.IN,
            description = "The path to the Gadget-2 info output file.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_timings_file():
        function = LegacyFunctionSpecification()
        function.addParameter('timings_file', dtype='string', direction=function.OUT,
            description = "The path to the Gadget-2 timings output file.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_timings_file():
        function = LegacyFunctionSpecification()
        function.addParameter('timings_file', dtype='string', direction=function.IN,
            description = "The path to the Gadget-2 timings output file.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_cpu_file():
        function = LegacyFunctionSpecification()
        function.addParameter('cpu_file', dtype='string', direction=function.OUT,
            description = "The path to the Gadget-2 cpu statistics output file.")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_cpu_file():
        function = LegacyFunctionSpecification()
        function.addParameter('cpu_file', dtype='string', direction=function.IN,
            description = "The path to the Gadget-2 cpu statistics output file.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_time_limit_cpu():
        function = LegacyFunctionSpecification()
        function.addParameter('time_limit_cpu', dtype='d', direction=function.OUT,
            description = "The cpu-time limit. Gadget2 will stop once 85% of this (wall-clock) time has passed.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_time_limit_cpu():
        function = LegacyFunctionSpecification()
        function.addParameter('time_limit_cpu', dtype='d', direction=function.IN,
            description = "The cpu-time limit. Gadget2 will stop once 85% of this (wall-clock) time has passed.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_comoving_integration_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('comoving_integration_flag', dtype='bool', direction=function.OUT,
            description = "Flag to do a cosmological run with comoving coordinates.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_comoving_integration_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('comoving_integration_flag', dtype='bool', direction=function.IN,
            description = "Flag to do a cosmological run with comoving coordinates.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_type_of_timestep_criterion():
        function = LegacyFunctionSpecification()
        function.addParameter('type_of_timestep_criterion', dtype='i', direction=function.OUT,
            description = "Timestep criterion to use. Can only be zero: timestep proportional to acceleration^-0.5")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_type_of_timestep_criterion():
        function = LegacyFunctionSpecification()
        function.addParameter('type_of_timestep_criterion', dtype='i', direction=function.IN,
            description = "Timestep criterion to use. Can only be zero: timestep proportional to acceleration^-0.5")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_time_begin():
        function = LegacyFunctionSpecification()
        function.addParameter('time_begin', dtype='d', direction=function.OUT,
            description = "The time at the start of the run.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_time_begin():
        function = LegacyFunctionSpecification()
        function.addParameter('time_begin', dtype='d', direction=function.IN,
            description = "The time at the start of the run.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_time_max():
        function = LegacyFunctionSpecification()
        function.addParameter('time_max', dtype='d', direction=function.OUT,
            description = "The time at the end of the run.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_time_max():
        function = LegacyFunctionSpecification()
        function.addParameter('time_max', dtype='d', direction=function.IN,
            description = "The time at the end of the run.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_omega_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_zero', dtype='d', direction=function.OUT,
            description = "Cosmological matter density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_omega_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_zero', dtype='d', direction=function.IN,
            description = "Cosmological matter density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_omega_lambda():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_lambda', dtype='d', direction=function.OUT,
            description = "Cosmological vacuum energy density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_omega_lambda():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_lambda', dtype='d', direction=function.IN,
            description = "Cosmological vacuum energy density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_omega_baryon():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_baryon', dtype='d', direction=function.OUT,
            description = "Cosmological baryonic density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_omega_baryon():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_baryon', dtype='d', direction=function.IN,
            description = "Cosmological baryonic density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_hubble_param():
        function = LegacyFunctionSpecification()
        function.addParameter('hubble_param', dtype='d', direction=function.OUT,
            description = "The cosmological Hubble parameter.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_hubble_param():
        function = LegacyFunctionSpecification()
        function.addParameter('hubble_param', dtype='d', direction=function.IN,
            description = "The cosmological Hubble parameter.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_err_tol_int_accuracy():
        function = LegacyFunctionSpecification()
        function.addParameter('err_tol_int_accuracy', dtype='d', direction=function.OUT,
            description = "Accuracy parameter used in timestep criterion. Actual timesteps are proportional to err_tol_int_accuracy^0.5")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_err_tol_int_accuracy():
        function = LegacyFunctionSpecification()
        function.addParameter('err_tol_int_accuracy', dtype='d', direction=function.IN,
            description = "Accuracy parameter used in timestep criterion. Actual timesteps are proportional to err_tol_int_accuracy^0.5")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_max_size_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('max_size_timestep', dtype='d', direction=function.OUT,
            description = "The maximum size of the timestep a particle may take.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_max_size_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('max_size_timestep', dtype='d', direction=function.IN,
            description = "The maximum size of the timestep a particle may take.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_min_size_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('min_size_timestep', dtype='d', direction=function.OUT,
            description = "The minimum size of the timestep a particle may take.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_min_size_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('min_size_timestep', dtype='d', direction=function.IN,
            description = "The minimum size of the timestep a particle may take.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_tree_domain_update_frequency():
        function = LegacyFunctionSpecification()
        function.addParameter('tree_domain_update_frequency', dtype='d', direction=function.OUT,
            description = "The frequency with which the tree and domain decomposition are fully updated, in terms of (# force computations / # particles).")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_tree_domain_update_frequency():
        function = LegacyFunctionSpecification()
        function.addParameter('tree_domain_update_frequency', dtype='d', direction=function.IN,
            description = "The frequency with which the tree and domain decomposition are fully updated, in terms of (# force computations / # particles).")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_time_between_statistics():
        function = LegacyFunctionSpecification()
        function.addParameter('time_between_statistics', dtype='d', direction=function.OUT,
            description = "The time between statistics output written to the output files.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_time_between_statistics():
        function = LegacyFunctionSpecification()
        function.addParameter('time_between_statistics', dtype='d', direction=function.IN,
            description = "The time between statistics output written to the output files.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_min_gas_temp():
        function = LegacyFunctionSpecification()
        function.addParameter('min_gas_temp', dtype='d', direction=function.OUT,
            description = "The minimum temperature of gas particles.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_min_gas_temp():
        function = LegacyFunctionSpecification()
        function.addParameter('min_gas_temp', dtype='d', direction=function.IN,
            description = "The minimum temperature of gas particles.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_min_gas_hsmooth_fractional():
        function = LegacyFunctionSpecification()
        function.addParameter('min_gas_hsmooth_fractional', dtype='d', direction=function.OUT,
            description = "The minimum smoothing length of gas particles relative to their softening lengths.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_min_gas_hsmooth_fractional():
        function = LegacyFunctionSpecification()
        function.addParameter('min_gas_hsmooth_fractional', dtype='d', direction=function.IN,
            description = "The minimum smoothing length of gas particles relative to their softening lengths.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_softening_gas_max_phys():
        function = LegacyFunctionSpecification()
        function.addParameter('softening_gas_max_phys', dtype='d', direction=function.OUT,
            description = "The maximum physical softening of gas particles for comoving integrations.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_softening_gas_max_phys():
        function = LegacyFunctionSpecification()
        function.addParameter('softening_gas_max_phys', dtype='d', direction=function.IN,
            description = "The maximum physical softening of gas particles for comoving integrations.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_softening_halo_max_phys():
        function = LegacyFunctionSpecification()
        function.addParameter('softening_halo_max_phys', dtype='d', direction=function.OUT,
            description = "The maximum physical softening of dm particles for comoving integrations.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_softening_halo_max_phys():
        function = LegacyFunctionSpecification()
        function.addParameter('softening_halo_max_phys', dtype='d', direction=function.IN,
            description = "The maximum physical softening of dm particles for comoving integrations.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_periodic_boundaries_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_periodic_boundaries_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_interpret_kicks_as_feedback_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_interpret_kicks_as_feedback_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_interpret_heat_as_feedback_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_interpret_heat_as_feedback_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_box_size():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_box_size():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_hydro_state_at_point():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','rhoe']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i' 
        return function
    

class Gadget2Doc(object):

    def __get__(self, instance, owner):
        return instance.legacy_interface.__doc__+"\n\n"+instance.parameters.__doc__

class Gadget2(GravitationalDynamics):
    
    __doc__ = Gadget2Doc()
    
    def __init__(self, unit_converter = None, mode = 'normal', **options):
        self.mode = mode
        legacy_interface = Gadget2Interface(mode = mode, **options)
        if unit_converter is None:
            unit_converter = ConvertBetweenGenericAndSiUnits(
                3.085678e21 | units.cm,   # 1.0 kpc
                1.989e43 | units.g,       # 1.0e10 solar masses
                1e5 | units.cm / units.s) # 1 km/sec
    
        self.stopping_conditions = StoppingConditions(self)
    
        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            unit_converter,
            **options
        )
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        if self.mode == self.legacy_interface.MODE_PERIODIC_BOUNDARIES:
            self.parameters.periodic_boundaries_flag = True
            
        ensure_data_directory_exists(self.get_output_directory())
        
        self.parameters.gadget_output_directory = self.get_output_directory()
        # The code's units are read-only, and set here to ensure they always match with the unit_converter
        self.set_unit_mass(self.unit_converter.to_si(generic_unit_system.mass).value_in(units.g))
        self.set_unit_length(self.unit_converter.to_si(generic_unit_system.length).value_in(units.cm))
        self.set_unit_time(self.unit_converter.to_si(generic_unit_system.time).value_in(units.s))
        return result
    
    def define_properties(self, object):
        object.add_property("get_kinetic_energy")
        object.add_property("get_potential_energy")
        object.add_property("get_thermal_energy")
        object.add_property("get_total_radius")
        object.add_property("get_center_of_mass_position")
        object.add_property("get_center_of_mass_velocity")
        object.add_property("get_total_mass")
        object.add_property('get_time', public_name = "model_time")
    
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        object.add_method('EDIT', 'new_dm_particle')
        object.add_method('UPDATE', 'new_dm_particle')
        object.add_transition('RUN', 'UPDATE', 'new_dm_particle', False)
        object.add_method('EDIT', 'new_sph_particle')
        object.add_method('UPDATE', 'new_sph_particle')
        object.add_transition('RUN', 'UPDATE', 'new_sph_particle', False)
        object.add_method('RUN', 'get_state_sph')
        object.add_method('RUN', 'get_acceleration')
        object.add_method('RUN', 'get_internal_energy')
        object.add_method('RUN', 'get_smoothing_length')
        object.add_method('RUN', 'get_density')
        object.add_method('RUN', 'get_pressure')
        object.add_method('RUN', 'get_d_internal_energy_dt')
        object.add_method('RUN', 'get_n_neighbours')
        object.add_method('RUN', 'get_epsilon_dm_part')
        object.add_method('RUN', 'get_epsilon_gas_part')
        object.add_method('RUN', 'set_state')
        object.add_method('RUN', 'set_state_sph')
        object.add_method('RUN', 'set_mass')
        object.add_method('RUN', 'set_position')
        object.add_method('RUN', 'set_velocity')
        object.add_method('RUN', 'set_internal_energy')
        
        object.add_method('RUN', 'get_kinetic_energy')
        object.add_method('RUN', 'get_potential_energy')
        object.add_method('RUN', 'get_thermal_energy')
        object.add_method('RUN', 'get_total_radius')
        object.add_method('RUN', 'get_center_of_mass_position')
        object.add_method('RUN', 'get_center_of_mass_velocity')
        object.add_method('RUN', 'get_total_mass')
        object.add_method('RUN', 'get_hydro_state_at_point')
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_epsilon_squared", 
            "set_epsilon_squared",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0001 | generic_unit_system.length * generic_unit_system.length
        )
        
        object.add_method_parameter(
            "get_time_step", 
            None,
            "timestep", 
            "timestep for system, Gadget2 calculates this by itself, based on particle acceleration.", 
            default_value = 1.0 | generic_unit_system.time
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
            "n_smooth", 
            "The target number of SPH neighbours.", 
            default_value = 50 | units.none
        )
        
        object.add_method_parameter(
            "get_unit_mass", 
            None,
            "code_mass_unit", 
            "The code mass unit (in g/h, 1.989e43 g = 10^10 MSun standard).", 
            default_value = 1.989e43 | units.g
        )
        
        object.add_method_parameter(
            "get_unit_length", 
            None,
            "code_length_unit", 
            "The code length unit (in cm/h, 3.085678e21 cm = 1 kpc standard).", 
            default_value = 3.085678e21 | units.cm
        )
        
        object.add_method_parameter(
            "get_unit_time", 
            None,
            "code_time_unit", 
            "The code time unit (in s/h, default: 3.085678e16 s = (1 kpc) / (1 km/s) ~ 0.9778 Gyr).", 
            default_value = 3.085678e16 | units.s
        )
        
        object.add_method_parameter(
            "get_unit_velocity", 
            None,
            "code_velocity_unit", 
            "The code velocity unit (in cm/s, default: 1e5 cm/s = 1 km/s).", 
            default_value = 1e5 | units.cm / units.s
        )
        
        
        object.add_method_parameter(
            "get_bh_tol", 
            "set_bh_tol",
            "opening_angle", 
            "Opening angle, theta, for building the tree: between 0 and 1 (unitless, 0.5).", 
            default_value = 0.5 | units.none
        )
        
        object.add_method_parameter(
            "get_gdgtol", 
            "set_gdgtol",
            "gadget_cell_opening_constant", 
            "Gadget-cell-openings criterion parameter  (unitless, 0.005)", 
            default_value = 0.005 | units.none
        )
        
        object.add_method_parameter(
            "get_epsgas", 
            "set_epsgas",
            "gas_epsilon", 
            "The gas gravitational smoothing epsilon.", 
            default_value = 0.01 | generic_unit_system.length
        )
        
        object.add_method_parameter(
            "get_gamma", 
            None,
            "polytropic_index_gamma", 
            "gas polytropic index (1.6666667 or 1 for isothermal"
                "(read-only: makefile option ISOTHERM_EQS).", 
            default_value = (5.0/3) | units.none
        )
        
        object.add_method_parameter(
            "get_alpha", 
            "set_alpha",
            "artificial_viscosity_alpha", 
            "SPH artificial viscosity alpha parameter (0.5)", 
            default_value = 0.5 | units.none
        )
        
        object.add_method_parameter(
            "get_courant", 
            "set_courant",
            "courant", 
            "SPH courant condition parameter (0.3). Note that we follow conventional smoothing length "
                "definitions, implying a factor 2 difference with Gadget's CourantFac parameter", 
            default_value = 0.3 | units.none
        )
        
        object.add_method_parameter(
            "get_nsmtol", 
            "set_nsmtol",
            "n_smooth_tol", 
            "fractional tolerance in number of SPH neighbours", 
            default_value = 0.1 | units.none
        )
        
        object.add_method_parameter(
            "get_gadget_output_directory", 
            "set_gadget_output_directory",
            "gadget_output_directory", 
            "Name of the Gadget-2 OutputDir", 
            default_value = ""
        )
        
        object.add_method_parameter(
            "get_energy_file", 
            "set_energy_file",
            "energy_file", 
            "The path to the Gadget-2 energy statistics output file.", 
            default_value = "energy.txt"
        )
        
        object.add_method_parameter(
            "get_info_file", 
            "set_info_file",
            "info_file", 
            "The path to the Gadget-2 info output file.", 
            default_value = "info.txt"
        )
        
        object.add_method_parameter(
            "get_timings_file", 
            "set_timings_file",
            "timings_file", 
            "The path to the Gadget-2 timings output file.", 
            default_value = "timings.txt"
        )
        
        object.add_method_parameter(
            "get_cpu_file", 
            "set_cpu_file",
            "cpu_file", 
            "The path to the Gadget-2 cpu statistics output file.", 
            default_value = "cpu.txt"
        )
        
        object.add_method_parameter(
            "get_time_limit_cpu", 
            "set_time_limit_cpu",
            "time_limit_cpu", 
            "The cpu-time limit. Gadget2 will stop once 85% of this (wall-clock) time has passed.", 
            default_value = 36000 | units.s
        )
        
        object.add_boolean_parameter(
            "get_comoving_integration_flag",
            "set_comoving_integration_flag",
            "comoving_integration_flag",
            "Flag to do a cosmological run with comoving coordinates.",
            False
        )
        
        object.add_method_parameter(
            "get_type_of_timestep_criterion", 
            "set_type_of_timestep_criterion",
            "type_of_timestep_criterion", 
            "Timestep criterion to use. Can only be zero: timestep proportional to acceleration^-0.5", 
            default_value = 0 | units.none
        )
        
        object.add_method_parameter(
            "get_time_begin", 
            "set_time_begin",
            "time_begin", 
            "The time at the start of the run.", 
            default_value = 0.0 | generic_unit_system.time
        )
        
        object.add_method_parameter(
            "get_time_max", 
            "set_time_max",
            "time_max", 
            "The time at the end of the run.", 
            default_value = 100.0 | generic_unit_system.time
        )
        
        object.add_method_parameter(
            "get_omega_zero", 
            "set_omega_zero",
            "omega_zero", 
            "Cosmological matter density parameter in units of the critical density at z=0.", 
            default_value = 0.0 | units.none
        )
        
        object.add_method_parameter(
            "get_omega_lambda", 
            "set_omega_lambda",
            "omega_lambda", 
            "Cosmological vacuum energy density parameter in units of the critical density at z=0.", 
            default_value = 0.0 | units.none
        )
        
        object.add_method_parameter(
            "get_omega_baryon", 
            "set_omega_baryon",
            "omega_baryon", 
            "Cosmological baryonic density parameter in units of the critical density at z=0.", 
            default_value = 0.0 | units.none
        )
        
        object.add_method_parameter(
            "get_hubble_param", 
            "set_hubble_param",
            "hubble_param", 
            "The cosmological Hubble parameter.", 
            default_value = 0.7 | 100 * units.km / units.s / units.Mpc
        )
        
        object.add_method_parameter(
            "get_err_tol_int_accuracy", 
            "set_err_tol_int_accuracy",
            "timestep_accuracy_parameter", 
            "Accuracy parameter used in timestep criterion. Actual timesteps are proportional to err_tol_int_accuracy^0.5", 
            default_value = 0.025 | units.none
        )
        
        object.add_method_parameter(
            "get_max_size_timestep", 
            "set_max_size_timestep",
            "max_size_timestep", 
            "The maximum size of the timestep a particle may take.", 
            default_value = 0.01 | generic_unit_system.time
        )
        
        object.add_method_parameter(
            "get_min_size_timestep", 
            "set_min_size_timestep",
            "min_size_timestep", 
            "The minimum size of the timestep a particle may take.", 
            default_value = 0.0 | generic_unit_system.time
        )
        
        object.add_method_parameter(
            "get_tree_domain_update_frequency", 
            "set_tree_domain_update_frequency",
            "tree_domain_update_frequency", 
            "The frequency with which the tree and domain decomposition are fully updated, in terms of (# force computations / # particles).", 
            default_value = 0.05 | units.none
        )
        
        object.add_method_parameter(
            "get_time_between_statistics", 
            "set_time_between_statistics",
            "time_between_statistics", 
            "The time between statistics output written to the output files.", 
            default_value = 0.1 | generic_unit_system.time
        )
        
        object.add_method_parameter(
            "get_min_gas_temp", 
            "set_min_gas_temp",
            "min_gas_temp", 
            "The minimum temperature of gas particles.", 
            default_value = 0.0 | units.K
        )
        
        object.add_method_parameter(
            "get_min_gas_hsmooth_fractional", 
            "set_min_gas_hsmooth_fractional",
            "min_gas_hsmooth_fractional", 
            "The minimum smoothing length of gas particles relative to their softening lengths.", 
            default_value = 0.0 | units.none
        )
        
        object.add_method_parameter(
            "get_softening_gas_max_phys", 
            "set_softening_gas_max_phys",
            "softening_gas_max_phys", 
            "The maximum physical softening of gas particles for comoving integrations.", 
            default_value = 0.0 | generic_unit_system.length
        )
        
        object.add_method_parameter(
            "get_softening_halo_max_phys", 
            "set_softening_halo_max_phys",
            "softening_halo_max_phys", 
            "The maximum physical softening of dm particles for comoving integrations.", 
            default_value = 0.0 | generic_unit_system.length
        )
        
        object.add_method_parameter(
            "get_box_size", 
            "set_box_size",
            "periodic_box_size", 
            "The size of the box in case of periodic boundary conditions.", 
            default_value = 1.0 | generic_unit_system.length
        )
        
        object.add_boolean_parameter(
            "get_periodic_boundaries_flag",
            "set_periodic_boundaries_flag",
            "periodic_boundaries_flag",
            "Periodic boundaries flag. True means: use periodic boundary conditions",
            False
        )
        
        object.add_boolean_parameter(
            "get_interpret_kicks_as_feedback_flag",
            "set_interpret_kicks_as_feedback_flag",
            "interpret_kicks_as_feedback",
            "Flag telling Gadget2 whether to interpret external changes to particles' velocities as feedback (for timestepping).",
            False
        )
        
        object.add_boolean_parameter(
            "get_interpret_heat_as_feedback_flag",
            "set_interpret_heat_as_feedback_flag",
            "interpret_heat_as_feedback",
            "Flag telling Gadget2 whether to interpret external changes to particles' internal energy as feedback (for timestepping).",
            True
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
        object.add_getter('dm_particles', 'get_epsilon_dm_part', names = ('radius',))
        object.add_getter('dm_particles', 'get_epsilon_dm_part', names = ('epsilon',))
        
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
        object.add_getter('gas_particles', 'get_smoothing_length')
        object.add_getter('gas_particles', 'get_density', names = ('rho',))
        object.add_getter('gas_particles', 'get_density', names = ('density',))
        object.add_getter('gas_particles', 'get_pressure')
        object.add_getter('gas_particles', 'get_d_internal_energy_dt')
        object.add_getter('gas_particles', 'get_n_neighbours')
        object.add_getter('gas_particles', 'get_epsilon_gas_part', names = ('radius',))
        object.add_getter('gas_particles', 'get_epsilon_gas_part', names = ('epsilon',))

        self.stopping_conditions.define_particle_set(object)
    
    def define_errorcodes(self, object):
        object.add_errorcode(-1, 'Unspecified, other error.')
        object.add_errorcode(-2, 'Called function is not implemented.')
        object.add_errorcode(-3, 'A particle with the given index was not found.')
        object.add_errorcode(-4, 'Parameter check failed.')
        object.add_errorcode(-5, 'CPU-time limit reached.')
        object.add_errorcode(-6, "Can't evolve backwards in time.")
        object.add_errorcode(-7, "Can't evolve further than time_max.")
        object.add_errorcode(-8, "A particle was assigned a timestep of size zero. The code_time_unit used may be too large.")
    
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        object.add_method('evolve_model', (generic_unit_system.time,), ( object.ERROR_CODE, ))
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
            "get_acceleration",
            (
                object.INDEX,
            ),
            (
                generic_unit_system.acceleration,
                generic_unit_system.acceleration,
                generic_unit_system.acceleration,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_potential",
            (
                object.NO_UNIT,
            ),
            (
                generic_unit_system.potential,
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
            "get_smoothing_length",
            (object.INDEX,),
            (generic_unit_system.length, object.ERROR_CODE)
        )
        object.add_method(
            "get_density",
            (object.INDEX,),
            (generic_unit_system.density, object.ERROR_CODE)
        )
        object.add_method(
            "get_pressure",
            (object.INDEX,),
            (generic_unit_system.pressure, object.ERROR_CODE)
        )
        object.add_method(
            "get_d_internal_energy_dt",
            (object.INDEX,),
            (generic_unit_system.specific_energy / generic_unit_system.time, object.ERROR_CODE)
        )
        object.add_method(
            "get_n_neighbours",
            (object.INDEX,),
            (units.none, object.ERROR_CODE)
        )
        object.add_method(
            "get_epsilon_dm_part",
            (object.INDEX,),
            (generic_unit_system.length, object.ERROR_CODE)
        )
        object.add_method(
            "get_epsilon_gas_part",
            (object.INDEX,),
            (generic_unit_system.length, object.ERROR_CODE)
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
        
        object.add_method(
            'get_hydro_state_at_point',
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,
                generic_unit_system.speed, generic_unit_system.speed, generic_unit_system.speed),
            (generic_unit_system.density, generic_unit_system.momentum_density, generic_unit_system.momentum_density, 
                generic_unit_system.momentum_density, generic_unit_system.energy_density, object.ERROR_CODE)
        )
        
        object.add_method(
            "get_epsilon_squared",
            (),
            (generic_unit_system.length * generic_unit_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_epsilon_squared",
            (generic_unit_system.length * generic_unit_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time_step",
            (),
            (generic_unit_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_nsmooth",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_nsmooth",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_unit_mass",
            (),
            (units.g, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_unit_length",
            (),
            (units.cm, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_unit_time",
            (),
            (units.s, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_unit_velocity",
            (),
            (units.cm / units.s, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_bh_tol",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_bh_tol",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_gdgtol",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_gdgtol",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_epsgas",
            (),
            (generic_unit_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_epsgas",
            (generic_unit_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_gamma",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_alpha",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_alpha",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_courant",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_courant",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_nsmtol",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_nsmtol",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_gadget_output_directory",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_gadget_output_directory",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_energy_file",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_energy_file",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_info_file",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_info_file",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_timings_file",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_timings_file",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_cpu_file",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_cpu_file",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time_limit_cpu",
            (),
            (units.s, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time_limit_cpu",
            (units.s, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_type_of_timestep_criterion",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_type_of_timestep_criterion",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time_begin",
            (),
            (generic_unit_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time_begin",
            (generic_unit_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time_max",
            (),
            (generic_unit_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time_max",
            (generic_unit_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_omega_zero",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_omega_zero",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_omega_lambda",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_omega_lambda",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_omega_baryon",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_omega_baryon",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_hubble_param",
            (),
            (100 * units.km / units.s / units.Mpc, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_hubble_param",
            (100 * units.km / units.s / units.Mpc, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_err_tol_int_accuracy",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_err_tol_int_accuracy",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_max_size_timestep",
            (),
            (generic_unit_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_max_size_timestep",
            (generic_unit_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_min_size_timestep",
            (),
            (generic_unit_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_min_size_timestep",
            (generic_unit_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tree_domain_update_frequency",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_tree_domain_update_frequency",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time_between_statistics",
            (),
            (generic_unit_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time_between_statistics",
            (generic_unit_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_min_gas_temp",
            (),
            (units.K, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_min_gas_temp",
            (units.K, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_min_gas_hsmooth_fractional",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_min_gas_hsmooth_fractional",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_softening_gas_max_phys",
            (),
            (generic_unit_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_softening_gas_max_phys",
            (generic_unit_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_softening_halo_max_phys",
            (),
            (generic_unit_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_softening_halo_max_phys",
            (generic_unit_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_box_size",
            (),
            (generic_unit_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_box_size",
            (generic_unit_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_kinetic_energy",
            (),
            (generic_unit_system.energy, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_potential_energy",
            (),
            (generic_unit_system.energy, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_thermal_energy",
            (),
            (generic_unit_system.energy, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_total_radius",
            (),
            (generic_unit_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_center_of_mass_position",
            (),
            (generic_unit_system.length,generic_unit_system.length,generic_unit_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_center_of_mass_velocity",
            (),
            (generic_unit_system.speed,generic_unit_system.speed,generic_unit_system.speed, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_total_mass",
            (),
            (generic_unit_system.mass, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_time',
            (),
            (generic_unit_system.time, object.ERROR_CODE,)
        )
        
        self.stopping_conditions.define_methods(object)

