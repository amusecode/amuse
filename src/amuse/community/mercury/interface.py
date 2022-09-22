import os
import numpy

from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

from amuse.datamodel import Particles

class MercuryInterface(CodeInterface, CommonCodeInterface, CodeWithDataDirectories,
                       LiteratureReferencesMixIn, StoppingConditionInterface):
    """
    Mercury N-body integrator package, version 6.2.
    Mercury is a general-purpose N-body integration package for problems in 
    celestial mechanics. The standard symplectic (MVS) algorithm is described in 
    Wisdom & Holman (1991). The hybrid symplectic algorithm is described 
    in Chambers (1999).
    
    Relevant references:
        .. [#] Chambers J. E., 1999, MNRAS, 304, 793 [1999MNRAS.304..793C]
        .. [#] Wisdom J. & Holman M., 1991, AJ, 102, 1528 [1991AJ....102.1528W]
    """

    use_modules = ['StoppingConditions', 'AmuseInterface']

    def __init__(self, **args):
        CodeInterface.__init__(self, name_of_the_worker = 'mercury_worker',**args)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
                
    @legacy_function    
    def commit_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function


    @legacy_function
    def recommit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function    
    def evolve_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def synchronize_model():
        """
        After an evolve the particles may be at different simulation
        times. Synchronize the particles to a consistent stat
        at the current simulation time
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK

        """
        return function

    @legacy_function    
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_initial_timestep():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN, unit = units.day)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_initial_timestep():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT, unit = units.day)
        function.result_type = 'i'
        return function

    @legacy_function          
    def delete_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function    

    @legacy_function    
    def new_orbiter():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT)
        for x in ['mass','density','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        for x,default in zip(['lx','ly','lz','celimit'], [0,0,0,3]):
            function.addParameter(x, dtype='float64', direction=function.IN, default = default)
            
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_central_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        
        for x in ['mass']:
            function.addParameter(x, dtype='d', direction=function.IN)
            
        for x in ['radius','j2','j4','j6','lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.IN, default = 0)
            
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_orbiter_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','density','x','y','z','vx','vy','vz','lx','ly','lz','celimit']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_orbiter_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','density','x','y','z','vx','vy','vz','lx','ly','lz','celimit']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_particle_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','j2','j4','j6','lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_particle_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','j2','j4','j6','lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the position from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current x position of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current y position of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the position for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new x position of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new y position of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the position from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current x position of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current y position of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the position for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new x position of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new y position of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new z position of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function    
    def get_density():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('density', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_density():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('density', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_radius():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_radius():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_celimit():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('celimit', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_celimit():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('celimit', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_oblateness():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['j2','j4','j6']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_oblateness():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['j2','j4','j6']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_spin():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_spin():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_angularmomentum():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_angularmomentum():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_number_of_orbiters():
        function = LegacyFunctionSpecification()   
        function.addParameter('norbiters', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_energy_deviation():
        function = LegacyFunctionSpecification()   
        function.addParameter('delta_e', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('ek', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_potential_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('ep', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_total_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('e_tot', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_total_angular_momentum():
        function = LegacyFunctionSpecification()   
        function.addParameter('am_tot', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_begin_time():
        """
        Retrieve the model time to start the evolution at.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT,
            description = "The begin time", unit = units.day)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time was retrieved
        -2 - ERROR
            The code does not have support for querying the begin time
        """
        return function
    
    @legacy_function
    def set_begin_time():
        """
        Set the model time to start the evolution at. This is an offset for
        all further calculations in the code.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "The model time to start at", unit = units.day)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Time value was changed
        -2 - ERROR
            The code does not support setting the begin time
        """
        return function

    @legacy_function
    def get_integrator():
        """
        Retrieve the integrator (algor parameter)
        """
        function = LegacyFunctionSpecification()
        function.addParameter('integrator', dtype='i', direction=function.OUT,
            description = "integrator to use")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_integrator():
        """
        Set the model integrator (algor parameter)
        """
        function = LegacyFunctionSpecification()
        function.addParameter('integrator', dtype='i', direction=function.IN,
            description = "integrator to use")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            integrator was changed
        -1 - ERROR
            integrator not supported (wrong input)
        """
        return function
        
    @legacy_function
    def get_rmax():
        """
        Retrieve the maximal radius (rmax) -- heliocentric distance at 
        which objects are considered ejected.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('rmax', dtype='float64', direction=function.OUT,
            description = "heliocentric distance at which objects are considered ejected", unit=units.AU)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_rmax():
        """
        Set the maximal radius (rmax) -- heliocentric distance at 
        which objects are considered ejected.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('rmax', dtype='float64', direction=function.IN,
            description = "heliocentric distance at which objects are considered ejected", unit=units.AU)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_cefac():
        """
        Retrieve the parameter to set the changeover distance RCRIT of the hybrid integrator: 
        RCRIT = cefac * R_HILL (n_1 in section 4.2 of Chambers 1999)
        """
        function = LegacyFunctionSpecification()
        function.addParameter('cefac', dtype='float64', direction=function.OUT,
            description = "Hybrid integrator changeover radius RCRIT (in Hill radii)", unit=units.none)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_cefac():
        """
        Set the parameter to set the changeover distance RCRIT of the hybrid integrator: 
        RCRIT = cefac * R_HILL (n_1 in section 4.2 of Chambers 1999)
        """
        function = LegacyFunctionSpecification()
        function.addParameter('cefac', dtype='float64', direction=function.IN,
            description = "Hybrid integrator changeover radius RCRIT (in Hill radii)", unit=units.none)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_elements_file():
        function = LegacyFunctionSpecification()
        function.addParameter('elements_file', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    @legacy_function
    def set_elements_file():
        function = LegacyFunctionSpecification()
        function.addParameter('elements_file', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_close_encounters_file():
        function = LegacyFunctionSpecification()
        function.addParameter('close_encounters_file', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    @legacy_function
    def set_close_encounters_file():
        function = LegacyFunctionSpecification()
        function.addParameter('close_encounters_file', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_info_file():
        function = LegacyFunctionSpecification()
        function.addParameter('info_file', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    @legacy_function
    def set_info_file():
        function = LegacyFunctionSpecification()
        function.addParameter('info_file', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_bigbody_file():
        function = LegacyFunctionSpecification()
        function.addParameter('bigbody_file', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    @legacy_function
    def set_bigbody_file():
        function = LegacyFunctionSpecification()
        function.addParameter('bigbody_file', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_smallbody_file():
        function = LegacyFunctionSpecification()
        function.addParameter('smallbody_file', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    @legacy_function
    def set_smallbody_file():
        function = LegacyFunctionSpecification()
        function.addParameter('smallbody_file', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_integration_parameters_file():
        function = LegacyFunctionSpecification()
        function.addParameter('integration_parameters_file', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    @legacy_function
    def set_integration_parameters_file():
        function = LegacyFunctionSpecification()
        function.addParameter('integration_parameters_file', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_restart_file():
        function = LegacyFunctionSpecification()
        function.addParameter('restart_file', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    @legacy_function
    def set_restart_file():
        function = LegacyFunctionSpecification()
        function.addParameter('restart_file', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gravity_at_point():
        """
        Get the gravitational acceleration at the given points. To calculate the force on
        bodies at those points, multiply with the mass of the bodies
        """
        function = LegacyFunctionSpecification()
        for x in ['eps','x','y','z']:
            function.addParameter(
              x,
              dtype='float64',
              direction=function.IN,
              unit=units.AU
            )
        for x in ['ax','ay','az']:
            function.addParameter(
                x,
                dtype='float64',
                direction=function.OUT,
                unit=units.AU/units.day**2
            )
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_potential_at_point():
        """
        Determine the gravitational potential on any given point
        """
        function = LegacyFunctionSpecification()
        for x in ['eps','x','y','z']:
            function.addParameter(
                x,
                dtype='float64',
                direction=function.IN,
                unit=units.AU
            )
        for x in ['phi']:
            function.addParameter(
                x,
                dtype='float64',
                direction=function.OUT,
                unit=units.AU**2/units.day**2
            )
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

class MercuryDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_doc+"\n\n"+instance.parameters.__doc__

class MercuryWayWard(GravitationalDynamics):

    __doc__ = MercuryDoc()

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = MercuryInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def initialize_code(self):
        self.overridden().initialize_code()
        self.set_elements_file(os.devnull)
        self.set_close_encounters_file(os.devnull)
        self.set_info_file(os.devnull)
        self.set_bigbody_file(os.devnull)
        self.set_smallbody_file(os.devnull)
        self.set_integration_parameters_file(os.devnull)
        self.set_restart_file(os.devnull)

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_time",
            None,
            "time",
            "current simulation time", 
            default_value = 0.0 | units.day
        )
        
        handler.add_method_parameter(
            "get_initial_timestep",
            "set_initial_timestep",
            "timestep",
            "initial timestep", 
            default_value = 8.0 | units.day
        )
        
        handler.add_method_parameter(
            "get_integrator",
            "set_integrator",
            "integrator",
            "integrator to use", 
            default_value = 10
        )
        
        handler.add_method_parameter(
            "get_rmax",
            "set_rmax",
            "rmax",
            "heliocentric distance at for ejection", 
            default_value = 100. | units.AU
        )
        
        handler.add_method_parameter(
            "get_cefac",
            "set_cefac",
            "cefac",
            "Hybrid integrator changeover radius RCRIT (in Hill radii)", 
            default_value = 3. | units.none
        )

        handler.add_method_parameter(
            "get_elements_file",
            "set_elements_file",
            "elements_file",
            "outputfile for mercury data (orbital elements) [/dev/null]", 
            default_value = None
        )        

        handler.add_method_parameter(
            "get_close_encounters_file",
            "set_close_encounters_file",
            "close_encounters_file",
            "outputfile for mercury data (close encounters) [/dev/null]", 
            default_value = None
        )
        
        handler.add_method_parameter(
            "get_info_file",
            "set_info_file",
            "info_file",
            "outputfile for mercury data (info) [/dev/null]", 
            default_value = None
        )        

        handler.add_method_parameter(
            "get_bigbody_file",
            "set_bigbody_file",
            "bigbody_file",
            "dumpfile for mercury data (big bodies) [/dev/null]", 
            default_value = None
        )

        handler.add_method_parameter(
            "get_smallbody_file",
            "set_smallbody_file",
            "smallbody_file",
            "dumpfile for mercury data (small bodies) [/dev/null]", 
            default_value = None
        )

        handler.add_method_parameter(
            "get_integration_parameters_file",
            "set_integration_parameters_file",
            "integration_parameters_file",
            "dumpfile for mercury data (integration parameters) [/dev/null]", 
            default_value = None
        )

        handler.add_method_parameter(
            "get_restart_file",
            "set_restart_file",
            "restart_file",
            "dumpfile for mercury data (restart data - only mercury internal state) [/dev/null]", 
            default_value = None
        )
        
        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | units.day
        )

        self.stopping_conditions.define_parameters(handler)

    
    def define_properties(self, handler):
        handler.add_property("get_kinetic_energy")
        handler.add_property("get_potential_energy")
        handler.add_property("get_total_energy")
        handler.add_property("get_center_of_mass_position")
        handler.add_property("get_center_of_mass_velocity")
        handler.add_property("get_total_mass")

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        handler.add_method('EDIT', 'new_central_particle')
        handler.add_method('EDIT', 'new_orbiter')
        handler.add_method('UPDATE', 'new_orbiter')
        handler.add_transition('RUN', 'UPDATE', 'new_central_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'new_orbiter', False)
        self.stopping_conditions.define_state(handler)



    def define_particle_sets(self, handler):
        handler.define_super_set('particles', ['central_particle','orbiters'],
            index_to_default_set=0)

        handler.define_set('orbiters', 'id')
        handler.set_new('orbiters', 'new_orbiter')
        handler.set_delete('orbiters', 'delete_particle')
        handler.add_setter('orbiters', 'set_orbiter_state') 
        handler.add_getter('orbiters', 'get_orbiter_state') 
        handler.add_setter('orbiters', 'set_mass')
        handler.add_getter('orbiters', 'get_mass')
        handler.add_setter('orbiters', 'set_density')
        handler.add_getter('orbiters', 'get_density')
        handler.add_setter('orbiters', 'set_position')
        handler.add_getter('orbiters', 'get_position')
        handler.add_setter('orbiters', 'set_velocity')
        handler.add_getter('orbiters', 'get_velocity')
        handler.add_setter('orbiters', 'set_angularmomentum')
        handler.add_getter('orbiters', 'get_angularmomentum')
        handler.add_setter('orbiters', 'set_celimit')
        handler.add_getter('orbiters', 'get_celimit')

        handler.define_set('central_particle', 'id')
        handler.set_new('central_particle', 'new_central_particle')
        handler.add_setter('central_particle', 'set_central_particle_state')
        handler.add_getter('central_particle', 'get_central_particle_state')
        handler.add_setter('central_particle', 'set_central_mass')
        handler.add_getter('central_particle', 'get_central_mass')
        handler.add_setter('central_particle', 'set_central_radius')
        handler.add_getter('central_particle', 'get_central_radius')
        handler.add_setter('central_particle', 'set_central_oblateness')
        handler.add_getter('central_particle', 'get_central_oblateness')
        handler.add_setter('central_particle', 'set_central_spin')
        handler.add_getter('central_particle', 'get_central_spin')

        #GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)

    def define_methods(self, handler):
        #GravitationalDynamics.define_methods(self, handler)
        handler.add_method('evolve_model', (units.day), ( handler.ERROR_CODE, ))
        handler.add_method(
            'new_orbiter',
            (
                units.MSun,
                units.g/units.cm**3,
                units.AU,
                units.AU,
                units.AU,
                units.AUd,
                units.AUd,
                units.AUd,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                handler.NO_UNIT
            ),
            (
                handler.INDEX, 
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            'get_orbiter_state',
            (
                handler.INDEX,
            ),
            (
                units.MSun,
                units.g/units.cm**3,
                units.AU,
                units.AU,
                units.AU,
                units.AUd,
                units.AUd,
                units.AUd,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                handler.NO_UNIT,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            'get_central_particle_state',
            (
                handler.INDEX,
            ),
            (
                units.MSun,
                units.AU,
                units.AU**2,
                units.AU**4,
                units.AU**6,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            'set_orbiter_state',
            (
                handler.INDEX,
                units.MSun,
                units.g/units.cm**3,
                units.AU,
                units.AU,
                units.AU,
                units.AUd,
                units.AUd,
                units.AUd,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                handler.NO_UNIT
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            'new_central_particle',
            (
                units.MSun,
                units.AU,
                units.AU**2,
                units.AU**4,
                units.AU**6,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day
            ),
            (
                handler.INDEX, 
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            'set_central_particle_state',
            (
                handler.INDEX,
                units.MSun,
                units.AU,
                units.AU**2,
                units.AU**4,
                units.AU**6,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day
            ),
            (
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "set_initial_timestep",
            (
                units.day,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_initial_timestep",
            (
            ),
            (
                units.day,
                handler.ERROR_CODE
            )
        )



        handler.add_method(
            "set_mass",
            (
                handler.INDEX,
                units.MSun,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_mass",
            (
                handler.INDEX,
            ),
            (
                units.MSun,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_central_mass",
            (
                handler.INDEX,
                units.MSun,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_central_mass",
            (
                handler.INDEX,
            ),
            (
                units.MSun,
                handler.ERROR_CODE
            )
        )

        #assuming celimit is RCEH, clouse-encounter limit
        #expressed in units Hill radius, I use unit none
        #see comments in: src/mercury_main.for
        handler.add_method(
            "set_celimit",
            (
                handler.INDEX,
                handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_celimit",
            (
                handler.INDEX,
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "set_position",
            (
                handler.INDEX,
                units.AU,
                units.AU,
                units.AU,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_position",
            (
                handler.INDEX,
            ),
            (
                units.AU,
                units.AU,
                units.AU,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_velocity",
            (
                handler.INDEX,
                units.AUd,
                units.AUd,
                units.AUd,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_velocity",
            (
                handler.INDEX,
            ),
            (
                units.AUd,
                units.AUd,
                units.AUd,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_angularmomentum",
            (
                handler.INDEX,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_angularmomentum",
            (
                handler.INDEX,
            ),
            (
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_central_spin",
            (
                handler.INDEX,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_central_spin",
            (
                handler.INDEX,
            ),
            (
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            'set_density',
            (
                handler.INDEX,
                units.g/units.cm**3
            ),
            (
                handler.ERROR_CODE
            )
        )    
        handler.add_method(
            'get_density',
            (
                handler.INDEX,
            ),
            (
                units.g/units.cm**3,
                handler.ERROR_CODE
            )
        )    
        handler.add_method(
            'get_central_oblateness',
            (
                handler.INDEX,
            ),
            (
                units.AU**2,
                units.AU**4,
                units.AU**6,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            'set_central_oblateness',
            (
                handler.INDEX,
                units.AU**2,
                units.AU**4,
                units.AU**6,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            'get_central_radius',
            (
                handler.INDEX,
            ),
            (
                units.AU,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            'set_central_radius',
            (
                handler.INDEX,
                units.AU,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        
        handler.add_method(
            "get_time",
            (),
            (units.day, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_kinetic_energy",
            (),
            (units.MSun*units.AU**2/units.day**2, handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_total_angular_momentum",
            (),
            (units.MSun*units.AU**2/units.day, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_potential_energy",
            (),
            (units.MSun*units.AU**2/units.day**2, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_total_energy",
            (),
            (units.MSun*units.AU**2/units.day**2, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_center_of_mass_position",
            (),
            (units.AU, units.AU, units.AU, handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_center_of_mass_velocity",
            (),
            (units.AUd, units.AUd, units.AUd, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_total_mass",
            (),
            (units.MSun, handler.ERROR_CODE,)
        )
        
        self.stopping_conditions.define_methods(handler)
 
class Mercury(MercuryWayWard):
    def __init__(self, *args, **kargs):
        MercuryWayWard.__init__(self, *args, **kargs)
        self._particles=Particles(0)
        self.model_time=0.|units.s
        self.particles_accessed=True
        self.committed=False

    def define_parameters(self, handler):
        MercuryWayWard.define_parameters(self,handler)

        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | units.AU * units.AU
        )

    def get_eps2(self):
        return 0|units.AU * units.AU

    def set_eps2(self, eps2):
       if eps2.number!=0:
            raise Exception("Mercury cannot use non-zero smoothing")

    def get_number_of_particles(self):
        return len(self.particles)

    def get_total_radius(self):
        max_dist2 = 0|units.AU * units.AU
        for p in self.particles:
            dist2 = p.x*p.x + p.y*p.y + p.z*p.z
            if dist2 > max_dist2:
                max_dist2 = dist2
        return ((max_dist2.number**0.5)|units.AU)

    @property
    def particles(self):
        if not self.particles_accessed:
            self.particles_accessed=True
        return self._particles

    @property
    def total_mass(self):
        return self.particles.total_mass()


    def get_center_of_mass_position(self):
        return self.particles.center_of_mass()

    def get_center_of_mass_velocity(self):
        return self.particles.center_of_mass_velocity()
        
    def commit_particles(self):
        N=len(self.particles)
        if N<=1:
            print("too few particles")
            return -11

        ic=self.particles.mass.argmax()
        self.central_particle=self.particles[ic]

        orbiters=(self.particles-self.central_particle).copy()

        maxmass=orbiters.mass.amax()
        if (maxmass/self.central_particle.mass) > 0.1:
            print("orbiters too massive")
            return -12

        orbiters.position=orbiters.position-self.central_particle.position
        orbiters.velocity=orbiters.velocity-self.central_particle.velocity

# note: lx,ly,lz, celimit?
        if not hasattr(orbiters,'density'):
            orbiters.density=orbiters.mass*3/(4*numpy.pi)/orbiters.radius**3

        self.overridden().central_particle.add_particle(self.central_particle)
                
        self.overridden().orbiters.add_particles(orbiters)
        
        self.overridden().commit_particles()

# commit_particles may already change the particle positions (I think it doesn't though)
        """
        com_position=self.particles.center_of_mass()
        com_velocity=self.particles.center_of_mass_velocity()

        self.particles.position=self.particles.position-self.central_particle.position
        self.particles.velocity=self.particles.velocity-self.central_particle.velocity

        channel=self.overridden().orbiters.new_channel_to(self.particles)

        channel.copy_attributes(["x","y","z","vx","vy","vz","mass"])
        
        self.particles.move_to_center()
        self.particles.position+=com_position
        self.particles.velocity+=com_velocity
        """
        self.model_time=self.overridden().get_time()
        self.committed=True
        self.particles_accessed=False
        return 0
        
        
    def recommit_particles(self):
        if not self.particles_accessed:
          return 0
        if self.overridden().central_particle[0] not in self.particles:
          raise Exception("you are not allowed to remove the central particle")
# these are allowed now
        if self.overridden().central_particle[0].mass != self.central_particle.mass:
          self.overridden().central_particle.mass=self.central_particle.mass
        if self.overridden().central_particle.radius != self.central_particle.radius:
          self.overridden().central_particle.radius=self.central_particle.radius

        orbiters=(self.particles-self.central_particle).copy()

        maxmass=orbiters.mass.amax()
        if (maxmass/self.central_particle.mass) > 0.1:
          print("orbiters too massive")
          return -12

        orbiters.position=orbiters.position-self.central_particle.position
        orbiters.velocity=orbiters.velocity-self.central_particle.velocity
        if not hasattr(orbiters,'density'):
          orbiters.density=orbiters.mass*3/(4*numpy.pi)/orbiters.radius**3
# note: lx,ly,lz, celimit?
   
        orbiters.synchronize_to(self.overridden().orbiters)                

        channel=orbiters.new_channel_to(self.overridden().orbiters)

        channel.copy_attributes(["x","y","z","vx","vy","vz","mass","density"])

        self.overridden().recommit_particles()

# recommit_particles could already change the particle positions (I think it doesn't though)
        """ 
        com_position=self.particles.center_of_mass()
        com_velocity=self.particles.center_of_mass_velocity()

        channel=self.overridden().orbiters.new_channel_to(self.particles)
        channel.copy_attributes(["x","y","z","vx","vy","vz","mass"])
        self.central_particle.position*=0
        self.central_particle.velocity*=0
        self.particles.move_to_center()
        self.particles.position+=com_position
        self.particles.velocity+=com_velocity
        """
        self.model_time=self.overridden().get_time()
        self.particles_accessed=False
        return 0
        
    def evolve_model(self, tend):
        if self.particles_accessed:
          if not self.committed:
            ret=self.commit_particles()
          else:  
            ret=self.recommit_particles()

        com_position=self.particles.center_of_mass()
        com_velocity=self.particles.center_of_mass_velocity()
        
        ret=self.overridden().evolve_model(tend)
        tmodel=self.overridden().get_time()
        com_position+=com_velocity*(tmodel-self.model_time)
        self.model_time=tmodel            
        channel=self.overridden().orbiters.new_channel_to(self.particles)
        channel.copy_attributes(["x","y","z","vx","vy","vz","mass"])
        self.central_particle.position*=0
        self.central_particle.velocity*=0
        self.particles.move_to_center()
        self.particles.position+=com_position
        self.particles.velocity+=com_velocity
        self.particles_accessed=False
        return ret

    def get_potential_at_point(self,eps,x,y,z):
        if self.particles_accessed:
          if not self.committed:
            ret=self.commit_particles()
          else:  
            ret=self.recommit_particles()
            
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z

        return self.overridden().get_potential_at_point(eps,xx,yy,zz)

    def get_gravity_at_point(self,eps,x,y,z):
        if self.particles_accessed:
          if not self.committed:
            ret=self.commit_particles()
          else:  
            ret=self.recommit_particles()
            
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z

        return self.overridden().get_gravity_at_point(eps,xx,yy,zz)

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        handler.add_method('EDIT', 'new_central_particle')
        handler.add_method('EDIT', 'new_orbiter')
        handler.add_method('UPDATE', 'new_orbiter')
        handler.add_transition('RUN', 'UPDATE', 'new_central_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'new_orbiter', False)
        handler.add_method('RUN', 'recommit_particles')

    def cleanup_code(self):
        self._particles=Particles(0)
        self.model_time=0.|units.s
        self.particles_accessed=True
        self.committed=False

