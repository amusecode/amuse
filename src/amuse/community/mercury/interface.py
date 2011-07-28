from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

class MercuryInterface(CodeInterface, LiteratureReferencesMixIn, StoppingConditionInterface):
    """
    Mercury N-body integrator package, version 6.2.
    Mercury is a general-purpose N-body integration package for problems in 
    celestial mechanics. The standard symplectic (MVS) algorithm is described in 
    Widsom & Holman (1991). The hybrid symplectic algorithm is described 
    in Chambers (1999).
    
    Relevant references:
        .. [#] Chambers J. E., 1999, MNRAS, 304, 793
        .. [#] Widsom J. & Holman M., 1991, AJ, 102, 1528
    """
    def __init__(self, **args):
        CodeInterface.__init__(self, name_of_the_worker = 'mercury_worker',**args)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function   
    def cleanup_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def commit_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def commit_parameters():
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
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_initial_timestep():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT)
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
        for x in ['mass','density','x','y','z','vx','vy','vz','Lx','Ly','Lz','celimit']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_central_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','radius','j2','j4','j6','Lx','Ly','Lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_orbiter_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','density','x','y','z','vx','vy','vz','Lx','Ly','Lz','celimit']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_orbiter_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','density','x','y','z','vx','vy','vz','Lx','Ly','Lz','celimit']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_particle_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','j2','j4','j6','Lx','Ly','Lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_particle_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','j2','j4','j6','Lx','Ly','Lz']:
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
        for x in ['Lx','Ly','Lz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_spin():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['Lx','Ly','Lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_angularmomentum():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['Lx','Ly','Lz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_angularmomentum():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['Lx','Ly','Lz']:
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


    def define_parameters(self, object):
        object.add_method_parameter(
            "get_time",
            None,
            "time",
            "current simulation time", 
            default_value = 0.0 | units.s
        )

    def define_parameters(self, object):
        object.add_method_parameter(
            "set_initial_timestep",
            "get_initial_timestep",
            "timestep",
            "current simulation time", 
            default_value = 0.0 | units.s
        )

        self.stopping_conditions.define_parameters(object)

    def define_properties(self, object):
        object.add_property("get_kinetic_energy")
        object.add_property("get_potential_energy")
        object.add_property("get_total_energy")
        object.add_property("get_center_of_mass_position")
        object.add_property("get_center_of_mass_velocity")
        object.add_property("get_total_mass")

    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_method('INITIALIZED', 'invoke_state_change')
        object.add_transition_to_method('END', 'cleanup_code')
        object.add_method('END', 'stop')

#    def define_state(self, object):
#        GravitationalDynamics.define_state(self, object)
#        object.add_method('EDIT', 'new_central_particle')
#        object.add_transition('RUN', 'UPDATE', 'new_central_particle', False)
#        object.add_method('EDIT', 'new_orbiter')
#        object.add_transition('RUN', 'UPDATE', 'new_orbiter', False)

    def define_particle_sets(self, object):
        object.define_super_set('particles', ['central_particle','orbiters'],
            index_to_default_set=0)

        object.define_set('orbiters', 'id')
        object.set_new('orbiters', 'new_orbiter')
        object.set_delete('orbiters', 'delete_particle')
        object.add_setter('orbiters', 'set_orbiter_state') 
        object.add_getter('orbiters', 'get_orbiter_state') 
        object.add_setter('orbiters', 'set_mass')
        object.add_getter('orbiters', 'get_mass')
        object.add_setter('orbiters', 'set_density')
        object.add_getter('orbiters', 'get_density')
        object.add_setter('orbiters', 'set_position')
        object.add_getter('orbiters', 'get_position')
        object.add_setter('orbiters', 'set_velocity')
        object.add_getter('orbiters', 'get_velocity')
        object.add_setter('orbiters', 'set_angularmomentum')
        object.add_getter('orbiters', 'get_angularmomentum')
        object.add_setter('orbiters', 'set_celimit')
        object.add_getter('orbiters', 'get_celimit')

        object.define_set('central_particle', 'id')
        object.set_new('central_particle', 'new_central_particle')
        object.add_setter('central_particle', 'set_central_particle_state')
        object.add_getter('central_particle', 'get_central_particle_state')
        object.add_setter('central_particle', 'set_central_mass')
        object.add_getter('central_particle', 'get_central_mass')
        object.add_setter('central_particle', 'set_central_radius')
        object.add_getter('central_particle', 'get_central_radius')
        object.add_setter('central_particle', 'set_central_oblateness')
        object.add_getter('central_particle', 'get_central_oblateness')
        object.add_setter('central_particle', 'set_central_spin')
        object.add_getter('central_particle', 'get_central_spin')

        #GravitationalDynamics.define_particle_sets(self, object)

    def define_methods(self, object):
        #GravitationalDynamics.define_methods(self, object)
        object.add_method('evolve_model', (units.day), ( object.ERROR_CODE, ))
        object.add_method(
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
                units.none
            ),
            (
                object.INDEX, 
                object.ERROR_CODE
            )
        )
        object.add_method(
            'get_orbiter_state',
            (
                object.INDEX,
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
                units.none,
                object.ERROR_CODE
            )
        )
        object.add_method(
            'get_central_particle_state',
            (
                object.INDEX,
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
                object.ERROR_CODE
            )
        )

        object.add_method(
            'set_orbiter_state',
            (
                object.INDEX,
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
                units.none
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
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
                object.INDEX, 
                object.ERROR_CODE
            )
        )

        object.add_method(
            'set_central_particle_state',
            (
                object.INDEX,
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
                object.ERROR_CODE
            )
        )

        object.add_method(
            "set_mass",
            (
                object.INDEX,
                units.MSun,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_mass",
            (
                object.INDEX,
            ),
            (
                units.MSun,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_central_mass",
            (
                object.INDEX,
                units.MSun,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_central_mass",
            (
                object.INDEX,
            ),
            (
                units.MSun,
                object.ERROR_CODE
            )
        )

        #assuming celimit is RCEH, clouse-encounter limit
        #expressed in units Hill radius, I use unit none
        #see comments in: src/mercury_main.for
        object.add_method(
            "set_celimit",
            (
                object.INDEX,
                units.none,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_celimit",
            (
                object.INDEX,
            ),
            (
                units.none,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "set_position",
            (
                object.INDEX,
                units.AU,
                units.AU,
                units.AU,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_position",
            (
                object.INDEX,
            ),
            (
                units.AU,
                units.AU,
                units.AU,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_velocity",
            (
                object.INDEX,
                units.AUd,
                units.AUd,
                units.AUd,
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
                units.AUd,
                units.AUd,
                units.AUd,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_angularmomentum",
            (
                object.INDEX,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_angularmomentum",
            (
                object.INDEX,
            ),
            (
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_central_spin",
            (
                object.INDEX,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_central_spin",
            (
                object.INDEX,
            ),
            (
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                units.MSun * units.AU**2/units.day,
                object.ERROR_CODE
            )
        )

        object.add_method(
            'set_density',
            (
                object.INDEX,
                units.g/units.cm**3
            ),
            (
                object.ERROR_CODE
            )
        )    
        object.add_method(
            'get_density',
            (
                object.INDEX,
            ),
            (
                units.g/units.cm**3,
                object.ERROR_CODE
            )
        )    
        object.add_method(
            'get_central_oblateness',
            (
                object.INDEX,
            ),
            (
                units.AU**2,
                units.AU**4,
                units.AU**6,
                object.ERROR_CODE
            )
        )
        object.add_method(
            'set_central_oblateness',
            (
                object.INDEX,
                units.AU**2,
                units.AU**4,
                units.AU**6,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            'get_central_radius',
            (
                object.INDEX,
            ),
            (
                units.AU,
                object.ERROR_CODE
            )
        )
        object.add_method(
            'set_central_radius',
            (
                object.INDEX,
                units.AU,
            ),
            (
                object.ERROR_CODE,
            )
        )
        
        object.add_method(
            "get_time",
            (),
            (units.s, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_kinetic_energy",
            (),
            (units.MSun*units.AU**2/units.day**2, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_potential_energy",
            (),
            (units.MSun*units.AU**2/units.day**2, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_total_energy",
            (),
            (units.MSun*units.AU**2/units.day**2, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_center_of_mass_position",
            (),
            (units.AU, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_center_of_mass_velocity",
            (),
            (units.AUd, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_total_mass",
            (),
            (units.MSun, object.ERROR_CODE,)
        )
        
        self.stopping_conditions.define_methods(object)
 
