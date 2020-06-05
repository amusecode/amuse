from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class ReboundInterface(CodeInterface,
                       LiteratureReferencesMixIn,
                       GravitationalDynamicsInterface,
                       StoppingConditionInterface,
                       #SinglePointGravityFieldInterface
    ):
    """
    REBOUND - An open-source multi-purpose N-body code
    
    .. [#] Rein, H., Liu, S.F., *Astronomy and Astrophysics* , **Volume 537**, A128 (2012)
    
    For different integrators, cite:
    ... IAS15:  Rein, H., Spiegel, D.S., *MNRAS* , **Volume 446**, Issue 2, p.1424-1437 (2015)
    ... WHFast: Rein, H., Tamayo, D., *MNRAS* , **Volume 452**, Issue 1, p.376-388 (2015)
    ... Hermes: Silburt, A., et al., in prep.
    ... SEI:    Rein, H., Tremaine, S., *MNRAS* , **Volume 415**, Issue 4, p.3168-3176 (2011)
    ... JANUS:  Rein, H., Tamayo, D., *arXiv* , 1704.07715 (2017)
        
    """
    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="rebound_worker",
                                 **options)
        LiteratureReferencesMixIn.__init__(self)
        
    @legacy_function
    def new_particle():
        """
        Define a new particle in the stellar dynamics code. The particle is initialized with the provided
        mass, radius, position and velocity. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )

        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle", default = 0)
        function.addParameter('subset', dtype='int32', direction=function.IN, description = "The subset index of the particle (defaults to 0, use new_subset for higher indices)", default = 0)
        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            particle was created and added to the model
        -1 - ERROR
            particle could not be created"""
        return function
        
    def delete_particle(self, index_of_the_particle, code_index=0):
        return self._delete_particle(index_of_the_particle, code_index)

    @legacy_function
    def _delete_particle():
        """
        Delete a particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN, description ="Index of the particle")
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            particle was deleted
        -1 - ERROR
            particle not deleted"""
        return function

    @legacy_function
    def _set_integrator():
        function = LegacyFunctionSpecification()      
        function.addParameter('integrator_name', dtype='i', direction=function.IN)
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
        

    @legacy_function
    def _get_integrator():
        function = LegacyFunctionSpecification()      
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('integrator_name', dtype='i', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
    

    INTEGRATORS = {"ias15": 0, "whfast": 1, "sei": 2, "leapfrog": 4, "hermes": 5, "whfast-helio": 6, "none": 7, "janus": 8}
    def set_integrator(self, name, code_index = 0 ):
        return self._set_integrator(self.INTEGRATORS[name], code_index)
    
    def get_integrator(self, code_index = 0):
        value, error = self._get_integrator(code_index)
        for key, index in self.INTEGRATORS.items():
            if value == index:
                return key
        return "none"

    @legacy_function
    def _set_solver():
        function = LegacyFunctionSpecification()      
        function.addParameter('solver_name', dtype='i', direction=function.IN)
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
        

    @legacy_function
    def _get_solver():
        function = LegacyFunctionSpecification()      
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('solver_name', dtype='i', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
    

    SOLVERS = {"none": 0, "basic": 1, "compensated": 2, "tree": 3}
    def set_solver(self, name, code_index = 0 ):
        return self._set_solver(self.SOLVERS[name], code_index)
    
    def get_solver(self, code_index = 0):
        value, error = self._get_solver(code_index)
        for key, index in self.SOLVERS.items():
            if value == index:
                return key
        return "none"

    @legacy_function
    def get_opening_angle2():
        function = LegacyFunctionSpecification()
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('opening_angle2', dtype='float64', direction=function.OUT,
                description = "theta, the opening angle for building the tree: between 0 and 1")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_opening_angle2():
        function = LegacyFunctionSpecification()
        function.addParameter('opening_angle2', dtype='float64', direction=function.IN,
                description = "theta, the opening angle for building the tree: between 0 and 1")
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eps2():
        function = LegacyFunctionSpecification()
        """
        Get epsilon^2, a softening parameter for gravitational potentials with point particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('epsilon_squared', dtype='float64', direction=function.OUT,
            description = "epsilon^2, a softening parameter for gravitational potentials with point particles",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_eps2():
        """
        Set epsilon^2, a softening parameter for gravitational potentials with point particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.IN,
            description = "epsilon^2, a softening parameter for gravitational potentials with point particles",
            unit = nbody_system.length * nbody_system.length)
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def _set_boundary():
        function = LegacyFunctionSpecification()      
        function.addParameter('boundary_name', dtype='i', direction=function.IN)
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
        

    @legacy_function
    def _get_boundary():
        function = LegacyFunctionSpecification()      
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('boundary_name', dtype='i', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
    

    BOUNDARIES = {"none": 0, "open": 1, "periodic": 2, "shear": 3}
    def set_boundary(self, name, code_index = 0 ):
        return self._set_boundary(self.BOUNDARIES[name], code_index)
    
    def get_boundary(self, code_index = 0):
        value, error = self._get_boundary(code_index)
        for key, index in self.BOUNDARIES.items():
            if value == index:
                return key
        return "none"

    @legacy_function
    def get_boundary_size():
        function = LegacyFunctionSpecification()
        """
        Get the size of the boundaries.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('boundary_size', dtype='float64', direction=function.OUT,
            description = "boundary size",
            unit = nbody_system.length)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_boundary_size():
        """
        Set size of the boundaries.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('boundary_size', dtype='float64', direction=function.IN,
            description = "boundary size",
            unit = nbody_system.length)
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function




    @legacy_function
    def set_time_step():
        """
        Update timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64', direction=function.IN,
            description = "timestep")
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            timestep was changed
        """
        return function


    @legacy_function
    def get_potential_energy():
        """
        Retrieve the current potential energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('potential_energy', dtype='float64', direction=function.OUT,
            description = "The potential energy of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the potential energy was set
        -1 - ERROR
            Kinetic potential could not be provided
        """
        return function
        

    @legacy_function
    def get_kinetic_energy():
        """
        Retrieve the current kinetic energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('kinetic_energy', dtype='float64', direction=function.OUT,
            description = "The kinetic energy of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the kinetic energy was set
        -1 - ERROR
            Kinetic energy could not be provided
        """
        return function


    @legacy_function
    def evolve_model():
        """
        Evolve the model until the given time, or until a stopping condition is set.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "Model time to evolve the code to. The model will be "
                "evolved until this time is reached exactly or just after.")
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound (default -1, evolve all systems)", default = -1)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_time():
        """
        Retrieve the model time. This time should be close to the end time specified
        in the evolve code. Or, when a collision was detected, it will be the
        model time of the collision.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('time', dtype='float64', direction=function.OUT,
            description = "The current model time")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time was retrieved
        -1 - ERROR
            The code does not have support for querying the time
        """
        return function
        


    @legacy_function
    def get_time_step():
        """
        Retrieve the model timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('code_index', dtype='int32', direction=function.IN, description = "Index of the code in rebound", default = 0)
        function.addParameter('time_step', dtype='float64', direction=function.OUT,
            description = "The current model timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time step was retrieved
        -1 - ERROR
            The code does not have support for querying the time
        """
        return function


    @legacy_function
    def new_subset():
        """
        Create a new particle subset (and corresponding code). This subset will evolve seperately from others.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_subset', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created subset
            """
            )

        function.addParameter('time_offset', dtype='float64', direction=function.IN, description = "Time of the system (defaults to the current model time)", default = -1)
        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            code was created
        -1 - ERROR
            code could not be created"""
        return function
        
        
    @legacy_function
    def get_state():
        """
        Retrieve the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, radius, position and velocity) is returned.
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
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
        function.addParameter('subset', dtype='int32', direction=function.OUT, description = "The current subset of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_subset():
        """
        Retrieve the subset index of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the subset of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('subset', dtype='int32', direction=function.OUT, description = "The current subset of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retreived
        -1 - ERROR
            particle could not be found
        """
        return function


    @legacy_function
    def stop_subset():
        """
        Stop a subset code
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_subset', dtype='int32', direction=function.IN, description =
            """
            An index assigned to an existing subset
            """
            )

        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            subset evolving was stopped
        -1 - ERROR
            subset evolving was already stopped"""
        return function
        
        

    @legacy_function
    def set_subset():
        """
        Retrieve the subset index of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the subset of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('subset', dtype='int32', direction=function.IN, description = "The new subset of the particle, as this is actually read only this will fail if changed!")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retreived
        -1 - ERROR
            particle could not be found
        """
        return function
    


class Rebound(GravitationalDynamics, GravityFieldCode):

    __interface__ = ReboundInterface


    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = self.__interface__(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    
    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        #GravityFieldCode.define_state(self, handler)
        self.stopping_conditions.define_state(handler)
        
        handler.add_method('EDIT', 'new_subset')
        handler.add_method('RUN', 'new_subset')
        


    def define_parameters(self, handler):
        self.stopping_conditions.define_parameters(handler)
        GravitationalDynamics.define_parameters(self, handler)
        
        handler.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            default_value = 0.0001 | nbody_system.time
        )


        handler.add_method_parameter(
            "get_integrator",
            "set_integrator",
            "integrator",
            "name of the integrator to use ({0})".format(sorted(self.INTEGRATORS.keys())), 
            default_value = "ias15"
        )


        handler.add_method_parameter(
            "get_solver",
            "set_solver",
            "solver",
            "name of the gravity solver to use ({0})".format(sorted(self.SOLVERS.keys())), 
            default_value = "compensated"
        )


        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )


        handler.add_method_parameter(
            "get_opening_angle2",
            "set_opening_angle2",
            "opening_angle2",
            "opening angle, theta, for building the tree in case of tree solver: between 0 and 1", 
            default_value = 0.5
        )

        handler.add_method_parameter(
            "get_boundary",
            "set_boundary",
            "boundary",
            "name of the boundary type to use ({0}) (required for tree solver)".format(sorted(self.BOUNDARIES.keys())), 
            default_value = "none"
        )

        handler.add_method_parameter(
            "get_boundary_size",
            "set_boundary_size",
            "boundary_size",
            "size of the boundaries, if the type is not none", 
            default_value = 1.0 | nbody_system.length
        )



    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        handler.add_method(
            "new_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                handler.NO_UNIT,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )                                                                  
        
        handler.add_method(
            "get_potential_energy",
            (handler.INDEX,),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_kinetic_energy",
            (handler.INDEX,),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, handler.ERROR_CODE,)
        )
        handler.add_method(
            'evolve_model',
            (
                nbody_system.time,
                handler.INDEX
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            'get_time',
            (handler.INDEX,),
            (nbody_system.time, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_time_step",
            (handler.INDEX,),
            (nbody_system.time, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_time_step",
            (nbody_system.time, handler.INDEX,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_state",
            (
                handler.NO_UNIT,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                handler.NO_UNIT,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_subset",
            (
                handler.NO_UNIT,
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_subset",
            (
                handler.NO_UNIT,
                handler.NO_UNIT,
            ),
            (
                
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            'new_subset',
            (
                nbody_system.time,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        self.stopping_conditions.define_methods(handler)
    


    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        
        self.stopping_conditions.define_particle_set(handler)
        
        handler.add_getter('particles', 'get_subset')
        handler.add_setter('particles', 'set_subset')

