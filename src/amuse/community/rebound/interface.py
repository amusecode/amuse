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
    

    INTEGRATORS = {"ias15": 0, "whfast": 1, "sei": 2, "wh": 3, "leapfrog": 4, "hybrid": 5, "none": 6}
    def set_integrator(self, name, code_index = 0 ):
        return self._set_integrator(self.INTEGRATORS[name], code_index)
    
    def get_integrator(self, code_index = 0):
        value, error = self._get_integrator(code_index)
        for key, index in self.INTEGRATORS.iteritems():
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
        print "set_solver name: %s code_index: %i"%(name,code_index)
        print self.SOLVERS[name]
        return self._set_solver(self.SOLVERS[name], code_index)
    
    def get_solver(self, code_index = 0):
        print "get_solver code_index: %i"%code_index
        value, error = self._get_solver(code_index)
        for key, index in self.SOLVERS.iteritems():
            if value == index:
                return key, error
        return "none", error

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
    __so_module__ = 'rebound_cython'


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

    
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        #GravityFieldCode.define_state(self, object)
        self.stopping_conditions.define_state(object)
        
        object.add_method('EDIT', 'new_subset')
        object.add_method('RUN', 'new_subset')
        


    def define_parameters(self, object):
        self.stopping_conditions.define_parameters(object)
        GravitationalDynamics.define_parameters(self, object)
        
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            default_value = 0.0001 | nbody_system.time
        )


        object.add_method_parameter(
            "get_integrator",
            "set_integrator",
            "integrator",
            "name of the integrator to use ({0})".format(sorted(self.INTEGRATORS.keys())), 
            default_value = "ias15"
        )


        object.add_method_parameter(
            "get_solver",
            "set_solver",
            "solver",
            "name of the gravity solver to use ({0})".format(sorted(self.SOLVERS.keys())), 
            default_value = "compensated"
        )


        object.add_method_parameter(
            "get_opening_angle2",
            "set_opening_angle2",
            "opening_angle2",
            "opening angle, theta, for building the tree in case of tree solver: between 0 and 1", 
            default_value = 0.5
        )


    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        object.add_method(
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
                object.NO_UNIT,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )                                                                  
        
        object.add_method(
            "get_potential_energy",
            (object.INDEX,),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, object.ERROR_CODE,)
        )


        object.add_method(
            "get_kinetic_energy",
            (object.INDEX,),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, object.ERROR_CODE,)
        )
        object.add_method(
            'evolve_model',
            (
                nbody_system.time,
                object.INDEX
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            'get_time',
            (object.INDEX,),
            (nbody_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time_step",
            (object.INDEX,),
            (nbody_system.time, object.ERROR_CODE,)
        )

        object.add_method(
            "set_time_step",
            (nbody_system.time, object.INDEX,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_state",
            (
                object.NO_UNIT,
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
                object.NO_UNIT,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_subset",
            (
                object.NO_UNIT,
            ),
            (
                object.NO_UNIT,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_subset",
            (
                object.NO_UNIT,
                object.NO_UNIT,
            ),
            (
                
                object.ERROR_CODE,
            )
        )
        self.stopping_conditions.define_methods(object)
    


    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        
        self.stopping_conditions.define_particle_set(object)
        
        object.add_getter('particles', 'get_subset')
        object.add_setter('particles', 'set_subset')

