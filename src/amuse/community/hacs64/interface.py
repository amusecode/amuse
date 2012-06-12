from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

# *** This script, together with the defaults in
# *** GravitationalDynamicsInterface, will be used to generate both
# *** the header file interface.h and the stub interface.cc.  Since
# *** interface.cc has been  hand-coded to implement the details,
# *** MAKE SURE TO SAVE IT SOMEWHERE, as build.py can overwrite it!

class Hacs64Interface(CodeInterface, GravitationalDynamicsInterface, StoppingConditionInterface):
    """
    HACS64, GPU-accelerated Hermite Ahmad-Cohen Scheme
    6th order irregular & 4th order regular step
    """

    # Interface specification.

    include_headers = ['interface.h', 'stopcond.h']

    MODE_GPU = 'gpu'
    MODE_CPU = 'cpu'

    def __init__(self, mode = MODE_CPU, **options):
        CodeInterface.__init__(
            self,
            name_of_the_worker=self.name_of_the_muse_worker(mode),
            **options
        )
        
    # Inheritance from GravitationalDynamicsInterface means that
    # functions in the standard interface don't need to be defined.
    # See interface.py.2 for a laboriously hand-coded version written
    # before I discovered this fact!    (Steve McMillan, 10/10)

    # Additional functions defined here will be reflected in
    # interface.h and must be provided in interface.cc in order for
    # ph4_worker to build.

    # The following functions aren't defined in the default interface:
    
    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_CPU:
            return 'hacs64_worker'
        elif mode == self.MODE_GPU:
            return 'hacs64_worker_gpu'
        else:
            return 'hacs64_worker'
    
    @legacy_function
    def get_number_of_particles_updated():
        """
        Return the number of particles added or deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_id_of_updated_particle():
        """
        Return the number of particles added or deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_update', dtype='int32',
                              direction=function.IN, 
                 description = 'index in the updated particles list')
        function.addParameter('index_of_particle', dtype='int32',
                              direction=function.OUT)
        function.addParameter('kind_of_update', dtype='int32',
                              direction=function.OUT,
                 description = 'kind of update (2, addition), (1, deletion)')
        function.can_handle_array = True
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_nmax():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('nmax', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nmax():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('nmax', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_dtmax():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dtmax', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dtmax():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dtmax', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_h2max():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('h2max', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_h2max():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('h2max', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eta_reg():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta_reg', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta_reg():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta_reg', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_eta_irr():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta_irr', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta_irr():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta_irr', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

        

class Hacs64(GravitationalDynamics):

    # The actual module.

    def __init__(self, convert_nbody = None, **keyword_arguments):

        legacy_interface = Hacs64Interface(**keyword_arguments)
        
        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(self,
                                       legacy_interface,
                                       convert_nbody,
                                       **keyword_arguments)

    def define_parameters(self, object):

        # Set/get parameters specific to the module, not part of the
        # standard interface.  Accessors used here must be defined
        # above and reflected in interface.cc.  Python access is
        # (e.g.)
        #
        #        ph4.parameters.timestep_parameter = xxx

        object.add_method_parameter(
            "get_nmax",                        # getter name in interface.cc
            "set_nmax",                        # setter name in interface.cc
            "nmax",                            # python parameter name
            "maximal number of particles",     # description
            default_value = -1    # default
        )

        object.add_method_parameter(
            "get_dtmax",                        # getter name in interface.cc
            "set_dtmax",                        # setter name in interface.cc
            "dtmax",                            # python parameter name
            "maximal timestep",     # description
            default_value = 0.0625 | nbody_system.time    # default
        )

        object.add_method_parameter(
            "get_h2max",                        # getter name in interface.cc
            "set_h2max",                        # setter name in interface.cc
            "h2max",                            # python parameter name
            "maximal neighbour radius squared", # description
            default_value = 0.5 | nbody_system.length*nbody_system.length   # default
        )

        object.add_method_parameter(
            "get_eta_reg",                        # getter name in interface.cc
            "set_eta_reg",                        # setter name in interface.cc
            "eta_reg",        # python parameter name
            "regular timestep parameter",        # description
            default_value = 0.1                # default
        )
        
        object.add_method_parameter(
            "get_eta_irr",                         # getter name in interface.cc
            "set_eta_irr",                         # setter name in interface.cc
            "eta_irr",              # python parameter name
            "irregular timestep parameter",        # description
            default_value = 0.6
        )

        object.add_method_parameter(
            "get_eps2",                        # already defined in standard interface
            "set_eps2",                        # already defined in standard interface
            "eps2",                             # python parameter name
            "smoothing parameter for gravity calculations", 
            default_value = 0.0  | nbody_system.length * nbody_system.length
        )
        
        self.stopping_conditions.define_parameters(object)

    def update_particle_set(self):
        """
        update the particle set after changes in the code
        
        this implementation needs to move to the
        amuse.support.data.incode_storage module, as
        it uses a lot of internal methods and info!
        
        """
        number_of_updated_particles = self.get_number_of_particles_updated()
        
        if number_of_updated_particles == 0:
            return
        
        indices_in_update_list = range(number_of_updated_particles)
        particle_indices, updates = self.get_id_of_updated_particle(indices_in_update_list)
        
        incode_storage = self.particles._private.attribute_storage
        
        indices_to_remove = []
        indices_to_add = []
        for index, status in zip(particle_indices, updates):
            if status == 1:                        # deletion
                indices_to_remove.append(index)
            elif status == 2:                        # addition
                indices_to_add.append(index)

        print ''
        print "indices_to_remove:", indices_to_remove
        print "indices_to_add:", indices_to_add

        if len(indices_to_remove) > 0:
            incode_storage._remove_indices(indices_to_remove)
        if len(indices_to_add) > 0:
            incode_storage._add_indices(indices_to_add)
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        # Similarly, we can add module-specific methods, if desired.
        # See hermite0/interface.py for examples.

        object.add_method("get_nmax", (),
            (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_nmax", (object.NO_UNIT, ),
            (object.ERROR_CODE,))
        
        object.add_method("get_dtmax", (),
            (nbody_system.time, object.ERROR_CODE,))
        object.add_method("set_dtmax", (nbody_system.time, ),
            (object.ERROR_CODE,))

        object.add_method("get_h2max", (),
            (nbody_system.length*nbody_system.length, object.ERROR_CODE,))
        object.add_method("set_h2max", (nbody_system.length*nbody_system.length, ),
            (object.ERROR_CODE,))

        object.add_method("get_eta_irr", (),
            (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_eta_irr", (object.NO_UNIT, ),
            (object.ERROR_CODE,))

        object.add_method("get_eta_reg", (),
            (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_eta_reg", (object.NO_UNIT, ),
            (object.ERROR_CODE,))

        object.add_method("get_eps2", (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,))
        object.add_method("set_eps2", (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,))
        
        self.stopping_conditions.define_methods(object)
    

    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        
        self.stopping_conditions.define_particle_set(object)

