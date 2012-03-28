from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics 
from amuse.community.interface.gd import GravitationalDynamicsInterface

class HiGPUsInterface(CodeInterface, GravitationalDynamicsInterface):
    
    include_headers = ['worker_code.h']
	    
    def __init__(self, **keyword_arguments):
	CodeInterface.__init__(self, name_of_the_worker="higpus_worker_gpu", **keyword_arguments)
    
    @legacy_function
    def echo_int():
	function = LegacyFunctionSpecification()  
      	function.addParameter('int_in', dtype='int32', direction=function.IN)
      	function.addParameter('int_out', dtype='int32', direction=function.OUT)
      	function.result_type = 'int32'
      	function.can_handle_array = True
      	return function

    @legacy_function
    def new_particle():
      	function = LegacyFunctionSpecification()
      	function.can_handle_array = True
      	function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
      	function.addParameter('mass', dtype='float64', direction=function.IN,
		 description = "The mass of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN,
                 description = "The radius of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN,
                 description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN,
                 description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN,
                 description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN,
                 description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN,
                 description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN,
                 description = "The initial velocity vector of the particle")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('dt', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time_begin():
	function = LegacyFunctionSpecification()
        function.addParameter('time_begin', dtype='float64', direction=function.IN, description = "start time.")
	function.result_type = 'int32'
        return function

    @legacy_function
    def get_time_begin():
	function = LegacyFunctionSpecification()
        function.addParameter('time_begin', dtype='float64', direction=function.OUT, description = "start time.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_eta6():
	function = LegacyFunctionSpecification()
        function.addParameter('eta6', dtype='float64', direction=function.IN, description = "eta parameter of time steps.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta6():
	function = LegacyFunctionSpecification()
        function.addParameter('eta6', dtype='float64', direction=function.OUT, description = "eta parameter of time steps.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eta4():
        function = LegacyFunctionSpecification()
        function.addParameter('eta4', dtype='float64', direction=function.IN, description = "eta parameter of time steps.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta4():
        function = LegacyFunctionSpecification()
        function.addParameter('eta4', dtype='float64', direction=function.OUT, description = "eta parameter of time steps.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_Plummer_core():
        function = LegacyFunctionSpecification()
        function.addParameter('Plummer_core', dtype='float64', direction=function.IN, description = "radius parameter for plummer potential.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_Plummer_core():
        function = LegacyFunctionSpecification()
        function.addParameter('Plummer_core', dtype='float64', direction=function.OUT, description = "radius parameter for plummer potential.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_Plummer_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('Plummer_mass', dtype='float64', direction=function.IN, description = "analytical mass of the plummer potential.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_Plummer_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('Plummer_mass', dtype='float64', direction=function.OUT, description = "analytical mass of the plummer potential.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eps():
	function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.IN, description = "softening parameter.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eps():
	function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.OUT, description = "softening parameter.")
        function.result_type = 'int32'
        return function
   
    @legacy_function
    def set_number_of_GPU():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_GPU', dtype='int32', direction=function.IN, description = "number of GPU.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_GPU():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_GPU', dtype='int32', direction=function.OUT, description = "number of GPU.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_number_of_Threads():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_Threads', dtype='int32', direction=function.IN, description = "number of Threads per block.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_Threads():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_Threads', dtype='int32', direction=function.OUT, description = "number of Threads per block.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_number_of_Print():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_Print', dtype='int32', direction=function.IN, description = "number of total file to print.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_Print():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_Print', dtype='int32', direction=function.OUT, description = "number of total file to print.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.IN, description = "number of particles.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT, description = "number of particles.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_DTPrint():
        function = LegacyFunctionSpecification()
        function.addParameter('DTPrint', dtype='float64', direction=function.IN, description = "unit time of snapshot.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_DTPrint():
	function = LegacyFunctionSpecification()
        function.addParameter('DTPrint', dtype='float64', direction=function.OUT, description = "unit time of snapshot.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_max_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('max_time_step', dtype='float64', direction=function.IN, description = "exponent of the maximum time step used.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_max_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('max_time_step', dtype='float64', direction=function.OUT, description = "maximum time step used.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_min_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('min_time_step', dtype='float64', direction=function.IN, description = "exponent of the minimum time step used.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_min_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('min_time_step', dtype='float64', direction=function.OUT, description = "minimum time step used.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_gpu_name():
        function = LegacyFunctionSpecification()
        function.addParameter('gpu_name', dtype='string', direction=function.IN, description = "name of the GPU used.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gpu_name():
        function = LegacyFunctionSpecification()
        function.addParameter('gpu_name', dtype='string', direction=function.OUT, description = "name of the GPU used.")
        function.result_type = 'int32'
        return function



class HiGPUs(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **keyword_arguments):
        legacy_interface = HiGPUsInterface(**keyword_arguments)

        GravitationalDynamics.__init__(self,  
                                       legacy_interface,
                                       convert_nbody,
                                       **keyword_arguments)


    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eta4",			          # getter name in interface.cc
            "set_eta4",			          # setter name in interface.cc
            "eta_4",                              # python parameter name
            "timestep parameter",                 # description
            default_value = 0.01 | units.none)

        object.add_method_parameter(
            "get_eta6",                           # getter name in interface.cc
            "set_eta6",                           # setter name in interface.cc
            "eta_6",                              # python parameter name
            "timestep parameter",                 # description
            default_value = 0.4 | units.none)

        object.add_method_parameter(
            "get_time_begin",                     # getter name in interface.cc
            "set_time_begin",                     # setter name in interface.cc
            "start_time",                         # python parameter name
            "start time",                         # description
            default_value = 0.0 | nbody_system.time
        )

	object.add_method_parameter(
            "get_Plummer_core",                   # getter name in interface.cc
            "set_Plummer_core",                   # setter name in interface.cc
            "r_core_plummer",                     # python parameter name
            "radius of Plummer potential",        # description
            default_value = 0.0 | nbody_system.length
        )

        object.add_method_parameter(
            "get_Plummer_mass",                   # getter name in interface.cc
            "set_Plummer_mass",                   # setter name in interface.cc
            "mass_plummer",                       # python parameter name
            "mass of Plummer potential",          # description
            default_value = 0.0 | nbody_system.mass
        )

        object.add_method_parameter(
            "get_eps",		                  # getter name in interface.cc
            "set_eps",          	          # setter name in interface.cc
            "softening",                          # python parameter name
            "softening",           	          # description
            default_value = 0.001 | nbody_system.length
        )

	object.add_method_parameter(
            "get_number_of_Threads",              # getter name in interface.cc
            "set_number_of_Threads",              # setter name in interface.cc
            "Threads",                            # python parameter name
	    "Threads per block",                  # description
            default_value = 128 | units.none
        )

        object.add_method_parameter(
            "get_number_of_particles",            # getter name in interface.cc
            "set_number_of_particles",            # setter name in interface.cc
            "N",                                  # python parameter name
            "number of particles",                # description
            default_value = 8192 | units.none
        )


        object.add_method_parameter(
            "get_number_of_Print",                # getter name in interface.cc
            "set_number_of_Print",                # setter name in interface.cc
            "n_Print",                            # python parameter name
	    "start number to print file",         # description
            default_value = 1000000 | units.none
        )

        object.add_method_parameter(
            "get_DTPrint",                        # getter name in interface.cc
            "set_DTPrint",                        # setter name in interface.cc
            "dt_Print",                           # python parameter name
	    "time for snapshot",                  # description
            default_value = 1.0 | nbody_system.time
        )

        object.add_method_parameter(
            "get_max_time_step",                  # getter name in interface.cc
            "set_max_time_step",                  # setter name in interface.cc
            "max_step",                           # python parameter name
	    "power of 2 for maximum time step",   # description
            default_value = -3.0 | nbody_system.time
        )

	object.add_method_parameter(
            "get_min_time_step",                  # getter name in interface.cc
            "set_min_time_step",                  # setter name in interface.cc
            "min_step",                           # python parameter name
            "power of 2 for minmum time step",    # description
            default_value = -30.0 | nbody_system.time
        )

	object.add_method_parameter(
            "get_gpu_name",                       # getter name in interface.cc
            "set_gpu_name",                       # setter name in interface.cc
            "gpu_name",                           # python parameter name
            "gpu name",                           # description
            default_value = "GeForce GTX 480" | units.none
        )

        object.add_method_parameter(
            "get_number_of_GPU",                  # getter name in interface.cc
            "set_number_of_GPU",                  # setter name in interface.cc
            "n_gpu",                              # python parameter name
            "number of gpus per node",            # description
            default_value = 2 | units.none
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
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
	        nbody_system.speed,
           ),
            (
               	object.INDEX,
                object.ERROR_CODE
            )
	)
        object.add_method(
            "set_time_begin",
            (
             	nbody_system.time
            ),
            (
             	object.ERROR_CODE
            )
	)
        object.add_method(
            "get_time_begin",
            (),
            (
             	nbody_system.time,
                object.ERROR_CODE
            )
	)
        object.add_method(
            "set_Plummer_core",
            (
               	nbody_system.length
            ),
            (
             	object.ERROR_CODE
            )
        )

        object.add_method(
            "get_Plummer_core",
            (),
            (
               	nbody_system.length,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "set_Plummer_mass",
            (
               	nbody_system.mass
            ),
            (
             	object.ERROR_CODE
            )
        )

        object.add_method(
            "get_Plummer_mass",
            (),
            (
               	nbody_system.mass,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "set_eta6",
            (
                units.none
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_eta6",
            (),
            (
                units.none,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "set_eta4",
            (
                units.none
            ),
            (
             	object.ERROR_CODE
            )
        )

        object.add_method(
            "get_eta4",
            (),
            (
                units.none,
                object.ERROR_CODE
            )
	)


        object.add_method(
            "set_eps",
            (
                nbody_system.length
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
	    "get_eps",
            (),
	    (
             	nbody_system.length,
                object.ERROR_CODE
            )
	)


	object.add_method(
            "set_number_of_GPU",
            (
             	units.none
            ),
            (
             	object.ERROR_CODE
            )
	)

	object.add_method(
            "get_number_of_GPU",
            (),
            (
             	units.none,
                object.ERROR_CODE
            )
	)

        object.add_method(
            "set_number_of_particles",
            (
                units.none
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_number_of_particles",
            (),
            (
                units.none,
                object.ERROR_CODE
            )
        )


        object.add_method(
            "set_number_of_Threads",
            (
                units.none
            ),
            (
             	object.ERROR_CODE
            )
        )

        object.add_method(
            "get_number_of_Threads",
            (),
            (
             	units.none,
                object.ERROR_CODE
            )
	)

	object.add_method(
            "set_number_of_Print",
            (
             	units.none
            ),
            (
             	object.ERROR_CODE
            )
	)

	object.add_method(
            "get_number_of_Print",
            (),
            (
             	units.none,
                object.ERROR_CODE
            )
	)

        object.add_method(
            "set_DTPrint",
            (
                nbody_system.time
            ),
            (
             	object.ERROR_CODE
            )
	)

	object.add_method(
            "get_DTPrint",
            (),
            (
             	nbody_system.time,
                object.ERROR_CODE
            )
	)

        object.add_method(
            "set_max_time_step",
            (
             	nbody_system.time
	    ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_max_time_step",
            (),
            (
                nbody_system.time,
	        object.ERROR_CODE
            )
	)

	object.add_method(
            "set_min_time_step",
            (
             	nbody_system.time
            ),
            (
             	object.ERROR_CODE
            )
	)

        object.add_method(
            "get_min_time_step",
            (),
            (
             	nbody_system.time,
                object.ERROR_CODE
            )
	)

	object.add_method(
            "set_gpu_name",
            (
                units.none
            ),
            (
             	object.ERROR_CODE
            )
	)

	object.add_method(
            "get_gpu_name",
            (),
            (
             	units.none,
                object.ERROR_CODE
            )
	)


    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)


