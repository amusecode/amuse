import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option


class SimpleXInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
    SimpleX(2) is a method for radiative transfer on an unstructured Delaunay 
    grid. The grid samples the medium through which photons are transported in 
    an optimal way for fast radiative transfer calculations.
    
    The relevant references are:
        .. [#] Paardekooper, J.-P., Kruip, C.J.H., Icke, V. 2010, A&A, 515, A79 (SimpleX2)
        .. [#] Ritzerveld, J., & Icke, V. 2006, Phys. Rev. E, 74, 26704 (SimpleX)
    """
    include_headers=['worker.h']
    
    def __init__(self, **kwargs):
        CodeInterface.__init__(self, name_of_the_worker = "simplex_worker", **kwargs)
        LiteratureReferencesMixIn.__init__(self)
    
    @option(type="string")
    def data_directory(self):
        """
        The root name of the directory for the SimpleX
        application data files.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'simplex', 'input')
    
    @option(type="string")
    def output_directory(self):
        """
        The root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'simplex', 'output')
    
    @legacy_function
    def set_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('output_path', dtype='string', direction=function.IN,
            description = "Name of the output directory")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def commit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function

    @legacy_function
    def recommit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        for x in ['x','y','z','rho','flux','xion','u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for x in ['x','y','z','rho','flux','xion','u']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
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
    def get_flux():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the flux from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('flux', dtype='float64', direction=function.OUT, description = "The current flux of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function
    
    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the density from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('rho', dtype='float64', direction=function.OUT, description = "The current density of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the internal energy from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('u', dtype='float64', direction=function.OUT, description = "The current internal energy of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    
    @legacy_function
    def get_ionisation():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the ionisation from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('xion', dtype='float64', direction=function.OUT, description = "The current ionisation of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function
    
    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for x in ['x','y','z','rho','flux','xion','u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
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
    def set_flux():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the flux for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('flux', dtype='float64', direction=function.IN, description = "The new flux of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
        return function
    
    @legacy_function
    def set_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the density for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('rho', dtype='float64', direction=function.IN, description = "The new density of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the internal energy for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('u', dtype='float64', direction=function.IN, description = "The new internal energy of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
        return function

    
    @legacy_function
    def set_ionisation():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the ionisation for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('xion', dtype='float64', direction=function.IN, description = "The new ionisation of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            new value was set
        -1 - ERROR
            particle could not be found
        """
        return function
    
    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN)
#        function.addParameter('synchronize', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
   
    
    @legacy_function
    def set_box_size_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('box_size', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_box_size_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('box_size', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_hilbert_order_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('hilbert_order', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_hilbert_order_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('hilbert_order', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    
    @legacy_function
    def get_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_source_Teff():
        function = LegacyFunctionSpecification()
        function.addParameter('sourceTeff', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_source_Teff():
        function = LegacyFunctionSpecification()
        function.addParameter('sourceTeff', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_number_frequency_bins():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_freq_bins', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_frequency_bins():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_freq_bins', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_thermal_evolution():
        function = LegacyFunctionSpecification()
        function.addParameter('thermal_evolution_flag', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_thermal_evolution():
        function = LegacyFunctionSpecification()
        function.addParameter('thermal_evolution_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    
    

    def invoke_state_change2(self):
        pass
    
    def synchronize_model(self):
        pass
    


class GlSimpleXInterface(SimpleXInterface):
    include_headers=['worker_gl.h']
    
    def __init__(self, **kwargs):
        CodeInterface.__init__(self,debug_with_gdb=False,
           name_of_the_worker = 'simplex_worker_gl', **kwargs)
    
    @legacy_function
    def start_viewer():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def stop_viewer():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    

class SimpleXDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_interface.__doc__+"\n\n"+instance.parameters.__doc__


class SimpleX(CommonCode):
    
    __doc__ = SimpleXDoc()
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, SimpleXInterface(**options))
        self.set_output_directory(self.output_directory)
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_timestep_parameter",
            "set_timestep_parameter", 
            "timestep", 
            "timestep for radiative transfer sweeps", 
            default_value = 0.05 | units.Myr
        )

        object.add_method_parameter(
            "get_hilbert_order_parameter",
            "set_hilbert_order_parameter", 
            "hilbert_order", 
            "hilbert_order for domain decomposition", 
            default_value = 1 | units.none
        )

        object.add_method_parameter(
            "get_number_frequency_bins",
            "set_number_frequency_bins", 
            "number_of_freq_bins", 
            "the number of bins of frequency", 
            default_value = 1 | units.none
        )


        object.add_method_parameter(
            "get_box_size_parameter",
            "set_box_size_parameter", 
            "box_size", 
            "boxsize for radiative transfer particle distribution", 
            default_value = 13200. | units.parsec
        )


    def define_methods(self, object):
        CommonCode.define_methods(self, object)
        object.add_method('evolve_model', (units.Myr,), ( object.ERROR_CODE, ))
        object.add_method(
            "new_particle",
            (
                units.parsec,
                units.parsec,
                units.parsec,
                units.amu / units.cm**3,
                1.0e48 / units.s,
                units.none,
                units.cm**2 / units.s**2
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
                units.parsec,
                units.parsec,
                units.parsec,
                units.amu / units.cm**3,
                1.0e48 / units.s,
                units.none,
                units.cm**2 / units.s**2,
                object.ERROR_CODE
                )
            )
        object.add_method(
            "set_state",
            (
                object.NO_UNIT,
                units.parsec,
                units.parsec,
                units.parsec,
                units.amu / units.cm**3,
                1.0e48 / units.s,
                units.none,
                units.cm**2 / units.s**2,
                ),
            (
                object.ERROR_CODE
                )
            )
        object.add_method(
            "get_internal_energy",
            (
                object.NO_UNIT,
                ),
            (
                units.cm**2 / units.s**2,
                object.ERROR_CODE
                )
            )
        object.add_method(
            "set_internal_energy",
            (
                object.NO_UNIT,
                units.cm**2 / units.s**2,
                ),
            (
                object.ERROR_CODE
                )
            )
        object.add_method(
            "set_position",
            (
                object.NO_UNIT,
                units.parsec,
                units.parsec,
                units.parsec,
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
                units.parsec,
                units.parsec,
                units.parsec,
                object.ERROR_CODE
                )
            )
        object.add_method(
            "set_density",
            (
                object.NO_UNIT,
                units.amu / units.cm**3,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_density",
            (
                object.NO_UNIT,
            ),
            (
                units.amu / units.cm**3,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_flux",
            (
                object.NO_UNIT,
                1.0e48 / units.s,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_flux",
            (
                object.NO_UNIT,
            ),
            (
                1.0e48 / units.s,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_ionisation",
            (
                object.NO_UNIT,
                units.none
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_ionisation",
            (
                object.NO_UNIT,
            ),
            (
                units.none,
                object.ERROR_CODE
            )
        )
        object.add_method(
            'commit_particles',
            (),
            (object.ERROR_CODE)
        )
        object.add_method(
            'recommit_particles',
            (),
            (object.ERROR_CODE)
        )
    
        object.add_method(
            "get_timestep_parameter",
            (),
            (units.Myr, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_timestep_parameter",
            (units.Myr, ),
            (object.ERROR_CODE,)
        )
    
        object.add_method(
            "get_hilbert_order_parameter",
            (),
            (units.none, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_hilbert_order_parameter",
            (units.none, ),
            (object.ERROR_CODE,)
        )
    
        object.add_method(
            "get_box_size_parameter",
            (),
            (units.parsec, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_box_size_parameter",
            (units.parsec, ),
            (object.ERROR_CODE,)
        )
    
    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')
        object.add_setter('particles', 'set_state')
        object.add_getter('particles', 'get_state')
        object.add_setter('particles', 'set_position')
        object.add_getter('particles', 'get_position')
        object.add_setter('particles', 'set_density')
        object.add_getter('particles', 'get_density')
        object.add_setter('particles', 'set_flux')
        object.add_getter('particles', 'get_flux')
        object.add_setter('particles', 'set_ionisation')
        object.add_getter('particles', 'get_ionisation')
        object.add_setter('particles', 'set_internal_energy')
        object.add_getter('particles', 'get_internal_energy')

    
    def define_state(self, object):
        CommonCode.define_state(self, object)
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','PARAMETER_CHANGE_A','invoke_state_change2')
        object.add_transition('EDIT','PARAMETER_CHANGE_B','invoke_state_change2')
        object.add_transition('PARAMETER_CHANGE_A','RUN','recommit_parameters')
        object.add_transition('PARAMETER_CHANGE_B','EDIT','recommit_parameters')
        object.add_method('EDIT', 'new_particle')
        object.add_method('EDIT', 'delete_particle')
        object.add_transition('EDIT', 'RUN', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_particle', False)
        object.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        object.add_method('EVOLVED', 'evolve_model')
        object.add_transition('EVOLVED','RUN', 'synchronize_model')
        object.add_method('RUN', 'synchronize_model')
        object.add_method('RUN', 'get_state')
        object.add_method('RUN', 'get_density')
        object.add_method('RUN', 'get_position')
        object.add_method('RUN', 'get_flux')
        object.add_method('RUN', 'get_ionisation')
        object.add_method('RUN', 'get_internal_energy')


    

