import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option

from amuse.datamodel import Particles

class SimpleXInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
    SimpleX(2.5) is a method for radiative transfer on an unstructured Delaunay 
    grid. The grid samples the medium through which photons are transported in 
    an optimal way for fast radiative transfer calculations.
    
    The relevant references are:
        .. [#] Kruip, C.J.H., Ph. D. thesis, University of Leiden (2011)
        .. [#] Paardekooper, J.-P., Ph. D. thesis, University of Leiden (2010)
        .. [#] Paardekooper, J.-P., Kruip, C.J.H., Icke, V. 2010, A&A, 515, A79 (SimpleX2)
        .. [#] Ritzerveld, J., & Icke, V. 2006, Phys. Rev. E, 74, 26704 (SimpleX)
    """
    include_headers=['worker_code.h']
    
    def __init__(self, **kwargs):
        CodeInterface.__init__(self, name_of_the_worker = "simplex_worker", **kwargs)
        LiteratureReferencesMixIn.__init__(self)
    
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
        The root name of the directory for the SimpleX
        application data files.
        """
        return os.path.join(self.input_data_root_directory, 'simplex', 'input')
    
    @option(type="string")
    def output_directory(self):
        """
        The root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, 'simplex', 'output')
    
    @legacy_function
    def set_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('output_path', dtype='string', direction=function.IN,
            description = "Name of the output directory")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_data_directory():
        """
        Update the path to the SimpleX database.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('data_directory', dtype='string', direction=function.IN,
            description = "Name of the SimpleX data directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function
    
    @legacy_function
    def get_data_directory():
        """
        Retrieve the path to the SimpleX database currently used.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('data_directory', dtype='string', direction=function.OUT,
            description = "Name of the SimpleX data directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Value was retrieved
        -1 - ERROR
            Could not retrieve value
        """
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
        function.addParameter('metallicity', dtype='float64',direction=function.IN, default=0.0 )
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
        function.addParameter('metallicity', dtype='float64', direction=function.OUT)    
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
    def get_mean_intensity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the mean intensity from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('Js', dtype='float64', direction=function.OUT, description = "The current mean intensity of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function


    @legacy_function
    def get_diffuse_intensity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the diffuse intensity from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('Jd', dtype='float64', direction=function.OUT, description = "The current diffuse intensity of the particle")
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
    def get_dinternal_energy_dt():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the d/dt internal energy from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('du_dt', dtype='float64', direction=function.OUT, description = "The current d/dt internal energy of the particle")
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
    def get_metallicity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the metallicity from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('metallicity', dtype='float64', direction=function.OUT, description = "The current metallicity of the particle")
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
        function.addParameter('metallicity',dtype='float64',direction=function.IN,default=0.0)    
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
    def set_dinternal_energy_dt():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the d/dt internal energy for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('du_dt', dtype='float64', direction=function.IN, description = "The new d/dt internal energy of the particle")
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
    def set_metallicity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to set the metallicity for. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('metallicity', dtype='float64', direction=function.IN, description = "The new metallicity of the particle")
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
    def set_box_size():
        function = LegacyFunctionSpecification()
        function.addParameter('box_size', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_box_size():
        function = LegacyFunctionSpecification()
        function.addParameter('box_size', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_hilbert_order():
        function = LegacyFunctionSpecification()
        function.addParameter('hilbert_order', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_hilbert_order():
        function = LegacyFunctionSpecification()
        function.addParameter('hilbert_order', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    
    @legacy_function
    def get_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
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
    def set_number_of_border_sites():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_border_sites', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_border_sites():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_border_sites', dtype='int32', direction=function.OUT)
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

    @legacy_function
    def set_blackbody_spectrum():
        function = LegacyFunctionSpecification()
        function.addParameter('blackbody_spectrum_flag', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_blackbody_spectrum():
        function = LegacyFunctionSpecification()
        function.addParameter('blackbody_spectrum_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_metal_cooling():
        function = LegacyFunctionSpecification()
        function.addParameter('metal_cooling_flag', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_metal_cooling():
        function = LegacyFunctionSpecification()
        function.addParameter('metal_cooling_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    
    @legacy_function
    def set_recombination_radiation():
        function = LegacyFunctionSpecification()
        function.addParameter('recombination_radiation_flag', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_recombination_radiation():
        function = LegacyFunctionSpecification()
        function.addParameter('recombination_radiation_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    
    @legacy_function
    def set_collisional_ionization():
        function = LegacyFunctionSpecification()
        function.addParameter('collisional_ionisation_flag', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_collisional_ionization():
        function = LegacyFunctionSpecification()
        function.addParameter('collisional_ionisation_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    def synchronize_model(self):
        pass
        
    def auto_go_to_run(self):
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
        self.set_data_directory(self.data_directory)

    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")
    
    def define_parameters(self, object):
        
        object.add_method_parameter(
            "get_timestep",
            "set_timestep", 
            "timestep", 
            "timestep for radiative transfer sweeps", 
            default_value = 0.05 | units.Myr
        )
        
        object.add_method_parameter(
            "get_source_Teff",
            "set_source_Teff", 
            "source_effective_T", 
            "effective temperature for sources", 
            default_value = 1.e5 | units.K
        )


        object.add_method_parameter(
            "get_hilbert_order",
            "set_hilbert_order", 
            "hilbert_order", 
            "hilbert_order for domain decomposition", 
            default_value = 1
        )

        object.add_method_parameter(
            "get_number_frequency_bins",
            "set_number_frequency_bins", 
            "number_of_freq_bins", 
            "the number of bins of frequency", 
            default_value = 1
        )

        object.add_method_parameter(
            "get_thermal_evolution",
            "set_thermal_evolution", 
            "thermal_evolution_flag", 
            "solve full thermal evolution if 1", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_metal_cooling",
            "set_metal_cooling", 
            "metal_cooling_flag", 
            "include cooling from metals if 1, not if zero", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_recombination_radiation",
            "set_recombination_radiation", 
            "recombination_radiation_flag", 
            "include recombination radiation if 1, not if zero", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_blackbody_spectrum",
            "set_blackbody_spectrum", 
            "blackbody_spectrum_flag", 
            "monochromatic if 0, blackbody_spectrum if 1", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_collisional_ionization",
            "set_collisional_ionization", 
            "collisional_ionization_flag", 
            "not use collisional ionization if 0, do if 1", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_box_size",
            "set_box_size", 
            "box_size", 
            "boxsize for radiative transfer particle distribution", 
            default_value = 13200. | units.parsec
        )

        object.add_method_parameter(
            "get_data_directory", 
            "set_data_directory",
            "simplex_data_directory", 
            "Name of the SimpleX data directory", 
            default_value = "."
        )
        object.add_method_parameter(
            "get_number_of_border_sites", 
            "set_number_of_border_sites",
            "number_of_border_sites", 
            "Number of border sites to generate on the boundary, needs to be a large number. These sites will be placed randomly outside the problem domain and will culled to create a closed boundary", 
            default_value = 25000
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
                object.NO_UNIT,
                units.cm**2 / units.s**2,
                object.NO_UNIT
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
                object.NO_UNIT,
                units.cm**2 / units.s**2,
                object.NO_UNIT,
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
                object.NO_UNIT,
                units.cm**2 / units.s**2,
                object.NO_UNIT
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
            "get_dinternal_energy_dt",
            (
                object.NO_UNIT,
                ),
            (
                units.cm**2 / units.s**3,
                object.ERROR_CODE
                )
            )
        object.add_method(
            "set_dinternal_energy_dt",
            (
                object.NO_UNIT,
                units.cm**2 / units.s**3,
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
            "get_mean_intensity",
            (
                object.NO_UNIT,
            ),
            (
                1.0e48 / units.s,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_diffuse_intensity",
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
                object.NO_UNIT
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
                object.NO_UNIT,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_metallicity",
            (
                object.NO_UNIT,
                object.NO_UNIT
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_metallicity",
            (
                object.NO_UNIT,
            ),
            (
                object.NO_UNIT,
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
            "get_timestep",
            (),
            (units.Myr, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_timestep",
            (units.Myr, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_source_Teff",
            (),
            (units.K, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_source_Teff",
            (units.K, ),
            (object.ERROR_CODE,)
        )

    
        object.add_method(
            "get_hilbert_order",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_hilbert_order",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_blackbody_spectrum",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_blackbody_spectrum",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )


        object.add_method(
            "get_thermal_evolution",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_thermal_evolution",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_collisional_ionization",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_collisional_ionization",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_number_frequency_bins",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_number_frequency_bins",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )


        object.add_method(
            "get_metal_cooling",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_metal_cooling",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        
        object.add_method(
            "get_recombination_radiation",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_recombination_radiation",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
    
        object.add_method(
            "get_box_size",
            (),
            (units.parsec, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_box_size",
            (units.parsec, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_time',
            (),
            (units.s,object.ERROR_CODE,)
        )
        object.add_method(
            'set_time',
            (units.s,),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_data_directory",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_data_directory",
            (object.NO_UNIT,),
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
        object.add_getter('particles', 'get_mean_intensity')
        object.add_getter('particles', 'get_diffuse_intensity')
        object.add_setter('particles', 'set_ionisation')
        object.add_getter('particles', 'get_ionisation')
        object.add_setter('particles', 'set_metallicity')
        object.add_getter('particles', 'get_metallicity')
        object.add_setter('particles', 'set_internal_energy')
        object.add_getter('particles', 'get_internal_energy')
        object.add_setter('particles', 'set_dinternal_energy_dt')
        object.add_getter('particles', 'get_dinternal_energy_dt')

    
    def define_state(self, object):
        CommonCode.define_state(self, object)
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter')
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter')
        object.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        object.add_method('UPDATE','before_get_parameter')
        object.add_method('EVOLVED','before_get_parameter')
        
        
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
        object.add_method('RUN', 'get_mean_intensity')
        object.add_method('RUN', 'get_diffuse_intensity')
        object.add_method('RUN', 'get_ionisation')
        object.add_method('RUN', 'get_internal_energy')
        object.add_method('RUN', 'set_dinternal_energy_dt')
        object.add_method('RUN', 'get_dinternal_energy_dt')
        object.add_method('UPDATE', 'set_dinternal_energy_dt')
        object.add_method('UPDATE', 'get_dinternal_energy_dt')

        object.add_method('INITIALIZED', 'set_hilbert_order')
        object.add_method('INITIALIZED', 'set_box_size')
        object.add_method('INITIALIZED', 'set_timestep')
        object.add_method('INITIALIZED', 'set_source_Teff')
        object.add_method('INITIALIZED', 'set_number_frequency_bins')
        object.add_method('INITIALIZED', 'set_thermal_evolution')
        object.add_method('INITIALIZED', 'set_blackbody_spectrum')
        object.add_method('INITIALIZED', 'set_metal_cooling')
        object.add_method('INITIALIZED', 'set_recombination_radiation')
        object.add_method('INITIALIZED', 'set_collisional_ionization')


    
class SimpleXSplitSet(SimpleX):
    
    def __init__(self,**options):
        SimpleX.__init__(self,**options)
        self.gas_particles=Particles()
        self.src_particles=Particles()

    def commit_particles(self):
        
        sites=self.gas_particles.copy()
        sites.flux=0. | units.s**-1

        for p in self.src_particles:
          nearest=sites.find_closest_particle_to(p.x,p.y,p.z)
          nearest.flux+=p.luminosity

        self.particles.add_particles(sites)
        
        self.simplex_to_gas_channel=self.particles.new_channel_to(self.gas_particles)

        self.overridden().commit_particles()
        
        if hasattr(sites,"du_dt"):
          attributes=["du_dt"]
          channel=sites.new_channel_to(self.particles)
          channel.copy_attributes(attributes)

        del sites
              
    def recommit_particles(self):  

        sites=self.gas_particles.copy()
        sites.flux=0. | units.s**-1

        for p in self.src_particles:
          nearest=sites.find_closest_particle_to(p.x,p.y,p.z)
          nearest.flux+=p.luminosity

#        sites.synchronize_to(self.particles)
        add_set=sites.difference(self.particles)
        remove_set=self.particles.difference(sites)

        if len(remove_set)>0: self.particles.remove_particles(remove_set)
        if len(add_set)>0: self.particles.add_particles(add_set)
        self.overridden().recommit_particles() 

        channel=sites.new_channel_to(self.particles)
        attributes=["x","y","z","rho","xion","u","flux"]
        if hasattr(sites,"metallicity"):
           attributes.append("metallicity") 
        if hasattr(sites,"du_dt"):
           attributes.append("du_dt") 
        channel.copy_attributes(attributes)
        del sites

        self.overridden().recommit_particles()

    def evolve_model(self,tend):
        self.overridden().evolve_model(tend)
        self.simplex_to_gas_channel.copy_attributes(["xion","u","metallicity"])
        
    def define_state(self, object):
        CommonCode.define_state(self, object)
        
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        object.add_method('UPDATE','before_get_parameter')
        object.add_method('EVOLVED','before_get_parameter')
        object.add_method('RUNCOMMIT','before_get_parameter')
        
        
        
        object.add_method('EDIT', 'new_particle')
        object.add_method('EDIT', 'delete_particle')
        object.add_transition('EDIT', 'RUNCOMMIT', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_particle', False)
        object.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        object.add_transition('RUN','RUNCOMMIT', 'recommit_particles')
        object.add_transition('RUNCOMMIT','RUN', 'auto_go_to_run')
        object.add_transition('RUNCOMMIT', 'EVOLVED', 'evolve_model', False)
#        object.add_method('EVOLVED', 'evolve_model')
        object.define_state('RUNCOMMIT')
        object.add_transition('EVOLVED','RUN', 'synchronize_model')
        object.add_method('RUN', 'synchronize_model')
        object.add_method('RUN', 'get_state')
        object.add_method('RUN', 'get_density')
        object.add_method('RUN', 'get_position')
        object.add_method('RUN', 'get_flux')
        object.add_method('RUN', 'get_ionisation')
        object.add_method('RUN', 'get_internal_energy')
        object.add_method('RUN', 'set_dinternal_energy_dt')
        object.add_method('RUN', 'get_dinternal_energy_dt')
        object.add_method('UPDATE', 'set_dinternal_energy_dt')
        object.add_method('UPDATE', 'get_dinternal_energy_dt')

        object.add_method('INITIALIZED', 'set_hilbert_order')
        object.add_method('INITIALIZED', 'set_box_size')
        object.add_method('INITIALIZED', 'set_timestep')
        object.add_method('INITIALIZED', 'set_source_Teff')
        object.add_method('INITIALIZED', 'set_number_frequency_bins')
        object.add_method('INITIALIZED', 'set_thermal_evolution')
        object.add_method('INITIALIZED', 'set_blackbody_spectrum')
        object.add_method('INITIALIZED', 'set_metal_cooling')
        object.add_method('INITIALIZED', 'set_recombination_radiation')
        object.add_method('INITIALIZED', 'set_collisional_ionization')
