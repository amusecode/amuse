from amuse.legacy import *
from amuse.legacy.interface.common import CommonCodeInterface
from amuse.support.core import OrderedDictionary

import tempfile
import os
import time

class MocassinInterface(LegacyInterface, CommonCodeInterface):
    MOCASSIN_VERSION = '2.02.66'
    use_modules = ['mocassin_interface',]
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, name_of_the_worker="mocassin_worker", **keyword_arguments)
    
    def get_default_input_directory(self):
        return (os.path.join(os.path.join(os.path.dirname(__file__), 'src'), 'mocassin.{0}'.format(self.MOCASSIN_VERSION))) + os.sep
        
    def setup_abundancies(self):
        fd, name = tempfile.mkstemp()
        print name
        if len(name) > 1024:
            raise Exception("Error in filename length of abundancies file, maximum length is 1024")
        with os.fdopen(fd, 'w') as f:
            for atomname, value in self.abundancies_table().iteritems():
                f.write("{0} !{1}\n".format(value, atomname))
        self.set_abundancies_filename(1,name)
            
        
    
    def abundancies_table(self):
        result = OrderedDictionary()
        result['H'] = 1.     
        result['He'] = 0.1    
        result['Li'] = 0.     
        result['Be'] = 0.     
        result['B'] = 0.     
        result['C'] = 2.2e-4 
        result['N'] = 4.e-5  
        result['O'] = 3.3e-4 
        result['F'] = 0.             
        result['Ne'] = 5.e-5     
        result['Na'] = 0.             
        result['Mg'] = 0.     
        result['Al'] = 0.        
        result['Si'] = 0.     
        result['P'] = 0.        
        result['S'] = 9.e-6     
        result['Cl'] = 0.        
        result['Ar'] = 0.     
        result['K'] = 0.        
        result['Ca'] = 0.        
        result['Sc'] = 0.        
        result['Ti'] = 0.        
        result['V'] = 0.        
        result['Cr'] = 0.        
        result['Mn'] = 0.        
        result['Fe'] = 0.     
        result['Co'] = 0.        
        result['Ni'] = 0.        
        result['Cu'] = 0.
        result['Zn'] = 0.
        return result

    @legacy_function
    def setup_mesh():
        function = LegacyFunctionSpecification() 
        function.addParameter('nmeshx', dtype='int32', direction=function.IN) 
        function.addParameter('nmeshy', dtype='int32', direction=function.IN) 
        function.addParameter('nmeshz', dtype='int32', direction=function.IN) 
        function.addParameter('xlength', dtype='float64', direction=function.IN) 
        function.addParameter('ylength', dtype='float64', direction=function.IN) 
        function.addParameter('zlength', dtype='float64', direction=function.IN) 
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1) 
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def setup_auto_convergence():
        function = LegacyFunctionSpecification() 
        function.addParameter('convergence_level_increase', dtype='float64', direction=function.IN) 
        function.addParameter('number_of_photons_increase', dtype='float64', direction=function.IN) 
        function.addParameter('maximum_number_of_photons', dtype='int32', direction=function.IN) 
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def has_auto_convergence():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='bool', direction=function.OUT) 
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def uset_auto_convergence():
        function = LegacyFunctionSpecification() 
        function.result_type = 'int32'
        return function
    
    @legacy_function    
    def get_position_of_index():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='int32', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function    
    def set_abundancies_filename():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('filename', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_abundancies_filename():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('filename', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def set_input_directory():
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_input_directory():
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function
        

    @legacy_function    
    def get_constant_hydrogen_density():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function    
    def set_constant_hydrogen_density():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function    
    def has_constant_hydrogen_density():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    
    @legacy_function    
    def redirect_outputs_to():
        function = LegacyFunctionSpecification()  
        function.addParameter('stdoutstring', dtype='s', direction=function.IN)
        function.addParameter('stderrstring', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_total_number_of_photons():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_total_number_of_photons():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_total_number_of_points_in_frequency_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_total_number_of_points_in_frequency_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
        
    @legacy_function
    def set_initial_nebular_temperature():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_initial_nebular_temperature():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def commit_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def commit_grid():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def iterate():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_symmetricXYZ():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_symmetricXYZ():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_maximum_number_of_monte_carlo_iterations():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_maximum_number_of_monte_carlo_iterations():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_minimum_convergence_level():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_minimum_convergence_level():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_high_limit_of_the_frequency_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_high_limit_of_the_frequency_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_low_limit_of_the_frequency_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_low_limit_of_the_frequency_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_inner_radius_of_the_ionised_region():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_inner_radius_of_the_ionised_region():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_outer_radius_of_the_ionised_region():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_outer_radius_of_the_ionised_region():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_convergence_limit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_convergence_limit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_number_of_ionisation_stages():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_ionisation_stages():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def define_stars():
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('temperature', dtype='float64', direction=function.IN)
        function.addParameter('luminosity', dtype='float64', direction=function.IN)
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.must_handle_array=True
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_write_snapshot_every_iteration():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_write_snapshot_every_iteration():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_number_of_elements_used():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_grid_electron_temperature():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
        
        function.addParameter('electron_temperature', dtype='float64', direction=function.OUT)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_grid_electron_density():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
        
        function.addParameter('electron_density', dtype='float64', direction=function.OUT)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_grid_active():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
        
        function.addParameter('is_active', dtype='bool', direction=function.OUT)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_emit_rate_of_photons():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_emit_rate_of_photons():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
class Mocassin(CodeInterface):

    def __init__(self):
        CodeInterface.__init__(self,  MocassinInterface())
    

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_get_abundancies_filename",
            "set_get_abundancies_filename", 
            "get_abundancies_filename", 
            "<fill>", 
            units.none, 
        )
    
    
        object.add_method_parameter(
            "get_get_constant_hydrogen_density",
            "set_get_constant_hydrogen_density", 
            "get_constant_hydrogen_density", 
            "<fill>", 
            units.none, 
            100.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_convergence_limit",
            "set_get_convergence_limit", 
            "get_convergence_limit", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_emit_rate_of_photons",
            "set_get_emit_rate_of_photons", 
            "get_emit_rate_of_photons", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_high_limit_of_the_frequency_mesh",
            "set_get_high_limit_of_the_frequency_mesh", 
            "get_high_limit_of_the_frequency_mesh", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_initial_nebular_temperature",
            "set_get_initial_nebular_temperature", 
            "get_initial_nebular_temperature", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_inner_radius_of_the_ionised_region",
            "set_get_inner_radius_of_the_ionised_region", 
            "get_inner_radius_of_the_ionised_region", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_input_directory",
            "set_get_input_directory", 
            "get_input_directory", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_low_limit_of_the_frequency_mesh",
            "set_get_low_limit_of_the_frequency_mesh", 
            "get_low_limit_of_the_frequency_mesh", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_maximum_number_of_monte_carlo_iterations",
            "set_get_maximum_number_of_monte_carlo_iterations", 
            "get_maximum_number_of_monte_carlo_iterations", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_minimum_convergence_level",
            "set_get_minimum_convergence_level", 
            "get_minimum_convergence_level", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_number_of_ionisation_stages",
            "set_get_number_of_ionisation_stages", 
            "get_number_of_ionisation_stages", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_outer_radius_of_the_ionised_region",
            "set_get_outer_radius_of_the_ionised_region", 
            "get_outer_radius_of_the_ionised_region", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_symmetricXYZ",
            "set_get_symmetricXYZ", 
            "get_symmetricXYZ", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_total_number_of_photons",
            "set_get_total_number_of_photons", 
            "get_total_number_of_photons", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_total_number_of_points_in_frequency_mesh",
            "set_get_total_number_of_points_in_frequency_mesh", 
            "get_total_number_of_points_in_frequency_mesh", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_get_write_snapshot_every_iteration",
            "set_get_write_snapshot_every_iteration", 
            "get_write_snapshot_every_iteration", 
            "<fill>", 
            units.none, 
            0.0 | units.none
        )
    
    
