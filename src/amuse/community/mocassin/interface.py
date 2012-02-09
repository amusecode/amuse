from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface
from amuse.support.core import OrderedDictionary
from amuse.units import derivedsi
from amuse.support.options import option

import tempfile
import os
import time

class MocassinInterface(CodeInterface, CommonCodeInterface):
    MOCASSIN_VERSION = '2.02.69'
    use_modules = ['mocassin_interface',]
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="mocassin_worker", **keyword_arguments)
        self._abundancies_table = None
        
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
        
    def get_default_input_directory(self):
        return (os.path.join(os.path.join(os.path.dirname(__file__), 'src'), 'mocassin.{0}'.format(self.MOCASSIN_VERSION))) + os.sep
        
    def get_default_output_directory(self):
        return os.path.join(self.input_data_root_directory, 'mocassin', 'output')+ os.sep
        
    def setup_abundancies(self):
        fd, name = tempfile.mkstemp()
        #print name
        if len(name) > 1024:
            raise Exception("Error in filename length of abundancies file, maximum length is 1024")
        with os.fdopen(fd, 'w') as f:
            for atomname, value in self.abundancies_table().iteritems():
                f.write("{0} !{1}\n".format(value, atomname))
        self.set_abundancies_filename(name, 1)
            
        
    
    def abundancies_table(self):
        if self._abundancies_table  is None:
            self._abundancies_table = self.default_abundancies_table()
        return self._abundancies_table
        
    def default_abundancies_table(self):
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
    def set_random_seed():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='int32', direction=function.IN) 
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
        function.addParameter('filename', dtype='s', direction=function.IN)
        function.addParameter('index', dtype='int32', direction=function.IN, default=1)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_abundancies_filename():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN, default=1)
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
    def set_output_directory():
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='s', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_output_directory():
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='s', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function    
    def get_percentage_converged():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
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
    def set_has_constant_hydrogen_density():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
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
    def step():
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
    def get_max_indices():
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
        for parametername in ['ni','nj','nk']:
            function.addParameter(parametername, dtype='int32', direction=function.OUT)
        function.can_handle_array = True
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
        
    @legacy_function
    def get_grid_hydrogen_density():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
        
        function.addParameter('hydrogen_density', dtype='float64', direction=function.OUT)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_grid_hydrogen_density():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        
        function.addParameter('hydrogen_density', dtype='float64', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
            
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
        
        function.addParameter('density', dtype='float64', direction=function.OUT)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_grid_electron_density():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        
        function.addParameter('density', dtype='float64', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function
        
    
    @legacy_function
    def get_grid_ion_density():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k','elemen', 'bin']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
        
        function.addParameter('value', dtype='float64', direction=function.OUT)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_grid_ion_density():
        function = LegacyFunctionSpecification()  
        for parametername in ['i','j','k','elemen', 'bin']:
            function.addParameter(parametername, dtype='int32', direction=function.IN)
        
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN, default = 1)
            
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.must_handle_array = True
        function.result_type = 'int32'
        return function
        
    
        
#class Parameters(object):
#    emit_rate_of_photons = MethodParameter(name = "emit_rate_of_photons", unit = units.seconds, default = )
#
#class Methods(object):
#    get_grid_electron_density = MethodWithUnit(name = "get_grid_electron_density", i = INDEX, j = INDEX, k = INDEX, is_active = NOUNIT )

units.ryd = 13.6 * derivedsi.eV

class Mocassin(InCodeComponentImplementation):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  MocassinInterface(**options), **options)
        ensure_data_directory_exists(self.get_default_output_directory())
        
    def get_index_range_inclusive(self, index_of_grid = 1):
        ni, nj, nk = self.get_max_indices(index_of_grid)
        return (1, ni, 1, nj, 1, nk)
        
    def get_index_range_inclusive_ion_density_grid(self, index_of_grid = 1):
        ni, nj, nk = self.get_max_indices(index_of_grid)
        nstages = self.get_number_of_ionisation_stages()
        nbins = self.get_total_number_of_points_in_frequency_mesh()
        return (1, ni, 1, nj, 1, nk, 1, nstages, 1, nbins)
        
    def define_methods(self, object):
        
        object.add_method(
            'commit_grid',
            (),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'commit_particles',
            (),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'iterate',
            (),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'step',
            (),
            (object.ERROR_CODE,)
        )
        
        
        object.add_method(
            'get_percentage_converged',
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        object.add_method(
            'get_position_of_index',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (units.cm, units.cm, units.cm, object.ERROR_CODE,)
        )
        object.add_method(
            'get_max_indices',
            (object.INDEX),
            (object.NO_UNIT, object.NO_UNIT, object.NO_UNIT, object.ERROR_CODE,)
            
        )
        object.add_method(
            'get_grid_electron_temperature',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (units.K, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_grid_hydrogen_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (units.cm**-3, object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_grid_hydrogen_density',
            (object.INDEX, object.INDEX, object.INDEX, units.cm**-3 , object.INDEX),
            (object.ERROR_CODE,)
        )
        
                
        object.add_method(
            'get_grid_electron_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (units.cm**-3, object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_grid_electron_density',
            (object.INDEX, object.INDEX, object.INDEX, units.cm**-3 , object.INDEX),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_ion_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_grid_ion_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX, units.none , object.INDEX),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_grid_active',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            'define_stars',
            (units.cm, units.cm, units.cm, units.K, 1e36 * units.erg * (units.s**-1)),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_abundancies_filename",
            (object.NO_UNIT,  ),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_abundancies_filename",
            (object.NO_UNIT, object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_constant_hydrogen_density",
            (),
            (1.0/units.cm**3, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_constant_hydrogen_density",
            (1.0/units.cm**3, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_convergence_limit",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_convergence_limit",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_emit_rate_of_photons",
            (),
            (1e36 / units.s, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_emit_rate_of_photons",
            (1e36 / units.s, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_high_limit_of_the_frequency_mesh",
            (),
            (units.ryd, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_high_limit_of_the_frequency_mesh",
            (units.ryd, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_initial_nebular_temperature",
            (),
            (units.K, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_initial_nebular_temperature",
            (units.K, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_inner_radius_of_the_ionised_region",
            (),
            (units.cm, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_inner_radius_of_the_ionised_region",
            (units.cm, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_input_directory",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_input_directory",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_output_directory",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_output_directory",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_low_limit_of_the_frequency_mesh",
            (),
            (units.ryd, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_low_limit_of_the_frequency_mesh",
            (units.ryd, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_maximum_number_of_monte_carlo_iterations",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_maximum_number_of_monte_carlo_iterations",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_minimum_convergence_level",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_minimum_convergence_level",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_number_of_ionisation_stages",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_number_of_ionisation_stages",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        
        object.add_method(
            "get_symmetricXYZ",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_symmetricXYZ",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_total_number_of_photons",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_total_number_of_photons",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_total_number_of_points_in_frequency_mesh",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_total_number_of_points_in_frequency_mesh",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_write_snapshot_every_iteration",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_write_snapshot_every_iteration",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'setup_mesh',
            (units.none,  units.none, units.none, units.cm, units.cm, units.cm, units.none,),
            (object.ERROR_CODE,)
        )
        
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_abundancies_filename",
            "set_abundancies_filename", 
            "abundancies_filename", 
            "<fill>", 
            default_value = ""
        )
    
    
        object.add_method_parameter(
            "get_constant_hydrogen_density",
            "set_constant_hydrogen_density", 
            "constant_hydrogen_density", 
            "<fill>", 
            default_value = 100.0 | (1.0/units.cm**3)
        )
    
    
        object.add_method_parameter(
            "get_convergence_limit",
            "set_convergence_limit", 
            "convergence_limit", 
            "<fill>", 
            default_value = 0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_emit_rate_of_photons",
            "set_emit_rate_of_photons", 
            "emit_rate_of_photons", 
            "This is the number of hydrogen-ionizing photons emitted by the source per unit time. Only used when a single star is modelled", 
            default_value = 0.0 | (1e36 / units.s)
        )
    
    
        object.add_method_parameter(
            "get_high_limit_of_the_frequency_mesh",
            "set_high_limit_of_the_frequency_mesh", 
            "high_limit_of_the_frequency_mesh", 
            "<fill>", 
            default_value = 15. | units.ryd
        )
    
    
        object.add_method_parameter(
            "get_initial_nebular_temperature",
            "set_initial_nebular_temperature", 
            "initial_nebular_temperature", 
            "Initial guess for the nebular temperature. ", 
            default_value = 10000.0 | units.K
        )
    
    
        object.add_method_parameter(
            "get_inner_radius_of_the_ionised_region",
            "set_inner_radius_of_the_ionised_region", 
            "inner_radius_of_the_ionised_region", 
            "Inner radius of the ionised region", 
            default_value = 0.0 | units.cm
        )
    
    
        object.add_method_parameter(
            "get_input_directory",
            "set_input_directory", 
            "input_directory", 
            "<fill>", 
            default_value = 0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_low_limit_of_the_frequency_mesh",
            "set_low_limit_of_the_frequency_mesh", 
            "low_limit_of_the_frequency_mesh", 
            "<fill>", 
            default_value = 1.001e-5 | units.ryd
        )
    
    
        object.add_method_parameter(
            "get_maximum_number_of_monte_carlo_iterations",
            "set_maximum_number_of_monte_carlo_iterations", 
            "maximum_number_of_monte_carlo_iterations", 
            "<fill>", 
            default_value = 0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_minimum_convergence_level",
            "set_minimum_convergence_level", 
            "minimum_convergence_level", 
            "<fill>", 
            default_value = 0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_number_of_ionisation_stages",
            "set_number_of_ionisation_stages", 
            "number_of_ionisation_stages", 
            "<fill>", 
            default_value = 6 | units.none
        )
    
    
        
    
        object.add_method_parameter(
            "get_symmetricXYZ",
            "set_symmetricXYZ", 
            "symmetricXYZ", 
            "If true assumes model is symetric in the X, Y and Z axes", 
            default_value = 0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_total_number_of_photons",
            "set_total_number_of_photons", 
            "total_number_of_photons", 
            "Total number of photons to start the iteration with", 
            default_value = 0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_total_number_of_points_in_frequency_mesh",
            "set_total_number_of_points_in_frequency_mesh", 
            "total_number_of_points_in_frequency_mesh", 
            "<fill>", 
            default_value = 0.0 | units.none
        )
    
    
        object.add_method_parameter(
            "get_write_snapshot_every_iteration",
            "set_write_snapshot_every_iteration", 
            "write_snapshot_every_iteration", 
            "If True will write the data to an output directory after every monte carlo iteration", 
            default_value = 0.0 | units.none
        )
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshx",
            "nx", 
            "number of cells in the x direction", 
            10 | units.none,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshy",
            "ny", 
            "number of cells in the y direction", 
            10 | units.none,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshz",
            "nz", 
            "number of cells in the z direction", 
            10 | units.none,
        )
        
        object.add_caching_parameter(
            "setup_mesh", 
            "xlength",
            "length_x", 
            "length of model in the x direction", 
            2e19 | units.cm,
        )
        object.add_caching_parameter(
            "setup_mesh", 
            "ylength",
            "length_y", 
            "length of model in the x direction", 
            2e19 | units.cm,
        )
        object.add_caching_parameter(
            "setup_mesh", 
            "zlength",
            "length_z", 
            "length of model in the z direction", 
            2e19 | units.cm,
        )
        
        object.add_vector_parameter(
            "mesh_size",
            "number of cells in the x, y and z directions",
            ("nx", "ny", "nz")
        )
        
        object.add_vector_parameter(
            "mesh_length",
            "length of the model in the x, y and z directions",
            ("length_x", "length_y", "length_z")
        )
        
    def commit_parameters(self):
        self.setup_abundancies()
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
        
    def define_particle_sets(self, object):
        object.define_grid('grid')
        object.set_grid_range('grid', 'get_index_range_inclusive')
        object.add_getter('grid', 'get_position_of_index', names=('x','y','z'))
        object.add_getter('grid', 'get_grid_electron_temperature', names=('electron_temperature',))
        object.add_getter('grid', 'get_grid_electron_density', names=('electron_density',))
        object.add_setter('grid', 'set_grid_electron_density', names=('electron_density',))
        object.add_getter('grid', 'get_grid_hydrogen_density', names=('hydrogen_density',))
        object.add_setter('grid', 'set_grid_hydrogen_density', names=('hydrogen_density',))
        
        object.define_extra_keywords('grid', {'index_of_grid':1})
        
        
        object.define_grid('ion_density_grid')
        object.set_grid_range('ion_density_grid', 'get_index_range_inclusive_ion_density_grid')
        object.add_getter('ion_density_grid', 'get_grid_ion_density', names=('density',))
        object.add_setter('ion_density_grid', 'set_grid_ion_density', names=('density',))
        object.define_extra_keywords('ion_density_grid', {'index_of_grid':1})
        
        object.define_inmemory_set('particles')
        
    def commit_particles(self):
        self.define_stars(
            self.particles.x,
            self.particles.y,
            self.particles.z,
            self.particles.temperature,
            self.particles.luminosity
        )
        self.overridden().commit_particles()
        
    
