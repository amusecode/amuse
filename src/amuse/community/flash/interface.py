from amuse.community import *
from amuse.community.interface.hydro import HydrodynamicsInterface
from amuse.support.options import OptionalAttributes, option
from amuse.units import generic_unit_system
from amuse.units import si
from amuse.community.interface.common import CommonCode
import numpy as np
import os

length = generic_unit_system.length
time = generic_unit_system.time
mass = generic_unit_system.mass
speed = generic_unit_system.speed
density = generic_unit_system.density
momentum =  generic_unit_system.momentum_density
energy =  generic_unit_system.energy_density
potential_energy =  generic_unit_system.energy
magnetic_field = generic_unit_system.mass / generic_unit_system.current / generic_unit_system.time ** 2
acc = generic_unit_system.acceleration
potential = generic_unit_system.potential

class FlashInterface(CodeInterface, HydrodynamicsInterface):

    use_modules = ['flash_run']

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="flash_worker", **keyword_arguments)
    

    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function

    @legacy_function
    def cleanup_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_grid_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid', 'nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx', 'rhovy', 'rhovz', 'rhoen']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid', 'nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx', 'rhovy', 'rhovz', 'rhoen']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_momentum_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid', 'nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rhovx', 'rhovy', 'rhovz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_momentum_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rhovx', 'rhovy', 'rhovz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function
        
    @legacy_function
    def get_grid_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_energy_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rhoen', dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_energy_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rhoen', dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rho', dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rho', dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function
        
    @legacy_function
    def get_potential():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k', 'index_of_grid']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('potential', dtype='d', direction=function.OUT)
        function.result_type='int32'
        return function
        
    @legacy_function
    def get_potential_at_point():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['eps','x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('gpot', dtype='d', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function
        
    @legacy_function
    def get_gravity_at_point():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['eps','x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['gax','gay','gaz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_number_of_grids():
        function = LegacyFunctionSpecification()
        function.addParameter('nproc', dtype='int32', direction=function.IN)
        function.addParameter('n', dtype='int32', direction=function.OUT)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_range():
        function = LegacyFunctionSpecification()
        for x in ['nx','ny','nz']:
            function.addParameter(x, dtype='int32', direction=function.OUT)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN)
        function.addParameter('nproc', dtype='int32', direction=function.IN)
        function.result_type='int32'
        return function

    @legacy_function    
    def get_position_of_index():
        """
        Retrieves the x, y and z position of the center of
        the cell with coordinates i, j, k in the grid specified
        by the index_of_grid
        """
        function = LegacyFunctionSpecification()  
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('Proc_ID', dtype='i', direction=function.IN, default = 0)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_index_of_position():
        """
        Retrieves the i,j and k index of the grid cell containing the
        given x, y and z position, the index of the grid and the local
        processor number on which this grid resides.
        """
        function = LegacyFunctionSpecification()  
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.OUT)
        function.addParameter('index_of_grid', dtype='i', direction=function.OUT)
        function.addParameter('Proc_ID', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_leaf_indices():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('dummy', dtype='i', direction=function.IN)
        function.addParameter('ind', dtype='i', direction=function.OUT)
        function.addParameter('ret_cnt', dtype='i', direction=function.OUT)
        function.addParameter('num_of_blks', dtype='i', direction=function.OUT)
        function.addParameter('nparts',dtype='i', direction=function.LENGTH)
        function.result_type='i'
        return function
        
    @legacy_function
    def get_max_refinement():
        function = LegacyFunctionSpecification()
        function.addParameter('max_refine', dtype='int32', direction=function.OUT)
        function.result_type='i'
        return function        

    @legacy_function
    def set_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def get_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def set_end_time():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN, default=0.0)
        function.result_type='i'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def get_end_time():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def set_max_num_steps():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def get_max_num_steps():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def set_begin_iter_step():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN, default=1)
        function.result_type='i'
        return function

    @legacy_function
    def get_begin_iter_step():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type='i'
        return function
        
    @legacy_function
    def initialize_restart():
		function = LegacyFunctionSpecification()
		function.result_type='i'
		return function
		
    @legacy_function
    def get_restart():
		function = LegacyFunctionSpecification()
		function.addParameter('value', dtype='b', direction=function.OUT)
		function.result_type='i'
		return function
		
    @legacy_function
    def set_restart():
		function = LegacyFunctionSpecification()
		function.addParameter('value', dtype='b', direction=function.IN)
		function.result_type='i'
		return function
        
    @legacy_function    
    def get_hydro_state_at_point():
        function = LegacyFunctionSpecification()  
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN, default = 0)
        for x in ['rho','rhovx','rhovy','rhovz','rhoen']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i' 
        return function
    
    # This needs to look like "get_grid_density" etc etc.    
    @legacy_function
    def get_cell_volume():
        function = LegacyFunctionSpecification()
        for x in ['block','i','j','k']:
            function.addParameter(x, dtype='int32', direction=function.IN)
        function.addParameter('vol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_number_of_procs():
        function = LegacyFunctionSpecification()
        function.addParameter('n', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_all_local_num_grids():
        function = LegacyFunctionSpecification()
        function.must_handle_array=True
        function.addParameter('num_grids_array', dtype='i', direction=function.INOUT)
        function.addParameter('nprocs', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function
    
    # WORK IN PROGRESS!!!!    
    #@legacy_function
    #def get_data_all_local_blks():
        #function = LegacyFunctionSpecification()
        #function.must_handle_array = True
        #function.addParameter('data_array', dtype='d', direction=function.INOUT)
        #function.addParameter('numcells', dtype='i', direction=function.LENGTH)
        #function.result_type = 'i'
        #return function
        
    @legacy_function
    def get_1blk_cell_coords():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('axis', dtype='i', direction=function.IN)
        function.addParameter('blockID', dtype='i', direction=function.IN)
        function.addParameter('limits', dtype='i', direction=function.IN)
        function.addParameter('coords', dtype='d', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def kick_grid():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def kick_block():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['accel_x', 'accel_y', 'accel_z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('blockID', dtype='i', direction=function.IN)
        function.addParameter('block_arr', dtype='i', direction=function.IN)
        function.addParameter('limits', dtype='i', direction=function.IN)
        function.addParameter('dt', dtype='d', direction=function.IN)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type='i'
        return function
        
########################################
# Particle Stuff
########################################
    
    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('n', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_accel_gas_on_particles():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('eps', dtype='d', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['ax','ay','az']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_num_part_prop():
        function = LegacyFunctionSpecification()
        function.addParameter('n', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_particle_position_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN,unit=NO_UNIT)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=length)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_particle_position():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=length)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_particle_velocity():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=speed)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_velocity_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=speed)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def make_sink():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('tags', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_particle_mass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_particle_mass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_particle_gpot():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('n', dtype='int32', direction=function.IN)
        function.addParameter('gpot', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_particle_gpot():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('gpot', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function    
        
    @legacy_function
    def get_particle_tags():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('n', dtype='int32', direction=function.IN)
        function.addParameter('tags', dtype='d', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_proc():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('n', dtype='int32', direction=function.IN)
        function.addParameter('procs', dtype='int32', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_particle_block():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('n', dtype='int32', direction=function.IN)
        function.addParameter('blocks', dtype='int32', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_new_tags():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('new_tags_length', dtype = 'i', direction=function.IN)
        function.addParameter('new_tags_array', dtype = 'i', direction=function.OUT)
        function.addParameter('nparts', dtype = 'i', direction=function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_number_of_new_tags():
        function = LegacyFunctionSpecification()
        function.addParameter('new_tag_num', dtype = 'int32', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def clear_new_tags():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function
        
    @legacy_function
    def particles_gather():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function
        
#    @legacy_function
#    def make_particle_tree():
#        function = LegacyFunctionSpecification()
#        function.result_type = 'i'
#        return function
        
        
########################
# IO Stuff
########################

    @legacy_function
    def write_chpt():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

###############################################
# Default implemenation made by build.py - Josh
###############################################

#class Flash(InCodeComponentImplementation):

#    def __init__(self, unit_converter = None, **options):
#        InCodeComponentImplementation.__init__(self,  FlashInterface(**options), **options)

#        object.set_converter(self.unit_converter.as_converter_from_si_to_generic())

#####################################################
# Attempt to copy amrvac's class implemenation, working so far! - Josh
#####################################################

class Flash(CommonCode):
    
        
        
    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        self.stopping_conditions = StoppingConditions(self)
        
        CommonCode.__init__(self,  FlashInterface(**options), **options)
        
#        self.set_parameters_filename(self.default_parameters_filename)
        
    def define_converter(self, object):
        if self.unit_converter is None:
            return
        
        object.set_converter(self.unit_converter.as_converter_from_si_to_generic())


    def get_index_range_inclusive(self, index_of_grid = 1, nproc=0):
        nx, ny, nz = self.get_grid_range(index_of_grid, nproc)
        
        return (0, nx-1, 0, ny-1, 0, nz-1)
        
    def get_particle_position(self, tags):
        
        [x, y, z] = self.get_particle_position_array(tags)
        
        pos_array = np.array([x.value_in(units.cm), y.value_in(units.cm), 
                  z.value_in(units.cm)]).transpose() | units.cm
        
        if (hasattr(x,"__len__") == False):
        
            pos_array = pos_array.flatten()

        return pos_array

    def get_particle_velocity(self, tags):
        
        [x, y, z] = self.get_particle_velocity_array(tags)
        
        vel_array = np.array([x.value_in(units.cm / units.s), y.value_in(units.cm / units.s), 
                  z.value_in(units.cm / units.s)]).transpose() | units.cm / units.s
        
        if (hasattr(x,"__len__") == False):
        
            vel_array = vel_array.flatten()

        return vel_array

    def define_methods(self, object):
        
        #length = units.cm
        #time = units.s
        #mass = units.g
        #speed = units.cm*units.s**-1
        #density = units.g*units.cm**-3
        #momentum =  density*speed
        #energy =  units.cm**2*units.g*units.s**-2
        #potential_energy =  energy
        #magnetic_field = units.g*0.1*units.C**-1*units.s**-1
        #acc = units.cm*units.s**-2

### These two are included in CommonCode

        #object.add_method(
            #'initialize_code',
            #(),
            #(object.ERROR_CODE)
        #)
                
        
        #object.add_method(
            #'cleanup_code',
            #(),
            #(object.ERROR_CODE)
        #)
        
        object.add_method(
            'evolve_model',
            (time,),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_position_of_index',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length, length, length, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_index_of_position',
            (length, length, length),
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX,
             object.INDEX, object.ERROR_CODE,)
        )

        object.add_method(
            "get_max_refinement",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_grid_state',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density, momentum, momentum, momentum, energy,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_state',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            density, momentum, momentum, momentum, energy),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_energy_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            ( energy,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_energy_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
             energy),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            density),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_momentum_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (momentum, momentum, momentum, 
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_grid_momentum_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            momentum, momentum, momentum), 
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_grid_velocity',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (speed, speed, speed, 
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_velocity',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            speed, speed, speed), 
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_potential',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (potential, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_potential_at_point',
            (length, length, length, length),
            (potential, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_gravity_at_point',
            (length, length, length, length),
            (acc, acc, acc, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_cell_volume',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length**3, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_number_of_procs',
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_all_local_num_grids',
            (object.INDEX),
            (object.INDEX, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_data_all_local_blks',
            (object.NO_UNIT),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_1blk_cell_coords',
            (object.NO_UNIT,object.NO_UNIT,object.NO_UNIT),
            (length, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_leaf_indices',
            (object.NO_UNIT),
            (object.NO_UNIT, object.INDEX, object.INDEX, object.ERROR_CODE)
        )
        
        object.add_method(
            'kick_grid',
            (time),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'kick_block',
            (acc, acc, acc, object.INDEX, object.INDEX, object.NO_UNIT, time),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            "get_timestep",
            (),
            (time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_timestep",
            (time),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_end_time",
            (),
            (time, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_end_time",
            (time, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_time',
            (),
            (time, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_max_num_steps',
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_max_num_steps',
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_restart',
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_restart',
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_hydro_state_at_point',
            (length, length, length,
                speed, speed, speed),
            (density, momentum, momentum, 
                momentum, energy, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_number_of_particles',
            (),
            (object.NO_UNIT,object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_accel_gas_on_particles',
            (length,length,length,length),
            (acc,acc,acc, object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_num_part_prop',
            (),
            (object.NO_UNIT,object.ERROR_CODE)
        )
        
        ### I'm implementing this with my own defined function
        ### so that the structure of the array return looks
        ### right.
        
        #object.add_method(
            #'get_particle_position',
            #(object.NO_UNIT),
            #(length, length, length, object.ERROR_CODE) 
        #)
        
        object.add_method(
            'set_particle_position',
            (object.NO_UNIT, length, length, length),
            (object.ERROR_CODE)
        )
        
        ### Same as above.
        
        #object.add_method(
            #'get_particle_velocity',
            #(object.NO_UNIT),
            #(speed, speed, speed, object.ERROR_CODE)
        #)
        
        object.add_method(
            'set_particle_velocity',
            (object.NO_UNIT, speed, speed, speed),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'make_sink',
            (length, length, length),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_mass',
            (object.INDEX),
            (mass, object.ERROR_CODE)
        )
        
        object.add_method(
            'set_particle_mass',
            (object.INDEX, mass),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_gpot',
            (object.INDEX),
            (potential, object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_gpot',
            (object.NO_UNIT, potential),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_tags',
            (object.INDEX),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_proc',
            (object.INDEX),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_block',
            (object.INDEX),
            (object.INDEX, object.ERROR_CODE)
        )
        
        object.add_method(
            'write_chpt',
            (),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'get_new_tags',
            (object.NO_UNIT),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_number_of_new_tags',
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'clear_new_tags',
            (),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'particles_gather',
            (),
            (object.ERROR_CODE)
        )
        
#        object.add_method(
#            'make_particle_tree',
#            (),
#            (object.ERROR_CODE)
#        )

    def specify_grid(self, definition, index_of_grid = 1, nproc=0):
        definition.set_grid_range('get_index_range_inclusive')
        
        definition.add_getter('get_position_of_index', names=('x','y','z'))
        
        definition.add_getter('get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        
        definition.add_getter('get_grid_density', names=('rho'))
        definition.add_setter('set_grid_density', names=('rho'))
        
#       if self.mode == self.MODE_SCALAR:
#           definition.add_getter('get_grid_scalar', names=('scalar',))
#           definition.add_setter('set_grid_scalar', names=('scalar',))
            
        definition.add_getter('get_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        definition.add_setter('set_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        
        #definition.add_getter('get_grid_velocity', names=('vx','vy','vz'))
        #definition.add_setter('set_grid_velocity', names=('vx','vy','vz'))
        
        definition.add_getter('get_grid_energy_density', names=('energy',))
        definition.add_setter('set_grid_energy_density', names=('energy',))
        
        
#       definition.add_getter('get_grid_gravitational_potential', names=('gravitational_potential',))
#       definition.add_getter('get_grid_gravitational_acceleration', names=('gravitational_acceleration_x','gravitational_acceleration_y','gravitational_acceleration_z',))
        
#        definition.add_getter('get_grid_acceleration', names=('ax','ay','az'))
#        definition.add_setter('set_grid_acceleration', names=('ax','ay','az'))
        
        definition.define_extra_keywords({'index_of_grid':index_of_grid,'nproc':nproc})

    @property
    def grid(self):
        return self._create_new_grid(self.specify_grid, index_of_grid = 1, nproc=0)



    # Define an object that returns a list of all the blocks in the simulation.
    # This iterates over all processors and then loops over the blocks on the local processors.
    def itergrids(self):
        m = self.get_number_of_procs()
        
        for x in range(m): # Loop over processors.
            n = self.get_number_of_grids(x)
            #n = max(num_grids)
            #print "N =",n, "X =",x
            for y in range(1, n+1): # Loop over blocks.
                yield self._create_new_grid(self.specify_grid, index_of_grid = y, nproc=x)


    def define_state(self, object): 
        CommonCode.define_state(self, object)   
        object.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        
        object.add_transition('EDIT', 'RUN', 'initialize_grid')
        object.add_method('RUN', 'evolve_model')
        object.add_method('RUN', 'get_hydro_state_at_point')
        
        for state in ['EDIT', 'RUN']:
            for methodname in [
                    'get_grid_state',
                    'set_grid_state',
                    'get_potential_at_point',
                    'get_potential',
                    'get_gravity_at_point',
#                    'set_potential',
                    'get_grid_density',
                    'set_grid_density',
                    'set_grid_energy_density',
                    'get_grid_energy_density',
                    'get_grid_momentum_density',
                    'set_grid_momentum_density',
                    'get_grid_velocity',
                    'set_grid_velocity', 
                    'get_position_of_index',
                    'get_index_of_position',
                    'get_max_refinement',
#                    'set_grid_scalar',
#                    'get_grid_scalar',
                    'get_number_of_grids',
                    'get_index_range_inclusive',
                    'get_cell_volume',
                    'get_number_of_procs',
                    'get_all_local_num_grids',
                    'get_data_all_local_blks',
                    'get_1blk_cell_coords',
                    'get_leaf_indices',
                    'kick_grid',
                    'kick_block',
#                    'get_boundary_state',
#                    'set_boundary_state',
#                    'get_boundary_position_if_index',
#                    'get_boundary_index_range_inclusive',
                    'get_timestep',
                    'set_timestep',
                    'get_end_time',
                    'set_end_time',
                    'get_time',
                    'get_max_num_steps',
                    'set_max_num_steps',
                    'get_restart',
                    'get_number_of_particles',
                    'get_accel_gas_on_particles',
                    'get_particle_position',
                    'set_particle_position',
                    'get_particle_velocity',
                    'set_particle_velocity',
                    'make_sink',
                    'get_particle_mass',
                    'set_particle_mass',
                    'set_particle_gpot',
                    'get_particle_gpot',
                    'get_particle_tags',
                    'get_particle_proc',
                    'get_particle_block',
                    'get_new_tags',
                    'get_number_of_new_tags',
                    'clear_new_tags',
                    'particles_gather',
                    #'make_particle_tree',
                    'write_chpt'
                ]:
                object.add_method(state, methodname)

