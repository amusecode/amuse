"""
Hydrodynamics Interface Defintion
"""

from amuse.support.codes.core import legacy_function, LegacyFunctionSpecification
from amuse.support.interface import CodeInterface
from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_converter
from amuse.community.interface import common

class HydrodynamicsInterface(common.CommonCodeInterface):

    
    def get_index_range_inclusive(self):
        """
        Returns the min and max values of indices in each
        direction. The range is inclusive, the min index
        and max index both exist and can be queried.
        The total number of cells in one direction
        is max - min + 1.
        
        The returntype is a list of 3 tuples, each tuples
        contains the minimum and maximum value in the
        index range.
        
        Fo C/C++ codes the returned values will usually
        be: ((0, nmeshx-1), (0, nmeshy-1), (0, nmeshz-1))
        
        
        Fo Fortran codes the returned values will usually
        be: ((1, nmeshx), (1, nmeshy), (1, nmeshz))
        """
        pass
        
    
    @legacy_function    
    def get_position_of_index():
        """
        Retrieves the x, y and z position of the center of
        the cell with coordinates i, j, k in the grid specified
        by level and domain.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['level','domain']:
            function.addParameter(x, dtype='i', direction=function.IN, default = 0)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_index_of_position():
        """
        Retrieves the i,j and k index of the grid cell containing the
        given x, y and z position. The cell is looked up
        in the grid specified by level and domain.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['level','domain']:
            function.addParameter(x, dtype='i', direction=function.IN, default = 0)
        for x in ['i','j','k']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
     
     
    
    def setup_mesh(self, nmeshx, nmeshy, nmeshz, xlength, ylength, zlength):
        """
        Sets the number of mesh cells in each direction, and the
        length of the grid in each direction.
        
        :argument nmeshx: number of mesh cells in the x direction
        :argument nmeshy: number of mesh cells in the y direction
        :argument nmeshz: number of mesh cells in the z direction
        :argument xlength: total length of the grid in the x direction
        :argument ylength: total length of the grid in the y direction
        :argument zlength: total length of the grid in the z direction
        """
        
        
    def set_boundary(self, xbound1, xbound2, ybound1, ybound2, zbound1, zbound2):
        """
        Sets the boundary conditions on the grid. Boundaries can be:
        "reflective", "periodic".
        
        :argument xbound1: inner or left boundary in the x direction
        :argument xbound2: outer or right boundary in the x direction
        :argument ybound1: inner or front boundary in the y direction
        :argument ybound2: outer or back boundary in the y direction
        :argument zbound1: inner or bottom boundary in the z direction
        :argument zbound1: outer or top boundary in the z direction
        """
    
        
    @legacy_function    
    def initialize_grid():
        """
        Perform accounting before evolving the model. This method will
        be called after setting the parameters and filling the grid points but
        just before evolving the system
        """
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_timestep():
        """
        Returns the timestep taken by the code.
        """
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.OUT) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_timestep():
        """
        Sets the timestep taken by the code (not implemented by
        all codes)
        """
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.IN) 
        function.result_type = 'i'
        return function
        
    
    @legacy_function
    def set_has_external_gravitational_potential():
        """
        When True enables the script to set and external gravitational
        potential. 
        
        note::
            Not every hydrodynamics code supports an external gravitational
            potential
        """
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='int32', direction=function.IN) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_has_external_gravitational_potential():
        """
        Returns true if an external gravitational potential is enabled.
        """
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='int32', direction=function.OUT) 
        function.result_type = 'i'
        return function
    
        
    @legacy_function
    def get_time():
        """
        Returns the current model time.
        """
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.OUT) 
        function.result_type = 'i'
        return function
        
        
        
    @legacy_function
    def evolve():
        """
        Evolve the model until the given end time (or just before).
        """
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.IN) 
        function.result_type = 'i'
        return function
        
        
    def get_index_range_for_potential(self):
        """
        Returns the min and max values of indices in each
        direction for the potential field, this
        range is 1 cell larger than the normal grid
        in all directions"""
        pass  
    
    @legacy_function    
    def set_potential():
        """
        Sets the gravitational potential on the given gridpoint (only
        for codes supporting an external gravitational potential).
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('potential', dtype='d', direction=function.IN)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_potential():
        """
        Retrieves the gravitational potential on the given gridpoint (only
        for codes supporting an external gravitational potential).
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('potential', dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_interpolated_gravitational_potential():
        """
        Return the interpolated gravitational potential.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('potential', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    

    @legacy_function
    def get_density():
        """
        Retreives the densitity at the given grid-point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['level','domain']:
            function.addParameter(x, dtype='i', direction=function.IN, default = 0)
        for x in ['rho',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    

    @legacy_function
    def get_energy_density():
        """
        Retreives the energy densitity at the given grid-point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['level','domain']:
            function.addParameter(x, dtype='i', direction=function.IN, default = 0)
        for x in ['en',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    

    @legacy_function
    def get_momentum_density():
        """
        Retreives the momentum densitity at the given grid-point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['level','domain']:
            function.addParameter(x, dtype='i', direction=function.IN, default = 0)
        for x in ['rhovx', 'rhovy', 'rhovz',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    

    @legacy_function
    def get_grid_state():
        """
        Retreives the densitity, momentum denisity and energy
        density at the given grid-point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['level','domain']:
            function.addParameter(x, dtype='i', direction=function.IN, default = 0)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        
        return function
    
    

    @legacy_function
    def set_grid_state():
        """
        Sets the densitity, momentum denisity and energy
        density at the given grid-point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['level','domain']:
            function.addParameter(x, dtype='i', direction=function.IN, default = 0)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    
