"""
Magnetohydrodynamics Interface Defintion
"""

from amuse.community.interface.hydro import HydrodynamicsInterface
from amuse.support.interface import InCodeComponentImplementation
from amuse.units import nbody_system
from amuse.units import generic_unit_converter
from amuse.community.interface import common

from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification
class MagnetohydrodynamicsInterface(HydrodynamicsInterface):
    
        
    @legacy_function
    def get_grid_magnetic_field():
        """
        Retreives the densitity at the given grid-point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['B1i', 'B2i', 'B3i']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    
    @legacy_function
    def set_grid_magnetic_field():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['B1i', 'B2i', 'B3i']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        

