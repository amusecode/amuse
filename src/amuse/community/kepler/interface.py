from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.options import option
from amuse.support.units import units
import os.path

class keplerInterface(CodeInterface, CommonCodeInterface):
    """
    Kepler structure functions, imported from Starlab.  Initialize an
    orbit from mass, pos, and vel, or mass, semi-major axis and
    eccentricity, and allow the user to manipulate the resulting
    structure.
    """

    # Interface specification.

    include_headers = ['interface.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(self,
                               name_of_the_worker = "kepler_worker",
                               **options)
    
    @legacy_function
    def initialize_from_dyn():
        """
        Initialize a new kepler system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.IN)
        function.addParameter('vy', dtype='float64', direction=function.IN)
        function.addParameter('vz', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            new kepler was created
        -1 - ERROR
            kepler could not be created"""
        return function
    
    @legacy_function
    def initialize_from_elements():
        """
        Initialize a new kepler system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.addParameter('semi', dtype='float64', direction=function.IN)
        function.addParameter('ecc', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            new kepler was created
        -1 - ERROR
            kepler could not be created"""
        return function
    
    @legacy_function
    def evolve_to_time():
        """
        Evolve the kepler system to the specified time.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            evolve to time OK
        -1 - ERROR
            could not evolve to time"""
        return function

    @legacy_function
    def advance_to_radius():
        """
        Evolve the kepler system to the specified time.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            evolve to time OK
        -1 - ERROR
            could not evolve to time"""
        return function

class kepler(CommonCode):

    # The actual module.

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,
                                               keplerInterface(**options),
                                               **options)
