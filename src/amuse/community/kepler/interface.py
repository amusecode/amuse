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
        function.addParameter('time', dtype='float64', direction=function.IN,
                              default = 0)
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
        function.addParameter('time', dtype='float64', direction=function.IN,
                              default = 0)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            new kepler was created
        -1 - ERROR
            kepler could not be created"""
        return function
    
    @legacy_function
    def transform_to_time():
        """
        Transform the kepler system to the specified time.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            transform to time OK
        -1 - ERROR
            could not transform to time"""
        return function

    @legacy_function
    def advance_to_radius():
        """
        Evolve the kepler system forward in time to the specified radius.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('radius', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            advance to radius OK
        -1 - ERROR
            could not advance to radius"""
        return function

    @legacy_function
    def return_to_radius():
        """
        Evolve the kepler system backward in time to the specified radius.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('radius', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            return to radius OK
        -1 - ERROR
            could not return to radius"""
        return function

    @legacy_function
    def advance_to_periastron():
        """
        Evolve the kepler system forward in time to the next periastron.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            advance to periastron OK
        -1 - ERROR
            could not advance to periastron"""
        return function

    @legacy_function
    def advance_to_apastron():
        """
        Evolve the kepler system forward in time to the next apastron.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            advance to apastron OK
        -1 - ERROR
            could not advance to apastron"""
        return function

    @legacy_function
    def return_to_periastron():
        """
        Evolve the kepler system backward in time to the previous periastron.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            return to periastron OK
        -1 - ERROR
            could not return to periastron"""
        return function

    @legacy_function
    def return_to_apastron():
        """
        Evolve the kepler system backward in time to the previous apastron.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            return to apastron OK
        -1 - ERROR
            could not return to apastron"""
        return function

    @legacy_function
    def get_elements():
        """
        Return the orbital elements (a,e,M) of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('semi', dtype='float64', direction=function.OUT)
        function.addParameter('ecc', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get elements OK
        -1 - ERROR
            could not get elements"""
        return function

    @legacy_function
    def get_separation():
        """
        Return the current separation vector (x,y,z) of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get separation OK
        -1 - ERROR
            could not get separation"""
        return function

class kepler(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,
                               keplerInterface(**options),
                               **options)

    def define_methods(self, object):
        CommonCode.define_methods(self, object)

        # Turn interface functions into methods.

        object.add_method(
            "initialize_from_dyn",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.time
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
            "initialize_from_elements",
            (
                nbody_system.mass,
                nbody_system.length,
                units.none,
                nbody_system.time
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method("transform_to_time",
                          (nbody_system.time),
                          (object.ERROR_CODE))
        object.add_method("advance_to_radius",
                          (nbody_system.length),
                          (object.ERROR_CODE))
        object.add_method("return_to_radius",
                          (nbody_system.length),
                          (object.ERROR_CODE))

        object.add_method("advance_to_periastron",
                          (),
                          (object.ERROR_CODE))
        object.add_method("advance_to_apastron",
                          (),
                          (object.ERROR_CODE))
        object.add_method("return_to_periastron",
                          (),
                          (object.ERROR_CODE))
        object.add_method("return_to_apastron",
                          (),
                          (object.ERROR_CODE))

        object.add_method("get_elements",
                          (),
                          (
                              nbody_system.length,
                              units.none,
                              object.ERROR_CODE
                          ))

        object.add_method("get_separation",
                          (),
                          (
                              nbody_system.length,
                              nbody_system.length,
                              nbody_system.length,
                              object.ERROR_CODE
                          ))
