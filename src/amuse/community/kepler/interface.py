from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.options import option
from amuse.units import units
import os.path

class KeplerInterface(CodeInterface,
                      CommonCodeInterface):
    """
    Kepler orbit manipulation functions, imported from Starlab.
    Initialize an orbit from mass, pos, and vel, or mass, semi-major
    axis and eccentricity, and allow the user to manipulate the
    resulting structure.  Most Starlab functionality is currently
    exposed.
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
        function.addParameter('mean_anomaly',
                              dtype='float64', direction=function.IN,
                              default = 0)
        function.addParameter('time', dtype='float64', direction=function.IN,
                              default = 0)
        function.addParameter('periastron',
                              dtype='float64', direction=function.IN,
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
    def get_total_mass():
        """
        Return the total mass (remind the user) of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('mass', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get mass OK
        -1 - ERROR
            could not get mass"""
        return function

    @legacy_function
    def get_time():
        """
        Return the current time of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get time OK
        -1 - ERROR
            could not get time"""
        return function

    @legacy_function
    def get_period():
        """
        Return the periodof the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('period', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get period OK
        -1 - ERROR
            could not get period"""
        return function

    @legacy_function
    def get_elements():
        """
        Return the orbital elements (a,e) of the system.
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
    def get_integrals():
        """
        Return the total energy and angular momentum of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('energy', dtype='float64', direction=function.OUT)
        function.addParameter('angular_momentum',
                              dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get integrals OK
        -1 - ERROR
            could not get integrals"""
        return function

    @legacy_function
    def get_separation_vector():
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
            get separation vector OK
        -1 - ERROR
            could not get separation vector"""
        return function

    @legacy_function
    def get_separation():
        """
        Return the current separation r of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('r', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get separation OK
        -1 - ERROR
            could not get separation"""
        return function

    @legacy_function
    def set_periastron():
        """
        Set the current periastron of the system (initialization only).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('peri', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            set periastron OK
        -1 - ERROR
            could not set periastron"""
        return function

    @legacy_function
    def get_periastron():
        """
        Return the current periastron of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('peri', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get periastron OK
        -1 - ERROR
            could not get periastron"""
        return function

    @legacy_function
    def get_apastron():
        """
        Return the current apastron of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('apo', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get apastron OK
        -1 - ERROR
            could not get apastron"""
        return function

    @legacy_function
    def get_velocity_vector():
        """
        Return the current relative velocity vector (x,y,z) of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('vx', dtype='float64', direction=function.OUT)
        function.addParameter('vy', dtype='float64', direction=function.OUT)
        function.addParameter('vz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get velocity vector OK
        -1 - ERROR
            could not get velocity vector"""
        return function

    @legacy_function
    def get_angles():
        """
        Return the current mean and true anomalies of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('true_anomaly',
                              dtype='float64', direction=function.OUT)
        function.addParameter('mean_anomaly',
                              dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get angles OK
        -1 - ERROR
            could not get angles"""
        return function

    @legacy_function
    def set_longitudinal_unit_vector():
        """
        Set the longitudinal unit vector of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('vx', dtype='float64', direction=function.IN)
        function.addParameter('vy', dtype='float64', direction=function.IN)
        function.addParameter('vz', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            set vector OK
        -1 - ERROR
            could not set vector"""
        return function

    @legacy_function
    def set_normal_unit_vector():
        """
        Set the normal unit vector of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('vx', dtype='float64', direction=function.IN)
        function.addParameter('vy', dtype='float64', direction=function.IN)
        function.addParameter('vz', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            set vector OK
        -1 - ERROR
            could not set vector"""
        return function

    @legacy_function
    def get_longitudinal_unit_vector():
        """
        Return the longitudinal unit vector of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('vx', dtype='float64', direction=function.OUT)
        function.addParameter('vy', dtype='float64', direction=function.OUT)
        function.addParameter('vz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get vector OK
        -1 - ERROR
            could not get vector"""
        return function

    @legacy_function
    def get_transverse_unit_vector():
        """
        Return the transverse unit vector of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('vx', dtype='float64', direction=function.OUT)
        function.addParameter('vy', dtype='float64', direction=function.OUT)
        function.addParameter('vz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get vector OK
        -1 - ERROR
            could not get vector"""
        return function

    @legacy_function
    def set_transverse_unit_vector():
        """
        Set the transverse unit vector of the system (tangent on longitudal uv).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('vx', dtype='float64', direction=function.IN)
        function.addParameter('vy', dtype='float64', direction=function.IN)
        function.addParameter('vz', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            set vector OK
        -1 - ERROR
            could not set vector"""
        return function

    @legacy_function
    def get_normal_unit_vector():
        """
        Return the normal unit vector of the system.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('vx', dtype='float64', direction=function.OUT)
        function.addParameter('vy', dtype='float64', direction=function.OUT)
        function.addParameter('vz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            get vector OK
        -1 - ERROR
            could not get vector"""
        return function

class Kepler(CommonCode):

    def __init__(self, unit_converter = None,  **options):
        self.unit_converter = unit_converter
        
        CommonCode.__init__(self,
                               KeplerInterface(**options),
                               **options)

    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
    
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
                object.NO_UNIT,
                object.NO_UNIT,
                nbody_system.time,
                nbody_system.length
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

        object.add_method("get_total_mass",
                          (),
                          (
                              nbody_system.mass,
                              object.ERROR_CODE
                          ))

        object.add_method("get_time",
                          (),
                          (
                              nbody_system.time,
                              object.ERROR_CODE
                          ))

        object.add_method("get_period",
                          (),
                          (
                              nbody_system.time,
                              object.ERROR_CODE
                          ))

        object.add_method("get_elements",
                          (),
                          (
                              nbody_system.length,
                              object.NO_UNIT,
                              object.ERROR_CODE
                          ))

        object.add_method("get_integrals",
                          (),
                          (
                              nbody_system.speed*nbody_system.speed,
                              nbody_system.length*nbody_system.speed,
                              object.ERROR_CODE
                          ))

        object.add_method("get_separation_vector",
                          (),
                          (
                              nbody_system.length,
                              nbody_system.length,
                              nbody_system.length,
                              object.ERROR_CODE
                          ))

        object.add_method("get_separation",
                          (),
                          (
                              nbody_system.length,
                              object.ERROR_CODE
                          ))

        object.add_method("set_periastron",
                          (nbody_system.length,),
                          (
                              object.ERROR_CODE
                          ))

        object.add_method("get_periastron",
                          (),
                          (
                              nbody_system.length,
                              object.ERROR_CODE
                          ))

        object.add_method("get_apastron",
                          (),
                          (
                              nbody_system.length,
                              object.ERROR_CODE
                          ))

        object.add_method("get_velocity_vector",
                          (),
                          (
                              nbody_system.speed,
                              nbody_system.speed,
                              nbody_system.speed,
                              object.ERROR_CODE
                          ))

        object.add_method("get_angles",
                          (),
                          (
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.ERROR_CODE
                          ))

        object.add_method("set_longitudinal_unit_vector",
                          (
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.NO_UNIT
                          ),
                          (
                              object.ERROR_CODE
                          ))

        object.add_method("set_transverse_unit_vector",
                          (
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.NO_UNIT
                          ),
                          (
                              object.ERROR_CODE
                          ))

        object.add_method("set_normal_unit_vector",
                          (
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.NO_UNIT
                          ),
                          (
                              object.ERROR_CODE
                          ))

        object.add_method("get_longitudinal_unit_vector",
                          (),
                          (
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.ERROR_CODE
                          ))

        object.add_method("get_transverse_unit_vector",
                          (),
                          (
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.ERROR_CODE
                          ))

        object.add_method("get_normal_unit_vector",
                          (),
                          (
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.NO_UNIT,
                              object.ERROR_CODE
                          ))
