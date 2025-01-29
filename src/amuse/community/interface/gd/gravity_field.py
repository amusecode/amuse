"""
Stellar Dynamics Interface Definition
"""

from amuse.units import nbody_system

from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification


class SinglePointGravityFieldInterface:
    """
    Codes implementing the gravity field interface provide functions to
    calculate the force and potential energy fields at any point.
    """

    @legacy_function
    def get_gravity_at_point():
        """
        Get the gravitational acceleration at the given points. To calculate
        the force on bodies at those points, multiply with the mass of the
        bodies
        """
        function = LegacyFunctionSpecification()
        for x in ["eps", "x", "y", "z"]:
            function.addParameter(
                x, dtype="float64", direction=function.IN, unit=nbody_system.length
            )
        for x in ["ax", "ay", "az"]:
            function.addParameter(
                x,
                dtype="float64",
                direction=function.OUT,
                unit=nbody_system.acceleration,
            )
        function.result_type = "int32"
        function.can_handle_array = True
        return function

    @legacy_function
    def get_potential_at_point():
        """
        Determine the gravitational potential on any given point
        """
        function = LegacyFunctionSpecification()
        for x in ["eps", "x", "y", "z"]:
            function.addParameter(
                x, dtype="float64", direction=function.IN, unit=nbody_system.length
            )
        for x in ["phi"]:
            function.addParameter(
                x, dtype="float64", direction=function.OUT, unit=nbody_system.potential
            )
        function.result_type = "int32"
        function.can_handle_array = True
        return function


class GravityFieldInterface:
    """
    Codes implementing the gravity field interface provide functions to
    calculate the force and potential energy fields at any point.
    """

    @legacy_function
    def get_gravity_at_point():
        """
        Get the gravitational acceleration at the given points. To calculate
        the force on bodies at those points, multiply with the mass of the
        bodies
        """
        function = LegacyFunctionSpecification()
        for x in ["eps", "x", "y", "z"]:
            function.addParameter(
                x, dtype="float64", direction=function.IN, unit=nbody_system.length
            )
        for x in ["ax", "ay", "az"]:
            function.addParameter(
                x,
                dtype="float64",
                direction=function.OUT,
                unit=nbody_system.acceleration,
            )
        function.addParameter("npoints", dtype="i", direction=function.LENGTH)
        function.result_type = "int32"
        function.must_handle_array = True
        return function

    @legacy_function
    def get_potential_at_point():
        """
        Determine the gravitational potential on any given point
        """
        function = LegacyFunctionSpecification()
        for x in ["eps", "x", "y", "z"]:
            function.addParameter(
                x, dtype="float64", direction=function.IN, unit=nbody_system.length
            )
        for x in ["phi"]:
            function.addParameter(
                x, dtype="float64", direction=function.OUT, unit=nbody_system.potential
            )
        function.addParameter("npoints", dtype="i", direction=function.LENGTH)
        function.result_type = "int32"
        function.must_handle_array = True
        return function


class GravityFieldCode:

    def define_state(self, handler):
        handler.add_method("RUN", "get_gravity_at_point")
        handler.add_method("RUN", "get_potential_at_point")
