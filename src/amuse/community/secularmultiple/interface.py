from amuse.community import *
from amuse.units import units,constants

### units used in the legacy code ###
### numerical values are defined in src/types.h ###
unit_l = units.AU
unit_m = units.MSun
unit_t = 1.0e6*units.yr
unit_h = unit_m*unit_l**2/unit_t ### specific angular momentum
unit_e = unit_m*unit_l**2/(unit_t**2) ### energy

print_name = 'SecularMultiple'

class SecularMultipleInterface(CodeInterface,    LiteratureReferencesMixIn):
    """
    SecularMultiple -- by Adrian Hamers, based on 2016MNRAS.459.2827H

    June 2016

    A code to compute the secular (orbit-averaged) gravitational dynamics of hierarchical multiple systems composed of nested binary orbits (simplex-type systems). with any configuration and any number of bodies. A particle can repesent a binary (`is_binary = True') or a body (`is_binary = False'). The structure of the system is determined by linking to other particles with the attributes child1 and child2. Tidal interactions and relativistic corrections are included in an ad hoc fashion (tides: treating the companion as a single body, even if it is not; relativistic terms: only including binary-binary interactions).

        .. [#] Hamers & Portegies Zwart, 2016, MNRAS 459, 2827

    """
    include_headers = ['interface.h','src/types.h','src/evolve.h','src/ODE_system.h']

    def __init__(self, **options):
#         CodeInterface.__init__(self, name_of_the_worker="secularmultiple_worker", **options)
         CodeInterface.__init__(self, **options)
         LiteratureReferencesMixIn.__init__(self)


    #######################
    ### basic interface ###
    #######################

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.OUT, unit=INDEX)
        function.addParameter('is_binary',              dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_children():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('child1',                 dtype='int32',      direction=function.IN,  unit=LINK('particles'))
        function.addParameter('child2',                 dtype='int32',      direction=function.IN,  unit=LINK('particles'))
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_children():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('child1',                 dtype='int32',      direction=function.OUT, unit=LINK('particles'))
        function.addParameter('child2',                 dtype='int32',      direction=function.OUT, unit=LINK('particles'))
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('mass',                   dtype='float64',    direction=function.IN,  unit=unit_m)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('mass',                   dtype='float64',    direction=function.OUT, unit=unit_m)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('radius',                 dtype='float64',    direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('radius',                 dtype='float64',    direction=function.OUT, unit=unit_l)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_level():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('level',                  dtype='int32',      direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    ####################
    ### spin vectors ###
    ####################

    @legacy_function
    def set_spin_vector():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('spin_vec_x',             dtype='float64',    direction=function.IN,  unit=1.0/unit_t)
        function.addParameter('spin_vec_y',             dtype='float64',    direction=function.IN,  unit=1.0/unit_t)
        function.addParameter('spin_vec_z',             dtype='float64',    direction=function.IN,  unit=1.0/unit_t)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_spin_vector():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('spin_vec_x',             dtype='float64',    direction=function.OUT, unit=1.0/unit_t)
        function.addParameter('spin_vec_y',             dtype='float64',    direction=function.OUT, unit=1.0/unit_t)
        function.addParameter('spin_vec_z',             dtype='float64',    direction=function.OUT, unit=1.0/unit_t)
        function.result_type = 'int32'
        return function


    ################################
    ### orbital vectors/elements ###
    ################################

    @legacy_function
    def set_orbital_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('e_vec_x',                dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('e_vec_y',                dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('e_vec_z',                dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('h_vec_x',                dtype='float64',    direction=function.IN,  unit=unit_h)
        function.addParameter('h_vec_y',                dtype='float64',    direction=function.IN,  unit=unit_h)
        function.addParameter('h_vec_z',                dtype='float64',    direction=function.IN,  unit=unit_h)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_orbital_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('e_vec_x',                dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('e_vec_y',                dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('e_vec_z',                dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('h_vec_x',                dtype='float64',    direction=function.OUT, unit=unit_h)
        function.addParameter('h_vec_y',                dtype='float64',    direction=function.OUT, unit=unit_h)
        function.addParameter('h_vec_z',                dtype='float64',    direction=function.OUT, unit=unit_h)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_orbital_elements():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',      dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('semimajor_axis',             dtype='float64',    direction=function.IN,  unit=unit_l)
        function.addParameter('eccentricity',               dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('inclination',                dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('argument_of_pericenter',     dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('longitude_of_ascending_node',dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_orbital_elements():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',      dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('semimajor_axis',             dtype='float64',    direction=function.OUT, unit=unit_l)
        function.addParameter('eccentricity',               dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('inclination',                dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('argument_of_pericenter',     dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('longitude_of_ascending_node',dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_inclination_relative_to_parent():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('inclination_relative_to_parent', dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_de_dt():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('de_dt',                  dtype='float64',    direction=function.OUT, unit=1.0/unit_t)
        function.result_type = 'int32'
        return function


    ################
    ### PN terms ###
    ################

    @legacy_function
    def set_include_pairwise_1PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',      dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_1PN_terms', dtype='bool',       direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_pairwise_1PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',      dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_1PN_terms', dtype='bool',       direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_include_pairwise_25PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',      dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_25PN_terms',dtype='bool',       direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_pairwise_25PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',      dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_25PN_terms',dtype='bool',       direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    #############
    ### tides ###
    #############

    @legacy_function
    def set_tides_method():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('tides_method',                   dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tides_method():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('tides_method',                   dtype='int32',      direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_tidal_friction_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_friction_terms',   dtype='bool',       direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_tidal_friction_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_friction_terms',   dtype='bool',       direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_tidal_bulges_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                dtype='int32', direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_bulges_precession_terms',dtype='bool',  direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_tidal_bulges_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                 dtype='int32', direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_bulges_precession_terms', dtype='bool',  direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_rotation_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',              dtype='int32', direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_rotation_precession_terms',  dtype='bool',  direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_rotation_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',              dtype='int32', direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_rotation_precession_terms',  dtype='bool',  direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_minimum_eccentricity_for_tidal_precession():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                     dtype='int32',  direction=function.IN, unit=NO_UNIT)
        function.addParameter('minimum_eccentricity_for_tidal_precession', dtype='float64',direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_minimum_eccentricity_for_tidal_precession():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                     dtype='int32', direction=function.IN, unit=NO_UNIT)
        function.addParameter('minimum_eccentricity_for_tidal_precession', dtype='float64',direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tides_apsidal_motion_constant():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('tides_apsidal_motion_constant',          dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_tides_apsidal_motion_constant():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,   unit=NO_UNIT)
        function.addParameter('tides_apsidal_motion_constant',          dtype='float64',    direction=function.OUT,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tides_gyration_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('tides_gyration_radius',                  dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_tides_gyration_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,   unit=NO_UNIT)
        function.addParameter('tides_gyration_radius',                  dtype='float64',    direction=function.OUT,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tides_viscous_time_scale():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('tides_viscous_time_scale',               dtype='float64',    direction=function.IN,  unit=unit_t)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_tides_viscous_time_scale():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,   unit=NO_UNIT)
        function.addParameter('tides_viscous_time_scale',               dtype='float64',    direction=function.OUT,  unit=unit_t)
        function.result_type = 'int32'
        return function

    ####################
    ### root finding ###
    ####################

    ### secular breakdown ###
    @legacy_function
    def set_check_for_secular_breakdown():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_secular_breakdown',    dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_secular_breakdown():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_secular_breakdown',    dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    ### dynamical instablity ###
    @legacy_function
    def set_check_for_dynamical_instability():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_dynamical_instability',dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_dynamical_instability():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_dynamical_instability',dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dynamical_instability_criterion():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_criterion',dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dynamical_instability_criterion():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_criterion',dtype='int32',      direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dynamical_instability_central_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                   dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_central_particle',  dtype='int32',      direction=function.IN,  unit=LINK('particles'))
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dynamical_instability_central_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                   dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_central_particle',  dtype='int32',      direction=function.OUT, unit=LINK('particles'))
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dynamical_instability_K_parameter():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                   dtype='float64',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_K_parameter',       dtype='float64',      direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dynamical_instability_K_parameter():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                   dtype='float64',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_K_parameter',       dtype='float64',      direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    ### physical collision / orbit crossing ###
    @legacy_function
    def set_check_for_physical_collision_or_orbit_crossing():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_physical_collision_or_orbit_crossing', dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_physical_collision_or_orbit_crossing():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_physical_collision_or_orbit_crossing', dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    ### minimum periapse distance reached ###
    @legacy_function
    def set_check_for_minimum_periapse_distance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_minimum_periapse_distance',            dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_minimum_periapse_distance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_minimum_periapse_distance',            dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_minimum_periapse_distance_value():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_minimum_periapse_distance_value',      dtype='float64',    direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_minimum_periapse_distance_value():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_minimum_periapse_distance_value',      dtype='float64',    direction=function.OUT, unit=unit_l)
        function.result_type = 'int32'
        return function


    ### RLOF at pericentre ###
    @legacy_function
    def set_check_for_RLOF_at_pericentre():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_RLOF_at_pericentre',   dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_RLOF_at_pericentre():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_RLOF_at_pericentre',   dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_RLOF_at_pericentre_use_sepinsky_fit():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                              dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_RLOF_at_pericentre_use_sepinsky_fit',      dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_RLOF_at_pericentre_use_sepinsky_fit():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('check_for_RLOF_at_pericentre_use_sepinsky_fit',  dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    ### retrieve root finding state ###
    @legacy_function
    def set_root_finding_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                              dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('secular_breakdown_has_occurred',                     dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_has_occurred',                 dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.addParameter('physical_collision_or_orbit_crossing_has_occurred',  dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.addParameter('minimum_periapse_distance_has_occurred',             dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.addParameter('RLOF_at_pericentre_has_occurred',                    dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_root_finding_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                              dtype='int32',  direction=function.IN,  unit=NO_UNIT)
        function.addParameter('secular_breakdown_has_occurred',                     dtype='bool',   direction=function.OUT, unit=NO_UNIT)
        function.addParameter('dynamical_instability_has_occurred',                 dtype='bool',   direction=function.OUT, unit=NO_UNIT)
        function.addParameter('physical_collision_or_orbit_crossing_has_occurred',  dtype='bool',   direction=function.OUT, unit=NO_UNIT)
        function.addParameter('minimum_periapse_distance_has_occurred',             dtype='bool',   direction=function.OUT, unit=NO_UNIT)
        function.addParameter('RLOF_at_pericentre_has_occurred',                    dtype='bool',   direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    ########################
    ### evolve interface ###
    ########################

    @legacy_function
    def evolve_interface():
        function = LegacyFunctionSpecification()
        function.addParameter('start_time',             dtype='float64',    direction=function.IN,  unit=unit_t)
        function.addParameter('time_step',              dtype='float64',    direction=function.IN,  unit=unit_t)
        function.addParameter('output_time',            dtype='float64',    direction=function.OUT, unit=unit_t)
        function.addParameter('hamiltonian',            dtype='float64',    direction=function.OUT, unit=unit_e)
        function.addParameter('flag',                   dtype='int32',      direction=function.OUT, unit=NO_UNIT)
        function.addParameter('error_code',             dtype='int32',      direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def determine_binary_parents_levels_and_masses_interface():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function


    #######################
    ### code parameters ###
    #######################


    ##################
    ### tolerances ###
    ##################

    @legacy_function
    def get_relative_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_relative_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_absolute_tolerance_eccentricity_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_absolute_tolerance_eccentricity_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_absolute_tolerance_angular_momentum_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_absolute_tolerance_angular_momentum_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_absolute_tolerance_spin_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_absolute_tolerance_spin_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_compute_orbital_elements_with_respect_to_total_angular_momentum_vector():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_compute_orbital_elements_with_respect_to_total_angular_momentum_vector():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    #############
    ### terms ###
    #############

    @legacy_function
    def get_include_quadrupole_order_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_quadrupole_order_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_octupole_order_binary_pair_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_octupole_order_binary_pair_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_octupole_order_binary_triplet_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_octupole_order_binary_triplet_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_hexadecupole_order_binary_pair_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_hexadecupole_order_binary_pair_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_dotriacontupole_order_binary_pair_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_dotriacontupole_order_binary_pair_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='bool',       direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

class SecularMultiple(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  SecularMultipleInterface(**options), **options)
        self.model_time = 0.0 | units.Myr
        self.particles_committed = False

        self.initial_hamiltonian = 0.0 | unit_e
        self.hamiltonian = 0.0 | unit_e
        self.flag = 0
        self.error_code = 0

#        self.verbose = True
#        self.debug = True

    def define_particle_sets(self,object):

        object.define_set('particles', 'index_of_the_particle')
        object.set_new   ('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')

        object.add_setter('particles', 'set_children', names = ('child1','child2') )
        object.add_getter('particles', 'get_children', names = ('child1','child2') )
        object.add_setter('particles', 'set_mass')
        object.add_getter('particles', 'get_mass')
        object.add_setter('particles', 'set_radius')
        object.add_getter('particles', 'get_radius')
        object.add_getter('particles', 'get_level')


        object.add_setter('particles', 'set_spin_vector')
        object.add_getter('particles', 'get_spin_vector')
        object.add_setter('particles', 'set_orbital_vectors')
        object.add_getter('particles', 'get_orbital_vectors')
        object.add_setter('particles', 'set_orbital_elements', names = ('semimajor_axis','eccentricity','inclination','argument_of_pericenter','longitude_of_ascending_node') )
        object.add_getter('particles', 'get_orbital_elements', names = ('semimajor_axis','eccentricity','inclination','argument_of_pericenter','longitude_of_ascending_node') )
        object.add_getter('particles', 'get_inclination_relative_to_parent')
        object.add_getter('particles', 'get_de_dt')

        object.add_setter('particles', 'set_include_pairwise_1PN_terms')
        object.add_getter('particles', 'get_include_pairwise_1PN_terms')
        object.add_setter('particles', 'set_include_pairwise_25PN_terms')
        object.add_getter('particles', 'get_include_pairwise_25PN_terms')

        object.add_setter('particles', 'set_tides_method')
        object.add_getter('particles', 'get_tides_method')
        object.add_setter('particles', 'set_include_tidal_friction_terms')
        object.add_getter('particles', 'get_include_tidal_friction_terms')
        object.add_setter('particles', 'set_include_tidal_bulges_precession_terms')
        object.add_getter('particles', 'get_include_tidal_bulges_precession_terms')
        object.add_setter('particles', 'set_include_rotation_precession_terms')
        object.add_getter('particles', 'get_include_rotation_precession_terms')
        object.add_setter('particles', 'set_minimum_eccentricity_for_tidal_precession')
        object.add_getter('particles', 'get_minimum_eccentricity_for_tidal_precession')

        object.add_setter('particles', 'set_tides_apsidal_motion_constant')
        object.add_getter('particles', 'get_tides_apsidal_motion_constant')
        object.add_setter('particles', 'set_tides_gyration_radius')
        object.add_getter('particles', 'get_tides_gyration_radius')
        object.add_setter('particles', 'set_tides_viscous_time_scale')
        object.add_getter('particles', 'get_tides_viscous_time_scale')

        object.add_setter('particles', 'set_check_for_dynamical_instability')
        object.add_getter('particles', 'get_check_for_dynamical_instability')
        object.add_setter('particles', 'set_dynamical_instability_criterion')
        object.add_getter('particles', 'get_dynamical_instability_criterion')
        object.add_setter('particles', 'set_dynamical_instability_central_particle')
        object.add_getter('particles', 'get_dynamical_instability_central_particle')
        object.add_setter('particles', 'set_dynamical_instability_K_parameter')
        object.add_getter('particles', 'get_dynamical_instability_K_parameter')

        object.add_setter('particles', 'set_check_for_physical_collision_or_orbit_crossing')
        object.add_getter('particles', 'get_check_for_physical_collision_or_orbit_crossing')

        object.add_setter('particles', 'set_check_for_minimum_periapse_distance')
        object.add_getter('particles', 'get_check_for_minimum_periapse_distance')
        object.add_setter('particles', 'set_check_for_minimum_periapse_distance_value')
        object.add_getter('particles', 'get_check_for_minimum_periapse_distance_value')

        object.add_setter('particles', 'set_check_for_RLOF_at_pericentre')
        object.add_getter('particles', 'get_check_for_RLOF_at_pericentre')
        object.add_setter('particles', 'set_check_for_RLOF_at_pericentre_use_sepinsky_fit')
        object.add_getter('particles', 'get_check_for_RLOF_at_pericentre_use_sepinsky_fit')

        object.add_setter('particles', 'set_root_finding_state')
        object.add_getter('particles', 'get_root_finding_state')

    def define_parameters(self, object):

        object.add_method_parameter(
            "get_relative_tolerance",
            "set_relative_tolerance",
            "relative_tolerance",
            "relative_tolerance",
            default_value = 1.0e-16
        )
        object.add_method_parameter(
            "get_absolute_tolerance_eccentricity_vectors",
            "set_absolute_tolerance_eccentricity_vectors",
            "absolute_tolerance_eccentricity_vectors",
            "absolute_tolerance_eccentricity_vectors",
            default_value = 1.0e-14
        )
        object.add_method_parameter(
            "get_absolute_tolerance_angular_momentum_vectors",
            "set_absolute_tolerance_angular_momentum_vectors",
            "absolute_tolerance_angular_momentum_vectors",
            "absolute_tolerance_angular_momentum_vectors",
            default_value = 1.0e-2 | unit_h
        )
        object.add_method_parameter(
            "get_absolute_tolerance_spin_vectors",
            "set_absolute_tolerance_spin_vectors",
            "absolute_tolerance_spin_vectors",
            "absolute_tolerance_spin_vectors",
            default_value = 1.0e4 | (1.0/unit_t)
        )
        object.add_method_parameter(
            "get_compute_orbital_elements_with_respect_to_total_angular_momentum_vector",
            "set_compute_orbital_elements_with_respect_to_total_angular_momentum_vector",
            "compute_orbital_elements_with_respect_to_total_angular_momentum_vector",
            "compute_orbital_elements_with_respect_to_total_angular_momentum_vector",
            default_value = False
        )
        object.add_method_parameter(
            "get_include_quadrupole_order_terms",
            "set_include_quadrupole_order_terms",
            "include_quadrupole_order_terms",
            "include_quadrupole_order_terms",
            default_value = True
        )
        object.add_method_parameter(
            "get_include_octupole_order_binary_pair_terms",
            "set_include_octupole_order_binary_pair_terms",
            "include_octupole_order_binary_pair_terms",
            "include_octupole_order_binary_pair_terms",
            default_value = True
        )
        object.add_method_parameter(
            "get_include_octupole_order_binary_triplet_terms",
            "set_include_octupole_order_binary_triplet_terms",
            "include_octupole_order_binary_triplet_terms",
            "include_octupole_order_binary_triplet_terms",
            default_value = False
        )
        object.add_method_parameter(
            "get_include_hexadecupole_order_binary_pair_terms",
            "set_include_hexadecupole_order_binary_pair_terms",
            "include_hexadecupole_order_binary_pair_terms",
            "include_hexadecupole_order_binary_pair_terms",
            default_value = False
        )
        object.add_method_parameter(
            "get_include_dotriacontupole_order_binary_pair_terms",
            "set_include_dotriacontupole_order_binary_pair_terms",
            "include_dotriacontupole_order_binary_pair_terms",
            "include_dotriacontupole_order_binary_pair_terms",
            default_value = False
        )


    def commit_particles(self):
        print print_name,' -- committing particles'
        particles = self.particles

        if len(particles) == 0:
            print print_name,' -- no particles have been added -- exiting'
            exit(-1)

        particles.add_vector_attribute("spin_vec",["spin_vec_x","spin_vec_y","spin_vec_z"])
        particles.add_vector_attribute("e_vec",["e_vec_x","e_vec_y","e_vec_z"])
        particles.add_vector_attribute("h_vec",["h_vec_x","h_vec_y","h_vec_z"])

        ### evaluate the initial hamiltonian ###
        time_step = 0.0 | units.Myr
        end_time,self.initial_hamiltonian,flag,error_code = self.evolve_interface(self.model_time,time_step)

        self.particles_committed = True

    def evolve_model(self,end_time):
        if end_time is None:
            print print_name,' -- end time not specified in evolve_model! exiting'
            exit(-1)
        if self.particles_committed == False:
            self.commit_particles()

        ### integrate system of ODEs ###
        start_time = self.model_time
        time_step = end_time - start_time
        end_time,self.hamiltonian,flag,error_code = self.evolve_interface(start_time,time_step)

        ### compute energy error ###
        self.relative_energy_error = numpy.fabs( (self.initial_hamiltonian - self.hamiltonian)/self.initial_hamiltonian )

        ### update model time ###
        self.model_time = end_time

        if (flag==99):
            print print_name,' -- error occurred during ODE integration'
            print print_name,' -- error code is ',error_code
        self.flag = flag
        self.error_code = error_code

    def determine_binary_parents_levels_and_masses(self):
        self.determine_binary_parents_levels_and_masses_interface()
