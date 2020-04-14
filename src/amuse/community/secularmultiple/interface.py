from amuse.community import *
from amuse.units import units,constants

### units used in the legacy code ###
### numerical values are defined in src/types.h ###
unit_l = units.AU
unit_m = units.MSun
unit_t = 1.0e6*units.yr
unit_h = unit_m*unit_l**2/unit_t ### specific angular momentum
unit_e = unit_m*unit_l**2/(unit_t**2) ### energy
unit_lum = unit_e/unit_t

print_name = 'SecularMultipleVar'

class SecularMultipleInterface(CodeInterface):
    """
    SecularMultiple -- by Adrian Hamers, based on 2016MNRAS.459.2827H
    
    A code to compute the secular (orbit-averaged) gravitational dynamics of hierarchical multiple systems composed of nested binary orbits (simplex-type systems). with any configuration and any number of bodies. A particle can repesent a binary (`is_binary = True') or a body (`is_binary = False'). The structure of the system is determined by linking to other particles with the attributes child1 and child2. Tidal interactions and relativistic corrections are included in an ad hoc fashion (tides: treating the companion as a single body, even if it is not; relativistic terms: only including binary-binary interactions).
    
    November 2017: Updates for external perturbations (flybys & supernovae), detailed in Hamers (2018, in prep)
    
    """
    include_headers = ['interface.h','src/types.h','src/evolve.h','src/ODE_system.h']

    def __init__(self, **options):
#         CodeInterface.__init__(self, name_of_the_worker="secularmultiple_worker", **options)
         CodeInterface.__init__(self, **options)

    #######################
    ### basic interface ###
    #######################

    ### particles ###
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
    def set_mass_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('mass_dot_external',      dtype='float64',    direction=function.IN,  unit=unit_m/unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_mass_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('mass_dot_external',      dtype='float64',    direction=function.OUT, unit=unit_m/unit_t)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('radius',                   dtype='float64',    direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('radius',                   dtype='float64',    direction=function.OUT, unit=unit_l)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_radius_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('radius_dot_external',    dtype='float64',    direction=function.IN,  unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_radius_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('radius_dot_external',    dtype='float64',    direction=function.OUT, unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_radius_ddot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('radius_ddot_external',    dtype='float64',    direction=function.IN,  unit=unit_l/(unit_t**2))
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_radius_ddot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('radius_ddot_external',    dtype='float64',    direction=function.OUT, unit=unit_l/(unit_t**2))
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

    @legacy_function
    def set_stellar_type():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('stellar_type',           dtype='int32',      direction=function.IN,  unit=units.stellar_type)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_stellar_type():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('stellar_type',           dtype='int32',      direction=function.OUT, unit=units.stellar_type)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_true_anomaly():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('true_anomaly',           dtype='float64',      direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_true_anomaly():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('true_anomaly',           dtype='float64',      direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_sample_orbital_phases_randomly():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',              dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('sample_orbital_phases_randomly',     dtype='bool',      direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_sample_orbital_phases_randomly():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',              dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('sample_orbital_phases_randomly',     dtype='bool',      direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function




    #################################################
    ### user-specified instantaneous perturbation ###
    #################################################
    
    @legacy_function
    def set_instantaneous_perturbation_delta_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                      dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('instantaneous_perturbation_delta_mass',      dtype='float64',      direction=function.IN,  unit=unit_m)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_instantaneous_perturbation_delta_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                      dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('instantaneous_perturbation_delta_mass',      dtype='float64',      direction=function.OUT, unit=unit_m)
        function.result_type = 'int32'
        return function
    
    
    @legacy_function
    def set_instantaneous_perturbation_delta_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                      dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('instantaneous_perturbation_delta_position_x',dtype='float64',      direction=function.IN,  unit=unit_l)
        function.addParameter('instantaneous_perturbation_delta_position_y',dtype='float64',      direction=function.IN,  unit=unit_l)
        function.addParameter('instantaneous_perturbation_delta_position_z',dtype='float64',      direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_instantaneous_perturbation_delta_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                      dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('instantaneous_perturbation_delta_position_x',dtype='float64',      direction=function.OUT, unit=unit_l)
        function.addParameter('instantaneous_perturbation_delta_position_y',dtype='float64',      direction=function.OUT, unit=unit_l)
        function.addParameter('instantaneous_perturbation_delta_position_z',dtype='float64',      direction=function.OUT, unit=unit_l)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_instantaneous_perturbation_delta_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                      dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('instantaneous_perturbation_delta_velocity_x',dtype='float64',      direction=function.IN,  unit=unit_l/unit_t)
        function.addParameter('instantaneous_perturbation_delta_velocity_y',dtype='float64',      direction=function.IN,  unit=unit_l/unit_t)
        function.addParameter('instantaneous_perturbation_delta_velocity_z',dtype='float64',      direction=function.IN,  unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_instantaneous_perturbation_delta_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                      dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('instantaneous_perturbation_delta_velocity_x',dtype='float64',      direction=function.OUT, unit=unit_l/unit_t)
        function.addParameter('instantaneous_perturbation_delta_velocity_y',dtype='float64',      direction=function.OUT, unit=unit_l/unit_t)
        function.addParameter('instantaneous_perturbation_delta_velocity_z',dtype='float64',      direction=function.OUT, unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function




    ##########################
    ### external particles ###
    ##########################
    
    @legacy_function
    def new_external_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.OUT, unit=INDEX)
        function.addParameter('mass',                   dtype='float64',    direction=function.IN,  unit=unit_m)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def delete_external_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.result_type = 'int32'
        return function
        
        
    @legacy_function
    def set_external_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('mass',                   dtype='float64',    direction=function.IN,  unit=unit_m)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('mass',                   dtype='float64',    direction=function.OUT, unit=unit_m)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_external_path():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('path',                   dtype='int32',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_path():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('path',                   dtype='int32',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_external_mode():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('mode',                   dtype='int32',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_mode():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('mode',                   dtype='int32',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
        
    @legacy_function
    def set_external_t_ref():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('t_ref',                  dtype='float64',    direction=function.IN,  unit=unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_t_ref():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('t_ref',                  dtype='float64',    direction=function.OUT, unit=unit_t)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_external_t_passed():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('t_passed',               dtype='float64',    direction=function.IN,  unit=unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_t_passed():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('t_passed',               dtype='float64',    direction=function.OUT, unit=unit_t)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_external_r0_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('r0_vec_x',               dtype='float64',    direction=function.IN,  unit=unit_l)
        function.addParameter('r0_vec_y',               dtype='float64',    direction=function.IN,  unit=unit_l)
        function.addParameter('r0_vec_z',               dtype='float64',    direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_r0_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('r0_vec_x',               dtype='float64',    direction=function.OUT, unit=unit_l)
        function.addParameter('r0_vec_y',               dtype='float64',    direction=function.OUT, unit=unit_l)
        function.addParameter('r0_vec_z',               dtype='float64',    direction=function.OUT, unit=unit_l)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_external_rdot_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('rdot_vec_x',             dtype='float64',    direction=function.IN,  unit=unit_l/unit_t)
        function.addParameter('rdot_vec_y',             dtype='float64',    direction=function.IN,  unit=unit_l/unit_t)
        function.addParameter('rdot_vec_z',             dtype='float64',    direction=function.IN,  unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_rdot_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('rdot_vec_x',             dtype='float64',    direction=function.OUT, unit=unit_l/unit_t)
        function.addParameter('rdot_vec_y',             dtype='float64',    direction=function.OUT, unit=unit_l/unit_t)
        function.addParameter('rdot_vec_z',             dtype='float64',    direction=function.OUT, unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_external_periapse_distance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('periapse_distance',               dtype='float64',    direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_periapse_distance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('periapse_distance',               dtype='float64',    direction=function.OUT, unit=unit_l)
        function.result_type = 'int32'
        return function
        
        
    @legacy_function
    def set_external_eccentricity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('eccentricity',               dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_eccentricity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('eccentricity',               dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function
                
        
    @legacy_function
    def set_external_e_hat_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('e_hat_vec_x',             dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('e_hat_vec_y',             dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('e_hat_vec_z',             dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_e_hat_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('e_hat_vec_x',             dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('e_hat_vec_y',             dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('e_hat_vec_z',             dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_external_h_hat_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('h_hat_vec_x',             dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('h_hat_vec_y',             dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('h_hat_vec_z',             dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_external_h_hat_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('h_hat_vec_x',             dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('h_hat_vec_y',             dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('h_hat_vec_z',             dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_external_r_vectors():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_external_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)        
        function.addParameter('r_vec_x',             dtype='float64',    direction=function.OUT, unit=unit_l)
        function.addParameter('r_vec_y',             dtype='float64',    direction=function.OUT, unit=unit_l)
        function.addParameter('r_vec_z',             dtype='float64',    direction=function.OUT, unit=unit_l)
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

    @legacy_function
    def set_spin_vector_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('spin_vec_x_dot_external',dtype='float64',    direction=function.IN,  unit=1.0/(unit_t**2))
        function.addParameter('spin_vec_y_dot_external',dtype='float64',    direction=function.IN,  unit=1.0/(unit_t**2))
        function.addParameter('spin_vec_z_dot_external',dtype='float64',    direction=function.IN,  unit=1.0/(unit_t**2))
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_spin_vector_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('spin_vec_x_dot_external',dtype='float64',    direction=function.OUT, unit=1.0/(unit_t**2))
        function.addParameter('spin_vec_y_dot_external',dtype='float64',    direction=function.OUT, unit=1.0/(unit_t**2))
        function.addParameter('spin_vec_z_dot_external',dtype='float64',    direction=function.OUT, unit=1.0/(unit_t**2))
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
    def set_orbital_vectors_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('e_vec_x_dot_external',   dtype='float64',    direction=function.IN,  unit=1.0/unit_t)
        function.addParameter('e_vec_y_dot_external',   dtype='float64',    direction=function.IN,  unit=1.0/unit_t)
        function.addParameter('e_vec_z_dot_external',   dtype='float64',    direction=function.IN,  unit=1.0/unit_t)
        function.addParameter('h_vec_x_dot_external',   dtype='float64',    direction=function.IN,  unit=unit_h/unit_t)
        function.addParameter('h_vec_y_dot_external',   dtype='float64',    direction=function.IN,  unit=unit_h/unit_t)
        function.addParameter('h_vec_z_dot_external',   dtype='float64',    direction=function.IN,  unit=unit_h/unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_orbital_vectors_dot_external():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('e_vec_x_dot_external',   dtype='float64',    direction=function.OUT, unit=1.0/unit_t)
        function.addParameter('e_vec_y_dot_external',   dtype='float64',    direction=function.OUT, unit=1.0/unit_t)
        function.addParameter('e_vec_z_dot_external',   dtype='float64',    direction=function.OUT, unit=1.0/unit_t)
        function.addParameter('h_vec_x_dot_external',   dtype='float64',    direction=function.OUT, unit=unit_h/unit_t)
        function.addParameter('h_vec_y_dot_external',   dtype='float64',    direction=function.OUT, unit=unit_h/unit_t)
        function.addParameter('h_vec_z_dot_external',   dtype='float64',    direction=function.OUT, unit=unit_h/unit_t)
        function.result_type = 'int32'
        return function            


    @legacy_function
    def set_orbital_elements():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('semimajor_axis',         dtype='float64',    direction=function.IN,  unit=unit_l)
        function.addParameter('eccentricity',           dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('inclination',            dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('argument_of_pericenter', dtype='float64',    direction=function.IN,  unit=NO_UNIT)
        function.addParameter('longitude_of_ascending_node', dtype='float64',direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_orbital_elements():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('semimajor_axis',         dtype='float64',    direction=function.OUT, unit=unit_l)
        function.addParameter('eccentricity',           dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('inclination',            dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('argument_of_pericenter', dtype='float64',    direction=function.OUT, unit=NO_UNIT)
        function.addParameter('longitude_of_ascending_node', dtype='float64',direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function 

        
    @legacy_function
    def get_inclination_relative_to_parent():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('inclination_relative_to_parent',dtype='float64',direction=function.OUT,unit=NO_UNIT)
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



    @legacy_function
    def set_position_vector():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('position_x',      dtype='float64',    direction=function.IN,  unit=unit_l)
        function.addParameter('position_y',      dtype='float64',    direction=function.IN,  unit=unit_l)
        function.addParameter('position_z',      dtype='float64',    direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_position_vector():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('position_x',      dtype='float64',    direction=function.OUT,  unit=unit_l)
        function.addParameter('position_y',      dtype='float64',    direction=function.OUT,  unit=unit_l)
        function.addParameter('position_z',      dtype='float64',    direction=function.OUT,  unit=unit_l)
        function.result_type = 'int32'
        return function    


    @legacy_function
    def set_velocity_vector():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('velocity_x',      dtype='float64',    direction=function.IN,  unit=unit_l/unit_t)
        function.addParameter('velocity_y',      dtype='float64',    direction=function.IN,  unit=unit_l/unit_t)
        function.addParameter('velocity_z',      dtype='float64',    direction=function.IN,  unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_velocity_vector():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN,  unit=INDEX)
        function.addParameter('velocity_x',      dtype='float64',    direction=function.OUT,  unit=unit_l/unit_t)
        function.addParameter('velocity_y',      dtype='float64',    direction=function.OUT,  unit=unit_l/unit_t)
        function.addParameter('velocity_z',      dtype='float64',    direction=function.OUT,  unit=unit_l/unit_t)
        function.result_type = 'int32'
        return function    




    ################
    ### PN terms ###
    ################
            
    @legacy_function
    def set_include_pairwise_1PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_1PN_terms', dtype='bool',   direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_pairwise_1PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_1PN_terms', dtype='bool',   direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_include_pairwise_25PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_25PN_terms', dtype='bool',  direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_pairwise_25PN_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_pairwise_25PN_terms', dtype='bool',  direction=function.OUT,unit=NO_UNIT)
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
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_friction_terms', dtype='bool',  direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_tidal_friction_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_friction_terms', dtype='bool',  direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_tidal_bulges_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_bulges_precession_terms', dtype='bool',  direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_tidal_bulges_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_tidal_bulges_precession_terms', dtype='bool',  direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_rotation_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_rotation_precession_terms', dtype='bool',  direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_rotation_precession_terms():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('include_rotation_precession_terms', dtype='bool',  direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_minimum_eccentricity_for_tidal_precession():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('minimum_eccentricity_for_tidal_precession', dtype='float64',  direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_minimum_eccentricity_for_tidal_precession():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('minimum_eccentricity_for_tidal_precession', dtype='float64',  direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    ### physical parameters ###
            
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

    @legacy_function
    def set_tides_viscous_time_scale_prescription():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('tides_viscous_time_scale_prescription',  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_tides_viscous_time_scale_prescription():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN, unit=NO_UNIT)
        function.addParameter('tides_viscous_time_scale_prescription',  dtype='int32',      direction=function.OUT,unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_convective_envelope_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('convective_envelope_mass',               dtype='float64',    direction=function.IN,  unit=unit_m)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_convective_envelope_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,   unit=NO_UNIT)
        function.addParameter('convective_envelope_mass',               dtype='float64',    direction=function.OUT,  unit=unit_m)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_convective_envelope_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('convective_envelope_radius',             dtype='float64',    direction=function.IN,  unit=unit_l)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_convective_envelope_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,   unit=NO_UNIT)
        function.addParameter('convective_envelope_radius',             dtype='float64',    direction=function.OUT,  unit=unit_l)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_luminosity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('luminosity',                             dtype='float64',    direction=function.IN,  unit=unit_lum)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_luminosity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                  dtype='int32',      direction=function.IN,   unit=NO_UNIT)
        function.addParameter('luminosity',                             dtype='float64',    direction=function.OUT,  unit=unit_lum)
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
        function.addParameter('index_of_the_particle',                              dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_central_particle',             dtype='int32',      direction=function.IN,  unit=LINK('particles'))
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dynamical_instability_central_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                              dtype='int32',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_central_particle',             dtype='int32',      direction=function.OUT, unit=LINK('particles'))
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_dynamical_instability_K_parameter():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='float64',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_K_parameter',              dtype='float64',      direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dynamical_instability_K_parameter():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle',                          dtype='float64',      direction=function.IN,  unit=NO_UNIT)
        function.addParameter('dynamical_instability_K_parameter',              dtype='float64',      direction=function.OUT, unit=NO_UNIT)
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

    @legacy_function
    def apply_external_perturbation_assuming_integrated_orbits_interface():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function


    @legacy_function
    def apply_user_specified_instantaneous_perturbation_interface():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_positions_and_velocities_interface():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
        
    #######################
    ### code parameters ###
    #######################
    
    @legacy_function
    def get_orbital_phases_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='int32',    direction=function.OUT, unit=NO_UNIT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_orbital_phases_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('value',                  dtype='int32',    direction=function.IN,  unit=NO_UNIT)
        function.result_type = 'int32'
        return function

        
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
       
    def define_particle_sets(self,handler):
       
        handler.define_set('particles', 'index_of_the_particle')
        handler.set_new   ('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        
        handler.add_setter('particles', 'set_children', names = ('child1','child2') )
        handler.add_getter('particles', 'get_children', names = ('child1','child2') )
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass')
        handler.add_setter('particles', 'set_mass_dot_external')
        handler.add_getter('particles', 'get_mass_dot_external')
        handler.add_setter('particles', 'set_radius')
        handler.add_getter('particles', 'get_radius')
        handler.add_setter('particles', 'set_radius_dot_external')
        handler.add_getter('particles', 'get_radius_dot_external')
        handler.add_setter('particles', 'set_radius_ddot_external')
        handler.add_getter('particles', 'get_radius_ddot_external')

        handler.add_getter('particles', 'get_level')
        
        handler.add_setter('particles', 'set_stellar_type')
        handler.add_getter('particles', 'get_stellar_type')

        handler.add_setter('particles', 'set_true_anomaly')
        handler.add_getter('particles', 'get_true_anomaly')

        handler.add_setter('particles', 'set_sample_orbital_phases_randomly')
        handler.add_getter('particles', 'get_sample_orbital_phases_randomly')

        handler.add_setter('particles', 'set_instantaneous_perturbation_delta_mass')
        handler.add_getter('particles', 'get_instantaneous_perturbation_delta_mass')
        handler.add_setter('particles', 'set_instantaneous_perturbation_delta_position')
        handler.add_getter('particles', 'get_instantaneous_perturbation_delta_position')
        handler.add_setter('particles', 'set_instantaneous_perturbation_delta_velocity')
        handler.add_getter('particles', 'get_instantaneous_perturbation_delta_velocity')

        handler.add_setter('particles', 'set_spin_vector')
        handler.add_getter('particles', 'get_spin_vector')
        handler.add_setter('particles', 'set_spin_vector_dot_external')
        handler.add_getter('particles', 'get_spin_vector_dot_external')
        handler.add_setter('particles', 'set_orbital_vectors')
        handler.add_getter('particles', 'get_orbital_vectors')
        handler.add_setter('particles', 'set_orbital_vectors_dot_external')
        handler.add_getter('particles', 'get_orbital_vectors_dot_external')
        handler.add_setter('particles', 'set_orbital_elements', names = ('semimajor_axis','eccentricity','inclination','argument_of_pericenter','longitude_of_ascending_node') )
        handler.add_getter('particles', 'get_orbital_elements', names = ('semimajor_axis','eccentricity','inclination','argument_of_pericenter','longitude_of_ascending_node') )
        handler.add_getter('particles', 'get_inclination_relative_to_parent')
        handler.add_getter('particles', 'get_de_dt')

        handler.add_setter('particles', 'set_position_vector')
        handler.add_getter('particles', 'get_position_vector')
        handler.add_setter('particles', 'set_velocity_vector')
        handler.add_getter('particles', 'get_velocity_vector')

        handler.add_setter('particles', 'set_include_pairwise_1PN_terms')
        handler.add_getter('particles', 'get_include_pairwise_1PN_terms')
        handler.add_setter('particles', 'set_include_pairwise_25PN_terms')
        handler.add_getter('particles', 'get_include_pairwise_25PN_terms')

        handler.add_setter('particles', 'set_tides_method')
        handler.add_getter('particles', 'get_tides_method')
        handler.add_setter('particles', 'set_include_tidal_friction_terms')
        handler.add_getter('particles', 'get_include_tidal_friction_terms')
        handler.add_setter('particles', 'set_include_tidal_bulges_precession_terms')
        handler.add_getter('particles', 'get_include_tidal_bulges_precession_terms')
        handler.add_setter('particles', 'set_include_rotation_precession_terms')
        handler.add_getter('particles', 'get_include_rotation_precession_terms')
        handler.add_setter('particles', 'set_minimum_eccentricity_for_tidal_precession')
        handler.add_getter('particles', 'get_minimum_eccentricity_for_tidal_precession')

        handler.add_setter('particles', 'set_tides_apsidal_motion_constant')
        handler.add_getter('particles', 'get_tides_apsidal_motion_constant')
        handler.add_setter('particles', 'set_tides_gyration_radius')
        handler.add_getter('particles', 'get_tides_gyration_radius')
        handler.add_setter('particles', 'set_tides_viscous_time_scale')
        handler.add_getter('particles', 'get_tides_viscous_time_scale')
        handler.add_setter('particles', 'set_tides_viscous_time_scale_prescription')
        handler.add_getter('particles', 'get_tides_viscous_time_scale_prescription')
        handler.add_setter('particles', 'set_convective_envelope_mass')
        handler.add_getter('particles', 'get_convective_envelope_mass')
        handler.add_setter('particles', 'set_convective_envelope_radius')
        handler.add_getter('particles', 'get_convective_envelope_radius')
        handler.add_setter('particles', 'set_luminosity')
        handler.add_getter('particles', 'get_luminosity')

        handler.add_setter('particles', 'set_check_for_secular_breakdown')
        handler.add_getter('particles', 'get_check_for_secular_breakdown')
        handler.add_setter('particles', 'set_check_for_dynamical_instability')
        handler.add_getter('particles', 'get_check_for_dynamical_instability')
        handler.add_setter('particles', 'set_dynamical_instability_criterion')
        handler.add_getter('particles', 'get_dynamical_instability_criterion')
        handler.add_setter('particles', 'set_dynamical_instability_central_particle')
        handler.add_getter('particles', 'get_dynamical_instability_central_particle')
        handler.add_setter('particles', 'set_dynamical_instability_K_parameter')
        handler.add_getter('particles', 'get_dynamical_instability_K_parameter')

        handler.add_setter('particles', 'set_check_for_physical_collision_or_orbit_crossing')
        handler.add_getter('particles', 'get_check_for_physical_collision_or_orbit_crossing')
        
        handler.add_setter('particles', 'set_check_for_minimum_periapse_distance')
        handler.add_getter('particles', 'get_check_for_minimum_periapse_distance')
        handler.add_setter('particles', 'set_check_for_minimum_periapse_distance_value')
        handler.add_getter('particles', 'get_check_for_minimum_periapse_distance_value')
        
        handler.add_setter('particles', 'set_check_for_RLOF_at_pericentre')
        handler.add_getter('particles', 'get_check_for_RLOF_at_pericentre')
        handler.add_setter('particles', 'set_check_for_RLOF_at_pericentre_use_sepinsky_fit')
        handler.add_getter('particles', 'get_check_for_RLOF_at_pericentre_use_sepinsky_fit')
        
        handler.add_setter('particles', 'set_root_finding_state')
        handler.add_getter('particles', 'get_root_finding_state')


        handler.define_set('external_particles', 'index_of_the_external_particle')
        handler.set_new   ('external_particles', 'new_external_particle')
        handler.set_delete('external_particles', 'delete_external_particle')
        
        handler.add_setter('external_particles', 'set_external_mass')
        handler.add_getter('external_particles', 'get_external_mass')
        handler.add_setter('external_particles', 'set_external_path')
        handler.add_getter('external_particles', 'get_external_path')
        handler.add_setter('external_particles', 'set_external_mode')
        handler.add_getter('external_particles', 'get_external_mode')

        handler.add_setter('external_particles', 'set_external_t_ref')
        handler.add_getter('external_particles', 'get_external_t_ref')
        handler.add_setter('external_particles', 'set_external_t_passed')
        handler.add_getter('external_particles', 'get_external_t_passed')

        handler.add_setter('external_particles', 'set_external_r0_vectors')
        handler.add_getter('external_particles', 'get_external_r0_vectors')        
        handler.add_setter('external_particles', 'set_external_rdot_vectors', names = ('rdot_vec_x','rdot_vec_y','rdot_vec_z') )
        handler.add_getter('external_particles', 'get_external_rdot_vectors', names = ('rdot_vec_x','rdot_vec_y','rdot_vec_z') )
        
        handler.add_setter('external_particles', 'set_external_eccentricity')
        handler.add_getter('external_particles', 'get_external_eccentricity')
        handler.add_setter('external_particles', 'set_external_periapse_distance')
        handler.add_getter('external_particles', 'get_external_periapse_distance')

        handler.add_setter('external_particles', 'set_external_e_hat_vectors')
        handler.add_getter('external_particles', 'get_external_e_hat_vectors')        
        handler.add_setter('external_particles', 'set_external_h_hat_vectors')
        handler.add_getter('external_particles', 'get_external_h_hat_vectors')        
        
        handler.add_getter('external_particles', 'get_external_r_vectors')        
        
    def define_parameters(self, handler):
        
        handler.add_method_parameter(
            "get_relative_tolerance",
            "set_relative_tolerance",
            "relative_tolerance",
            "relative_tolerance",
            default_value = 1.0e-16
        )
        handler.add_method_parameter(
            "get_absolute_tolerance_eccentricity_vectors",
            "set_absolute_tolerance_eccentricity_vectors",
            "absolute_tolerance_eccentricity_vectors",
            "absolute_tolerance_eccentricity_vectors",
            default_value = 1.0e-14
        )
        handler.add_method_parameter(
            "get_include_quadrupole_order_terms",
            "set_include_quadrupole_order_terms",
            "include_quadrupole_order_terms",
            "include_quadrupole_order_terms",
            default_value = True
        )
        handler.add_method_parameter(
            "get_include_octupole_order_binary_pair_terms",
            "set_include_octupole_order_binary_pair_terms",
            "include_octupole_order_binary_pair_terms",
            "include_octupole_order_binary_pair_terms",
            default_value = True
        )
        handler.add_method_parameter(
            "get_include_octupole_order_binary_triplet_terms",
            "set_include_octupole_order_binary_triplet_terms",
            "include_octupole_order_binary_triplet_terms",
            "include_octupole_order_binary_triplet_terms",
            default_value = False
        )
        handler.add_method_parameter(
            "get_include_hexadecupole_order_binary_pair_terms",
            "set_include_hexadecupole_order_binary_pair_terms",
            "include_hexadecupole_order_binary_pair_terms",
            "include_hexadecupole_order_binary_pair_terms",
            default_value = False
        )
        handler.add_method_parameter(
            "get_include_dotriacontupole_order_binary_pair_terms",
            "set_include_dotriacontupole_order_binary_pair_terms",
            "include_dotriacontupole_order_binary_pair_terms",
            "include_dotriacontupole_order_binary_pair_terms",
            default_value = False
        )
        handler.add_method_parameter(
            "get_orbital_phases_random_seed",
            "set_orbital_phases_random_seed",
            "orbital_phases_random_seed",
            "orbital_phases_random_seed",
            default_value = 0
        )


    def define_methods(self, handler):
        pass

    def before_get_parameter(self):
        """
        Called everytime just before a parameter is retrieved in using::
            instance.parameter.name
        """
        pass
        
    def before_set_parameter(self):
        """
        Called everytime just before a parameter is updated in using::
            instance.parameter.name = newvalue
        """
        pass

    def commit_particles(self):
        print(print_name,' -- committing particles')
        particles = self.particles

        if len(particles) == 0:
            print(print_name,' -- no particles have been added -- exiting')
            exit(-1)

        particles.add_vector_attribute("spin_vec",["spin_vec_x","spin_vec_y","spin_vec_z"])
        particles.add_vector_attribute("e_vec",["e_vec_x","e_vec_y","e_vec_z"])
        particles.add_vector_attribute("h_vec",["h_vec_x","h_vec_y","h_vec_z"])
        
        self.external_particles.add_vector_attribute("r_vec",["r_vec_x","r_vec_y","r_vec_z"])

        ### evaluate the initial hamiltonian ###
        time_step = 0.0 | units.Myr 
        end_time,self.initial_hamiltonian,flag,error_code = self.evolve_interface(self.model_time,time_step)


        self.particles_committed = True
        
    def evolve_model(self,end_time):
        if end_time is None:
            print(print_name,' -- end time not specified in evolve_model! exiting')
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
            print(print_name,' -- error occurred during ODE integration')
            print(print_name,' -- error code is ',error_code)
        self.flag = flag
        self.error_code = error_code

    def determine_binary_parents_levels_and_masses(self):
        self.determine_binary_parents_levels_and_masses_interface()

    def apply_external_perturbation_assuming_integrated_orbits(self):
        self.apply_external_perturbation_assuming_integrated_orbits_interface()

    def apply_user_specified_instantaneous_perturbation(self):
        self.apply_user_specified_instantaneous_perturbation_interface()

    def set_positions_and_velocities(self):
        self.set_positions_and_velocities_interface()


Secularmultiple = SecularMultiple
