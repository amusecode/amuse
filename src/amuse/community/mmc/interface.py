from amuse.community import *
from amuse.support.options import OptionalAttributes, option
import os
class mmcInterface(LegacyInterface):
    
    use_modules = ['MMC']
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, **keyword_arguments)

    @option(type="string")
    def data_directory(self):
        """
        The root name of the directory for the EVTwin
        application data files. This directory should contain the
        zams data and init.run and init.dat.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'mmc')

    @legacy_function
    def nonstandard_init():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def commit_parameters():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mmc_data_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('data_directory', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT)
        function.addParameter('r', dtype='float64', direction=function.OUT)
        function.addParameter('vr', dtype='float64', direction=function.OUT)


        function.addParameter('vt', dtype='float64', direction=function.OUT)
        function.addParameter('x1', dtype='float64', direction=function.OUT)
        function.addParameter('x2', dtype='float64', direction=function.OUT)
        function.addParameter('x3', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function    

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()  
        function.addParameter('Ek', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def initial_run():
        function = LegacyFunctionSpecification()  
        function.addParameter('res', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def run():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function    
    
    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function    

    @legacy_function
    def set_nitesc():
        function = LegacyFunctionSpecification()
        function.addParameter('iteration_tidal_radius', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iflagbh():
        function = LegacyFunctionSpecification()
        function.addParameter('natal_kicks_bh', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iflagns():
        function = LegacyFunctionSpecification()
        function.addParameter('natal_kicks_ns', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ikroupa():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_bins_parameters', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_zini():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_metalicity', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ybmin():
        function = LegacyFunctionSpecification()
        function.addParameter('mult_bmin0', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ytau():
        function = LegacyFunctionSpecification()
        function.addParameter('mult_tau0', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_xtau():
        function = LegacyFunctionSpecification()
        function.addParameter('call_mloss', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tcrevo():
        function = LegacyFunctionSpecification()
        function.addParameter('critical_time_step_dtte0_to_dtte', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dtte_naught():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_mloss_call_tph_lt_tcr', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dtte():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_mloss_call', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dttp():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_profile_output', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rplum():
        function = LegacyFunctionSpecification()
        function.addParameter('rsplum_scale_radius_plummer_model', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_xtid():
        function = LegacyFunctionSpecification()
        function.addParameter('coeff_front_tidal_energy', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_gamma():
        function = LegacyFunctionSpecification()
        function.addParameter('param_in_coulomb_log', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tau_naught():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_complete_cluster_model', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bmax():
        function = LegacyFunctionSpecification()
        function.addParameter('max_val_of_sin_betasqr', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bmin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_val_of_sin_betasqr', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_w_naught():
        function = LegacyFunctionSpecification()
        function.addParameter('king_model_parameter', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_zmbar():
        function = LegacyFunctionSpecification()
        function.addParameter('total_mass_cluster', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rbar():
        function = LegacyFunctionSpecification()
        function.addParameter('tidal_radius', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_qvir():
        function = LegacyFunctionSpecification()
        function.addParameter('virial_ratio', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_amax():
        function = LegacyFunctionSpecification()
        function.addParameter('max_semi_major_ax_of_bins', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_amin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_semi_major_ax_of_bins', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_body_fracb():
        function = LegacyFunctionSpecification()
        function.addParameter('primordial_bin_fraction', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_body_n():
        function = LegacyFunctionSpecification()
        function.addParameter('min_particle_mass_before_scaling', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_body_1():
        function = LegacyFunctionSpecification()
        function.addParameter('max_particle_mass_before_scaling', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_brakem():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_in_which_IMF_breaks', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_alphah():
        function = LegacyFunctionSpecification()
        function.addParameter('pwr_law_index_ht_breake_mass', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_alphal():
        function = LegacyFunctionSpecification()
        function.addParameter('pwr_law_index_lt_breake_mass', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_qe():
        function = LegacyFunctionSpecification()
        function.addParameter('energy_tolerance', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tcomp():
        function = LegacyFunctionSpecification()
        function.addParameter('max_comp_time_in_hours', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tcrit():
        function = LegacyFunctionSpecification()
        function.addParameter('termination_time_units_crossing_time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iexch():
        function = LegacyFunctionSpecification()
        function.addParameter('exchange_mode', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ib3f():
        function = LegacyFunctionSpecification()
        function.addParameter('spitzer_or_heggie', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iprint():
        function = LegacyFunctionSpecification()
        function.addParameter('diag_info', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_imodel():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_model', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ntwo():
        function = LegacyFunctionSpecification()
        function.addParameter('max_index_od_two', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nminzo():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_stars_in_a_zone', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nzonc():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_zones_in_the_core', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nz_naught():
        function = LegacyFunctionSpecification()
        function.addParameter('numb_of_stars_in_each_zone_at_t0', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nmin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_stars_to_calc_c_parms', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ncor():
        function = LegacyFunctionSpecification()
        function.addParameter('numb_of_stars_to_calc_c_parms', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_istart():
        function = LegacyFunctionSpecification()
        function.addParameter('start_or_restart', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nt():
        function = LegacyFunctionSpecification()
        function.addParameter('tot_numb_of_objs', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_irun():
        function = LegacyFunctionSpecification()
        function.addParameter('init_sequence_of_rnd_numbs', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
        
    
    
class mmc(InCodeComponentImplementation):

    def __init__(self):
        InCodeComponentImplementation.__init__(self,  NearestNeighborInterface())
    
