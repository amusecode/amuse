from amuse.community import *
from amuse.community.mmc.amuselib.interface import supportInterface
from amuse.support.options import OptionalAttributes, option
from amuse.ext.polarsupport import PolarSupport
import os
import numpy as np

class mmcInterface(CodeInterface, PolarSupport):

    use_modules = ['MMC']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, **keyword_arguments)

    @option(type="string")
    def data_directory(self):
        """
        The root name of the directory for the mmc
        application data files.
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
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass_', dtype='float64', direction=function.IN)
        function.addParameter('r_', dtype='float64', direction=function.IN)
        function.addParameter('vr_', dtype='float64', direction=function.IN)
        function.addParameter('vt_', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def commit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function
    
    @legacy_function
    def recommit_particles():
        """
        Let the code perform initialization actions
        after the number of particles have been updated
        or particle attributes have been updated from
        the script.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function
    
    @legacy_function
    def test_sort_routine():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('n', dtype='int32', direction=function.LENGTH)
        function.addParameter('aa', dtype='float64', direction=function.INOUT)
        function.addParameter('bb', dtype='int32', direction=function.INOUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_crossing_time():
        function = LegacyFunctionSpecification()
        function.addParameter('tcr', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_timet():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def run_a_while():
        function = LegacyFunctionSpecification()
        function.addParameter('iphase', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('time_end', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def amuse_input():
        function = LegacyFunctionSpecification()
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
        function.addParameter('n_', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_irun():
        function = LegacyFunctionSpecification()
        function.addParameter('init_sequence_of_rnd_numbs', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iseed():
        function = LegacyFunctionSpecification()
        function.addParameter('init_sequence_of_rnd_numbs', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_flagr():
        function = LegacyFunctionSpecification()
        function.addParameter('lagrangeradii', dtype='float64', direction=function.IN)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_flagr():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('flagr_', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def get_nlagra():
        function = LegacyFunctionSpecification()
        function.addParameter('nlagrange', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_irun():
        function = LegacyFunctionSpecification()
        function.addParameter('init_sequence_of_rnd_numbs', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nt():
        function = LegacyFunctionSpecification()
        function.addParameter('tot_numb_of_objs', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nt0():
        function = LegacyFunctionSpecification()
        function.addParameter('tot_numb_of_objs', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nt00():
        function = LegacyFunctionSpecification()
        function.addParameter('tot_numb_of_objs', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nt():
        function = LegacyFunctionSpecification()
        function.addParameter('tot_numb_of_objs', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_istart():
        function = LegacyFunctionSpecification()
        function.addParameter('start_or_restart', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_istart():
        function = LegacyFunctionSpecification()
        function.addParameter('start_or_restart', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ncor():
        function = LegacyFunctionSpecification()
        function.addParameter('numb_of_stars_to_calc_c_parms', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ncor():
        function = LegacyFunctionSpecification()
        function.addParameter('numb_of_stars_to_calc_c_parms', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nmin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_stars_to_calc_c_parms', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nmin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_stars_to_calc_c_parms', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nz0():
        function = LegacyFunctionSpecification()
        function.addParameter('numb_of_stars_in_each_zone_at_t0', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nz0():
        function = LegacyFunctionSpecification()
        function.addParameter('numb_of_stars_in_each_zone_at_t0', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nzonc():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_zones_in_the_core', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nzonc():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_zones_in_the_core', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nminzo():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_stars_in_a_zone', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nminzo():
        function = LegacyFunctionSpecification()
        function.addParameter('min_numb_of_stars_in_a_zone', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ntwo():
        function = LegacyFunctionSpecification()
        function.addParameter('max_index_od_two', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ntwo():
        function = LegacyFunctionSpecification()
        function.addParameter('max_index_od_two', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_imodel():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_model', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_imodel():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_model', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iprint():
        function = LegacyFunctionSpecification()
        function.addParameter('diag_info', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_iprint():
        function = LegacyFunctionSpecification()
        function.addParameter('diag_info', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ib3f():
        function = LegacyFunctionSpecification()
        function.addParameter('spitzer_or_heggie', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ib3f():
        function = LegacyFunctionSpecification()
        function.addParameter('spitzer_or_heggie', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iexch():
        function = LegacyFunctionSpecification()
        function.addParameter('exchange_mode', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_iexch():
        function = LegacyFunctionSpecification()
        function.addParameter('exchange_mode', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tcrit():
        function = LegacyFunctionSpecification()
        function.addParameter('termination_time_units_crossing_time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tcrit():
        function = LegacyFunctionSpecification()
        function.addParameter('termination_time_units_crossing_time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tcomp():
        function = LegacyFunctionSpecification()
        function.addParameter('max_comp_time_in_hours', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tcomp():
        function = LegacyFunctionSpecification()
        function.addParameter('max_comp_time_in_hours', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_qe():
        function = LegacyFunctionSpecification()
        function.addParameter('energy_tolerance', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_qe():
        function = LegacyFunctionSpecification()
        function.addParameter('energy_tolerance', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_alphal():
        function = LegacyFunctionSpecification()
        function.addParameter('pwr_law_index_lt_breake_mass', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_alphal():
        function = LegacyFunctionSpecification()
        function.addParameter('pwr_law_index_lt_breake_mass', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_alphah():
        function = LegacyFunctionSpecification()
        function.addParameter('pwr_law_index_ht_breake_mass', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_alphah():
        function = LegacyFunctionSpecification()
        function.addParameter('pwr_law_index_ht_breake_mass', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_brakem():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_in_which_IMF_breaks', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_brakem():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_in_which_IMF_breaks', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_body1():
        function = LegacyFunctionSpecification()
        function.addParameter('max_particle_mass_before_scaling', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_body1():
        function = LegacyFunctionSpecification()
        function.addParameter('max_particle_mass_before_scaling', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bodyn():
        function = LegacyFunctionSpecification()
        function.addParameter('min_particle_mass_before_scaling', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_bodyn():
        function = LegacyFunctionSpecification()
        function.addParameter('min_particle_mass_before_scaling', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_fracb():
        function = LegacyFunctionSpecification()
        function.addParameter('primordial_bin_fraction', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_fracb():
        function = LegacyFunctionSpecification()
        function.addParameter('primordial_bin_fraction', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_amin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_semi_major_ax_of_bins', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_amin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_semi_major_ax_of_bins', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_amax():
        function = LegacyFunctionSpecification()
        function.addParameter('max_semi_major_ax_of_bins', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_amax():
        function = LegacyFunctionSpecification()
        function.addParameter('max_semi_major_ax_of_bins', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_qvir():
        function = LegacyFunctionSpecification()
        function.addParameter('virial_ratio', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_qvir():
        function = LegacyFunctionSpecification()
        function.addParameter('virial_ratio', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rbar():
        function = LegacyFunctionSpecification()
        function.addParameter('tidal_radius', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_rbar():
        function = LegacyFunctionSpecification()
        function.addParameter('tidal_radius', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_zmbar():
        function = LegacyFunctionSpecification()
        function.addParameter('total_mass_cluster', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_zmbar():
        function = LegacyFunctionSpecification()
        function.addParameter('total_mass_cluster', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_w0():
        function = LegacyFunctionSpecification()
        function.addParameter('king_model_parameter', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_w0():
        function = LegacyFunctionSpecification()
        function.addParameter('king_model_parameter', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bmin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_val_of_sin_betasqr', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bmin0():
        function = LegacyFunctionSpecification()
        function.addParameter('min_val_of_sin_betasqr', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_bmin():
        function = LegacyFunctionSpecification()
        function.addParameter('min_val_of_sin_betasqr', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bmax():
        function = LegacyFunctionSpecification()
        function.addParameter('max_val_of_sin_betasqr', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_bmax():
        function = LegacyFunctionSpecification()
        function.addParameter('max_val_of_sin_betasqr', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tau0():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_complete_cluster_model', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tau0():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_complete_cluster_model', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_gamma():
        function = LegacyFunctionSpecification()
        function.addParameter('param_in_coulomb_log', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gamma():
        function = LegacyFunctionSpecification()
        function.addParameter('param_in_coulomb_log', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_xtid():
        function = LegacyFunctionSpecification()
        function.addParameter('coeff_front_tidal_energy', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_xtid():
        function = LegacyFunctionSpecification()
        function.addParameter('coeff_front_tidal_energy', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rplum():
        function = LegacyFunctionSpecification()
        function.addParameter('rsplum_scale_radius_plummer_model', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_rplum():
        function = LegacyFunctionSpecification()
        function.addParameter('rsplum_scale_radius_plummer_model', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dttp():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_profile_output', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dttp():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_profile_output', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dtte():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_mloss_call', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dtte():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_mloss_call', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dtte0():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_mloss_call_tph_lt_tcr', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dtte0():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step_for_mloss_call_tph_lt_tcr', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tcrevo():
        function = LegacyFunctionSpecification()
        function.addParameter('critical_time_step_dtte0_to_dtte', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tcrevo():
        function = LegacyFunctionSpecification()
        function.addParameter('critical_time_step_dtte0_to_dtte', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_xtau():
        function = LegacyFunctionSpecification()
        function.addParameter('call_mloss', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_xtau():
        function = LegacyFunctionSpecification()
        function.addParameter('call_mloss', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ytau():
        function = LegacyFunctionSpecification()
        function.addParameter('mult_tau0', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ytau():
        function = LegacyFunctionSpecification()
        function.addParameter('mult_tau0', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ybmin():
        function = LegacyFunctionSpecification()
        function.addParameter('mult_bmin0', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ybmin():
        function = LegacyFunctionSpecification()
        function.addParameter('mult_bmin0', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_zini():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_metalicity', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_zini():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_metalicity', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ikroupa():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_bins_parameters', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ikroupa():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_bins_parameters', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iflagns():
        function = LegacyFunctionSpecification()
        function.addParameter('natal_kicks_ns', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_iflagns():
        function = LegacyFunctionSpecification()
        function.addParameter('natal_kicks_ns', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iflagbh():
        function = LegacyFunctionSpecification()
        function.addParameter('natal_kicks_bh', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_iflagbh():
        function = LegacyFunctionSpecification()
        function.addParameter('natal_kicks_bh', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nitesc():
        function = LegacyFunctionSpecification()
        function.addParameter('iteration_tidal_radius', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nitesc():
        function = LegacyFunctionSpecification()
        function.addParameter('iteration_tidal_radius', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_rc():
        function = LegacyFunctionSpecification()
        function.addParameter('coreradius', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nc():
        function = LegacyFunctionSpecification()
        function.addParameter('numberdensity', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def call_zone():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def call_relaxt():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def call_amuse_output():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    def get_positions_3d_incode(self, particles):
        #this is an experiment!
        r = self.get_state(particles).r
        i = supportInterface()
        x = np.zeros(len(r))
        y = np.zeros(len(r))
        z = np.zeros(len(r))
        R = i.rnd_points_on_sphere(x, y, z)
        #R = i.many_points_on_sphere(x, y, z)
        i.stop()
        return {'x':np.array(r*R.x), 'y':np.array(r*R.y),'z':np.array(r*R.z)}

    def get_positions_3d(self, particles):
        r = self.get_state(particles).r
        return self.position_to_cartesian(r)

class mmc(InCodeComponentImplementation):

    def __init__(self):
        InCodeComponentImplementation.__init__(self,  mmcInterface())
