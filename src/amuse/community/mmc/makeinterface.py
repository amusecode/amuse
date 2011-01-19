funcs =\
[['irun','','init_sequence_of_rnd_numbs'],
['nt','','tot_numb_of_objs'],
['istart','','start_or_restart'],
['ncor','','numb_of_stars_to_calc_c_parms'],
['nmin','','min_numb_of_stars_to_calc_c_parms'],
['nz_naught','','numb_of_stars_in_each_zone_at_t0'],
['nzonc','','min_numb_of_zones_in_the_core'],
['nminzo','','min_numb_of_stars_in_a_zone'],
['ntwo','','max_index_od_two'],
['imodel','','initial_model'],
['iprint','','diag_info'],
['ib3f','','spitzer_or_heggie'],
['iexch','','exchange_mode'],
['tcrit','DOUBLE PRECISION','termination_time_units_crossing_time'],
['tcomp','DOUBLE PRECISION','max_comp_time_in_hours'],
['qe', 'DOUBLE PRECISION','energy_tolerance'],
['alphal','','pwr_law_index_lt_breake_mass'],
['alphah', '','pwr_law_index_ht_breake_mass'],
['brakem','DOUBLE PRECISION','mass_in_which_IMF_breaks'],
['body_1','DOUBLE PRECISION','max_particle_mass_before_scaling'],
['body_n','DOUBLE PRECISION','min_particle_mass_before_scaling'],
['body_fracb','DOUBLE PRECISION','primordial_bin_fraction'],
['amin','DOUBLE PRECISION','min_semi_major_ax_of_bins'],
['amax','DOUBLE PRECISION','max_semi_major_ax_of_bins'],
['qvir','DOUBLE PRECISION','virial_ratio'],
['rbar','DOUBLE PRECISION','tidal_radius'],
['zmbar','DOUBLE PRECISION','total_mass_cluster'],
['w_naught','DOUBLE PRECISION','king_model_parameter'],
['bmin','DOUBLE PRECISION','min_val_of_sin_betasqr'],
['bmax','DOUBLE PRECISION','max_val_of_sin_betasqr'],
['tau_naught','DOUBLE PRECISION','time_step_for_complete_cluster_model'],
['gamma','DOUBLE PRECISION','param_in_coulomb_log'],
['xtid','DOUBLE PRECISION','coeff_front_tidal_energy'],
['rplum','DOUBLE PRECISION','rsplum_scale_radius_plummer_model'],
['dttp','DOUBLE PRECISION','time_step_for_profile_output'],
['dtte','DOUBLE PRECISION','time_step_for_mloss_call'],
['dtte_naught','DOUBLE PRECISION','time_step_for_mloss_call_tph_lt_tcr'],
['tcrevo','DOUBLE PRECISION','critical_time_step_dtte0_to_dtte'],
['xtau','DOUBLE PRECISION','call_mloss'],
['ytau','DOUBLE PRECISION','mult_tau0'],
['ybmin','DOUBLE PRECISION','mult_bmin0'],
['zini','DOUBLE PRECISION','initial_metalicity'],
['ikroupa','','initial_bins_parameters'],
['iflagns','','natal_kicks_ns'],
['iflagbh','','natal_kicks_bh'],
['nitesc','','iteration_tidal_radius']]

codestringfort = \
"""
MODULE MMC

CONTAINS

FUNCTION nonstandard_init()
  INTEGER :: nonstandard_init
  INTEGER :: init_sequence
  INTEGER :: res
  ! read initial parameters
  PRINT*,'calling input'
  res = init_sequence()
  nonstandard_init = res
  PRINT*,'init done'
END FUNCTION

FUNCTION set_mmc_data_directory(data_directory) 
  INTEGER :: set_mmc_data_directory
  CHARACTER(len=200) :: data_directory
  CALL amuse_set_mmc_data_directory(data_directory)
  set_mmc_data_directory = 0
END FUNCTION

FUNCTION get_time(time_)
  INTEGER :: get_time
  COMMON /SYSTEM/ TIME
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: time_
  time_ = time
  get_time = 0
END FUNCTION

FUNCTION set_time(time_)
  INTEGER :: set_time
  COMMON /SYSTEM/ TIME
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: time_
  time = time_
  set_time = 0
END FUNCTION


FUNCTION get_kinetic_energy(Ek)
  IMPLICIT NONE
  INTEGER :: res
  INTEGER :: total_kinetic_energy
  INTEGER :: get_kinetic_energy
  DOUBLE PRECISION :: Ek

  res = total_kinetic_energy(Ek)
  get_kinetic_energy = 0
END FUNCTION


FUNCTION new_particle( index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: new_particle
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  new_particle=0
END FUNCTION

FUNCTION delete_particle( index_of_the_particle)
  IMPLICIT NONE
  INTEGER :: delete_particle
  INTEGER :: index_of_the_particle
  delete_particle=0
END FUNCTION

FUNCTION set_state( index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: set_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  set_state=0
END FUNCTION

FUNCTION get_state( index_of_the_particle, mass_, r_, vr_, vt_,x1,x2,x3)
  INTEGER :: get_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: mass_, r_, vr_, vt_, x1,x2,x3
  
  COMMON /BODY/ BODY(1000),VR(1000),VT(1000), &
                U(1000),XESCP(1000),XESCT(1000),VRR(1000),R(1000)

  COMMON /POSVEL/ X(1000,3), XDOT(1000,3)
  REAL*8  BODY, VR, VT, R
  REAL*8  X,XDOT

  mass_ = BODY(index_of_the_particle)
  r_ = R(index_of_the_particle)
  vr_ = VR(index_of_the_particle)
  vt_ = VT(index_of_the_particle)
  x1 = X(index_of_the_particle,1)
  y1 = X(index_of_the_particle,2)
  z1 = X(index_of_the_particle,3)

  get_state=0
END FUNCTION

FUNCTION get_number_of_particles( VALUE)
  IMPLICIT NONE
  INTEGER :: get_number_of_particles
  INTEGER :: VALUE
  get_number_of_particles=0
END FUNCTION

FUNCTION internal__redirect_outputs( stdoutfile, stderrfile)
  IMPLICIT NONE
  INTEGER :: internal__redirect_outputs
  CHARACTER(LEN=*) :: stdoutfile, stderrfile
  internal__redirect_outputs=0
END FUNCTION

FUNCTION run( )
  IMPLICIT NONE
  INTEGER :: run
  run=0
END FUNCTION

FUNCTION commit_parameters( )
  IMPLICIT NONE
  INTEGER commit_parameters
  commit_parameters=0
END FUNCTION

END MODULE
"""

codestringpyth = \
"""
from amuse.community import *
from amuse.support.options import OptionalAttributes, option
import os
class mmcInterface(LegacyInterface):
    
    use_modules = ['MMC']
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, **keyword_arguments)

    @option(type="string")
    def data_directory(self):
        \"""
        The root name of the directory for the EVTwin
        application data files. This directory should contain the
        zams data and init.run and init.dat.
        \"""
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
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
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
        
    
    
class mmc(CodeInterface):

    def __init__(self):
        CodeInterface.__init__(self,  NearestNeighborInterface())
    

"""

if __name__ == '__main__':
    f = open('interface.f90','w')

    f.write(codestringfort)

    for O in funcs:
        i=O[0]
        tp = O[1]
        f.write( "FUNCTION set_{0}(tmp_)\n".format(i))
        f.write( "  COMMON /IPARAM/ {0}\n".format(i))


        if tp == '':
            f.write("  INTEGER {0}\n".format(i))
            f.write("  INTEGER tmp_\n")
        else:
            f.write( "  {0} {1}\n".format(tp,i))
            f.write( "  {0} tmp_\n".format(tp))

        f.write( "  INTEGER set_{0}\n".format(i))
        f.write( "  {0} = tmp_\n".format(i))
        f.write( "  set_{0} = 0\n".format(i))
        f.write( "END FUNCTION set_{0}\n".format(i))
        f.write('\n')

    f.close()
    
    print codestringpyth
