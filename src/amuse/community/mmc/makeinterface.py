"""
   this is a tmp script, after all fctns are in place etc..
   I'll remove this
"""

funcs =\
[['irun','','init_sequence_of_rnd_numbs'],
['nt','','tot_numb_of_objs'],
['nt_naught','','num'],
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
['tau0','DOUBLE PRECISION','time_step_for_complete_cluster_model'],
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

FUNCTION evolve(time_end)
  INTEGER :: evolve
  INTEGER :: evolve_src
  INTEGER :: res
  DOUBLE PRECISION time_end
  res = evolve_src(time_end)
  evolve = res
END FUNCTION  

FUNCTION get_time(time_)
  INTEGER :: get_time
  INTEGER :: get_time_src
  INTEGER :: res
  DOUBLE PRECISION :: time_
  res = get_time_src(time_)
  get_time = res
END FUNCTION

FUNCTION get_timet(time_)
  INTEGER :: get_timet
  INTEGER :: get_timet_src
  INTEGER :: res
  DOUBLE PRECISION :: time_
  res = get_timet_src()
  get_timet = res
END FUNCTION

FUNCTION get_crossing_time(time_)
  INTEGER :: get_crossing_time
  INTEGER :: get_crossing_time_src
  INTEGER :: res
  DOUBLE PRECISION :: time_
  res = get_crossing_time_src()
  get_crossing_time = res
END FUNCTION

FUNCTION set_time(time_)
  INTEGER :: set_time
  COMMON /SYSTEM/ TIME
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: time_
  time = time_
  set_time = 0
END FUNCTION

FUNCTION get_number_of_particles(n_)
  INTEGER :: get_number_of_particles
  INTEGER :: get_number_of_particles_, res
  INTEGER :: n_
  res = get_number_of_particles_src(n_)
  print*, "n call = ",n_
  get_number_of_particles=0
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

FUNCTION set_state( index_of_the_particle, mass_, r_, vr_, vt_)
  IMPLICIT NONE
  INTEGER :: set_state, set_state_src
  INTEGER :: index_of_the_particle
  INTEGER :: res
  DOUBLE PRECISION :: mass_, r_, vr_, vt_
  res = set_state_src(index_of_the_particle, mass_, r_, vr_, vt_)
  set_state=0
END FUNCTION

FUNCTION get_state( index_of_the_particle, mass_, r_, vr_, vt_)
  INTEGER :: get_state
  INTEGER :: get_state_src,res
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: mass_, r_, vr_, vt_
  
  res = get_state_src(index_of_the_particle, mass_,r_,vr_,vt_)
  get_state=0

END FUNCTION

FUNCTION get_total_kinetic_energy(T)
  INTEGER :: get_total_kinetic_energy
  INTEGER :: get_total_kinetic_energy_src
  INTEGER :: res
  DOUBLE PRECISION :: T

  res = get_total_kinetic_energy_src(T)
  get_total_kinetic_energy = 0

END FUNCTION

FUNCTION get_total_potential_energy(V)
  INTEGER :: get_total_potential_energy
  INTEGER :: get_total_potential_energy_src
  INTEGER :: res
  DOUBLE PRECISION :: V

  res = get_total_potential_energy_src(T)
  get_total_potential_energy = 0

END FUNCTION

FUNCTION commit_particles()
  CALL commit_particles_src()
END FUNCTION 

FUNCTION recommit_particles()
  CALL recommit_particles_src()
END FUNCTION 

FUNCTION test_sort_routine(n, aa, bb)
  INTEGER n
  REAL*8 aa
  INTEGER bb
  INTEGER test_sort_routine
  INTEGER res
  dimension aa(n), bb(n)

  res = test_sort(n, aa, bb)
  test_sort_routine = 0
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

FUNCTION amuse_input()
  INTEGER amuse_input
  call amuse_input_src()
  amuse_input = 0
END FUNCTION

FUNCTION commit_parameters( )
  IMPLICIT NONE
  INTEGER commit_parameters
  commit_parameters=0
END FUNCTION

"""

if __name__ == '__main__':
    f = open('interface.f90','w')
    g = open('src/amuse_interface2.f','w')
    h = open('interface_specs.tmp','w')
    
    f.write(codestringfort)

    for O in funcs:
        i=O[0]
        tp = O[1]
        nm = O[2]
        
        f.write( "FUNCTION set_{0}(tmp_)\n".format(i))
        f.write( "  INTEGER :: set_{0}\n".format(i))
        if tp == '':
            #f.write("  INTEGER {0}\n".format(i))
            f.write("  INTEGER :: tmp_\n")
        else:
            #f.write( "  {0} {1}\n".format(tp,i))
            f.write( "  {0} :: tmp_\n".format(tp))
        f.write( "  INTEGER :: res\n")
        f.write( "  INTEGER :: set_{0}_src\n".format(i))
        f.write( "  res = set_{0}_src(tmp_)\n".format(i))
        f.write( "  set_{0} = res\n".format(i))
        f.write( "END FUNCTION\n".format(i))
        f.write('\n')

        f.write( "FUNCTION get_{0}(tmp_)\n".format(i))
        f.write( "  INTEGER :: get_{0}\n".format(i))
        if tp == '':
            f.write("  INTEGER :: tmp_\n")
        else:
            f.write( "  {0} :: tmp_\n".format(tp))
        f.write( "  INTEGER :: res\n")
        f.write( "  INTEGER :: get_{0}_src\n".format(i))
        f.write( "  res = get_{0}_src(tmp_)\n".format(i))
        f.write( "  get_{0} = res\n".format(i))
        f.write( "END FUNCTION\n".format(i))
        f.write('\n')

        g.write( "      FUNCTION set_{0}_src(tmp_)\n".format(i))
        g.write( "      INCLUDE 'common.h'\n")
        if tp == '':
            #g.write("      INTEGER {0}\n".format(i))
            g.write("      INTEGER tmp_\n")
        else:
            #g.write( "      {0} {1}\n".format(tp,i))
            g.write( "      {0} tmp_\n".format(tp))
        g.write( "      INTEGER set_{0}_src\n".format(i))
        g.write( "      {0} = tmp_\n".format(i))
        g.write( "      set_{0}_src = 0\n".format(i))
        g.write( "      END FUNCTION\n".format(i))
        g.write('\n')

        g.write( "      FUNCTION get_{0}_src(tmp_)\n".format(i))
        g.write( "      INCLUDE 'common.h'\n")
        if tp == '':
            g.write("      INTEGER tmp_\n")
        else:
            g.write( "      {0} tmp_\n".format(tp))
        g.write( "      INTEGER get_{0}_src\n".format(i))
        g.write( "      tmp_ = {0}\n".format(i))
        g.write( "      get_{0}_src = 0\n".format(i))
        g.write( "      END FUNCTION\n".format(i))
        g.write('\n')


        if tp == '':
            typus = 'int32'
        else:
            typus = 'float64'

        h.write("@legacy_function\n")
        h.write("def set_{0}():\n".format(i))
        h.write("    function = LegacyFunctionSpecification()\n")
        h.write("    function.addParameter('{0}', dtype='{1}', direction=function.IN)\n".format(nm, typus))
        h.write("    function.result_type = 'int32'\n")
        h.write("    return function\n\n")

        h.write("@legacy_function\n")
        h.write("def get_{0}():\n".format(i))
        h.write("    function = LegacyFunctionSpecification()\n")
        h.write("    function.addParameter('{0}', dtype='{1}', direction=function.OUT)\n".format(nm,typus))
        h.write("    function.result_type = 'int32'\n")
        h.write("    return function\n\n")

    f.write("END MODULE\n")

    f.close()
    g.close()
    h.close()
