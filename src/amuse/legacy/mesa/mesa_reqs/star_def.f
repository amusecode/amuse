! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   you should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module star_def

      use const_def, only: ln10
      use chem_def, only: num_chem_elements
      use net_def, only: num_categories
      use utils_def, only: integer_dict
      
      implicit none
      
            
      ! wind schemes
      integer, parameter :: no_automatic_wind = 0
      integer, parameter :: Reimers_wind = 1
      integer, parameter :: Blocker_wind = 2
      integer, parameter :: de_Jager_wind = 3
      integer, parameter :: Dutch_wind = 4
      integer, parameter :: Mattsson_wind = 5


      ! mesh cell types
      integer, parameter :: split_type = -1
      integer, parameter :: unchanged_type = 0
      integer, parameter :: merged_type = 1
      
  
      ! there are several options for how eps_grav is evaluated.
      integer, parameter :: eps_grav_Tds_formula = 1 ! from K&W, eqn 4.27
         ! convert eps_grav = - T Ds/Dt into the following:
         ! eps_grav = -cp T ((1 - grad_ad chiT) DlnT/Dt - grad_ad chiRho Dlnd/Dt)
      integer, parameter :: eps_grav_dE_formula = 2
         ! E is internal energy from equation of state.
         ! eps_grav = -DE/Dt + P/rho Dlnd/Dt [e.g., see K&W, 4.25]
         ! NOTE: this form seems to have worse numerical problems than the
         ! Tds form when timesteps get small.
      
      
      type history_node
         integer :: nvals, n_ivals
         double precision, pointer :: vals(:) ! (nvals)
         integer, pointer :: ivals(:) ! (n_ivals)
         double precision :: time ! list is in order of decreasing time
         type (history_node), pointer :: next ! null() for end of list
      end type history_node


      integer, parameter :: star_def_version = 14
      integer, parameter :: max_generations = 3
      integer, parameter :: strlen = 256
      integer, parameter :: net_name_len = strlen
      integer, parameter :: name_len = 32           
      integer, parameter :: max_num_mixing_regions = 100


      type star_info
         
         include "star_controls.dek"
         include "star_data.dek"
   
         ! handles
         integer :: eos_handle
         integer :: kap_handle
         integer :: net_handle
         integer :: burn_and_mix_handle
         integer :: hydro_handle
         
         ! private         
         include "private_controls.dek"
         integer :: retry_cnt
         integer :: dbg_control

         ! bookkeeping
         integer :: id
         logical :: in_use
         
      end type star_info
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! values for result_reason
      integer, parameter :: result_reason_normal = 1
      integer, parameter :: dt_is_zero = 2 
         ! indicates that t+dt == t, perhaps because of round-off with tiny dt.
      integer, parameter :: nonzero_ierr = 3 ! some routine returned with ierr /= 0
      integer, parameter :: hydro_failed_to_converge = 4
      integer, parameter :: hydro_err_too_large = 5
      integer, parameter :: burn_and_mix_failed = 6
      integer, parameter :: diffusion_failed = 7
      integer, parameter :: too_many_steps_for_burn = 8
      integer, parameter :: too_many_steps_for_diffusion = 9
      integer, parameter :: too_many_steps_for_hydro = 10
      integer, parameter :: adjust_mesh_failed = 11
      integer, parameter :: adjust_mass_failed = 12
      integer, parameter :: core_dump_model_number = 13
      integer, parameter :: timestep_limits = 14
      integer, parameter :: variable_change_limits = 15
      integer, parameter :: bad_lnd_prediction = 16
      integer, parameter :: bad_lnT_prediction = 17
      integer, parameter :: op_split_failed_to_converge = 18
      
      integer, parameter :: num_reasons = 18
      
      character (len=32) :: result_reason_str(num_reasons)


      integer, parameter :: max_star_handles = 1000 ! this can be increased as necessary
      type (star_info), target, save :: star_handles(max_star_handles) 
         ! gfortran requires "save" here. strange, but true.
         
         
      character (len=256) :: data_dir_for_mesa
         
         
      ! debugging storage
      integer, parameter :: max_ndbg = 9

      
      integer, parameter :: max_nvar = 200      


      ! phases of evolution for logs and profiles
      integer, parameter :: phase_starting = 0
      integer, parameter :: phase_early_main_seq = 1
      integer, parameter :: phase_mid_main_seq = 2
      integer, parameter :: phase_wait_for_he = 3
      integer, parameter :: phase_he_ignition_over = 4
      integer, parameter :: phase_he_igniting = 5
      integer, parameter :: phase_helium_burning = 6
      integer, parameter :: phase_carbon_burning = 7
      
      ! time_step limit identifiers
      integer, parameter :: Tlim_struc = 1
      integer, parameter :: Tlim_num_jacobians = 2
      integer, parameter :: Tlim_num_diff_solver_steps = 3
      integer, parameter :: Tlim_dX = 4
      integer, parameter :: Tlim_dH = 5
      integer, parameter :: Tlim_dHe = 6
      integer, parameter :: Tlim_dH_div_H = 7
      integer, parameter :: Tlim_dHe_div_He = 8
      integer, parameter :: Tlim_dX_div_X = 9
      integer, parameter :: Tlim_dL_div_L = 10
      integer, parameter :: Tlim_dlgP = 11
      integer, parameter :: Tlim_dlgRho = 12
      integer, parameter :: Tlim_dlgT = 13
      integer, parameter :: Tlim_dlgR = 14
      integer, parameter :: Tlim_dlgL_nuc_cat = 15
      integer, parameter :: Tlim_dlgL_H = 16
      integer, parameter :: Tlim_dlgL_He = 17
      integer, parameter :: Tlim_d_deltaR = 18
      integer, parameter :: Tlim_dlgL_z = 19
      integer, parameter :: Tlim_dlgL_nuc = 20
      integer, parameter :: Tlim_dlgTeff = 21
      integer, parameter :: Tlim_dxdt_nuc = 22
      integer, parameter :: Tlim_dlgRho_cntr = 23
      integer, parameter :: Tlim_dlgT_cntr = 24
      integer, parameter :: Tlim_lg_XH_cntr = 25
      integer, parameter :: Tlim_dmstar = 26
      integer, parameter :: Tlim_lgL = 27
      integer, parameter :: Tlim_max_timestep = 28
      integer, parameter :: Tlim_timestep_hold = 29
      integer, parameter :: Tlim_lg_XHe_cntr = 30
      integer, parameter :: Tlim_h_bdy = 31
      integer, parameter :: Tlim_he_bdy = 32
      integer, parameter :: Tlim_h1_czb = 33
      
      integer, parameter :: numTlim = 33
      
      character (len=14) :: dt_why_str(numTlim) ! indicates the reson for the timestep choice
      character (len=100) :: dt_why_long_str(numTlim)
      
   
   ! log column options


      integer, parameter :: l_model_number = 1
      integer, parameter :: l_star_age = l_model_number + 1
      integer, parameter :: l_star_mass = l_star_age + 1
      integer, parameter :: l_star_mdot = l_star_mass + 1
      integer, parameter :: l_log_abs_mdot = l_star_mdot + 1
      integer, parameter :: l_time_step = l_log_abs_mdot + 1
      integer, parameter :: l_num_zones = l_time_step + 1
      integer, parameter :: l_conv_mx1_top = l_num_zones + 1
      integer, parameter :: l_conv_mx1_bot = l_conv_mx1_top + 1
      integer, parameter :: l_conv_mx2_top = l_conv_mx1_bot + 1
      integer, parameter :: l_conv_mx2_bot = l_conv_mx2_top + 1
      integer, parameter :: l_mx1_top = l_conv_mx2_bot + 1
      integer, parameter :: l_mx1_bot = l_mx1_top + 1
      integer, parameter :: l_mx2_top = l_mx1_bot + 1
      integer, parameter :: l_mx2_bot = l_mx2_top + 1
      integer, parameter :: l_mixing_regions = l_mx2_bot + 1
      integer, parameter :: l_epsnuc_M_1 = l_mixing_regions + 1
      integer, parameter :: l_epsnuc_M_2 = l_epsnuc_M_1 + 1
      integer, parameter :: l_epsnuc_M_3 = l_epsnuc_M_2 + 1
      integer, parameter :: l_epsnuc_M_4 = l_epsnuc_M_3 + 1
      integer, parameter :: l_epsnuc_M_5 = l_epsnuc_M_4 + 1
      integer, parameter :: l_epsnuc_M_6 = l_epsnuc_M_5 + 1
      integer, parameter :: l_epsnuc_M_7 = l_epsnuc_M_6 + 1
      integer, parameter :: l_epsnuc_M_8 = l_epsnuc_M_7 + 1
      integer, parameter :: l_h1_boundary_mass = l_epsnuc_M_8 + 1
      integer, parameter :: l_he4_boundary_mass = l_h1_boundary_mass + 1
      integer, parameter :: l_power_h_burn = l_he4_boundary_mass + 1
      integer, parameter :: l_power_he_burn = l_power_h_burn + 1
      integer, parameter :: l_log_center_T = l_power_he_burn + 1
      integer, parameter :: l_log_center_Rho = l_log_center_T + 1
      integer, parameter :: l_v_div_csound_surf = l_log_center_Rho + 1
      integer, parameter :: l_surface_accel_div_grav = l_v_div_csound_surf + 1
      integer, parameter :: l_log_dt = l_surface_accel_div_grav + 1
      integer, parameter :: l_log_LH = l_log_dt + 1
      integer, parameter :: l_log_LHe = l_log_LH + 1
      integer, parameter :: l_log_L = l_log_LHe + 1
      integer, parameter :: l_log_R = l_log_L + 1
      integer, parameter :: l_log_Teff = l_log_R + 1
      integer, parameter :: l_log_g = l_log_Teff + 1
      integer, parameter :: l_log_L_div_Ledd = l_log_g + 1
      integer, parameter :: l_num_retries = l_log_L_div_Ledd + 1
      integer, parameter :: l_num_backups = l_num_retries + 1
      integer, parameter :: l_h1_czb_mass = l_num_backups + 1
      integer, parameter :: l_surf_c12_minus_o16 = l_h1_czb_mass + 1

      integer, parameter :: l_log_center_P = l_surf_c12_minus_o16 + 1
      integer, parameter :: l_center_degeneracy = l_log_center_P + 1
      integer, parameter :: l_center_gamma = l_center_degeneracy + 1
      integer, parameter :: l_h1_boundary_radius = l_center_gamma + 1
      integer, parameter :: l_h1_boundary_lgT = l_h1_boundary_radius + 1
      integer, parameter :: l_h1_boundary_lgRho = l_h1_boundary_lgT + 1
      integer, parameter :: l_h1_boundary_L = l_h1_boundary_lgRho + 1
      integer, parameter :: l_h1_boundary_v = l_h1_boundary_L + 1
      integer, parameter :: l_he4_boundary_radius = l_h1_boundary_v + 1
      integer, parameter :: l_he4_boundary_lgT = l_he4_boundary_radius + 1
      integer, parameter :: l_he4_boundary_lgRho = l_he4_boundary_lgT + 1
      integer, parameter :: l_he4_boundary_L = l_he4_boundary_lgRho + 1
      integer, parameter :: l_he4_boundary_v = l_he4_boundary_L + 1
      integer, parameter :: l_c12_boundary_mass = l_he4_boundary_v + 1
      integer, parameter :: l_c12_boundary_radius = l_c12_boundary_mass + 1
      integer, parameter :: l_c12_boundary_lgT = l_c12_boundary_radius + 1
      integer, parameter :: l_c12_boundary_lgRho = l_c12_boundary_lgT + 1
      integer, parameter :: l_c12_boundary_L = l_c12_boundary_lgRho + 1
      integer, parameter :: l_c12_boundary_v = l_c12_boundary_L + 1
      integer, parameter :: l_envelope_mass = l_c12_boundary_v + 1
      integer, parameter :: l_envelope_fraction_left = l_envelope_mass + 1
      
      integer, parameter :: l_tau10_mass = l_envelope_fraction_left + 1
      integer, parameter :: l_tau10_radius = l_tau10_mass + 1
      integer, parameter :: l_tau10_lgP = l_tau10_radius + 1
      integer, parameter :: l_tau10_lgT = l_tau10_lgP + 1
      integer, parameter :: l_tau10_lgRho = l_tau10_lgT + 1
      integer, parameter :: l_tau10_L = l_tau10_lgRho + 1
      integer, parameter :: l_tau100_mass = l_tau10_L + 1
      integer, parameter :: l_tau100_radius = l_tau100_mass + 1
      integer, parameter :: l_tau100_lgP = l_tau100_radius + 1
      integer, parameter :: l_tau100_lgT = l_tau100_lgP + 1
      integer, parameter :: l_tau100_lgRho = l_tau100_lgT + 1
      integer, parameter :: l_tau100_L = l_tau100_lgRho + 1
      integer, parameter :: l_dynamic_timescale = l_tau100_L + 1
      integer, parameter :: l_kh_timescale = l_dynamic_timescale + 1
      integer, parameter :: l_nuc_timescale = l_kh_timescale + 1 
      
      integer, parameter :: l_eps_h_max = l_nuc_timescale + 1
      integer, parameter :: l_eps_h_max_lgT = l_eps_h_max + 1
      integer, parameter :: l_eps_h_max_lgRho = l_eps_h_max_lgT + 1
      integer, parameter :: l_eps_h_max_m = l_eps_h_max_lgRho + 1
      integer, parameter :: l_eps_h_max_lgR = l_eps_h_max_m + 1
      integer, parameter :: l_eps_h_max_lgP = l_eps_h_max_lgR + 1
      integer, parameter :: l_eps_h_max_opacity = l_eps_h_max_lgP + 1
      
      integer, parameter :: l_eps_he_max = l_eps_h_max_opacity + 1
      integer, parameter :: l_eps_he_max_lgT = l_eps_he_max + 1
      integer, parameter :: l_eps_he_max_lgRho = l_eps_he_max_lgT + 1
      integer, parameter :: l_eps_he_max_m = l_eps_he_max_lgRho + 1
      integer, parameter :: l_eps_he_max_lgR = l_eps_he_max_m + 1
      integer, parameter :: l_eps_he_max_lgP = l_eps_he_max_lgR + 1
      integer, parameter :: l_eps_he_max_opacity = l_eps_he_max_lgP + 1
      
      integer, parameter :: l_eps_z_max = l_eps_he_max_opacity + 1
      integer, parameter :: l_eps_z_max_lgT = l_eps_z_max + 1
      integer, parameter :: l_eps_z_max_lgRho = l_eps_z_max_lgT + 1
      integer, parameter :: l_eps_z_max_m = l_eps_z_max_lgRho + 1
      integer, parameter :: l_eps_z_max_lgR = l_eps_z_max_m + 1
      integer, parameter :: l_eps_z_max_lgP = l_eps_z_max_lgR + 1
      integer, parameter :: l_eps_z_max_opacity = l_eps_z_max_lgP + 1
      
      
      
      integer, parameter :: l_col_id_max = l_eps_z_max_opacity
      
      integer, parameter :: maxlen_log_column_name = 32
      character (len=maxlen_log_column_name) :: log_column_name(l_col_id_max)
      type (integer_dict), pointer :: log_column_names_dict
      
      
      
   ! profile column options

      ! abundances -- names from chem_Name array in chem_def
      ! nuclear reaction eps - names from reaction_Name array in rates_def
      ! nuclear reaction category eps - names from category_name array in net_def
      ! items specified in another file - include 'filename'
      ! star_data items from the following list

      integer, parameter :: p_zone = 1
      integer, parameter :: p_luminosity = p_zone + 1
      integer, parameter :: p_logL = p_luminosity + 1
      integer, parameter :: p_velocity = p_logL + 1
      integer, parameter :: p_radius = p_velocity + 1
      integer, parameter :: p_logR = p_radius + 1
      integer, parameter :: p_q = p_logR + 1
      integer, parameter :: p_dq = p_q + 1
      
      integer, parameter :: p_pressure_beta = p_dq + 1
      integer, parameter :: p_mass = p_pressure_beta + 1
      integer, parameter :: p_mmid = p_mass + 1
      integer, parameter :: p_xm = p_mmid + 1
      integer, parameter :: p_logxq = p_xm + 1
      integer, parameter :: p_logdq = p_logxq + 1
      integer, parameter :: p_dq_ratio = p_logdq + 1
      integer, parameter :: p_tau = p_dq_ratio + 1
      integer, parameter :: p_log_opacity = p_tau + 1
      integer, parameter :: p_energy = p_log_opacity + 1
      integer, parameter :: p_logM = p_energy + 1
      integer, parameter :: p_logtau = p_logM + 1
      integer, parameter :: p_temperature = p_logtau + 1
      integer, parameter :: p_logT = p_temperature + 1
      integer, parameter :: p_rho = p_logT + 1
      integer, parameter :: p_logRho = p_rho + 1
      integer, parameter :: p_pgas = p_logRho+ 1
      integer, parameter :: p_logPgas = p_pgas + 1
      integer, parameter :: p_prad = p_logPgas + 1
      integer, parameter :: p_pressure = p_prad + 1
      integer, parameter :: p_logP = p_pressure + 1
      integer, parameter :: p_logE = p_logP + 1
      integer, parameter :: p_grada = p_logE + 1
      integer, parameter :: p_dE_dRho = p_grada + 1
      integer, parameter :: p_cv = p_dE_dRho + 1
      integer, parameter :: p_cp = p_cv + 1
      integer, parameter :: p_logS = p_cp + 1
      integer, parameter :: p_gamma1 = p_logS + 1
      integer, parameter :: p_gamma3 = p_gamma1 + 1
      integer, parameter :: p_eta = p_gamma3 + 1
      integer, parameter :: p_theta_e = p_eta + 1
      integer, parameter :: p_gam = p_theta_e + 1
      integer, parameter :: p_mu = p_gam + 1
      integer, parameter :: p_v_div_r = p_mu + 1
      integer, parameter :: p_v_div_csound = p_v_div_r + 1
      integer, parameter :: p_csound = p_v_div_csound + 1
      integer, parameter :: p_scale_height = p_csound + 1
      integer, parameter :: p_gradr_sub_grada = p_scale_height + 1
      integer, parameter :: p_entropy = p_gradr_sub_grada + 1
      integer, parameter :: p_free_e = p_entropy + 1
      integer, parameter :: p_logfree_e = p_free_e + 1
      integer, parameter :: p_chiRho = p_logfree_e + 1
      integer, parameter :: p_chiT = p_chiRho + 1
      integer, parameter :: p_abar = p_chiT + 1
      integer, parameter :: p_zbar = p_abar + 1
      integer, parameter :: p_z2bar = p_zbar + 1
      integer, parameter :: p_ye = p_z2bar + 1
      integer, parameter :: p_opacity = p_ye + 1
      integer, parameter :: p_eps_nuc = p_opacity + 1
      integer, parameter :: p_non_nuc_neu = p_eps_nuc + 1
      integer, parameter :: p_nonnucneu_plas = p_non_nuc_neu + 1
      integer, parameter :: p_nonnucneu_brem = p_nonnucneu_plas + 1
      integer, parameter :: p_nonnucneu_phot = p_nonnucneu_brem + 1
      integer, parameter :: p_nonnucneu_pair = p_nonnucneu_phot + 1
      integer, parameter :: p_extra_heat = p_nonnucneu_pair + 1
      integer, parameter :: p_logPvisc = p_extra_heat + 1
      integer, parameter :: p_eps_grav = p_logPvisc + 1
      integer, parameter :: p_mlt_mixing_length = p_eps_grav + 1
      integer, parameter :: p_log_cdc = p_mlt_mixing_length + 1
      integer, parameter :: p_log_cdc_Eulerian = p_log_cdc + 1
      integer, parameter :: p_log_conv_vel = p_log_cdc_Eulerian + 1
      integer, parameter :: p_conv_vel_div_csound = p_log_conv_vel + 1
      integer, parameter :: p_conv_mixing_type = p_conv_vel_div_csound + 1

      integer, parameter :: p_use_gradr_for_gradT = p_conv_mixing_type + 1
      integer, parameter :: p_log_mlt_cdc = p_use_gradr_for_gradT + 1
      integer, parameter :: p_log_mlt_cdc_Eulerian = p_log_mlt_cdc + 1
      integer, parameter :: p_pressure_scale_height = p_log_mlt_cdc_Eulerian + 1

      integer, parameter :: p_gradT = p_pressure_scale_height + 1
      integer, parameter :: p_gradr = p_gradT + 1
      integer, parameter :: p_dlnR_dm = p_gradr + 1
      integer, parameter :: p_dlnP_dm = p_dlnR_dm + 1
      integer, parameter :: p_dlnT_dm = p_dlnP_dm + 1
      integer, parameter :: p_dL_dm = p_dlnT_dm + 1
      integer, parameter :: p_dqdot = p_dL_dm + 1
      integer, parameter :: p_qdot = p_dqdot + 1
      integer, parameter :: p_qdot_times_dt = p_qdot + 1
      integer, parameter :: p_accel_div_grav = p_qdot_times_dt + 1
      
      integer, parameter :: p_dlnd_dt = p_accel_div_grav + 1
      integer, parameter :: p_dlnT_dt = p_dlnd_dt + 1
      integer, parameter :: p_dlnR_dt = p_dlnT_dt + 1
      integer, parameter :: p_dv_dt = p_dlnR_dt + 1
      integer, parameter :: p_cno_div_z = p_dv_dt + 1

      integer, parameter :: p_delta_r = p_cno_div_z + 1
      integer, parameter :: p_delta_v = p_delta_r + 1
      integer, parameter :: p_dt_dv_div_dr = p_delta_v + 1

      integer, parameter :: p_dlnH1_dlnP = p_dt_dv_div_dr + 1
      integer, parameter :: p_dlnHe3_dlnP = p_dlnH1_dlnP + 1
      integer, parameter :: p_dlnHe4_dlnP = p_dlnHe3_dlnP + 1
      integer, parameter :: p_dlnC12_dlnP = p_dlnHe4_dlnP + 1
      integer, parameter :: p_dlnC13_dlnP = p_dlnC12_dlnP + 1
      integer, parameter :: p_dlnN14_dlnP = p_dlnC13_dlnP + 1
      integer, parameter :: p_dlnO16_dlnP = p_dlnN14_dlnP + 1
      integer, parameter :: p_dlnNe20_dlnP = p_dlnO16_dlnP + 1
      integer, parameter :: p_dlnMg24_dlnP = p_dlnNe20_dlnP + 1
      integer, parameter :: p_dlnSi28_dlnP = p_dlnMg24_dlnP + 1

      integer, parameter :: p_dlog_pp_dlogP = p_dlnSi28_dlnP + 1
      integer, parameter :: p_dlog_cno_dlogP = p_dlog_pp_dlogP + 1
      integer, parameter :: p_dlog_3alf_dlogP = p_dlog_cno_dlogP + 1
         
      integer, parameter :: p_dlog_burn_c_dlogP = p_dlog_3alf_dlogP + 1
      integer, parameter :: p_dlog_burn_n_dlogP = p_dlog_burn_c_dlogP + 1
      integer, parameter :: p_dlog_burn_o_dlogP = p_dlog_burn_n_dlogP + 1
         
      integer, parameter :: p_dlog_burn_ne_dlogP = p_dlog_burn_o_dlogP + 1
      integer, parameter :: p_dlog_burn_na_dlogP = p_dlog_burn_ne_dlogP + 1
      integer, parameter :: p_dlog_burn_mg_dlogP = p_dlog_burn_na_dlogP + 1
         
      integer, parameter :: p_dlog_cc_dlogP = p_dlog_burn_mg_dlogP + 1
      integer, parameter :: p_dlog_co_dlogP = p_dlog_cc_dlogP + 1
      integer, parameter :: p_dlog_oo_dlogP = p_dlog_co_dlogP + 1
         
      integer, parameter :: p_dlog_burn_si_dlogP = p_dlog_oo_dlogP + 1
      integer, parameter :: p_dlog_burn_s_dlogP = p_dlog_burn_si_dlogP + 1
      integer, parameter :: p_dlog_burn_ar_dlogP = p_dlog_burn_s_dlogP + 1
      integer, parameter :: p_dlog_burn_ca_dlogP = p_dlog_burn_ar_dlogP + 1
      integer, parameter :: p_dlog_burn_ti_dlogP = p_dlog_burn_ca_dlogP + 1
      integer, parameter :: p_dlog_burn_cr_dlogP = p_dlog_burn_ti_dlogP + 1
      integer, parameter :: p_dlog_burn_fe_dlogP = p_dlog_burn_cr_dlogP + 1
         
      integer, parameter :: p_dlog_pnhe4_dlogP = p_dlog_burn_fe_dlogP + 1
      integer, parameter :: p_dlog_photo_dlogP = p_dlog_pnhe4_dlogP + 1
      integer, parameter :: p_dlog_other_dlogP = p_dlog_photo_dlogP + 1

      integer, parameter :: p_binding_energy = p_dlog_other_dlogP + 1

      integer, parameter :: p_brunt_N2 = p_binding_energy + 1
      integer, parameter :: p_brunt_Astar = p_brunt_N2 + 1
      integer, parameter :: p_brunt_B = p_brunt_Astar + 1
      integer, parameter :: p_chiY = p_brunt_B + 1
      integer, parameter :: p_dlnY_dlnP = p_chiY + 1
      integer, parameter :: p_logQ = p_dlnY_dlnP + 1
      integer, parameter :: p_cs_at_cell_bdy = p_logQ + 1

      integer, parameter :: p_col_id_max = p_cs_at_cell_bdy
      
      integer, parameter :: maxlen_profile_column_name = 32
      character (len=maxlen_profile_column_name) :: profile_column_name(p_col_id_max)
      type (integer_dict), pointer :: profile_column_names_dict


      contains
      
      
      subroutine result_reason_init         
         result_reason_str(result_reason_normal) = 'normal'
         result_reason_str(dt_is_zero) = 'dt_is_zero'
         result_reason_str(nonzero_ierr) = 'nonzero_ierr'
         result_reason_str(hydro_failed_to_converge) = 'hydro_failed_to_converge'
         result_reason_str(hydro_err_too_large) = 'hydro_err_too_large'
         result_reason_str(burn_and_mix_failed) = 'burn_and_mix_failed'
         result_reason_str(diffusion_failed) = 'element_diffusion_failed'
         result_reason_str(too_many_steps_for_burn) = 'too_many_steps_for_burn'
         result_reason_str(too_many_steps_for_diffusion) = 'too_many_steps_for_diffusion'
         result_reason_str(too_many_steps_for_hydro) = 'too_many_steps_for_hydro'
         result_reason_str(adjust_mesh_failed) = 'adjust_mesh_failed'
         result_reason_str(adjust_mass_failed) = 'adjust_mass_failed'
         result_reason_str(core_dump_model_number) = 'core_dump_model_number'
         result_reason_str(timestep_limits) = 'timestep_limits'
         result_reason_str(variable_change_limits) = 'variable_change_limits'
         result_reason_str(bad_lnd_prediction) = 'bad_lnd_prediction'
         result_reason_str(bad_lnT_prediction) = 'bad_lnT_prediction'
         result_reason_str(op_split_failed_to_converge) = 'op_split_failed_to_converge'
      end subroutine result_reason_init
      
      
      subroutine log_column_names_init
         use utils_lib, only: integer_dict_define
         integer :: i, cnt, ierr
         
         cnt = 0
         log_column_name(:) = ''

         log_column_name(l_model_number) = 'model_number'
         log_column_name(l_star_age) = 'star_age'
         log_column_name(l_star_mass) = 'star_mass'
         log_column_name(l_star_mdot) = 'star_mdot'
         log_column_name(l_log_abs_mdot) = 'log_abs_mdot'
         log_column_name(l_time_step) = 'time_step'
         log_column_name(l_num_zones) = 'num_zones'
         log_column_name(l_conv_mx1_top) = 'conv_mx1_top'
         log_column_name(l_conv_mx1_bot) = 'conv_mx1_bot'
         log_column_name(l_conv_mx2_top) = 'conv_mx2_top'
         log_column_name(l_conv_mx2_bot) = 'conv_mx2_bot'
         log_column_name(l_mx1_top) = 'mx1_top'
         log_column_name(l_mx1_bot) = 'mx1_bot'
         log_column_name(l_mx2_top) = 'mx2_top'
         log_column_name(l_mx2_bot) = 'mx2_bot'
         log_column_name(l_mixing_regions) = 'mixing_regions'
         log_column_name(l_epsnuc_M_1) = 'epsnuc_M_1'
         log_column_name(l_epsnuc_M_2) = 'epsnuc_M_2'
         log_column_name(l_epsnuc_M_3) = 'epsnuc_M_3'
         log_column_name(l_epsnuc_M_4) = 'epsnuc_M_4'
         log_column_name(l_epsnuc_M_5) = 'epsnuc_M_5'
         log_column_name(l_epsnuc_M_6) = 'epsnuc_M_6'
         log_column_name(l_epsnuc_M_7) = 'epsnuc_M_7'
         log_column_name(l_epsnuc_M_8) = 'epsnuc_M_8'
         log_column_name(l_h1_boundary_mass) = 'h1_boundary_mass'
         log_column_name(l_he4_boundary_mass) = 'he4_boundary_mass'
         log_column_name(l_power_h_burn) = 'power_h_burn'
         log_column_name(l_power_he_burn) = 'power_he_burn'
         log_column_name(l_log_center_T) = 'log_center_T'
         log_column_name(l_log_center_Rho) = 'log_center_Rho'
         log_column_name(l_v_div_csound_surf) = 'v_div_csound_surf'
         log_column_name(l_surface_accel_div_grav) = 'surface_accel_div_grav'
         log_column_name(l_log_dt) = 'log_dt'
         log_column_name(l_log_LH) = 'log_LH'
         log_column_name(l_log_LHe) = 'log_LHe'
         log_column_name(l_log_L) = 'log_L'
         log_column_name(l_log_R) = 'log_R'
         log_column_name(l_log_Teff) = 'log_Teff'
         log_column_name(l_log_g) = 'log_g'
         log_column_name(l_log_L_div_Ledd) = 'log_L_div_Ledd'
         log_column_name(l_num_retries) = 'num_retries'
         log_column_name(l_num_backups) = 'num_backups'
         log_column_name(l_h1_czb_mass) = 'h1_czb_mass'
         log_column_name(l_surf_c12_minus_o16) = 'surf_c12_minus_o16'

         log_column_name(l_log_center_P) = 'log_center_P' ! log10(center pressure in dynes/cm^2)
         log_column_name(l_center_degeneracy) = 'center_degeneracy' ! the electron chemical potential in units of k*T
         log_column_name(l_center_gamma) = 'center_gamma' ! plasma interaction parameter
         log_column_name(l_h1_boundary_radius) = 'h1_boundary_radius'
         log_column_name(l_h1_boundary_lgT) = 'h1_boundary_lgT'
         log_column_name(l_h1_boundary_lgRho) = 'h1_boundary_lgRho'
         log_column_name(l_h1_boundary_L) = 'h1_boundary_L'
         log_column_name(l_h1_boundary_v) = 'h1_boundary_v'
         log_column_name(l_he4_boundary_radius) = 'he4_boundary_radius'
         log_column_name(l_he4_boundary_lgT) = 'he4_boundary_lgT'
         log_column_name(l_he4_boundary_lgRho) = 'he4_boundary_lgRho'
         log_column_name(l_he4_boundary_L) = 'he4_boundary_L'
         log_column_name(l_he4_boundary_v) = 'he4_boundary_v'
         log_column_name(l_c12_boundary_mass) = 'c12_boundary_mass'
         log_column_name(l_c12_boundary_radius) = 'c12_boundary_radius'
         log_column_name(l_c12_boundary_lgT) = 'c12_boundary_lgT'
         log_column_name(l_c12_boundary_lgRho) = 'c12_boundary_lgRho'
         log_column_name(l_c12_boundary_L) = 'c12_boundary_L'
         log_column_name(l_c12_boundary_v) = 'c12_boundary_v'
         log_column_name(l_envelope_mass) = 'envelope_mass'
         log_column_name(l_envelope_fraction_left) = 'envelope_fraction_left'
         
         log_column_name(l_tau10_mass) = 'tau10_mass' ! mass in solar units where optical depth = 10
         log_column_name(l_tau10_radius) = 'tau10_radius' ! radius in solar units where optical depth = 10
         log_column_name(l_tau10_lgP) = 'tau10_lgP' ! estimate for log10(P) at tau = 10
         log_column_name(l_tau10_lgT) = 'tau10_lgT' ! estimate for log10(T) at tau = 10
         log_column_name(l_tau10_lgRho) = 'tau10_lgRho' ! estimate for log10(density) at tau = 10
         log_column_name(l_tau10_L) = 'tau10_L' ! estimate for L/Lsun at tau = 10
         log_column_name(l_tau100_mass) = 'tau100_mass' ! location in solar units where optical depth = 100
         log_column_name(l_tau100_radius) = 'tau100_radius' ! location in solar units where optical depth = 100
         log_column_name(l_tau100_lgP) = 'tau100_lgP' ! estimates for values at tau = 100
         log_column_name(l_tau100_lgT) = 'tau100_lgT'
         log_column_name(l_tau100_lgRho) = 'tau100_lgRho'
         log_column_name(l_tau100_L) = 'tau100_L'
         log_column_name(l_dynamic_timescale) = 'dynamic_timescale'
         log_column_name(l_kh_timescale) = 'kh_timescale'
         log_column_name(l_nuc_timescale) = 'nuc_timescale' 
         
         log_column_name(l_eps_h_max) = 'eps_h_max' ! erg/g/s
         log_column_name(l_eps_h_max_lgT) = 'eps_h_max_lgT' ! log10 temperature at location of max burn
         log_column_name(l_eps_h_max_lgRho) = 'eps_h_max_lgRho' ! log10 density at location of max burn
         log_column_name(l_eps_h_max_m) = 'eps_h_max_m' ! mass coordinate at location of max burn (Msun units)
         log_column_name(l_eps_h_max_lgR) = 'eps_h_max_lgR'
         log_column_name(l_eps_h_max_lgP) = 'eps_h_max_lgP'
         log_column_name(l_eps_h_max_opacity) = 'eps_h_max_opacity'
         
         log_column_name(l_eps_he_max) = 'eps_he_max' ! erg/g/s
         log_column_name(l_eps_he_max_lgT) = 'eps_he_max_lgT' ! log10 temperature at location of max burn
         log_column_name(l_eps_he_max_lgRho) = 'eps_he_max_lgRho' ! log10 density at location of max burn
         log_column_name(l_eps_he_max_m) = 'eps_he_max_m' ! mass coordinate at location of max burn (Msun units)
         log_column_name(l_eps_he_max_lgR) = 'eps_he_max_lgR'
         log_column_name(l_eps_he_max_lgP) = 'eps_he_max_lgP'
         log_column_name(l_eps_he_max_opacity) = 'eps_he_max_opacity'
         
         log_column_name(l_eps_z_max) = 'eps_z_max' ! erg/g/s
         log_column_name(l_eps_z_max_lgT) = 'eps_z_max_lgT' ! log10 temperature at location of max burn
         log_column_name(l_eps_z_max_lgRho) = 'eps_z_max_lgRho' ! log10 density at location of max burn
         log_column_name(l_eps_z_max_m) = 'eps_z_max_m' ! mass coordinate at location of max burn (Msun units)      
         log_column_name(l_eps_z_max_lgR) = 'eps_z_max_lgR'
         log_column_name(l_eps_z_max_lgP) = 'eps_z_max_lgP'
         log_column_name(l_eps_z_max_opacity) = 'eps_z_max_opacity'
                  
         cnt = 0
         do i=1,l_col_id_max
            if (len_trim(log_column_name(i)) == 0) then
               write(*,*) 'missing name for log column id', i
               if (i > 1) write(*,*) 'following ' // trim(log_column_name(i-1))
               write(*,*) 
               cnt = cnt+1
            end if
         end do

         if (cnt > 0) stop 'log_column_names_init'
         
         nullify(log_column_names_dict)
         do i=1,l_col_id_max
            call integer_dict_define(log_column_names_dict, log_column_name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: log_column_names_init failed in integer_dict_define'
               stop 1
            end if
         end do

      end subroutine log_column_names_init         
            
      
		integer function do_get_log_id(cname)
         use net_def
         use utils_lib
			character (len=*), intent(in)  :: cname 
			! returns id for the log column if there is a matching name
			! returns 0 otherwise.
			integer :: ierr, value
			call integer_dict_lookup(log_column_names_dict, cname, value, ierr)
			if (ierr /= 0) value = 0
			do_get_log_id = value
		end function do_get_log_id
         
      
      
      
      subroutine profile_column_names_init
         use utils_lib, only: integer_dict_define
         integer :: i, cnt, ierr
         
         cnt = 0
         profile_column_name(:) = ''

         profile_column_name(p_zone) = 'zone'
         profile_column_name(p_luminosity) = 'luminosity'
         profile_column_name(p_logL) = 'logL'
         profile_column_name(p_velocity) = 'velocity'
         profile_column_name(p_radius) = 'radius'
         profile_column_name(p_logR) = 'logR'
         profile_column_name(p_q) = 'q'
         profile_column_name(p_dq) = 'dq'
         profile_column_name(p_logtau) = 'logtau'
         profile_column_name(p_pressure_beta) = 'pressure_beta'
         profile_column_name(p_mass) = 'mass'
         profile_column_name(p_mmid) = 'mmid'
         profile_column_name(p_xm) = 'xm'
         profile_column_name(p_logxq) = 'logxq'
         profile_column_name(p_logdq) = 'logdq'
         profile_column_name(p_dq_ratio) = 'dq_ratio'
         profile_column_name(p_tau) = 'tau'
         profile_column_name(p_log_opacity) = 'log_opacity'
         profile_column_name(p_energy) = 'energy'
         profile_column_name(p_logM) = 'logM'
         profile_column_name(p_temperature) = 'temperature'
         profile_column_name(p_logT) = 'logT'
         profile_column_name(p_rho) = 'rho'
         profile_column_name(p_logRho) = 'logRho'
         profile_column_name(p_pgas) = 'pgas'
         profile_column_name(p_logPgas) = 'logPgas'
         profile_column_name(p_prad) = 'prad'
         profile_column_name(p_pressure) = 'pressure'
         profile_column_name(p_logP) = 'logP'
         profile_column_name(p_logE) = 'logE'
         profile_column_name(p_grada) = 'grada'
         profile_column_name(p_dE_dRho) = 'dE_dRho'
         profile_column_name(p_cv) = 'cv'
         profile_column_name(p_cp) = 'cp'
         profile_column_name(p_logS) = 'logS'
         profile_column_name(p_gamma1) = 'gamma1'
         profile_column_name(p_gamma3) = 'gamma3'
         profile_column_name(p_eta) = 'eta'
         profile_column_name(p_theta_e) = 'theta_e'
         profile_column_name(p_gam) = 'gam'
         profile_column_name(p_mu) = 'mu'
         profile_column_name(p_v_div_r) = 'v_div_r'
         profile_column_name(p_v_div_csound) = 'v_div_csound'
         profile_column_name(p_csound) = 'csound'
         profile_column_name(p_gradr_sub_grada) = 'gradr_sub_grada'
         profile_column_name(p_scale_height) = 'scale_height'
         profile_column_name(p_entropy) = 'entropy'
         profile_column_name(p_free_e) = 'free_e'
         profile_column_name(p_logfree_e) = 'logfree_e'
         profile_column_name(p_chiRho) = 'chiRho'
         profile_column_name(p_chiT) = 'chiT'
         profile_column_name(p_abar) = 'abar'
         profile_column_name(p_zbar) = 'zbar'
         profile_column_name(p_z2bar) = 'z2bar'
         profile_column_name(p_ye) = 'ye'
         profile_column_name(p_opacity) = 'opacity'
         profile_column_name(p_eps_nuc) = 'eps_nuc'
         profile_column_name(p_non_nuc_neu) = 'non_nuc_neu'
         profile_column_name(p_nonnucneu_plas) = 'nonnucneu_plas'
         profile_column_name(p_nonnucneu_brem) = 'nonnucneu_brem'
         profile_column_name(p_nonnucneu_phot) = 'nonnucneu_phot'
         profile_column_name(p_nonnucneu_pair) = 'nonnucneu_pair'
         profile_column_name(p_extra_heat) = 'extra_heat'
         profile_column_name(p_logPvisc) = 'logPvisc'
         profile_column_name(p_eps_grav) = 'eps_grav'
         profile_column_name(p_mlt_mixing_length) = 'mlt_mixing_length'
         profile_column_name(p_log_cdc) = 'log_cdc'
         profile_column_name(p_log_cdc_Eulerian) = 'log_cdc_Eulerian'
         profile_column_name(p_log_conv_vel) = 'log_conv_vel'
         profile_column_name(p_conv_vel_div_csound) = 'conv_vel_div_csound'
         profile_column_name(p_conv_mixing_type) = 'conv_mixing_type'
         profile_column_name(p_use_gradr_for_gradT) = 'use_gradr_for_gradT'
         profile_column_name(p_log_mlt_cdc) = 'log_mlt_cdc'
         profile_column_name(p_log_mlt_cdc_Eulerian) = 'log_mlt_cdc_Eulerian'
         profile_column_name(p_pressure_scale_height) = 'pressure_scale_height'
         profile_column_name(p_gradT) = 'gradT'
         profile_column_name(p_gradr) = 'gradr'
         profile_column_name(p_dlnR_dm) = 'dlnR_dm'
         profile_column_name(p_dlnP_dm) = 'dlnP_dm'
         profile_column_name(p_dlnT_dm) = 'dlnT_dm'
         profile_column_name(p_dL_dm) = 'dL_dm'
         profile_column_name(p_dqdot) = 'dqdot'
         profile_column_name(p_qdot) = 'qdot'
         profile_column_name(p_qdot_times_dt) = 'qdot_times_dt'
         profile_column_name(p_accel_div_grav) = 'accel_div_grav'
         profile_column_name(p_dlnd_dt) = 'dlnd_dt'
         profile_column_name(p_dlnT_dt) = 'dlnT_dt'
         profile_column_name(p_dlnR_dt) = 'dlnR_dt'
         profile_column_name(p_dv_dt) = 'dv_dt'
         profile_column_name(p_cno_div_z) = 'cno_div_z'

         profile_column_name(p_delta_r) = 'delta_r'
         profile_column_name(p_delta_v) = 'delta_v'
         profile_column_name(p_dt_dv_div_dr) = 'dt_dv_div_dr'
         
         profile_column_name(p_dlnH1_dlnP) = 'dlog_h1_dlogP'
         profile_column_name(p_dlnHe3_dlnP) = 'dlog_he3_dlogP'
         profile_column_name(p_dlnHe4_dlnP) = 'dlog_he4_dlogP'
         profile_column_name(p_dlnC12_dlnP) = 'dlog_c12_dlogP'
         profile_column_name(p_dlnC13_dlnP) = 'dlog_c13_dlogP'
         profile_column_name(p_dlnN14_dlnP) = 'dlog_n14_dlogP'
         profile_column_name(p_dlnO16_dlnP) = 'dlog_o16_dlogP'
         profile_column_name(p_dlnNe20_dlnP) = 'dlog_ne20_dlogP'
         profile_column_name(p_dlnMg24_dlnP) = 'dlog_mg24_dlogP'
         profile_column_name(p_dlnSi28_dlnP) = 'dlog_si28_dlogP'

         profile_column_name(p_dlog_pp_dlogP) = 'dlog_pp_dlogP'
         profile_column_name(p_dlog_cno_dlogP) = 'dlog_cno_dlogP'
         profile_column_name(p_dlog_3alf_dlogP) = 'dlog_3alf_dlogP'
         
         profile_column_name(p_dlog_burn_c_dlogP) = 'dlog_burn_c_dlogP'
         profile_column_name(p_dlog_burn_n_dlogP) = 'dlog_burn_n_dlogP'
         profile_column_name(p_dlog_burn_o_dlogP) = 'dlog_burn_o_dlogP'
         
         profile_column_name(p_dlog_burn_ne_dlogP) = 'dlog_burn_ne_dlogP'
         profile_column_name(p_dlog_burn_na_dlogP) = 'dlog_burn_na_dlogP'
         profile_column_name(p_dlog_burn_mg_dlogP) = 'dlog_burn_mg_dlogP'
         
         profile_column_name(p_dlog_cc_dlogP) = 'dlog_cc_dlogP'
         profile_column_name(p_dlog_co_dlogP) = 'dlog_co_dlogP'
         profile_column_name(p_dlog_oo_dlogP) = 'dlog_oo_dlogP'
         
         profile_column_name(p_dlog_burn_si_dlogP) = 'dlog_burn_si_dlogP'
         profile_column_name(p_dlog_burn_s_dlogP) = 'dlog_burn_s_dlogP'
         profile_column_name(p_dlog_burn_ar_dlogP) = 'dlog_burn_ar_dlogP'
         profile_column_name(p_dlog_burn_ca_dlogP) = 'dlog_burn_ca_dlogP'
         profile_column_name(p_dlog_burn_ti_dlogP) = 'dlog_burn_ti_dlogP'
         profile_column_name(p_dlog_burn_cr_dlogP) = 'dlog_burn_cr_dlogP'
         profile_column_name(p_dlog_burn_fe_dlogP) = 'dlog_burn_fe_dlogP'
         
         profile_column_name(p_dlog_pnhe4_dlogP) = 'dlog_pnhe4_dlogP'
         profile_column_name(p_dlog_photo_dlogP) = 'dlog_photo_dlogP'
         profile_column_name(p_dlog_other_dlogP) = 'dlog_other_dlogP'

         profile_column_name(p_binding_energy) = 'binding_energy'

         profile_column_name(p_brunt_N2) = 'brunt_N2'
         profile_column_name(p_brunt_Astar) = 'brunt_Astar'
         profile_column_name(p_brunt_B) = 'brunt_B'
         profile_column_name(p_chiY) = 'chiY'
         profile_column_name(p_dlnY_dlnP) = 'dlnY_dlnP'
         profile_column_name(p_cs_at_cell_bdy) = 'cs_at_cell_bdy'
         profile_column_name(p_logQ) = 'logQ'
         
         cnt = 0
         do i=1,p_col_id_max
            if (len_trim(profile_column_name(i)) == 0) then
               write(*,*) 'missing name for profile column id', i
               if (i > 1) write(*,*) 'following ' // trim(profile_column_name(i-1))
               write(*,*) 
               cnt = cnt+1
            end if
         end do

         if (cnt > 0) stop 'profile_column_names_init'
         
         nullify(profile_column_names_dict)
         do i=1,p_col_id_max
            call integer_dict_define(profile_column_names_dict, profile_column_name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: profile_column_names_init failed in integer_dict_define'
               stop 1
            end if
         end do

      end subroutine profile_column_names_init         
            
      
		integer function do_get_profile_id(cname)
         use net_def
         use utils_lib
			character (len=*), intent(in)  :: cname 
			! returns id for the profile column if there is a matching name
			! returns 0 otherwise.
			integer :: ierr, value
			call integer_dict_lookup(profile_column_names_dict, cname, value, ierr)
			if (ierr /= 0) value = 0
			do_get_profile_id = value
		end function do_get_profile_id
      
      
      end module star_def

