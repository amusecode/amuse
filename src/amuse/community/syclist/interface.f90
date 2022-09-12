function initialize_code()
  use DataStructure, only: verbose
  use VariousParameters, only: Std_Path,Std_Path_tables
  use VariousParameters, only: grid, table_format, PMS, star_number, i_metallicity, &
          fixed_metallicity, IMF_type, m_IMF_inf, m_IMF_sup, ivdist, om_ivdist, &
          iangle, Fixed_AoV_latitude, inoise, sigma_mv, sigma_bv, binary_prob, &
          Colour_Calibration_mode, grav_dark, limb_dark
  use Population_Mode, only: Pop_Mass_Beam_Number, Pop_Omega_Beam_Number, N_Time_step
  implicit none
  integer :: initialize_code

  ! Set default values.
  ! Not all of these are used in AMUSE but we set them just in case for now.
  verbose = .true.  ! Control output with 'redirection' option in AMUSE
  ! Root path for Syclist data. Should move to AMUSE data dir at some point.
  Std_Path = './src/SYCLIST'  
  Std_Path_tables = trim(Std_Path)//'/tables/'
  Std_Path = trim(Std_Path)//'/inputs/'

  !grid = 'BeGrids'
  grid = 'Grids2012'
  table_format = 1
  PMS = .false.
  star_number = 1
  i_metallicity = 0
  fixed_metallicity = 0.014d0
  IMF_type = 1
  m_IMF_inf = 0.8d0
  m_IMF_sup = 120.d0
  ivdist = 2
  om_ivdist = 0.4d0
  iangle = 0
  Fixed_AoV_latitude = 22.d0
  inoise = 0
  sigma_mv = 0.01d0
  sigma_bv = 0.0025d0
  binary_prob = 0.d0
  Colour_Calibration_mode = 2
  grav_dark = 2
  limb_dark = 0

  Pop_Mass_Beam_Number = 2000
  Pop_Omega_Beam_Number = 100
  N_Time_step = 500

  initialize_code = 0
end function

function set_models_path(path)
  use VariousParameters, only: Std_Path
  implicit none
  character(256) :: path
  integer :: set_models_path
  Std_Path = path
  set_models_path = 0
end function

function get_models_path(path)
  use VariousParameters, only: Std_Path
  implicit none
  character(256) :: path
  integer :: get_models_path
  path = Std_Path
  get_models_path = 0
end function

function evolve_star(&
                age, initial_mass, metallicity, omega,&
                mass, luminosity, temperature,&
                ab_surf_H1, ab_surf_He4, ab_surf_C12, ab_surf_C13, ab_surf_N14,&
                ab_surf_O16, ab_surf_O17, ab_surf_O18, ab_surf_Ne20, ab_surf_Ne22,&
                ab_surf_Al26, MccMt, temperature_wr, Md, rhoc, ab_core_H1, ab_core_He4,&
                ab_core_C12, ab_core_C13, ab_core_N14, ab_core_O16, ab_core_O17, ab_core_Ne20,&
                ab_core_Ne22, ab_core_Al26, omega_surf, omega_core, RpReq, MdMd0,&
                v_crit1, v_crit2, v_eq, OmOmcr, Gamma_Ed, MdotMech, Ltot,&
                !StarType, &
                n&
                )
  use DataStructure, only: Table_Line_Number
  use DataStructure, only: &
      i_mass,i_logL,i_logTeff_corr,i_H1_Surf,i_He4_surf,i_C12_surf,i_C13_surf,i_N14_surf, &
      i_O16_surf,i_O17_surf,i_O18_surf,i_Ne20_surf,i_Ne22_surf,i_Al26_surf,i_logTeff,i_Mdot, &
      i_Omega_surf,i_oblat,i_v_crit1,i_v_crit2,i_v_equa,i_Omega_Omcrit,i_Gamma_Ed,i_Mdot_mec, &
      i_PolarRadius,i_polar_gravity,Table_Line_Number,Data_Number,i_MV_noisy,i_BV_noisy,&
      i_MBol,i_MV,i_UB,i_BV,i_B2V1,i_VK,i_VR,i_VI,i_JK,i_HK,i_BC,i_logL_gd,i_logTeff_gd,i_logL_lgd, &
      i_logTeff_lgd,i_mean_gravity,i_GV,i_GbpV,i_GrpV,i_Gflag,i_crossing_nb_F,i_crossing_nb_1O,i_P_F, &
      i_Pdot_P_F,i_omi_omr_F,i_P_1O,i_Pdot_P_1O,i_omi_omr_1O, &
      i_al26_cen, i_c12_cen, i_c13_cen, i_h1_cen, i_he4_cen, i_l_tot, i_mcc, i_mdot_enhencement, &
      i_n14_cen, i_ne20_cen, i_ne22_cen, i_o16_cen, i_o17_cen, i_omega_cen, i_rhoc
  use InterpolationLoop, only:Initialise_Position_and_factor
  use interpolmod, only: All_Positions_and_factors,Make_InterpolatedModel,Make_TimeModel
  use Additional_Data, only: Compute_Additional,Compute_Additional_Single
  use VariousParameters, only:&
          Comp_Mode,age_log,init_AoV,&
          Star_Z,Star_mass,Star_omega,Star_AoV,Fixed_AoV,Current_Number,&
          Z_Number, Z_List
  use VariousParameters, only: inoise,iangle,star_number
  use LoopVariables, only:&
          Z_Position,Z_factor,omega_Position,omega_factor,&
          mass_Position,mass_factor,Interpolated_Model,CurrentTime_Model
  implicit none
  integer :: n, i
  double precision :: age(n), initial_mass(n), metallicity(n), omega(n) !in
  double precision :: mass(n), luminosity(n), temperature(n),&
          ab_surf_H1(n), ab_surf_He4(n), ab_surf_C12(n), ab_surf_C13(n), ab_surf_N14(n),&
          ab_surf_O16(n), ab_surf_O17(n), ab_surf_O18(n), ab_surf_Ne20(n), ab_surf_Ne22(n),&
          ab_surf_Al26(n), MccMt(n), temperature_wr(n), Md(n), rhoc(n), ab_core_H1(n), ab_core_He4(n),&
          ab_core_C12(n), ab_core_C13(n), ab_core_N14(n), ab_core_O16(n), ab_core_O17(n), ab_core_Ne20(n),&
          ab_core_Ne22(n), ab_core_Al26(n), omega_surf(n), omega_core(n), RpReq(n), MdMd0(n),&
          v_crit1(n), v_crit2(n), v_eq(n), OmOmcr(n), Gamma_Ed(n), MdotMech(n), Ltot(n)
  integer :: evolve_star
  !integer :: StarType(n)

  ! StarTypes:
  ! positives:
  ! stars that can be calculated with Syclist
  ! negatives:
  ! stars that fall outside the scope
  ! Numbers are from AMUSE
  !  "deeply or fully convective low mass MS star",  # 0
  !  "Main Sequence star",  # 1
  !  "Hertzsprung Gap",  # 2
  !  "First Giant Branch",  # 3
  !  "Core Helium Burning",  # 4
  !  "First Asymptotic Giant Branch",  # 5
  !  "Second Asymptotic Giant Branch",  # 6
  !  "Main Sequence Naked Helium star",  # 7
  !  "Hertzsprung Gap Naked Helium star",  # 8
  !  "Giant Branch Naked Helium star",  # 9
  !  "Helium White Dwarf",  # 10
  !  "Carbon/Oxygen White Dwarf",  # 11
  !  "Oxygen/Neon White Dwarf",  # 12
  !  "Neutron Star",  # 13
  !  "Black Hole",  # 14
  !  "Massless Supernova",  # 15
  !  "Unknown stellar type",  # 16
  !  "Pre-main-sequence Star",  # 17
  !  "Planet",  # 18
  !  "Brown Dwarf",  # 19

  !Comp_Mode = 1  ! stellar cluster
  !Comp_Mode = 2  ! isochrones
  !Comp_Mode = 3  ! population as a function of time
  Comp_Mode = 4  ! single star
  
  i = 1  
  Current_Number = 1

  !write(*,*) metallicity(i), initial_mass(i), omega(i), log10(age(i)), n

  do while (i <= n)
  !StarType = 1  ! Main sequence star as default

  ! Make sure metallicity is in the range of models, and adjust if needed.
  Star_Z = metallicity(i)
  if (Star_Z < Z_List(1)) then
          Star_Z = Z_List(1)
  else if (Star_Z > Z_List(Z_Number)) then
          Star_Z = Z_List(Z_Number)
  end if
  metallicity(i) = Star_Z

  Star_mass = initial_mass(i)
  Star_omega = omega(i)
  Star_AoV = Fixed_AoV
  age_log = log10(age(i))  ! age in (julian)year?

  !write(*,*) Star_mass, Star_Z, Star_omega, age_log

  !call Amuse_SingleModelMode
  inoise = 0
  if (iangle /= 3) then
          iangle = 0
  endif

  ! Number of stars to compute
  star_number = Table_Line_Number
  
  call init_AoV
  ! if (iangle == 4) then
  !   call init_angle_external
  ! endif

  call Initialise_Position_and_factor(&
          Z_Position,Z_factor,mass_Position,mass_factor,omega_Position, &
          omega_factor)

  call All_Positions_and_factors(&
          Star_Z,Z_Position,Z_factor,Star_mass,mass_Position(:),mass_factor(:), &
          Star_omega,omega_Position(:,:),omega_factor(:,:))

  ! Perform the interpolation
  call Make_InterpolatedModel(&
          Z_Position,Z_factor,mass_Position,mass_factor,omega_Position, &
          omega_factor,Interpolated_Model)

  call Make_TimeModel(Interpolated_Model,age_log,CurrentTime_Model(Current_Number))
  call Compute_Additional(CurrentTime_Model(Current_Number))
  
  !! Computation of the additional quantities.
  !Call Compute_Additional_Single(&
  !        Interpolated_Model,Star_AoV,CurrentTime_Model,Table_Line_Number)

  !call WriteResults(table_format)
  mass(i) = CurrentTime_Model(Current_Number)%Data_Line(i_mass)
  luminosity(i) = 10**(CurrentTime_Model(Current_Number)%Data_Line(i_logL))
  temperature(i) = 10**(CurrentTime_Model(Current_Number)%Data_Line(i_logTeff))
  ab_surf_H1(i) = CurrentTime_Model(Current_Number)%Data_Line(i_H1_Surf)
  ab_surf_He4(i) = CurrentTime_Model(Current_Number)%Data_Line(i_He4_Surf)
  ab_surf_C12(i) = CurrentTime_Model(Current_Number)%Data_Line(i_C12_Surf)
  ab_surf_C13(i) = CurrentTime_Model(Current_Number)%Data_Line(i_C13_Surf)
  ab_surf_N14(i) = CurrentTime_Model(Current_Number)%Data_Line(i_N14_Surf)
  ab_surf_O16(i) = CurrentTime_Model(Current_Number)%Data_Line(i_O16_Surf)
  ab_surf_O17(i) = CurrentTime_Model(Current_Number)%Data_Line(i_O17_Surf)
  ab_surf_O18(i) = CurrentTime_Model(Current_Number)%Data_Line(i_O18_Surf)
  ab_surf_Ne20(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Ne20_Surf)
  ab_surf_Ne22(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Ne22_Surf)
  ab_surf_Al26(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Al26_Surf)
  MccMt(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Mcc)
  temperature_wr(i) = CurrentTime_Model(Current_Number)%Data_Line(i_logTeff_corr)
  Md(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Mdot)
  rhoc(i) = CurrentTime_Model(Current_Number)%Data_Line(i_rhoc)
  ab_core_H1(i) = CurrentTime_Model(Current_Number)%Data_Line(i_H1_cen)
  ab_core_He4(i) = CurrentTime_Model(Current_Number)%Data_Line(i_He4_cen)
  ab_core_C12(i) = CurrentTime_Model(Current_Number)%Data_Line(i_C12_cen)
  ab_core_C13(i) = CurrentTime_Model(Current_Number)%Data_Line(i_C13_cen)
  ab_core_N14(i) = CurrentTime_Model(Current_Number)%Data_Line(i_N14_cen)
  ab_core_O16(i) = CurrentTime_Model(Current_Number)%Data_Line(i_O16_cen)
  ab_core_O17(i) = CurrentTime_Model(Current_Number)%Data_Line(i_O17_cen)
  ab_core_Ne20(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Ne20_cen)
  ab_core_Ne22(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Ne22_cen)
  ab_core_Al26(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Al26_cen)
  omega_surf(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Omega_surf)
  omega_core(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Omega_cen)
  RpReq(i) = CurrentTime_Model(Current_Number)%Data_Line(i_oblat) ! ?
  MdMd0(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Mdot_enhencement) ! ?
  v_crit1(i) = CurrentTime_Model(Current_Number)%Data_Line(i_v_crit1)
  v_crit2(i) = CurrentTime_Model(Current_Number)%Data_Line(i_v_crit2)
  v_eq(i) = CurrentTime_Model(Current_Number)%Data_Line(i_v_equa)
  OmOmcr(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Omega_Omcrit)
  Gamma_Ed(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Gamma_Ed)
  MdotMech(i) = CurrentTime_Model(Current_Number)%Data_Line(i_Mdot_mec)
  Ltot(i) = CurrentTime_Model(Current_Number)%Data_Line(i_L_tot)
  i = i + 1
  end do
  evolve_star=0
end function

function evolve_for(index_of_the_star, delta_t)
  implicit none
  integer :: index_of_the_star
  double precision :: delta_t
  integer :: evolve_for
  evolve_for=0
end function

function recommit_particles()
  implicit none
  integer :: recommit_particles
  recommit_particles=0
end function

function cleanup_code()
  implicit none
  integer :: cleanup_code
  cleanup_code=0
end function

function recommit_parameters()
  implicit none
  integer :: recommit_parameters
  recommit_parameters=0
end function

function evolve_one_step(index_of_the_star)
  implicit none
  integer :: index_of_the_star
  integer :: evolve_one_step
  evolve_one_step=0
end function

function commit_parameters()
  use ReadData, only: &
          init_Huang,init_HG,init_Correction,init_VcritOmega, init_SurfaceOmega, init_Correct_fact
  use InOut, only: InitialiseData
  use InterpolationLoop, only: Initialise
  use VariousParameters, only: Comp_Mode

  implicit none
  integer :: commit_parameters

  call init_Huang
  call init_HG
  call init_Correction
  call init_VcritOmega
  call init_SurfaceOmega
  call init_Correct_fact
  Comp_Mode = 4

  ! Various initialisations
  call Initialise

  ! Searching, opening, sorting and reading evolution files. These files must be listed (with the whole
  ! access path) in a file called "ListFile.txt".
  call InitialiseData

  commit_parameters=0
end function


