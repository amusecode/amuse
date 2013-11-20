subroutine rad_init
  use cooling_module
  use data_common
  use radiation_commons
  implicit none
  integer::aton_init_gpu

  call aton_get_grid_size( &
       & grid_size_x, &
       & grid_size_y, &
       & grid_size_z, &
       & boundary_size)
  num_cpu_x = 2**levelmin / grid_size_x
  num_cpu_y = 2**levelmin / grid_size_y
  num_cpu_z = 2**levelmin / grid_size_z

  if (rad_aton_version.eq.'gpu') then
    if (aton_init_gpu(allow_gpu_overload).eq.0) then
      write(*,*) "aton_init_gpu error"
      call clean_stop
    end if
    call aton_gpu_malloc(rad_num_sources)
  end if

  call rad_cpu_malloc()

  cpu_e=1e-7
  cpu_f=0.0
  cpu_x=0.0
  cpu_t=0.0
  cpu_d=0.0

  ! These are initialized each time step because they depends on dx.
  cpu_s=0.0
  cpu_spos=0
  cpu_spos=0
  cpu_spos=0

  call init_rad_boundary()

  if(nrestart>0)then
     call restore_radiation()
  end if

  call observe_init()

end subroutine rad_init

subroutine rad_cpu_malloc()
  use data_common
  use radiation_commons
  implicit none
  integer::ncellx,ncelly,ncellz,nbnd;
  integer::nmem,nmemf;

  call aton_get_grid_size(ncellx, ncelly, ncellz, nbnd);
  nmem = (ncellx + 2*nbnd)*(ncelly + 2*nbnd)*(ncellz + 2*nbnd);
  nmemf = 3*nmem;

  allocate(cpu_e(0:nmem-1))
  allocate(cpu_d(0:nmem-1))
  allocate(cpu_t(0:nmem-1))
  allocate(cpu_x(0:nmem-1))
  allocate(cpu_photon_source(0:nmem-1))

  allocate(cpu_f(0:nmemf-1))

  allocate(cpu_s(0:rad_num_sources-1))
  allocate(cpu_spos(0:3*rad_num_sources-1))
end subroutine rad_cpu_malloc

subroutine clean_radiation
  implicit none
  call clean_rad_boundary()
  call observe_stop()
end subroutine clean_radiation

subroutine read_radiation_params(file_desc)
  use radiation_commons
  implicit none
  integer::file_desc

  namelist/radiation_params/&
       & allow_gpu_overload,radiation_feedback,&
       & rad_max_time_step, &
       & rad_light_speed_factor, &
       & rad_num_sources,rad_source_x,rad_source_y, &
       & rad_source_z,rad_source_rate, &
       & rad_flux_x,rad_flux_y,rad_flux_z,rad_density, &
       & rad_flux_at_x_min_boundary, &
       & rad_boundary_condition, &
       & rad_escape_fraction, &
       & rad_aton_version

  rewind(file_desc)
  read(file_desc,NML=radiation_params)

  if (rad_num_sources.gt.1) then
     write(*,*) "Sorry, rad_num_sources > 1 is not supported yet."
     ! NOTE(tstranex): It should be easy to fix this if needed.
     stop
  end if

  if ((rad_aton_version.ne.'gpu').and.(rad_aton_version.ne.'cpu')) then
     write(*,*) "Invalid value for rad_aton_version."
     stop
  end if
  write(*,*) "Using ATON version: ", rad_aton_version
  if (rad_aton_version.eq.'cpu') then
    write(*,*) "Warning: CPU version of ATON is experimental."
  end if

  c_light = 2.99792458e8 * rad_light_speed_factor

  call rad_init()

end subroutine
