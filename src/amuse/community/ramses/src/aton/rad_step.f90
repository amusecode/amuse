! Radiation module for RAMSES.
! It uses the ATON library.
! 2010
! Timothy Stranex <timothy@stranex.com>

! FIXME: Rename data_common to rad_arrays and move to its own file.
module data_common
  ! Scalar quantities.
  ! cpu_t: Temperature [Kelvin]
  ! cpu_d: Gas number density [atoms / m^3]
  ! cpu_x: Ionization fraction [dimensionless]
  ! cpu_e: photon number density [photons / m^3]
  ! cpu_photon_source: photon source field [photons / m^3 / s]
  real(kind=8),dimension(:),allocatable:: cpu_e,cpu_d,cpu_t,cpu_x
  real(kind=8),dimension(:),allocatable:: cpu_photon_source

  ! 3-vector quantities
  ! cpu_f: photon flux [photons / m^2 / s]
  real(kind=8),dimension(:),allocatable:: cpu_f

  ! Emission rates of the sources. [photons / m^3 / s]
  ! length is 0:nsrc-1
  real(kind=8),dimension(:),allocatable:: cpu_s
  ! Source positions.
  ! length is 0:nsrc-1.
  integer,dimension(:),allocatable::cpu_spos
end module data_common


subroutine rad_step(time_step_user)
  use amr_parameters
  use radiation_commons
  use cooling_module
  use timing
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  real(dp)::time_step_user  ! Hydro time step in user units.
  real(dp)::time_step_s     ! Hydro time step in seconds.
  real(kind=8)::dx          ! Lattice spacing in metres.

  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::hubblet  ! H [1/seconds]

  real(kind=8)::deltat  ! Radiation time step in seconds
  integer::num_steps    ! Number of radiation steps.
  integer::i            ! Radiation step counter.
  
  integer::ierr

  ! 1.0 corresponds to the maximum Courant condition.
  real(kind=8)::cfl=0.9

  
  call timer_start(total_timer)
  call timer_stop(ramses_timer)
  
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  myid = myid + 1  ! We need a 1-based id.
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  time_step_s = time_step_user * scale_t
  dx=boxlen*scale_l/100 / (2**levelmin)
  hubblet=h0*3.2408608e-20*HsurH0(1.0/aexp-1.,omega_m,omega_l,1.-omega_m-omega_l)

  ! Choose the radiation time step.
  deltat=cfl*dx/c_light/3.
  num_steps=time_step_s/deltat
  num_steps=max(num_steps, 1)
  
  if (myid.eq.1) then
     write(*,*) 'Starting ATON radiation step:'
     write(*,*) '  time step = ', time_step_s, ' (seconds)'
     write(*,*) '  dx = ', dx, ' (metres)'
     write(*,*) '  hubblet = ', hubblet, ' (1/seconds)'
     write(*,*) '  aexp = ', aexp
     write(*,*) '  J0min = ', J0min
     write(*,*) '  J0simple(aexp) = ', J0simple(aexp)
     write(*,*) '  deltat = ', deltat, ' (seconds)'
     write(*,*) '  Running ', num_steps, ' radiation steps.'
  end if

  call compute_Srad()

  call timer_start(mpi_timer)
  call start_mpi(dx)
  call timer_stop(mpi_timer)
  call timer_inc_count(mpi_timer, -1)


  call fill_cpu_field_from_hydro(dx)

  if (rad_aton_version.eq.'gpu') then
     call aton_gpu_loop(num_steps,dx,deltat)
  else
     call aton_cpu_loop(num_steps,dx,deltat)
  end if
    
  call timer_start(mpi_timer)
  call fill_hydro_from_cpu_field(dx)
  call end_mpi()
  call timer_stop(mpi_timer)

  ! Observations
  call observe_level(levelmin)
  call observe_output()

  ! Timings
  call timer_stop(total_timer)
  !call rad_output_timing()
  call timer_start(ramses_timer)

end subroutine rad_step

subroutine aton_cpu_loop(num_steps,dx,deltat)
  use radiation_commons
  use data_common
  implicit none
  integer::num_steps,i
  real(kind=8)::deltat
  real(kind=8)::dx

  !call aton_validate(myid, c_light, cpu_e, cpu_d, cpu_t, cpu_x, cpu_f)

  i=0
  do while(i.lt.num_steps)
     if(mod(i,1).eq.0) then
        !write(*,*)'.'
     end if

     call timer_start(boundary_memory_timer)
     call cpu_pack_boundary_values(cpu_e, cpu_f, boundary_send)
     call timer_stop(boundary_memory_timer)
     call timer_inc_count(boundary_memory_timer, -1)

     call timer_start(boundary_timer)
     call cpu_mpi_boundary()
     call timer_stop(boundary_timer)
     call timer_inc_count(boundary_timer, -1)

     call timer_start(boundary_memory_timer)
     call cpu_unpack_boundary_values(boundary_recv, cpu_e, cpu_f)
     call timer_stop(boundary_memory_timer)
     call timer_inc_count(boundary_memory_timer, -1)

     call timer_start(aton_timer)
     call aton_cpu_rad( &
          & myid, &
          & c_light,dx,deltat,rad_num_sources,fudgecool,aexp,0d0, &
          & cpu_e, cpu_d, cpu_t, cpu_x, cpu_photon_source, cpu_f, &
          & cpu_spos, cpu_s)
     call timer_stop(aton_timer)
     call timer_inc_count(aton_timer, -1)

     i=i+1
  end do

  !call aton_validate(myid, c_light, cpu_e, cpu_d, cpu_t, cpu_x, cpu_f)

end subroutine aton_cpu_loop

subroutine aton_gpu_loop(num_steps,dx,deltat)
  use radiation_commons
  use data_common
  implicit none
  integer::num_steps,i
  real(kind=8)::deltat
  real(kind=8)::dx

  call timer_start(full_memory_timer)
  call aton_cpu_to_gpu_full(cpu_e,cpu_f,cpu_x,cpu_t,cpu_d,cpu_s,cpu_spos,rad_num_sources,cpu_photon_source)
  call timer_stop(full_memory_timer)
  call timer_inc_count(full_memory_timer, -1)

  ! For debugging:
  !call aton_debug_dump(cpu_e, cpu_f, cpu_x, cpu_t, cpu_s, cpu_spos, cpu_d, rad_num_sources, cpu_photon_source, 0.0, 100+myid)
  !call aton_validate(10000+myid, c_light, cpu_e, cpu_d, cpu_t, cpu_x, cpu_f)

  !num_steps = 0
  
  i=0
  do while(i.lt.num_steps)
     if(mod(i,100).eq.0) then
        !write(*,*)'.'
     end if

     call timer_start(boundary_memory_timer)
     call aton_gpu_to_cpu_boundary(boundary_send)
     call timer_stop(boundary_memory_timer)
     call timer_inc_count(boundary_memory_timer, -1)

     call timer_start(boundary_timer)
     call cpu_mpi_boundary()
     call timer_stop(boundary_timer)
     call timer_inc_count(boundary_timer, -1)

     call timer_start(boundary_memory_timer)
     call aton_cpu_to_gpu_boundary(boundary_recv)
     call timer_stop(boundary_memory_timer)
     call timer_inc_count(boundary_memory_timer, -1)

     call timer_start(aton_timer)
     call aton_gpu_step(c_light,dx,deltat,rad_num_sources,fudgecool,aexp,0d0)
     call timer_stop(aton_timer)
     call timer_inc_count(aton_timer, -1)

     i=i+1
  end do

  call timer_inc_count(boundary_memory_timer, 1)
  call timer_inc_count(boundary_timer, 1)
  call timer_inc_count(aton_timer, 1)

  call timer_start(full_memory_timer)
  call aton_gpu_to_cpu_full(cpu_e,cpu_f,cpu_x,cpu_t,cpu_d,cpu_s,cpu_spos,rad_num_sources)
  call timer_stop(full_memory_timer)

  ! For debugging:
  !call aton_debug_dump(cpu_e, cpu_f, cpu_x, cpu_t, cpu_s, cpu_spos, cpu_d, rad_num_sources, cpu_photon_source, time_step_s, 200+myid)
  !call aton_validate(20000+myid, c_light, cpu_e, cpu_d, cpu_t, cpu_x, cpu_f)

end subroutine aton_gpu_loop

! Return the maximum total time step (sum of all substeps) for ATON in
! user units.
! -1 is returned if there is no maximum.
function aton_time_step()
  use radiation_commons
  implicit none
  real(dp)::aton_time_step
  aton_time_step = rad_max_time_step
end function

subroutine start_mpi(dx)
  use amr_commons
  use radiation_commons
  use data_common
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  real(kind=8)::dx

  integer::i,j,k
  integer::icpu
  integer::idx

  integer,dimension(1:nvector)::ind_tree,lev_tree
  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer,dimension(1:nvector)::cc

  allocate(sendbuf(1:ncpu))
  allocate(recvbuf(1:ncpu))
  allocate(icount(1:ncpu))
  allocate(receive_request(1:ncpu))
  allocate(send_request(1:ncpu))
  allocate(request_status(1:MPI_STATUS_SIZE,1:ncpu))

  if(ncpu.ne.num_cpu_x*num_cpu_y*num_cpu_z)then
     if(myid==1)then
        write(*,*)'Cartesian partition incorrect'
        write(*,*)' ncpu=',ncpu
        write(*,*)' cpu dimensions:',num_cpu_x,num_cpu_y,num_cpu_z
        write(*,*)'You need to either change the number of MPI nodes or'
        write(*,*)'change NCELLX, NCELLY, NCELLZ in ATON.'
     endif
     call clean_stop
  endif

  ! Find our cpu cell coordinates.
  my_k=(myid-1)/(num_cpu_x*num_cpu_y)
  my_j=(myid-1-my_k*num_cpu_x*num_cpu_y)/num_cpu_x
  my_i=(myid-1-my_j*num_cpu_x-my_k*num_cpu_x*num_cpu_y)

  ! Calculate corner coordinates of our subdomain.
  ! TODO(tstranex): Need to shift for boundaries.
  my_xmin=dble(my_i)/dble(num_cpu_x)
  my_ymin=dble(my_j)/dble(num_cpu_y)
  my_zmin=dble(my_k)/dble(num_cpu_z)
  my_xmax=dble(my_i+1)/dble(num_cpu_x)
  my_ymax=dble(my_j+1)/dble(num_cpu_y)
  my_zmax=dble(my_k+1)/dble(num_cpu_z)

  ! Calculate the sendbuf size for each cpu.
  allocate(cpu_rad(1:grid_size_x,1:grid_size_y,1:grid_size_z))
  sendbuf=0
  do i=1,grid_size_x
  do j=1,grid_size_y
  do k=1,grid_size_z
     xx_dp(1,1) = my_xmin + (dble(i)-0.5)/dble(2**levelmin)
     xx_dp(1,2) = my_ymin + (dble(j)-0.5)/dble(2**levelmin)
     xx_dp(1,3) = my_zmin + (dble(k)-0.5)/dble(2**levelmin)
     call cmp_cpumap(xx_dp,cc,1)
     cpu_rad(i,j,k)=cc(1)
     if(cc(1).ne.myid)then
        sendbuf(cc(1))=sendbuf(cc(1))+1
     endif
  end do
  end do
  end do

  ! Allocate send and receive buffers.
  do icpu=1,ncpu
     if(sendbuf(icpu).gt.0)then
        allocate(emission(icpu,1)%u_radiation(1:sendbuf(icpu),1:7))
     end if
  end do
  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)
  do icpu=1,ncpu
     if(recvbuf(icpu)>0)then
        allocate(reception(icpu,1)%u_radiation(1:recvbuf(icpu),1:7))
     end if
  end do

  ! Fill the send buffer with coordinates we need from others.
  icount=0
  do i=1,grid_size_x
  do j=1,grid_size_y
  do k=1,grid_size_z
     icpu = cpu_rad(i,j,k)
     if (icpu.eq.myid) cycle
     icount(icpu) = icount(icpu) + 1

     xx_dp(1,1) = my_xmin + (dble(i)-0.5)/dble(2**levelmin)
     xx_dp(1,2) = my_ymin + (dble(j)-0.5)/dble(2**levelmin)
     xx_dp(1,3) = my_zmin + (dble(k)-0.5)/dble(2**levelmin)

     emission(icpu,1)%u_radiation(icount(icpu),1)=xx_dp(1,1)
     emission(icpu,1)%u_radiation(icount(icpu),2)=xx_dp(1,2)
     emission(icpu,1)%u_radiation(icount(icpu),3)=xx_dp(1,3)
     emission(icpu,1)%u_radiation(icount(icpu),4)=0.0
     emission(icpu,1)%u_radiation(icount(icpu),5)=0.0
     emission(icpu,1)%u_radiation(icount(icpu),6)=0.0
     emission(icpu,1)%u_radiation(icount(icpu),7)=0.0
  end do
  end do
  end do

  ! Receive coordinates requests (non-blocking).
  n_receive_request=0
  do icpu=1,ncpu
     if (recvbuf(icpu).eq.0) cycle
     n_receive_request = n_receive_request + 1
     call MPI_IRECV( &
          & reception(icpu,1)%u_radiation, &
          & recvbuf(icpu)*7, &
          & MPI_DOUBLE_PRECISION, &
          & icpu-1, &
          & tag, &
          & MPI_COMM_WORLD, &
          & receive_request(n_receive_request), &
          & info)
  end do

  ! Send coordinate requests (non-blocking).
  n_send_request=0
  do icpu=1,ncpu
     if (sendbuf(icpu).eq.0) cycle
     n_send_request = n_send_request + 1
     call MPI_ISEND( &
          & emission(icpu,1)%u_radiation, &
          & sendbuf(icpu)*7, &
          & MPI_DOUBLE_PRECISION, &
          & icpu-1, &
          & tag, &
          & MPI_COMM_WORLD, &
          & send_request(n_send_request), &
          & info)
  end do

  call MPI_WAITALL(n_receive_request,receive_request,request_status,info)
  call MPI_WAITALL(n_send_request,send_request,request_status,info)

  ! Now we know which values the other cpus want. Lookup the requested values.
  do icpu=1,ncpu
     if (recvbuf(icpu).eq.0) cycle
     do i=1,recvbuf(icpu)
        xx_dp(1,1)=reception(icpu,1)%u_radiation(i,1)
        xx_dp(1,2)=reception(icpu,1)%u_radiation(i,2)
        xx_dp(1,3)=reception(icpu,1)%u_radiation(i,3)
        call aton_get_cell_index(ind_tree,lev_tree,xx_dp,levelmin,1)
        call get_rad_quantities_from_cell( &
             & ind_tree(1), &
             & dx, &
             & reception(icpu,1)%u_radiation(i,4), &
             & reception(icpu,1)%u_radiation(i,5), &
             & reception(icpu,1)%u_radiation(i,6), &
             & reception(icpu,1)%u_radiation(i,7))
     end do
  end do

  ! Receive our requested values from other cpus (non-blocking).
  n_receive_request=0
  do icpu=1,ncpu
     if (sendbuf(icpu).eq.0) cycle
     n_receive_request = n_receive_request + 1
     call MPI_IRECV( &
          & emission(icpu,1)%u_radiation, &
          & sendbuf(icpu)*7, &
          & MPI_DOUBLE_PRECISION, &
          & icpu-1, &
          & tag, &
          & MPI_COMM_WORLD, &
          & receive_request(n_receive_request), &
          & info)
  end do

  ! Send values to other cpus (non-blocking).
  n_send_request=0
  do icpu=1,ncpu
     if (recvbuf(icpu).eq.0) cycle
     n_send_request = n_send_request + 1
     call MPI_ISEND( &
          & reception(icpu,1)%u_radiation, &
          & recvbuf(icpu)*7, &
          & MPI_DOUBLE_PRECISION, &
          & icpu-1, &
          & tag, &
          & MPI_COMM_WORLD, &
          & send_request(n_send_request), &
          & info)
  end do

  call MPI_WAITALL(n_receive_request,receive_request,request_status,info)
  call MPI_WAITALL(n_send_request,send_request,request_status,info)

end subroutine


subroutine end_mpi()
  use radiation_commons
  implicit none

  integer::icpu

  deallocate(cpu_rad)

  ! Deallocate send and receive buffers.
  do icpu=1,ncpu
     if (sendbuf(icpu).gt.0)then
        deallocate(emission(icpu,1)%u_radiation)
     end if
     if (recvbuf(icpu).gt.0)then
        deallocate(reception(icpu,1)%u_radiation)
     end if
  end do

  deallocate(sendbuf)
  deallocate(recvbuf)
  deallocate(icount)
  deallocate(receive_request)
  deallocate(send_request)
  deallocate(request_status)

end subroutine


subroutine get_rad_quantities_from_cell(cell_index,dx,density,xion,temperature,radsource)
  use amr_commons
  use cooling_module
  use hydro_commons
  use radiation_commons
  implicit none

  integer::cell_index
  real(kind=8)::dx
  real(dp)::density,xion,temperature,radsource
  integer::idim
  real(dp)::ekk,T2
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Density rho.
  ! cudaton needs units in atoms/m^3
  density = uold(cell_index,1)*scale_d*0.76/mH*1.0e6
  if (density.eq.0d0) then
     write(*,*)'Warning: zero density'
  end if

  ! Ionization fraction x.
  xion = uold(cell_index,ixion)/uold(cell_index,1)

  ! Temperature.
  ! cudaton needs units in Kelvin
  ! TODO(tstranex): Make this more readable.
  T2=uold(cell_index,5)
  ekk=0.0
  do idim=1,ndim
     ekk=ekk+0.5*uold(cell_index,idim+1)**2/uold(cell_index,1)
  end do
  T2=(gamma-1.0)*(T2-ekk)
  T2=T2/uold(cell_index,1)*scale_T2
  temperature = T2 / (1 + xion)

  radsource = Srad(cell_index)/dx/dx/dx
end subroutine


recursive subroutine update_rad_quantities_in_cell(cell_index,photon_density,xion,temperature)
  use amr_commons
  use cooling_module
  use hydro_commons
  use radiation_commons
  implicit none

  integer::cell_index
  real(dp)::photon_density,xion,temperature
  integer::idim
  real(dp)::ekk,T2
  real(dp)::N0
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  integer::ichild

  ! N [m^-3] = 4*pi/hc 10^4 J0 [erg s^-1 cm^-2 ster^-1 Hz^-1]
  N0=6.326e22
  if (rad_aton_version.eq.'gpu2') then
     N0 = N0 * 2.99792458e8
  end if

  ! For reduced speed of light:
  ! c * N is independent of the speed of light so we need to scale N0 by c when c != 1.
  N0 = N0 / rad_light_speed_factor

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Ionization fraction x.
  uold(cell_index,ixion) = xion * uold(cell_index,1)

  ! Photon density.
  Erad(cell_index) = photon_density / N0

  ! Temperature.
  ! We only use the temperature from ATON if ramses cooling is switched
  ! off.
  ! TODO(tstranex): Make this more readable.
  if(.not.cooling)then
     ekk=0.0
     do idim=1,ndim
        ekk=ekk+0.5*uold(cell_index,idim+1)**2/uold(cell_index,1)
     end do
     T2 = temperature * (1 + xion)
     T2 = T2/scale_T2*uold(cell_index,1)
     T2 = T2/(gamma-1.0) + ekk
     uold(cell_index,5) = T2
  end if

  ! Recursively fill the subtree.
  if (son(cell_index).ne.0) then
     do ichild=0,7
        call update_rad_quantities_in_cell( &
             & son(cell_index) + ncoarse + ichild*ngridmax, &
             & photon_density, &
             & xion, &
             & temperature)
     end do
  end if

end subroutine

subroutine fill_cpu_field_from_hydro(dx)
  use data_common
  use amr_commons
  use hydro_commons
  use cooling_module
  use radiation_commons
  implicit none

  real(kind=8)::dx  ! Lattice spacing in metres.

  integer::i,j,k
  integer::idx
  integer::icpu
  integer::aton_cell_index

  integer::i_min,j_min,k_min

  real(dp)::density_dp, xion_dp, temperature_dp, radsource_dp
  real(dp)::ekk,T2
  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer,dimension(1:nvector)::ind_tree,lev_tree

  i_min = my_i*(2**levelmin)/num_cpu_x
  j_min = my_j*(2**levelmin)/num_cpu_y
  k_min = my_k*(2**levelmin)/num_cpu_z

  if(rad_num_sources.gt.1)then
     if (myid.eq.1) then
        write(*,*) "rad_num_sources > 1 is not supported"
     end if
     call clean_stop
  endif

  if (rad_num_sources.eq.1) then
     cpu_spos(0) = rad_source_x*2**levelmin - i_min
     cpu_spos(1) = rad_source_y*2**levelmin - j_min
     cpu_spos(2) = rad_source_z*2**levelmin - k_min
     if ((cpu_spos(0).lt.0).or.(cpu_spos(0).ge.grid_size_x).or. &
          & (cpu_spos(1).lt.0).or.(cpu_spos(1).ge.grid_size_y).or. &
          & (cpu_spos(2).lt.0).or.(cpu_spos(2).ge.grid_size_z)) then
        cpu_s(0)=0.0
     else
        cpu_s(0)=rad_source_rate/dx/dx/dx
     end if
  endif

  icount=0
  do i=1,grid_size_x
  do j=1,grid_size_y
  do k=1,grid_size_z
     idx = aton_cell_index(i-1, j-1, k-1)
     icpu = cpu_rad(i,j,k)

     if (icpu.eq.myid) then
        xx_dp(1,1) = my_xmin + (dble(i)-0.5)/dble(2**levelmin)
        xx_dp(1,2) = my_ymin + (dble(j)-0.5)/dble(2**levelmin)
        xx_dp(1,3) = my_zmin + (dble(k)-0.5)/dble(2**levelmin)
        call aton_get_cell_index(ind_tree,lev_tree,xx_dp,levelmin,1)
        call get_rad_quantities_from_cell( &
             & ind_tree(1), &
             & dx, &
             & density_dp, &
             & xion_dp, &
             & temperature_dp, &
             & radsource_dp)
        ! cpu_* are single precision.
        cpu_d(idx) = density_dp
        cpu_x(idx) = xion_dp
        cpu_t(idx) = temperature_dp
        cpu_photon_source(idx) = radsource_dp
     else
        icount(icpu) = icount(icpu) + 1
        cpu_d(idx) = emission(icpu,1)%u_radiation(icount(icpu),4)
        cpu_x(idx) = emission(icpu,1)%u_radiation(icount(icpu),5)
        cpu_t(idx) = emission(icpu,1)%u_radiation(icount(icpu),6)
        cpu_photon_source(idx) = emission(icpu,1)%u_radiation(icount(icpu),7)
     end if

  end do
  end do
  end do

  ! Convert the photon density and flux from comoving to physical units.
  ! Both N and F scale as a^-3.
  ! See Eq. A3 from Gnedin and Ostriker, 1997, ApJ, 486.
  cpu_e = cpu_e / (aexp**3)
  cpu_f = cpu_f / (aexp**3)

  if (.not.radiation_feedback) then
     cpu_t = 1d4
  end if

end subroutine fill_cpu_field_from_hydro


subroutine fill_hydro_from_cpu_field(dx)
  use data_common
  use amr_commons
  use hydro_commons
  use cooling_module
  use radiation_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  real(kind=8)::dx  ! Lattice spacing in metres.

  integer::i,j,k
  integer::idx
  integer::icpu
  integer::aton_cell_index

  real(dp)::photon_density_dp, xion_dp, temperature_dp

  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer,dimension(1:nvector)::ind_tree,lev_tree

  icount=0
  do i=1,grid_size_x
  do j=1,grid_size_y
  do k=1,grid_size_z
     idx = aton_cell_index(i-1, j-1, k-1)
     icpu = cpu_rad(i,j,k)
     if (icpu.eq.myid) cycle

     icount(icpu) = icount(icpu) + 1
     emission(icpu,1)%u_radiation(icount(icpu),4) = cpu_e(idx)
     emission(icpu,1)%u_radiation(icount(icpu),5) = cpu_x(idx)
     emission(icpu,1)%u_radiation(icount(icpu),6) = cpu_t(idx)
     emission(icpu,1)%u_radiation(icount(icpu),7) = 0.0  ! Unused.

  end do
  end do
  end do

  ! Receive updated values from other cpus (non-blocking).
  n_receive_request=0
  do icpu=1,ncpu
     if (recvbuf(icpu).eq.0) cycle
     n_receive_request = n_receive_request + 1
     call MPI_IRECV( &
          & reception(icpu,1)%u_radiation, &
          & recvbuf(icpu)*6, &
          & MPI_DOUBLE_PRECISION, &
          & icpu-1, &
          & tag, &
          & MPI_COMM_WORLD, &
          & receive_request(n_receive_request), &
          & info)
  end do

  ! Send updated values to other cpus (non-blocking).
  n_send_request=0
  do icpu=1,ncpu
     if (sendbuf(icpu).eq.0) cycle
     n_send_request = n_send_request + 1
     call MPI_ISEND( &
          & emission(icpu,1)%u_radiation, &
          & sendbuf(icpu)*6, &
          & MPI_DOUBLE_PRECISION, &
          & icpu-1, &
          & tag, &
          & MPI_COMM_WORLD, &
          & send_request(n_send_request), &
          & info)
  end do

  call MPI_WAITALL(n_receive_request,receive_request,request_status,info)
  call MPI_WAITALL(n_send_request,send_request,request_status,info)


  ! Need to update our cells from radiation values computed on other cpus.
  do icpu=1,ncpu
     if (recvbuf(icpu).eq.0) cycle

     do i=1,recvbuf(icpu)
        xx_dp(1,1)=reception(icpu,1)%u_radiation(i,1)
        xx_dp(1,2)=reception(icpu,1)%u_radiation(i,2)
        xx_dp(1,3)=reception(icpu,1)%u_radiation(i,3)
        call aton_get_cell_index(ind_tree,lev_tree,xx_dp,levelmin,1)

        call update_rad_quantities_in_cell( &
             & ind_tree(1), &
             & reception(icpu,1)%u_radiation(i,4), &
             & reception(icpu,1)%u_radiation(i,5), &
             & reception(icpu,1)%u_radiation(i,6))

     end do
  end do

  ! Update our cells.
  do i=1,grid_size_x
  do j=1,grid_size_y
  do k=1,grid_size_z
     idx = aton_cell_index(i-1, j-1, k-1)
     icpu = cpu_rad(i,j,k)
     if (icpu.ne.myid) cycle

     xx_dp(1,1) = my_xmin + (dble(i)-0.5)/dble(2**levelmin)
     xx_dp(1,2) = my_ymin + (dble(j)-0.5)/dble(2**levelmin)
     xx_dp(1,3) = my_zmin + (dble(k)-0.5)/dble(2**levelmin)
     call aton_get_cell_index(ind_tree,lev_tree,xx_dp,levelmin,1)

     photon_density_dp = cpu_e(idx)
     xion_dp = cpu_x(idx)
     temperature_dp = cpu_t(idx)
     call update_rad_quantities_in_cell( &
          & ind_tree(1), &
          & photon_density_dp, &
          & xion_dp, &
          & temperature_dp)

  end do
  end do
  end do

  ! Convert the photon density and flux from physical to comoving units.
  cpu_e = cpu_e * (aexp**3)
  cpu_f = cpu_f * (aexp**3)

end subroutine

! Copied from Romain's radiation.f90
! Returns the index of the smallest AMR cell containing the given point.
subroutine aton_get_cell_index(cell_index,cell_levl,xpart,ilevel,np)
  use amr_commons
  implicit none
  integer::np,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart
  ! This function returns the index of the cell, at maximum level
  ! ilevel, in which the input particle sits
  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     if (myid.eq.1) then
        write(*,*)"nx=ny=nz != 1,3 is not supported."
     end if
     call clean_stop
  end if

  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,np
     xx = xpart(i,1) + (nx-1)/2.0
     yy = xpart(i,2) + (ny-1)/2.0
     zz = xpart(i,3) + (nz-1)/2.0
     igrid=igrid0
     do j=1,ilevel
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=ilevel
  end do
end subroutine aton_get_cell_index

subroutine rad_output_timing()
  use radiation_commons
  implicit none

  ! TODO: Use MPI to max the timers over all tasks.
  if (myid.eq.1) then
    write(*,*) 'Radiation timers:'
    write(*,*) ' rad_total_timer', total_timer%sum, total_timer%count
    write(*,*) ' rad_full_memory_timer', full_memory_timer%sum, full_memory_timer%count
    write(*,*) ' rad_boundary_memory_timer', boundary_memory_timer%sum, boundary_memory_timer%count
    write(*,*) ' rad_boundary_timer', boundary_timer%sum, boundary_timer%count
    write(*,*) ' rad_aton_timer', aton_timer%sum, aton_timer%count
    write(*,*) ' rad_mpi_timer', mpi_timer%sum, mpi_timer%count
    write(*,*) ' rad_ramses_timer', ramses_timer%sum, ramses_timer%count
  end if

  call timer_init(total_timer)
  call timer_init(ramses_timer)
  call timer_init(full_memory_timer)
  call timer_init(boundary_memory_timer)
  call timer_init(boundary_timer)
  call timer_init(aton_timer)
  call timer_init(mpi_timer)

end subroutine rad_output_timing
