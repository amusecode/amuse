! Call once per process.
subroutine init_rad_boundary()
  use radiation_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::lattice_size

  if (rad_boundary_condition.eq.0) then
     if (myid.eq.1) then
        write(*,*)"Using zero-gradient boundary conditions for radiation."
     end if
  else if (rad_boundary_condition.eq.1) then
     if (myid.eq.1) then
        write(*,*)"Using periodic boundary conditions for radiation."
     end if
  else
     if (myid.eq.1) then
        write(*,*)"Invalid rad_boundary_condition:",rad_boundary_condition
     end if
     call clean_stop
  end if

  lattice_size = MAX(grid_size_x, MAX(grid_size_y, grid_size_z))
  boundary_memory_size = 4*lattice_size**2
  allocate(boundary_send(1:boundary_memory_size,1:6))
  allocate(boundary_recv(1:boundary_memory_size,1:6))
  allocate(boundary_request_status(1:MPI_STATUS_SIZE,1:6))

  boundary_send = 0.0
  boundary_recv = 0.0
end subroutine

! Call once per process.
subroutine clean_rad_boundary()
  use radiation_commons
  implicit none
  
  deallocate(boundary_send)
  deallocate(boundary_recv)
  deallocate(boundary_request_status)
end subroutine

! Return the MPI rank for the processes with the given cpu grid coordinates.
! The result is in the range [1, ncpu].
! If periodic=true, then out of bound grid coordinates are wrapped back inside.
! If periodic=false and the grid coordinates are out of bounds,
! then 0 is returned.
function calc_cpuid(ii,jj,kk,periodic)
  use radiation_commons
  implicit none
  integer::calc_cpuid
  integer::ii,jj,kk
  logical::periodic
  integer::i,j,k

  if (periodic) then
     i = mod(ii + num_cpu_x, num_cpu_x)
     j = mod(jj + num_cpu_y, num_cpu_y)
     k = mod(kk + num_cpu_z, num_cpu_z)
  else
     i = ii
     j = jj
     k = kk
  end if

  if ((i.lt.0).or.(i.ge.num_cpu_x).or. &
       & (j.lt.0).or.(j.ge.num_cpu_y).or. &
       & (k.lt.0).or.(k.ge.num_cpu_z)) then
     calc_cpuid = 0
  else
     calc_cpuid = 1 + i + j*num_cpu_x + k*num_cpu_x*num_cpu_y
  endif
end function

subroutine impose_flux_at_x_min_boundary()
  use radiation_commons
  implicit none
  integer::calc_cpuid
  integer::i,j
  integer::bound_index
  integer::lattice_size,stride
  logical::periodic

  ! Periodic boundary conditions:
  periodic = (rad_boundary_condition.eq.1)

  if (calc_cpuid(my_i-1, my_j, my_k, periodic).ne.0) then
     ! This is an internal boundary.
     return
  end if

  lattice_size = MAX(grid_size_x, MAX(grid_size_y, grid_size_z))
  stride = lattice_size**2

  do i=0,lattice_size-1
     do j=0,lattice_size-1
        bound_index = 1 + i + lattice_size*j
        boundary_recv(bound_index + 0*stride, 1) = rad_density
        boundary_recv(bound_index + 1*stride, 1) = rad_flux_x
        boundary_recv(bound_index + 2*stride, 1) = rad_flux_y
        boundary_recv(bound_index + 3*stride, 1) = rad_flux_z
     end do
  end do

end subroutine

subroutine cpu_mpi_boundary()
  use radiation_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer,dimension(6)::adjacent_cpu
  integer,dimension(6)::boundary_receive_request,boundary_send_request

  integer::calc_cpuid

  integer::i,j,k
  integer::global_i,global_j,global_k
  integer::i_min,j_min,k_min,i_max,j_max,k_max
  integer::icpu
  integer::global_size
  logical::periodic
  integer::tag_recv,tag_send

  ! Periodic boundary conditions:
  periodic = (rad_boundary_condition.eq.1)

  adjacent_cpu(1) = calc_cpuid(my_i-1, my_j, my_k, periodic)
  adjacent_cpu(2) = calc_cpuid(my_i+1, my_j, my_k, periodic)
  adjacent_cpu(3) = calc_cpuid(my_i, my_j-1, my_k, periodic)
  adjacent_cpu(4) = calc_cpuid(my_i, my_j+1, my_k, periodic)
  adjacent_cpu(5) = calc_cpuid(my_i, my_j, my_k-1, periodic)
  adjacent_cpu(6) = calc_cpuid(my_i, my_j, my_k+1, periodic)

  ! Receive boundary values from other cpus (non-blocking).
  n_receive_request=0
  do i=1,6
     if (adjacent_cpu(i).eq.0) cycle
     n_receive_request = n_receive_request + 1
     tag_recv = i
     call MPI_IRECV( &
          & boundary_recv(:,i), &
          & boundary_memory_size, &
          & MPI_REAL8, &
          & adjacent_cpu(i)-1, &
          & tag_recv, &
          & MPI_COMM_WORLD, &
          & boundary_receive_request(n_receive_request), &
          & info)
  end do

  ! Send boundary values to other cpus (non-blocking).
  n_send_request=0
  do i=1,6
     if (adjacent_cpu(i).eq.0) cycle
     n_send_request = n_send_request + 1

     ! Sometimes we need to get multiple boundary values from the same cpu.
     ! This can happen when using periodic boundary conditions in a small
     ! simulation (e.g. with only two CPUs).
     ! In this case, we need to use tags to distinguish the particular boundary.
     if (i == 1) tag_send = 2
     if (i == 2) tag_send = 1
     if (i == 3) tag_send = 4
     if (i == 4) tag_send = 3
     if (i == 5) tag_send = 6
     if (i == 6) tag_send = 5

     call MPI_ISEND( &
          & boundary_send(:,i), &
          & boundary_memory_size, &
          & MPI_REAL8, &
          & adjacent_cpu(i)-1, &
          & tag_send, &
          & MPI_COMM_WORLD, &
          & boundary_send_request(n_send_request), &
          & info)
  end do

  call MPI_WAITALL(n_receive_request, &
       & boundary_receive_request, &
       & boundary_request_status, &
       info)
  call MPI_WAITALL(n_send_request, &
       & boundary_send_request, &
       & boundary_request_status, &
       info)

  ! Impose boundary conditions on external boundaries.
  do i=1,6
     if (adjacent_cpu(i).eq.0) then
        ! If there is no adjacent cell, use zero-gradient boundary conditions.
        boundary_recv(:,i) = boundary_send(:,i)
     end if
  end do

  if (rad_flux_at_x_min_boundary) then
     call impose_flux_at_x_min_boundary()
  end if

end subroutine
