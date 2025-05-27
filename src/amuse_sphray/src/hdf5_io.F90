!> \file hdf5_io.f90

!> \brief Wrappers for HDF5 I/O 
!!
!! 
!< 

module hdf5_io_mod
use myf03_mod
use hdf5
implicit none
private


public :: rw_attr


interface rw_attr
   module procedure rw_attr_int,  rw_attr_int_arr, &
                    rw_attr_real, rw_attr_real_arr
end interface


integer, parameter :: hdf5_id_t = HID_T
integer, parameter :: hdf5_size_t = HSIZE_T

contains



! read/write a scalar, integer, group attribute
!------------------------------------------------------------------
subroutine rw_attr_int(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  integer(i4b) :: buf
  character(1) :: rorw
 
  integer(HSIZE_T), dimension(1), parameter :: dims = [1]
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_f(H5S_SCALAR_F, space_id, ierr) 
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_int'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_int


! read/write a vector, integer, group attribute
!------------------------------------------------------------------
subroutine rw_attr_int_arr(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  integer(i4b) :: buf(0:5)
  character(1) :: rorw

  integer(HSIZE_T), dimension(1), parameter :: dims = [6]
  integer(i4b), parameter :: rank = 1
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_simple_f(rank, dims, space_id, ierr)
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_int_arr'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_int_arr



! read/write a scalar, real, group attribute
!------------------------------------------------------------------
subroutine rw_attr_real(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  real(r8b) :: buf
  character(1) :: rorw
 
  integer(HSIZE_T), dimension(1), parameter :: dims = [1]
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_f(H5S_SCALAR_F, space_id, ierr) 
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_real'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_real


! read/write a vector, real, group attribute
!------------------------------------------------------------------
subroutine rw_attr_real_arr(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  real(r8b) :: buf(0:5)
  character(1) :: rorw

  integer(HSIZE_T), dimension(1), parameter :: dims = [6]
  integer(i4b), parameter :: rank = 1
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_simple_f(rank, dims, space_id, ierr)
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_real_arr'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_real_arr


end module hdf5_io_mod
