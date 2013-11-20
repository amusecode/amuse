! Copyright (c) 2010-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

!TODO: faire attention au communicateur

module hints
  use mpi

  implicit none

  private
  public :: hints_init
  public :: hints_finalize
  public :: hints_print
  public :: hints_put
  public :: hints_set
  public :: hints_reset

  integer,parameter,public :: MAX_HINTS_NUMBER = 28

  integer,parameter,public :: MAX_HINT_VALUELEN = 255

  integer,parameter :: NO = 0, YES = 1, UNKNOWN = -1
  integer,parameter :: INT = 0, STRING = 2, BOOL = 3, OTHER = 4

  logical,save :: is_initialized = .false.

  integer :: ierr

  type mpi_hint
    character(len=MPI_MAX_INFO_KEY) :: key
    integer :: value_type
    integer :: available = UNKNOWN
    logical :: used = .false.
    character(len=MAX_HINT_VALUELEN) :: current_value
    character(len=MAX_HINT_VALUELEN) :: default_value
  end type mpi_hint

  type(mpi_hint),dimension(:),allocatable :: hints_list



  contains

  subroutine hints_init(comm)
    integer,intent(in) :: comm

    if(is_initialized) then
      print *,'Error: hints_init already called'
      return
    endif

    call hints_build_standard_list
    call hints_get_list(comm)

    is_initialized = .true.

  end subroutine



  subroutine hints_finalize
    if(.not.is_initialized) then
      print *,'Error: hints_init never called'
      return
    endif

    deallocate(hints_list)

  end subroutine hints_finalize



  subroutine hints_build_standard_list
    logical,save :: standard_list_is_initialized = .false.

    integer :: idx

    if(standard_list_is_initialized) then
      print *,'Error: hints_build_standard_list already called'
      return
    endif

    allocate(hints_list(MAX_HINTS_NUMBER))

    idx = 1

    hints_list(idx)%key = "romio_cb_read"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "romio_cb_write"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "romio_cb_fr_types"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "romio_cb_fr_alignment"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "romio_cb_alltoall"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "romio_cb_pfr"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "romio_cb_ds_threshold"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "cb_buffer_size"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "cb_nodes"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "cb_config_list"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "romio_no_indep_rw"
    hints_list(idx)%value_type = BOOL
    idx = idx + 1

    hints_list(idx)%key = "ind_rd_buffer_size"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "ind_wr_buffer_size"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "romio_ds_read"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "romio_ds_write"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "filename"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    hints_list(idx)%key = "file_perm"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "IBM_io_buffer_size"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "IBM_largeblock_io"
    hints_list(idx)%value_type = BOOL
    idx = idx + 1

    hints_list(idx)%key = "IBM_sparse_access"
    hints_list(idx)%value_type = BOOL
    idx = idx + 1

    hints_list(idx)%key = "striping_unit"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "striping_factor"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "romio_lustre_start_iodevice"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "direct_read"
    hints_list(idx)%value_type = BOOL
    idx = idx + 1

    hints_list(idx)%key = "direct_write"
    hints_list(idx)%value_type = BOOL
    idx = idx + 1

    hints_list(idx)%key = "romio_lustre_co_ratio"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "romio_lustre_coll_threshold"
    hints_list(idx)%value_type = INT
    idx = idx + 1

    hints_list(idx)%key = "romio_lustre_ds_in_coll"
    hints_list(idx)%value_type = STRING
    idx = idx + 1

    standard_list_is_initialized = .true.
  end subroutine hints_build_standard_list



  subroutine hints_get_list(comm)
    integer,intent(in) :: comm

    logical::flag,match
    integer::i,j,infos,nkeys,sz
    integer::fid,mode
    character(len=MPI_MAX_INFO_VAL)::value
    character(len=MPI_MAX_INFO_KEY)::key

    mode = MPI_MODE_RDWR+MPI_MODE_CREATE+MPI_MODE_DELETE_ON_CLOSE
    call MPI_File_open(comm,"test_hint_get_list123456",mode,MPI_INFO_NULL,fid,ierr)
    if (ierr/=MPI_SUCCESS) then
      print *,'Error in hint_get_list'
      return
    end if

    call mpi_file_get_info(fid,infos,ierr)
    call mpi_info_get_nkeys(infos,nkeys,ierr)
    do i=1,nkeys
      match = .false.
      call mpi_info_get_nthkey(infos,i-1,key,ierr)
      call mpi_info_get_valuelen(infos,key,sz,flag,ierr)
      call mpi_info_get(infos,key,MAX_HINT_VALUELEN,value,flag,ierr)
      do j=1,MAX_HINTS_NUMBER
        if (trim(key)==trim(hints_list(j)%key)) then
          match = .true.
          hints_list(j)%available = YES
          if(sz>MAX_HINT_VALUELEN) then
            print *,'Error: MAX_HINT_VALUELEN is too small'
          else
            hints_list(j)%default_value = trim(value)
            hints_list(j)%current_value = trim(value)
          end if
          continue
        end if
      end do
      if(.not.match) then
        print *,'Error:',trim(key),' unknown'
      endif
    end do

    !Predefined-key availability is set to NO if not found
    do i = 1,MAX_HINTS_NUMBER
      if(hints_list(i)%available==UNKNOWN) hints_list(i)%available=NO
    end do

    call mpi_info_free(infos,ierr)

    call MPI_File_close(fid,ierr)

  end subroutine hints_get_list



  subroutine hints_put(key,value)
    character(len=*),intent(in)  :: key
    character(len=*),intent(in) :: value

    integer :: i
    logical :: match

    if(len(key)>MPI_MAX_INFO_KEY) then
      print *,'Error: key too long in hints_set'
      return
    end if
    if(len(value)>MAX_HINT_VALUELEN) then
      print *,'Error: value of key too long in hints_set'
      return
    end if

    match = .false.
    do i = 1,MAX_HINTS_NUMBER
      if(trim(key)==trim(hints_list(i)%key)) then
        match = .true.
        hints_list(i)%current_value = value
        hints_list(i)%used = .true.
      end if
    end do

    if(.not.match) print *,'Error: key not in the list (',trim(key),')'

  end subroutine hints_put



  subroutine hints_set(infos)
    integer,intent(inout)::infos

    integer :: i

    do i = 1,MAX_HINTS_NUMBER
      if(hints_list(i)%used) then
        call MPI_Info_set(infos,trim(hints_list(i)%key),trim(hints_list(i)%current_value),ierr)
      end if
    end do
  end subroutine hints_set


  subroutine hints_reset
    integer :: i

    if(.not.is_initialized) then
      print *,'Error: hints_init never called'
      return
    endif

    do i = 1,MAX_HINTS_NUMBER
      if(hints_list(i)%used) then
        hints_list(i)%current_value = hints_list(i)%default_value
      end if
    end do
  end subroutine hints_reset


  subroutine hints_print(comm)
    integer,intent(in) :: comm

    integer           :: i, rank
    character(len=16) :: value_type

    if(.not.is_initialized) then
      print *,'Error: hints_init never called'
      return
    endif

    call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

    if(rank==0) then
      write(*,'("--------------------------------------------------------------------------------------")')
      write(*,'("         MPI hints")')
      write(*,'("--------------------------------------------------------------------------------------")')
      write(*,'(" Key                      |   Avail | Used  | Type    | Current value | Default value ")')
      write(*,'("--------------------------------------------------------------------------------------")')
      do i = 1,MAX_HINTS_NUMBER
        if(hints_list(i)%value_type==INT) then
          value_type='integer'
        else if(hints_list(i)%value_type==STRING) then
          value_type='string'
        else if(hints_list(i)%value_type==BOOL) then
          value_type='bool'
        else if(hints_list(i)%value_type==OTHER) then
          value_type='other'
        else
          print *,'Error in hints_print: unknown value_type in hint'
          value_type=''
        end if

        if(hints_list(i)%available==YES) then
          write(*,'(" ",A24," | ",A7," | ",L5," | ",A7," | ",A13," | ",A13)') &
               hints_list(i)%key,"yes",hints_list(i)%used,value_type, &
               hints_list(i)%current_value,hints_list(i)%default_value
        else if(hints_list(i)%available==NO) then
          write(*,'(" ",A24," | ",A7," | ",L5," | ",A7," |               |")') &
               hints_list(i)%key,"no",hints_list(i)%used,value_type
        else if(hints_list(i)%available==UNKNOWN) then
          write(*,'(" ",A24," | ",A7," | ",L5," | ",A7," |               |")') &
               hints_list(i)%key,"unknown",hints_list(i)%used,value_type
        else
          print *,'Error in hints_print: unknown available flag in hint'
        end if
      end do
      write(*,'("--------------------------------------------------------------------------------------")')
    end if
  end subroutine hints_print
end module hints
