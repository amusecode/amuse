

SUBROUTINE run_loop
  INCLUDE 'mpif.h'
  integer :: rank, parent, ioerror
  integer :: must_run_loop
  integer i
  integer mpiStatus(MPI_STATUS_SIZE,4)
  
  integer header(5)
  
  integer :: tag_in, tag_out
  
  integer :: len_in, len_out
  
  integer, DIMENSION(:), ALLOCATABLE ::integers_in
  integer, DIMENSION(:), ALLOCATABLE ::integers_out
  integer :: number_of_integers_out, number_of_integers_in
  real*4, DIMENSION(:), ALLOCATABLE ::floats_in
  real*4, DIMENSION(:), ALLOCATABLE ::floats_out
  integer :: number_of_floats_out, number_of_floats_in
  real*8, DIMENSION(:), ALLOCATABLE ::doubles_in
  real*8, DIMENSION(:), ALLOCATABLE ::doubles_out
  integer :: number_of_doubles_out, number_of_doubles_in
  
  ALLOCATE(integers_in(1275000))
  ALLOCATE(integers_out(1275000))
  ALLOCATE(floats_in(1275000))
  ALLOCATE(floats_out(1275000))
  ALLOCATE(doubles_in(1275000))
  ALLOCATE(doubles_out(1275000))
  
  call MPI_COMM_GET_PARENT(parent, ioerror)
  call MPI_COMM_RANK(parent, rank, mpierror)
  
  must_run_loop = 1
  
  do while (must_run_loop .eq. 1)
    call MPI_BCast(header, 5, MPI_INTEGER, 0, parent,&
      ioerror)
    
    tag_in = header(1)
    
    len_in = header(2)
    number_of_doubles_in =  header(3)
    number_of_integers_in =  header(4)
    number_of_floats_in =  header(5)
    
    tag_out = tag_in
    len_out = len_in
    number_of_doubles_out = 0
    number_of_integers_out = 0
    number_of_floats_out = 0
    
    if (number_of_doubles_in .gt. 0) then
      call MPI_BCast(doubles_in, number_of_doubles_in * len_in, &
        MPI_DOUBLE_PRECISION, 0, parent,&
        ioError);
    end if
    if (number_of_integers_in .gt. 0) then
      call MPI_BCast(integers_in, number_of_integers_in * len_in, &
        MPI_INTEGER, 0, parent,&
        ioError);
    end if
    if (number_of_floats_in .gt. 0) then
      call MPI_BCast(floats_in, number_of_floats_in * len_in, &
        MPI_SINGLE_PRECISION, 0, parent,&
        ioError);
    end if
    
    SELECT CASE (tag_in)
      CASE(0)
        must_run_loop = 0
      CASE(670175614)
        CALL initialize( &
          doubles_in(1) ,&
          doubles_in(2) ,&
          doubles_in(3) ,&
          doubles_in(4) ,&
          doubles_in(5) ,&
          integers_in(1) ,&
          integers_in(2) ,&
          integers_in(3) ,&
          integers_in(4) ,&
          integers_in(5) ,&
          doubles_in(6) ,&
          doubles_in(7) ,&
          doubles_in(8) ,&
          integers_out(1) &
        )
        number_of_integers_out = 1
        
      CASE(1024680297)
        CALL get_time_step( &
          integers_in(1) ,&
          doubles_in(1) ,&
          doubles_in(2) ,&
          doubles_in(3) ,&
          doubles_in(4) ,&
          doubles_in(5) ,&
          doubles_out(1) &
        )
        number_of_doubles_out = 1
        
      CASE(1658568341)
        CALL evolve0( &
          integers_in(1) ,&
          doubles_in(1) ,&
          doubles_in(2) ,&
          doubles_in(3) ,&
          doubles_in(4) ,&
          doubles_in(5) ,&
          doubles_in(6) ,&
          doubles_in(7) ,&
          doubles_in(8) ,&
          doubles_in(9) ,&
          doubles_in(10) ,&
          doubles_in(11) ,&
          doubles_in(12) ,&
          doubles_in(13) &
        )
        integers_out(1) = integers_in(1)
        doubles_out(1) = doubles_in(1)
        doubles_out(2) = doubles_in(2)
        doubles_out(3) = doubles_in(3)
        doubles_out(4) = doubles_in(4)
        doubles_out(5) = doubles_in(5)
        doubles_out(6) = doubles_in(6)
        doubles_out(7) = doubles_in(7)
        doubles_out(8) = doubles_in(8)
        doubles_out(9) = doubles_in(9)
        doubles_out(10) = doubles_in(10)
        doubles_out(11) = doubles_in(11)
        doubles_out(12) = doubles_in(12)
        doubles_out(13) = doubles_in(13)
        number_of_integers_out = 1
        number_of_doubles_out = 13
        
      CASE DEFAULT
        tag_out = -1
    END SELECT
    
    header(1) = tag_out
    header(2) = len_out
    header(3) = number_of_doubles_out
    header(4) = number_of_integers_out
    header(5) = number_of_floats_out
    
    call MPI_SEND(header, 5, MPI_INTEGER, 0, 999, &
      parent, mpierror);
    
    if (number_of_doubles_out .gt. 0) then
      call MPI_SEND(doubles_out, number_of_doubles_out * len_out, &
        MPI_DOUBLE_PRECISION, 0, 999, &
        parent, mpierror);
    end if
    if (number_of_integers_out .gt. 0) then
      call MPI_SEND(integers_out, number_of_integers_out * len_out, &
        MPI_INTEGER, 0, 999, &
        parent, mpierror);
    end if
    if (number_of_floats_out .gt. 0) then
      call MPI_SEND(floats_out, number_of_floats_out * len_out, &
        MPI_SINGLE_PRECISION, 0, 999, &
        parent, mpierror);
    end if
  end do
  
  DEALLOCATE(integers_in)
  DEALLOCATE(integers_out)
  DEALLOCATE(floats_in)
  DEALLOCATE(floats_out)
  DEALLOCATE(doubles_in)
  DEALLOCATE(doubles_out)
  return
end subroutine

program muse_worker
  INCLUDE 'mpif.h'
  call MPI_INIT(mpierror)
  
  call run_loop()
  
  call MPI_FINALIZE(mpierror)
end program muse_worker
