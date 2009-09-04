

SUBROUTINE run_loop
  INCLUDE 'mpif.h'
  integer :: rank, parent, ioerror
  integer :: must_run_loop
  integer mpiStatus(MPI_STATUS_SIZE,4)
  
  integer header(4)
  
  integer :: tag_in, tag_out
  integer integers_in(255)
  integer integers_out(255)
  integer :: number_of_integers_out, number_of_integers_in
  real*8 doubles_in(255)
  real*8 doubles_out(255)
  integer :: number_of_doubles_out, number_of_doubles_in
  real*4 floats_in(255)
  real*4 floats_out(255)
  integer :: number_of_floats_out_out, number_of_floats_out_in
  
  call MPI_COMM_GET_PARENT(parent, ioerror)
  call MPI_COMM_RANK(parent, rank, mpierror)
  
  must_run_loop = 1
  
  do while (must_run_loop .eq. 1)
    call MPI_RECV(header, 4, MPI_INTEGER, 0, 0, parent,&
      mpiStatus, ioerror)
    
    tag_in = header(1)
    number_of_doubles_in =  header(2)
    number_of_integers_in =  header(3)
    number_of_floats_in =  header(4)
    
    tag_out = tag_in
    number_of_doubles_out = 0
    number_of_integers_out = 0
    number_of_floats_out = 0
    
    if (number_of_doubles_in .gt. 0) then
      call MPI_RECV(doubles_in, number_of_doubles_in, &
        MPI_DOUBLE_PRECISION, 0, 0, parent,&
        mpiStatus, ioError);
    end if
    if (number_of_integers_in .gt. 0) then
      call MPI_RECV(integers_in, number_of_integers_in, &
        MPI_INTEGER, 0, 0, parent,&
        mpiStatus, ioError);
    end if
    if (number_of_floats_in .gt. 0) then
      call MPI_RECV(floats_in, number_of_floats_in, &
        MPI_SINGLE_PRECISION, 0, 0, parent,&
        mpiStatus, ioError);
    end if
    
    SELECT CASE (tag_in)
      CASE(0)
        must_run_loop = 0
      CASE(1)
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
        
      CASE(2)
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
        
      CASE(3)
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
        
      CASE DEFAULT
        tag_out = -1
    END SELECT
    
    header(1) = tag_out
    header(2) = number_of_doubles_out
    header(3) = number_of_integers_out
    
    call MPI_SEND(header, 4, MPI_INTEGER, 0, 999, &
      parent, mpierror);
    
    if (number_of_doubles_out .gt. 0) then
      call MPI_SEND(doubles_out, number_of_doubles_out, &
        MPI_DOUBLE_PRECISION, 0, 999, &
        parent, mpierror);
    end if
    if (number_of_integers_out .gt. 0) then
      call MPI_SEND(integers_out, number_of_integers_out, &
        MPI_INTEGER, 0, 999, &
        parent, mpierror);
    end if
    if (number_of_floats_out .gt. 0) then
      call MPI_SEND(floats_out, number_of_floats_out, &
        MPI_SINGLE_PRECISION, 0, 999, &
        parent, mpierror);
    end if
  end do
  return
end subroutine

program muse_worker
  INCLUDE 'mpif.h'
  call MPI_INIT(mpierror)
  
  call run_loop()
  
  call MPI_FINALIZE(mpierror)
end program muse_worker
