SUBROUTINE run_loop()

    INCLUDE 'mpif.h'

    integer :: rank, parent, ioerror
    integer :: must_run_loop
    integer mpiStatus(MPI_STATUS_SIZE,4)

    integer header(3)
    
    real*8 doubles_in(20)
    real*8 doubles_out(20)
    integer integers_in(20)
    integer integers_out(20)
    integer :: number_of_doubles_in, number_of_integers_in, tag_in
    integer :: number_of_doubles_out, number_of_integers_out, tag_out
    
    call MPI_COMM_GET_PARENT(parent, ioerror)
    call MPI_COMM_RANK(parent, rank, mpierror)

    must_run_loop = 1

    do while (must_run_loop .eq. 1)
        call MPI_RECV(header, 3, MPI_INTEGER, 0, rank, parent,&
                mpiStatus, ioerror)
        tag_in = header(1)
        number_of_doubles_in = header(2)
        number_of_integers_in = header(3)

    
        tag_out = tag_in
        number_of_doubles_out = 0
        number_of_integers_out = 0
        
        if (number_of_doubles_in .gt. 0) then
            call MPI_RECV(doubles_in, number_of_doubles_in,&
                 MPI_DOUBLE_PRECISION, 0,&
                 rank, parent,&
                 mpiStatus, ioerror)
        end if
      
        if (number_of_integers_in .gt. 0) then
            call MPI_RECV(integers_in, number_of_integers_in,&
                MPI_INTEGER, 0,&
                rank, parent,&
                mpiStatus, ioerror)
        end if
      
        
        !write (*,*) 'number_of_doubles_in=',number_of_doubles_in
      
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

          END SELECT
      
      
        header(1) = tag_out
        header(2) = number_of_doubles_out
        header(3) = number_of_integers_out
     
        call MPI_SEND(header, 3, MPI_INTEGER, 0, 999, parent,&
           ioerror)
        if (number_of_doubles_out .gt. 0) then
            call MPI_SEND(doubles_out, number_of_doubles_out, &
                MPI_DOUBLE_PRECISION, 0, &
                999, parent, &
                ioerror)
        end if
      
        if (number_of_integers_out .gt. 0) then
            call MPI_SEND(integers_out, number_of_integers_out, &
                MPI_INTEGER, 0, &
                999, parent, &
                ioerror)
        end if
    enddo

    return
END

program muse_worker

    INCLUDE 'mpif.h'
    
    write (*,*) 'starting'
    call MPI_INIT(mpierror)
    call run_loop()
    call MPI_FINALIZE(mpierror)
    write (*,*) 'stopping'
end program muse_worker
    
