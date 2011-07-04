module forsockets
    integer HEADER_FLAGS, HEADER_CALL_ID, HEADER_FUNCTION_ID, HEADER_CALL_COUNT, & 
        HEADER_INTEGER_COUNT, HEADER_LONG_COUNT, HEADER_FLOAT_COUNT, & 
        HEADER_DOUBLE_COUNT, HEADER_BOOLEAN_COUNT, HEADER_STRING_COUNT, & 
        HEADER_SIZE

    parameter (HEADER_FLAGS=1, HEADER_CALL_ID=2, HEADER_FUNCTION_ID=3, & 
        HEADER_CALL_COUNT=4, HEADER_INTEGER_COUNT=5, HEADER_LONG_COUNT=6, & 
        HEADER_FLOAT_COUNT=7, HEADER_DOUBLE_COUNT=8, & 
        HEADER_BOOLEAN_COUNT=9, HEADER_STRING_COUNT=10, & 
        HEADER_SIZE=10)

    interface
        subroutine receive_integers & 
            (ints, length) & 
            bind(c, name='forsockets_receive_integers')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: ints
            integer (c_int32_t), value :: length
        end subroutine receive_integers

        subroutine receive_longs & 
            (longs, length) & 
            bind(c, name='forsockets_receive_longs')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: longs
            integer (c_int32_t), value :: length
        end subroutine receive_longs

        subroutine receive_floats & 
            (floats, length) & 
            bind(c, name='forsockets_receive_floats')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: floats
            integer (c_int32_t), value :: length
        end subroutine receive_floats

        subroutine receive_doubles & 
            (doubles, length) & 
            bind(c, name='forsockets_receive_doubles')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: doubles
            integer (c_int32_t), value :: length
        end subroutine receive_doubles

        subroutine receive_booleans & 
            (booleans, length) & 
            bind(c, name='forsockets_receive_booleans')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: booleans
            integer (c_int32_t), value :: length
        end subroutine receive_booleans

        subroutine receive_string & 
            (string, length) & 
            bind(c, name='forsockets_receive_string')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: string
            integer (c_int32_t), value :: length
        end subroutine receive_string

        subroutine send_integers & 
            (ints, length) & 
            bind(c, name='forsockets_send_integers')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: ints
            integer (c_int32_t), value :: length
        end subroutine send_integers

        subroutine send_longs & 
            (longs, length) & 
            bind(c, name='forsockets_send_longs')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: longs
            integer (c_int32_t), value :: length
        end subroutine send_longs

        subroutine send_floats & 
            (floats, length) & 
            bind(c, name='forsockets_send_floats')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: floats
            integer (c_int32_t), value :: length
        end subroutine send_floats

        subroutine send_doubles & 
            (doubles, length) & 
            bind(c, name='forsockets_send_doubles')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: doubles
            integer (c_int32_t), value :: length
        end subroutine send_doubles

        subroutine send_booleans & 
            (booleans, length) & 
            bind(c, name='forsockets_send_booleans')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: booleans
            integer (c_int32_t), value :: length
        end subroutine send_booleans

        subroutine send_string & 
            (string, length) & 
            bind(c, name='forsockets_send_string')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: string
            integer (c_int32_t), value :: length
        end subroutine send_string

        subroutine forsockets_init & 
            (port) & 
            bind(c, name='forsockets_init')
            use iso_c_binding
            implicit none
            integer (c_int32_t), value :: port
        end subroutine forsockets_init

        subroutine forsockets_close & 
            () & 
            bind(c, name='forsockets_close')
            use iso_c_binding
            implicit none
        end subroutine forsockets_close

    end interface
end module forsockets



function internal__redirect_outputs(stdoutfile, stderrfile)
    use iso_c_binding
    
    implicit none
    
    character(kind=c_char, len = *) , intent(in) :: stdoutfile, stderrfile
    integer(c_int32_t) :: internal__redirect_outputs
    
    print*, 'NOT redirecting output to', stdoutfile, ' and ', stderrfile
    call flush()
    
    internal__redirect_outputs = 0
end function

subroutine run_loop
  use iso_c_binding
  use forsockets
  
  implicit none
  
  INCLUDE 'mpif.h'
  
  
  
  integer (c_int32_t) internal__redirect_outputs
  
  
  integer (c_int32_t) :: max_length = 255, MAX_STRING_LENGTH = 256
  logical :: must_run_loop, error
  integer i, length
  character (c_char), allocatable, target :: characters_in(:) * 256
  character (c_char), allocatable, target :: characters_out(:) * 256
  
  integer (c_int32_t), target :: header_in(HEADER_SIZE)
  integer (c_int32_t), target :: header_out(HEADER_SIZE)
  
  
  integer (c_int32_t), allocatable, target :: integers_in(:)
  integer (c_int32_t), allocatable, target :: integers_out(:)
  
  integer (c_int64_t), allocatable, target :: longs_in(:)
  integer (c_int64_t), allocatable, target :: longs_out(:)
  
  real (c_float), allocatable, target :: floats_in(:)
  real (c_float), allocatable, target :: floats_out(:)
  
  real (c_double), allocatable, target :: doubles_in(:)
  real (c_double), allocatable, target :: doubles_out(:)
  
  logical (c_bool), allocatable, target :: booleans_in(:)
  logical (c_bool), allocatable, target :: booleans_out(:)
  
  integer (c_int32_t), allocatable, target :: strings_in(:)
  integer (c_int32_t), allocatable, target :: strings_out(:)
  
  
  allocate(integers_in( max_length *7))
  allocate(integers_out( max_length * 2))
  allocate(doubles_in( max_length *26))
  allocate(doubles_out( max_length * 26))
  allocate(strings_in( max_length *2))
  
  must_run_loop = .true.
  
  do while (must_run_loop)
    
    call receive_integers(c_loc(header_in), HEADER_SIZE)
    
    length = header_in(HEADER_CALL_COUNT)
    
    if (length .gt. max_length) then
      max_length = length + 255;
      deallocate(integers_in)
      deallocate(integers_out)
      deallocate(strings_in)
      deallocate(doubles_in)
      deallocate(doubles_out)
      allocate(integers_in( max_length *7))
      allocate(integers_out( max_length * 2))
      allocate(doubles_in( max_length *26))
      allocate(doubles_out( max_length * 26))
      allocate(strings_in( max_length *2))
    end if
    
    if (header_in(HEADER_INTEGER_COUNT) .gt. 0) then
      call receive_integers(c_loc(integers_in), header_in(HEADER_INTEGER_COUNT))
    end if
    
    if (header_in(HEADER_LONG_COUNT) .gt. 0) then
      call receive_longs(c_loc(longs_in), header_in(HEADER_LONG_COUNT))
    end if
    
    if (header_in(HEADER_FLOAT_COUNT) .gt. 0) then
      call receive_floats(c_loc(floats_in), header_in(HEADER_FLOAT_COUNT))
    end if
    
    if (header_in(HEADER_DOUBLE_COUNT) .gt. 0) then
      call receive_doubles(c_loc(doubles_in), header_in(HEADER_DOUBLE_COUNT))
    end if
    
    if (header_in(HEADER_BOOLEAN_COUNT) .gt. 0) then
      call receive_booleans(c_loc(booleans_in), header_in(HEADER_BOOLEAN_COUNT))
    end if
    
    if (header_in(HEADER_STRING_COUNT) .gt. 0) then
      
      
      call receive_integers(c_loc(strings_in), header_in(HEADER_STRING_COUNT))

      !print*, 'received string header:', strings_in
      !call flush()


      do i = 1, header_in(HEADER_STRING_COUNT), 1
        if (strings_in(i) .gt. MAX_STRING_LENGTH) then
            print*, 'error! cannot receive strings exeeding length ', MAX_STRING_LENGTH
        end if
      end do

      !space for all strings in this one call
      allocate(characters_in(header_in(HEADER_STRING_COUNT)))

      do i = 1, header_in(HEADER_STRING_COUNT), 1
          characters_in(i) = ' '
          
          call receive_string(c_loc(characters_in(i)), strings_in(i))

          print*, 'received string:', characters_in(i)

          call flush()

      end do

    end if
    
    
    header_out = 0
    header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
    header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
    header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)
    
    error = .false.
    
    select case (header_in(HEADER_FUNCTION_ID))
      case(0)
        must_run_loop = .false.
      CASE(670175614)
        header_out(HEADER_INTEGER_COUNT) = 1 * length
        CALL initialize( &
          doubles_in(1) ,&
          doubles_in(2) ,&
          doubles_in(3) ,&
          doubles_in(4) ,&
          doubles_in(5) ,&
          doubles_in(6) ,&
          integers_in(1) ,&
          integers_in(2) ,&
          integers_in(3) ,&
          integers_in(4) ,&
          integers_in(5) ,&
          integers_in(6) ,&
          doubles_in(7) ,&
          integers_in(7) ,&
          doubles_in(8) ,&
          doubles_in(9) ,&
          doubles_in(10) ,&
          doubles_in(11) ,&
          doubles_in(12) ,&
          doubles_in(13) ,&
          doubles_in(14) ,&
          doubles_in(15) ,&
          doubles_in(16) ,&
          doubles_in(17) ,&
          integers_out(1) &
        )
        
      
      CASE(1024680297)
        header_out(HEADER_DOUBLE_COUNT) = 1 * length
        do i = 1, length, 1
          CALL get_time_step( &
            integers_in(i) ,&
            integers_in(( 1 * length) + i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * length) + i) ,&
            doubles_in(( 2 * length) + i) ,&
            doubles_in(( 3 * length) + i) ,&
            doubles_in(( 4 * length) + i) ,&
            doubles_in(( 5 * length) + i) ,&
            doubles_in(( 6 * length) + i) ,&
            doubles_in(( 7 * length) + i) ,&
            doubles_in(( 8 * length) + i) ,&
            doubles_out(i) &
          )
        end do
        
      
      CASE(1141573512)
        header_out(HEADER_INTEGER_COUNT) = 1 * length
        integers_out(1) = internal__redirect_outputs( &
          trim(characters_in(1)) ,&
          trim(characters_in(2)) &
        )
        
      
      CASE(1713719459)
        header_out(HEADER_INTEGER_COUNT) = 2 * length
        header_out(HEADER_DOUBLE_COUNT) = 26 * length
        do i = 1, length, 1
          CALL evolve_binary( &
            integers_in(i) ,&
            integers_in(( 1 * length) + i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * length) + i) ,&
            doubles_in(( 2 * length) + i) ,&
            doubles_in(( 3 * length) + i) ,&
            doubles_in(( 4 * length) + i) ,&
            doubles_in(( 5 * length) + i) ,&
            doubles_in(( 6 * length) + i) ,&
            doubles_in(( 7 * length) + i) ,&
            doubles_in(( 8 * length) + i) ,&
            doubles_in(( 9 * length) + i) ,&
            doubles_in(( 10 * length) + i) ,&
            doubles_in(( 11 * length) + i) ,&
            doubles_in(( 12 * length) + i) ,&
            doubles_in(( 13 * length) + i) ,&
            doubles_in(( 14 * length) + i) ,&
            doubles_in(( 15 * length) + i) ,&
            doubles_in(( 16 * length) + i) ,&
            doubles_in(( 17 * length) + i) ,&
            doubles_in(( 18 * length) + i) ,&
            doubles_in(( 19 * length) + i) ,&
            doubles_in(( 20 * length) + i) ,&
            doubles_in(( 21 * length) + i) ,&
            doubles_in(( 22 * length) + i) ,&
            doubles_in(( 23 * length) + i) ,&
            doubles_in(( 24 * length) + i) ,&
            doubles_in(( 25 * length) + i) &
          )
          integers_out(i) = integers_in(i)
          integers_out(( 1 * length) + i) = integers_in(( 1 * length) + i)
          doubles_out(i) = doubles_in(i)
          doubles_out(( 1 * length) + i) = doubles_in(( 1 * length) + i)
          doubles_out(( 2 * length) + i) = doubles_in(( 2 * length) + i)
          doubles_out(( 3 * length) + i) = doubles_in(( 3 * length) + i)
          doubles_out(( 4 * length) + i) = doubles_in(( 4 * length) + i)
          doubles_out(( 5 * length) + i) = doubles_in(( 5 * length) + i)
          doubles_out(( 6 * length) + i) = doubles_in(( 6 * length) + i)
          doubles_out(( 7 * length) + i) = doubles_in(( 7 * length) + i)
          doubles_out(( 8 * length) + i) = doubles_in(( 8 * length) + i)
          doubles_out(( 9 * length) + i) = doubles_in(( 9 * length) + i)
          doubles_out(( 10 * length) + i) = doubles_in(( 10 * length) + i)
          doubles_out(( 11 * length) + i) = doubles_in(( 11 * length) + i)
          doubles_out(( 12 * length) + i) = doubles_in(( 12 * length) + i)
          doubles_out(( 13 * length) + i) = doubles_in(( 13 * length) + i)
          doubles_out(( 14 * length) + i) = doubles_in(( 14 * length) + i)
          doubles_out(( 15 * length) + i) = doubles_in(( 15 * length) + i)
          doubles_out(( 16 * length) + i) = doubles_in(( 16 * length) + i)
          doubles_out(( 17 * length) + i) = doubles_in(( 17 * length) + i)
          doubles_out(( 18 * length) + i) = doubles_in(( 18 * length) + i)
          doubles_out(( 19 * length) + i) = doubles_in(( 19 * length) + i)
          doubles_out(( 20 * length) + i) = doubles_in(( 20 * length) + i)
          doubles_out(( 21 * length) + i) = doubles_in(( 21 * length) + i)
          doubles_out(( 22 * length) + i) = doubles_in(( 22 * length) + i)
          doubles_out(( 23 * length) + i) = doubles_in(( 23 * length) + i)
          doubles_out(( 24 * length) + i) = doubles_in(( 24 * length) + i)
          doubles_out(( 25 * length) + i) = doubles_in(( 25 * length) + i)
        end do
        
      
      case default
        error = .true.
    end select
    
    if (header_in(HEADER_STRING_COUNT) .gt. 0) then
      
      deallocate(characters_in)
    
    end if
    
    !print*, 'sending header', header_out
    
    !call flush()
    
    call send_integers(c_loc(header_out), HEADER_SIZE)
    
    if (header_out(HEADER_INTEGER_COUNT) .gt. 0) then
      call send_integers(c_loc(integers_out), header_out(HEADER_INTEGER_COUNT))
    end if
    
    if (header_out(HEADER_LONG_COUNT) .gt. 0) then
      call send_longs(c_loc(longs_out), header_out(HEADER_LONG_COUNT))
    end if
    
    if (header_out(HEADER_FLOAT_COUNT) .gt. 0) then
      call send_floats(c_loc(floats_out), header_out(HEADER_FLOAT_COUNT))
    end if
    
    if (header_out(HEADER_DOUBLE_COUNT) .gt. 0) then
      call send_doubles(c_loc(doubles_out), header_out(HEADER_DOUBLE_COUNT))
    end if
    
    if (header_out(HEADER_BOOLEAN_COUNT) .gt. 0) then
      call send_booleans(c_loc(booleans_out), header_out(HEADER_BOOLEAN_COUNT))
    end if
    
    if (header_out(HEADER_STRING_COUNT) .gt. 0) then
      
      !figure out length of all strings
      do i = 1, header_out(HEADER_STRING_COUNT), 1
          strings_out(i) = len_trim(characters_out(i))
      end do

      !send string header
      call send_integers(c_loc(strings_out), header_out(HEADER_STRING_COUNT))

      do i = 1, header_out(HEADER_STRING_COUNT), 1
          print*, 'sending string', characters_out(i)
          call flush()
          call send_string(c_loc(characters_out(i)), strings_out(i))
      end do     

      
      deallocate(characters_out)
    end if
    
  end do
  
  deallocate(integers_in)
  deallocate(integers_out)
  deallocate(strings_in)
  deallocate(doubles_in)
  deallocate(doubles_out)
  return
end subroutine


  program amuse_worker
    use iso_c_binding
    use forsockets
    
    implicit none
    
    include 'mpif.h'
    integer :: provided,ioerror, port
    character(len=32) :: port_string

    call mpi_init_thread(mpi_thread_multiple, provided, ioerror)

    call get_command_argument(1, port_string)

    read (port_string,*) port

    call forsockets_init(port)
    
    call run_loop()
    
    call mpi_finalize(ioerror)

    call forsockets_close()

  end program amuse_worker
