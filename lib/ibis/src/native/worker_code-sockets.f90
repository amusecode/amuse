module forsockets
	integer HEADER_FLAGS, HEADER_CALL_ID, HEADER_FUNCTION_ID, HEADER_LENGTH, &
		HEADER_INT_COUNT, HEADER_LONG_COUNT, HEADER_FLOAT_COUNT, &
		HEADER_DOUBLE_COUNT, HEADER_BOOLEAN_COUNT, HEADER_STRING_COUNT, &
		HEADER_SIZE

	parameter (HEADER_FLAGS=1, HEADER_CALL_ID=2, HEADER_FUNCTION_ID = 3, &
    	HEADER_LENGTH = 4, HEADER_INT_COUNT = 5, HEADER_LONG_COUNT = 6, &
    	HEADER_FLOAT_COUNT = 7, HEADER_DOUBLE_COUNT = 8, &
    	HEADER_BOOLEAN_COUNT = 9, HEADER_STRING_COUNT = 10, &
    	HEADER_SIZE = 10)

	interface
		integer (c_int32_t) function receive_ints &
			(ints, length) &
			bind(c, name='forsockets_receive_ints')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: ints
			integer (c_int32_t), value :: length
		end function receive_ints

		integer (c_int32_t) function receive_longs &
			(longs, length) &
			bind(c, name='forsockets_receive_longs')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: longs
			integer (c_int32_t), value :: length
		end function receive_longs

		integer (c_int32_t) function receive_floats &
			(floats, length) &
			bind(c, name='forsockets_receive_floats')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: floats
			integer (c_int32_t), value :: length
		end function receive_floats

		integer (c_int32_t) function receive_doubles &
			(doubles, length) &
			bind(c, name='forsockets_receive_doubles')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: doubles
			integer (c_int32_t), value :: length
		end function receive_doubles

		integer (c_int32_t) function receive_booleans &
			(booleans, length) &
			bind(c, name='forsockets_receive_booleans')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: booleans
			integer (c_int32_t), value :: length
		end function receive_booleans

		integer (c_int32_t) function receive_string &
			(string, length) &
			bind(c, name='forsockets_receive_string')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: string
			integer (c_int32_t), value :: length
		end function receive_string

		integer (c_int32_t) function send_ints &
			(ints, length) &
			bind(c, name='forsockets_send_ints')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: ints
			integer (c_int32_t), value :: length
		end function send_ints

		integer (c_int32_t) function send_longs &
			(longs, length) &
			bind(c, name='forsockets_send_longs')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: longs
			integer (c_int32_t), value :: length
		end function send_longs

		integer (c_int32_t) function send_floats &
			(floats, length) &
			bind(c, name='forsockets_send_floats')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: floats
			integer (c_int32_t), value :: length
		end function send_floats

		integer (c_int32_t) function send_doubles &
			(doubles, length) &
			bind(c, name='forsockets_send_doubles')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: doubles
			integer (c_int32_t), value :: length
		end function send_doubles

		integer (c_int32_t) function send_booleans &
			(booleans, length) &
			bind(c, name='forsockets_send_booleans')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: booleans
			integer (c_int32_t), value :: length
		end function send_booleans

		integer (c_int32_t) function send_string &
			(string, length) &
			bind(c, name='forsockets_send_string')
			use iso_c_binding
			implicit none
			type (c_ptr), value :: string
			integer (c_int32_t), value :: length
		end function send_string

		integer (c_int32_t) function forsockets_init &
			(port) &
			bind(c, name='forsockets_init')
			use iso_c_binding
			implicit none
			integer (c_int32_t), value :: port
		end function forsockets_init

		integer (c_int32_t) function forsockets_close &
			() &
			bind(c, name='forsockets_close')
			use iso_c_binding
			implicit none
		end function forsockets_close

	end interface
end module forsockets


subroutine run_loop
  use iso_c_binding
  use forsockets
  use iso_fortran_env

  implicit none

  integer (c_int32_t) :: ioerror, maxlen = 1
  integer :: must_run_loop
  integer i, str_len, offset, max_string_length,error, length

!  character (c_char), allocatable, target :: characters(:,:)

!  character (len=100000) :: characters
!
!  character (len=100000) :: output_characters
  character (c_char), target :: input_characters(100)
  
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

  integer (c_int32_t) success
  
 ! allocate(characters(100,100))

  allocate(integers_in( maxlen *4))
  allocate(integers_out( maxlen * 1))
  allocate(longs_in( maxlen *4))
  allocate(longs_out( maxlen * 1))
  allocate(floats_in( maxlen *4))
  allocate(floats_out( maxlen * 1))
  allocate(doubles_in( maxlen *13))
  allocate(doubles_out( maxlen * 13))
  allocate(booleans_in( maxlen * 13))
  allocate(booleans_out( maxlen * 13))
  allocate(strings_in( maxlen *2))
  allocate(strings_out( maxlen *2))
  
  must_run_loop = 1

  do while (must_run_loop .eq. 1)
    success = receive_ints(c_loc(header_in), header_size)

    length = header_in(HEADER_LENGTH)

!	print*, 'received header:'
!
!    print*, header_in
!    call flush()

!    call mpi_bcast(header, 8, mpi_integer, 0, parent,&
!      ioerror)
    
    if (length .gt. maxlen) then
      maxlen = length + 255;
      deallocate(integers_in)
      deallocate(integers_out)
      deallocate(strings_in)
      deallocate(doubles_in)
      deallocate(doubles_out)
      allocate(integers_in( maxlen *4))
      allocate(integers_out( maxlen * 1))
      allocate(strings_in( maxlen *2))
      allocate(doubles_in( maxlen *13))
      allocate(doubles_out( maxlen * 13))
    end if

    if (header_in(HEADER_INT_COUNT) .gt. 0) then
      success = receive_ints(c_loc(integers_in), header_in(HEADER_INT_COUNT))
    end if
    if (header_in(HEADER_LONG_COUNT) .gt. 0) then
      success = receive_longs(c_loc(longs_in), header_in(HEADER_LONG_COUNT))
    end if
    if (header_in(HEADER_FLOAT_COUNT) .gt. 0) then
      success = receive_floats(c_loc(floats_in), header_in(HEADER_FLOAT_COUNT))
    end if
    if (header_in(HEADER_DOUBLE_COUNT) .gt. 0) then
      success = receive_doubles(c_loc(doubles_in), header_in(HEADER_DOUBLE_COUNT))
    end if
    if (header_in(HEADER_BOOLEAN_COUNT) .gt. 0) then
      success = receive_booleans(c_loc(booleans_in), header_in(HEADER_BOOLEAN_COUNT))
    end if

    if (header_in(HEADER_STRING_COUNT) .gt. 0) then
      success = receive_ints(c_loc(strings_in), header_in(HEADER_STRING_COUNT))

!        print*, 'received string header:'
!    	print*, strings_in
!    	call flush()

	  max_string_length = 0
      do i = 1, header_in(HEADER_STRING_COUNT), 1
        if (strings_in(i) .gt. max_string_length) then
        	max_string_length = strings_in(i)
        end if
      end do

	  !space for all strings in this one call
	  !allocate(characters(header_in(HEADER_STRING_COUNT), max_string_length))

      do i = 1, header_in(HEADER_STRING_COUNT), 1
      	success = receive_string(c_loc(input_characters), strings_in(i))

!      	print*, 'received string:'
!
!    	call flush()

      end do

    end if

	header_out = 0

! 	header_out(HEADER_FLAGS) = 0
    header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
    header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
    header_out(HEADER_LENGTH) = header_in(HEADER_LENGTH)
!	header_out(HEADER_INT_COUNT) = 0
!	header_out(HEADER_LONG_COUNT) = 0
!	header_out(HEADER_FLOAT_COUNT) = 0
!	header_out(HEADER_DOUBLE_COUNT) = 0
!	header_out(HEADER_BOOLEAN_COUNT) = 0
!	header_out(HEADER_STRING_COUNT) = 0

    error = 0

	!call exit(23)

	print*, "performing function", header_in(HEADER_FUNCTION_ID)
	call flush()

    select case (header_in(HEADER_FUNCTION_ID))
      case(0)
        must_run_loop = 0
      case(670175614)

        call initialize( &
          doubles_in(1) ,&
          doubles_in(2) ,&
          doubles_in(3) ,&
          doubles_in(4) ,&
          doubles_in(5) ,&
          integers_in(1) ,&
          integers_in(2) ,&
          integers_in(3) ,&
          integers_in(4) ,&
          doubles_in(6) ,&
          doubles_in(7) ,&
          doubles_in(8) ,&
          doubles_in(9) ,&
          integers_out(1) &
        )
        header_out(HEADER_INT_COUNT) = 1

        print*, 'initialized, result ='
        print*,   integers_out(1)
      
      case(1024680297)
        do i = 1, header_in(HEADER_LENGTH), 1
          call get_time_step( &
            integers_in(i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 2 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 3 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 4 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_out(i) &
          )
        end do

        header_out(HEADER_DOUBLE_COUNT) = 1 * header_in(HEADER_LENGTH)
      
      case(1141573512)
        integers_out(1) = 0
        header_out(HEADER_INT_COUNT) = 1
      
      case(1658568341)
        do i = 1, header_in(HEADER_LENGTH), 1
          call evolve0( &
            integers_in(i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 2 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 3 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 4 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 5 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 6 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 7 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 8 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 9 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 10 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 11 * header_in(HEADER_LENGTH)) + i) ,&
            doubles_in(( 12 * header_in(HEADER_LENGTH)) + i) &
          )
          integers_out(i) = integers_in(i)
          doubles_out(i) = doubles_in(i)
          doubles_out(( 1 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 1 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 2 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 2 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 3 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 3 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 4 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 4 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 5 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 5 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 6 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 6 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 7 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 7 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 8 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 8 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 9 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 9 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 10 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 10 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 11 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 11 * header_in(HEADER_LENGTH)) + i)
          doubles_out(( 12 * header_in(HEADER_LENGTH)) + i) = doubles_in(( 12 * header_in(HEADER_LENGTH)) + i)
        end do
        header_out(HEADER_INT_COUNT) = 1 * header_in(HEADER_LENGTH)
        header_out(HEADER_DOUBLE_COUNT) = 13 * header_in(HEADER_LENGTH)
      
      case default
        error = 1
    end select
    
    print*, 'sending reply header', header_out
!    print*, 'sending doubles'
!    print*, doubles_out(1)

    call flush()

    success = send_ints(c_loc(header_out), HEADER_SIZE)

  	if (header_out(HEADER_INT_COUNT) .gt. 0) then
      success = send_ints(c_loc(integers_out), header_out(HEADER_INT_COUNT))
    end if

  	if (header_out(HEADER_LONG_COUNT) .gt. 0) then
    	success = send_longs(c_loc(longs_out), header_out(HEADER_LONG_COUNT))
  	end if

  	if (header_out(HEADER_FLOAT_COUNT) .gt. 0) then
    	success = send_floats(c_loc(floats_out), header_out(HEADER_FLOAT_COUNT))
  	end if

  	if (header_out(HEADER_DOUBLE_COUNT) .gt. 0) then
    	success = send_doubles(c_loc(doubles_out), header_out(HEADER_DOUBLE_COUNT))
  	end if

  	if (header_out(HEADER_BOOLEAN_COUNT) .gt. 0) then
    	success = send_booleans(c_loc(booleans_out), header_out(HEADER_BOOLEAN_COUNT))
  	end if

  	if (header_out(HEADER_STRING_COUNT) .gt. 0) then
  	  print*, 'error! cannot send strings'
  	  call flush()


  	end if

    if (header_in(HEADER_STRING_COUNT) .gt. 0) then
    	!deallocate(characters)
    end if


    print*, 'done sending reply'
    call flush()

    end do

    deallocate(integers_in)
    deallocate(integers_out)
    deallocate(strings_in)
    deallocate(doubles_in)
    deallocate(doubles_out)
    return
  end subroutine
  
  program muse_worker
    use iso_c_binding
  	use forsockets
    include 'mpif.h'
    integer :: provided,ioerror, port, success
    character(len=32) :: port_string

    call mpi_init_thread(mpi_thread_multiple, provided, ioerror)

    call get_command_argument(1, port_string)

    read (port_string,*) port

    success = forsockets_init(port)
    
    call run_loop()
    
    call mpi_finalize(ioerror)

    success = forsockets_close()

  end program muse_worker
