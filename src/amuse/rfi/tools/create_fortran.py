from amuse.support.core import late
from amuse.support import exceptions

from amuse import config

from amuse.rfi.tools.create_code import GenerateASourcecodeString
from amuse.rfi.tools.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.rfi.tools.create_code import DTypeSpec
from amuse.rfi.tools.create_code import dtypes
from amuse.rfi.tools.create_code import DTypeToSpecDictionary
from amuse.rfi.tools import create_definition
from amuse.rfi.core import LegacyFunctionSpecification



dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('integers_in','integers_out','HEADER_INTEGER_COUNT', 'integer', 'integer'),
    'int64' : DTypeSpec('longs_in', 'longs_out', 'HEADER_LONG_COUNT', 'integer*8', 'long'),
    'float32' : DTypeSpec('floats_in', 'floats_out', 'HEADER_FLOAT_COUNT', 'real*4', 'float'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out', 'HEADER_DOUBLE_COUNT', 'real*8', 'double'),
    'bool' : DTypeSpec('booleans_in', 'booleans_out', 'HEADER_BOOLEAN_COUNT', 'logical', 'boolean'),
    'string' : DTypeSpec('strings_in', 'strings_out', 'HEADER_STRING_COUNT', 'integer*4', 'integer'),
})

CONSTANTS_STRING = """
  integer HEADER_FLAGS, HEADER_CALL_ID, HEADER_FUNCTION_ID, HEADER_CALL_COUNT, & 
        HEADER_INTEGER_COUNT, HEADER_LONG_COUNT, HEADER_FLOAT_COUNT, & 
        HEADER_DOUBLE_COUNT, HEADER_BOOLEAN_COUNT, HEADER_STRING_COUNT, & 
        HEADER_SIZE

  parameter (HEADER_FLAGS=1, HEADER_CALL_ID=2, HEADER_FUNCTION_ID=3, & 
        HEADER_CALL_COUNT=4, HEADER_INTEGER_COUNT=5, HEADER_LONG_COUNT=6, & 
        HEADER_FLOAT_COUNT=7, HEADER_DOUBLE_COUNT=8, & 
        HEADER_BOOLEAN_COUNT=9, HEADER_STRING_COUNT=10, & 
        HEADER_SIZE=10)
"""

ARRAY_DEFINES_STRING = """
  integer*4, target :: header_in(HEADER_SIZE)
  integer*4, target :: header_out(HEADER_SIZE)
  
  integer*4, allocatable, target :: integers_in(:)
  integer*4, allocatable, target :: integers_out(:)
  
  integer*8, allocatable, target :: longs_in(:)
  integer*8, allocatable, target :: longs_out(:)
  
  real*4, allocatable, target :: floats_in(:)
  real*4, allocatable, target :: floats_out(:)
  
  real*8, allocatable, target :: doubles_in(:)
  real*8, allocatable, target :: doubles_out(:)
  
  logical, allocatable, target :: booleans_in(:)
  logical, allocatable, target :: booleans_out(:)
  
  integer*4, allocatable, target :: string_sizes_in(:)
  integer*4, allocatable, target :: string_sizes_out(:)
  
  character (len=256), allocatable, target :: strings_in(:)
  character (len=256), allocatable, target :: strings_out(:)
  
  character (len=100000) :: characters_in
  character (len=100000) :: characters_out
"""

ISO_ARRAY_DEFINES_STRING = """
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
  
  logical, allocatable, target :: booleans_in(:)
  logical, allocatable, target :: booleans_out(:)
  
  integer (c_int32_t), allocatable, target :: string_sizes_in(:)
  integer (c_int32_t), allocatable, target :: string_sizes_out(:)

  character (c_char), allocatable, target :: strings_in(:) * 256
  character (c_char), allocatable, target :: strings_out(:) * 256
  
  character (kind=c_char, len=100000), target :: characters_in
  character (kind=c_char, len=100000), target :: characters_out
"""

MODULE_GLOBALS_STRING = """
  integer, save :: polling_interval = 0
"""

POLLING_FUNCTIONS_STRING = """
    FUNCTION internal__get_message_polling_interval(outval)
        INTEGER,intent(out) :: outval
        INTEGER :: internal__get_message_polling_interval
        outval = polling_interval
        internal__get_message_polling_interval = 0
    END FUNCTION
    FUNCTION internal__set_message_polling_interval(inval)
        INTEGER,intent(in) :: inval
        INTEGER :: internal__set_message_polling_interval
        polling_interval = inval
        internal__set_message_polling_interval = 0
    END FUNCTION
"""

RECV_HEADER_SLEEP_STRING = """
    SUBROUTINE mpi_recv_header(parent, ioerror)
        use iso_c_binding
        implicit none
        
        INCLUDE 'mpif.h'
      
        integer,intent(in) :: parent
        integer,intent(inout) :: ioerror
        integer :: request_status(MPI_STATUS_SIZE),header_request
        logical is_finished
        
        INTERFACE
          INTEGER (C_INT) FUNCTION usleep(useconds) bind(C)
          !SUBROUTINE usleep(useconds) bind(C)
            use iso_c_binding
            implicit none
            INTEGER(c_int32_t), value  :: useconds
          END
        END INTERFACE
        
        call MPI_Irecv(header_in, HEADER_SIZE, MPI_INTEGER, 0, 989, parent, header_request, ioerror)
        if(polling_interval.GT.0) then
            is_finished = .false.
            call MPI_Test(header_request, is_finished, request_status, ioerror)
            DO WHILE(.NOT. is_finished)
                ioerror =  usleep(int(polling_interval, c_int32_t))
                call MPI_Test(header_request, is_finished, request_status, ioerror)
            END DO
            call MPI_Wait(header_request, request_status, ioerror)
        else
            call MPI_Wait(header_request, request_status, ioerror)
        endif
    END SUBROUTINE
"""

RECV_HEADER_WAIT_STRING = """
    SUBROUTINE mpi_recv_header(parent, ioerror)
        implicit none
        INCLUDE 'mpif.h'
        integer,intent(in) :: parent
        integer,intent(inout) :: ioerror
        integer :: request_status(MPI_STATUS_SIZE),header_request
        
        call MPI_Irecv(header_in, HEADER_SIZE, MPI_INTEGER, 0, 989, parent, header_request, ioerror)
        call MPI_Wait(header_request, request_status, ioerror)
    END SUBROUTINE
"""
RUN_LOOP_MPI_STRING = """
    SUBROUTINE run_loop_mpi
      implicit none
      
      INCLUDE 'mpif.h'
      
      integer :: provided
      integer :: rank, parent, ioerror, max_call_count = 255
      integer :: must_run_loop, maximum_size, total_string_length
      integer i, offset, call_count
      
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, ioerror)
      
      ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
      ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
      ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
      ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
      ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
      ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
      ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
      ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
      ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
      ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
      
      ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
      !ensure there is at least one string to return an error code in
      ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      
      call MPI_COMM_GET_PARENT(parent, ioerror)
      call MPI_COMM_RANK(parent, rank, ioerror)
      
      must_run_loop = 1
      
      do while (must_run_loop .eq. 1)
        call mpi_recv_header(parent, ioerror)
        
        !print*, 'fortran: got header ', header_in
        
        call_count = header_in(HEADER_CALL_COUNT)
        
        IF (call_count .gt. max_call_count) THEN
          max_call_count = call_count + 255;
          DEALLOCATE(integers_in)
          DEALLOCATE(integers_out)
          DEALLOCATE(longs_in)
          DEALLOCATE(longs_out)
          DEALLOCATE(floats_in)
          DEALLOCATE(floats_out)
          DEALLOCATE(doubles_in)
          DEALLOCATE(doubles_out)
          DEALLOCATE(booleans_in)
          DEALLOCATE(booleans_out)
          DEALLOCATE(string_sizes_in)
          DEALLOCATE(string_sizes_out)
          DEALLOCATE(strings_in)
          DEALLOCATE(strings_out)
          ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
          ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
          ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
          ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
          ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
          ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
          ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
          ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
          ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
          ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
        END IF
    
        if (header_in(HEADER_INTEGER_COUNT) .gt. 0) then
          call MPI_BCast(integers_in, header_in(HEADER_INTEGER_COUNT), MPI_INTEGER, 0, parent, ioError);
        end if
        if (header_in(HEADER_LONG_COUNT) .gt. 0) then
          call MPI_BCast(longs_in, header_in(HEADER_LONG_COUNT), MPI_INTEGER8, 0, parent, ioError);
        end if
        if (header_in(HEADER_FLOAT_COUNT) .gt. 0) then
          call MPI_BCast(floats_in, header_in(HEADER_FLOAT_COUNT), MPI_REAL, 0, parent, ioError);
        end if
        if (header_in(HEADER_DOUBLE_COUNT) .gt. 0) then
          call MPI_BCast(doubles_in, header_in(HEADER_DOUBLE_COUNT), MPI_REAL8, 0, parent, ioError);
        end if
        if (header_in(HEADER_BOOLEAN_COUNT) .gt. 0) then
          call MPI_BCast(booleans_in, header_in(HEADER_BOOLEAN_COUNT), MPI_LOGICAL, 0, parent, ioError);
        end if
        if (header_in(HEADER_STRING_COUNT) .gt. 0) then
          strings_in = ' '
          call MPI_BCast(string_sizes_in, header_in(HEADER_STRING_COUNT), MPI_INTEGER, 0, parent, ioError);

          maximum_size = 0
          total_string_length = 0
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              total_string_length = total_string_length + string_sizes_in(i) + 1
              if (string_sizes_in(i) .gt. maximum_size) then
                maximum_size = string_sizes_in(i)
              end if
          end do
          
          call MPI_BCast(characters_in, total_string_length, MPI_CHARACTER, 0, parent, ioError);
          
          offset = 1
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              strings_in(i) = ' '
              strings_in(i)  = characters_in(offset : (offset + string_sizes_in(i)))
              strings_in(i)((string_sizes_in(i) + 1):(string_sizes_in(i) + 1)) = ' ' 
              offset = offset + string_sizes_in(i) + 1
              !print*, 'fortran: strings_in(i) ', i, strings_in(i) , ' of length ', string_sizes_in(i), &
              !' actually of size ', len_trim(strings_in(i))

          end do
          
        end if
        
        header_out = 0
        header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
        header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
        header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)
        
        strings_out = ' '
        
        must_run_loop = handle_call()
        
        !print*, 'fortran: sending header ', header_out
    
        if (rank .eq. 0 ) then
    
          call MPI_SEND(header_out, HEADER_SIZE, MPI_INTEGER, 0, 999, parent, ioerror);
    
          if (header_out(HEADER_INTEGER_COUNT) .gt. 0) then
            call MPI_SEND(integers_out,  header_out(HEADER_INTEGER_COUNT), MPI_INTEGER, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_LONG_COUNT) .gt. 0) then
            call MPI_SEND(longs_out,  header_out(HEADER_LONG_COUNT), MPI_INTEGER8, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_FLOAT_COUNT) .gt. 0) then
            call MPI_SEND(floats_out,  header_out(HEADER_FLOAT_COUNT), MPI_REAL, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_DOUBLE_COUNT) .gt. 0) then
            call MPI_SEND(doubles_out, header_out(HEADER_DOUBLE_COUNT), MPI_REAL8, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_BOOLEAN_COUNT) .gt. 0) then
            call MPI_SEND(booleans_out, header_out(HEADER_BOOLEAN_COUNT), MPI_LOGICAL, 0, 999, parent, ioerror)
          end if
       
          if (header_out(HEADER_STRING_COUNT) .gt. 0) then
          
            offset = 1
            do i = 1, header_out(HEADER_STRING_COUNT),1
              
              string_sizes_out(i) = len_trim(strings_out(i))
              
              !print*, 'fortran: sending strings, strings_out(i) ', i, strings_out(i) , ' of length ', string_sizes_out(i), &
              !' actually of size ', len_trim(strings_out(i))
              
              characters_out(offset:offset+string_sizes_out(i)) = strings_out(i)
              offset = offset + string_sizes_out(i) + 1
              characters_out(offset:offset) = char(0)
            end do
            
            call MPI_SEND(string_sizes_out, header_out(HEADER_STRING_COUNT), MPI_INTEGER, 0, 999, parent, ioerror)
            call MPI_SEND(characters_out, offset -1, MPI_CHARACTER, 0, 999, parent, ioerror)
          end if
        end if
      end do
    
      DEALLOCATE(integers_in)
      DEALLOCATE(integers_out)
      DEALLOCATE(longs_in)
      DEALLOCATE(longs_out)
      DEALLOCATE(floats_in)
      DEALLOCATE(floats_out)
      DEALLOCATE(doubles_in)
      DEALLOCATE(doubles_out)
      DEALLOCATE(booleans_in)
      DEALLOCATE(booleans_out)
      DEALLOCATE(string_sizes_in)
      DEALLOCATE(string_sizes_out)
      DEALLOCATE(strings_in)
      DEALLOCATE(strings_out)

      call MPI_COMM_DISCONNECT(parent, ioerror)
      call MPI_FINALIZE(ioerror)
      return
    end subroutine
"""

RUN_LOOP_SOCKETS_STRING = """
    SUBROUTINE run_loop_sockets
      use iso_c_binding
      use FortranSocketsInterface
      
      implicit none
    
      integer :: max_call_count = 255
      integer :: must_run_loop
      integer i, call_count, port
      character(len=32) :: port_string
      logical (c_bool), allocatable, target :: c_booleans_in(:)
      logical (c_bool), allocatable, target :: c_booleans_out(:)
      
      ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
      ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
      ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
      ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
      ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
      ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
      ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
      ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
      ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      
      ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      
      ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
      ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
      
      !ensure there is at least one string to return an error code in
      ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      ALLOCATE(string_sizes_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      
      call get_command_argument(1, port_string)

      read (port_string,*) port

      call forsockets_init(port)
      
      must_run_loop = 1
      
      do while (must_run_loop .eq. 1)
        call receive_integers(c_loc(header_in), HEADER_SIZE)
        
        !print*, 'fortran sockets: got header ', header_in
        
        call_count = header_in(HEADER_CALL_COUNT)
        
        IF (call_count .gt. max_call_count) THEN
          max_call_count = call_count + 255;
          DEALLOCATE(integers_in)
          DEALLOCATE(integers_out)
          DEALLOCATE(longs_in)
          DEALLOCATE(longs_out)
          DEALLOCATE(floats_in)
          DEALLOCATE(floats_out)
          DEALLOCATE(doubles_in)
          DEALLOCATE(doubles_out)
          DEALLOCATE(booleans_in)
          DEALLOCATE(booleans_out)
          DEALLOCATE(c_booleans_in)
          DEALLOCATE(c_booleans_out)
          DEALLOCATE(string_sizes_in)
          DEALLOCATE(string_sizes_out)
          DEALLOCATE(strings_in)
          DEALLOCATE(strings_out)
          ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
          ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
          ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
          ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
          ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
          ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
          ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
          ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
          ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
          ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
        END IF
    
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
          call receive_booleans(c_loc(c_booleans_in), header_in(HEADER_BOOLEAN_COUNT))
          do i = 1, header_in(HEADER_BOOLEAN_COUNT), 1
              booleans_in(i) = logical(c_booleans_in(i))
          end do
        end if
        if (header_in(HEADER_STRING_COUNT) .gt. 0) then
          strings_in = ' '
          call receive_integers(c_loc(string_sizes_in), header_in(HEADER_STRING_COUNT))

          do i = 1, header_in(HEADER_STRING_COUNT), 1
              strings_in(i) = ' '
              call receive_string(c_loc(strings_in(i)), string_sizes_in(i))
          end do
          
        end if
        
        header_out = 0
        header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
        header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
        header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)
        
        strings_out = ' '
        
        must_run_loop = handle_call()
        
        !print*, 'fortran: sending header ', header_out
    
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
          do i = 1, header_out(HEADER_BOOLEAN_COUNT), 1
              c_booleans_out(i) = logical(booleans_out(i), c_bool) 
          end do

          call send_booleans(c_loc(c_booleans_out), header_out(HEADER_BOOLEAN_COUNT))
        end if
   
        if (header_out(HEADER_STRING_COUNT) .gt. 0) then
      
          do i = 1, header_out(HEADER_STRING_COUNT),1
            string_sizes_out(i) = len_trim(strings_out(i))
            
            !print*, 'fortran: sending strings, strings_out(i) ', i, strings_out(i) , ' of length ', string_sizes_out(i), &
            !'actually of size ', len_trim(strings_out(i))
          end do
        
          call send_integers(c_loc(string_sizes_out), header_out(HEADER_STRING_COUNT))
          
          do i = 1, header_out(HEADER_STRING_COUNT),1
            call send_string(c_loc(strings_out(i)), string_sizes_out(i))
          end do
        end if
      end do
    
      DEALLOCATE(integers_in)
      DEALLOCATE(integers_out)
      DEALLOCATE(longs_in)
      DEALLOCATE(longs_out)
      DEALLOCATE(floats_in)
      DEALLOCATE(floats_out)
      DEALLOCATE(doubles_in)
      DEALLOCATE(doubles_out)
      DEALLOCATE(booleans_in)
      DEALLOCATE(booleans_out)
      DEALLOCATE(c_booleans_in)
      DEALLOCATE(c_booleans_out)
      DEALLOCATE(string_sizes_in)
      DEALLOCATE(string_sizes_out)
      DEALLOCATE(strings_in)
      DEALLOCATE(strings_out)

      call forsockets_close()
      return
    end subroutine
"""

EMPTY_RUN_LOOP_SOCKETS_STRING = """
    subroutine run_loop_sockets
    
      print*, 'fortran: sockets channel not supported in this worker'
      return
    
    end subroutine
"""

RUN_LOOP_SOCKETS_MPI_STRING = """
    SUBROUTINE run_loop_sockets_mpi
      use iso_c_binding
      use FortranSocketsInterface

      implicit none
      
      include 'mpif.h'
      
      integer :: provided
      integer :: max_call_count = 255
      integer :: must_run_loop
      integer i, call_count, port, rank, ioerror
      character(len=32) :: port_string
      logical (c_bool), allocatable, target :: c_booleans_in(:)
      logical (c_bool), allocatable, target :: c_booleans_out(:)

      
      ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
      ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
      ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
      ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
      ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
      ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
      ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
      ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
      ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      
      ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
      ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
      
      !ensure there is at least one string to return an error code in
      ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      ALLOCATE(string_sizes_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      
      call mpi_init_thread(mpi_thread_multiple, provided, ioerror)
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ioerror)

      if (rank .eq. 0) then
        call get_command_argument(1, port_string)

        read (port_string,*) port

        call forsockets_init(port)
      end if
      
      must_run_loop = 1
      
      do while (must_run_loop .eq. 1)
        if (rank .eq. 0) then
          call receive_integers(c_loc(header_in), HEADER_SIZE)
        end if
        call MPI_BCast(header_in, HEADER_SIZE , MPI_INTEGER, 0, MPI_COMM_WORLD, ioerror)
        
        !print*, 'fortran sockets mpi: got header ', header_in
        
        call_count = header_in(HEADER_CALL_COUNT)
        
        IF (call_count .gt. max_call_count) THEN
          max_call_count = call_count + 255;
          DEALLOCATE(integers_in)
          DEALLOCATE(integers_out)
          DEALLOCATE(longs_in)
          DEALLOCATE(longs_out)
          DEALLOCATE(floats_in)
          DEALLOCATE(floats_out)
          DEALLOCATE(doubles_in)
          DEALLOCATE(doubles_out)
          DEALLOCATE(booleans_in)
          DEALLOCATE(booleans_out)
          DEALLOCATE(c_booleans_in)
          DEALLOCATE(c_booleans_out)
          DEALLOCATE(string_sizes_in)
          DEALLOCATE(string_sizes_out)
          DEALLOCATE(strings_in)
          DEALLOCATE(strings_out)
          ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
          ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
          ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
          ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
          ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
          ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
          ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
          ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
          ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
          ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
        END IF
    
        if (header_in(HEADER_INTEGER_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_integers(c_loc(integers_in), header_in(HEADER_INTEGER_COUNT))
          end if
          call MPI_BCast(integers_in, header_in(HEADER_INTEGER_COUNT), MPI_INTEGER, 0, MPI_COMM_WORLD, ioError);
        end if
        
        if (header_in(HEADER_LONG_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_longs(c_loc(longs_in), header_in(HEADER_LONG_COUNT))
          end if
            call MPI_BCast(longs_in, header_in(HEADER_LONG_COUNT), MPI_INTEGER8, 0, MPI_COMM_WORLD, ioError);
        end if
        
        if (header_in(HEADER_FLOAT_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_floats(c_loc(floats_in), header_in(HEADER_FLOAT_COUNT))
          end if
          call MPI_BCast(floats_in,  header_in(HEADER_FLOAT_COUNT), MPI_REAL, 0, MPI_COMM_WORLD, ioerror)
        end if
        
        if (header_in(HEADER_DOUBLE_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_doubles(c_loc(doubles_in), header_in(HEADER_DOUBLE_COUNT))
          end if
          call MPI_BCast(doubles_in, header_in(HEADER_DOUBLE_COUNT), MPI_REAL8, 0, MPI_COMM_WORLD, ioerror)
        end if
        
        if (header_in(HEADER_BOOLEAN_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_booleans(c_loc(c_booleans_in), header_in(HEADER_BOOLEAN_COUNT))
            do i = 1, header_in(HEADER_BOOLEAN_COUNT), 1
              booleans_in(i) = logical(c_booleans_in(i))
            end do
          end if
          call MPI_BCast(booleans_in, header_in(HEADER_BOOLEAN_COUNT), MPI_LOGICAL, 0, MPI_COMM_WORLD, ioerror)
        end if
        
        if (header_in(HEADER_STRING_COUNT) .gt. 0) then

          strings_in = ' '
            
          if (rank .eq. 0) then
            call receive_integers(c_loc(string_sizes_in), header_in(HEADER_STRING_COUNT))
          end if
            
          call MPI_BCast(string_sizes_in, header_in(HEADER_STRING_COUNT), MPI_INTEGER, 0, MPI_COMM_WORLD, ioError);

          do i = 1, header_in(HEADER_STRING_COUNT), 1
            strings_in(i) = ' '
            if (rank .eq. 0) then
              call receive_string(c_loc(strings_in(i)), string_sizes_in(i))
            end if
            call MPI_BCast(strings_in(i), string_sizes_in(i), MPI_CHARACTER, 0, MPI_COMM_WORLD, ioError);
          end do
        end if
        
        header_out = 0
        header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
        header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
        header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)
        
        strings_out = ' '
        
        must_run_loop = handle_call()
        
        call MPI_Barrier(MPI_COMM_WORLD, ioerror)
        
        if (rank .eq. 0) then
        
          !print*, 'fortran: sending header ', header_out
    
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
            do i = 1, header_out(HEADER_BOOLEAN_COUNT), 1
              c_booleans_out(i) = logical(booleans_out(i), c_bool)
              !print*, 'fortran sockets mpi: sending boolean', booleans_out(i) , i, ' send as ', c_booleans_out(i) 
            end do
        
            call send_booleans(c_loc(c_booleans_out), header_out(HEADER_BOOLEAN_COUNT))
          end if
   
          if (header_out(HEADER_STRING_COUNT) .gt. 0) then
      
            do i = 1, header_out(HEADER_STRING_COUNT),1
              string_sizes_out(i) = len_trim(strings_out(i))
            
              !print*, 'fortran: sending strings, strings_out(i) ', i, strings_out(i) , ' of length ', string_sizes_out(i), &
              !'actually of size ', len_trim(strings_out(i))
            end do
        
            call send_integers(c_loc(string_sizes_out), header_out(HEADER_STRING_COUNT))
          
            do i = 1, header_out(HEADER_STRING_COUNT),1
              call send_string(c_loc(strings_out(i)), string_sizes_out(i))
            end do
          end if
        end if
      end do
    
      DEALLOCATE(integers_in)
      DEALLOCATE(integers_out)
      DEALLOCATE(longs_in)
      DEALLOCATE(longs_out)
      DEALLOCATE(floats_in)
      DEALLOCATE(floats_out)
      DEALLOCATE(doubles_in)
      DEALLOCATE(doubles_out)
      DEALLOCATE(booleans_in)
      DEALLOCATE(booleans_out)
      DEALLOCATE(string_sizes_in)
      DEALLOCATE(string_sizes_out)
      DEALLOCATE(strings_in)
      DEALLOCATE(strings_out)

      if (rank .eq. 0) then
        call forsockets_close()
      end if
      
      call MPI_FINALIZE(ioerror)
      return
    end subroutine
"""

EMPTY_RUN_LOOP_SOCKETS_MPI_STRING = """
    subroutine run_loop_sockets_mpi
    
      print*, 'fortran: sockets channel not supported in this worker'
      return
    
    end subroutine
"""

MAIN_STRING = """
  integer :: count
  logical :: use_mpi
  character(len=32) :: use_mpi_string
 
  count = command_argument_count()
  
  use_mpi = NEEDS_MPI

  if (count .eq. 0) then
    call run_loop_mpi()
  else
    if (count .eq. 2) then
      call get_command_argument(2, use_mpi_string)
      
      if (use_mpi_string .eq. 'no_mpi') then
        use_mpi = .false.
      else if (use_mpi_string .eq. 'mpi') then
        use_mpi = .true.
      else
        use_mpi = NEEDS_MPI
      end if
    end if
  
    if (use_mpi) then
      call run_loop_sockets_mpi()
    else 
      call run_loop_sockets()
    end if
  end if
"""
        
class GenerateAFortranStringOfAFunctionSpecification(GenerateASourcecodeString):
    MAX_STRING_LEN = 256
    
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
    
    @late
    def underscore_functions_from_specification_classes(self):
        return []
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
        
    def index_string(self, index, must_copy_in_to_out = False):
        if self.specification.must_handle_array and not must_copy_in_to_out:
            if index == 0:
                return '1'
            else:
                return '( %d * call_count) + 1' % (index )
        elif self.specification.can_handle_array or (self.specification.must_handle_array and must_copy_in_to_out):
            if index == 0:
                return 'i'
            else:
                if index == -1:
                    return "i - 1"
                else:
                    return '( %d * call_count) + i' % index
        else:
            return index + 1
            
    def start(self):        
        self.specification.prepare_output_parameters()
         
        self.output_casestmt_start()
        self.out.indent()
        
        #self.output_lines_before_with_clear_out_variables()
        #self.output_lines_before_with_clear_input_variables()
        
        if self.specification.must_handle_array:
            pass
        elif self.specification.can_handle_array:
            self.out.lf() + 'do i = 1, call_count, 1'
            self.out.indent()
        
        #self.output_lines_before_with_inout_variables()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_lines_with_inout_variables()
        
        
        if self.specification.must_handle_array:
            if not self.specification.result_type is None:
                spec = self.dtype_to_spec[self.specification.result_type]
                self.out.lf() + 'DO i = 2, call_count'
                self.out.indent()
                self.out.lf() + spec.output_var_name + '(i)' + ' = ' + spec.output_var_name + '(1)'
                self.out.dedent()
                self.out.lf() + 'END DO'
        elif self.specification.can_handle_array:
            self.out.dedent()
            self.out.lf() + 'end do'
            
        self.output_lines_with_number_of_outputs()
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
        
    def output_function_parameters(self):
        self.out.indent()
        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
                self.out + ' &'
            else:
                self.out + ' ,&'
                
            if parameter.direction == LegacyFunctionSpecification.IN:
#                if parameter.datatype == 'string':
#                    self.out.n() + 'input_characters('
#                    self.out  + '( (' + self.index_string(parameter.input_index) + ')* ' + self.MAX_STRING_LEN + ')'
#                    self.out  + ':' + '(((' + self.index_string(parameter.input_index) + ')* ' + self.MAX_STRING_LEN + ') +'
#                    self.out  +  '(' + spec.input_var_name + '(' + self.index_string(parameter.input_index) + ')' + '-' 
#                    self.out  + 'get_offset(' + self.index_string(parameter.input_index) + ' - 1 , '+spec.input_var_name +') ))'
#                    self.out  + ')'
#                else:
                if parameter.datatype == 'string':
                    self.out.n() + 'strings_in(' + self.index_string(parameter.input_index) + ')'
                else:
                    self.out.n() + spec.input_var_name 
                    self.out + '(' + self.index_string(parameter.input_index) + ')'
            if parameter.direction == LegacyFunctionSpecification.INOUT:
#                if parameter.datatype == 'string':
#                    self.out.n() + 'output_characters('
#                    self.out  + '((' + self.index_string(parameter.output_index) + ')* ' + self.MAX_STRING_LEN + ')'
#                    self.out  + ':' + '(((' + self.index_string(parameter.output_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
#                    self.out  + ')'
#                else:
#                if parameter.datatype == 'string':
#                    self.out.n() + spec.input_var_name 
#                    self.out + '(' + self.index_string(parameter.input_index) + ', :)'
#                else:
                    self.out.n() + spec.input_var_name 
                    self.out + '(' + self.index_string(parameter.input_index) + ')'
            elif parameter.direction == LegacyFunctionSpecification.OUT:
#                if parameter.datatype == 'string':
#                    self.out.n() + 'output_characters('
#                    self.out  + '((' + self.index_string(parameter.output_index) + ')* ' + self.MAX_STRING_LEN + ')'
#                    self.out  + ':' + '(((' + self.index_string(parameter.output_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
#                    self.out  + ')'
#                else:
#                if parameter.datatype == 'string':
#                    self.out.n() + spec.output_var_name
#                    self.out + '(' + self.index_string(parameter.output_index) + ')(1:50)'
#                else:
                    self.out.n() + spec.output_var_name
                    self.out + '(' + self.index_string(parameter.output_index) + ')'
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out.n() + 'call_count'
                
        self.out.dedent()
        
    def output_lines_with_inout_variables(self):
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'DO i = 1, call_count'
                    self.out.indent() 
                    
                self.out.n() + spec.output_var_name 
                self.out + '(' + self.index_string(parameter.output_index, must_copy_in_to_out = True)  + ')' 
                self.out + ' = ' 
                self.out + spec.input_var_name + '(' + self.index_string(parameter.input_index, must_copy_in_to_out = True) + ')'
        
                if self.specification.must_handle_array:
                    self.out.dedent() 
                    self.out.lf() + 'END DO'
    
    def output_lines_before_with_clear_out_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.is_output():
                if parameter.datatype == 'string': 
                    self.out.lf() + 'output_characters = "x"'  
                    return
     
    def output_lines_before_with_clear_input_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.is_input():
                if parameter.datatype == 'string': 
                    self.out.lf() + 'input_characters = "x"'  
                    return
     
                
                    
    def output_lines_before_with_inout_variables(self):
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            
            if parameter.direction == LegacyFunctionSpecification.IN:
                if parameter.datatype == 'string':
                    self.out.n() + 'input_characters('
                    self.out  + '( (' + self.index_string(parameter.input_index) + ')* ' + self.MAX_STRING_LEN + ')'
                    self.out  + ':' + '(((' + self.index_string(parameter.input_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
                    self.out  + ') = &'
                    self.out.lf()
                    self.out + 'characters('
                    self.out + 'get_offset(' + self.index_string(parameter.input_index) + ' - 1 , '+spec.input_var_name +')'
                    self.out  + ':' + spec.input_var_name + '(' + self.index_string(parameter.input_index) + ')'
                    self.out  + ')' 
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string':
                    self.out.n() + 'output_characters('
                    self.out  + '( (' + self.index_string(parameter.output_index) + ')* ' + self.MAX_STRING_LEN + ')'
                    self.out  + ':' + '(((' + self.index_string(parameter.output_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
                    self.out  + ') = &'
                    self.out.lf()
                    self.out + 'characters('
                    self.out + 'get_offset(' + self.index_string(parameter.input_index) + ' - 1 , '+spec.input_var_name +')'
                    self.out  + ':' + spec.input_var_name + '(' + self.index_string(parameter.input_index) + ')'
                    self.out  + ')' 
                    
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for parameter in self.specification.output_parameters:
            count = dtype_to_count.get(parameter.datatype, 0)
            dtype_to_count[parameter.datatype] = count + 1
                
        if not self.specification.result_type is None:
            count = dtype_to_count.get(self.specification.result_type, 0)
            dtype_to_count[self.specification.result_type] = count + 1
            
        for dtype in dtype_to_count:       
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_count[dtype]
            self.out.n() + 'header_out(' + spec.counter_name + ') = ' + count + ' * call_count'
            pass
            
    def output_function_end(self):
        self.out + ' &'
        self.out.n() + ')'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
#            if self.specification.result_type == 'string':
#                self.out + 'output_characters('
#                self.out  + '( (' + self.index_string(0) + ')* ' + self.MAX_STRING_LEN + ')'
#                self.out  + ':' + '(((' + self.index_string(0) + ')+1)*' + self.MAX_STRING_LEN + '-1)'
#                self.out  + ') = &'
#                self.out.lf()
#            else:
            self.out + spec.output_var_name
            self.out + '(' + self.index_string(0) + ')' + ' = '
        else:    
            self.out + 'CALL ' 
        self.out +  self.specification.name
        if self.must_add_underscore_to_function(self.specification):
            self.out + '_'
        self.out + '('
        
    def output_casestmt_start(self):
        self.out + 'CASE(' + self.specification.id + ')'
        
    def output_casestmt_end(self):
        self.out.n() 
        
    def must_add_underscore_to_function(self, x):
           
        for cls in self.underscore_functions_from_specification_classes:
            if hasattr(cls, x.name):
                return True
        
        return False
        
        
class GenerateAFortranSourcecodeStringFromASpecificationClass(GenerateASourcecodeStringFromASpecificationClass):
    MAX_STRING_LEN = 256

    @late
    def dtype_to_spec(self):
        return dtype_to_spec 
   
    @late
    def number_of_types(self):
        return len(self.dtype_to_spec)
        
    @late
    def length_of_the_header(self):
        return 2 + self.number_of_types
        
    @late
    def underscore_functions_from_specification_classes(self):
        return []
        
    def output_sourcecode_for_function(self):
        result = GenerateAFortranStringOfAFunctionSpecification()
        result.underscore_functions_from_specification_classes = self.underscore_functions_from_specification_classes
        return result
    
    def output_needs_mpi(self):
        self.out.lf() + 'logical NEEDS_MPI'
        
        if hasattr(self, 'needs_mpi'):
            if self.needs_mpi:
                self.out.lf() + 'parameter (NEEDS_MPI=.true.)'
            else:
                self.out.lf() + 'parameter (NEEDS_MPI=.false.)'
        else:
            self.out.lf() + 'parameter (NEEDS_MPI=.true.)'
                
        self.out.lf().lf()
   
    def start(self):
        self.use_iso_c_bindings = config.compilers.fc_iso_c_bindings

        self.out + 'program amuse_worker_program'
        self.out.indent()
        
        self.output_modules()
        
        if self.use_iso_c_bindings:    
            self.out.n() + 'use iso_c_binding'
        
        self.out.n() + 'implicit none'

        self.out.n() + CONSTANTS_STRING
        
        self.output_needs_mpi()

        self.output_maximum_constants()

        self.out.lf().lf() + MODULE_GLOBALS_STRING
            
        if self.use_iso_c_bindings:
            self.out.n() + ISO_ARRAY_DEFINES_STRING
        else:
            self.out.n() + ARRAY_DEFINES_STRING
            
        
        self.out.lf().lf() + MAIN_STRING
        
        self.out.lf().lf() + 'CONTAINS'
        
        self.out + POLLING_FUNCTIONS_STRING
            
        if self.must_generate_mpi:
            
            if self.use_iso_c_bindings:    
                self.out + RECV_HEADER_SLEEP_STRING
            else:
                self.out + RECV_HEADER_WAIT_STRING
                
            self.out + RUN_LOOP_MPI_STRING

        if self.use_iso_c_bindings:
            self.out.n() + RUN_LOOP_SOCKETS_STRING
            self.out.n() + RUN_LOOP_SOCKETS_MPI_STRING
        else:
            self.out.n() + EMPTY_RUN_LOOP_SOCKETS_STRING
            self.out.n() + EMPTY_RUN_LOOP_SOCKETS_MPI_STRING
        
        self.output_handle_call()

        self.out.dedent()
        self.out.n() + 'end program amuse_worker_program'

        self._result = self.out.string

    def output_mpi_include(self):
        self.out.n() + "INCLUDE 'mpif.h'"
        
    def output_modules(self):
        self.out.n()
        if hasattr(self.specification_class, 'use_modules'):
            for x in self.specification_class.use_modules:
                self.out.n() + 'use ' + x 
                
    def must_include_declaration_of_function(self, x):
        if x.specification.name.startswith("internal__"):
            return False
        
        return True
        
        
    def output_declarations_for_the_functions(self):
        if not hasattr(self.specification_class, 'use_modules'):
            for x in self.interface_functions:
                if not self.must_include_declaration_of_function(x):
                    continue
                    
                specification = x.specification
                if specification.id == 0:
                    continue
                if specification.result_type is None:
                    continue
                if specification.result_type == 'string':
                    type = 'CHARACTER(len=255)'
                else:
                    spec = self.dtype_to_spec[specification.result_type]
                    type = spec.type
                self.out.lf() +  type + ' :: ' + specification.name
                
                if self.must_add_underscore_to_function(x):
                    self.out + '_'
        
    
    def must_add_underscore_to_function(self, x):
           
        for cls in self.underscore_functions_from_specification_classes:
            if hasattr(cls, x.specification.name):
                return True
        
        return False
         
    def output_handle_call(self):
        
        self.out.lf() + 'integer function handle_call()'
        self.out.indent().n()
        self.out.lf() + 'implicit none'

        
        self.output_declarations_for_the_functions()
        
        self.out.lf() + 'integer i, call_count'
        self.out.lf() + 'call_count = header_in(HEADER_CALL_COUNT)'
        self.out.lf() + 'handle_call = 1'
        self.out.lf() + 'SELECT CASE (header_in(HEADER_FUNCTION_ID))'
        self.out.indent().n()
        self.out.lf() + 'CASE(0)'
        self.out.indent().lf()+'handle_call = 0'
        self.out.dedent()
        
        self.output_sourcecode_for_functions()

        self.out.lf() + 'CASE DEFAULT'
        self.out.indent()
        self.out.lf() + 'header_out(HEADER_STRING_COUNT) = 1'
        self.out.lf() + 'header_out(HEADER_FLAGS) = IOR(header_out(HEADER_FLAGS), 256) '
        self.out.lf() + "strings_out(1) = 'error, illegal function id'" 
        self.out.dedent()
        
        self.out.dedent().n() + 'END SELECT'

        self.out.n() + 'return'
        self.out.dedent()
        self.out.n() + 'end function'
        
    def output_maximum_constants(self):
                
        self.out.lf() + 'integer MAX_INTEGERS_IN, MAX_INTEGERS_OUT, MAX_LONGS_IN, MAX_LONGS_OUT, &'
        self.out.lf() + 'MAX_FLOATS_IN, MAX_FLOATS_OUT, MAX_DOUBLES_IN,MAX_DOUBLES_OUT, &'
        self.out.lf() + 'MAX_BOOLEANS_IN,MAX_BOOLEANS_OUT, MAX_STRINGS_IN, MAX_STRINGS_OUT'
        self.out.lf()

        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            maximum = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype,0)

            self.out.n() + 'parameter (MAX_' + dtype_spec.input_var_name.upper() + '=' + maximum + ')'
            
            maximum =self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
            
            self.out.n() + 'parameter (MAX_' + dtype_spec.output_var_name.upper() + '=' + maximum + ')'
    

class GenerateAFortranStubStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def dtype_to_spec(self):
        return dtype_to_spec
  
    @late
    def ignore_functions_from_specification_classes(self):
        return []
        
    @late
    def underscore_functions_from_specification_classes(self):
        return []
        
    def output_sourcecode_for_function(self):
        result = create_definition.CreateFortranStub()
        result.output_definition_only = False
        return result
        
    def start(self):  
        
        self.output_modules()
        
        self.out.lf()
        
        self.output_sourcecode_for_functions()
        
        self.out.lf()
        
        self._result = self.out.string
        
    
    def must_include_interface_function_in_output(self, x):
        if x.specification.name.startswith("internal__"):
            return False
            
        for cls in self.ignore_functions_from_specification_classes:
            if hasattr(cls, x.specification.name):
                return False
        
        return True
        
    def output_modules(self):
        self.out.n()
        if hasattr(self.specification_class, 'use_modules'):
            for x in self.specification_class.use_modules:
                self.out.n() + 'use ' + x 
        
    
        
        

        
       
    
        
        
        
