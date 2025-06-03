function echo_int(int_in, int_out)
    implicit none
    integer :: int_in, int_out
    integer :: echo_int
    
    int_out = int_in
    
    echo_int = 0
end function

function echo_double(double_in, double_out)
    implicit none
    DOUBLE PRECISION :: double_in, double_out
    integer :: echo_double
    
    double_out = double_in
    
    echo_double = 0
end function

function echo_float(float_in, float_out)
    implicit none
    REAL(kind=4) :: float_in, float_out
    integer :: echo_float
    
    float_out = float_in
    
    echo_float = 0
end function

function echo_string(string_in, string_out)
    implicit none
    character(len=*) :: string_in, string_out
    integer :: echo_string
    
    string_out = string_in
    
    echo_string = 0
end function

function echo_strings(string_inout1, string_inout2)
    implicit none
    character(len=*) :: string_inout1, string_inout2
    integer :: echo_strings
    
    string_inout1(1:1) = 'A'
    string_inout2(1:1) = 'B'
    
    
    echo_strings = 0
end function


function return_string(string_in)
    implicit none
    character(len=*) :: string_in, return_string
    
    return_string = string_in
end function


function hello_string(string_out)
    implicit none
    character(len=*) :: string_out
    integer :: hello_string
    
    string_out = 'hello'
    
    hello_string = 0
end function

function print_string(string_in)
    implicit none
    character(len=*) :: string_in
    integer :: print_string
    
    write (*,*) string_in
    
    print_string = 0
end function

function print_error_string(string_in)
    implicit none
    character(len=*) :: string_in
    integer :: print_error_string
    
    write (0,*) string_in
    
    print_error_string = 0
end function

function echo_string_fixed_len(string_in, string_out)
    implicit none
    character(len=30) :: string_in, string_out
    integer :: echo_string_fixed_len
    
    string_out = string_in
    
    echo_string_fixed_len = 0
end function

function echo_array_with_result(int_in, int_out, N)
    implicit none
    integer, intent(in) :: N
    integer :: int_in(N), int_out(N)
    integer :: echo_array_with_result,  i
    
    do i = 1, N
     int_out(i) = int_in(i)
    end do
    
    echo_array_with_result = -1
end function



function echo_inout_array_with_result(inout, N) 
    implicit none
    integer, intent(in) :: N
    integer :: inout(N)
    integer :: echo_inout_array_with_result,  i
    
    do i = 1, N
     inout(i) = inout(i) + 10
    end do
    
    echo_inout_array_with_result = 11;
end function


function echo_logical(input, output)
    implicit none
    logical :: input, output
    integer :: echo_logical
    
    output = input
    print *, "INPUT=", input
    
    echo_logical = 0
end function

function echo_logical2(input, output, n)
    implicit none
    logical :: input(n), output(n)
    integer :: echo_logical2, n,i
    
    output(1:n)=.FALSE.
    do i=1,n
      if(input(i)) then
        output(i) = .TRUE.
      endif
    enddo
    
    echo_logical2 = 0
end function


function get_element_status(ind,x) result(ret)
  integer :: ind,ret
  character(len=*) :: x
  x="dry"
  ret=0
end function

