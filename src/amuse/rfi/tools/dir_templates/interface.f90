module {interface_module}


contains

  function echo_int(input, output)
      integer :: echo
      integer :: echo_int
      integer ::  input, output
      output = echo(input)
      echo_int = 0
  end function

end module

