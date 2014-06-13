!!     Helper function file_exists; placed in its own module because all FORTRAN 
!     compilers are apparently broken and can't find it if its in the twin_library 
!     module

module file_exists_module
contains
  
  !     Function returning TRUE or FALSE, depending on whether a file exists or not
  function file_exists(filename)
    use real_kind
    
    implicit none
    logical :: file_exists
    logical :: status
    character(*) :: filename
    
    inquire(file=filename, exist=status)
    file_exists = status
  end function file_exists
  
end module file_exists_module


