! Helper function file_exists; placed in its own module because all FORTRAN 
!  compilers are apparently broken and can't find it if its in the twin_library 
!  module
      module file_exists_module
      contains
! Function returning TRUE or FALSE, depending on whether a file exists or not
      FUNCTION FILE_EXISTS(FILENAME)
      IMPLICIT NONE
      LOGICAL :: FILE_EXISTS
      LOGICAL :: STATUS
      CHARACTER(*) :: FILENAME
      INQUIRE(FILE=FILENAME, EXIST=STATUS)
      FILE_EXISTS = STATUS
      RETURN
      END FUNCTION
      end module


