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
      CHARACTER :: delimiter
      INQUIRE(FILE=FILENAME, EXIST=STATUS)
      IF (.NOT. STATUS) THEN
! delimiter = CALL get_environment_variable('DELIMITER',delimiter)
        delimiter ='/'
        INQUIRE(file=FILENAME//delimiter//'.',EXIST=STATUS)
      END IF
      FILE_EXISTS = STATUS
      RETURN
      END FUNCTION
      end module


