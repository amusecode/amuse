module string_manipulation

  IMPLICIT NONE
    
contains
  
!------------------------------------------------------------------------------
!
!  Program name:
!
!    convert_case
!
!  Purpose:
!
!    This program shows how to use DO and CASE control constructs. 
!    It reads in a character string, changes the case of letters and 
!    writes out new string.
!
!  Note:
!
!    Need to know the difference in the collation sequence of the upper 
!    and lower case characters - use position of upper and lower case A: 
!
!                        IACHAR('A') - IACHAR('a')
!
!------------------------------------------------------------------------------

  subroutine convert_case (string,direction)
    
    CHARACTER (LEN = *),intent(inout) :: string
    integer,intent(in) :: direction

    INTEGER :: i, upper_to_lower, len_string
    
    upper_to_lower = IACHAR("a") - IACHAR("A")
    
    ! Find length of string excluding trailing blanks
    len_string = LEN_TRIM(string)
    
    DO i = 1, len_string
       
       SELECT CASE (string(i:i))
          
          
       CASE ('A':'Z')            ! Upper case character found:
          if (direction == 0) &
               string(i:i) = ACHAR(IACHAR(string(i:i)) + upper_to_lower)
          
       CASE ('a':'z')            ! Lower case character found:
          if (direction == 1) &
               string(i:i) = ACHAR(IACHAR(string(i:i)) - upper_to_lower)
          !  No change for any other characters
       END SELECT
       
    END DO
    
  END subroutine  convert_case

end module string_manipulation
