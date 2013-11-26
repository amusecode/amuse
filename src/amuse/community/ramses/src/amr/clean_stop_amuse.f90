module amuse_error_code
    integer::error_code = 0
end module amuse_error_code


subroutine clean_stop
    use amuse_error_code
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    write(*,*) "'clean_stop' was called, but we let AMUSE handle it"
    error_code = -1
end subroutine clean_stop
