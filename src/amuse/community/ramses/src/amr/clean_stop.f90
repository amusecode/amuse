subroutine clean_stop
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::info
#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif
  stop
end subroutine clean_stop
