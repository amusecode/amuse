! Module for timing passages of code.
!
! Example:
!
!  use timing
!  implicit none
!  timer_state::my_timer
!
!  call timer_init(my_timer)
!  call timer_start(my_timer)
!   ... do stuff ...
!  call timer_stop(my_timer)
!  call timer_output(<file>, my_timer, "my_timer")

module timing
  implicit none

  type timer_state
     integer::count
     real(kind=8)::start
     real(kind=8)::sum
  end type timer_state

end module timing

subroutine timer_init(timer)
  use timing
  implicit none
  type(timer_state)::timer

  timer%count = 0
  timer%start = 0.0
  timer%sum = 0.0
end subroutine

subroutine timer_start(timer)
  use timing
  implicit none
  include 'mpif.h'
  type(timer_state)::timer
  integer::info

  timer%start = MPI_WTIME(info)
end subroutine

subroutine timer_stop(timer)
  use timing
  implicit none
  include 'mpif.h'
  type(timer_state)::timer
  real(kind=8)::end
  integer::info

  if (timer%start.le.0.0) return

  end = MPI_WTIME(info)
  timer%sum = timer%sum + end - timer%start
  timer%count = timer%count + 1
end subroutine

subroutine timer_inc_count(timer, delta)
  use timing
  implicit none
  include 'mpif.h'
  type(timer_state)::timer
  integer::delta

  if (timer%start.le.0.0) return
  timer%count = timer%count + delta
end subroutine
