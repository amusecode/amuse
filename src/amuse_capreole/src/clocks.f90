module clocks

  use precision, only: dp
  use file_admin, only: log_unit
  use my_mpi, only: rank

  implicit none

  ! Start and end time for CPU report
  real :: cputime1 !< Start time for CPU report
  real :: cputime2 !< End time for CPU report
  real(kind=dp) :: cpu_seconds=0.0
  integer :: cpu_hours=0
  integer :: cpu_minutes=0
  
  ! Wall clock time variables
  integer :: cntr1 !< Start time wall clock
  integer :: cntr2 !< End time wall clock
  integer :: countspersec !< counts per second (for wall clock time)
  real(kind=dp) :: clock_seconds=0.0
  integer :: clock_hours=0
  integer :: clock_minutes=0
  
contains

  !=======================================================================

  subroutine setup_clocks
    
    call setup_cpuclock()
    call setup_wallclock()
    
  end subroutine setup_clocks
  
  !=======================================================================

  subroutine setup_cpuclock
    
    ! Initialize cpu timer
    call cpu_time(cputime1)
    
  end subroutine setup_cpuclock
  
  !=======================================================================

  subroutine setup_wallclock
    
    ! Initialize wall cock timer
    call system_clock(cntr1)
    
  end subroutine setup_wallclock
  
  !=======================================================================

  subroutine update_clocks
    
    call update_cpuclock
    call update_wallclock
    
  end subroutine update_clocks
  
  !=======================================================================

  subroutine update_cpuclock
    
    ! Find out intermediate CPU time (to avoid overflowing the counter)
    call cpu_time(cputime2)
    cpu_seconds=cpu_seconds+real(cputime2-cputime1)
    cputime1=cputime2
    cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
    cpu_seconds = MOD ( cpu_seconds , 60.0d0 )
    cpu_hours = cpu_hours + cpu_minutes / 60
    cpu_minutes = MOD ( cpu_minutes , 60 )
    
  end subroutine update_cpuclock
  
  !=======================================================================

  subroutine update_wallclock
    
    call system_clock(cntr2,countspersec)
    clock_seconds=clock_seconds+real(cntr2-cntr1)/real(countspersec)
    cntr1=cntr2
    clock_minutes = clock_minutes + int(clock_seconds) / 60
    clock_seconds = MOD ( clock_seconds , 60.0d0 )
    clock_hours = clock_hours + clock_minutes / 60
    clock_minutes = MOD ( clock_minutes , 60 )
    
  end subroutine update_wallclock
  
  !=======================================================================

  subroutine report_clocks
    
    call report_cpuclock
    call report_wallclock
    
  end subroutine report_clocks
  
  !=======================================================================

  subroutine report_cpuclock
    
    call update_cpuclock ()
    if (rank == 0) then
       write(log_unit,*) "CPU time: ",cpu_hours,' hours',cpu_minutes,' minutes', &
            cpu_seconds,' seconds.'
    endif

  end subroutine report_cpuclock
  
  !=======================================================================

  subroutine report_wallclock
    
    call update_wallclock ()
    if (rank == 0) then
       write(log_unit,*) "Wall clock time: ",clock_hours,' hours', &
            clock_minutes,' minutes',clock_seconds,' seconds.'
    endif

  end subroutine report_wallclock
  
end module clocks
