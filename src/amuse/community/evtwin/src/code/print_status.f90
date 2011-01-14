!> \file print_status.f90  Print the exit status when stopping the code

!> \brief Print a comprehensible exit message to standard output
!! 
!! The message includes the reason for termination and the the run time of the code.
!! This run time is the wall-clock time, not the CPU time 
!! (which are different if he code does not get 100% of the CPU)
subroutine print_exit_message(jo1,jo2,time0)
  implicit none
  integer, intent(in) :: jo1,jo2
  real, intent(in) :: time0
  integer :: i,jo
  real :: time
  
  do i = 1,2
     jo = jo1
     if(i.eq.2) jo = jo2
     
     if(jo.ne.-1) then
        select case(jo1)
        case(-2) 
           write(6,'(A,I2,A)', advance='no')'  Requested mesh too large for star',i,'  -'
        case(0)
           write(6,'(A,I2,A)', advance='no')'  Finished required timesteps for star',i,'  -'
        case(1)
           write(6,'(A,I2,A)', advance='no')'  Failed; backup, reduce timestep for star',i,'  -'
        case(2,12,22,32)
           write(6,'(A,I2,A)', advance='no')'  Time step reduced below limit for star',i,'  -'
        case(3)
           write(6,'(A)', advance='no')'  Star 2 evoled beyond last star-1 model -'
        case(4)
           write(6,'(A)', advance='no')'  Radius of star 1 exceeds Roche lobe by limit -'
        case(5)
           write(6,'(A,I2,A)', advance='no')'  Age greater than limit for star',i,'  -'
        case(6)
           write(6,'(A,I2,A)', advance='no')'  C-burning exceeds limit for star',i,'  -'
        case(7)
           write(6,'(A)', advance='no')'  Radius of star 2 exceeds Roche lobe by limit -'
        case(8)
           write(6,'(A,I2,A)', advance='no')'  He flash for star',i,'  -'
        case(9)
           write(6,'(A,I2,A)', advance='no')'  Massive (>1.2Mo) degenerate C/O core for star',i,'  -'
        case(10)
           write(6,'(A,I2,A)', advance='no')'  |M1dot| exceeds limit for star',i,'  -'
        case(11)
           write(6,'(A)', advance='no')'  Impermissible FDT for star 2 -'
        case(14)
           write(6,'(A,I2,A)', advance='no')'  Funny composition distribution for star',i,'  -'
        case(15)
           write(6,'(A)', advance='no')'  Terminated by hand -'
        case(16)
           write(6,'(A,I2,A)', advance='no')'  ZAHB did not converge for star',i,'  -'
        case(17)
           write(6,'(A,I2,A)', advance='no')'  Nucleosynthesis did not converge for star',i,'  -'
        case(51)
           write(6,'(A,I2,A)', advance='no')'  End of MS (core H below limit) for star',i,'  -'
        case(52)
           write(6,'(A,I2,A)', advance='no')'  Radius exceeds limit for star',i,'  -'
        case(53)
           write(6,'(A,I2,A)', advance='no')'  Convergence to target model reached minimum for star',i,'  -'
        end select
     end if
     
  end do
  
  ! The run time reported here is wall-clock time, not CPU time!
  call cpu_time(time)
  write(6,'(A,F8.1,A)')'  Evolution done, run (wall) time:', time-time0, 's.'
  
end subroutine print_exit_message

