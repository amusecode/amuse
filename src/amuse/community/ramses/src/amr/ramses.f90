program ramses
  implicit none  

  ! Read run parameters
  call read_params

  ! Start time integration
  call adaptive_loop

end program ramses

