program twin
   use real_kind
   use twinlib
  
   implicit none
   real(double), parameter :: mass = 1.0
   real(double), parameter :: age = 0.0
   integer, parameter :: max_id = 2
   integer :: status
   integer :: star_id(max_id)
   integer :: iter, i
  
   print *, 'Welcome to the TWIN library test program'
  
   ! Initialise TWIN library
   status = initialise_twin('', 20, 0.018d0, .true.)
   print *, 'TWIN initialisation finished with errnum', status

   if (status /= 0) stop

   if (have_zams_library()) then
      print *, 'ZAMS library available for this metallicity'
   else
      print *, 'No ZAMS library is available, models will be constructed from pre-MS models.'
   end if
   
  !     Load some ZAMS models
  do i=1, max_id
  !do i=max_id, 1, -1
     print *, 'loading ZAMS model with mass ', 4.-(2*i-1)*mass
     star_id(i) = new_zams_star(4.-(2*i-1)*mass, age)
     print *, 'New star with id ', star_id(i)
  end do
  print *, 'star IDs:', star_id(1:max_id)
  !
  do iter=1, 1200
     !        Evolve all stars
     do i=1, max_id
        status = evolve_one_timestep(star_id(i))
        print *, 'Star', i, 'Timestep', iter, 'completed with code', status
        if (status /= 0) stop
     end do

     !        Write output: star id, time, log T, log L, M, R
     do i=1, max_id
        write (6, '(1X,1P,I5, I3, I3, 5E15.7)') iter, star_id(i), stellar_type_of(star_id(i)), &
             age_of(star_id(i)),  log10(temperature_of(star_id(i))), log10(luminosity_of(star_id(i))),   &
             mass_of(star_id(i)), radius_of(star_id(i))
     end do
     flush(6)
  end do
  
end program twin

