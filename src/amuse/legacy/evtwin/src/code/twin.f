      program twin
      use twin_library
      implicit none
      double precision, parameter :: mass = 1.0
      double precision, parameter :: age = 0.0
      integer, parameter :: max_id = 4
      integer :: status
      integer :: star_id(max_id)
      integer :: iter, id

      print *, 'Welcome to the TWIN library test program'
      
!     Initialise TWIN library: files are under the current directory, Z=0.02
      status = initialise_twin('.', 20, '02');
      print *, 'TWIN initialisation finished with errnum', status
      
!     Load some ZAMS models
      do id=1, max_id
         print *, 'loading ZAMS model with mass ', (2*id-1)*mass
         star_id(id) = load_zams_star((2*id-1)*mass, age)
         call swap_out_star(star_id(id))
         print *, 'New star with id ', star_id(id)
      end do
      print *, 'star IDs:', star_id(1:max_id)
      
      do iter=1, 1200
!        Evolve all stars
         do id=1, max_id
            call select_star(star_id(id))
            status = twin_evolve()
            print *, 'Star', id, 'Timestep', iter, 'completed with code', status
            if (status /= 0) stop
         end do

!        Flush cached data
         call flush_star()

!        Write output: star id, time, log T, log L, M, R
         do id=1, max_id
            write (6, '(I5, I3, 5E14.7)') iter, star_id(id), get_age(star_id(id)),
     &         log10(get_temperature(star_id(id))), log10(get_luminosity(star_id(id))), 
     &         get_mass(star_id(id)), get_radius(star_id(id))
         end do
         call flush(6)
      end do
      
      end program
