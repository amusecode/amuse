! ***********************************************************************
!
!   Copyright (C) 2008  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module create_zams
      
      use const_def
      use alert_lib
      use star_def
      use star_lib
      use net_def
      use chem_def
      use run_star_support

      implicit none


         
      ! controls
      character (len=256) :: zams_name
      double precision :: create_z, mlo, mhi, dmass
      namelist /create_zams_job/ &
         zams_name, create_z, mlo, mhi, dmass
      

      contains


      
      subroutine AMUSE_do_create_zams( &
            metallicity_in, zams_outfile_in, AMUSE_inlist_path, &
            dmass_in, mlo_in, mhi_in, ierr)
         use mtx_lib, only: lapack_decsol
         use num_def, only: square_matrix_type
         use num_lib
         use utils_lib, only:alloc_iounit, free_iounit
!         use star_lib, only: alloc_star, star_setup
!         use star_private_def, only: star_info, get_star_ptr
         type (star_info), pointer :: s
         character (len=*) :: zams_outfile_in, AMUSE_inlist_path
         double precision, intent(in) :: metallicity_in
         double precision, intent(in) :: dmass_in, mlo_in, mhi_in
         integer, intent(out) :: ierr
         
         integer :: io_ms_mod, io_ms_index      
         double precision :: init_m
         integer :: i, j, k, n, id, result, result_reason
!         integer :: id
         character (len=256) :: ms_file
         character (len=1024) :: line
         logical :: first_try, okay

         1 format(a40, 1pe26.16)
         2 format(a40, i6, 1pe26.16)
         3 format(a15, 2x, f15.6)
         14 format(a40, e24.14)
         
         ierr = 0

! Don't read the ZAMS controls from file...
!         call read_zams_controls(s, zams_inlist, ierr)
!         if (failed('read_zams_controls')) return
! ... but set them here:
         create_z = metallicity_in ! 2d-2
         zams_name = zams_outfile_in ! 'z2m2'
         dmass = dmass_in
         mlo = mlo_in
         mhi = mhi_in

         okay = .true.

         io_ms_mod = alloc_iounit(ierr)
         if (failed('alloc_iounit')) return
         
         io_ms_index = alloc_iounit(ierr)
         if (failed('alloc_iounit')) return

         do j=1, 10
            write(*, *)
         end do

         ms_file = trim(zams_name) // '.data'
         write(*, *) 'creating ' // trim(ms_file)
         open(unit=io_ms_index, file=trim(ms_file), action='write', status='replace')

         ms_file = trim(zams_name) // '_mod.data'
         open(unit=io_ms_mod, file=trim(ms_file), action='write', status='replace')
         n = (mhi-mlo)/dmass + 1
         
         write(*,1) 'mlo', mlo
         write(*,1) 'mhi', mhi
         
         mass_loop: do i=1, n
         
            id = alloc_star(ierr)
            if (failed('alloc_star')) return
            s => star_handles(id)
!            call get_star_ptr(id, s, ierr)
!            if (failed('get_star_ptr')) return
            call star_setup(id, AMUSE_inlist_path, ierr)
            if (failed('star_setup')) return
      
            init_m = 10**(mlo+(i-1)*dmass)
            
            if (init_m > 1) s% mesh_delta_coeff = 0.3
            if (init_m > 80) s% mesh_delta_coeff = 0.2
            
            if (init_m > 10**mhi) exit
            do j=1, 10
               write(*, *)
            end do

            s% initial_z = create_z
            s% initial_mass = init_m
            write(*, 14) 'do ' // trim(zams_name), s% initial_mass

            if (i==1) call write_index_head
            
            call star_create_pre_ms_model( &
               id, pre_ms_T_c, pre_ms_guess_rho_c, pre_ms_d_log10_P, ierr)
            if (failed('star_create_pre_ms_model')) exit
            
            call evolve_to_zams(s, id, ierr)
            if (failed('evolve_to_zams')) exit 

            call write_model(id, create_z, io_ms_mod, io_ms_index, ierr)             
            if (failed('write_model')) exit  

         end do mass_loop
         
         11 format(3x, f15.8, i15)
         write(io_ms_index, 11) -1d0, -1 ! marks end of index
         write(io_ms_index, *) ! blank line at end of index
         
         close ( io_ms_mod )
         open(unit=io_ms_mod, file=trim(ms_file), action='read', status='old', iostat=ierr)
         if (failed('open mods to read')) return
         
         do 
            read(io_ms_mod, fmt='(a)', iostat=ierr) line
            if (ierr /= 0) then
               ierr = 0; exit
            end if
            write(io_ms_index, fmt='(a)') trim(line)
         end do
         
         close ( io_ms_mod )
         close ( io_ms_index )
         
         call free_iounit(io_ms_mod)
         call free_iounit(io_ms_index)
         
         call free_star(id, ierr)
         if (failed('free_star')) return 
         
         call star_shutdown
         
         write(*, *)
         if (okay) then
            write(*, '(a)') 'finished create main sequence'
         else
            write(*, '(a)') 'failed during attempt to create main sequence'
         end if
         write(*, *)
         
         contains
         
         subroutine write_index_head
            use chem_def
            integer :: i, time_vals(8)
            character (len=10) :: date_str, time_str, zone_str
            type (star_info), pointer :: s
            1 format(a32, 2x, 1pe26.14)
            2 format(a32, 2x, i9)
            3 format(a32, 3x, a8)
            4 format(a32, 3x, a)
            call star_ptr(id, s, ierr)
            if (ierr /= 0) then
               write(*, *) 'write_model: star_ptr failed'
               return
            end if
            write(io_ms_index, '(a,/)') '          1 -- mesa/star zams'
            call date_and_time(date_str, time_str, zone_str, time_vals)
            ! write property list
            write(io_ms_index, 3) 'year_month_day_when_created', date_str(1:8)
            write(io_ms_index, 4) 'net_name', "'basic.net'"
            write(io_ms_index, 2) 'species', num_isos_for_Basic
            write(io_ms_index, 1) 'initial_z', create_z
            write(io_ms_index, *) ! blank line for end of property list
            write(io_ms_index, '(a)') '          M/Msun           n_shells'
         end subroutine write_index_head
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *)
               write(*, *) trim(str) // ' ierr', ierr
!               write(*, '(a)') trim(alert_message)
               okay = .false.
               !stop 1
            end if
         end function failed
         
         subroutine dump_initial_model 
            character (len=256) :: fname
            fname = 'initial_model.data'
            write(*, *) 'dump initial model to ' // trim(fname)
            call write_internals(s% id, fname, ierr)
            
            stop 'dump_initial_model'
            
         end subroutine dump_initial_model

      end subroutine AMUSE_do_create_zams
      
      
      
      subroutine do_create_zams( &
            s, zams_inlist, log_columns_file_in, profile_columns_file_in, ierr)
         use mtx_lib, only: lapack_decsol
         use num_def, only: square_matrix_type
         use num_lib
         use utils_lib, only:alloc_iounit, free_iounit
         type (star_info), pointer :: s
         character (len=*) :: zams_inlist, log_columns_file_in, profile_columns_file_in
         integer, intent(out) :: ierr
         
         integer :: io_ms_mod, io_ms_index      
         double precision :: init_m
         integer :: i, j, k, n, id, result, result_reason
         character (len=256) :: ms_file
         character (len=1024) :: line
         logical :: first_try, okay

         1 format(a40, 1pe26.16)
         2 format(a40, i6, 1pe26.16)
         3 format(a15, 2x, f15.6)
         14 format(a40, e24.14)
         
         ierr = 0
         id = s% id

         call read_zams_controls(s, zams_inlist, ierr)
         if (failed('read_zams_controls')) return

         okay = .true.

         if (.false.) then ! DEBUGGING
            s% log_cnt = 5
            s% photostep = 50
            s% profile_interval = 50
            s% use_graboske_et_al_screening = .true.
            s% T_mix_limit = 1d4

            s% mesh_delta_coeff = 0.5

            dmass = 0.1d0
            mlo = log10(0.42600914858318d+02)
            mhi = log10(61d0)
         
            s% de_Jager_wind_eta = 0
            s% Reimers_wind_eta = 0
            s% Blocker_wind_eta = 0
         
            s% net_name = 'basic.net'
            s% species = num_isos_for_Basic
            s% v_flag = .false.

         end if

         ierr = 0
         io_ms_mod = alloc_iounit(ierr)
         if (failed('alloc_iounit')) return
         
         io_ms_index = alloc_iounit(ierr)
         if (failed('alloc_iounit')) return

         do j=1, 10
            write(*, *)
         end do

         ms_file = trim(zams_name) // '.data'
         write(*, *) 'creating ' // trim(ms_file)
         open(unit=io_ms_index, file=trim(ms_file), action='write', status='replace')

         ms_file = trim(zams_name) // '_mod.data'
         open(unit=io_ms_mod, file=trim(ms_file), action='write', status='replace')
         n = (mhi-mlo)/dmass + 1
         
         write(*,1) 'mlo', mlo
         write(*,1) 'mhi', mhi
         
         mass_loop: do i=1, n
         
            init_m = 10**(mlo+(i-1)*dmass)
            
            if (init_m > 1) s% mesh_delta_coeff = 0.3
            if (init_m > 80) s% mesh_delta_coeff = 0.2
            
            if (init_m > 10**mhi) exit
            do j=1, 10
               write(*, *)
            end do

            s% initial_z = create_z
            s% initial_mass = init_m
            write(*, 14) 'do ' // trim(zams_name), s% initial_mass

            if (i==1) call write_index_head
            
            call star_create_pre_ms_model( &
               id, pre_ms_T_c, pre_ms_guess_rho_c, pre_ms_d_log10_P, ierr)
            if (failed('star_create_pre_ms_model')) exit
            
            call evolve_to_zams(s, id, ierr)
            if (failed('evolve_to_zams')) exit 

            call write_model(id, create_z, io_ms_mod, io_ms_index, ierr)             
            if (failed('write_model')) exit  

         end do mass_loop
         
         11 format(3x, f15.8, i15)
         write(io_ms_index, 11) -1d0, -1 ! marks end of index
         write(io_ms_index, *) ! blank line at end of index
         
         close ( io_ms_mod )
         open(unit=io_ms_mod, file=trim(ms_file), action='read', status='old', iostat=ierr)
         if (failed('open mods to read')) return
         
         do 
            read(io_ms_mod, fmt='(a)', iostat=ierr) line
            if (ierr /= 0) then
               ierr = 0; exit
            end if
            write(io_ms_index, fmt='(a)') trim(line)
         end do
         
         close ( io_ms_mod )
         close ( io_ms_index )
         
         call free_iounit(io_ms_mod)
         call free_iounit(io_ms_index)
         
         call free_star(id, ierr)
         if (failed('free_star')) return 
         
         call star_shutdown
         
         write(*, *)
         if (okay) then
            write(*, '(a)') 'finished create main sequence'
         else
            write(*, '(a)') 'failed during attempt to create main sequence'
         end if
         write(*, *)
         
         contains
         
         subroutine write_index_head
            use chem_def
            integer :: i, time_vals(8)
            character (len=10) :: date_str, time_str, zone_str
            type (star_info), pointer :: s
            1 format(a32, 2x, 1pe26.14)
            2 format(a32, 2x, i9)
            3 format(a32, 3x, a8)
            4 format(a32, 3x, a)
            call star_ptr(id, s, ierr)
            if (ierr /= 0) then
               write(*, *) 'write_model: star_ptr failed'
               return
            end if
            write(io_ms_index, '(a,/)') '          1 -- mesa/star zams'
            call date_and_time(date_str, time_str, zone_str, time_vals)
            ! write property list
            write(io_ms_index, 3) 'year_month_day_when_created', date_str(1:8)
            write(io_ms_index, 4) 'net_name', "'basic.net'"
            write(io_ms_index, 2) 'species', num_isos_for_Basic
            write(io_ms_index, 1) 'initial_z', create_z
            write(io_ms_index, *) ! blank line for end of property list
            write(io_ms_index, '(a)') '          M/Msun           n_shells'
         end subroutine write_index_head
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *)
               write(*, *) trim(str) // ' ierr', ierr
               write(*, '(a)') trim(alert_message)
               okay = .false.
               !stop 1
            end if
         end function failed
         
         subroutine dump_initial_model 
            character (len=256) :: fname
            fname = 'initial_model.data'
            write(*, *) 'dump initial model to ' // trim(fname)
            call write_internals(s% id, fname, ierr)
            
            stop 'dump_initial_model'
            
         end subroutine dump_initial_model

      end subroutine do_create_zams

      
      subroutine write_model(id, create_z, io_ms_mod, io_ms_index, ierr)
         use chem_def
         double precision, intent(in) :: create_z
         integer, intent(in) :: id, io_ms_mod, io_ms_index
         integer, intent(out) :: ierr
         
         integer :: k, j, species, nz
         type (star_info), pointer :: s
         1 format(a32, 2x, 1pe26.16)
         2 format(a32, 2x, i9)
         11 format(3x, f15.8, i15)
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*, *) 'write_model: star_ptr failed'
            return
         end if
         
         species = s% species
         nz = s% nz
         
         write(io_ms_index, 11) s% star_mass, nz

         ! write property list
         write(io_ms_mod, 1) 'M/Msun', s% star_mass
         write(io_ms_mod, 2) 'n_shells', nz
         write(io_ms_mod, *) ! blank line for end of property list

         write(io_ms_mod, fmt='(7x, a9, 1x, 99(a24, 1x))', advance='no') &
            'lnd', 'lnT', 'lnR', 'L', 'dq'
         do j=1, species
            write(io_ms_mod, fmt='(a24, 1x)', advance='no') trim(chem_isos% name(s% chem_id(j)))
         end do
         write(io_ms_mod, *)
         do k=1, nz
            write(io_ms_mod, fmt='(i5, 1x)', advance='no') k
            write(io_ms_mod, fmt='(99(1pe24.16, 1x))', advance='no') &
                  s% lnd(k), s% lnT(k), s% lnR(k), s% L(k), s% dq(k)
            do j=1, species
               write(io_ms_mod, fmt='(1pe24.16, 1x)', advance='no') s% xa(j, k)
            end do
            write(io_ms_mod, *)
         end do      
         write(io_ms_mod, *)  
      
      end subroutine write_model
      
      
      subroutine evolve_to_zams(s, id, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer, parameter :: lipar=0, lrpar=0
         integer :: ipar(lipar)
         double precision :: rpar(lrpar)
         call star_evolve_to_check_point( &
            id, before_evolve_to_zams, evolve_to_zams_check_model, lipar, ipar, lrpar, rpar, ierr)
      end subroutine evolve_to_zams


      subroutine before_evolve_to_zams(s, id, lipar, ipar, lrpar, rpar, ierr)
         use star_def, only:star_info
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout) :: ipar(lipar)
         double precision, intent(inout) :: rpar(lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         !s% matrix_type = 3 ! use block tridiagonal for first step
      end subroutine before_evolve_to_zams
      
      
      integer function evolve_to_zams_check_model(s, id, lipar, ipar, lrpar, rpar)
         use star_def, only:star_info
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout) :: ipar(lipar)
         double precision, intent(inout) :: rpar(lrpar)
         evolve_to_zams_check_model = bare_bones_check_model(id) 
         if (evolve_to_zams_check_model /= keep_going) return
         if (s% L_nuc_burn_total >= s% L(1)/Lsun) evolve_to_zams_check_model = terminate
         s% matrix_type = 2
      end function evolve_to_zams_check_model


      subroutine read_zams_controls(s, zams_inlist, ierr)
         use utils_lib
         type (star_info), pointer :: s
         character (len=*), intent(in) :: zams_inlist
         integer, intent(out) :: ierr

         character (len=256) :: filename, message
         integer :: unit
         
         11 format(a30, f16.6)
         
         ierr = 0
         
         ! set defaults
         create_z = 2d-2
         zams_name = 'z2m2'
         dmass = 0.1d0
         mlo = 0
         mhi = mlo + 2*dmass

         unit=alloc_iounit(ierr)
         if (ierr /= 0) return

         filename = zams_inlist
         open(unit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(message, *) 'Failed to open control namelist file ', trim(filename)
            call alert(ierr, message)
            write(*,*) trim(message)
         else
            read(unit, nml=create_zams_job, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(message, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') trim(message)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(unit=unit, file=trim(filename), action='read',  &
                  delim='quote', status='old', iostat=ierr)
               read(unit, nml=create_zams_job)
               close(unit)
               call alert(ierr, message)
            end if  
         end if
         call free_iounit(unit)

      end subroutine read_zams_controls


      end module create_zams
     
     

     
