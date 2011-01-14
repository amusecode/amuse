module filenames
   character, save :: evpath*500                   ! Path to the evolution code
   integer, parameter :: n_inp_files = 14          ! Number of input files
   integer, parameter :: n_general_out_files = 6   ! General output (always written)
   integer, parameter :: n_star_out_files = 3      ! Output for each star
   integer, parameter :: n_nucsyn_out_files = 3    ! Nucleosynthesis output for each star
   integer, parameter :: n_out_files = n_general_out_files + 2*(n_star_out_files + n_nucsyn_out_files)

   character(len=500) :: basename
   character(len=500) :: inputfilenames(n_inp_files)

   character, parameter :: suffix(n_out_files)*8 = (/                         &
      '.io12   ','.log    ','.out    ','.last1  ', '.last2  ', '.mod    ',    &
      '.out1   ','.plt1   ','.mdl1   ','.nucout1', '.nucplt1', '.nucmdl1',    &
      '.out2   ','.plt2   ','.mdl2   ','.nucout2', '.nucplt2', '.nucmdl2'/)

   ! Default metallicity (string), used for constructing filenames
   character, save :: zstr*8  = '02'

   ! Standard units for output files
   integer, parameter :: units(n_out_files) = (/ 3, 8, 9,13,14,15,      &  ! General
                                                 1,31,33,35,37,39,      &  ! Star 1
                                                 2,32,34,36,38,40/)        ! Star 2

   ! Standard units for input files
   integer, parameter :: inputunits(n_inp_files)     = (/12, 24, 16, 19, 20, 21, 26, 63, 22, 23, 41, 42, 43, 44/)

   ! Whether a file is required (>0), optional (0) or unneeded (-1)
   integer, parameter :: input_required(n_inp_files) = (/ 0,  0,  1,  1,  0, 1,  1, -1,  1,  1,  0,  1, 1, 1/)

   ! Number of output files, single stars and binaries
   integer, parameter :: nout(2) = (/ n_general_out_files+n_star_out_files+n_nucsyn_out_files, n_out_files /)

contains

   !> \brief Opens the default files directly by name, unless the associated fort.n exists
   subroutine set_default_filenames
   use file_exists_module
   implicit none
   ! Open input files; this replaces the symbolic links to fort.nnn files
   ! This is actually not the best way to do this (better would be to specify
   ! the filenames from a configuration file), but it'll do.
   ! Opens the file directly by name, unless the associated fort.nnn file exists
   inputfilenames(1)  = trim(evpath)//"/input/zahb"//trim(zstr)//".mod"
   inputfilenames(2)  = trim(evpath)//"/input/zahb"//".dat"
   inputfilenames(3)  = trim(evpath)//"/input/zams/zams"//trim(zstr)//".mod"
   inputfilenames(4)  = trim(evpath)//"/input/zams/zams"//trim(zstr)//".mas"
   inputfilenames(5)  = trim(evpath)//"/input/metals/z"//trim(zstr)//"/phys.z"//trim(zstr)
   inputfilenames(6)  = trim(evpath)//"/input/lt2ubv.dat"
   inputfilenames(7)  = trim(evpath)//"/input/nucdata.dat"
   inputfilenames(8)  = trim(evpath)//"/input/mutate.dat"
   inputfilenames(9)  = 'init.dat'
   inputfilenames(10) = 'init.run'
   inputfilenames(11) = trim(evpath)//"/input/COtables/COtables_z"//trim(zstr)
   inputfilenames(12) = trim(evpath)//"/input/physinfo.dat"
   inputfilenames(13) = trim(evpath)//"/input/rates.dat"
   inputfilenames(14) = trim(evpath)//"/input/nrates.dat"

   ! If init.dat does not exist, try name.dat; likewise for init.run
   if ( (.not. file_exists(inputfilenames(9))) .and. (file_exists(trim(basename)//".dat")) )  &
        inputfilenames(9)=trim(basename)//".dat"

   if ( (.not. file_exists(inputfilenames(10))) .and. (file_exists(trim(basename)//".run")) )  &
        inputfilenames(10)=trim(basename)//".run"
   end subroutine set_default_filenames


   !> \brief Check whether I/O files exist and act depending on whether they're required or not
   subroutine assert_input_files_exist
     use file_exists_module
     implicit none
     integer :: i
     character(len=500) :: fname
     
     do i=1, n_inp_files
        ! If a fort.n file exists for an input file, then use that in preference of the named file
        ! This is for backward compatibility
        write (fname, '("fort.",i2)') inputunits(i)
        if (file_exists(inputfilenames(i)) .and. .not. file_exists(fname)) then
           open(unit = inputunits(i), action="read", file=inputfilenames(i))
        end if
        
        ! Check if files exist and act appropriately
        ! If the file is required (INPUT_REQUIRED>0), abort on error
        ! If the file is optional (INPUT_REQUIRED==0), give a warning
        ! If the file is probably unneeded (INPUT_REQUIRED<0), do nothing
        if (.not. (file_exists(inputfilenames(i)) .or. file_exists(fname))) then
           if (input_required(i) > 0) then
              write (0, *) 'Required input file ', trim(inputfilenames(i)), ' (', trim(fname), ') not found'
              stop
           else if (input_required(i) == 0) then
              write (0, *) 'Warning: input file ', trim(inputfilenames(i)), ' (', trim(fname), ') not found'
           end if
        end if
     end do
   end subroutine assert_input_files_exist



   !> \brief Open standard output files
   !! 
   !! Only use named files if no fort.nnn file for the same unit exists
   !!
   !! Files opened (unit: extension):
   !! - 3: io12, 
   !! - 8: log, 
   !! - 9: out, 
   !! - 13: last1, 
   !! - 14: last2, 
   !! - 15: mod,
   !! - 29: mas
   !<
   subroutine open_standard_output_files
     use file_exists_module
     implicit none
     character :: fname*500
     integer :: i
     
     do i=1, n_general_out_files
        if (units(i)<10) then
           write (fname, '("fort.",i1)') units(i)
        else
           write (fname, '("fort.",i2)') units(i)
        end if
        if (.not. file_exists(fname)) then
           open(unit = units(i), file=trim(basename)//suffix(i))
        end if
     end do
     open(unit = 29, file=trim(basename)//'.mas')
   end subroutine open_standard_output_files




   !> \brief Open output files for star
   !! 
   !! Only use named files if no fort.nnn file for the same unit exists
   !! 
   !! Files opened (unit: extension):
   !! -  1: out1, 
   !! - 31: plt1, 
   !! - 33: mdl1, 
   !! -  2: out2, 
   !! - 32: plt2, 
   !! - 34: mdl2
   !!
   !! \param Jstar  Binary member (1 or 2)
   !<
   subroutine open_star_output_files(Jstar)
     use file_exists_module
     implicit none
     integer, intent(in) :: Jstar
     character :: fname*500
     integer :: i, j
     
     do j=1, n_star_out_files
        i = n_general_out_files + (Jstar-1) * (n_star_out_files+n_nucsyn_out_files) + j
        if (units(i)<10) then
           write (fname, '("fort.",i1)') units(i)
        else
           write (fname, '("fort.",i2)') units(i)
        end if
        if (.not. file_exists(fname)) then
           open(unit = units(i), file=trim(basename)//suffix(i))
        end if
     end do
   end subroutine open_star_output_files



   !> \brief Open output files for nucleosynthesis
   !! 
   !! Only use named files if no fort.nnn file for the same unit exists
   !!
   !! Files opened (unit: extension):
   !! - 35: nucout1, 
   !! - 37: nucplt1, 
   !! - 39: nucmdl1, 
   !! - 36: nucout2, 
   !! - 38: nucplt2, 
   !! - 40: nucmdl2, 
   !!
   !! \param Jstar  Binary member (1 or 2)
   !<
   subroutine open_nucsyn_output_files(Jstar)
     use file_exists_module
     implicit none
     integer, intent(in) :: Jstar
     character :: fname*500
     integer :: i, j
     
     do j=1, n_nucsyn_out_files
        i = n_general_out_files + (Jstar-1) * (n_star_out_files+n_nucsyn_out_files) + n_star_out_files + j
        if (units(i)<10) then
           write (fname, '("fort.",i1)') units(i)
        else
           write (fname, '("fort.",i2)') units(i)
        end if
        if (.not. file_exists(fname)) then
           open(unit = units(i), file=trim(basename)//suffix(i))
        end if
     end do
   end subroutine open_nucsyn_output_files
   
end module filenames
