!> module save_init
!!
!! Saves the value of kpt (the maximum number of timesteps in the model) from init.run as s_kpt.
!!
!! This variable used to be in the common block saveinit
module save_init
   use real_kind

   implicit none
   integer, save :: s_kpt  !Saved value of kpt
end module save_init






!> Module init_run
!!
!! Contains:
!! - input settings from init.run  (formerly stored in the common block cinit_run)
!! - read_init_run() to read an init.run file from disc
!<
module init_run

   use real_kind

   implicit none

   integer :: ip1               ! Input unit for initial model for star 1 (13-16)
   integer :: im1               ! Model number for initial model for star 1
   integer :: ip2               ! Input unit for initial model for star 2 (13-16)
   integer :: im2               ! Model number for initial model for star 2

   integer :: kpt               ! Maximum number of timesteps for each component
   integer :: kp                ! Number of timesteps before switching to the other component in a  non-TWIN binary

   integer :: kml               ! Number of iterations to perform in the primary-mass loop
   integer :: kql               ! Number of iterations to perform in the mass-ratio loop
   integer :: kxl               ! Number of iterations to perform in the orbital-period loop

   integer :: kr                ! Used to set nitial rotation of the star (see user manual)
   integer :: jmx               ! Starting model number (see user manual)


   real(double) :: ml1                ! log(M) for the primary in the first iteration in the primary-mass loop
   real(double) :: dml                ! dlog(M) to add to log(M) in subsequent iterations

   real(double) :: ql1                ! log(q) for the mass ratio in the first iteration in the mass-ratio loop
   real(double) :: dql                ! dlog(q) to add to log(q) in subsequent iterations

   real(double) :: xl1                ! log(P) for the orbital period in the first iteration in the orbital-period loop
   real(double) :: dxl                ! dlog(P) to add to log(P) in subsequent iterations

   real(double) :: rot                ! Initial rotation (see user manual)
   real(double) :: ex                 ! Initial exentricity


   character(len=500) :: startfile    ! File name to read starting model from (*1)
   character(len=500) :: startfile2   ! File name to read starting model from (*2)


contains



   ! ------------------------------------------------------------------------------
   ! READ_INIT_RUN
   !  Read the initial conditions for the current run from unit IR
   ! ------------------------------------------------------------------------------
   !  Input:
   !     IR    -  The FORTRAN unit from which to read the input (typically 23)
   !  Output:
   !     Various options in various common blocks (it's a bit of a mess)
   ! ------------------------------------------------------------------------------
   subroutine read_init_run(ir)
   use real_kind
   use mesh
   use mesh_enc
   use settings
   use file_exists_module
   use current_model_properties
   use starting_values0
   use save_init
   use stopping_conditions
   use polytrope, only: NMP

   implicit none
   integer, intent(in) :: ir

   integer :: ioerror
   character(len=500) :: name
   real(double) :: sm, dty, age, per, bms, ecc, p, enc
   ! Temporary (dummy) values used for importing merger models.
   real(double) :: t1sm, t1dty, t1age, t1per, t1bms, t1ecc, t1p, t1enc
   integer :: kh1


   ! EG: Namelist for namelist I/O of init.run (fort.23)
   namelist /init_run/ isb, ktw, ip1, im1, ip2, im2, kpt, kp, ml1, dml, kml,  &
         ql1, dql, kql, xl1, dxl, kxl, rot, kr, ex, &
         sm, dty, age, per, bms, ecc, p, enc, jmx, &
         uc, start_model_only, composition_only, startfile, startfile2, &
         start_with_rigid_rotation, stop_at_target_mass

   ! first try to use namelist i/o; if that fails use old-fashioned io
   startfile = ''
   startfile2 = ''
   start_model_only = .true.
   composition_only = .false.
   read (ir, nml=init_run, iostat=ioerror)
   rewind (ir)

   !> \todo This is a temporary solution, let's see how we fix this properly later
   t0sm = sm
   t0dty = dty
   t0age = age
   t0per = per
   t0bms = bms
   t0ecc = ecc
   t0p = p
   t0enc = enc

   if (ioerror /= 0) then
      read (ir, *) isb, ktw, ip1, im1, ip2, im2, kpt, kp, ml1, dml, kml,   &
                   ql1, dql, kql, xl1, dxl, kxl, rot, kr, ex, &
                   t0sm, t0dty, t0age, t0per, t0bms, t0ecc, t0p, t0enc, jmx, uc
      rewind (ir)
   end if

   ! If a starting model was named in the input file for star 1, use that instead.
   if ( len(trim(startfile))>0 ) then
      if (.not. file_exists(startfile)) then
         inquire(unit=ir, name=name)
         write (0, *) 'start file "', trim(startfile), '" not found. check setting in ', trim(name), '.'
         stop
      else
         ip1 = 62
         open(unit = ip1, action="read", file=startfile)
         read (ip1, *) t1sm, t1dty, t1age, t1per, t1bms, t1ecc, t1p, t1enc, kh1
         rewind(ip1)
         max_nm = max(max_nm, kh1)
      end if
   end if

   ! If a starting model was named in the input file for star 2, use that instead.
   if ( len(trim(startfile2))>0 ) then
      if (.not. file_exists(startfile2)) then
         inquire(unit=ir, name=name)
         write (0, *) 'start file "', trim(startfile2), '" not found. check setting in ', trim(name), '.'
         stop
      else
         ip2 = 61
         open(unit = ip2, action="read", file=startfile2)
         read (ip2, *) t1sm, t1dty, t1age, t1per, t1bms, t1ecc, t1p, t1enc, kh1
         rewind(ip2)
         max_nm = max(max_nm, kh1)
      end if
   end if

   ! If starting from a polytrope make sure array sizes are adequate
   if (ip1 == 0 .or. ip2 == 0) max_nm = max(max_nm, NMP)

   ! EG: store some information we may later need so we don't have to reread the file
   s_kpt = kpt

   ! Sanity check: if not doing a binary set ktw=1
   if (isb == 1) ktw = 1

   ! EG: Do some magic for ISB == -1, -2 meaning we want to do a mutation rather than
   !  a regular evolution run. A mutation run is a run where we start with a
   !  normal star and end up with something weird that we want to evolve next.
   !  Useful for getting merger remnants, for instance, to evolve
   if (isb == -1) then

      ! Make sure all parts of the code that need to know know we're mutating
      mutate = .true.

      ! Keep the type of mutation we want to do
      mutation_mode = isb

      ! We really do want the code to calculate things for one star
      isb = 1

      ! Read in target model global data
      tprof_in = ip1
      read (tprof_in, *) t1sm, t1dty, t1age, t1per, t1bms, t1ecc, t1p, t1enc, kh1
      rewind (tprof_in)

      max_nm = max(max_nm, kh1)

      ! To start our mutation run, we first need a suitable model from the
      !  main sequence library, so we have to setup some variables to the
      !  mass that we want.
      ! If the mass and period specified in the init.run file are negative, then
      ! we adopt the values from the input model (as expected), otherwise we
      ! override those. This is mainly useful for ignoring rotation and to
      ! add "missing mass" from a partial SPH model.
      if (t0sm < 0.0) t0sm = t1sm
      if (t0p < 0.0) t0p = t1p
      ml1 = log10(t1sm)
      ip1 = 16
   end if
   end subroutine read_init_run



end module init_run



