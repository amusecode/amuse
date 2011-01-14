!> \file ev_main.f90  Contains the main routine
!! 
!! \mainpage Program EV
!! An incarnation of Peter Eggleton's TWIN binary evolution code.
!!
!! \par
!! This version of the code can be run as a standalone program, meaning
!! it doesn't need symbolic links or environment variables, although it will
!! continue to use those if found.
!!
!! \par Usage:
!!   ev name [metalicity] [evpath]
!!
!! \par
!! - Stores output files with the basename "name".
!! - The metalicity can be specified on the command line as a pure number, or
!! through the environment variable Z. If specified on the command line, it
!! overrides the environment variable. If not specified at all, the metalicity
!! defaults to solar (0.02).
!! - The evpath parameter identifies the location of the code datafiles, which are
!! expected in evpath/input/. Defaults to the current directory if not specified.
!! Settings are read from files init.run and init.dat in the current directory.
!! If these don't exist, "name".dat and "name".run will be used instead.
!! An alternative starting model can be specified (by name) using the STARTFILE
!! parameter in the .run file.
!!
!! \par Examples:
!! - Run the code, saving the output file as "sun":
!!   - ev sun
!! - As above, but change the metalicity to 0.0001:
!!   - ev sun 0001
!!
!! \see 
!! - The working of this programme is described in a file `writeup.tex'
!! - http://www.phys.ualberta.ca/~sluys/index.php?title=ev
!<


!> \brief The main routine
program ev_main
   use real_kind
   use mesh
   use mesh_enc
   use settings
   use constants
   use file_exists_module
   use init_run
   use svn_version
   use current_model_properties
   use starting_values0
   use filenames

   implicit none

   !Local variables:
   real(double) :: sm1, dty1, age1, per1, bms1, ecc1, p1, enc1

   integer :: i,j,kdm,jop
   integer :: kh1,kp1,jmod1,jb1,jn1,jf1,kh2,kp2,jmod2,jb2,jn2,jf2,jm1,jm2,k
   integer :: kpp,knt,jo1,jip,jc1,ksv,jo3,jo2,jc2,system
   real(double) :: mlo,dm,mhi,tn,perc,dmm

   real(double) :: sm2, dty2, age2, per2, bms2, ecc2, p2, enc2
   real(double) :: h1(nvar,nm), h2(nvar,nm)
   real(double) :: hn1(50, nm), hn2(50, nm)
   real(double) :: dh1(nvar,nm), dh2(nvar,nm)

   !Create different filenames for different loops:
   integer :: nml,nql,nxl,status
   character :: fname2*6,tmpstr*80
   integer :: dateval(8)
   
   real :: time0
   
   ! Time the start of execution:
   call cpu_time(time0)
   
   write (6, *) svn_version_string(:)
   if(debug.ge.1) write(6,'(A,I5)')' Running in debug mode',debug

   ! No arguments specified, abort and complain.
   ! This is a change of behaviour for the code
   if (command_argument_count() == 0) then
      write(6,'(A)')' '
      write(6,'(A)')' Calculate the evolution of single or binary stars.'
      write(6,'(A)')' Originally developed by Peter Eggleton (1971, 2002)'
      write(6,'(A)')' '
      write(6,'(A)')' Usage:'
      write(6,'(A)')'   ev name [metalicity] [evpath]'
      write(6,'(A)')' '
      write(6,'(A)')' Output is stored in files with the basename "name".'
      write(6,'(A)')' '
      write(6,'(A)')' The metalicity can be specified on the command line as a pure number, or'
      write(6,'(A)')' through the environment variable Z. If specified on the command line, it'
      write(6,'(A)')' overrides the environment variable. If not specified at all, the metalicity'
      write(6,'(A)')' defaults to solar (0.02).'
      write(6,'(A)')' '
      write(6,'(A)')' The evpath parameter identifies the location of the code datafiles, which are'
      write(6,'(A)')' expected in evpath/input/. Defaults to the current directory if not specified.'
      write(6,'(A)')' Settings are read from files init.run and init.dat in the current directory.'
      write(6,'(A)')' If these don'//"'"//'t exist, "name".dat and "name".run will be used instead.'
      stop
   end if

   ! Check for command-line arguments and determine the output filename
   call get_command_argument(1,basename)

   ! Setup metalicity and evolution code directory.
   ! In order of precedence: command line option, environment variable, default value
   if (command_argument_count() >= 2) then
      call get_command_argument(2,zstr)
   else
      call get_environment_variable("z", zstr)
   end if
   if (len(trim(zstr)) == 0) zstr = "02"
   if (command_argument_count() >= 3) then
      call get_command_argument(3,evpath)
   else
      call get_environment_variable("evpath", evpath)
   end if
   if (len(trim(evpath)) == 0) evpath = "."

   ! Default names for input files
   call set_default_filenames

   ! Check whether all files that are required exist
   call assert_input_files_exist


   ! We need to have a fort.11
   write (11, *) 0
   close (11)

   ! We need a scratch file for output we don't want or need
   if (.not. file_exists('fort.25')) open (unit=25, action='write', status = 'scratch')

   call setsup
   call load_nucleosynthesis_rates(inputunits(13), inputunits(14))
92 format (1x, 1p, 40d23.15, 0p)
95 format (1x, 1p, 8d23.15, 0p, 6i6)
   read (19, *) mlo, dm, mhi, kdm
   ! Read data for a triple set of loops, in mass, mass ratio, and period
   call read_init_run(23)

   call open_standard_output_files
   if (mutate) open(unit = 65, file=trim(basename)//'.pmutate')
   do i=1, isb
      call open_star_output_files(i)
      !> \todo FIXME: move this call to later, wftare we load init.dat and know whether we actually need nucleosynthesis output.
      call open_nucsyn_output_files(i)
   end do

   ! Write out termination (JO) codes on fort.8 for later reference
   write (8,*) ' -2 -- BEGINN -- requested mesh too large'
   write (8,*) ' -1 -- STAR12 -- no timesteps required'
   write (8,*) '  0 -- STAR12 -- finished required timesteps OK'
   write (8,*) '  1 -- SOLVER -- failed; backup, reduce timestep'
   ! JO = 1 won't stop the code, but it may lead to JO = 2
   write (8,*) '  2 -- BACKUP -- tstep reduced below limit; quit'
   write (8,*) '  3 -- NEXTDT -- *2 evolving beyond last *1 model'
   write (8,*) '  4 -- PRINTB -- *1 rstar exceeds rlobe by limit'
   write (8,*) '  5 -- PRINTB -- age greater than limit'
   write (8,*) '  6 -- PRINTB -- C-burning exceeds limit'
   write (8,*) '  7 -- PRINTB -- *2 radius exceeds rlobe by limit'
   write (8,*) '  8 -- PRINTB -- close to He flash'
   write (8,*) '  9 -- PRINTB -- massive (>1.2 msun) deg. C/O core'
   write (8,*) ' 10 -- PRINTB -- |M1dot| exceeds limit'
   write (8,*) ' 11 -- NEXTDT -- impermissible FDT for *2'
   write (8,*) ' 14 -- PRINTB -- funny compos. distribution'
   write (8,*) ' 15 -- STAR12 -- terminated by hand'
   write (8,*) ' 16 -- MAIN   -- ZAHB didnt converge'
   write (8,*) ' 17 -- BEGINN -- Nucleosynthesis didn'//"'"//'t converge'
   write (8,*) ' 51 -- PRINTB -- end of MS (core H below limit)'
   write (8,*) ' 52 -- PRINTB -- Radius exceeds limit'
   write (8,*) ' 53 -- PRINTB -- Convergence to target model reached minimum'
   write (8,*) ' 12, 22, 32 -- as 2, core H, He, C non-zero, resp.'
   write (8,*) ''
   write (8,*) 'CODE*1 C*2  log(Mi)   log(Qi)   log(Pi)  Modl  JMOD 27-JOP:'
   jop = 14
   if ( isb == 1 ) ktw = 1
   if ( isb == 1.or.ktw == 2 ) kp = kpt
   if ( ip1.le.15.or.ip2.le.14 ) kml = 1
   if ( ip1.le.14.or.ip2.le.15 ) kql = 1
   if ( ip1.le.14.or.ip2.le.14 ) kxl = 1
   if ( ip1.le.14 ) ip2 = ip1

   ! Create a list of model properties for this run, including date and pwd:
   open(unit=50, status='unknown', form='formatted', file=trim(basename)//'.list')
   write (50, *) svn_version_string(:)
   write(50,'(a1)')' '
   ! Get date and time for list file
   call date_and_time(tmpstr, tmpstr, tmpstr, dateval)
   write(50, '(i4,"-",i2,"-",i2, " ", i2,":",i2,":",i2, " ut",a5)') dateval(1:3), dateval(5:7), tmpstr
   ! Write hostname
   call get_environment_variable("HOSTNAME", tmpstr);
   write(50,'(a80)') tmpstr
   ! Write current working directory
   call get_environment_variable("PWD", tmpstr);
   write(50,'(a80)') tmpstr
   ! Write names of files used
   do i=1, n_inp_files
      if (file_exists(inputfilenames(i)) .and. .not. file_exists(basename)) then
         write (50, *) "Using ",trim(inputfilenames(i))," for unit", inputunits(i)
      end if
   end do
   if ( len(trim(startfile))>0 ) then
      write (50, *) "Using initial model (*1) ",trim(startfile)," for unit", ip1
   end if
   if ( len(trim(startfile2))>0 ) then
      write (50, *) "Using initial model (*2) ",trim(startfile2)," for unit", ip2
   end if
   write(50,'(a1)')' '
   write(50,'(16x,a4,4(8x,a3))')'file','m1i','m2i',' qi',' pi'
   close(50)
   

   !BEGIN LOOP
   ! Cycle on *1's mass.
   do nml = 1, kml
      ml = ml1 + (nml-1)*dml
      
      ! Locate the *1 model, optionally in the ZAMS (fort.16) library
      if ( ip1==16 .and. len(trim(zamsfile))>0 ) then
         ! Override the ZAMS file from which the starting model is to be obtained; Necessary for stars > 200Msun
         ! FIXME: what exactly is this doing that cannot be done by specifying the startfile???
         open(unit=50, file=trim(zamsfile))
         im1 = 1
         call load_star_model(50,im1, h1, dh1, hn1, sm1,dty1,age1,per1,bms1,ecc1,p1,enc1,kh1,kp1,jmod1,jb1,jn1,jf1)
         close(50)
      else
         if ( ip1 == 16 ) then
            !If no mass loop and setting M1 by hand, determine ZAMS-library mass automatically (ignore im1 in init.run):
            if(kml <= 1 .and. t0sm > 0) ml = log10(t0sm)
            im1 = (ml - mlo)/dm + 1.501_dbl
         end if
         call load_star_model(ip1,im1, h1, dh1, hn1, sm1,dty1,age1,per1,bms1,ecc1,p1,enc1,kh1,kp1,jmod1,jb1,jn1,jf1)
      end if
      if ( jmx == 0.and.ip1 == 16 ) sm1 = 10.0_dbl**ml

      ! Initialise some parameters relating to *1
      ! TN = nuclear timescale, PERC = log period for RLOF on ZAMS
      tn = 1.0d10*h1(4, 1)/h1(8, 1) * clsn/cmsn
      perc = sqrt( exp(3.0 * h1(7, 1))/(8.157_dbl*h1(4, 1)) * cmsn/crsn**3 )


      ! Cycle on mass ratio.
      do nql = 1, kql
         ql = ql1 + (nql-1)*dql
         sm2 = sm1/10.0_dbl**ql
         if ( isb == 2 ) then
            ! Locate the *2 model, optionally in the ZAMS (fort.16) library
            if ( ip2 == 16 ) im2 = (ml - ql - mlo)/dm + 1.501_dbl
            call load_star_model(ip2,im2, h2, dh2, hn2, sm2,dty2,age2,per2,bms2,ecc2,p2,enc2,kh2,kp2,jmod2,jb2,jn2,jf2)
            if ( jmx == 0.and.ip1 == 16 ) sm2 = sm1/10.0_dbl**ql

            ! Store total binary mass in case we're setting up the binary from previous models
            if ( ip1 /= 16 .and. ip2 /= 16 ) then
               bms1 = sm1 + sm2
               bms2 = sm1 + sm2
            end if
         end if


         ! Cycle on initial period
         do nxl = 1, kxl
            xl = xl1 + (nxl-1)*dxl

            !BEGIN ACTUAL LOOP
            !Make a list of models in each file
            write(fname2,'(3i2.2)')nml,nql,nxl
            open(unit=50, status='unknown', position='append',form='formatted', file=trim(basename)//'.list')
            if(kml*kql*kxl == 1) then
               write(50,'(a20,4(2x,es9.3))')trim(basename),10**ml,10**(ml-ql),10**ql,10**xl
            else
               write(50,'(a20,4(2x,es9.3))')trim(basename)//fname2,10**ml,10**(ml-ql),10**ql,10**xl
            end if
            close(50)


            ! For ZAMS models, replace SMn, ..., JMn by correct values
            if ( ip1 == 16 ) then
               age1 = 0.0_dbl
               per1 = 0.0_dbl
               bms1 = 0.0_dbl
               ecc1 = 0.0_dbl
               p1   = h1(13,1)
               enc1 = 0.0_dbl
               dty1 = 1.0d-4*tn
               bms1 = sm1 + sm2
               per1 = perc*10.0_dbl**xl
               if ( kr == 1 ) p1 = perc*10.0_dbl**rot
               if ( kr == 2 ) p1 = dmax1(perc, per1*10.0_dbl**rot)
               if ( kr >= 3 ) p1 = per1 !SdM almost sync. rotation
               ecc1 = ex
               jm1 = 0
            end if
            if ( ip2 == 16 ) then
               dty2 = dty1
               age2 = age1
               per2 = per1
               bms2 = bms1
               ecc2 = ecc1
               p2   = p1
               enc2 = enc1
               jm2 = 0
            end if

            ! Optionally replace SM, ..., JM by values AM, DTY, .., JMX from fort.23
            if ( kr == 3 ) then
               t0per = t0per*((1.0_dbl - t0ecc**2)/(1.0_dbl - ex**2))**1.5_dbl
               t0ecc = ex
            end if

            if (t0sm >= 0.0_dbl)  sm1  = t0sm
            if (t0dty >= 0.0_dbl) dty1 = t0dty
            if (t0age >= 0.0_dbl) age1 = t0age
            if (t0per >= 0.0_dbl) per1 = t0per
            if (t0bms >= 0.0_dbl) bms1 = t0bms
            if (t0ecc >= 0.0_dbl) ecc1 = t0ecc
            if (t0p >= 0.0_dbl)   p1   = t0p
            if (t0enc >= 0.0_dbl) enc1 = t0enc

            if (t0sm >= 0.0_dbl)  sm2  = t0sm    ! Note: this is overruled a couple of lines below
            if (t0dty >= 0.0_dbl) dty2 = t0dty
            if (t0age >= 0.0_dbl) age2 = t0age
            if (t0per >= 0.0_dbl) per2 = t0per
            if (t0bms >= 0.0_dbl) bms2 = t0bms
            if (t0ecc >= 0.0_dbl) ecc2 = t0ecc
            if (t0p >= 0.0_dbl)   p2   = t0p
            if (t0enc >= 0.0_dbl) enc2 = t0enc

            if ( jmx >= 0 ) then  ! Use desired model number from init.run
               jm1 = jmx
               jm2 = jmx
            else
               jm1 = jmod1        ! Take number from existing model 
               if(isb==2) jm2 = jmod2
            end if
            
            sm2 = bms2 - sm1                    ! M2 = binary mass - M1
            if ( ktw == 1 ) p2 = 1.0d6          ! Set P_rot to a large number for non-TWIN secondaries

            ! Store the *1 and *2 initial models in one file (fort.14 = last2)
            write (jop, 95)  sm1, dty1, age1, per1, bms1, ecc1, p1, enc1, kh1, kp, jm1, 1, jn1, jf1
            do k = 1, kh1        ! Model data
               write (jop, 92) (h1(j, k), j = 1, jn1)
            end do
            if (iand(jf1, 4) == 4) then
               do k = 1, kh1     ! Derivative information
                  write (jop, 92) (dh1(j, k), j = 1, jn1)
               end do
            end if
            if (iand(jf1, 8) == 8) then
               do k = 1, kh1     ! Nucleosynthesis abundances
                  write (jop, '(1x, 1p, 50d23.15, 0p)') (hn1(j, k), j = 1, 50)
               end do
            end if
            if ( isb /= 1 ) then
               kpp = 5000
               if ( ktw == 2 ) kpp = kp
               write (jop, 95)  sm2, dty2, age2, per2, bms2, ecc2, p2, enc2, kh2, kpp, jm2, 2, jn2, jf2
               do k = 1, kh2        ! Model data
                  write (jop, 92) (h2(j, k), j = 1, jn2)
               end do
               if (iand(jf2, 4) == 4) then
                  do k = 1, kh2     ! Derivative information
                     write (jop, 92) (dh2(j, k), j = 1, jn2)
                  end do
               end if
               if (iand(jf2, 8) == 8) then
                  do k = 1, kh2     ! Nucleosynthesis abundances
                     write (jop, '(1x, 1p, 50d23.15, 0p)') (hn2(j, k), j = 1, 50)
                  end do
               end if
            end if
            call flush (jop)
            rewind (jop)

            knt = 0
            jo1 = -1
            do
               jip = jop
               if (mutate) then
                  ! Do mutation run of the code
                  call ms2bss ( jop, jo1 )
                  if (jo1 /= -1) exit
                  ! exit the code after constructing a starting model if that's all we want
                  if (start_model_only) then
                     jo1 = 0
                     exit
                  end if
               end if

               ! Do alternately KP steps of *1, then sufficient steps of *2 for *2 to
               ! catch up, then back to *1 ..., until either KPT steps are taken, or
               ! one star terminates. The `termination code' JO is possibly altered
               ! by various subroutines: see above.
               call star12 ( jo1, jc1, jop, jip, ksv, 22 )
               rewind (22)

               ! He flash (JO=8): replace last model by ZAHB model, then continue
               if ( jo1 == 8 ) then
                  jo3 = -1
                  call fgb2hb ( jop, jo3 )
                  jo1 = jo3   ! Pass on ZAHB flag to STAR12
                  if ( jo3 == 13 ) cycle
                  jo1 = 16
               end if

               jm1 = jmod
               rewind (1)
               call pruner ( 1, 3, isb )  !Unit 1,3 = out1, io12
               close (3)
               close (1)
               ! Next line was: unit 1 = out1 - FIXME: does this only happen for star 1? - this goes wrong for loops
               open(unit = 1, file=trim(basename)//'.out1',status='old',position='append')
               jo2 = -1
               jmod = 0

               ! Cycle for *2, but only if doing a binary in normal (ie, non-TWIN) mode
               if ( isb == 2 .and. ktw == 1 ) then
                  jip = jop
                  call star12 ( jo2, jc2, jop, jip, ksv, 22 )
                  call flush (3)
                  call flush (jop)
                  call flush (27 - jop)
                  rewind (3)
                  rewind (jop)
                  rewind (27 - jop)
                  rewind (22)
                  knt = knt + kp

                  ! Loop back if *2 has caught up with *1, etc
                  if ( knt < kpt .and. jo1 == 0 .and. jo2 == 3 ) then
                     jop = 27 - jop
                     cycle
                  end if
               end if

               exit
            end do  !do


            ! Write final answer on unit 9; initial values and termination code on unit 8
            do i=1, isb
               rewind (i)
               call pruner ( i, 9, isb )
               close (i)
               if(i.eq.1) open(unit = 1, file=trim(basename)//'.out1',status='old',position='append')  !unit=1, file=out1
               if(i.eq.2) open(unit = 2, file=trim(basename)//'.out2',status='old',position='append')  !unit=2, file=out2
            end do
            call flush ( 9 )
            rewind (22)
            rewind (jop)
            rewind (27 - jop)
            if ( jo1 == 2 ) jo1 = jo1 + 10*jc1
            if ( jo2 == 2 ) jo2 = jo2 + 10*jc2
            write (8,'(2i5, 3f10.5, 7i5)') jo1, jo2, ml, ql, xl, jm1, jmod, 27 - jop
            call flush ( 8 )

            ! Write description of the mass grid to file.mas. Really only useful for ZAMS libraries though
            dmm = -log10(1.0_dbl - cmi*t0dty*csy)*ksv
            write (29,'(3f10.5, i5)') ml, dmm, ml + dmm*kpt/ksv, ksv/kt4
            call flush (29)

            !If the code was started for more than one model, rename the files for each model
            !  write(fname2,'(3i2.2)')nml,nql,nxl
            !  if(kml*kql*kxl.gt.1) then
            !  do i=1,nout(isb)
            !    close(units(i))
            !    status = system('mv -f '//trim(basename)//suffix(i)//' '//trim(basename)//fname2//suffix(i))
            !  end do
            !  end if

            !If the code was started for more than one model, create a directory per model and move the files
            if(kml*kql*kxl.gt.1) then
               status = system('mkdir '//trim(basename)//fname2)
               do i=1,nout(isb)
                  close(units(i))
                  status = system('mv -f '//trim(basename)//suffix(i)//' '//trim(basename)//fname2)
               end do
            end if
            
         end do  ! nxl - orbital period
         
      end do  ! nql - mass ratio
      
   end do  ! nml - primary mass


   !End of a run loop
   call print_exit_message(jo1,jo2,time0)
   
end program ev_main

