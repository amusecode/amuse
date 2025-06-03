program make_merger_startup_model
   use real_kind
   use twinlib
   use import

   implicit none
   character(len=500) :: inputfilename, evpath, zstr, outputfilename, basename
   integer :: status, id, IK
   integer :: number_of_meshpoints
   integer, parameter :: number_of_variables = 13
   character :: tmpstr*(80)
   real(double), allocatable :: model(:,:)
   real(double) :: Z

   print *, 'Stand-alone merger startup model creation program for TWIN'
   print *, 'Written by Evert Glebbeek 2007,2008,2009'

   if (command_argument_count() == 0) then
      print *, "Error: need an input filename"
      print *, "Usage: mkmergermod file"
      print *, "The specified file should contain the number of mass shells"
      print *, "on the first line, followed by the individual mass shells"
      print *, "starting at the surface and down to the centre."
      print *, ""
      print *, "Order of variables stored at each meshpoint:"
      print *, "Mass coordinate [Msun], Radius [Rsun], log density [cgs],"
      print *, "log pressure [cgs], XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE"
      stop
   end if
   call get_command_argument(1, inputfilename);

   ! Read metalicity environment variable
   call getenv("Z", zstr);
   if (len(trim(zstr)) == 0) then
      print *, "Warning: Metalicity Z not set!"
      print *, "Set the environment variable Z to the desired metalicity"
      print *, "Using default Z=02"
      zstr = "02"
   end if

   tmpstr = "0."//zstr
   read (tmpstr, *) Z

   ! Global path to the evolution code   
   evpath=""
   call getenv("evpath", evpath);
   if (len(trim(evpath)) == 0) then
      print *, "Warning: TWIN path evpath not set!"
      print *, "Set the environment variable evpath to point to the stars/ directory"
      print *, "Using default evpath=."
      evpath = "./"
   end if

   IK = len(trim(inputfilename))
   basename = inputfilename(1:ik-4);

   ! The init.run we read in here is just a stub really...
   status = initialise_twin(trim(evpath), 3, Z, verb=.false.)

   print *, 'TWIN initialisation finished with errnum', status
   if (status/=0) stop

   ! Read target model
   open(unit = 66, file=inputfilename, action = 'read', status = 'old')
   read (66, *) number_of_meshpoints

   ! Allocate enough memory to store the model
   allocate(model(number_of_variables, number_of_meshpoints))

   ! Read in all values
   do ik = 1, number_of_meshpoints
      read (66, *) model(1:number_of_variables, ik)
   end do
   close (66)

   ! Convert model to TWIN input file
   status = import_stellar_merger(id, number_of_meshpoints, number_of_variables, model, 0.0d0, .false.)

   ! Store the output model
   outputfilename = trim(basename)//'.pmutate'
   print *, 'Writing output to ', trim(outputfilename)
   status = write_star_to_file(id, outputfilename)

end program make_merger_startup_model
