      program make_merger_startup_model
      use twin_library
      use fudge_control
      use mesh_enc
      use constants
      use settings
      implicit none
      character*500 :: inputfilename, evpath, zstr, outputfilename, basename;
      double precision :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      integer :: KH, KTW, KW(260)
      integer :: ISB, IP1, IM1, IP2, IM2, KPT, KP, KML, KQL, KXL, KR, JMX
      double precision :: ML1, DML, QL1, DQL, XL1, DXL, ROT, EX
      double precision :: AX, SM, DTY, AGE, PER, BMS, ECC, P, ENC
      character*500 :: STARTFILE
      double precision :: mass, t, dt, entropy_max, target_diffsqr
      double precision :: composition_model(NVAR,NM)
      integer :: status, iter, id, IK, IJ
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!      COMMON /CINIT_RUN/ ISB, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
!     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, JMX, STARTFILE
!      COMMON /T0 / SM, DTY, AGE, PER, BMS, ECC, P, ENC
      integer :: number_of_meshpoints
      integer, parameter :: number_of_variables = 13
      double precision, allocatable :: model(:)
      
      print *, 'Stand-alone merger startup model creation program for TWIN'
      print *, 'Written by Evert Glebbeek 2007,2008,2009'

      if (iargc() == 0) then
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
      call getarg(1, inputfilename);

      ! Read metalicity environment variable
      call getenv("Z", zstr);
      if (len(trim(zstr)) == 0) then
         print *, "Warning: Metalicity Z not set!"
         print *, "Set the environment variable Z to the desired metalicity"
         print *, "Using default Z=02"
         zstr = "02"
      end if

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
      call set_init_run_name(trim(evpath)//'/stars_standards/init.run')
      call set_init_dat_name(trim(evpath)//'/input/mutate.dat')
      status = initialise_twin(trim(evpath), 20, trim(zstr));
      print *, 'TWIN initialisation finished with errnum', status
      if (status/=0) stop

      ! Read target model
      open(unit = 66, name=inputfilename, action = 'read', status = 'old')
      read (66, *) number_of_meshpoints

      ! Allocate enough memory to store the model
      allocate(model(number_of_meshpoints*number_of_variables))

      ! Read in all values
      read (66, *) model(1:number_of_meshpoints*number_of_variables)
      close (66)

      ! Convert model to TWIN input file
      call import_stellar_merger(number_of_meshpoints, number_of_variables,
     &                            model, 0.0d0)

      ! Store the output model
      outputfilename = trim(basename)//'.pmutate'
      print *, 'Writing output to ', trim(outputfilename)
      call dump_twin_model(id, outputfilename);

      end program
