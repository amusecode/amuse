subroutine read_params
  use amr_commons
  use pm_parameters
  use poisson_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,ierr,levelmax
  character(LEN=80)::infile
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  logical::nml_ok

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/run_params/cosmo,pic,poisson,hydro,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering,static,geom
  namelist/output_params/noutput,foutput,fbackup,aout,tout,output_mode
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen
  namelist/poisson_params/epsilon,gravity_type,gravity_params
  namelist/halo_parameters/Mvir,overdensity
  namelist/jet_parameters/r_jet,h_jet,e_jet,t_jet

  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Carefull with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif
  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  if(myid==1)then
  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 2.0                            '
  write(*,*)'       written by Romain Teyssier (CEA/DSM/DAPNIA/SAP)         '
  write(*,*)'                     (c) CEA 1999-2005                         '
  write(*,*)' '
  write(*,'(" Working with nproc = ",I4," for ndim = ",I1)')ncpu,ndim
  write(*,*)' '

  ! Read namelist filename from command line argument
  narg = iargc()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: hydro3d input.nml'
     write(*,*)'File input.nml should contain a parameter namelist'
     call clean_stop
  END IF
  CALL getarg(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)then
        write(*,*)'File '//TRIM(infile)//' does not exist'
     endif
     call clean_stop
  end if

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue
  rewind(1)
  read(1,NML=halo_parameters)
  rewind(1)
  read(1,NML=jet_parameters)

  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nml_ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=ngridtot/int(ncpu,kind=8)
     endif
  end if
  if(npartmax==0)then
     npartmax=nparttot/int(ncpu,kind=8)
  endif
  if(myid>1)verbose=.false.
  if(pic.and.(.not.poisson))then
     poisson=.true.
  endif

  call read_hydro_params(nml_ok)

  close(1)

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if
  
  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
  end do
     
  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params

