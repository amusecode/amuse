subroutine init_radiation
  use amr_commons
  use hydro_commons
  use cooling_module, ONLY: force_j0_one
  use radiation_commons, ONLY: Erad,Srad
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer::nvar_expected
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_hydro'
  
  !------------------------------------------
  ! Allocate radiation arrays on the AMR grid
  !------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(Erad(1:ncell))
  allocate(Srad(1:ncell))
  Erad=0.0  ! Photon density = 0 initially.
  Srad=0.0  ! No initial sources.
  
  ! Force constant UV bkg
  force_j0_one=.true.

  !-----------------------------
  ! For a restart, read rad file
  !-----------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/rad_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2

     nvar_expected = 1

     if(nvar2.ne.nvar_expected)then
        write(*,*)'File rad.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar_expected
        call clean_stop
     end if

     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File rad.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 read(ilun)xx
                 do i=1,ncache
                    Erad(ind_grid(i)+iskip)=xx(i)
                 end do
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)
#ifndef WITHOUTMPI
     if(debug)write(*,*)'rad.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'RAD backup files read completed'
  end if

end subroutine init_radiation




