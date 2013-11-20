subroutine init_poisson
  use pm_commons
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_poisson'

  !------------------------------------------------------
  ! Allocate cell centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(rho (1:ncell))
  allocate(phi (1:ncell))
  allocate(phi_old (1:ncell))
  allocate(f   (1:ncell,1:3))
  rho=0.0D0; phi=0.0D0; f=0.0D0
  if(cic_levelmax>0)then
     allocate(rho_top(1:ncell))
     rho_top=0d0
  endif

  !------------------------------------------------------
  ! Allocate multigrid variables
  !------------------------------------------------------
  ! Allocate communicators for coarser multigrid levels
  allocate(active_mg    (1:ncpu,1:nlevelmax-1))
  allocate(emission_mg  (1:ncpu,1:nlevelmax-1))
  do ilevel=1,nlevelmax-1
     do i=1,ncpu
        active_mg   (i,ilevel)%ngrid=0
        active_mg   (i,ilevel)%npart=0
        emission_mg (i,ilevel)%ngrid=0
        emission_mg (i,ilevel)%npart=0
     end do
  end do
  allocate(safe_mode(1:nlevelmax))
  safe_mode = .false.

  !--------------------------------
  ! For a restart, read poisson file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/grav_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     if(ndim2.ne.ndim)then
        write(*,*)'File poisson.tmp is not compatible'
        write(*,*)'Found   =',ndim2
        write(*,*)'Expected=',ndim
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
              write(*,*)'File poisson.tmp is not compatible'
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
                 ! Read potential
                 read(ilun)xx
                 do i=1,ncache
                    phi(ind_grid(i)+iskip)=xx(i)
                 end do
                 ! Read force
                 do ivar=1,ndim
                    read(ilun)xx
                    do i=1,ncache
                       f(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)
#ifndef WITHOUTMPI
     if(debug)write(*,*)'poisson.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'POISSON backup files read completed'
  end if

end subroutine init_poisson


