subroutine backup_poisson(filename)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,info
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering backup_poisson'

  ilun=ncpu+myid+10
     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              ! Write potential
              do i=1,ncache
                 xdp(i)=phi(ind_grid(i)+iskip)
              end do
              write(ilun)xdp
              ! Write force
              do ivar=1,ndim
                 do i=1,ncache
                    xdp(i)=f(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
     
end subroutine backup_poisson





