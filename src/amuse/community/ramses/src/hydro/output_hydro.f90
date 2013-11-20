subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10
     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
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
              do ivar=1,nvar
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                    end do
                 else if(ivar>=2.and.ivar<=ndim+1)then ! Write velocity field
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                    end do                    
                 else if(ivar==ndim+2)then ! Write pressure
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ndim+2)
                       xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,2)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#if NDIM>1
                       xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,3)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NDIM>2
                       xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,4)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
                       xdp(i)=(gamma-1d0)*xdp(i)
                    end do                                        
                 else ! Write passive scalars if any
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                 endif
                 write(ilun)xdp
              end do
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
     
end subroutine backup_hydro





