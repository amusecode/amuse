subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  real(dp)::d,u,v,w,A,B,C,e
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
  write(ilun)nvar+3
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
              do ivar=1,4
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                    end do
                 else ! Write velocity field
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ivar)/uold(ind_grid(i)+iskip,1)
                    end do
                 endif
                 write(ilun)xdp
              end do
              do ivar=6,8 ! Write left B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
              do ivar=nvar+1,nvar+3 ! Write right B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
              do i=1,ncache ! Write pressure
                 d=uold(ind_grid(i)+iskip,1)
                 u=uold(ind_grid(i)+iskip,2)/d
                 v=uold(ind_grid(i)+iskip,3)/d
                 w=uold(ind_grid(i)+iskip,4)/d
                 A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                 B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                 C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                 e=uold(ind_grid(i)+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
                 xdp(i)=(gamma-1d0)*e
              end do
              write(ilun)xdp
#if NVAR > 8
              do ivar=9,nvar ! Write passive scalars if any
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)/uold(ind_grid(i)+iskip,1)
                 end do
                 write(ilun)xdp
              end do
#endif
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
     
end subroutine backup_hydro





