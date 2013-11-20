subroutine init_hydro
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  real(dp)::d,u,v,w,A,B,C,e
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_hydro'
  
  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar+3))
  allocate(unew(1:ncell,1:nvar+3))
  uold=0.0d0; unew=0.0d0
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(nvar2.ne.(nvar+3))then
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar+3
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
              write(*,*)'File hydro.tmp is not compatible'
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
                 ! Loop over conservative variables
                 do ivar=1,4
                    read(ilun)xx
                    if(ivar==1)then ! Read density
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else  ! Read velocity field
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*uold(ind_grid(i)+iskip,1)
                       end do
                    end if
                 end do
                 do ivar=6,8 ! Read left B field
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
                 do ivar=nvar+1,nvar+3 ! Read right B field
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
                 read(ilun)xx ! Read pressure
                 do i=1,ncache
                    e=xx(i)/(gamma-1d0)
                    d=uold(ind_grid(i)+iskip,1)
                    u=uold(ind_grid(i)+iskip,2)/d
                    v=uold(ind_grid(i)+iskip,3)/d
                    w=uold(ind_grid(i)+iskip,4)/d
                    A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                    B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                    C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                    uold(ind_grid(i)+iskip,5)=e+0.5*d*(u**2+v**2+w**2)+0.5*(A**2+B**2+C**2)
                 end do
#if NVAR > 8
                 do ivar=9,nvar ! Read passive scalars if any
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)*uold(ind_grid(i)+iskip,1)
                    end do
                 end do
#endif
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)
#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'

  end if

end subroutine init_hydro




