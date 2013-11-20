subroutine init_hydro
  use amr_commons
  use hydro_commons
#ifdef RT      
  use rt_parameters,only: convert_birth_times
#endif
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
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_hydro'
  
  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
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
     if(.not.(neq_chem.or.rt) .and. nvar2.ne.nvar)then
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar
        call clean_stop
     end if
#ifdef RT
     if((neq_chem.or.rt).and.nvar2.lt.nvar)then ! OK to add ionization fraction vars
        ! Convert birth times for RT postprocessing:
        if(rt.and.static) convert_birth_times=.true.
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found nvar2  =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        if(myid==1) write(*,*)'..so only reading first ',nvar2, &
                  'variables and setting the rest to zero'
     end if
     if((neq_chem.or.rt).and.nvar2.gt.nvar)then ! Not OK to drop variables 
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found   =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        call clean_stop
     end if
#endif
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
#ifndef RT
                 do ivar=1,nvar
#else
                 do ivar=1,min(nvar,nvar2)
#endif
                    read(ilun)xx
                    if(ivar==1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*uold(ind_grid(i)+iskip,1)
                       end do
                    else if(ivar==ndim+2)then
                       do i=1,ncache
                          xx(i)=xx(i)/(gamma-1d0)
                          if (uold(ind_grid(i)+iskip,1)>0.)then                          
                             xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,2)**2/uold(ind_grid(i)+iskip,1)
#if NDIM>1
                             xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,3)**2/uold(ind_grid(i)+iskip,1)
#endif
#if NDIM>2
                             xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,4)**2/uold(ind_grid(i)+iskip,1)
#endif
                             else if(uold(ind_grid(i)+iskip,2) /= 0.)then 
                                write(*,*)'Problem in init_hydro with zero or negative density'
                                call clean_stop
#if NDIM>1
                             else if(uold(ind_grid(i)+iskip,3) /= 0.)then 
                                write(*,*)'Problem in init_hydro with zero or negative density'
                                call clean_stop
#endif
#if NDIM>2
                             else if(uold(ind_grid(i)+iskip,4) /= 0.)then 
                                write(*,*)'Problem in init_hydro with zero or negative density'
                                call clean_stop
#endif
                             end if
                             uold(ind_grid(i)+iskip,ivar)=xx(i)
                       end do
                    else
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*uold(ind_grid(i)+iskip,1)
                       end do
                    endif
                 end do
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




