!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_frame()
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  
  integer::dummy_io,info
  integer,parameter::tag=100

  character(len=5) :: istep_str
  character(len=100) :: moviedir, moviecmd, moviefile
  character(len=100) :: moviefile1,moviefile2,moviefile3
  
  integer::icell,ncache,iskip,ngrid,nlevelmax_frame
  integer::ilun,nx_loc,ipout,npout,npart_out,ind,ix,iy,iz
  integer::imin,imax,jmin,jmax,ii,jj
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp)::scale,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::xcen,ycen,zcen,delx,dely,delz
  real(dp)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame
  real(dp)::xleft,xright,yleft,yright,zleft,zright
  real(dp)::xxleft,xxright,yyleft,yyright,zzleft,zzright
  real(dp)::dx_frame,dy_frame,dx,dx_loc
  real(dp)::dx_cell,dy_cell,dz_cell,dvol
  real(kind=8)::cell_value
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  logical,dimension(1:nvector)::ok
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(kind=8),dimension(:,:,:),allocatable::data_frame,data_frame_all
  real(kind=4),dimension(:,:),allocatable::data_single
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: observer(3),thetay,thetaz,theta,phi,temp,ekk
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel
  integer::i,ig,ip,npart1
  integer::nalloc1,nalloc2

  integer,dimension(1:nvector),save::ind_part
  logical::opened
  opened=.false.
  
#if NDIM > 1

  ! Update counter
  imov=imov+1
  if(imov>imovout)return

  ! Determine the filename, dir, etc
  if(myid==1)write(*,*)'Computing and dumping movie frame'

  ! Determine the filename, dir, etc
  if(myid==1)write(*,*)'Computing and dumping movie frame'

  call title(imov, istep_str)
  moviedir = 'movie/'
  moviecmd = 'mkdir -p '//trim(moviedir)
  if(myid==1) write(*,*) "Writing frame ", istep_str
#ifdef NOSYSTEM
  if(myid==1)call PXFMKDIR(TRIM(moviedir),LEN(TRIM(moviedir)),O'755',info)  
#else
  if(myid==1)call system(moviecmd)
#endif

  moviefile = trim(moviedir)//'info_'//trim(istep_str)//'.txt'
  if(myid==1)call output_info(moviefile)

  moviefile1 = trim(moviedir)//'dens_'//trim(istep_str)//'.map'
  moviefile2 = trim(moviedir)//'temp_'//trim(istep_str)//'.map'
  moviefile3 = trim(moviedir)//'metal_'//trim(istep_str)//'.map'

  if(levelmax_frame==0)then
     nlevelmax_frame=nlevelmax
  else
     nlevelmax_frame=levelmax_frame
  endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Compute frame boundaries
  xcen=xcentre_frame(1)+xcentre_frame(2)*aexp+xcentre_frame(3)*aexp**2+xcentre_frame(4)*aexp**3
  ycen=ycentre_frame(1)+ycentre_frame(2)*aexp+ycentre_frame(3)*aexp**2+ycentre_frame(4)*aexp**3
  zcen=zcentre_frame(1)+zcentre_frame(2)*aexp+zcentre_frame(3)*aexp**2+zcentre_frame(4)*aexp**3
  delx=deltax_frame(1)+deltax_frame(2)/aexp !+deltax_frame(3)*aexp**2+deltax_frame(4)*aexp**3  !Essentially comoving or physical
  dely=deltay_frame(1)+deltay_frame(2)/aexp !+deltay_frame(3)*aexp**2+deltay_frame(4)*aexp**3
  delz=deltaz_frame(1)+deltaz_frame(2)/aexp !+deltaz_frame(3)*aexp**2+deltaz_frame(4)*aexp**3
  xleft_frame=xcen-delx/2.
  xright_frame=xcen+delx/2.
  yleft_frame=ycen-dely/2.
  yright_frame=ycen+dely/2.
  zleft_frame=zcen-delz/2.
  zright_frame=zcen+delz/2.
  
  ! Allocate image
  allocate(data_frame(1:nx_frame,1:ny_frame,1:4))
  data_frame=0d0
  dx_frame=delx/dble(nx_frame)
  dy_frame=dely/dble(ny_frame)

  ! Loop over levels
  do ilevel=levelmin,nlevelmax_frame

     ! Mesh size at level ilevel in coarse cell units
     dx=0.5D0**ilevel
     
     ! Set position of cell centres relative to grid centre
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do
  
     dx_loc=dx*scale
     ncache=active(ilevel)%ngrid

     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather cell centre positions
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do
           
           ! Check if cell is to be considered
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0.or.ilevel==nlevelmax_frame
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Check if the cell intersect the domain
                 xleft=xx(i,1)-dx_loc/2.
                 xright=xx(i,1)+dx_loc/2.
                 yleft=xx(i,2)-dx_loc/2.
                 yright=xx(i,2)+dx_loc/2.
#if NDIM>2                 
                 zleft=xx(i,3)-dx_loc/2.
                 zright=xx(i,3)+dx_loc/2.
                 if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                      & yright.lt.yleft_frame.or.yleft.ge.yright_frame.or.&
                      & zright.lt.zleft_frame.or.zleft.ge.zright_frame)cycle
#else
                 if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                      & yright.lt.yleft_frame.or.yleft.ge.yright_frame)cycle
#endif
                 ! Compute map indices for the cell
                 if(xleft>xleft_frame)then
                    imin=min(int((xleft-xleft_frame)/dx_frame)+1,nx_frame)
                 else
                    imin=1
                 endif
                 imax=min(int((xright-xleft_frame)/dx_frame)+1,nx_frame)
                 if(yleft>yleft_frame)then
                    jmin=min(int((yleft-yleft_frame)/dy_frame)+1,ny_frame)
                 else
                    jmin=1
                 endif
                 jmax=min(int((yright-yleft_frame)/dy_frame)+1,ny_frame)
                 
                 ! Fill up map with projected mass
#if NDIM>2                 
                 dz_cell=min(zright_frame,zright)-max(zleft_frame,zleft)
#endif
                 do ii=imin,imax
                    xxleft=xleft_frame+dble(ii-1)*dx_frame
                    xxright=xxleft+dx_frame
                    dx_cell=min(xxright,xright)-max(xxleft,xleft)
                    do jj=jmin,jmax
                       yyleft=yleft_frame+dble(jj-1)*dy_frame
                       yyright=yyleft+dy_frame
                       dy_cell=min(yyright,yright)-max(yyleft,yleft)
                       ! Intersection volume
                       dvol=dx_cell*dy_cell
#if NDIM>2                 
                       dvol=dvol*dz_cell
#endif
                       data_frame(ii,jj,1)=data_frame(ii,jj,1)+dvol*uold(ind_cell(i),1)
                       data_frame(ii,jj,2)=data_frame(ii,jj,2)+dvol*uold(ind_cell(i),1)**2

                       !Get temperature
                       ekk=0.0d0
                       do idim=1,3
                          ekk=ekk+0.5*uold(ind_cell(i),idim+1)**2/uold(ind_cell(i),1)
                       enddo
                       temp=(gamma-1.0)*(uold(ind_cell(i),5)-ekk) !pressure
                       temp=temp/uold(ind_cell(i),1)*scale_T2 !temperature in K

                       data_frame(ii,jj,3)=data_frame(ii,jj,3)+dvol*uold(ind_cell(i),1)*temp !mass weighted temperature

                       if(metal)then
                       data_frame(ii,jj,4)=data_frame(ii,jj,4)+dvol*uold(ind_cell(i),6)
                       endif
                    end do
                 end do
              end if
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over levels

  ! Convert into mass weighted
!  do ii=1,nx_frame
!     do jj=1,ny_frame
!        data_frame(ii,jj,2)=data_frame(ii,jj,2)/data_frame(ii,jj,1)
!        data_frame(ii,jj,3)=data_frame(ii,jj,3)/data_frame(ii,jj,1)
!        if(metal)then
!        data_frame(ii,jj,4)=data_frame(ii,jj,4)/data_frame(ii,jj,1)
!        endif
!     end do
!  end do
#ifndef WITHOUTMPI
  allocate(data_frame_all(1:nx_frame,1:ny_frame,1:4))
  call MPI_ALLREDUCE(data_frame,data_frame_all,nx_frame*ny_frame*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  data_frame=data_frame_all
  deallocate(data_frame_all)
#endif
  ! Convert into mass weighted                                                                                                         
  do ii=1,nx_frame
     do jj=1,ny_frame
        data_frame(ii,jj,2)=data_frame(ii,jj,2)/data_frame(ii,jj,1)
        data_frame(ii,jj,3)=data_frame(ii,jj,3)/data_frame(ii,jj,1)
        if(metal)then
        data_frame(ii,jj,4)=data_frame(ii,jj,4)/data_frame(ii,jj,1)
        endif
     end do
  end do

!     write(*,*) 'testing1', data_frame(100,100,1),data_frame(100,100,2)

  if(myid==1)then
     ilun=10
     allocate(data_single(1:nx_frame,1:ny_frame))
     ! Output mass weighted density
     open(ilun,file=TRIM(moviefile1),form='unformatted')
     data_single=data_frame(:,:,2)
     rewind(ilun)  
     if(tendmov>0)then
        write(ilun)t,delx,dely,delz
     else
        write(ilun)aexp,delx,dely,delz
     endif
     write(ilun)nx_frame,ny_frame
     write(ilun)data_single
     close(ilun)
     ! Output mass weighted temperature
     open(ilun,file=TRIM(moviefile2),form='unformatted')
     data_single=data_frame(:,:,3)
!     write(*,*) 'testing', data_single(100,100)
     rewind(ilun)  
     if(tendmov>0)then
        write(ilun)t,delx,dely,delz
     else
        write(ilun)aexp,delx,dely,delz
     endif
     write(ilun)nx_frame,ny_frame
     write(ilun)data_single
     close(ilun)
     ! Output mass weighted metal fraction
     if(metal)then
        open(ilun,file=TRIM(moviefile3),form='unformatted')
        data_single=data_frame(:,:,4)
        rewind(ilun)  
        if(tendmov>0)then
           write(ilun)t,delx,dely,delz
        else
           write(ilun)aexp,delx,dely,delz
        endif
        write(ilun)nx_frame,ny_frame
        write(ilun)data_single
        close(ilun)
     endif
     deallocate(data_single)
  endif

  deallocate(data_frame)
#endif

end subroutine output_frame



