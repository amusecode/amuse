!################################################################
!################################################################
!################################################################
!################################################################
subroutine rt_init_flow  
  use amr_commons
  use rt_hydro_commons, ONLY: nrtvar, rtuold
  implicit none

  integer::ilevel,ivar
  
  if(verbose)write(*,*)'Entering init_flow'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call rt_init_flow_fine(ilevel)
     call rt_upload_fine(ilevel)
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)
  end do
  if(verbose)write(*,*)'Complete init_flow'

end subroutine rt_init_flow
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rt_init_flow_fine(ilevel)
  use amr_commons
  use rt_hydro_commons
  use rt_cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  
  integer::i,icell,igrid,ncache,iskip,ngrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,info
  integer ,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy,vz,ek,ei,pp,xx1,xx2,xx3,dx_loc,scale,xval
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector),save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nrtvar),save::uu

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file1,ok_file2,ok_file3,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid

  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------
  filename=TRIM(initfile(ilevel))//'/ic_d'
  INQUIRE(file=filename,exist=ok_file1)
  if(multiple)then
     filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.00001'
     INQUIRE(file=filename,exist=ok_file2)
  else
     filename=TRIM(initfile(ilevel))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok_file2)
  endif
  ok_file = ok_file1 .or. ok_file2
  if(ok_file)then

     !-------------------------------------------------------------------------
     ! First step: compute level boundaries in terms of initial condition array
     !-------------------------------------------------------------------------
     if(ncache>0)then
     i1_min=n1(ilevel)+1; i1_max=0
     i2_min=n2(ilevel)+1; i2_max=0
     i3_min=n3(ilevel)+1; i3_max=0
     do ind=1,twotondim           
        do i=1,ncache
           igrid=active(ilevel)%igrid(i)
           xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
           xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
           xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
           xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
           xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
           xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
           i1_min=MIN(i1_min,int(xx1)+1)
           i1_max=MAX(i1_max,int(xx1)+1)
           i2_min=MIN(i2_min,int(xx2)+1)
           i2_max=MAX(i2_max,int(xx2)+1)
           i3_min=MIN(i3_min,int(xx3)+1)
           i3_max=MAX(i3_max,int(xx3)+1)
        end do
     end do
     error=.false.
     if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
     if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
     if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
     if(error) then
        write(*,*)'Some grid are outside initial conditions sub-volume'
        write(*,*)'for ilevel=',ilevel
        write(*,*)i1_min,i1_max
        write(*,*)i2_min,i2_max
        write(*,*)i3_min,i3_max
        write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
        call clean_stop
     end if
     endif

     !-----------------------------------------
     ! Second step: read initial condition file
     !-----------------------------------------
     ! Allocate initial conditions array
     if(ncache>0)allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
     allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
     ! Loop over input variables
     do ivar=1,nrtvar
        call title(ivar,ncharvar)
        filename=TRIM(initfile(ilevel))//'/ic_rt_'//TRIM(ncharvar)

        INQUIRE(file=filename,exist=ok_file3)
        if(ok_file3)then
           ! Reading the existing file   
           if(myid==1)write(*,*)'Reading file '//TRIM(filename)
           if(multiple)then
              ilun=ncpu+myid+10
              open(ilun,file=filename,form='unformatted')
              rewind ilun
              read(ilun) ! skip first line
              do i3=1,n3(ilevel)
                 read(ilun) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                 if(i3.ge.i3_min.and.i3.le.i3_max)then
                    init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                         & init_plane(i1_min:i1_max,i2_min:i2_max)
                 end if
              end do
              close(ilun)
           else
              if(myid==1)then
                 open(10,file=filename,form='unformatted')
                 rewind 10
                 read(10) ! skip first line
              endif
              do i3=1,n3(ilevel)
                 if(myid==1)then
                    read(10) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                 else
                    init_plane=0.0
                 endif
                 buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                 call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
                 if(ncache>0)then
                    if(i3.ge.i3_min.and.i3.le.i3_max)then
                       init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane(i1_min:i1_max,i2_min:i2_max)
                    end if
                 endif
              end do
              if(myid==1)close(10)
           endif
        else
           ! If file doesn't exist, initialize variable to default value 
           ! In most cases, this is zero (you can change that if necessary)
           if(myid==1)write(*,*)'File '//TRIM(filename)//' not found'
           if(myid==1)write(*,*)'Initialize corresponding variable to default value'
           init_array=0d0
        endif

        if(ncache>0)then

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ncache
              igrid=active(ilevel)%igrid(i)
              icell=igrid+iskip
              xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
              xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
              xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
              xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
              xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
              xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
              i1=int(xx1)+1
              i1=int(xx1)+1
              i2=int(xx2)+1
              i2=int(xx2)+1
              i3=int(xx3)+1
              i3=int(xx3)+1
              ! Scatter to corresponding primitive variable
              rtuold(icell,ivar)=init_array(i1,i2,i3)
           end do
        end do
        ! End loop over cells
        endif
     end do
     ! End loop over input variables

     ! Deallocate initial conditions array
     if(ncache>0)deallocate(init_array)
     deallocate(init_plane) 

  !-------------------------------------------------------
  ! Compute initial conditions from subroutine condinit
  !-------------------------------------------------------
  else

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

           ! Call initial condition routine
           call rt_condinit(xx,uu,dx_loc,ngrid)
           ! Scatter variables
           do ivar=1,nrtvar
              do i=1,ngrid
                 rtuold(ind_cell(i),ivar)=uu(i,ivar)
              end do
           end do
        end do
        ! End loop over cells
     end do
     ! End loop over grids

  end if

111 format('   Entering rt_init_flow_fine for level ',I2)

end subroutine rt_init_flow_fine
!################################################################
!################################################################
!################################################################
!################################################################

!************************************************************************
SUBROUTINE rt_region_condinit(x,uu,dx,nn)

! Initialize RT regions, as defined in the namelist setup file.
!
! x  =>  ncells*ndim: positions of grid cells
! uu <=  ncells*nrtvar: cell variables
! dx =>  real cell width
! nn =>  int number of cells
!------------------------------------------------------------------------
  use rt_parameters
  implicit none
  integer ::nn
  real(dp)::dx,dx_cgs
  real(dp),dimension(1:nvector,1:nrtvar)::uu
  real(dp),dimension(1:nvector,1:ndim)  ::x
  integer::i,k,group_ind
  real(dp)::vol,r,xn,yn,zn,en
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_np,scale_fp
!------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)
  dx_cgs=dx*scale_l
  ! Set some (tiny) default values in case n_region=0
  do i=1,nGroups ! Starting indices in uold and unew of each photon group
     uu(1:nn,iGroups(i))=smallNp
  end do

  ! Loop over RT regions
  do k=1,rt_nregion
     
     if (rt_n_region(k).le.0.0) rt_n_region(k)=smallnp
     group_ind = iGroups(rt_reg_group(k))
     if(rt_reg_group(k) .le. 0 .or. rt_reg_group(k) .gt. nGroups) cycle
     ! For "square" regions only:
     if(rt_region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=rt_exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-rt_reg_x_center(k))/rt_reg_length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-rt_reg_y_center(k))/rt_reg_length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-rt_reg_z_center(k))/rt_reg_length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(rt_exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region, inject value
           if(r .lt. 1.0)then
              uu(i,group_ind)=rt_n_region(k)
              uu(i,group_ind+1)=rt_u_region(k) * rt_c  
#if NDIM>1 
              uu(i,group_ind+2)=rt_v_region(k) * rt_c
#endif
#if NDIM>2
              uu(i,group_ind+3)=rt_w_region(k) * rt_c
#endif
           end if
        end do
     end if
     
     ! For "point" regions only:
     if(rt_region_type(k) .eq. 'point')then
        ! Volume element
        vol=dx_cgs**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-rt_reg_x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-rt_reg_y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-rt_reg_z_center(k))/dx,0.0_dp)
#endif
           r=xn*yn*zn
           if(r .gt. 0.) then
              ! If cell lies within CIC cloud, inject value
              ! Convert photon number to photon number density
              uu(i,group_ind) = rt_n_region(k)/scale_Np *r/vol 
              uu(i,group_ind+1) = rt_u_region(k)/scale_Np*r/vol*rt_c
#if NDIM>1
              uu(i,group_ind+2) = rt_v_region(k)/scale_Np*r/vol*rt_c
#endif
#if NDIM>2
              uu(i,group_ind+3) = rt_w_region(k)/scale_Np *r/vol*rt_c
#endif
           endif
        end do
     end if
  end do

  return
END SUBROUTINE rt_region_condinit




