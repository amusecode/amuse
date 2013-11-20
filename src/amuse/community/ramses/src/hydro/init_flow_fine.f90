!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow  
  use amr_commons
  use hydro_commons, ONLY: nvar, uold
  implicit none

  integer::ilevel,ivar
  
  if(verbose)write(*,*)'Entering init_flow'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call init_flow_fine(ilevel)
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  end do
  if(verbose)write(*,*)'Complete init_flow'

end subroutine init_flow
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  
  integer::i,icell,igrid,ncache,iskip,ngrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,info,nvar_in
  integer ,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy,vz,ek,ei,pp,xx1,xx2,xx3,dx_loc,scale,xval
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ,save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu

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
     do ivar=1,nvar
        if(cosmo)then
           ! Read baryons initial overdensity and displacement at a=aexp
           if(multiple)then
              call title(myid,nchar)
              if(ivar==1)filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
              if(ivar==2)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
              if(ivar==3)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
              if(ivar==4)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
              if(ivar==5)filename=TRIM(initfile(ilevel))//'/dir_tempb/ic_tempb.'//TRIM(nchar)
           else
              if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_deltab'
              if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_velcx'
              if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_velcy'
              if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_velcz'
              if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_tempb'
           endif
        else
           ! Read primitive variables
           if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_d'
           if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_u'
           if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_v'
           if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_w'
           if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_p'
        endif
        call title(ivar,ncharvar)
        if(ivar>5)then
           call title(ivar-5,ncharvar)
           filename=TRIM(initfile(ilevel))//'/ic_pvar_'//TRIM(ncharvar)
        endif

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
                 if(ncache>0)then
                    if(i3.ge.i3_min.and.i3.le.i3_max)then
                       init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane(i1_min:i1_max,i2_min:i2_max)
                    end if
                 endif
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
           if(ncache>0)then
              init_array=0d0
              ! Default value for metals
              if(cosmo.and.ivar==imetal.and.metal)init_array=z_ave*0.02 ! from solar units
              ! Default value for ionization fraction
              xval=sqrt(omega_m)/(h0/100.*omega_b) ! From the book of Peebles p. 173
              if(cosmo.and.ivar==ixion.and.aton)init_array=1.2d-5*xval
           endif
        endif

        if(ncache>0)then

        ! For cosmo runs, rescale initial conditions to code units
        if(cosmo)then
           ! Compute approximate average temperature in K
           if(.not. cooling)T2_start=1.356d-2/aexp**2
           if(ivar==1)init_array=(1.0+dfact(ilevel)*init_array)*omega_b/omega_m
           if(ivar==2)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==3)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==4)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==ndim+2)init_array=(1.0+init_array)*T2_start/scale_T2
        endif

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
              uold(icell,ivar)=init_array(i1,i2,i3)
           end do
        end do
        ! End loop over cells
        endif
     end do
     ! End loop over input variables

     ! Deallocate initial conditions array
     if(ncache>0)deallocate(init_array)
     deallocate(init_plane) 

     !----------------------------------------------------------------
     ! For cosmology runs: compute pressure, prevent negative density
     !----------------------------------------------------------------
     if(cosmo)then
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
              ! Prevent negative density
              do i=1,ngrid
                 rr=max(uold(ind_cell(i),1),0.1*omega_b/omega_m)
                 uold(ind_cell(i),1)=rr
              end do
              ! Compute pressure from temperature and density
              do i=1,ngrid
                 uold(ind_cell(i),ndim+2)=uold(ind_cell(i),1)*uold(ind_cell(i),ndim+2)
              end do
           end do
           ! End loop over cells
        end do
        ! End loop over grids
     end if

     !---------------------------------------------------
     ! Third step: compute initial conservative variables
     !---------------------------------------------------
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        vy=0.0
        vz=0.0
        ! Loop over cells
        do ind=1,twotondim
           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Compute total energy density
           do i=1,ngrid
              rr=uold(ind_cell(i),1)
              vx=uold(ind_cell(i),2)
#if NDIM>1
              vy=uold(ind_cell(i),3)
#endif
#if NDIM>2
              vz=uold(ind_cell(i),4)
#endif
              pp=uold(ind_cell(i),ndim+2)
              ek=0.5d0*(vx**2+vy**2+vz**2)
              ei=pp/(gamma-1.0)
              vv(i)=ei+rr*ek
           end do
           ! Scatter to corresponding conservative variable
           do i=1,ngrid
              uold(ind_cell(i),ndim+2)=vv(i)
           end do
           ! Compute momentum density
           do ivar=1,ndim
              do i=1,ngrid
                 rr=uold(ind_cell(i),1)
                 vx=uold(ind_cell(i),ivar+1)
                 vv(i)=rr*vx
              end do
              ! Scatter to corresponding conservative variable
              do i=1,ngrid
                 uold(ind_cell(i),ivar+1)=vv(i)
              end do
           end do
#if NVAR > NDIM + 2
           ! Compute passive variable density
           do ivar=ndim+3,nvar
              do i=1,ngrid
                 rr=uold(ind_cell(i),1)
                 uold(ind_cell(i),ivar)=rr*uold(ind_cell(i),ivar)
              end do
           enddo
#endif
        end do
        ! End loop over cells
        
     end do
     ! End loop over grids

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
           call condinit(xx,uu,dx_loc,ngrid)
           ! Scatter variables
           do ivar=1,nvar
              do i=1,ngrid
                 uold(ind_cell(i),ivar)=uu(i,ivar)
              end do
           end do
        end do
        ! End loop over cells
     end do
     ! End loop over grids

  end if
  
111 format('   Entering init_flow_fine for level ',I2)

end subroutine init_flow_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine region_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar)::q
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,ivar,k
  real(dp)::vol,r,xn,yn,zn,en

  ! Set some (tiny) default values in case n_region=0
  q(1:nn,1)=smallr
  q(1:nn,2)=0.0d0
#if NDIM>1
  q(1:nn,3)=0.0d0
#endif
#if NDIM>2
  q(1:nn,4)=0.0d0
#endif
  q(1:nn,ndim+2)=smallr*smallc**2/gamma
#if NVAR > NDIM + 2
  do ivar=ndim+3,nvar
     q(1:nn,ivar)=0.0d0
  end do
#endif

  ! Loop over initial conditions regions
  do k=1,nregion
     
     ! For "square" regions only:
     if(region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-x_center(k))/length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-y_center(k))/length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-z_center(k))/length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region,
           ! REPLACE primitive variables by region values
           if(r<1.0)then
              q(i,1)=d_region(k)
              q(i,2)=u_region(k)
#if NDIM>1
              q(i,3)=v_region(k)
#endif
#if NDIM>2
              q(i,4)=w_region(k)
#endif
              q(i,ndim+2)=p_region(k)
#if NVAR>NDIM+2
              do ivar=ndim+3,nvar
                 q(i,ivar)=var_region(k,ivar-ndim-2)
              end do
#endif
           end if
        end do
     end if
     
     ! For "point" regions only:
     if(region_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-z_center(k))/dx,0.0_dp)
#endif
           r=xn*yn*zn
           ! If cell lies within CIC cloud, 
           ! ADD to primitive variables the region values
           q(i,1)=q(i,1)+d_region(k)*r/vol
           q(i,2)=q(i,2)+u_region(k)*r
#if NDIM>1
           q(i,3)=q(i,3)+v_region(k)*r
#endif
#if NDIM>2
           q(i,4)=q(i,4)+w_region(k)*r
#endif
           q(i,ndim+2)=q(i,ndim+2)+p_region(k)*r/vol
#if NVAR>NDIM+2
           do ivar=ndim+3,nvar
              q(i,ivar)=var_region(k,ivar-ndim-2)
           end do
#endif
        end do
     end if
  end do

  return
end subroutine region_condinit
