!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  integer::ilevel,icycle,nsubdiff,ivar,i,iskip,ind,info
  real(dp)::dx,scale,dx_loc,dtdiff,norm
  real(dp),dimension(1:3)::skip_loc

  ! Determine minimum mesh size
  do ilevel=levelmin,nlevelmax
     if(numbtot(1,ilevel)>0)dx=0.5D0**ilevel
  end do
  
  ! Rescaling factors
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(icoarse_max-icoarse_min+1)/boxlen
  dx_loc=dx/scale

  dtdiff=0.05*dx_loc**2/eta_mag
  nsubdiff=dtnew(levelmin)/dtdiff
  nsubdiff=nsubdiff+1
  dtdiff=-dtnew(levelmin)/dble(nsubdiff)*eta_mag

  do icycle=1,nsubdiff
     if(myid==1)write(*,*)icycle,nsubdiff,dtdiff
     do ilevel=levelmin,nlevelmax
        call set_unew(ilevel)
     end do
     do ilevel=nlevelmax,levelmin,-1
        call diffusion_fine(ilevel,dtdiff)
#ifndef WITHOUTMPI
        call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        call make_virtual_reverse_dp(unew(1,6),ilevel)
        call make_virtual_reverse_dp(unew(1,7),ilevel)
        call make_virtual_reverse_dp(unew(1,8),ilevel)
        call make_virtual_reverse_dp(unew(1,nvar+1),ilevel)
        call make_virtual_reverse_dp(unew(1,nvar+2),ilevel)
        call make_virtual_reverse_dp(unew(1,nvar+3),ilevel)
        call set_uold(ilevel)
        call upload_fine(ilevel)
        call make_virtual_fine_dp(uold(1,6),ilevel)
        call make_virtual_fine_dp(uold(1,7),ilevel)
        call make_virtual_fine_dp(uold(1,8),ilevel)
        call make_virtual_fine_dp(uold(1,nvar+1),ilevel)
        call make_virtual_fine_dp(uold(1,nvar+2),ilevel)
        call make_virtual_fine_dp(uold(1,nvar+3),ilevel)
     end do
  end do

end subroutine diffusion
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion_fine(ilevel,dtdiff)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dtdiff

  integer::i,ivar,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call diffine1(ind_grid,ngrid,dtdiff,ilevel)
  end do

111 format('   Entering godunov_fine for level ',i2)

end subroutine diffusion_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffine1(ind_grid,ncache,dtdiff,ilevel)
  use amr_commons
  use hydro_commons
 implicit none
  integer::ilevel,ncache
  real(dp)::dtdiff
  integer,dimension(1:ncache)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first MHD variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It then computes
  ! the current at cell edges. Finally, currents are corrected from finer level
  ! and boundary currents are stored in buffer regions. Updated 
  ! conservative variables are stored in array unew(:).
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:6   ),save::B1
  integer ,dimension(1:nvector,0:twondim)         ,save::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:6   ),save::B2
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::v2
  real(dp),dimension(1:nvector,1:ndim),save::vv,xx

  logical ,dimension(1:nvector,-1:4,-1:4,-1:4),save::ok
  logical ,dimension(1:nvector,-1:4,-1:4,-1:4),save::buffer
  real(dp),dimension(1:nvector, 0:4,-1:4,-1:4),save::Bx
  real(dp),dimension(1:nvector,-1:4, 0:4,-1:4),save::By
  real(dp),dimension(1:nvector,-1:4,-1:4, 0:4),save::Bz
  real(dp),dimension(1:nvector, 1:2, 1:3, 1:3),save::emfx
  real(dp),dimension(1:nvector, 1:3, 1:2, 1:3),save::emfy
  real(dp),dimension(1:nvector, 1:3, 1:3, 1:2),save::emfz
  real(dp),dimension(1:nvector),save :: dB

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,igrid
  logical,dimension(1:nvector),save::exist_nbor

  real(dp),dimension(1:3)::skip_loc
  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::ind,ix,iy,iz
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::ind_father1,ind_father2,ind_father3
  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::interpol_type_old,ivar1,ivar2,ivar3,ivar4,ivar5,ivar6
  real(dp)::dx,dflux,weight,dflux_x,dflux_y,dflux_z,scale,dx_loc

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel  

  ! Rescaling factors
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(icoarse_max-icoarse_min+1)/boxlen
  dx_loc=dx/scale
  ivar1=6; ivar2=7; ivar3=8
  ivar4=nvar+1; ivar5=nvar+2; ivar6=nvar+3

  ! Gather 3^ndim neighboring father cells
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=0,2
  do j1=0,2
  do i1=0,2
     
     ! Check if neighboring grid exists
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        exist_nbor(i)=igrid_nbor(i)>0
     end do
     
     ! If not, interpolate variables from parent cells
     nbuffer=0
     do i=1,ncache
        if(.not. exist_nbor(i))then
           nbuffer=nbuffer+1
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
           igrid     (nbuffer)=son(nbors_father_cells(i,ind_father))
        end if
     end do

     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do i=1,nbuffer
              B1(i,j,1)=uold(ibuffer_father(i,j),ivar1)
              B1(i,j,2)=uold(ibuffer_father(i,j),ivar2)
              B1(i,j,3)=uold(ibuffer_father(i,j),ivar3)
              B1(i,j,4)=uold(ibuffer_father(i,j),ivar4)
              B1(i,j,5)=uold(ibuffer_father(i,j),ivar5)
              B1(i,j,6)=uold(ibuffer_father(i,j),ivar6)
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_mag(B1,ind1,B2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=0,1
     do j2=0,1
     do i2=0,1

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+igrid_nbor(i)
        end do
        
        i3=1+2*(i1-1)+i2
        j3=1+2*(j1-1)+j2
        k3=1+2*(k1-1)+k2
        
        ! Gather MHD variables
        if(i3>=0)then
           ibuffer=0
           do i=1,ncache
              if(exist_nbor(i))then
                 Bx(i,i3,j3,k3)=uold(ind_cell(i),ivar1)
              else
                 ibuffer=ibuffer+1
                 Bx(i,i3,j3,k3)=B2(ibuffer,ind_son,1)
              end if
           end do
        endif
        if(j3>=0)then
           ibuffer=0
           do i=1,ncache
              if(exist_nbor(i))then
                 By(i,i3,j3,k3)=uold(ind_cell(i),ivar2)
              else
                 ibuffer=ibuffer+1
                 By(i,i3,j3,k3)=B2(ibuffer,ind_son,2)
              end if
           end do
        endif
        if(k3>=0)then
           ibuffer=0
           do i=1,ncache
              if(exist_nbor(i))then
                 Bz(i,i3,j3,k3)=uold(ind_cell(i),ivar3)
              else
                 ibuffer=ibuffer+1
                 Bz(i,i3,j3,k3)=B2(ibuffer,ind_son,3)
              end if
           end do
        endif
        
        ! Gather refinement flag
        do i=1,ncache
           if(exist_nbor(i))then
              ok(i,i3,j3,k3)=son(ind_cell(i))>0
              buffer(i,i3,j3,k3)=.false.
           else
              ok(i,i3,j3,k3)=.false.
              buffer(i,i3,j3,k3)=.true.
           end if
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !----------------
  ! Compute current
  !----------------
  emfx=0.0d0; emfy=0.0d0; emfz=0.0d0
  call cmp_current(Bx,By,Bz,emfx,emfy,emfz,buffer,2,2,2,ncache,dx_loc,dx_loc,dx_loc)

  !-------------------------------------------------
  ! Reset current along direction x at refined edges
  !-------------------------------------------------
  do k3=1,3
  do j3=1,3
  do i3=1,2
     do i=1,ncache
        if(ok(i,i3,j3  ,k3  ) .or. ok(i,i3,j3  ,k3-1) .or.  &
         & ok(i,i3,j3-1,k3  ) .or. ok(i,i3,j3-1,k3-1))then
           emfx(i,i3,j3,k3)=0.0d0
        end if
     end do
  end do
  end do
  end do
  !-------------------------------------------------
  ! Reset current along direction y at refined edges
  !-------------------------------------------------
  do k3=1,3
  do j3=1,2
  do i3=1,3
     do i=1,ncache
        if(ok(i,i3  ,j3,k3  ) .or. ok(i,i3  ,j3,k3-1) .or.  &
         & ok(i,i3-1,j3,k3  ) .or. ok(i,i3-1,j3,k3-1))then
           emfy(i,i3,j3,k3)=0.0d0
        end if
     end do
  end do
  end do
  end do
  !-------------------------------------------------
  ! Reset current along direction z at refined edges
  !-------------------------------------------------
  do k3=1,2
  do j3=1,3
  do i3=1,3
     do i=1,ncache
        if(ok(i,i3  ,j3  ,k3) .or. ok(i,i3  ,j3-1,k3) .or.  &
         & ok(i,i3-1,j3  ,k3) .or. ok(i,i3-1,j3-1,k3))then
           emfz(i,i3,j3,k3)=0.0d0
        end if
     end do
  end do
  end do
  end do

  !------------------------------------
  ! Conservative update at level ilevel
  !------------------------------------
  do k3=1,2
  do j3=1,2
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update Bx using constraint transport
     do i=1,ncache
        dflux_x=( emfy(i,i3,j3,k3)-emfy(i,i3,j3,k3+1) ) &
          &    -( emfz(i,i3,j3,k3)-emfz(i,i3,j3+1,k3) )
        unew(ind_cell(i),ivar1)=unew(ind_cell(i),ivar1)+dflux_x*dtdiff/dx_loc
        dflux_x=( emfy(i,i3+1,j3,k3)-emfy(i,i3+1,j3,k3+1) ) &
          &    -( emfz(i,i3+1,j3,k3)-emfz(i,i3+1,j3+1,k3) )  
        unew(ind_cell(i),ivar4)=unew(ind_cell(i),ivar4)+dflux_x*dtdiff/dx_loc
     end do
     ! Update By using constraint transport
     do i=1,ncache
        dflux_y=( emfz(i,i3,j3,k3)-emfz(i,i3+1,j3,k3) ) &
          &    -( emfx(i,i3,j3,k3)-emfx(i,i3,j3,k3+1) )
        unew(ind_cell(i),ivar2)=unew(ind_cell(i),ivar2)+dflux_y*dtdiff/dx_loc
        dflux_y=( emfz(i,i3,j3+1,k3)-emfz(i,i3+1,j3+1,k3) ) &
          &    -( emfx(i,i3,j3+1,k3)-emfx(i,i3,j3+1,k3+1) )
        unew(ind_cell(i),ivar5)=unew(ind_cell(i),ivar5)+dflux_y*dtdiff/dx_loc
     end do
     ! Update Bz using constraint transport
     do i=1,ncache
        dflux_z=( emfx(i,i3,j3,k3)-emfx(i,i3,j3+1,k3) ) &
          &    -( emfy(i,i3,j3,k3)-emfy(i,i3+1,j3,k3) )
        unew(ind_cell(i),ivar3)=unew(ind_cell(i),ivar3)+dflux_z*dtdiff/dx_loc
        dflux_z=( emfx(i,i3,j3,k3+1)-emfx(i,i3,j3+1,k3+1) ) &
          &    -( emfy(i,i3,j3,k3+1)-emfy(i,i3+1,j3,k3+1) )
        unew(ind_cell(i),ivar6)=unew(ind_cell(i),ivar6)+dflux_z*dtdiff/dx_loc
     end do
  end do
  end do
  end do

  if(ilevel>levelmin)then

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  i1=1; j1=1; k1=1

  !--------------------------------------
  ! Deal with 4 EMFx edges
  !--------------------------------------

  ! Update coarse By and Bz using fine EMFx on Y=0 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
  ind_father2=1+(i1  )+3*(j1-1)+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,1)+emfx(i,2,1,1))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar2)=unew(ind_buffer1,ivar2)+dflux
     unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)+dflux
     unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)-dflux
     unew(ind_buffer3,ivar3)=unew(ind_buffer3,ivar3)-dflux
  end do

  ! Update coarse By and Bz using fine EMFx on Y=0 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
  ind_father2=1+(i1  )+3*(j1-1)+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,3)+emfx(i,2,1,3))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar6)=unew(ind_buffer1,ivar6)-dflux
     unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)-dflux
     unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)-dflux
     unew(ind_buffer3,ivar2)=unew(ind_buffer3,ivar2)-dflux
  end do

  ! Update coarse By and Bz using fine EMFx on Y=1 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
  ind_father2=1+(i1  )+3*(j1+1)+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,3)+emfx(i,2,3,3))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar5)=unew(ind_buffer1,ivar5)-dflux
     unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)-dflux
     unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)+dflux
     unew(ind_buffer3,ivar6)=unew(ind_buffer3,ivar6)+dflux
  end do

  ! Update coarse By and Bz using fine EMFx on Y=1 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
  ind_father2=1+(i1  )+3*(j1+1)+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,1)+emfx(i,2,3,1))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar3)=unew(ind_buffer1,ivar3)+dflux
     unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)+dflux
     unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)+dflux
     unew(ind_buffer3,ivar5)=unew(ind_buffer3,ivar5)+dflux
  end do

  !--------------------------------------
  ! Deal with 4 EMFy edges
  !--------------------------------------

  ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
  ind_father2=1+(i1-1)+3*(j1  )+9*(k1-1)
  ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,1)+emfy(i,1,2,1))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar1)=unew(ind_buffer1,ivar1)-dflux
     unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)-dflux
     unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)+dflux
     unew(ind_buffer3,ivar3)=unew(ind_buffer3,ivar3)+dflux
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=1 grid edge
  ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1  )+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,3)+emfy(i,1,2,3))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar6)=unew(ind_buffer1,ivar6)+dflux
     unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)+dflux
     unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)+dflux
     unew(ind_buffer3,ivar1)=unew(ind_buffer3,ivar1)+dflux
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
  ind_father2=1+(i1+1)+3*(j1  )+9*(k1+1)
  ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,3)+emfy(i,3,2,3))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar4)=unew(ind_buffer1,ivar4)+dflux
     unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)+dflux
     unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)-dflux
     unew(ind_buffer3,ivar6)=unew(ind_buffer3,ivar6)-dflux
  end do

  ! Update coarse Bx and Bz using fine EMFx on X=1 and Z=0 grid edge
  ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1  )+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,1)+emfy(i,3,2,1))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar3)=unew(ind_buffer1,ivar3)-dflux
     unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)-dflux
     unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)-dflux
     unew(ind_buffer3,ivar4)=unew(ind_buffer3,ivar4)-dflux
  end do

  !--------------------------------------
  ! Deal with 4 EMFz edges
  !--------------------------------------

  ! Update coarse Bx and By using fine EMFz on X=0 and Y=0 grid edge
  ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1-1)+9*(k1  )
  ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,1,1)+emfz(i,1,1,2))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar1)=unew(ind_buffer1,ivar1)+dflux
     unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)+dflux
     unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)-dflux
     unew(ind_buffer3,ivar2)=unew(ind_buffer3,ivar2)-dflux
  end do

  ! Update coarse Bx and By using fine EMFz on X=0 and Y=1 grid edge
  ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1+1)+9*(k1  )
  ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,3,1)+emfz(i,1,3,2))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar5)=unew(ind_buffer1,ivar5)-dflux
     unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)-dflux
     unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)-dflux
     unew(ind_buffer3,ivar1)=unew(ind_buffer3,ivar1)-dflux
  end do

  ! Update coarse Bx and By using fine EMFz on X=1 and Y=1 grid edge
  ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1+1)+9*(k1  )
  ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,3,1)+emfz(i,3,3,2))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar4)=unew(ind_buffer1,ivar4)-dflux
     unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)-dflux
     unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)+dflux
     unew(ind_buffer3,ivar5)=unew(ind_buffer3,ivar5)+dflux
  end do

  ! Update coarse Bx and By using fine EMFz on X=1 and Y=0 grid edge
  ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1-1)+9*(k1  )
  ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,1,1)+emfz(i,3,1,2))*0.25*weight*dtdiff/dx_loc
     unew(ind_buffer1,ivar2)=unew(ind_buffer1,ivar2)+dflux
     unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)+dflux
     unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)+dflux
     unew(ind_buffer3,ivar4)=unew(ind_buffer3,ivar4)+dflux
  end do

  endif

end subroutine diffine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_current(Bx,By,Bz,Ex_arete,Ey_arete,Ez_arete,buffer, &
     & Nx,Ny,Nz,ngrid,dx,dy,dz)
  use amr_parameters,ONLY:dp,verbose,nvector
  use hydro_parameters,ONLY:slope_type
  implicit none
  integer :: Nx,Ny,Nz,ngrid
  real(dp),dimension(1:nvector, 0:Nx+2,-1:Ny+2,-1:Nz+2) :: Bx
  real(dp),dimension(1:nvector,-1:Nx+2, 0:Ny+2,-1:Nz+2) :: By
  real(dp),dimension(1:nvector,-1:Nx+2,-1:Ny+2, 0:Nz+2) :: Bz
  logical ,dimension(1:nvector,-1:Nx+2,-1:Ny+2,-1:Nz+2) :: buffer
  real(dp)::dx,dy,dz
  real(dp),dimension(1:nvector, 1:Nx  , 1:Ny+1, 1:Nz+1) :: Ex_arete
  real(dp),dimension(1:nvector, 1:Nx+1, 1:Ny  , 1:Nz+1) :: Ey_arete
  real(dp),dimension(1:nvector, 1:Nx+1, 1:Ny+1, 1:Nz  ) :: Ez_arete
  !
  real(dp) :: dBx_arete_dy,dBx_arete_dz
  real(dp) :: dBy_arete_dx,dBy_arete_dz
  real(dp) :: dBz_arete_dx,dBz_arete_dy
  real(dp) :: dx_L,dx_R,dy_L,dy_R,dz_L,dz_R
  integer  :: ic,i,j,k,im1,jm1,km1
  integer  :: Nxp1,Nyp1,Nzp1

  Nxp1=Nx+1
  Nyp1=Ny+1
  Nzp1=Nz+1

  ! Aretes paralleles a l'axe des x
  do k=1,Nzp1
     km1=k-1
     do j=1,Nyp1
        jm1=j-1
        do i=1,Nx
           do ic=1,ngrid
              dBz_arete_dy=(Bz(ic,i,j,k)-Bz(ic,i,jm1,k))
              dBy_arete_dz=(By(ic,i,j,k)-By(ic,i,j,km1))
              Ex_arete(ic,i,j,k)=(dBz_arete_dy-dBy_arete_dz)/dx
           enddo
        enddo
     enddo
  enddo

  ! Aretes paralleles a l'axe des y
  do k=1,Nzp1
     km1=k-1
     do j=1,Ny
        do i=1,Nxp1
           im1=i-1
           do ic=1,ngrid
              dBx_arete_dz=(Bx(ic,i,j,k)-Bx(ic,i,j,km1))
              dBz_arete_dx=(Bz(ic,i,j,k)-Bz(ic,im1,j,k))
              Ey_arete(ic,i,j,k)=(dBx_arete_dz-dBz_arete_dx)/dx
           enddo
        enddo
     enddo
  enddo

  ! Aretes paralleles a l'axe des z
  do k=1,Nz
     do j=1,Nyp1
        jm1=j-1
        do i=1,Nxp1
           im1=i-1             
           do ic=1,ngrid
              dBy_arete_dx=(By(ic,i,j,k)-By(ic,im1,j,k))
              dBx_arete_dy=(Bx(ic,i,j,k)-Bx(ic,i,jm1,k))
              Ez_arete(ic,i,j,k)=(dBy_arete_dx-dBx_arete_dy)/dx
           enddo
        enddo
     enddo
  enddo

end subroutine cmp_current
!###########################################################
!###########################################################
!###########################################################
!###########################################################
