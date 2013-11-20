!###########################################################
!########################################################### 
!###########################################################
!###########################################################
subroutine upload_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the hydro variables.
  !----------------------------------------------------------------------
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_split
  integer,dimension(1:nvector),save::ind_unsplit,igrid_son
  integer ,dimension(1:nvector,0:twondim),save::igridn

  integer,dimension(1:3,1:2,1:8)::iii,jjj
  integer::ind_left,ind_right,neul=5
  integer::id1,id2,ig1,ig2,ih1,ih2
  integer::i,icpu,idim,ivar,ncache,igrid,ngrid,ind,iskip,nsplit,icell

  real(dp)::emag

  logical,dimension(1:nvector),save::ok,ok_leaf

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
 
  !------------------------------------------------------------
  ! Average down all MHD variables in split cells
  !------------------------------------------------------------
  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
 
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather split cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))>0
        end do
        
        ! Count split cells
        nsplit=0
        do i=1,ngrid
           if(ok(i))nsplit=nsplit+1
        end do
        
        ! Upload for selected cells
        if(nsplit>0)then
           icell=0
           do i=1,ngrid
              if(ok(i))then
                 icell=icell+1
                 ind_split(icell)=ind_cell(i)
              end if
           end do
           call upl(ind_split,nsplit)
        end if
        
     end do
     ! End loop over cells

  end do
  ! End loop over grids

  !-----------------------------------------------------------------
  ! Average down the magnetic field on each face of unsplit cells
  ! that is in contact with neighboring split cells
  !----------------------------------------------------------------
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left =nbor(ind_grid(i),2*idim-1)
           ind_right=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left )
           igridn(i,2*idim  )=son(ind_right)
        end do
     end do
 
    ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather unsplit cells
        do i=1,ngrid
           ok_leaf(i)=son(ind_cell(i))==0
        end do

        ! Loop over dimensions
        do idim=1,ndim

           ! Select unsplit cells with a refined left neighboring cell
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 ok(i)=ok_leaf(i).and.son(igridn(i,ig1)+ih1)>0
              else
                 ok(i)=.false.
              endif
           enddo

           ! Count selected cells
           nsplit=0
           do i=1,ngrid
              if(ok(i))nsplit=nsplit+1
           end do
        
           ! Upload for selected cells
           if(nsplit>0)then
              icell=0
              do i=1,ngrid
                 if(ok(i))then
                    icell=icell+1
                    ind_unsplit(icell)=ind_cell(i)
                    igrid_son  (icell)=son(igridn(i,ig1)+ih1)
                 end if
              end do
              if(interpol_var==1)then
                 ! Remove magnetic energy
                 do i=1,nsplit
                    emag=0.125d0*(uold(ind_unsplit(i),neul+idim)+ &
                         & uold(ind_unsplit(i),nvar+idim))**2
                    uold(ind_unsplit(i),neul)=uold(ind_unsplit(i),neul)-emag
                 end do
              endif
              call upl_left(ind_unsplit,igrid_son,idim,nsplit)
              if(interpol_var==1)then
                 ! Add magnetic energy
                 do i=1,nsplit
                    emag=0.125d0*(uold(ind_unsplit(i),neul+idim)+ &
                         & uold(ind_unsplit(i),nvar+idim))**2
                    uold(ind_unsplit(i),neul)=uold(ind_unsplit(i),neul)+emag
                 end do
              endif
           end if
        
           !  Select unsplit cells with a refined right neighboring cell
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 ok(i)=ok_leaf(i).and.son(igridn(i,ig2)+ih2)>0
              else
                 ok(i)=.false.
              endif
           enddo

           ! Count selected cells
           nsplit=0
           do i=1,ngrid
              if(ok(i))nsplit=nsplit+1
           end do
        
           ! Upload for selected cells
           if(nsplit>0)then
              icell=0
              do i=1,ngrid
                 if(ok(i))then
                    icell=icell+1
                    ind_unsplit(icell)=ind_cell(i)
                    igrid_son  (icell)=son(igridn(i,ig2)+ih2)
                 end if
              end do
              if(interpol_var==1)then
                 ! Remove magnetic energy
                 do i=1,nsplit
                    emag=0.125d0*(uold(ind_unsplit(i),neul+idim)+ &
                         & uold(ind_unsplit(i),nvar+idim))**2
                    uold(ind_unsplit(i),neul)=uold(ind_unsplit(i),neul)-emag
                 end do
              endif
              call upl_right(ind_unsplit,igrid_son,idim,nsplit)
              if(interpol_var==1)then
                 ! Add magnetic energy
                 do i=1,nsplit
                    emag=0.125d0*(uold(ind_unsplit(i),neul+idim)+ &
                         & uold(ind_unsplit(i),nvar+idim))**2
                    uold(ind_unsplit(i),neul)=uold(ind_unsplit(i),neul)+emag
                 end do
              endif
           end if
        
        end do
        ! End loop over dimensions

     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering upload_fine for level',i2)

end subroutine upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine upl(ind_cell,ncell)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ncell
  integer,dimension(1:nvector)::ind_cell
  !---------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the following variables:
  ! interpol_var=0: use rho, rho u and E
  ! interpol_tar=1: use rho, rho u and rho epsilon
  !---------------------------------------------------------------------
  integer ::ivar,i,idim,ind_son,iskip_son,ind,neul=5
  integer ,dimension(1:nvector),save::igrid_son,ind_cell_son
  real(dp),dimension(1:nvector),save::getx,ekin,emag
  integer,dimension(1:6,1:4)::hhh

  ! Get child oct index
  do i=1,ncell
     igrid_son(i)=son(ind_cell(i))
  end do

  !----------------------------------
  ! Loop over cell centered variables
  !----------------------------------
  do ivar=1,nvar
  if(ivar<=neul.or.ivar>neul+ndim)then

     ! Average conservative variable
     getx(1:ncell)=0.0d0
     do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        do i=1,ncell
           getx(i)=getx(i)+uold(ind_cell_son(i),ivar)
        end do
     end do
     
     ! Scatter result to cells
     do i=1,ncell
        uold(ind_cell(i),ivar)=getx(i)/dble(twotondim)
     end do

  end if
  end do
  ! End loop over cell centered variables

  ! Update cell centered magnetic field also in redundant array
#if NDIM==1
  do i=1,ncell
     uold(ind_cell(i),nvar+2)=uold(ind_cell(i),neul+2)
     uold(ind_cell(i),nvar+3)=uold(ind_cell(i),neul+3)
  end do
#endif
#if NDIM==2
  do i=1,ncell
     uold(ind_cell(i),nvar+3)=uold(ind_cell(i),neul+3)
  end do
#endif

  !----------------------------------
  ! Loop over face centered variables
  !----------------------------------
  hhh(1,1:4)=(/1,3,5,7/) 
  hhh(2,1:4)=(/2,4,6,8/) 
  hhh(3,1:4)=(/1,2,5,6/) 
  hhh(4,1:4)=(/3,4,7,8/) 
  hhh(5,1:4)=(/1,2,3,4/)
  hhh(6,1:4)=(/5,6,7,8/)

  ! Loop over dimensions
  do idim=1,ndim

     !----------------------
     ! Left B in parent cell
     !----------------------    
     getx(1:ncell)=0.0d0
     do ind=1,twotondim/2
        ind_son=hhh(2*idim-1,ind)
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        ! Update average
        do i=1,ncell
           getx(i)=getx(i)+uold(ind_cell_son(i),neul+idim)
        end do
     end do
     ! Scatter result to cells
     do i=1,ncell
        uold(ind_cell(i),neul+idim)=getx(i)/dble(twotondim/2)
     end do

     !-----------------------
     ! Right B in parent cell
     !-----------------------     
     getx(1:ncell)=0.0d0
     do ind=1,twotondim/2
        ind_son=hhh(2*idim,ind)
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        ! Update average
        do i=1,ncell
           getx(i)=getx(i)+uold(ind_cell_son(i),idim+nvar)
        end do
     end do
     ! Scatter result to cells
     do i=1,ncell
        uold(ind_cell(i),idim+nvar)=getx(i)/dble(twotondim/2)
     end do

  end do
  ! End loop over dimensions

  !------------------------
  ! Average internal energy
  !------------------------
  if(interpol_var==1)then

     getx(1:ncell)=0.0d0
     do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        ! Compute child kinetic energy
        ekin(1:ncell)=0.0d0
        do idim=1,3
           do i=1,ncell
              ekin(i)=ekin(i)+0.5d0*uold(ind_cell_son(i),1+idim)**2 &
                   &               /uold(ind_cell_son(i),1)
           end do
        end do
        ! Compute child magnetic energy
        emag(1:ncell)=0.0d0
        do idim=1,3
           do i=1,ncell
              emag(i)=emag(i)+0.125d0*(uold(ind_cell_son(i),neul+idim)+ &
                   &                   uold(ind_cell_son(i),nvar+idim))**2
           end do
        end do
        ! Update average
        do i=1,ncell
           getx(i)=getx(i)+uold(ind_cell_son(i),neul)-ekin(i)-emag(i)
        end do
     end do
        
     ! Compute new kinetic energy
     ekin(1:ncell)=0.0d0
     do idim=1,3
        do i=1,ncell
           ekin(i)=ekin(i)+0.5d0*uold(ind_cell(i),1+idim)**2 &
                &               /uold(ind_cell(i),1)
        end do
     end do
     ! Compute new magnetic energy
     emag(1:ncell)=0.0d0
     do idim=1,3
        do i=1,ncell
           emag(i)=emag(i)+0.125d0*(uold(ind_cell(i),neul+idim)+ &
                &                   uold(ind_cell(i),nvar+idim))**2
        end do
     end do
     
     ! Scatter result to cells
     do i=1,ncell
        uold(ind_cell(i),neul)=getx(i)/dble(twotondim)+ekin(i)+emag(i)
     end do
     
  endif

end subroutine upl
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine upl_left(ind_cell,igrid_son,idim,ncell)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ncell,idim
  integer,dimension(1:nvector)::ind_cell,igrid_son
  !---------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the magnetic field on cell faces
  !---------------------------------------------------------------------
  integer::i,ind_son,iskip_son,ind,neul=5
  integer ,dimension(1:nvector),save::ind_cell_son
  real(dp),dimension(1:nvector),save::getx
  integer,dimension(1:6,1:4)::hhh

  hhh(1,1:4)=(/1,3,5,7/)
  hhh(2,1:4)=(/2,4,6,8/)
  hhh(3,1:4)=(/1,2,5,6/)
  hhh(4,1:4)=(/3,4,7,8/)
  hhh(5,1:4)=(/1,2,3,4/)
  hhh(6,1:4)=(/5,6,7,8/)

  !--------------------------------------------------
  ! Left B in parent cell is computed as the average
  ! over right B in all left children cells
  !--------------------------------------------------
  getx(1:ncell)=0.0d0
  do ind=1,twotondim/2
     ind_son=hhh(2*idim,ind)
     iskip_son=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncell
        ind_cell_son(i)=iskip_son+igrid_son(i)
     end do
     ! Update average
     do i=1,ncell
        getx(i)=getx(i)+uold(ind_cell_son(i),idim+nvar)
     end do
  end do
  ! Scatter result to cells
  do i=1,ncell
     uold(ind_cell(i),idim+neul)=getx(i)/dble(twotondim/2)
  end do
  
end subroutine upl_left
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine upl_right(ind_cell,igrid_son,idim,ncell)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ncell,idim
  integer,dimension(1:nvector)::ind_cell,igrid_son
  !---------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the magnetic field on cell faces
  !---------------------------------------------------------------------
  integer::i,ind_son,iskip_son,ind,neul=5
  integer ,dimension(1:nvector),save::ind_cell_son
  real(dp),dimension(1:nvector),save::getx
  integer,dimension(1:6,1:4)::hhh

  hhh(1,1:4)=(/1,3,5,7/) 
  hhh(2,1:4)=(/2,4,6,8/) 
  hhh(3,1:4)=(/1,2,5,6/) 
  hhh(4,1:4)=(/3,4,7,8/) 
  hhh(5,1:4)=(/1,2,3,4/)
  hhh(6,1:4)=(/5,6,7,8/)

  !--------------------------------------------------
  ! Right B in parent cell is computed as the average 
  ! over left B in all right children cells
  !--------------------------------------------------
  getx(1:ncell)=0.0d0
  do ind=1,twotondim/2
     ind_son=hhh(2*idim-1,ind)
     iskip_son=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncell
        ind_cell_son(i)=iskip_son+igrid_son(i)
     end do
     ! Update average
     do i=1,ncell
        getx(i)=getx(i)+uold(ind_cell_son(i),idim+neul)
     end do
  end do
  ! Scatter result to cells
  do i=1,ncell
     uold(ind_cell(i),idim+nvar)=getx(i)/dble(twotondim/2)
  end do

end subroutine upl_right
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_hydro(u1,ind1,u2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3)::u1
  integer ,dimension(1:nvector,0:twondim)           ::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  ! interpol_var=0: rho, rho u and E
  ! interpol_var=1: rho, rho u and rho epsilon
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation without limiters
  !----------------------------------------------------------
  integer::i,j,ivar,idim,ind,ix,iy,iz,neul=5

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w
  real(dp),dimension(1:nvector),save::ekin,emag
  real(dp),dimension(1:nvector,0:twondim  ,1:6),save::B1
  real(dp),dimension(1:nvector,1:twotondim,1:6),save::B2

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! If necessary, convert father total energy into internal energy
  if(interpol_var==1)then
     do j=0,twondim
        ekin(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u1(i,j,idim+1)**2/u1(i,j,1)
           end do
        end do
        emag(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              emag(i)=emag(i)+0.125d0*(u1(i,j,idim+neul)+u1(i,j,idim+nvar))**2
           end do
        end do
        do i=1,nn
           u1(i,j,neul)=u1(i,j,neul)-ekin(i)-emag(i)
        end do
     end do
  end if


  !------------------------------------------------
  ! Loop over cell-centered interpolation variables
  !------------------------------------------------
  do ivar=1,nvar
  if(ivar<=neul.or.ivar>neul+ndim)then

     ! Load father variable
     do j=0,twondim
        do i=1,nn 
           a(i,j)=u1(i,j,ivar)
        end do
     end do

     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w,nn)
     if(interpol_type==2)call compute_limiter_central(a,w,nn)
     if(interpol_type==3)call compute_central(a,w,nn)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(1:nn,ind,ivar)=a(1:nn,0)
        do idim=1,ndim
           do i=1,nn
              u2(i,ind,ivar)=u2(i,ind,ivar)+w(i,idim)*xc(ind,idim)
           end do
        end do
     end do

  end if
  end do
  ! End loop over cell-centered variables

  ! Update cell centered magnetic field also in redundant array
#if NDIM<2
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,2+nvar)=u2(i,ind,2+neul)
     end do
  end do
#endif
#if NDIM<3
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,3+nvar)=u2(i,ind,3+neul)
     end do
  end do
#endif

  !------------------------------------------------
  ! Loop over face-centered interpolation variables
  !------------------------------------------------
  do j=0,twondim
     do i=1,nn
        B1(i,j,1)=u1(i,j,neul+1)
        B1(i,j,2)=u1(i,j,neul+2)
        B1(i,j,3)=u1(i,j,neul+3)
        B1(i,j,4)=u1(i,j,nvar+1)
        B1(i,j,5)=u1(i,j,nvar+2)
        B1(i,j,6)=u1(i,j,nvar+3)
     end do
  end do
  call interpol_mag(B1,ind1,B2,nn)
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,neul+1)=B2(i,ind,1)
        u2(i,ind,nvar+1)=B2(i,ind,4)
#if NDIM>1        
        u2(i,ind,neul+2)=B2(i,ind,2)
        u2(i,ind,nvar+2)=B2(i,ind,5)
#endif
#if NDIM>2
        u2(i,ind,neul+3)=B2(i,ind,3)
        u2(i,ind,nvar+3)=B2(i,ind,6)
#endif
     end do
  end do

  ! If necessary, convert children internal energy into total energy
  if(interpol_var==1)then
     do ind=1,twotondim
        ekin(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u2(i,ind,idim+1)**2/u2(i,ind,1)
           end do
        end do
        emag(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              emag(i)=emag(i)+0.125d0*(u2(i,ind,idim+neul)+u2(i,ind,idim+nvar))**2
           end do
        end do
        do i=1,nn
           u2(i,ind,neul)=u2(i,ind,neul)+ekin(i)+emag(i)
        end do
     end do
  end if

end subroutine interpol_hydro
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_minmod(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------
  ! MinMod slope
  !---------------
  integer::i,idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     do i=1,nn
        diff_left=0.5*(a(i,2*idim)-a(i,0))
        diff_right=0.5*(a(i,0)-a(i,2*idim-1))
        if(diff_left*diff_right<=0.0)then
           minmod=0.0
        else
           minmod=MIN(ABS(diff_left),ABS(diff_right)) &
                &   *diff_left/ABS(diff_left)
        end if
        w(i,idim)=minmod
     end do
  end do

end subroutine compute_limiter_minmod
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------
  ! MinMod slope
  !---------------
  integer::i,idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.5*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

end subroutine compute_central
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------------------
  ! Monotonized Central slope
  !---------------------------
  integer::i,j,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::xxc
  real(dp),dimension(1:nvector,1:twotondim),save::ac
  real(dp),dimension(1:nvector),save::corner,kernel,diff_corner,diff_kernel
  real(dp),dimension(1:nvector),save::max_limiter,min_limiter,limiter

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Second order central slope
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

  ! Compute corner interpolated values
  do ind=1,twotondim
     do i=1,nn
        ac(i,ind)=a(i,0)
     end do
  end do
  do idim=1,ndim
     do ind=1,twotondim
        xxc = xc(ind,idim)
        do i=1,nn
           corner(i)=ac(i,ind)+2.D0*w(i,idim)*xxc
        end do
        do i=1,nn
           ac(i,ind)=corner(i)
        end do
     end do
  end do

  ! Compute max of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MAX(corner(i),ac(i,j))
     end do
  end do

  ! Compute max of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MAX(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  max_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        max_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute min of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MIN(corner(i),ac(i,j))
     end do
  end do

  ! Compute min of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MIN(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  min_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        min_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute limiter
  do i=1,nn
     limiter(i)=MIN(min_limiter(i),max_limiter(i))
  end do

  ! Correct gradient with limiter
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=w(i,idim)*limiter(i)
     end do
  end do

end subroutine compute_limiter_central
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_mag(B1,ind1,B2,nn)
  use amr_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:6)::B1
  integer ,dimension(1:nvector,0:twondim)      ::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:6)::B2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are Bx, By and Bz.
  ! Divergence free is garanteed.
  ! The scheme is the one invented by Toth and Balsara.
  ! interpol_type=0: straight injection
  ! interpol_type=1: linear interpolation with MinMod slope
  ! interpol_type=2: linear interpolation with Monotonized Central slope
  ! interpol_type=3: linear interpolation without limiters
  !----------------------------------------------------------
  integer::i,j,k,ind,l,idim,imax,jmax,kmax
  real(dp),dimension(1:nvector,-1:1,0:1,0:1),save::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1),save::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1),save::w

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Compute interpolated fine B over coarse side faces
  call interpol_faces(B1,u,v,w,nn)
  
  ! Get fine B from refined faces, if any
  call copy_from_refined_faces(B1,ind1,u,v,w,nn)
 
  ! Compute interpolated fine B inside coarse cell.
  call cmp_central_faces(u,v,w,nn)

  ! Scatter results
  do i=0,imax
  do j=0,jmax
  do k=0,kmax
     ind=1+i+2*j+4*k
     do l=1,nn
        B2(l,ind,1)=u(l,i-1,j,k)
        B2(l,ind,2)=v(l,i,j-1,k)
        B2(l,ind,3)=w(l,i,j,k-1)
        B2(l,ind,4)=u(l,i,j,k)
        B2(l,ind,5)=v(l,i,j,k)
        B2(l,ind,6)=w(l,i,j,k)
     end do
  end do
  end do
  end do

end subroutine interpol_mag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_faces(b1,u,v,w,nn)
  use amr_commons
  use hydro_commons, ONLY: interpol_type
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim,1:6)::b1
  real(dp),dimension(1:nvector,-1:1,0:1,0:1)::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1)::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1)::w

  ! TVD interpolation from coarse faces
  integer::i,j,k,l,imax,jmax,kmax
  real(dp),dimension(1:nvector,0:4),save::b
  real(dp),dimension(1:nvector,1:2),save::s

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Left face along direction x (interpolate Bx)
  do l=1,nn
     b(l,0)=b1(l,0,1)
  end do
#if NDIM>1     
  do l=1,nn
     b(l,1)=b1(l,3,1)
     b(l,2)=b1(l,4,1)
  end do
#endif
#if NDIM>2
  do l=1,nn
     b(l,3)=b1(l,5,1)
     b(l,4)=b1(l,6,1)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,nn)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,nn)
#endif
  do j=0,jmax
  do k=0,kmax
     do l=1,nn
        u(l,-1,j,k)=b(l,0)+0.5*s(l,1)*(dble(j)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

  ! Right face along direction x (interpolate Bx)
  do l=1,nn
     b(l,0)=b1(l,0,4)
  end do
#if NDIM>1     
  do l=1,nn
     b(l,1)=b1(l,3,4)
     b(l,2)=b1(l,4,4)
  end do
#endif
#if NDIM>2
  do l=1,nn
     b(l,3)=b1(l,5,4)
     b(l,4)=b1(l,6,4)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,nn)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,nn)
#endif
  do j=0,jmax
  do k=0,kmax
     do l=1,nn
        u(l,+1,j,k)=b(l,0)+0.5*s(l,1)*(dble(j)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

#if NDIM>1
  ! Left face along direction y (interpolate By)
  do l=1,nn
     b(l,0)=b1(l,0,2)
  end do
  do l=1,nn
     b(l,1)=b1(l,1,2)
     b(l,2)=b1(l,2,2)
  end do
#if NDIM>2     
  do l=1,nn
     b(l,3)=b1(l,5,2)
     b(l,4)=b1(l,6,2)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,nn)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,nn)
#endif
  do i=0,imax
  do k=0,kmax
     do l=1,nn
        v(l,i,-1,k)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

  ! Right face along direction y (interpolate By)
  do l=1,nn
     b(l,0)=b1(l,0,5)
  end do
  do l=1,nn
     b(l,1)=b1(l,1,5)
     b(l,2)=b1(l,2,5)
  end do
#if NDIM>2     
  do l=1,nn
     b(l,3)=b1(l,5,5)
     b(l,4)=b1(l,6,5)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,nn)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,nn)
#endif
  do i=0,imax
  do k=0,kmax
     do l=1,nn
        v(l,i,+1,k)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do
#endif

#if NDIM>2
  ! Left face along direction z (interpolate Bz)
  do l=1,nn
     b(l,0)=b1(l,0,3)
     b(l,1)=b1(l,1,3)
     b(l,2)=b1(l,2,3)
     b(l,3)=b1(l,3,3)
     b(l,4)=b1(l,4,3)
  end do

  s(1:nn,1:2)=0.0
  if(interpol_type>0)call compute_2d_tvd(b,s,nn)
  do i=0,1
     do j=0,1
        do l=1,nn
           w(l,i,j,-1)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(j)-0.5)
        end do
     end do
  end do

  ! Right face along direction z (interpolate Bz)
  do l=1,nn
     b(l,0)=b1(l,0,6)
     b(l,1)=b1(l,1,6)
     b(l,2)=b1(l,2,6)
     b(l,3)=b1(l,3,6)
     b(l,4)=b1(l,4,6)
  end do

  s(1:nn,1:2)=0.0
  if(interpol_type>0)call compute_2d_tvd(b,s,nn)
  do i=0,1
     do j=0,1
        do l=1,nn
           w(l,i,j,+1)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(j)-0.5)
        end do
     end do
  end do
#endif

end subroutine interpol_faces
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine copy_from_refined_faces(b1,ind1,u,v,w,nn)
  use amr_commons
  use hydro_commons, ONLY: nvar,uold
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim,1:6)::b1
  integer ,dimension(1:nvector,0:twondim)::ind1
  real(dp),dimension(1:nvector,-1:1,0:1,0:1)::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1)::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1)::w

  ! TVD interpolation from coarse faces
  integer::i,j,k,l,ind,iskip,imax,jmax,kmax,neul=5

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Left face along direction x (interpolate Bx)
  do j=0,jmax
  do k=0,kmax
     ind=1+1+j*2+k*4
     iskip=ncoarse+(ind-1)*ngridmax
     do l=1,nn
        if(ind1(l,1)>0)then
           u(l,-1,j,k)=uold(iskip+ind1(l,1),nvar+1)
        end if
     end do
  end do
  end do

  ! Right face along direction x (interpolate Bx)
  do j=0,jmax
  do k=0,kmax
     ind=1+0+j*2+k*4
     iskip=ncoarse+(ind-1)*ngridmax
     do l=1,nn
        if(ind1(l,2)>0)then
           u(l,+1,j,k)=uold(iskip+ind1(l,2),neul+1)
        end if
     end do
  end do
  end do

#if NDIM>1
  ! Left face along direction y (interpolate By)
  do i=0,imax
  do k=0,kmax
     ind=1+i+1*2+k*4
     iskip=ncoarse+(ind-1)*ngridmax
     do l=1,nn
        if(ind1(l,3)>0)then
           v(l,i,-1,k)=uold(iskip+ind1(l,3),nvar+2)
        end if
     end do
  end do
  end do

  ! Right face along direction y (interpolate By)
  do i=0,imax
  do k=0,kmax
     ind=1+i+0*2+k*4
     iskip=ncoarse+(ind-1)*ngridmax
     do l=1,nn
        if(ind1(l,4)>0)then
           v(l,i,+1,k)=uold(iskip+ind1(l,4),neul+2)
        end if
     end do
  end do
  end do
#endif

#if NDIM>2
  ! Left face along direction z (interpolate Bz)
  do i=0,imax
  do j=0,kmax
     ind=1+i+j*2+1*4
     iskip=ncoarse+(ind-1)*ngridmax
     do l=1,nn
        if(ind1(l,5)>0)then
           w(l,i,j,-1)=uold(iskip+ind1(l,5),nvar+3)
        end if
     end do
  end do
  end do

  ! Right face along direction z (interpolate Bz)
  do i=0,imax
  do j=0,kmax
     ind=1+i+j*2+0*4
     iskip=ncoarse+(ind-1)*ngridmax
     do l=1,nn
        if(ind1(l,6)>0)then
           w(l,i,j,+1)=uold(iskip+ind1(l,6),neul+3)
        end if
     end do
  end do
  end do
#endif

end subroutine copy_from_refined_faces
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_central_faces(u,v,w,nn)
  use amr_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,-1:1,0:1,0:1)::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1)::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1)::w

  integer::i,j,k,l,ii,jj,kk,imax,jmax,kmax
  real(dp),dimension(1:nvector),save::UXX,VYY,WZZ,UXYZ,VXYZ,WXYZ

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  do l = 1,nn
     UXX (l)=0.0_dp
     VYY (l)=0.0_dp
     WZZ (l)=0.0_dp
     UXYZ(l)=0.0_dp
     VXYZ(l)=0.0_dp
     WXYZ(l)=0.0_dp
  end do

#if NDIM==2
  do i=0,imax
  do j=0,jmax
  do k=0,kmax
     ii=2*i-1
     jj=2*j-1
     do l = 1,nn
        UXX (l)=UXX (l)+(ii*jj*v(l,i,jj,k))*0.25
        VYY (l)=VYY (l)+(ii*jj*u(l,ii,j,k))*0.25
     enddo
  enddo
  enddo
  enddo
#endif

#if NDIM==3
  do i=0,imax
  do j=0,jmax
  do k=0,kmax
     ii=2*i-1
     jj=2*j-1
     kk=2*k-1
     do l = 1,nn
        UXX (l)=UXX (l)+(ii*jj*v(l,i,jj,k)+ii*kk*w(l,i,j,kk))*0.125
        VYY (l)=VYY (l)+(jj*kk*w(l,i,j,kk)+ii*jj*u(l,ii,j,k))*0.125
        WZZ (l)=WZZ (l)+(ii*kk*u(l,ii,j,k)+jj*kk*v(l,i,jj,k))*0.125
        UXYZ(l)=UXYZ(l)+(ii*jj*kk*u(l,ii,j,k))*0.125
        VXYZ(l)=VXYZ(l)+(ii*jj*kk*v(l,i,jj,k))*0.125
        WXYZ(l)=WXYZ(l)+(ii*jj*kk*w(l,i,j,kk))*0.125
     enddo
  enddo
  enddo
  enddo
#endif

#if NDIM==1
  ! Bx on central faces
  do j=0,jmax
  do k=0,kmax
     do l = 1,nn
        u(l,0,j,k)=0.5*(u(l,-1,j,k)+u(l,+1,j,k))
     enddo
  enddo
  enddo
#endif
#if NDIM==2
  ! Bx on central faces
  do j=0,jmax
  do k=0,kmax
     do l = 1,nn
        u(l,0,j,k)=0.5*(u(l,-1,j,k)+u(l,+1,j,k)) + UXX(l)
     enddo
  enddo
  enddo
  do i=0,imax
  do k=0,kmax
     do l = 1,nn
        v(l,i,0,k)=0.5*(v(l,i,-1,k)+v(l,i,+1,k)) + VYY(l)
     enddo
  enddo
  enddo
#endif

#if NDIM==3
  ! Bx on central faces
  do j=0,jmax
  do k=0,kmax
     do l = 1,nn
        u(l,0,j,k)=0.5*(u(l,-1,j,k)+u(l,+1,j,k)) + UXX(l)     &
             &      + (dble(k)-0.5)*VXYZ(l) + (dble(j)-0.5)*WXYZ(l)
     enddo
  enddo
  enddo
  do i=0,imax
  do k=0,kmax
     do l = 1,nn
        v(l,i,0,k)=0.5*(v(l,i,-1,k)+v(l,i,+1,k)) + VYY(l)     &
             &      + (dble(i)-0.5)*WXYZ(l) + (dble(k)-0.5)*UXYZ(l)
     enddo
  enddo
  enddo
  do i=0,imax
  do j=0,jmax
     do l = 1,nn
        w(l,i,j,0)=0.5*(w(l,i,j,-1)+w(l,i,j,+1)) + WZZ(l)     &
             &      + (dble(j)-0.5)*UXYZ(l) + (dble(i)-0.5)*VXYZ(l)
     enddo
  enddo
  enddo
#endif

end subroutine cmp_central_faces
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_2d_tvd(b,s,nn)
  use amr_commons, ONLY: nvector
  use hydro_commons, ONLY: interpol_type
  use const
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:4)::b
  real(dp),dimension(1:nvector,1:2)::s
  
  integer::i
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop

  if(interpol_type==3)then
     do i=1,nn
        dlft = half*(b(i,0) - b(i,1))
        drgt = half*(b(i,2) - b(i,0))
        s(i,1) = dlft+drgt
     end do
     do i=1,nn
        dlft = half*(b(i,0) - b(i,3))
        drgt = half*(b(i,4) - b(i,0))
        s(i,2) = dlft+drgt
     end do
     return
  endif

  do i=1,nn
     dlft = interpol_type*(b(i,0) - b(i,1))
     drgt = interpol_type*(b(i,2) - b(i,0))
     dcen = half*(dlft+drgt)/interpol_type
     dsgn = sign(one, dcen)
     slop = min(abs(dlft),abs(drgt))
     dlim = slop
     if((dlft*drgt)<=zero)dlim=zero
     s(i,1) = dsgn*min(dlim,abs(dcen))
  end do

  do i=1,nn
     dlft = interpol_type*(b(i,0) - b(i,3))
     drgt = interpol_type*(b(i,4) - b(i,0))
     dcen = half*(dlft+drgt)/interpol_type
     dsgn = sign(one, dcen)
     slop = min(abs(dlft),abs(drgt))
     dlim = slop
     if((dlft*drgt)<=zero)dlim=zero
     s(i,2) = dsgn*min(dlim,abs(dcen))
  end do


end subroutine compute_2d_tvd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_1d_tvd(b,s,nn)
  use amr_commons, ONLY: nvector
  use hydro_commons, ONLY: interpol_type
  use const
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:4)::b
  real(dp),dimension(1:nvector,1:2)::s
  
  integer::i
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop

  if(interpol_type==3)then
     do i=1,nn
        dlft = half*(b(i,0) - b(i,1))
        drgt = half*(b(i,2) - b(i,0))
        s(i,1) = dlft+drgt
     end do
     return
  endif
  do i=1,nn
     dlft = interpol_type*(b(i,0) - b(i,1))
     drgt = interpol_type*(b(i,2) - b(i,0))
     dcen = half*(dlft+drgt)/interpol_type
     dsgn = sign(one, dcen)
     slop = min(abs(dlft),abs(drgt))
     dlim = slop
     if((dlft*drgt)<=zero)dlim=zero
     s(i,1) = dsgn*min(dlim,abs(dcen))
  end do

end subroutine compute_1d_tvd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
