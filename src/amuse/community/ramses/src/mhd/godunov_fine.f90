!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call godfine1(ind_grid,ngrid,ilevel)
  end do

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip
  real(dp)::d,u,v,w,e,A,B,C

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar+3
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(pressure_fix)then
        do i=1,active(ilevel)%ngrid
           divu(active(ilevel)%igrid(i)+iskip) = 0.0
        end do
        do i=1,active(ilevel)%ngrid
           d=uold(active(ilevel)%igrid(i)+iskip,1)
           u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           A=0.5*(uold(active(ilevel)%igrid(i)+iskip,6)+uold(active(ilevel)%igrid(i)+iskip,nvar+1))
           B=0.5*(uold(active(ilevel)%igrid(i)+iskip,7)+uold(active(ilevel)%igrid(i)+iskip,nvar+2))
           C=0.5*(uold(active(ilevel)%igrid(i)+iskip,8)+uold(active(ilevel)%igrid(i)+iskip,nvar+3))
           e=uold(active(ilevel)%igrid(i)+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
           enew(active(ilevel)%igrid(i)+iskip)=e
        end do
     end if
  end do

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar+3
        do i=1,reception(icpu,ilevel)%ngrid
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
     end do
     if(pressure_fix)then
        do i=1,reception(icpu,ilevel)%ngrid
           divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
        end do
     end if
  end do
  end do

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold to its new value unew after the
  ! hydro step.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc
  real(dp)::scale,d,u,v,w,A,B,C
  real(dp)::e_mag,e_kin,e_cons,e_prim,e_trunc,div,dx,fact,d_old

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

  ! Add gravity source term at time t with half time step                                                                                    
  if(poisson)then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           d=unew(active(ilevel)%igrid(i)+iskip,1)
           u=unew(active(ilevel)%igrid(i)+iskip,2)/d
           v=unew(active(ilevel)%igrid(i)+iskip,3)/d
           w=unew(active(ilevel)%igrid(i)+iskip,4)/d
           e_kin=0.5*d*(u**2+v**2+w**2)
           e_prim=unew(active(ilevel)%igrid(i)+iskip,5)-e_kin
           d_old=uold(active(ilevel)%igrid(i)+iskip,1)
           fact=d_old/d*0.5*dtnew(ilevel)
           if(ndim>0)then
              u=u+f(active(ilevel)%igrid(i)+iskip,1)*fact
              unew(active(ilevel)%igrid(i)+iskip,2)=d*u
           endif
           if(ndim>1)then
              v=v+f(active(ilevel)%igrid(i)+iskip,2)*fact
              unew(active(ilevel)%igrid(i)+iskip,3)=d*v
           end if
           if(ndim>2)then
              w=w+f(active(ilevel)%igrid(i)+iskip,3)*fact
              unew(active(ilevel)%igrid(i)+iskip,4)=d*w
           endif
           e_kin=0.5*d*(u**2+v**2+w**2)
           unew(active(ilevel)%igrid(i)+iskip,5)=e_prim+e_kin
        end do
     end do
  end if

  ! Set uold to unew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar+3
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(pressure_fix)then
        fact=(gamma-1.0d0)
        do i=1,active(ilevel)%ngrid
           d=uold(active(ilevel)%igrid(i)+iskip,1)
           u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           A=0.5*(uold(active(ilevel)%igrid(i)+iskip,6)+uold(active(ilevel)%igrid(i)+iskip,nvar+1))
           B=0.5*(uold(active(ilevel)%igrid(i)+iskip,7)+uold(active(ilevel)%igrid(i)+iskip,nvar+2))
           C=0.5*(uold(active(ilevel)%igrid(i)+iskip,8)+uold(active(ilevel)%igrid(i)+iskip,nvar+3))
           e_kin=0.5*d*(u**2+v**2+w**2)
           e_mag=0.5*(A**2+B**2+C**2)
           e_cons=uold(active(ilevel)%igrid(i)+iskip,5)-e_kin-e_mag
           e_prim=enew(active(ilevel)%igrid(i)+iskip)*(1.0d0+fact* &
                & divu(active(ilevel)%igrid(i)+iskip))! Note: here divu=-div.u*dt
           div=abs(divu(active(ilevel)%igrid(i)+iskip))*dx/dtnew(ilevel)
           ! Estimate of the local truncation errors
           e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
           if(e_cons<e_trunc)then
              uold(active(ilevel)%igrid(i)+iskip,5)=e_prim+e_kin+e_mag
           end if
        end do
     end if
  end do

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1=0.0d0
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2=0.0d0

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,1:3,1:3,1:3),save::emfx,emfy,emfz
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::neul=5
  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dflux_x,dflux_y,dflux_z
  real(dp)::dx,scale,oneontwotondim
  real(dp)::dflux,weight

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     
     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do
        
        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if
        
        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call mag_unsplit(uloc,gloc,flux,emfx,emfy,emfz,tmp,dx,dx,dx,dtnew(ilevel),ncache)

  if(scheme.eq.'induction')then
  !---------------------------------
  ! Reset all Euler variables fluxes
  !---------------------------------
  do idim=1,ndim
  do ivar=1,neul
     do k3=k3min,k3max+1
     do j3=j3min,j3max+1
     do i3=i3min,i3max+1
        do i=1,ncache
           flux(i,i3,j3,k3,ivar,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
  end do

  endif

  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do
  !---------------------------------------------------------------
  ! Reset Euler fluxes for Bx
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,6,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#if NDIM>1
  !---------------------------------------------------------------
  ! Reset electromotive force along direction z at refined edges
  !---------------------------------------------------------------
  do k3=k3min,k3max
  do j3=1,3
  do i3=1,3
     do i=1,ncache
        if(ok(i,i3  ,j3  ,k3) .or. ok(i,i3  ,j3-1,k3) .or.  &
         & ok(i,i3-1,j3  ,k3) .or. ok(i,i3-1,j3-1,k3))then
           emfz(i,i3,j3,k3)=0.0d0
#if NDIM==2
           emfz(i,i3,j3,k3+1)=0.0d0
#endif
        end if
     end do
  end do
  end do
  end do
  !---------------------------------------------------------------
  ! Reset Euler fluxes for By
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,7,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#endif
#if NDIM>2
  !---------------------------------------------------------------
  ! Reset electromotive force along direction y at refined edges
  !---------------------------------------------------------------
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
  !---------------------------------------------------------------
  ! Reset electromotive force along direction x at refined edges
  !---------------------------------------------------------------
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
  !---------------------------------------------------------------
  ! Reset Euler fluxes for Bz
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,8,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#endif

  !-----------------------------------------------------
  ! Conservative update at level ilevel for Euler system
  !-----------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do ivar=1,nvar
           do i=1,ncache
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        do ivar=1,3
           do i=1,ncache
              unew(ind_cell(i),nvar+ivar)=unew(ind_cell(i),nvar+ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,neul+ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,neul+ivar,idim))
           end do
        end do

        if(pressure_fix)then
        ! Update velocity divergence
        do i=1,ncache
           divu(ind_cell(i))=divu(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
        end do
        ! Update internal energy
        do i=1,ncache
           enew(ind_cell(i))=enew(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
        end do
        end if
     end do
     end do
     end do
  end do

  !---------------------------------------------------------
  ! Conservative update at level ilevel for induction system
  !---------------------------------------------------------
  do k3=k3min,k3max
  do j3=j3min,j3max
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
        unew(ind_cell(i),neul+1)=unew(ind_cell(i),neul+1)+dflux_x
        dflux_x=( emfy(i,i3+1,j3,k3)-emfy(i,i3+1,j3,k3+1) ) &
          &    -( emfz(i,i3+1,j3,k3)-emfz(i,i3+1,j3+1,k3) )  
        unew(ind_cell(i),nvar+1)=unew(ind_cell(i),nvar+1)+dflux_x
     end do
  end do
  end do
  end do
#if NDIM>1  
  do k3=k3min,k3max
  do j3=1,2
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update By using constraint transport
     do i=1,ncache
        dflux_y=( emfz(i,i3,j3,k3)-emfz(i,i3+1,j3,k3) ) &
          &    -( emfx(i,i3,j3,k3)-emfx(i,i3,j3,k3+1) )
        unew(ind_cell(i),neul+2)=unew(ind_cell(i),neul+2)+dflux_y
        dflux_y=( emfz(i,i3,j3+1,k3)-emfz(i,i3+1,j3+1,k3) ) &
          &    -( emfx(i,i3,j3+1,k3)-emfx(i,i3,j3+1,k3+1) )
        unew(ind_cell(i),nvar+2)=unew(ind_cell(i),nvar+2)+dflux_y
     end do
  end do
  end do
  end do
#endif
#if NDIM>2
  do k3=1,2
  do j3=1,2
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update Bz using constraint transport
     do i=1,ncache
        dflux_z=( emfx(i,i3,j3,k3)-emfx(i,i3,j3+1,k3) ) &
          &    -( emfy(i,i3,j3,k3)-emfy(i,i3+1,j3,k3) )
        unew(ind_cell(i),neul+3)=unew(ind_cell(i),neul+3)+dflux_z
        dflux_z=( emfx(i,i3,j3,k3+1)-emfx(i,i3,j3+1,k3+1) ) &
          &    -( emfy(i,i3,j3,k3+1)-emfy(i,i3+1,j3,k3+1) )
        unew(ind_cell(i),nvar+3)=unew(ind_cell(i),nvar+3)+dflux_z
     end do
  end do
  end do
  end do
#endif

  if(ilevel>levelmin)then

  !-----------------------------------------------------------
  ! Conservative update at level ilevel-1 for the Euler system
  !-----------------------------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     
     !----------------------
     ! Left flux at boundary
     !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     do ivar=1,3
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),nvar+ivar)=unew(ind_buffer(i),nvar+ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,neul+ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if
     
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     do ivar=1,3
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              unew(ind_buffer(i),nvar+ivar)=unew(ind_buffer(i),nvar+ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,neul+ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if

  end do
  ! End loop over dimensions

#if NDIM>1  
  !---------------------------------------------------------------
  ! Conservative update at level ilevel-1 for the induction system
  !---------------------------------------------------------------
  i1=1; j1=0; k1=0
  if(ndim>1)j1=1
  if(ndim>2)k1=1

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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,1,1)+emfz(i,1,1,2))*0.25*weight
     unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)+dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)+dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)-dflux
     unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux*0.5
        unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,3,1)+emfz(i,1,3,2))*0.25*weight
     unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)-dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)-dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)-dflux
     unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
        unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,3,1)+emfz(i,3,3,2))*0.25*weight
     unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)-dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)+dflux
     unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)+dflux*0.5
        unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,1,1)+emfz(i,3,1,2))*0.25*weight
     unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)+dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)+dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)+dflux
     unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
        unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)-dflux*0.5
     endif
  end do
#endif
  !--------------------------------------
  ! Deal with 4 EMFx edges
  !--------------------------------------
#if NDIM>2
  ! Update coarse By and Bz using fine EMFx on Y=0 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
  ind_father2=1+(i1  )+3*(j1-1)+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,1)+emfx(i,2,1,1))*0.25*weight
     unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)+dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)+dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)-dflux
     unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux*0.5
        unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,3)+emfx(i,2,1,3))*0.25*weight
     unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)-dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)-dflux
     unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
        unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)+dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,3)+emfx(i,2,3,3))*0.25*weight
     unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)-dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)-dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)+dflux
     unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
        unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)-dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,1)+emfx(i,2,3,1))*0.25*weight
     unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)+dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)+dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)+dflux
     unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux*0.5
        unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,1)+emfy(i,1,2,1))*0.25*weight
     unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)-dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)-dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)+dflux
     unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux*0.5
        unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,3)+emfy(i,1,2,3))*0.25*weight
     unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)+dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)+dflux
     unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)-dflux*0.5
        unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux*0.5
     endif
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
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,3)+emfy(i,3,2,3))*0.25*weight
     unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)+dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)-dflux
     unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)-dflux*0.5
        unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)+dflux*0.5
     endif
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=0 grid edge
  ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1  )+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,1)+emfy(i,3,2,1))*0.25*weight
     unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)-dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)-dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)-dflux
     unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux*0.5
        unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)+dflux*0.5
     endif
  end do
#endif

  end if

end subroutine godfine1
