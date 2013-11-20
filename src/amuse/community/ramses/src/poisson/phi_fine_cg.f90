!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine phi_fine_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !=========================================================
  ! Iterative Poisson solver with Conjugate Gradient method 
  ! to solve A x = b
  ! r  : stored in f(i,1)
  ! p  : stored in f(i,2)
  ! A p: stored in f(i,3)
  ! x  : stored in phi(i)
  ! b  : stored in rho(i)
  !=========================================================
  integer::i,idim,info,ind,iter,iskip,itermax,nx_loc
  integer::idx
  real(dp)::error,error_ini
  real(dp)::dx2,fourpi,scale,oneoversix,fact,fact2
  real(dp)::r2_old,alpha_cg,beta_cg
  real(kind=8)::r2,pAp,rhs_norm,r2_all,pAp_all,rhs_norm_all

  if(gravity_type>0)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel


  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2
  fact2 = fact*fact

  !===============================
  ! Compute initial phi
  !===============================
   if(ilevel>levelmin)then
      call make_initial_phi(ilevel,icount)              ! Interpolate phi down
   else
      call make_multipole_phi(ilevel)            ! Fill up with simple initial guess
   endif
   call make_virtual_fine_dp(phi(1),ilevel)      ! Update boundaries
   call make_boundary_phi(ilevel)                ! Update physical boundaries

  !===============================
  ! Compute right-hand side norm
  !===============================
  rhs_norm=0.d0
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        idx=active(ilevel)%igrid(i)+iskip
        rhs_norm=rhs_norm+fact2*(rho(idx)-rho_tot)*(rho(idx)-rho_tot)
     end do
  end do
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       & MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
#endif
  rhs_norm=DSQRT(rhs_norm/dble(twotondim*numbtot(1,ilevel)))

  !==============================================
  ! Compute r = b - Ax and store it into f(i,1)
  ! Also set p = r and store it into f(i,2)
  !==============================================
  call cmp_residual_cg(ilevel,icount)

  !====================================
  ! Main iteration loop
  !====================================
  iter=0; itermax=10000
  error=1.0D0; error_ini=1.0D0
  do while(error>epsilon*error_ini.and.iter<itermax)

     iter=iter+1

     !====================================
     ! Compute residual norm
     !====================================
     r2=0.0d0
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           r2=r2+f(idx,1)*f(idx,1)
        end do
     end do
     ! Compute global norm
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(r2,r2_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
     r2=r2_all
#endif
     !====================================
     ! Compute beta factor
     !====================================
     if(iter==1)then
        beta_cg=0.
     else
        beta_cg=r2/r2_old
     end if
     r2_old=r2

     !====================================
     ! Recurrence on p
     !====================================
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           f(idx,2)=f(idx,1)+beta_cg*f(idx,2)
        end do
     end do
     ! Update boundaries
     call make_virtual_fine_dp(f(1,2),ilevel)

     !==============================================
     ! Compute z = Ap and store it into f(i,3)
     !==============================================
     call cmp_Ap_cg(ilevel)

     !====================================
     ! Compute p.Ap scalar product
     !====================================
     pAp=0.0d0
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           pAp=pAp+f(idx,2)*f(idx,3)
        end do
     end do
     ! Compute global sum
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(pAp,pAp_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
     pAp=pAp_all
#endif

     !====================================
     ! Compute alpha factor
     !====================================
     alpha_cg = r2/pAp

     !====================================
     ! Recurrence on x
     !====================================
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           phi(idx)=phi(idx)+alpha_cg*f(idx,2)
        end do
     end do

     !====================================
     ! Recurrence on r
     !====================================
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           f(idx,1)=f(idx,1)-alpha_cg*f(idx,3)
        end do
     end do

     ! Compute error
     error=DSQRT(r2/dble(twotondim*numbtot(1,ilevel)))
     if(iter==1)error_ini=error
     if(verbose)write(*,112)iter,error/rhs_norm,error/error_ini

  end do
  ! End main iteration loop

  if(myid==1)write(*,115)ilevel,iter,error/rhs_norm,error/error_ini
  if(iter >= itermax)then
     if(myid==1)write(*,*)'Poisson failed to converge...'
  end if

  ! Update boundaries
  call make_virtual_fine_dp(phi(1),ilevel)

111 format('   Entering phi_fine_cg for level ',I2)
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x))
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

end subroutine phi_fine_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in f(i,1).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip,nx_loc
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(dp)::dx2,fourpi,scale,oneoversix,fact
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right
  real(dp),dimension(1:nvector),save::residu

  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
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
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do
     
     ! Interpolate potential from upper level
     do idim=1,ndim
        call interpol_phi(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_phi(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do

     ! Loop over cells
     do ind=1,twotondim
        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=phi(igridn(i,ig1)+ih1)
              else
                 phig(i,idim)=phi_left(i,id1,idim)
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=phi(igridn(i,ig2)+ih2)
              else
                 phid(i,idim)=phi_right(i,id2,idim)
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute residual using 6 neighbors potential
        do i=1,ngrid
           residu(i)=phi(ind_cell(i))
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)-oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        do i=1,ngrid
           residu(i)=residu(i)+fact*(rho(ind_cell(i))-rho_tot)
        end do

        ! Store results in f(i,1)
        do i=1,ngrid
           f(ind_cell(i),1)=residu(i)
        end do

        ! Store results in f(i,2)
        do i=1,ngrid
           f(ind_cell(i),2)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cmp_residual_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Ap_cg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into f(i,3).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(dp)::oneoversix
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector),save::residu

  ! Set constants
  oneoversix=1.0D0/dble(twondim)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
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
           igridn(i,2*idim-1)=son(nbor(ind_grid(i),2*idim-1))
           igridn(i,2*idim  )=son(nbor(ind_grid(i),2*idim  ))
        end do
     end do
     
     ! Loop over cells
     do ind=1,twotondim

        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=f(igridn(i,ig1)+ih1,2)
              else
                 phig(i,idim)=0.
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=f(igridn(i,ig2)+ih2,2)
              else
                 phid(i,idim)=0.
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute Ap using neighbors potential
        do i=1,ngrid
           residu(i)=-f(ind_cell(i),2)
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)+oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        ! Store results in f(i,3)
        do i=1,ngrid
           f(ind_cell(i),3)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cmp_Ap_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_initial_phi(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !
  !
  !
  integer::igrid,ncache,i,ngrid,ind,iskip,idim,ibound
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,ind_cell_father
  real(dp),dimension(1:nvector,1:twotondim),save::phi_int

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
 
     if(ilevel==1)then
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              phi(ind_cell(i))=0.0d0
           end do
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=0.0
              end do
           end do
        end do
        ! End loop over cells
     else
        ! Compute father cell index
        do i=1,ngrid
           ind_cell_father(i)=father(ind_grid(i))
        end do
        
        ! Interpolate
        call interpol_phi(ind_cell_father,phi_int,ngrid,ilevel,icount)
        
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              phi(ind_cell(i))=phi_int(i,ind)
           end do
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=0.0
              end do
           end do
        end do
        ! End loop over cells
     end if

  end do
  ! End loop over grids

end subroutine make_initial_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_multipole_phi(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !
  !
  !
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,dx_loc,scale,fourpi,boxlen2,eps,r2
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector),save::rr,pp
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ndim),save::ff


  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  fourpi=4.D0*ACOS(-1.0D0)
  boxlen2=boxlen**2
  eps=dx_loc

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ! Gather nvector grids
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
        
        if(simple_boundary)then
           ! Compute cell center in code units
           do idim=1,ndim
               do i=1,ngrid
                  xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
               end do
           end do

           ! Rescale position from code units to user units
           rr(1:ngrid)=0.0d0
           do idim=1,ndim
               do i=1,ngrid
                  xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                  rr(i)=rr(i)+(xx(i,idim)-multipole(idim+1)/multipole(1))**2
               end do
           end do

           do i=1,ngrid
               rr(i)=max(eps,sqrt(rr(i)))       ! Cutoff
           end do

           if(ngrid>0) call phi_ana(rr,pp,ngrid)

           ! Scatter variables
           do i=1,ngrid
               phi(ind_cell(i))=pp(i)/scale
           end do

        else
           do i=1,ngrid
               phi(ind_cell(i))=0d0
           end do
        endif
        
        ! End loop over cells
     end do
     
  end do
  ! End loop over grids

end subroutine make_multipole_phi
