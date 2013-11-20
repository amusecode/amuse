!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine full_multigrid(icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  !
  integer::ilevel,icount
  logical::multigrid=.false.

  ! Exit routine if no self-gravity
  if(gravity_type>0)return

  ! Compute rho for all coarser levels using restriction operator
  do ilevel=levelmin-1,1,-1
     call restriction_fine(ilevel,multigrid)
  end do

  ! Compute potential and acceleration for all coarser levels
  do ilevel=1,levelmin-1
     call multigrid_coarse(ilevel,icount)
     call force_fine(ilevel,icount)
   end do

   ! Compute potential at levelmin
   call multigrid_coarse(levelmin,icount)

end subroutine full_multigrid
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine multigrid_coarse(ilevel,icount)
  use hydro_commons
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !--------------------------------------------------------
  ! Multigrid Poisson solver using Gauss-Seidel smoother
  ! with Red-Black ordering.
  !--------------------------------------------------------
  integer::i,idim,info,ind,iter,iterj,itermax,niter_jacobi,iskip,ibound,nx_loc
  logical::multigrid=.true.,redstep=.true.,blackstep=.false.
  real(kind=8)::dx2,oneoversix,fourpi,scale,fact,error_ini,floor,tms
  real(kind=8)::error,rhs_norm,error_all,rhs_norm_all,prec

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  
  ! Set local constants
  dx2=(0.5d0**ilevel)**2
  oneoversix=1.0D0/dble(twondim)
  niter_jacobi=2

  !-----------------------------------------------------------
  ! Interpolate potential from coarse to fine as a first guess
  ! and initialize array f(:,1:ndim) to zero
  !-----------------------------------------------------------
  call make_initial_phi(ilevel,icount)

  ! Divide by 4PI
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        phi(active(ilevel)%igrid(i)+iskip)=phi(active(ilevel)%igrid(i)+iskip)/fourpi
     end do
  end do
  ! Update boundaries for phi
  call make_virtual_fine_dp(phi(1),ilevel)
  ! Update physical boundaries for phi
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           phi(boundary(ibound,ilevel)%igrid(i)+iskip)=0d0
        end do
     end do
  end do

  ! Substract rho_tot to rho
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rho(active(ilevel)%igrid(i)+iskip)=rho(active(ilevel)%igrid(i)+iskip)-rho_tot
     end do
  end do

  !-----------------------------
  ! Compute right-hand side norm
  !-----------------------------
  rhs_norm=0.0d0; rhs_norm_all=0.0d0
  fact=(oneoversix*dx2)**2/dble(twotondim*numbtot(1,ilevel))
#ifndef NPRE
  prec=1d-10
#else
#if NPRE==4
  prec=1d-4
#else
  prec=1d-10
#endif
#endif
  floor=prec*sqrt(fact*dble(twotondim*numbtot(1,ilevel)))*rho_tot
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rhs_norm=rhs_norm+fact*rho(active(ilevel)%igrid(i)+iskip)*rho(active(ilevel)%igrid(i)+iskip)
     end do
  end do
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
#endif
  rhs_norm=SQRT(rhs_norm)
  
  !-------------------------------
  ! Upward Gauss-Seidel iterations
  !-------------------------------
  do iterj=1,niter_jacobi
     call gauss_seidel(ilevel,redstep)
     call make_virtual_fine_dp(phi(1),ilevel)
     call gauss_seidel(ilevel,blackstep)
     call make_virtual_fine_dp(phi(1),ilevel)
  end do

  !------------------------------------------
  ! Compute residual and store it in f(i,1)
  !------------------------------------------
  call cmp_residual_mg(ilevel)

  !----------------------------------------
  ! Compute residual norm
  !----------------------------------------
  error=0.0d0; error_all=0.0d0
  fact=(oneoversix*dx2)**2/dble(twotondim*numbtot(1,ilevel))
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        error=error+fact*f(active(ilevel)%igrid(i)+iskip,1)*f(active(ilevel)%igrid(i)+iskip,1)
     end do
  end do
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(error,error_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  error=error_all
#endif
  error=SQRT(error)

  !--------------------
  ! Main iteration loop
  !--------------------
  error_ini=error
  iter=0; itermax=100
  if(debug.and.myid==1)then
     write(*,*)'rhs_norm =',rhs_norm
     write(*,*)'error_ini=',error_ini
     write(*,*)'floor=    ',floor
  endif
  do while(error>epsilon*(error_ini+floor).and.iter<itermax)
  iter=iter+1
  
  !--------------------------------------
  ! Get correction from upper level solve
  !--------------------------------------
  if(ilevel>1)then
     call restriction_fine(ilevel-1,multigrid)    ! Apply restriction operator 
     call multigrid_iterator(ilevel-1)           ! Recursive call to multigrid
     call prolong(ilevel)! Correct current potential with prolongated solution
     call make_virtual_fine_dp(phi(1),ilevel)      ! Update boundaries for phi
  end if
  
  !---------------------------------
  ! Downward Gauss-Seidel iterations
  !---------------------------------
  do iterj=1,niter_jacobi
     call gauss_seidel(ilevel,redstep)
     call make_virtual_fine_dp(phi(1),ilevel)
     call gauss_seidel(ilevel,blackstep)
     call make_virtual_fine_dp(phi(1),ilevel)
  end do

  !------------------------------------------
  ! Compute residual and store it in f(i,1)
  !------------------------------------------
  call cmp_residual_mg(ilevel)

  !----------------------------------------
  ! Compute residual norm
  !----------------------------------------
  error=0.0d0; error_all=0.0d0
  fact=(oneoversix*dx2)**2/dble(twotondim*numbtot(1,ilevel))
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        error=error+fact*f(active(ilevel)%igrid(i)+iskip,1)*f(active(ilevel)%igrid(i)+iskip,1)
     end do
  end do
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(error,error_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  error=error_all
#endif
  error=SQRT(error)

  if(verbose)write(*,112)iter,error/(rhs_norm+floor),error/(error_ini+floor)
  if(debug.and.myid==1)write(*,112)iter,error/(rhs_norm+floor),error/(error_ini+floor)
  end do
  !---------------------------
  ! End main iteration loop
  !---------------------------
  if(myid==1)write(*,115)ilevel,iter,error/(rhs_norm+floor),error/(error_ini+floor)
  if(iter>=itermax)then
     if(myid==1)write(*,*)'Poisson failed to converge...'
  end if

  ! Multiply by 4PI
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        phi(active(ilevel)%igrid(i)+iskip)=phi(active(ilevel)%igrid(i)+iskip)*fourpi
     end do
  end do
  ! Update boundaries for phi
  call make_virtual_fine_dp(phi(1),ilevel)
  ! Update physical boundaries for phi
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           phi(boundary(ibound,ilevel)%igrid(i)+iskip)=0d0
        end do
     end do
  end do
  ! Add rho_tot to rho
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rho(active(ilevel)%igrid(i)+iskip)=rho(active(ilevel)%igrid(i)+iskip)+rho_tot
     end do
  end do

111 format('   Entering multigrid_coarse for level ',I2)
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x))
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))


end subroutine multigrid_coarse
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine multigrid_iterator(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------
  ! Multigrid Poisson solver using Gauss-Seidel smoother
  ! with Red-Black ordering
  !--------------------------------------------------------
  integer::i,ind,iter,niter_jacobi,iskip
  logical::multigrid=.true.,redstep=.true.,blackstep=.false.
  real(kind=8)::dx2,oneoversix,fact,tms
 
  if(numbtot(1,ilevel)==0)return

  ! Set local constants
  dx2=(0.5d0**ilevel)**2
  oneoversix=1.0D0/dble(twondim)
  niter_jacobi=3

  !---------------------------------------------
  ! Compute first guess as the diagonal solution
  !---------------------------------------------
  fact=oneoversix*dx2
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        phi(active(ilevel)%igrid(i)+iskip)=-fact*rho(active(ilevel)%igrid(i)+iskip)
     end do
  end do
  ! Update boundaries for phi
  call make_virtual_fine_dp(phi(1),ilevel)

  !-------------------------------
  ! Upward Gauss-Seidel iterations
  !-------------------------------
  do iter=1,niter_jacobi
     call gauss_seidel(ilevel,redstep)
     call make_virtual_fine_dp(phi(1),ilevel)
     call gauss_seidel(ilevel,blackstep)
     call make_virtual_fine_dp(phi(1),ilevel)
  end do

  !------------------------------------------
  ! Compute residual and store it in f(i,1)
  !------------------------------------------
  call cmp_residual_mg(ilevel)

  !--------------------------------------
  ! Get correction from upper level solve
  !--------------------------------------
  if(ilevel>1)then
     call restriction_fine(ilevel-1,multigrid)    ! Apply restriction operator
     call multigrid_iterator(ilevel-1)           ! Recursive call to multigrid
     call prolong(ilevel)! Correct current potential with prolongated solution
     call make_virtual_fine_dp(phi(1),ilevel)      ! Update boundaries for phi
  end if     

  !---------------------------------
  ! Downward Gauss-Seidel iterations
  !---------------------------------
  do iter=1,niter_jacobi
     call gauss_seidel(ilevel,redstep)
     call make_virtual_fine_dp(phi(1),ilevel)
     call gauss_seidel(ilevel,blackstep)
     call make_virtual_fine_dp(phi(1),ilevel)
  end do

end subroutine multigrid_iterator
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine gauss_seidel(ilevel,redstep)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel
  logical::redstep
  !------------------------------------------------------------------
  ! This routine computes one relaxation sweep 
  ! for one Gauss-Seidel iteration.
  !------------------------------------------------------------------
  integer::i,ind0,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(kind=8)::oneoversix,dx,dx2
  integer,dimension(1:3,1:4)::ired,iblack
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector),save::residu

  ! Set constants
  dx=0.5d0**ilevel
  dx2=dx*dx
  oneoversix=1.0D0/dble(twondim)

  ired  (1,1:4)=(/1,0,0,0/)
  iblack(1,1:4)=(/2,0,0,0/)
  ired  (2,1:4)=(/1,4,0,0/)
  iblack(2,1:4)=(/2,3,0,0/)
  ired  (3,1:4)=(/1,4,6,7/)
  iblack(3,1:4)=(/2,3,5,8/)
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
     
     ! Loop over red or black cells
     do ind0=1,twotondim/2
        if(redstep)then
           ind=ired  (ndim,ind0)        
        else
           ind=iblack(ndim,ind0)
        end if
        
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              phig(i,idim)=phi(igridn(i,ig1)+ih1)
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              phid(i,idim)=phi(igridn(i,ig2)+ih2)
           end do
        end do

        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute new potential using old neighbors potential
        do i=1,ngrid
           residu(i)=0.0d0
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)+oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        do i=1,ngrid
           residu(i)=residu(i)-oneoversix*dx2*rho(ind_cell(i))
        end do
        do i=1,ngrid
           phi(ind_cell(i))=residu(i)
        end do

     end do
     ! Loop over cells

  end do
  ! Loop over grids
  
end subroutine gauss_seidel
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_mg(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the residual r = b - A x
  ! and stores it in f(i,1)
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(kind=8)::oneoversix,dx,dx2
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(kind=8),dimension(1:nvector),save::residu

  ! Set constants
  dx=0.5d0**ilevel
  dx2=dx*dx
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
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              phig(i,idim)=phi(igridn(i,ig1)+ih1)
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              phid(i,idim)=phi(igridn(i,ig2)+ih2)
           end do
        end do
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Compute residual
        do i=1,ngrid
           residu(i)=phi(ind_cell(i))
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)-oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        do i=1,ngrid
           residu(i)=residu(i)+oneoversix*dx2*rho(ind_cell(i))
        end do
        do i=1,ngrid
           f(ind_cell(i),1)=residu(i)/oneoversix/dx2
        end do
     end do
     ! End loop over cells

  end do
  ! End loop over grids
  
end subroutine cmp_residual_mg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine restriction_fine(ilevel,multigrid)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  logical::multigrid
  !-------------------------------------------------------------------
  ! This routine compute array rho (source term for Poisson equation)
  ! by first reseting array rho to zero, then 
  ! by affecting the gas density to leaf cells, and finally
  ! by performing a restriction operation for split cells.
  ! For pure particle runs, the restriction is not necessary and the
  ! routine only set rho to zero. On the other hand, for the Multigrid
  ! solver, the restriction is necessary in any case.
  !-------------------------------------------------------------------
  integer ::ind,i,icpu,ncache,igrid,ngrid,iskip,info,ibound,nx_loc
  integer ::idim,nleaf,ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid, ind_cell, ind_leaf
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector),save::dd
  real(kind=8)::vol,dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Initialize density field to zero
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
           rho(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
        end do
     end do
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rho(active(ilevel)%igrid(i)+iskip)=0.0D0
     end do
  end do
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           rho(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0D0
        end do
     end do
  end do
  
  ! Perform a restriction over split cells (ilevel+1)
  if(ilevel<nlevelmax)then
     ncache=active(ilevel+1)%ngrid
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel+1)%igrid(igrid+i-1)
        end do
        call restrict(ind_grid,ngrid,ilevel+1,multigrid)
     end do
  end if
  ! Update boundaries
  call make_virtual_reverse_dp(rho(1),ilevel)
  call make_virtual_fine_dp   (rho(1),ilevel)

111 format('   Entering restriction_fine for level',i2)

end subroutine restriction_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine restrict(ind_grid,ngrid,ilevel,multigrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel
  logical::multigrid
  integer,dimension(1:nvector)::ind_grid
  !
  !
  integer ,dimension(1:nvector),save::ind_cell,ind_cell_father
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector),save::new_rho

  real(dp)::a,b,c,d,coeff
  real(dp),dimension(1:8)::bbb
  integer ,dimension(1:8,1:8)::ccc

  integer::i,ind_father,ind_average,ind,iskip
  
  a = 1.0D0/4.0D0**ndim
  b = 3.0D0*a
  c = 9.0D0*a
  d = 27.D0*a
  
  bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)
  bbb=bbb/dble(twotondim)

  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)
    
  ! Compute father cell index
  do i=1,ngrid
     ind_cell(i)=father(ind_grid(i))
  end do

  ! Gather 3x3x3 neighboring parent cells
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  ! Update residual for coarse grid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Loop over relevant parent cells
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff     =bbb(ind_average)
        do i=1,ngrid
           ind_cell_father(i)=nbors_father_cells(i,ind_father)
        end do
        ! Gather rho in temporary array
        do i=1,ngrid
           new_rho(i)=rho(ind_cell_father(i))
        end do
        ! Perform CIC projection
        if(multigrid)then
           do i=1,ngrid
              new_rho(i)=new_rho(i)+coeff*f(ind_cell(i),1)
           end do
        else
           do i=1,ngrid
              new_rho(i)=new_rho(i)+coeff*rho(ind_cell(i))
           end do
        end if
        ! Update array rho
        do i=1,ngrid
           rho(ind_cell_father(i))=new_rho(i)
        end do
        
     end do
  end do

end subroutine restrict
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine prolong(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ilevel
  ! This routine updates the current solution using the correction
  ! given by the prolongated solution at the finer level.
  integer::i,ind_father,ind_average,ind,iskip,ncache,igrid,ngrid

  real(dp)::a,b,c,d,coeff
  real(dp),dimension(1:8)::bbb
  integer,dimension(1:8,1:8)::ccc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector),save::new_rho

  ! Local constants
  a = 1.0D0/4.0D0**ndim
  b = 3.0D0*a
  c = 9.0D0*a
  d = 27.D0*a
  
  bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Compute father cell index
     do i=1,ngrid
        ind_cell(i)=father(ind_grid(i))
     end do
     
     ! Gather 3x3x3 neighboring parent cells
     call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)
     
     ! Update solution for fine grid cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        new_rho=0.0d0

        ! Loop over relevant parent cells
        do ind_average=1,twotondim
           ind_father=ccc(ind_average,ind)
           coeff     =bbb(ind_average)
           do i=1,ngrid
              new_rho(i)=new_rho(i)+coeff*phi(nbors_father_cells(i,ind_father))
           end do
        end do

        ! Correct potential
        do i=1,ngrid
           phi(ind_cell(i))=phi(ind_cell(i))+new_rho(i)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine prolong




