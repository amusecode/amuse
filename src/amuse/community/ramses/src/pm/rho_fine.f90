!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Particles that are not entirely in
  ! level ilevel contribute also to the level density field
  ! (boundary particles) using buffer grids.
  ! Array flag1, flag2 and phi are used as temporary work space.
  ! Array rho and cpu_map2 are stored with:
  ! - rho containing the Poisson source term
  ! - cpu_map2 containing the refinement map due to particle
  !   number density criterion (quasi Lagrangian mesh).
  !------------------------------------------------------------------
  integer::iskip,icpu,ind,i,info,nx_loc,ibound,idim
  real(dp)::dx,d_scale,scale,dx_loc,scalar
  real(dp)::d0,m_refine_loc,dx_min,vol_min,mstar,msnk,nISM,nCOM
  real(kind=8)::total,total_all,total2,total2_all,tms
  real(kind=8),dimension(2)::totals_in,totals_out
  logical::multigrid=.false.
  real(kind=8),dimension(1:ndim+1)::multipole_in,multipole_out

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  if(ilevel==levelmin)multipole=0d0

  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------
  if(ilevel==levelmin.or.icount>1)then
     do i=nlevelmax,ilevel,-1
        ! Compute mass multipole
        if(hydro)call multipole_fine(i)
        ! Perform CIC using pseudo-particle
        call cic_from_multipole(i)
        ! Update boundaries
        call make_virtual_reverse_dp(rho(1),i)
        call make_virtual_fine_dp   (rho(1),i)
     end do
  end if

  !--------------------------
  ! Initialize fields to zero
  !--------------------------
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        phi(active(ilevel)%igrid(i)+iskip)=0.0D0
     end do
     if(ilevel==cic_levelmax)then
        do i=1,active(ilevel)%ngrid
           rho_top(active(ilevel)%igrid(i)+iskip)=0.0D0
        end do
     endif
  end do
  if(cic_levelmax>0.and.ilevel>cic_levelmax)then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           rho_top(active(ilevel)%igrid(i)+iskip)=rho_top(father(active(ilevel)%igrid(i)))
           rho(active(ilevel)%igrid(i)+iskip)=rho(active(ilevel)%igrid(i)+iskip)+ &
                & rho_top(active(ilevel)%igrid(i)+iskip)
        end do
     end do
  endif

  !-------------------------------------------------------------------------
  ! Initialize "number density" field to baryon number density in array phi.
  !-------------------------------------------------------------------------
  if(m_refine(ilevel)>-1.0d0)then
     d_scale=max(mass_sph/dx_loc**ndim,smallr)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        if(hydro)then
           if(ivar_refine>0)then
              do i=1,active(ilevel)%ngrid
                 scalar=uold(active(ilevel)%igrid(i)+iskip,ivar_refine) &
                      & /uold(active(ilevel)%igrid(i)+iskip,1)
                 if(scalar>var_cut_refine)then
                    phi(active(ilevel)%igrid(i)+iskip)= &
                         & rho(active(ilevel)%igrid(i)+iskip)/d_scale
                 endif
              end do
           else
              do i=1,active(ilevel)%ngrid
                 phi(active(ilevel)%igrid(i)+iskip)= &
                      & rho(active(ilevel)%igrid(i)+iskip)/d_scale
              end do
           endif
        endif
     end do
  endif

  !-------------------------------------------------------
  ! Initialize rho and phi to zero in virtual boundaries
  !-------------------------------------------------------
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
           rho(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
           phi(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
        end do
        if(ilevel==cic_levelmax)then
           do i=1,reception(icpu,ilevel)%ngrid
              rho_top(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
           end do
        endif
     end do
  end do

  !---------------------------------------------------------
  ! Compute particle contribution to density field
  !---------------------------------------------------------
  ! Compute density due to current level particles
  if(pic)then
     call rho_from_current_level(ilevel)
  end if
  ! Update boudaries
  call make_virtual_reverse_dp(rho(1),ilevel)
  call make_virtual_fine_dp   (rho(1),ilevel)
  if(ilevel==cic_levelmax)then
     call make_virtual_reverse_dp(rho_top(1),ilevel)
  endif
  if(cic_levelmax>0.and.ilevel>=cic_levelmax)then
     call make_virtual_fine_dp   (rho_top(1),ilevel)
  endif
  if(m_refine(ilevel)>-1.0d0)then
     call make_virtual_reverse_dp(phi(1),ilevel)
     call make_virtual_fine_dp   (phi(1),ilevel)
  endif

  !--------------------------------------------------------------
  ! Compute multipole contribution from all cpus and set rho_tot
  !--------------------------------------------------------------
#ifndef WITHOUTMPI
  if(ilevel==levelmin)then
     multipole_in=multipole
     call MPI_ALLREDUCE(multipole_in,multipole_out,ndim+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     multipole=multipole_out
  endif
#endif
  if(nboundary==0)then
     rho_tot=multipole(1)/scale**ndim
     if(debug)write(*,*)'rho_average=',rho_tot
  else
     rho_tot=0d0
  endif

  !----------------------------------------------------
  ! Reset rho and phi in physical boundaries
  !----------------------------------------------------
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           phi(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
           rho(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
        end do
     end do
  end do

  !-----------------------------------------
  ! Compute quasi Lagrangian refinement map
  !-----------------------------------------
  if(m_refine(ilevel)>-1.0d0)then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(phi(active(ilevel)%igrid(i)+iskip)>=m_refine(ilevel))then
              cpu_map2(active(ilevel)%igrid(i)+iskip)=1
           else
              cpu_map2(active(ilevel)%igrid(i)+iskip)=0
           end if
        end do
     end do
     ! Update boundaries
     call make_virtual_fine_int(cpu_map2(1),ilevel)
  end if

!!$  do ind=1,twotondim
!!$     iskip=ncoarse+(ind-1)*ngridmax
!!$     do i=1,active(ilevel)%ngrid
!!$        print*,rho(active(ilevel)%igrid(i)+iskip),rho_tot
!!$     end do
!!$  end do
  

111 format('   Entering rho_fine for level ',I2)
  
end subroutine rho_fine
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_from_current_level(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme from particles that are not entirely in
  ! level ilevel (boundary particles).
  ! Arrays flag1 and flag2 are used as temporary work space.
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu
  integer::i,ig,ip,npart1
  real(dp)::dx

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim),save::x0
    
  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0   
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           
           ! Loop over particles
           do jpart=1,npart1
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              if(ip==nvector)then
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do i=1,ig
                       x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
                    end do
                 end do
                 do i=1,ig
                    ind_cell(i)=father(ind_grid(i))
                 end do
                 call cic_amr(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=nextp(ipart)  ! Go to next particle
           end do
           ! End loop over particles
           
        end if

        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ig
              x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
           end do
        end do
        do i=1,ig
           ind_cell(i)=father(ind_grid(i))
        end do
        call cic_amr(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
     end if

  end do
  ! End loop over cpus

end subroutine rho_from_current_level
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine cic_amr(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: mass_sph
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim)::x0
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  logical::error
  integer::j,ind,idim,nx_loc
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm
  real(dp),dimension(1:nvector),save::ttt=0d0
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim


  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Gather particle mass
  do j=1,np
     mmm(j)=mp(ind_part(j))
  end do

  if(ilevel==levelmin)then
     do j=1,np
        multipole(1)=multipole(1)+mp(ind_part(j))
     end do
     do idim=1,ndim
        do j=1,np
           multipole(idim+1)=multipole(idim+1)+mp(ind_part(j))*xp(ind_part(j),idim)
        end do
     end do
  end if

  ! Gather particle birth epoch
  if(star)then
     do j=1,np
        ttt(j)=tp(ind_part(j))
     end do
  endif

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in cic'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=dd(j,idim)
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif
        
  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icg(j,idim)=ig(j,idim)-2*igg(j,idim)
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
     icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif

  ! Compute parent cell adress
  do ind=1,twotondim
     do j=1,np
        indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
     end do
  end do

  ! Update mass density and number density fields
  do ind=1,twotondim

     do j=1,np
        ok(j)=igrid(j,ind)>0
     end do

     do j=1,np
        vol2(j)=mmm(j)*vol(j,ind)/vol_loc
     end do

     if(cic_levelmax==0.or.ilevel<=cic_levelmax)then
        do j=1,np
           if(ok(j))then
              rho(indp(j,ind))=rho(indp(j,ind))+vol2(j)
           end if
        end do
     else if(ilevel>cic_levelmax)then
        do j=1,np
           if(ok(j).and.ttt(j).ne.0d0)then
              rho(indp(j,ind))=rho(indp(j,ind))+vol2(j)
           end if
        end do
     endif

     if(ilevel==cic_levelmax)then
        do j=1,np
           if(ok(j).and.ttt(j)==0d0)then
              rho_top(indp(j,ind))=rho_top(indp(j,ind))+vol2(j)
           end if
        end do
     endif

     do j=1,np
        vol2(j)=vol(j,ind)
     end do

     ! Remove test particles for static runs
     if(static)then
        do j=1,np
           ok(j)=ok(j).and.mmm(j)>0.0
        end do
     endif

     ! Remove massive dark matter particle
     if(mass_cut_refine>0.0)then
        do j=1,np
           if(ttt(j)==0d0)then
              ok(j)=ok(j).and.mmm(j)<mass_cut_refine
           endif
        end do
     endif

     ! For low mass baryon particles
     if(star)then
        do j=1,np
           if(ttt(j).ne.0.0)then
              vol2(j)=vol2(j)*mmm(j)/mass_sph
           endif
        end do
     endif

     if(cic_levelmax==0.or.ilevel<cic_levelmax)then
        do j=1,np
           if(ok(j))then
              phi(indp(j,ind))=phi(indp(j,ind))+vol2(j)
           end if
        end do
     else if(ilevel>=cic_levelmax)then
        do j=1,np
           if(ok(j).and.ttt(j).ne.0d0)then
              phi(indp(j,ind))=phi(indp(j,ind))+vol2(j)
           end if
        end do
     endif

  end do

end subroutine cic_amr
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine multipole_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
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
  integer ::idim,nleaf,nsplit,ix,iy,iz,iskip_son,ind_son,ind_grid_son,ind_cell_son
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf,ind_split
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector),save::dd
  real(kind=8)::vol,dx,dx_loc,scale,vol_loc,mm
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Initialize fields to zero
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        unew(active(ilevel)%igrid(i)+iskip,1)=0.0D0
     end do
     do idim=1,ndim
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,idim+1)=0.0D0
        end do
     end do
  end do

  ! Compute mass multipoles in each cell
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ! Gather cell indices
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells and compute cell centers
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
              do idim=1,ndim
                 xx(nleaf,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
              end do
           end if
        end do
        
        ! Compute gas multipole for leaf cells only
        if(hydro)then
           do i=1,nleaf
              mm=max(uold(ind_leaf(i),1),smallr)*vol_loc
              unew(ind_leaf(i),1)=unew(ind_leaf(i),1)+mm
           end do
           do idim=1,ndim
              do i=1,nleaf
                 mm=max(uold(ind_leaf(i),1),smallr)*vol_loc
                 unew(ind_leaf(i),idim+1)=unew(ind_leaf(i),idim+1)+mm*xx(i,idim)
              end do
           end do
        endif

        ! Add analytical density profile for leaf cells only
        if(gravity_type < 0)then           
           ! Call user defined routine rho_ana
           call rho_ana(xx,dd,dx_loc,nleaf)
           ! Scatter results to array phi
           do i=1,nleaf
              unew(ind_leaf(i),1)=unew(ind_leaf(i),1)+dd(i)*vol_loc
           end do
           do idim=1,ndim
              do i=1,nleaf
                 mm=dd(i)*vol_loc
                 unew(ind_leaf(i),idim+1)=unew(ind_leaf(i),idim+1)+mm*xx(i,idim)
              end do
           end do           
        end if

        ! Gather split cells
        nsplit=0
        do i=1,ngrid
           if(son(ind_cell(i))>0)then
              nsplit=nsplit+1
              ind_split(nsplit)=ind_cell(i)
           end if
        end do

        ! Add children multipoles
        do ind_son=1,twotondim
           iskip_son=ncoarse+(ind_son-1)*ngridmax
           do i=1,nsplit
              ind_grid_son=son(ind_split(i))
              ind_cell_son=iskip_son+ind_grid_son
              unew(ind_split(i),1)=unew(ind_split(i),1)+unew(ind_cell_son,1)
           end do
           do idim=1,ndim
              do i=1,nsplit
                 ind_grid_son=son(ind_split(i))
                 ind_cell_son=iskip_son+ind_grid_son
                 unew(ind_split(i),idim+1)=unew(ind_split(i),idim+1)+unew(ind_cell_son,idim+1)
              end do
           end do
        end do

     end do
  enddo

  ! Update boundaries
  do idim=1,ndim+1
     call make_virtual_fine_dp(unew(1,idim),ilevel)
  end do

111 format('   Entering multipole_fine for level',i2)

end subroutine multipole_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_from_multipole(ilevel)
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
  integer ::ind,i,j,icpu,ncache,ngrid,iskip,info,ibound,nx_loc
  integer ::idim,nleaf,ix,iy,iz,igrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
  ! Reset rho in physical boundaries
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           rho(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
        end do
     end do
  end do
  
  if(hydro)then
     ! Perform a restriction over split cells (ilevel+1)
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        call cic_cell(ind_grid,ngrid,ilevel)
     end do
  end if

111 format('   Entering cic_from_multipole for level',i2)

end subroutine cic_from_multipole
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_cell(ind_grid,ngrid,ilevel)
  use amr_commons
  use poisson_commons
  use hydro_commons, ONLY: unew
  implicit none
  integer::ngrid,ilevel
  integer,dimension(1:nvector)::ind_grid
  !
  !
  integer::i,j,idim,ind_cell_son,iskip_son,np,ind_son,nx_loc,ind
  integer ,dimension(1:nvector),save::ind_cell,ind_cell_father
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector),save::new_rho
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm,ttt
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  logical::error
  
  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  np=ngrid

  ! Compute father cell index
  do i=1,ngrid
     ind_cell(i)=father(ind_grid(i))
  end do

  ! Gather 3x3x3 neighboring parent cells
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  ! Loop over grid cells
  do ind_son=1,twotondim
     iskip_son=ncoarse+(ind_son-1)*ngridmax

     ! Compute pseudo particle (centre of mass) position
     do idim=1,ndim
        do j=1,np
           ind_cell_son=iskip_son+ind_grid(j)
           x(j,idim)=unew(ind_cell_son,idim+1)/unew(ind_cell_son,1)
        end do
     end do
     
     ! Compute total multipole
     if(ilevel==levelmin)then
        do idim=1,ndim+1
           do j=1,np
              ind_cell_son=iskip_son+ind_grid(j)
              multipole(idim)=multipole(idim)+unew(ind_cell_son,idim)
           end do
        end do
     endif

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/scale+skip_loc(idim)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)-(xg(ind_grid(j),idim)-3d0*dx)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/dx
        end do
     end do
     
     ! Gather particle mass
     do j=1,np
        ind_cell_son=iskip_son+ind_grid(j)
        mmm(j)=unew(ind_cell_son,1)
     end do
     
     ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
     do idim=1,ndim
        do j=1,np
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=dd(j,idim)
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end do
     end do
     
     ! Check for illegal moves
     error=.false.
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
        end do
     end do
     if(error)then
        write(*,*)'problem in cic'
        do idim=1,ndim
           do j=1,np
              if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
                 write(*,*)x(j,1:ndim)
              endif
           end do
        end do
        stop
     end if

     ! Compute cloud volumes
#if NDIM==1
     do j=1,np
        vol(j,1)=dg(j,1)
        vol(j,2)=dd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        vol(j,1)=dg(j,1)*dg(j,2)
        vol(j,2)=dd(j,1)*dg(j,2)
        vol(j,3)=dg(j,1)*dd(j,2)
        vol(j,4)=dd(j,1)*dd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
        vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
        vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
        vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
        vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
        vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
        vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
        vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
     end do
#endif
     
     ! Compute parent grids
     do idim=1,ndim
        do j=1,np
           igg(j,idim)=ig(j,idim)/2
           igd(j,idim)=id(j,idim)/2
        end do
     end do
#if NDIM==1
     do j=1,np
        kg(j,1)=1+igg(j,1)
        kg(j,2)=1+igd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        kg(j,1)=1+igg(j,1)+3*igg(j,2)
        kg(j,2)=1+igd(j,1)+3*igg(j,2)
        kg(j,3)=1+igg(j,1)+3*igd(j,2)
        kg(j,4)=1+igd(j,1)+3*igd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
        kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
        kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
        kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
        kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
        kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
        kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
        kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
     end do
#endif
     do ind=1,twotondim
        do j=1,np
           igrid(j,ind)=son(nbors_father_cells(j,kg(j,ind)))
        end do
     end do
     
     ! Compute parent cell position
     do idim=1,ndim
        do j=1,np
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end do
     end do
#if NDIM==1
     do j=1,np
        icell(j,1)=1+icg(j,1)
        icell(j,2)=1+icd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end do
#endif
     
     ! Compute parent cell adress
     do ind=1,twotondim
        do j=1,np
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        end do
     end do
     
     ! Update mass density and number density fields
     do ind=1,twotondim
        do j=1,np
           ok(j)=igrid(j,ind)>0
        end do
        do j=1,np
           vol2(j)=mmm(j)*vol(j,ind)/vol_loc
        end do        
        do j=1,np
           if(ok(j))then
              rho(indp(j,ind))=rho(indp(j,ind))+vol2(j)
           end if
        end do
     end do
     
  end do
  ! End loop over grid cells

end subroutine cic_cell
!###########################################################
!###########################################################
!###########################################################
!###########################################################
