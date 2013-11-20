!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag
  use amr_commons
  implicit none
  integer::ilevel

  if(verbose)write(*,*)'Entering flag'
  do ilevel=nlevelmax-1,1,-1
     call flag_fine(ilevel,2)
  end do
  call flag_coarse
  if(verbose)write(*,*)'Complete flag'

end subroutine flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag_coarse
  use amr_commons
  implicit none
  !--------------------------------------------------------------
  ! This routine compute the refinement map at the coarse level.
  !--------------------------------------------------------------
  integer::ind,nxny,ix,iy,iz

  if(verbose)write(*,*)'  Entering flag_coarse'
  ! Constants
  nxny=nx*ny
  ! Reset flag1 array at coarse level
  flag1(0:ncoarse)=0
  ! Set flag1 to 1 at coarse level for inner cells only
  nflag=0
  do iz=kcoarse_min,kcoarse_max
     do iy=jcoarse_min,jcoarse_max
        do ix=icoarse_min,icoarse_max
           ind=1+ix+iy*nx+iz*nxny
           flag1(ind)=1
           nflag=nflag+1
        end do
     end do
  end do
  if(verbose)write(*,112)nflag  
  call make_virtual_coarse_int(flag1(1))
  if(simple_boundary)call make_boundary_coarse

112 format('   ==> Flag ',i6,' cells')

end subroutine flag_coarse
!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag_fine(ilevel,icount)
  use amr_commons
  implicit none
  integer::ilevel,icount
  !--------------------------------------------------------
  ! This routine builds the refinement map at level ilevel.
  !--------------------------------------------------------
  integer::iexpand

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Step 1: initialize refinement map to minimal refinement rules
  call init_flag(ilevel)
  if(verbose)write(*,*) '  ==> end step 1',nflag

  ! If ilevel < levelmin, exit routine
  if(ilevel<levelmin)return
  if(balance)return

  ! Step 2: make one cubic buffer around flagged cells,
  ! in order to enforce numerical rule.
  call smooth_fine(ilevel)
  if(verbose)write(*,*) '  ==> end step 2',nflag

  ! Step 3: if cell satisfies user-defined physical citeria,
  ! then flag cell for refinement.
  call userflag_fine(ilevel)    
  if(verbose)write(*,*) '  ==> end step 3',nflag

  ! Step 4: make nexpand cubic buffers around flagged cells.
  do iexpand=1,nexpand(ilevel)
     call smooth_fine(ilevel)
  end do
  if(verbose)write(*,*) '  ==> end step 4',nflag

  if(verbose)write(*,112)nflag

  ! In case of adaptive time step ONLY, check for refinement rules.
  if(ilevel>levelmin)then
     if(icount<nsubcycle(ilevel-1))then
        call ensure_ref_rules(ilevel)
     end if
  end if

111 format('   Entering flag_fine for level ',I2)
112 format('   ==> Flag ',i6,' cells')

end subroutine flag_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flag(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  !-------------------------------------------
  ! This routine initialize the refinement map
  ! to a minimal state in order to satisfy the
  ! refinement rules.
  !-------------------------------------------
  integer::i,ind,iskip

  ! Initialize flag1 to 0
  nflag=0
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        flag1(active(ilevel)%igrid(i)+iskip)=0
     end do
  end do

  ! If load balancing operations, flag only refined cells
  if(balance)then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip)>0)then
              flag1(active(ilevel)%igrid(i)+iskip)=1
              nflag=nflag+1
           end if
        end do
     end do
  else
     ! If cell is refined and contains a flagged son
     ! or a refined son, then flag cell for refinement.
     if(ilevel>=levelmin)then
        call test_flag(ilevel)
     else
        ! If ilevel < levelmin, set flag to 1 for all cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(ilevel)%ngrid
              flag1(active(ilevel)%igrid(i)+iskip)=1
           end do
           nflag=nflag+active(ilevel)%ngrid
        end do
     end if
  end if
  
  ! Update boundaries
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary)call make_boundary_flag(ilevel)

end subroutine init_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine test_flag(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets flag1 to 1 if cell is refined and 
  ! contains a flagged son or a refined son.
  ! This ensures that refinement rules are satisfied.
  !---------------------------------------------------------
  integer::i,ind_son,ind,iskip
  integer::iskip_son,ind_grid_son,ind_cell_son
  logical::ok

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     ! Test all refined cells
     do i=1,active(ilevel)%ngrid
        ! Gather child grid number
        ind_grid_son=son(active(ilevel)%igrid(i)+iskip)
        ! Test child if it exists
        ok=.false.
        if(ind_grid_son>0)then
           ! Loop over children cells
           do ind_son=1,twotondim
              iskip_son=ncoarse+(ind_son-1)*ngridmax
              ind_cell_son=iskip_son+ind_grid_son
              ok=(ok.or.(son  (ind_cell_son)> 0))
              ok=(ok.or.(flag1(ind_cell_son)==1))
           end do
        end if
        ! If ok, then flag1 cells.
        if(ok)then
           flag1(active(ilevel)%igrid(i)+iskip)=1
           nflag=nflag+1
        end if
     end do
  end do
  ! End loop over cells

end subroutine test_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine ensure_ref_rules(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------
  ! This routine determines if all grids at level ilevel are 
  ! surrounded by 26 neighboring grids, in order to enforce the 
  ! strict refinement rule. 
  ! Used in case of adaptive time steps only.
  !-----------------------------------------------------------------
  integer::i,ind,iskip,igrid,ngrid,ncache
  integer,dimension(1:nvector),save::ind_cell,ind_grid
  integer,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  logical,dimension(1:nvector),save::ok

  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Gather neighboring father cells (should be present anytime !)
     do i=1,ngrid
        ind_cell(i)=father(ind_grid(i))
     end do
     call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids &
          & ,ngrid,ilevel)
     
     do i=1,ngrid
        ok(i)=.true.
     end do

     do ind=1,threetondim
        do i=1,ngrid
           ind_cell(i)=nbors_father_cells(i,ind)
        end do
        do i=1,ngrid
           if(son(ind_cell(i))==0)ok(i)=.false.
        end do
     end do
     
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           if(.not.ok(i))flag1(ind_cell(i))=0
        end do
     end do

  end do

  ! Update boundaries
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary)call make_boundary_flag(ilevel)

end subroutine ensure_ref_rules 
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine userflag_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  integer::i,j,ncache,nok,ix,iy,iz,iskip
  integer::igrid,ind,idim,ngrid,ivar
  integer::nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  integer,dimension(1:nvector,1:twondim),save::indn

  logical,dimension(1:nvector),save::ok

  real(dp)::dx,dx_loc,scale
  real(dp)::d0,dx_min,vol_min,mstar,msnk,nISM,nCOM
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  
  logical::prevent_refine

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return

  ! Conversion factor from user units to cgs units                              
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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

  ! Do we prevent the whole level from refining ?
  prevent_refine=.false.
  
  ! Prevent from refining too much relative to available star mass res.      
  ! This translates into a constant physical res., except at very high z.    
  if(star)then
     ! Finest cell size                                                      
     dx_min=(0.5D0**nlevelmax)*scale
     vol_min=dx_min**ndim
     ! Typical ISM mass density from H/cc to code units                      
     nISM = n_star
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
     nISM = MAX(nCOM,nISM)
     d0   = nISM/scale_nH
     ! Star particle mass                                                    
     mstar=n_star/(scale_nH*aexp**3)*vol_min
     ! Test is designed so that nlevelmax is activated at aexp \simeq 0.8
     if(d0*(dx_loc/2.0)**ndim.lt.mstar/2d0)prevent_refine=.true.
  endif
  
  ! Prevent from refining too much relative to available sink mass res.      
  ! This translates into a constant physical res., except at very high z.    
  if(sink)then
     ! Finest cell size                                                      
     dx_min=(0.5D0**nlevelmax)*scale
     vol_min=dx_min**ndim
     ! Typical ISM mass density from H/cc to code units
     nISM = n_sink
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
     nISM = MAX(nCOM,nISM)
     d0   = nISM/scale_nH
     ! Sink particle mass                                                    
     msnk=n_sink/(scale_nH*aexp**3)*vol_min
     ! Test is designed so that nlevelmax is activated at aexp \simeq 0.8    
     if(d0*(dx_loc/2.0)**ndim.lt.msnk/2d0)prevent_refine=.true.
  endif
  
  if(prevent_refine)return

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over active grids
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

        ! Initialize refinement to false
        do i=1,ngrid
           ok(i)=.false.
        end do

        ! Apply purely local Lagrangian refinement criteria
        if(m_refine(ilevel)>-1.0d0)then
           call poisson_refine(ind_cell,ok,ngrid,ilevel)
           ! Apply geometry-based refinement criteria
           if(r_refine(ilevel)>-1.0)then
              ! Compute cell center in code units
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
              call geometry_refine(xx,ind_cell,ok,ngrid,ilevel)
           end if
        end if

        ! Count newly flagged cells
        nok=0
        do i=1,ngrid
           if(flag1(ind_cell(i))==0.and.ok(i))then
              nok=nok+1
           end if
        end do
        
        do i=1,ngrid
           if(ok(i))flag1(ind_cell(i))=1
        end do

        nflag=nflag+nok
     end do
     ! End loop over cells

  end do
  ! End loop over grids

  ! Do the same for hydro solver
  if(hydro)call hydro_flag(ilevel)

#ifdef RT
  ! Do the same for RT solver
  if(rt)call rt_hydro_flag(ilevel)
#endif

  ! Update boundaries
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary)call make_boundary_flag(ilevel)

end subroutine userflag_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine poisson_refine(ind_cell,ok,ncell,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  logical,dimension(1:nvector)::ok
  !-------------------------------------------------
  ! This routine sets flag1 to 1 if cell statisfy
  ! user-defined physical criterion for refinement.
  !-------------------------------------------------
  integer::i,nx_loc
  real(dp)::d_scale,scale,dx,dx_loc,vol_loc

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  if(poisson)then

     if(.not. init) then
        do i=1,ncell
           ok(i)=ok(i).or.(cpu_map2(ind_cell(i))==1)
        end do
     else
        if(ivar_refine==0)then
           do i=1,ncell
              ok(i)=ok(i).or.(cpu_map2(ind_cell(i))==1)
           end do
        else if(ivar_refine>0)then
           do i=1,ncell
              ok(i)=ok(i).or. &
                   & (uold(ind_cell(i),ivar_refine)/uold(ind_cell(i),1) &
                   & > var_cut_refine)
           end do
        else if(m_refine(ilevel)==0.0)then
           do i=1,ncell
              ok(i)=.true.
           end do
        endif
     endif

  else 

     if(hydro)then
        d_scale=mass_sph/vol_loc
        do i=1,ncell
           ok(i)=ok(i).or.(uold(ind_cell(i),1)>=m_refine(ilevel)*d_scale)
        end do
     endif

  end if

end subroutine poisson_refine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine geometry_refine(xx,ind_cell,ok,ncell,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:ndim)::xx
  logical ,dimension(1:nvector)::ok
  !-------------------------------------------------
  ! This routine sets flag1 to 1 if cell statisfy
  ! user-defined physical criterion for refinement.
  !-------------------------------------------------
  real(dp)::er,xr,yr,zr,rr,xn,yn,zn,r,aa,bb
  integer ::i

  ! Authorize refinement if cell lies within region,
  ! otherwise unmark cell (no refinement outside region)
  if(r_refine(ilevel)>-1.0)then
     er=exp_refine(ilevel) ! Exponent defining norm
     xr=x_refine  (ilevel) ! Region centre
     yr=y_refine  (ilevel)
     zr=z_refine  (ilevel)
     rr=r_refine  (ilevel) ! Region DIAMETER (beware !)
     aa=a_refine  (ilevel) ! Ellipticity (Y/X)
     bb=b_refine  (ilevel) ! Ellipticity (Z/X)
     do i=1,ncell
        xn=0.0d0; yn=0.0d0; zn=0.0d0
        xn=2.0d0*abs(xx(i,1)-xr)/rr
#if NDIM > 1
        yn=2.0d0*abs(xx(i,2)-yr)/(aa*rr)
#endif
#if NDIM >2
        zn=2.0d0*abs(xx(i,3)-zr)/(bb*rr)
#endif
        if(er<10)then
           r=(xn**er+yn**er+zn**er)**(1.0/er)
        else
           r=max(xn,yn,zn)
        end if
        ok(i)=ok(i).and.(r < 1.0)
     end do
  endif

end subroutine geometry_refine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine smooth_fine(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Dilatation operator.
  ! This routine makes one cell width cubic buffer around flag1 cells 
  ! at level ilevel by following these 3 steps:
  ! step 1: flag1 cells with at least 1 flag1 neighbors (if ndim > 0) 
  ! step 2: flag1 cells with at least 2 flag1 neighbors (if ndim > 1) 
  ! step 3: flag1 cells with at least 2 flag1 neighbors (if ndim > 2) 
  ! Array flag2 is used as temporary workspace.
  ! -------------------------------------------------------------------
  integer::ismooth
  integer::i,ncache,iskip,ngrid
  integer::igrid,ind
  integer,dimension(1:3)::n_nbor
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return

  n_nbor(1:3)=(/1,2,2/)
  flag1(0)=0
  ncache=active(ilevel)%ngrid
  ! Loop over steps
  do ismooth=1,ndim
     ! Initialize flag2 to 0
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              flag2(ind_cell(i))=0
           end do
        end do
     end do
     ! Count neighbors and set flag2 accordingly
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        call getnborgrids(ind_grid,igridn,ngrid)
        do ind=1,twotondim
           call count_nbors(igridn,ind,n_nbor(ismooth),ngrid)
        end do
     end do
     ! Set flag1=1 for cells with flag2=1
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              if(flag1(ind_cell(i))==1)flag2(ind_cell(i))=0
           end do
           do i=1,ngrid
              if(flag2(ind_cell(i))==1)then
                 flag1(ind_cell(i))=1
                 nflag=nflag+1
              end if
           end do
        end do
     end do
     ! Update boundaries
     call make_virtual_fine_int(flag1(1),ilevel)
     if(simple_boundary)call make_boundary_flag(ilevel)
  end do
  ! End loop over steps

end subroutine smooth_fine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine count_nbors(igridn,ind,n_nbor,nn)
  use amr_commons
  implicit none
  integer::ind,nn,n_nbor
  integer,dimension(1:nvector,0:twondim)::igridn
  !----------------------------------------------------
  ! This routine computes the number of neighbors 
  ! for cell ind in grid igridn(:,0) for which flag1=1.
  ! The user must provide the neighboring grids index
  ! stored in igridn(:,:) and the threshold n_nbor
  ! If the number of flag1 neighbors exceeds n_nbor, 
  ! then cell is marked with flag2=1
  !----------------------------------------------------
  integer::i,in,iskip
  integer,dimension(1:nvector),save::ind_cell,i_nbor
  integer,dimension(1:nvector,1:twondim),save::indn
  ! Compute cell number
  iskip=ncoarse+(ind-1)*ngridmax
  do i=1,nn
     ind_cell(i)=iskip+igridn(i,0)
  end do
  ! Gather neighbors
  call getnborcells(igridn,ind,indn,nn)
  ! Check if neighboring cell exists and count it as a flagged neighbor
  i_nbor(1:nn)=0
  do in=1,twondim
     do i=1,nn
        i_nbor(i)=i_nbor(i)+flag1(indn(i,in))
     end do
  end do
  ! flag2 cell if necessary
  do i=1,nn
     if(i_nbor(i)>=n_nbor)flag2(ind_cell(i))=1
  end do
end subroutine count_nbors
!############################################################
!############################################################
!############################################################
!############################################################
subroutine count_nbors2(igridn,ind,n_nbor,nn)
  use amr_commons
  implicit none
  integer::ind,nn,n_nbor
  integer,dimension(1:nvector,0:twondim)::igridn
  !----------------------------------------------------
  ! This routine computes the number of neighbors 
  ! for cell ind in grid igridn(:,0) for which flag2=1.
  ! The user must provide the neighboring grids index
  ! stored in igridn(:,:) and the threshold n_nbor
  ! If the number of flag2 neighbors exceeds n_nbor, 
  ! then cell is marked with flag1=1
  !----------------------------------------------------
  integer::i,in,iskip
  integer,dimension(1:nvector),save::ind_cell,i_nbor
  integer,dimension(1:nvector,1:twondim),save::indn
  ! Compute cell number
  iskip=ncoarse+(ind-1)*ngridmax
  do i=1,nn
     ind_cell(i)=iskip+igridn(i,0)
  end do
  ! Gather neighbors
  call getnborcells(igridn,ind,indn,nn)
  ! Check if neighboring cell exists and count it as a flagged neighbor
  i_nbor(1:nn)=0
  do in=1,twondim
     do i=1,nn
        i_nbor(i)=i_nbor(i)+flag2(indn(i,in))
     end do
  end do
  ! flag2 cell if necessary
  do i=1,nn
     if(i_nbor(i)>=n_nbor)flag1(ind_cell(i))=1
  end do
end subroutine count_nbors2
!############################################################
!############################################################
!############################################################
subroutine init_refmap
  use amr_commons
  implicit none

  integer::ilevel,ivar

  if(verbose)write(*,*)'Entering init_refmap'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call init_refmap_fine(ilevel)
     call make_virtual_fine_int(cpu_map2(1),ilevel)
  end do
  if(verbose)write(*,*)'Complete init_refmap'

end subroutine init_refmap
!############################################################
!############################################################
!############################################################
!############################################################
subroutine init_refmap_fine(ilevel)
  use amr_commons
  use hydro_commons
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
  real(dp)::dx,rr,vx,vy,vz,ek,ei,pp,xx1,xx2,xx3,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ,save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file1,ok_file2,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
  if(multiple)then
     filename=TRIM(initfile(ilevel))//'/dir_refmap/ic_refmap.00001'
     INQUIRE(file=filename,exist=ok_file2)
  else
     filename=TRIM(initfile(ilevel))//'/ic_refmap'
     INQUIRE(file=filename,exist=ok_file2)
  endif
  if (.not. ok_file2)then
     if(myid==1)write(*,*)'File '//TRIM(filename)//' not found'
     call clean_stop
  endif
     
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
           cpu_map2(icell)=int(init_array(i1,i2,i3))
        end do
     end do
     ! End loop over cells
  endif

  ! Deallocate initial conditions array
  if(ncache>0)deallocate(init_array)
  deallocate(init_plane) 
  
111 format('   Entering init_refmap_fine ',I2)

end subroutine init_refmap_fine
