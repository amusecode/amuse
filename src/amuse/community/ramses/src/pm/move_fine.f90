subroutine move_fine(ilevel)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel. 
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,info
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Update particles position and velocity
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     if(npart1>0)then        
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <---- Very important !!!
           next_part=nextp(ipart)
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig   
           if(ip==nvector)then
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
              ip=0
              ig=0
           end if
           ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)

111 format('   Entering move_fine for level ',I2)

end subroutine move_fine
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move1(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine computes the force on each particle by
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,length,dx_loc,scale,vol_loc,r2
  ! Grid-based arrays
  integer ,dimension(1:nvector),save::father_cell
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp,dd,dg
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
  vol_loc=dx_loc**3

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do
  
  ! Gather neighboring father cells (should be present anytime !)
  do i=1,ng
     father_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,&
       & ng,ilevel)

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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in move'
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

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=dd(j,idim)
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
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
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)   
     end if
  end do
#endif
        
  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
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
  
  ! Gather 3-force
  ff(1:np,1:ndim)=0.0D0
  if(tracer.and.hydro)then
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+uold(indp(j,ind),idim+1)*vol(j,ind)
           end do
        end do
     end do
  endif
  if(poisson)then
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
           end do
        end do
#ifdef OUTPUT_PARTICLE_POTENTIAL
        do j=1,np
           ptcl_phi(ind_part(j)) = phi(indp(j,ind))
        end do
#endif
     end do
  endif

  ! Update velocity
  do idim=1,ndim
     if(static.or.tracer)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dtnew(ilevel)
        end do
     endif
  end do

  ! For sink cloud particle only
  if(sink)then
     ! Overwrite cloud particle velocity with sink velocity
     do idim=1,ndim
        do j=1,np
           isink=-idp(ind_part(j))
           if(isink>0)then
              new_vp(j,idim)=vsnew(isink,idim,ilevel)
           end if
        end do
     end do
  end if

  ! Store velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! Update position
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)
        end do
     else
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)+new_vp(j,idim)*dtnew(ilevel)
        end do
     endif
  end do
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do

end subroutine move1
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
