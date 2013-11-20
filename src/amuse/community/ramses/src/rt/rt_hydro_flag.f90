! RT patch: Only change here is to possibly do refinement if rt_err_grad_n
!           or rt_err_grad_xH namelist parameters are set
! ________________________________________________________________________

subroutine rt_hydro_flag(ilevel)
  use amr_commons
  use rt_hydro_commons
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
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nrtvar),save::uug,uum,uud
  ! ----------------------------------------------------------------------
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(.not. rt_refine .or. aexp.lt.rt_refine_aexp) return

  ! Rescaling factors
  dx=0.5d0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  if( rt_err_grad_n==-1.0) return

  ! Loop over active grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring offsets
     call getnborgrids(ind_grid,igridn,ngrid)

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

        ! Gather neighboring cells
        call getnborcells(igridn,ind,indn,ngrid)

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do j=1,twondim
           do i=1,ngrid
              if(indn(i,j)==0)then
                 indn(i,j)=nbor(ind_grid(i),j)
              end if
           end do
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather hydro variables
           do ivar=1,nrtvar
              do i=1,ngrid
                 uug(i,ivar)=rtuold(indn(i,2*idim-1),ivar)
                 uum(i,ivar)=rtuold(ind_cell(i     ),ivar)
                 uud(i,ivar)=rtuold(indn(i,2*idim  ),ivar)
              end do
           end do
           call rt_hydro_refine(uug,uum,uud,ok,ngrid)
        end do
     
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

end subroutine rt_hydro_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################

