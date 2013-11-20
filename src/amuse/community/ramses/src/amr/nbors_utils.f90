!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get3cubefather(ind_cell_father,nbors_father_cells,&
     &                    nbors_father_grids,ncell,ilevel)
  use amr_commons 
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell_father
  integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector,1:twotondim)::nbors_father_grids
  !------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input father cell. According to the refinement rule, 
  ! they should be present anytime.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  integer,dimension(1:nvector),save::ix,iy,iz,iix,iiy,iiz
  integer,dimension(1:nvector),save::pos,ind_grid_father,ind_grid_ok
  integer,dimension(1:nvector,1:threetondim),save::nbors_father_ok
  integer,dimension(1:nvector,1:twotondim),save::nbors_grids_ok
  logical::oups

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...

     oups=.false.
     do i=1,ncell
        if(ind_cell_father(i)>ncoarse)oups=.true.
     end do
     if(oups)then
        write(*,*)'get3cubefather'
        write(*,*)'oupsssss !'
        call clean_stop
     endif

     do i=1,ncell
        iz(i)=(ind_cell_father(i)-1)/nxny
     end do
     do i=1,ncell
        iy(i)=(ind_cell_father(i)-1-iz(i)*nxny)/nx
     end do
     do i=1,ncell
        ix(i)=(ind_cell_father(i)-1-iy(i)*nx-iz(i)*nxny)
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=2
     j1min=0; j1max=0
     if(ndim > 1)j1max=2
     k1min=0; k1max=0
     if(ndim > 2)k1max=2

     ! Loop over 3^ndim neighboring father cells
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(i)=iz(i)+k1-1
              if(iiz(i) < 0   )iiz(i)=nz-1
              if(iiz(i) > nz-1)iiz(i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(i)=iy(i)+j1-1
                 if(iiy(i) < 0   )iiy(i)=ny-1
                 if(iiy(i) > ny-1)iiy(i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(i)=ix(i)+i1-1
                    if(iix(i) < 0   )iix(i)=nx-1
                    if(iix(i) > nx-1)iix(i)=0
                 end do
              end if
              ind_father=1+i1+3*j1+9*k1
              do i=1,ncell
                 nbors_father_cells(i,ind_father)=1 &
                      & +iix(i) &
                      & +iiy(i)*nx &
                      & +iiz(i)*nxny
              end do
           end do
        end do
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=1
     j1min=0; j1max=0
     if(ndim > 1)j1max=1
     k1min=0; k1max=0
     if(ndim > 2)k1max=1

     ! Loop over 2^ndim neighboring father grids
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(i)=iz(i)+2*k1-1
              if(iiz(i) < 0   )iiz(i)=nz-1
              if(iiz(i) > nz-1)iiz(i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(i)=iy(i)+2*j1-1
                 if(iiy(i) < 0   )iiy(i)=ny-1
                 if(iiy(i) > ny-1)iiy(i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(i)=ix(i)+2*i1-1
                    if(iix(i) < 0   )iix(i)=nx-1
                    if(iix(i) > nx-1)iix(i)=0
                 end do
              end if
              ind_father=1+i1+2*j1+4*k1
              do i=1,ncell
                 nbors_father_grids(i,ind_father)=1 &
                      & +(iix(i)/2) &
                      & +(iiy(i)/2)*(nx/2) &
                      & +(iiz(i)/2)*(nxny/4)
              end do
           end do
        end do
     end do

  else    ! else, more complicated...
     
     ! Get father cell position in the grid
     do i=1,ncell
        pos(i)=(ind_cell_father(i)-ncoarse-1)/ngridmax+1
     end do
     ! Get father grid
     do i=1,ncell
        ind_grid_father(i)=ind_cell_father(i)-ncoarse-(pos(i)-1)*ngridmax
     end do

     ! Loop over position
     do ind=1,twotondim

        ! Select father cells that sit at position ind
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              ind_grid_ok(iok)=ind_grid_father(i)
           end if
        end do

        if(iok>0)&
        & call get3cubepos(ind_grid_ok,ind,nbors_father_ok,nbors_grids_ok,iok)

        ! Store neighboring father cells for selected cells
        do j=1,threetondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 nbors_father_cells(i,j)=nbors_father_ok(iok,j)
              end if
           end do
        end do

        ! Store neighboring father grids for selected cells
        do j=1,twotondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 nbors_father_grids(i,j)=nbors_grids_ok(iok,j)
              end if
           end do
        end do

     end do

  end if

end subroutine get3cubefather
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get3cubepos(ind_grid,ind,nbors_father_cells,nbors_father_grids,ng)
  use amr_commons
  implicit none
  integer::ng,ind
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector,1:twotondim)::nbors_father_grids
  !--------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input cell at position ind in grid ind_grid. According to 
  ! the refinements rules and since the input cell is refined, 
  ! they should be present anytime.
  !--------------------------------------------------------------------
  integer::i,j,iskip
  integer::ii,iimin,iimax
  integer::jj,jjmin,jjmax
  integer::kk,kkmin,kkmax
  integer::icell,igrid,inbor
  integer,dimension(1:8)::iii=(/1,2,1,2,1,2,1,2/)
  integer,dimension(1:8)::jjj=(/3,3,4,4,3,3,4,4/)
  integer,dimension(1:8)::kkk=(/5,5,5,5,6,6,6,6/)
  integer,dimension(1:27,1:8,1:3)::lll,mmm
  integer,dimension(1:nvector),save::ind_grid1,ind_grid2,ind_grid3
  integer,dimension(1:nvector,1:twotondim),save::nbors_grids
  
  call getindices3cube(lll,mmm)

  iimin=0; iimax=0
  if(ndim>0)iimax=1
  jjmin=0; jjmax=0
  if(ndim>1)jjmax=1
  kkmin=0; kkmax=0
  if(ndim>2)kkmax=1

  do kk=kkmin,kkmax
     do i=1,ng
        ind_grid1(i)=ind_grid(i)
     end do
     if(kk>0)then
        inbor=kkk(ind)
        do i=1,ng
           ind_grid1(i)=son(nbor(ind_grid(i),inbor))
        end do
     end if

     do jj=jjmin,jjmax
        do i=1,ng
           ind_grid2(i)=ind_grid1(i)
        end do
        if(jj>0)then
           inbor=jjj(ind)
           do i=1,ng
              ind_grid2(i)=son(nbor(ind_grid1(i),inbor))
           end do
        end if
 
        do ii=iimin,iimax
           do i=1,ng
              ind_grid3(i)=ind_grid2(i)
           end do
           if(ii>0)then
              inbor=iii(ind)
              do i=1,ng
                 ind_grid3(i)=son(nbor(ind_grid2(i),inbor))
              end do
           end if

           inbor=1+ii+2*jj+4*kk
           do i=1,ng
              nbors_grids(i,inbor)=ind_grid3(i)
           end do     

        end do
     end do
  end do     

  do j=1,twotondim
     do i=1,ng
        nbors_father_grids(i,j)=nbors_grids(i,j)
     end do
  end do

  do j=1,threetondim
     igrid=lll(j,ind,ndim)
     icell=mmm(j,ind,ndim)
     iskip=ncoarse+(icell-1)*ngridmax
     do i=1,ng
        nbors_father_cells(i,j)=iskip+nbors_grids(i,igrid)
     end do
  end do

end subroutine get3cubepos
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getindices3cube(lll,mmm)
  implicit none
  integer,dimension(1:27,1:8,1:3)::lll,mmm
  
  lll=0; mmm=0
  ! -> ndim=1
  ! @ind =1
  lll(1:3,1,1)=(/2,1,1/)
  mmm(1:3,1,1)=(/2,1,2/)
  ! @ind =2
  lll(1:3,2,1)=(/1,1,2/)
  mmm(1:3,2,1)=(/1,2,1/)

  ! -> ndim=2
  ! @ind =1
  lll(1:9,1,2)=(/4,3,3,2,1,1,2,1,1/)
  mmm(1:9,1,2)=(/4,3,4,2,1,2,4,3,4/)
  ! @ind =2
  lll(1:9,2,2)=(/3,3,4,1,1,2,1,1,2/)
  mmm(1:9,2,2)=(/3,4,3,1,2,1,3,4,3/)
  ! @ind =3
  lll(1:9,3,2)=(/2,1,1,2,1,1,4,3,3/)
  mmm(1:9,3,2)=(/2,1,2,4,3,4,2,1,2/)
  ! @ind =4
  lll(1:9,4,2)=(/1,1,2,1,1,2,3,3,4/)
  mmm(1:9,4,2)=(/1,2,1,3,4,3,1,2,1/)

  ! -> ndim= 3
  ! @ind = 1
  lll(1:27,1,3)=(/8,7,7,6,5,5,6,5,5,4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1/)
  mmm(1:27,1,3)=(/8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8/)
  ! @ind = 2
  lll(1:27,2,3)=(/7,7,8,5,5,6,5,5,6,3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2/)
  mmm(1:27,2,3)=(/7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7/)
  ! @ind = 3
  lll(1:27,3,3)=(/6,5,5,6,5,5,8,7,7,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3/)
  mmm(1:27,3,3)=(/6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6/)
  ! @ind = 4
  lll(1:27,4,3)=(/5,5,6,5,5,6,7,7,8,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4/)
  mmm(1:27,4,3)=(/5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5/)
  ! @ind = 5
  lll(1:27,5,3)=(/4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,8,7,7,6,5,5,6,5,5/)
  mmm(1:27,5,3)=(/4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4/)
  ! @ind = 6
  lll(1:27,6,3)=(/3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,7,7,8,5,5,6,5,5,6/)
  mmm(1:27,6,3)=(/3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3/)
  ! @ind = 7
  lll(1:27,7,3)=(/2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3,6,5,5,6,5,5,8,7,7/)
  mmm(1:27,7,3)=(/2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2/)
  ! @ind = 8
  lll(1:27,8,3)=(/1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4,5,5,6,5,5,6,7,7,8/)
  mmm(1:27,8,3)=(/1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1/)

end subroutine getindices3cube
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborcells(igridn,ind,icelln,ng)
  use amr_commons
  implicit none
  integer::ng,ind
  integer,dimension(1:nvector,0:twondim)::igridn
  integer,dimension(1:nvector,1:twondim)::icelln
  !--------------------------------------------------------------
  ! This routine computes the index of 6-neighboring cells
  ! The user must provide igridn = index of the 6 neighboring
  ! grids and the cell's grid (see routine getnborgrids). 
  ! ind is the cell index in the grid.
  !--------------------------------------------------------------
  integer::i,in,ig,ih,iskip
  integer,dimension(1:8,1:6)::ggg,hhh

  ggg(1:8,1)=(/1,0,1,0,1,0,1,0/); hhh(1:8,1)=(/2,1,4,3,6,5,8,7/)
  ggg(1:8,2)=(/0,2,0,2,0,2,0,2/); hhh(1:8,2)=(/2,1,4,3,6,5,8,7/)
  ggg(1:8,3)=(/3,3,0,0,3,3,0,0/); hhh(1:8,3)=(/3,4,1,2,7,8,5,6/)
  ggg(1:8,4)=(/0,0,4,4,0,0,4,4/); hhh(1:8,4)=(/3,4,1,2,7,8,5,6/)
  ggg(1:8,5)=(/5,5,5,5,0,0,0,0/); hhh(1:8,5)=(/5,6,7,8,1,2,3,4/)
  ggg(1:8,6)=(/0,0,0,0,6,6,6,6/); hhh(1:8,6)=(/5,6,7,8,1,2,3,4/)

  ! Reset indices
  icelln(1:ng,1:twondim)=0
  ! Compute cell numbers
  do in=1,twondim
     ig=ggg(ind,in)
     ih=hhh(ind,in)
     iskip=ncoarse+(ih-1)*ngridmax
     do i=1,ng
        if(igridn(i,ig)>0)then
           icelln(i,in)=iskip+igridn(i,ig)
        end if
     end do
  end do

end subroutine getnborcells
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborfather(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell). 
  ! If for some reasons they don't exist, the routine returns 
  ! the neighboring father cells of the input cell.
  !-----------------------------------------------------------------
  integer::nxny,i,idim,j,iok,ind
  integer,dimension(1:3)::ibound,iskip1,iskip2
  integer,dimension(1:nvector,1:3),save::ix
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok

  nxny=nx*ny

  if(ilevel==1)then

     ibound(1)=nx-1
     iskip1(1)=1
     iskip2(1)=nx-1
     ibound(2)=ny-1
     iskip1(2)=nx
     iskip2(2)=(ny-1)*nx
     ibound(3)=nz-1
     iskip1(3)=nxny
     iskip2(3)=(nz-1)*nxny

     ! Get father cell
     do i=1,ncell
        ind_father(i,0)=ind_cell(i)
     end do

     do i=1,ncell
        ix(i,3)=(ind_father(i,0)-1)/nxny
     end do
     do i=1,ncell
        ix(i,2)=(ind_father(i,0)-1-ix(i,3)*nxny)/nx
     end do
     do i=1,ncell
        ix(i,1)=(ind_father(i,0)-1-ix(i,2)*nx-ix(i,3)*nxny)
     end do

     do idim=1,ndim
        do i=1,ncell
           if(ix(i,idim)>0)then
              ind_father(i,2*idim-1)=ind_father(i,0)-iskip1(idim)
           else
              ind_father(i,2*idim-1)=ind_father(i,0)+iskip2(idim)
           end if
        end do
        do i=1,ncell
           if(ix(i,idim)<ibound(idim))then
              ind_father(i,2*idim)=ind_father(i,0)+iskip1(idim)
           else
              ind_father(i,2*idim)=ind_father(i,0)-iskip2(idim)
           end if
        end do
     end do

  else

     ! Get father cell
     do i=1,ncell
        ind_father(i,0)=ind_cell(i)
     end do

     ! Get father cell position in the grid
     do i=1,ncell
        pos(i)=(ind_father(i,0)-ncoarse-1)/ngridmax+1
     end do
     
     ! Get father grid
     do i=1,ncell
        ind_grid_father(i)=ind_father(i,0)-ncoarse-(pos(i)-1)*ngridmax
     end do

     ! Get neighboring father grids
     call getnborgrids(ind_grid_father,igridn,ncell)

     ! Loop over position
     do ind=1,twotondim

        ! Select father cells that sit at position ind
        do j=0,twondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 igridn_ok(iok,j)=igridn(i,j)
              end if
           end do
        end do

        ! Get neighboring cells for selected cells
        if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)

        ! Update neighboring father cells for selected cells
        do j=1,twondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 if(icelln_ok(iok,j)>0)then
                    ind_father(i,j)=icelln_ok(iok,j)
                 else
                    ind_father(i,j)=nbor(ind_grid_father(i),j)
                 end if
              end if
           end do
        end do

     end do

  end if

end subroutine getnborfather
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborgrids(igrid,igridn,ngrid)
  use amr_commons
  implicit none
  integer::ngrid
  integer,dimension(1:nvector)::igrid
  integer,dimension(1:nvector,0:twondim)::igridn
  !---------------------------------------------------------
  ! This routine computes the index of the 6 neighboring 
  ! grids for grid igrid(:). The index for the central 
  ! grid is stored in igridn(:,0). If for some reasons
  ! the neighboring grids don't exist, then igridn(:,j) = 0.
  !---------------------------------------------------------
  integer::i,j

  ! Store central grid
  do i=1,ngrid
     igridn(i,0)=igrid(i)
  end do
  ! Store neighboring grids
  do j=1,twondim
     do i=1,ngrid
        igridn(i,j)=son(nbor(igrid(i),j))
     end do
  end do
    
end subroutine getnborgrids
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborgrids_check(igrid,igridn,ngrid)
  use amr_commons
  implicit none
  integer::ngrid
  integer,dimension(1:nvector)::igrid
  integer,dimension(1:nvector,0:twondim)::igridn
  !---------------------------------------------------------
  ! This routine does EXACTLY the same as getnborgrids
  ! but it checks if the neighbor of igrid exists in order
  ! to avoid son(0) leading to crash.
  !---------------------------------------------------------
  integer::i,j

  ! Store central grid
  do i=1,ngrid
     igridn(i,0)=igrid(i)
  end do
  ! Store neighboring grids
  do j=1,twondim
     do i=1,ngrid
        if (nbor(igrid(i),j)>0)igridn(i,j)=son(nbor(igrid(i),j))
     end do
  end do
    
end subroutine getnborgrids_check


