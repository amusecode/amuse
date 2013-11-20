!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchro_hydro_fine(ilevel,dteff)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dteff
  !-------------------------------------------------------------------
  ! Update velocity  from gravitational acceleration
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
           ind_cell(i)=ind_grid(i)+iskip
        end do
        call synchydrofine1(ind_cell,ngrid,dteff)
     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering synchro_hydro_fine for level',i2)

end subroutine synchro_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchydrofine1(ind_cell,ncell,dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell
  real(dp)::dteff
  integer,dimension(1:nvector)::ind_cell
  !-------------------------------------------------------------------
  ! Gravity update for hydro variables
  !-------------------------------------------------------------------
  integer::i,idim,neul=ndim+2,nndim=ndim
  real(dp),dimension(1:nvector),save::pp

  ! Compute internal + magnetic energy
  do i=1,ncell
     pp(i)=uold(ind_cell(i),neul)
  end do
  do idim=1,nndim
     do i=1,ncell
        pp(i)=pp(i)-0.5*uold(ind_cell(i),idim+1)**2/uold(ind_cell(i),1)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do
  
  ! Update momentum
  do idim=1,ndim
     do i=1,ncell
        pp(i)=uold(ind_cell(i),idim+1)+ &
             & uold(ind_cell(i),1)*f(ind_cell(i),idim)*dteff
     end do
     do i=1,ncell
        uold(ind_cell(i),idim+1)=pp(i)
     end do
  end do
  
  ! Update total energy
  do i=1,ncell
     pp(i)=uold(ind_cell(i),neul)
  end do
  do idim=1,nndim
     do i=1,ncell
        pp(i)=pp(i)+0.5*uold(ind_cell(i),idim+1)**2/uold(ind_cell(i),1)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do
  
end subroutine synchydrofine1
!################################################################
!################################################################
!################################################################
!################################################################

