subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if(cooling.and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,neul=5
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekin,emag,T2min,Zsolar

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do
     
     ! Compute kinetic energy
     do i=1,nleaf
        ekin(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           ekin(i)=ekin(i)+0.5d0*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     ! Compute magnetic energy
     do i=1,nleaf
        emag(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+neul)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do

     ! Compute temperature from polytrope
     do i=1,nleaf
        T2min(i) = T2_star*(1.0+(nH(i)*scale_d/n_star)**g_star)/scale_T2
     end do

     ! Compute minimal total energy from polytrope
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/(gamma-1.0) + ekin(i) + emag(i)
     end do

     ! Update total fluid energy
     do i=1,nleaf
        uold(ind_leaf(i),neul) = T2min(i)
     end do

  end do
  ! End loop over cells

end subroutine coolfine1



