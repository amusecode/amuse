! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-coarse-level related routines
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     potential            phi            active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------


! ------------------------------------------------------------------------
! Mask restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_mask_coarse(ifinelevel,allmasked)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel
   logical, intent(out) :: allmasked

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_amr, iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1
   allmasked=.true.

   ! Loop over coarse cells of the myid active comm
   do ind_c_cell=1,twotondim
      iskip_c_amr=ncoarse+(ind_c_cell-1)*ngridmax
      iskip_c_mg =(ind_c_cell-1)*active_mg(myid,icoarselevel)%ngrid

      ! Loop over coarse grids of myid
      do igrid_c_mg=1,active_mg(myid,icoarselevel)%ngrid
         igrid_c_amr=active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr=iskip_c_amr+igrid_c_amr
         icell_c_mg =iskip_c_mg +igrid_c_mg
         igrid_f_amr=son(icell_c_amr)
         cpu_amr=cpu_map(icell_c_amr)
         if(igrid_f_amr==0) then
            ! Cell is not refined
            ngpmask      = -1.0d0
         else
            ! Cell is refined
            ! Check if son grid is in MG hierarchy
            igrid_f_mg=lookup_mg(igrid_f_amr)
            if(igrid_f_mg<=0) then
               ! Child oct is not in multigrid hierarchy
               ngpmask=-1.0d0
            else
               ! Child oct is within MG hierarchy
               ! Loop over fine cells and gather ngpmask
               ngpmask=0.0d0
               do ind_f_cell=1,twotondim
                  ! Extract fine mask value in the corresponding MG comm
                  iskip_f_mg=(ind_f_cell-1)*active_mg(cpu_amr,ifinelevel)%ngrid
                  icell_f_mg=iskip_f_mg+igrid_f_mg
                  ngpmask=ngpmask+active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,4)
               end do
               ngpmask=ngpmask/dtwotondim
            end if
         end if
         ! Store cell mask
         active_mg(myid,icoarselevel)%u(icell_c_mg,4)=ngpmask
         allmasked=allmasked .and. (ngpmask<=0.0)
      end do
   end do

end subroutine restrict_mask_coarse

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_mask_coarse_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)
         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))
         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg
         ! Stack cell volume fraction in coarse cell
         ngpmask=(1d0+active_mg(myid,ifinelevel)%u(icell_f_mg,4))/2d0/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)+ngpmask
      end do
   end do

end subroutine restrict_mask_coarse_reverse

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_coarse(ilevel)
   ! Computes the residual for pure MG levels, and stores it into active_mg(myid,ilevel)%u(:,3)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2, phi_c, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: icell_mg, iskip_mg, igrid_nbor_mg, icell_nbor_mg
   integer  :: igrid_amr, iskip_amr, cpu_nbor_amr 
   integer  :: igshift, igrid_nbor_amr

   real(dp) :: dtwondim = (twondim)

   ! Set constants
   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active_mg(myid,ilevel)%ngrid

   ! Loop over cells myid
   do ind=1,twotondim
      iskip_mg  = (ind-1)*ngrid
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids myid
      do igrid_mg=1,ngrid
         igrid_amr = active_mg(myid,ilevel)%igrid(igrid_mg)
         icell_mg = igrid_mg + iskip_mg

         phi_c = active_mg(myid,ilevel)%u(icell_mg,1)
         nb_sum=0.0d0  ! Sum of phi on neighbors

         ! SCAN FLAG TEST
         if(.not. btest(active_mg(myid,ilevel)%f(icell_mg,1),0)) then ! NO SCAN
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if
                  igrid_nbor_mg = lookup_mg(igrid_nbor_amr) 
                  ! Add up
                  icell_nbor_mg = igrid_nbor_mg + &
                      (jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid
                  nb_sum = nb_sum + &
                      active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
               end do
            end do
         else ! PERFORM SCAN
            if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0) then
               active_mg(myid,ilevel)%u(icell_mg,3)=0.0
               cycle
            end if
            do idim=1,ndim
               do inbor=1,2
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell !
                     ! Virtual phi value on unrefnd neighbor cell : -phi_c/mask_c
                     ! (simulates mask=-1.0 for the nonexistent refined cell)
                     nb_sum = nb_sum - phi_c/active_mg(myid,ilevel)%u(icell_mg,4)
                  else
                     ! Fetch neighbor cell
                     igrid_nbor_mg  = lookup_mg(igrid_nbor_amr)
                     if(igrid_nbor_mg<=0) then
                        nb_sum=nb_sum-phi_c/active_mg(myid,ilevel)%u(icell_mg,4)
                        cycle
                     end if

                     icell_nbor_mg  = igrid_nbor_mg + &
                       (jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid
                     if(active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)<=0.0) then
                        ! Neighbor cell is masked : compute its virtual phi with the mask
                        nb_sum = nb_sum + phi_c * &
                           (active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)/active_mg(myid,ilevel)%u(icell_mg,4))
                     else
                        ! Neighbor cell is active, use its true potential
                        nb_sum = nb_sum + active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
                     end if
                  end if
               end do
            end do
         end if ! END SCAN TEST

         ! Store ***MINUS THE RESIDUAL***
         active_mg(myid,ilevel)%u(icell_mg,3) = &
          -oneoverdx2*( nb_sum - dtwondim*phi_c )+active_mg(myid,ilevel)%u(icell_mg,2)
      end do
   end do

end subroutine cmp_residual_mg_coarse

! ##################################################################
! ##################################################################

subroutine cmp_uvar_norm2_coarse(ivar, ilevel, norm2)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel, ivar
   real(dp), intent(out) :: norm2

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg, icell_mg, iskip_mg

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active_mg(myid,ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_mg = (ind-1)*ngrid
      ! Loop over active grids
      do igrid_mg=1,ngrid
         icell_mg = iskip_mg + igrid_mg
         if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0 .and. ivar/=4) cycle
         norm2 = norm2 + active_mg(myid,ilevel)%u(icell_mg,ivar)**2
      end do
   end do
   norm2 = dx2*norm2
end subroutine cmp_uvar_norm2_coarse

! ##################################################################
! ##################################################################

subroutine cmp_fvar_norm2_coarse(ivar, ilevel, norm2)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel, ivar
   real(dp), intent(out) :: norm2

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg, icell_mg, iskip_mg

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active_mg(myid,ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_mg = (ind-1)*ngrid
      ! Loop over active grids
      do igrid_mg=1,ngrid
         icell_mg = iskip_mg + igrid_mg
         if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0) cycle
         norm2 = norm2 + active_mg(myid,ilevel)%f(icell_mg,ivar)**2
      end do
   end do
   norm2 = dx2*norm2
end subroutine cmp_fvar_norm2_coarse

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_coarse(ilevel,safe,redstep)
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel
   logical, intent(in) :: safe
   logical, intent(in) :: redstep

   integer, dimension(1:3,1:2,1:8) :: iii, jjj
   integer, dimension(1:3,1:4)     :: ired, iblack

   real(dp) :: dx2, nb_sum, weight
   integer  :: ngrid
   integer  :: ind, ind0, igrid_mg, idim, inbor
   integer  :: igrid_amr, cpu_nbor_amr 
   integer  :: iskip_mg, igrid_nbor_mg, icell_mg, icell_nbor_mg
   integer  :: igshift, igrid_nbor_amr
   real(dp) :: dtwondim = (twondim)

   ! Set constants
   dx2  = (0.5d0**ilevel)**2

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

   ngrid=active_mg(myid,ilevel)%ngrid

   ! Loop over cells, with red/black ordering
   do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep
      if(redstep) then
         ind = ired  (ndim,ind0)
      else
         ind = iblack(ndim,ind0)
      end if

      iskip_mg  = (ind-1)*ngrid

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active_mg(myid,ilevel)%igrid(igrid_mg)
         icell_mg  = iskip_mg  + igrid_mg

         nb_sum=0.0d0                       ! Sum of phi on neighbors
         ! Read scan flag
         if(.not. btest(active_mg(myid,ilevel)%f(icell_mg,1),0)) then
            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
            ! Those cells are active, have all their neighbors active
            ! and all neighbors are in the AMR+MG trees
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid shift
                  igshift = iii(idim,inbor,ind)
                  ! Get neighbor grid
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if
                  ! Get neighbor cpu
                  igrid_nbor_mg  = lookup_mg(igrid_nbor_amr)
                  icell_nbor_mg  = igrid_nbor_mg + &
                      (jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid
                  nb_sum = nb_sum + &
                      active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
               end do
            end do
            ! Update the potential, solving for potential on icell_amr
            active_mg(myid,ilevel)%u(icell_mg,1)=(nb_sum-dx2*active_mg(myid,ilevel)%u(icell_mg,2))/dtwondim
         else
            ! Use the finer "solve" Gauss-Seidel near boundaries,
            ! with all necessary checks

            if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0) cycle
            if(safe .and. active_mg(myid,ilevel)%u(icell_mg,4)<1.0) cycle

            weight=0.0d0                       ! Central weight for "Solve G-S"

            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid shift
                  igshift = iii(idim,inbor,ind)

                  ! Get neighbor grid
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell, set mask=-1 on nonexistent neighbor cell
                     weight = weight - 1.0d0/active_mg(myid,ilevel)%u(icell_mg,4)
                  else
                     ! Fetch neighbor cell
                     igrid_nbor_mg  = lookup_mg(igrid_nbor_amr)
                     if(igrid_nbor_mg<=0) then
                        ! No MG neighbor
                        weight = weight - 1.0d0/active_mg(myid,ilevel)%u(icell_mg,4)
                     else
                        icell_nbor_mg  = igrid_nbor_mg  + (jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid

                        if(active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)<=0.0) then
                           ! Neighbor cell is masked
                           weight = weight + &
                                 active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)/active_mg(myid,ilevel)%u(icell_mg,4)
                        else
                           ! Neighbor cell is active, increment neighbor sum
                           nb_sum = nb_sum + active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
                        end if
                     end if
                  end if
               end do
            end do
            ! Update the potential, solving for potential on icell_amr
            active_mg(myid,ilevel)%u(icell_mg,1) = (nb_sum - dx2*active_mg(myid,ilevel)%u(icell_mg,2)) &
                     / (dtwondim - weight)
         end if
      end do
   end do
end subroutine gauss_seidel_mg_coarse

! ------------------------------------------------------------------------
! Residual restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_residual_coarse(ifinelevel)
   ! Restrict coarser (MG) residual at level ifinelevel using NGP into coarser residual at level
   ! ifinelevel-1
   ! Restricted residual is stored into the RHS at the coarser level
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   real(dp) :: val, w
   integer  :: icoarselevel, cpu_amr
   integer  :: ngrid_c, ind_c, iskip_c_amr, iskip_c_mg, igrid_c_amr, icell_c_amr, icell_c_mg, igrid_c_mg
   integer  :: ind_f, igrid_f_amr, igrid_f_mg, icell_f_mg

   icoarselevel=ifinelevel-1

   ! Loop over coarse MG cells
   ngrid_c=active_mg(myid,icoarselevel)%ngrid
   do ind_c=1,twotondim
      iskip_c_amr = ncoarse + (ind_c-1)*ngridmax
      iskip_c_mg  = (ind_c-1)*ngrid_c

      do igrid_c_mg=1,ngrid_c
         igrid_c_amr = active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr = igrid_c_amr + iskip_c_amr
         cpu_amr     = cpu_map(icell_c_amr)
         icell_c_mg  = igrid_c_mg  + iskip_c_mg

         ! Get AMR child grid
         igrid_f_amr = son(icell_c_amr)
         if(igrid_f_amr==0) then
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
            cycle
         end if

         ! Get child MG grid id
         igrid_f_mg = lookup_mg(igrid_f_amr)
         if(igrid_f_mg<=0) then
            ! Son grid is not in MG hierarchy
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
            cycle
         end if

         ! Loop over child (fine MG) cells
         val = 0.0d0
         w = 0d0
         do ind_f=1,twotondim
            icell_f_mg = igrid_f_mg + (ind_f-1)*active_mg(cpu_amr,ifinelevel)%ngrid

            if (active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,4)<=0.0) cycle
            val = val + active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,3)
            w = w + 1d0
         end do
         ! Store restricted residual into RHS of coarse level
         if(w>0) then
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = val/w
         else
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0d0
         end if
      end do
   end do
end subroutine restrict_residual_coarse



! ------------------------------------------------------------------------
! Residual restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_residual_coarse_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: res
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         ! Is fine cell masked?
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Is coarse cell masked?
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         res=active_mg(myid,ifinelevel)%u(icell_f_mg,3)/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)+res
      end do
   end do

end subroutine restrict_residual_coarse_reverse

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_coarse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr, ngrid_f, istart, nbatch
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg, iskip_f_mg, icell_f_mg
   integer  :: icoarselevel, ind_c, cpu_c_amr

   real(dp) :: a, b, c, d, coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector), save                :: igrid_f_amr, icell_amr, cpu_amr
   integer,  dimension(1:nvector,1:threetondim), save  :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim), save    :: nbors_father_grids
   real(dp), dimension(1:nvector), save                :: corr

   ! Local constants
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.D0*a
   icoarselevel=ifinelevel-1

   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Loop over fine grids by vector sweeps
   ngrid_f=active_mg(myid,ifinelevel)%ngrid
   do istart=1,ngrid_f,nvector

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active_mg(myid,ifinelevel)%igrid(istart+i-1)
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))
         cpu_amr(i)  =cpu_map(icell_amr(i))
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather(icell_amr,nbors_father_cells,nbors_father_grids,nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax
         iskip_f_mg  = (ind_f-1)*ngrid_f

         do i=1,nbatch
            ! Compute fine cell indices
            icell_amr(i) = iskip_f_amr + igrid_f_amr(i)
         end do
         corr=0.0d0

         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            do i=1,nbatch
               icell_f_mg  = iskip_f_mg + istart+i-1
               if(active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,4)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
               ind_c       = (icell_c_amr-ncoarse)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               cpu_c_amr   = cpu_map(father(igrid_c_amr))
               if(igrid_c_mg<=0) cycle

               icell_c_mg  = (ind_c-1)*active_mg(cpu_c_amr,icoarselevel)%ngrid + igrid_c_mg
               corr(i)=corr(i)+coeff*active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,1)
            end do
         end do

         ! Correct potential
         do i=1,nbatch
            icell_f_mg  = iskip_f_mg + istart+i-1
            active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1) = active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1) + corr(i)
         end do

      end do
      ! End loop over cells

   end do
   ! End loop over grids
end subroutine interpolate_and_correct_coarse


! ------------------------------------------------------------------------
! Flag setting
! ------------------------------------------------------------------------

subroutine set_scan_flag_coarse(ilevel)
   use amr_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel

   integer :: ind, ngrid, scan_flag
   integer :: igrid_mg, inbor, idim, igshift
   integer :: igrid_amr, igrid_nbor_amr, cpu_nbor_amr, icell_nbor_amr

   integer :: iskip_mg, icell_mg, igrid_nbor_mg, icell_nbor_mg

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid = active_mg(myid,ilevel)%ngrid
   if(ngrid==0) return

   ! Loop over cells and set coarse SCAN flag
   do ind=1,twotondim
      iskip_mg  = (ind-1)*ngrid
      do igrid_mg=1,ngrid
         igrid_amr = active_mg(myid,ilevel)%igrid(igrid_mg)
         icell_mg  = iskip_mg  + igrid_mg

         if(active_mg(myid,ilevel)%u(icell_mg,4)==1d0) then
            scan_flag=0       ! Init flag to 'no scan needed'
            scan_flag_loop: do inbor=1,2
               do idim=1,ndim
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     scan_flag=1
                     exit scan_flag_loop
                  else
                     igrid_nbor_mg = lookup_mg(igrid_nbor_amr)
                     if(igrid_nbor_mg<=0) then
                        scan_flag=1
                        exit scan_flag_loop
                     else
                        icell_nbor_mg  = igrid_nbor_mg  + &
                                 (jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid
                        if(active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)<=0.0) then
                           scan_flag=1
                           exit scan_flag_loop
                        end if
                     end if
                  end if
               end do
            end do scan_flag_loop
         else
            scan_flag=1
         end if
         active_mg(myid,ilevel)%f(icell_mg,1)=scan_flag
      end do
   end do
end subroutine set_scan_flag_coarse
