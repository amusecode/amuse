! KNOWN BUGS:
!  * modifying flag1 in bisection causes load balancing to malfunction
! POSSIBLE IMPROVEMENTS
!  * cell cost is fixed at 1 for now

module bisection
   use amr_parameters
   use amr_commons
   
   implicit none

contains

   subroutine cmp_bisection_cpumap(x,c,nn)
      ! This routine takes a ndim x nvector array as input, representing nvector space points,
      ! and returns the nvector-array of integers corresponding to the matching CPU id in
      ! the domain decomposition.

      ! Array of input coordinates
      real(dp), intent(in), dimension(:,:) :: x
      ! Array of cpu ids to output
      integer, intent(out), dimension(:)   :: c
      integer, intent(in) :: nn
      
      integer :: p, dir, id, cur, half

!      if(verbose) print *, 'entering cmp_bisection_cpumap'

      ! Loop on input points
      do p=1,nn
         ! Begin splitting along the first coordinate
         dir=1
         ! Go down binary tree starting from root
         cur=bisec_root
         do while(bisec_next(cur,1)>0) ! Keep exploring tree downwards till cur is a leaf
            ! Choose relevant half
            half=1; if(x(p,dir)>bisec_wall(cur))half=2
            ! Next node in the tree, in the matching branch
            cur=bisec_next(cur,half)
            ! Next direction
            dir=dir+1; if(dir>ndim)dir=1
         end do
         ! cur should be a leaf by now
         ! Save point cpu id into the output array
         c(p)=bisec_indx(cur)
      end do
 !     if(verbose) print *, 'done with cmp_bisection_cpumap'
   end subroutine cmp_bisection_cpumap


   ! MAIN BISECTION CREATION/UPDATING ROUTINE
   subroutine build_bisection(update)
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif

      logical, intent(in) :: update

      ! Tree-wide variables (needed between levels)
      integer,  dimension(1:nbinodes) :: tmp_imin, tmp_imax
      integer,  dimension(1:nbinodes) :: tmp_load
      real(dp), dimension(1:nbinodes,1:ndim) :: tmp_bxmin, tmp_bxmax

      ! Level-wide variables (needed within one level)
      logical,  dimension(1:nbileafnodes) :: skip
      real(dp), dimension(1:nbileafnodes) :: u_limit, l_limit
      real(dp), dimension(1:nbileafnodes) :: last_wall, best_wall, best_score, walls

      integer,  dimension(1:nbileafnodes) :: load1, myload, totload
      integer,  dimension(1:nbileafnodes) :: lncpu1, lncpu2

      logical :: all_skip, start_bisec
      integer :: load2, mytmp, tottmp
      real(dp) :: scale, nload1, nload2, score
      real(dp) :: mean, var, stdev
      real(dp), dimension(1:ndim) :: xmin, xmax
  
      integer :: nc, dir, i, lvl, ierr, iter
      integer :: lncpu, cpuid, child1, child2
      integer :: cur_levelstart, cur_cell

      scale=boxlen/dble(icoarse_max-icoarse_min+1)

      if(verbose) print *,'entering build_bisection with update = ',update

      ! TREE INIT
      bisec_root=1; bisec_next=0

      tmp_imin=0; tmp_imax=0
      tmp_imin(1)=1; tmp_imax(1)=ncpu
      tmp_bxmin(1,:)=0.0; tmp_bxmax(1,:)=scale

      l_limit(1)=0.0; u_limit(1)=1.0
      cur_levelstart=1; nc=0; dir=1; cur_cell=0

      if(update) then
         ! init histogram
         call init_bisection_histogram

         ! build top-level histogram in x dir to get total load for current cpu
         call build_bisection_histogram(0,dir)

         ! compute total load for comp. box
         mytmp=bisec_hist(1,bisec_nres)   ! total load for current cpu
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(mytmp,tottmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
         tottmp=mytmp
#endif
         tmp_load(1)=tottmp
      end if

      ! Loop through levels
      lvl=0
      level_loop: do
         
         ! Number of cells for the current level
         nc  = 2**lvl

         ! Rebuild histograms if needed, level 0 is already done by init
         if(update .and. lvl>0 .and. lvl<nbilevelmax) call build_bisection_histogram(lvl,dir)

         ! WALL-FINDING
         start_bisec=.true.
         skip=.false.

         iter=0
         dichotomy_loop: do ! main dichotomy loop, for levels 0..nbilevelmax-1
            if(lvl==nbilevelmax) exit    ! No need to bisect anything at very last level
            iter=iter+1
            
            all_skip=.true. ! init at true

            ! This loop sets the skip flag in some special cases
            do i=1,nc
               cur_cell = cur_levelstart + (i-1)      ! cell id
               if (tmp_imax(cur_cell)==0) then
                  ! this is an empty slot left by a leaf cell on an upper level
                  skip(i) = .true.
                  cycle
               end if

               ! calc number of left and right cpus
               lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
               if (lncpu==1) then
                  ! leaf cell, consider it done
                  skip(i) = .true.
                  cycle
               end if

               lncpu1(i) = lncpu/2
               lncpu2(i) = lncpu - lncpu1(i)

               all_skip = all_skip .and. skip(i)
            end do
            ! Check if dichotomy is over
            if (all_skip) exit

            ! skip flag set, start the init
            ! treat new tree creation separately
            build_from_scratch: if (.not. update) then
            do i=1,nc
               cur_cell = cur_levelstart + (i-1)
               lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
               if(skip(i)) cycle

               bisec_wall(cur_cell) = round_to_bisec_res( ( tmp_bxmin(cur_cell,dir)*lncpu2(i) &
                                                         + tmp_bxmax(cur_cell,dir)*lncpu1(i) )/lncpu )
               if(bisec_wall(cur_cell)==tmp_bxmin(cur_cell,dir) .or. bisec_wall(cur_cell)==tmp_bxmax(cur_cell,dir)) then
                  if(myid==1) print *,"Problem in bisection tree creation : insufficient resolution"
#ifndef WITHOUTMPI
                  call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
                  stop
               end if
            end do
            ! don't stay in the dichotomy loop
            exit dichotomy_loop
            end if build_from_scratch


            ! tree update, two cases : 1. init, 2. bisection
            if(start_bisec) then
               do i=1,nc
                  cur_cell = cur_levelstart + (i-1)
                  lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
                  if(skip(i)) cycle
                  ! check whether wall position is compatible with bounding box
                  if( bisec_wall(cur_cell)<=tmp_bxmin(cur_cell,dir) .or. bisec_wall(cur_cell)>=tmp_bxmax(cur_cell,dir) ) then
                     bisec_wall(cur_cell) = round_to_bisec_res( 0.5 * (tmp_bxmin(cur_cell,dir) + &
                                                                          tmp_bxmax(cur_cell,dir) ) )
                  end if

                  ! get local load for left subcell knowing current wall pos
                  myload(i) = bisec_hist( i , floor(bisec_wall(cur_cell)/bisec_res)+1 )
               end do
               ! sum the local left loads from every cpu into load1
#ifndef WITHOUTMPI
               call MPI_ALLREDUCE(myload,load1,nbileafnodes,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
               load1=myload
#endif

               do i=1,nc
                  cur_cell = cur_levelstart + (i-1)
                  if(skip(i)) cycle
                  
                  ! init best wall and init score to a terrible value
                  best_wall(i) = bisec_wall(cur_cell)
#if NPRE==4
                  best_score(i) = huge(1.0e0)
#endif
#if NPRE==8
                  best_score(i) = huge(1.0d0)
#endif

                  ! setup domain limits
                  l_limit(i) = tmp_bxmin(cur_cell, dir)
                  u_limit(i) = tmp_bxmax(cur_cell, dir)
               end do
               start_bisec=.false.
            else
               ! not starting new dichotomic stuff
               do i=1,nc
                  cur_cell = cur_levelstart + (i-1)
                  if(skip(i)) cycle

                  ! retrieve differential load
                  myload(i) = abs( bisec_hist( i , floor(max(bisec_wall(cur_cell),last_wall(i))/bisec_res)+1 ) &
                              - bisec_hist( i , floor(min(bisec_wall(cur_cell),last_wall(i))/bisec_res)+1 ) )
               end do

               ! sum up all differential loads into totload
#ifndef WITHOUTMPI
               call MPI_ALLREDUCE(myload,totload,nbileafnodes,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
               totload=myload
#endif

               do i=1,nc
                  ! transfer differential load
                  cur_cell = cur_levelstart + (i-1)
                  if(skip(i)) cycle
                  
                  if(bisec_wall(cur_cell)>last_wall(i)) then
                     load1(i) = load1(i) + totload(i)
                  else
                     load1(i) = load1(i) - totload(i)
                  end if

               end do
            end if

            do i=1,nc
               cur_cell = cur_levelstart + (i-1)
               lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
               if(skip(i)) cycle

               ! compute imbalance
               load2 = tmp_load(cur_cell)-load1(i)
               nload1 = dble(load1(i))/dble(lncpu1(i)); nload2 = dble(load2)/dble(lncpu2(i))  ! normalized loads
               score = abs(nload1-nload2)/(nload1+nload2)

               ! tolerance met ?
               if(score < bisec_tol) then
                  skip(i)=.true.
                  cycle
               end if

               ! tolerance is not met, proceed with bisection.
               ! if wall is the best one yet, store it
               if (score < best_score(i)) then
                   best_score(i) = score
                   best_wall(i)  = bisec_wall(cur_cell)
               end if
               ! compute new wall position for next bisection step
               if(nload1>nload2) then
                  ! move left
                  u_limit(i) = bisec_wall(cur_cell)
               else
                  ! move right
                  l_limit(i) = bisec_wall(cur_cell)
               end if

               ! wall pos for next step
               last_wall(i) = bisec_wall(cur_cell)
               bisec_wall(cur_cell) = round_to_bisec_res( 0.5 * (u_limit(i) + l_limit(i)) )
               ! check if we're at resolution limit for next bisection step
               if( abs(bisec_wall(cur_cell)-u_limit(i))<0.5*bisec_res &
                       .or. abs(bisec_wall(cur_cell)-l_limit(i))<0.5*bisec_res )  then
                  ! restore best wall, mark cell done and loop
                  bisec_wall(cur_cell) = best_wall(i)
                  skip(i) = .true.
                  cycle
               end if
            end do

         end do dichotomy_loop

         ! CHILDREN CREATION AND LEAF PROCESSING
         ! this is done at every level, including the very last one (for leaf processing)
         walls=0.0
         children_and_leaves: do i=1,nc
            cur_cell = cur_levelstart + (i-1)
            if (tmp_imax(cur_cell)==0) cycle

            ! Is current cell a leaf?
            if (tmp_imin(cur_cell)==tmp_imax(cur_cell)) then
               cpuid=tmp_imin(cur_cell)
               ! save cpu id
               bisec_indx(cur_cell)=cpuid
               ! save cpu bound box.
               ! for a first computation (update=false) this goes into bisec_cpubox_min
               ! for update=true, this goes into bisec_cpubox_min2, as the old boxes are still
               ! needed for load balancing (virtual_boundaries.f90)
               if(update) then
                  bisec_cpubox_min2(cpuid,:)=tmp_bxmin(cur_cell,:)
                  bisec_cpubox_max2(cpuid,:)=tmp_bxmax(cur_cell,:)
               else
                  bisec_cpubox_min(cpuid,:)=tmp_bxmin(cur_cell,:)
                  bisec_cpubox_max(cpuid,:)=tmp_bxmax(cur_cell,:)
               end if
               ! save cpu workload
               bisec_cpu_load(cpuid)=tmp_load(cur_cell)
               ! make sure has no child
               bisec_next(cur_cell,:)=0
               ! skip to next cpu
               cycle
            end if

            ! node ids of the two children
            child1 = cur_levelstart + nc + 2*i-2
            child2 = cur_levelstart + nc + 2*i-1
            ! create links in the tree
            bisec_next(cur_cell,1)=child1
            bisec_next(cur_cell,2)=child2
            ! store bounding boxes
            tmp_bxmin(child1,:)=tmp_bxmin(cur_cell,:)
            tmp_bxmax(child1,:)=tmp_bxmax(cur_cell,:)
            tmp_bxmin(child2,:)=tmp_bxmin(cur_cell,:)
            tmp_bxmax(child2,:)=tmp_bxmax(cur_cell,:)
            
            tmp_bxmax(child1,dir)=bisec_wall(cur_cell)
            tmp_bxmin(child2,dir)=bisec_wall(cur_cell)
            ! store index ranges
            lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
            tmp_imin(child1)=tmp_imin(cur_cell); tmp_imax(child1)=tmp_imin(cur_cell)+lncpu1(i)-1
            tmp_imin(child2)=tmp_imin(cur_cell)+lncpu1(i); tmp_imax(child2)=tmp_imax(cur_cell)
            ! store load of subcells
            tmp_load(child1)=load1(i)
            tmp_load(child2)=tmp_load(cur_cell)-load1(i)
            ! store walls for histogram updating
            walls(i)=bisec_wall(cur_cell)
         end do children_and_leaves


         if(lvl==nbilevelmax) exit    ! terminate level loop here if at bottom level

         ! split and sort for next histogram computation
         if(update) call splitsort_bisection_histogram(lvl,dir,walls)

         ! LEVEL CHANGE

         ! Advance level start cursor
         cur_levelstart = cur_levelstart + nc
         
         ! Alternate direction
         dir=dir+1; if(dir>ndim) dir=1

         ! Next level
         lvl = lvl + 1

      end do level_loop


      mean  = sum(dble(bisec_cpu_load))/ncpu
      var   = sum(dble(bisec_cpu_load)*dble(bisec_cpu_load))/ncpu
      stdev = sqrt(var-mean*mean)

      if(verbose .and. update) then
      print *,"Load balancing report ..."
      print *,"   Average CPU load     :", mean
      print *,"   Standard deviation   :", stdev
      print *,"   Balancing accuracy   :", stdev/mean
      print *,"   Requested tolerance  :", bisec_tol
      end if

      if(verbose) print *,'done with build_bisection'

   end subroutine build_bisection

   
   function round_to_bisec_res(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y

      y = dble(nint(x/bisec_res))*bisec_res
   end function round_to_bisec_res
    


   subroutine splitsort_bisection_histogram(lev,dir,walls)
      ! rearranges ind_to_part and updates ind_min and ind_max according to a split of domains in
      ! direction dir at wall positions walls
      ! this effectively creates two subdomains ready for histogram building

      real(dp), intent(in), dimension(:) :: walls
      integer, intent(in) :: dir, lev
   
      integer :: i, nc, lmost, rmost, tmp, nxny
      real(dp) :: xc_lmost, subcell_c, dx, scale

      integer :: ix,iy,iz
      integer :: icell, igrid, isubcell, nx_loc
      integer, dimension(1:3) :: iarray, icoarse_array

      if(verbose) print *,'entering splitsort_bisection_histogram'

      iarray=0; icoarse_array=(/ icoarse_min, jcoarse_min, kcoarse_min /)
      nx_loc=icoarse_max-icoarse_min+1; scale=boxlen/dble(nx_loc)

      new_hist_bounds=0

      nc=2**lev; nxny=nx*ny
      do i=1,nc
         lmost=bisec_hist_bounds(i); rmost=bisec_hist_bounds(i+1)-1       ! note the -1
         do while(rmost-lmost>0)
            ! compute the dir coordinate of center of leftmost cell
            ! corresponding cell id
            icell=bisec_ind_cell(lmost)
            dx = 0.5d0**cell_level(icell)
            if(icell<=ncoarse) then
               ! extract ix,iy,iz for coarse cells
               iz=(icell-1)/nxny
               iy=((icell-1)-iz*nxny)/nx
               ix=((icell-1)-iy*nx-iz*nxny)

               iarray = (/ ix, iy, iz /)
               xc_lmost = scale*(dble(iarray(dir))+0.5d0-dble(icoarse_array(dir))) ! TODO : check
            else
               ! refined cell : extract grid id and subcell id
               isubcell = ((icell-ncoarse)/ngridmax)+1
               igrid = icell - ncoarse - ngridmax*(isubcell-1)
               ! compute subcell center shift
               iz=(isubcell-1)/4
               iy=(isubcell-1-4*iz)/2
               ix=(isubcell-1-2*iy-4*iz)
   
               iarray = (/ ix, iy, iz /)
               subcell_c = (dble(iarray(dir))-0.5d0)*dx - dble(icoarse_array(dir)) 
               xc_lmost = scale*( xg(igrid,dir) + subcell_c )
            end if

            if(xc_lmost<walls(i)) then    ! NOTE : strict inequality
               lmost=lmost+1
            else
               ! swap lmost and rmost cells
               tmp=bisec_ind_cell(lmost); bisec_ind_cell(lmost)=bisec_ind_cell(rmost); bisec_ind_cell(rmost)=tmp
               rmost=rmost-1
            end if
         end do
         ! NOTE : rmost==lmost by now
         ! new histogram splitting bounds
         new_hist_bounds(2*i-1) = bisec_hist_bounds(i)
         new_hist_bounds(2*i  ) = lmost
         new_hist_bounds(2*i+1) = bisec_hist_bounds(i+1)
      end do
      bisec_hist_bounds=new_hist_bounds
   end subroutine splitsort_bisection_histogram


   subroutine init_bisection_histogram()
      ! This sets up bisec_ind_cell, bisec_hist_bounds ready to start a level-0 hist build
      ! by looping over all AMR cells

#ifndef WITHOUTMPI
      include 'mpif.h'
#endif

      integer::igrid,ncache,ngrid,ierr
      integer::ilevel,i,ind,idim

      integer::nc,ibcell,p,slot

      integer::nx_loc
      integer::icpu,ncell,ncell_loc
      integer::nxny,ix,iy,iz,iskip
  
      integer,dimension(1:nvector),save::ind_grid,ind_cell

      real(dp)::dx,scale
      real(dp),dimension(1:twotondim,1:3)::xc
      
      if(verbose) print *,'entering init_bisection_histogram'

      ! Local constants
      nxny=nx*ny
      nx_loc=icoarse_max-icoarse_min+1
      scale=boxlen/dble(nx_loc)

      ! Init ncell for current subroutine

      ! Init bisec_ind_cell (partitioned index space -> cell id mapping)
      bisec_ind_cell=0
      ncell=0
      ncell_loc=1
      dx=1.0*scale
      do iz=0,nz-1
      do iy=0,ny-1
      do ix=0,nx-1
         ind=1+ix+iy*nx+iz*nxny
         if(cpu_map(ind)==myid.and.son(ind)==0)then
            ncell=ncell+1
            flag1(ncell)=ind
            bisec_ind_cell(ncell)=ind
            cell_level(ncell)=0     ! 0 for coarse levels
         end if
      end do
      end do
      end do

      ! LOOP OVER ALL CPU AMR CELLS
      ! Loop over levels
      level_loop: do ilevel=1,nlevelmax
         ! Loop over cpus
         cpu_loop: do icpu=1,ncpu
               if(icpu==myid)then
               ncache=active(ilevel)%ngrid
            else
               ncache=reception(icpu,ilevel)%ngrid
            end if
            ! Loop over grids by vector sweeps
            grid_loop: do igrid=1,ncache,nvector
               ! Gather nvector grids
               ngrid=MIN(nvector,ncache-igrid+1)
               if(icpu==myid)then
                  do i=1,ngrid
                  ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
               end do
               else
               do i=1,ngrid
                  ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
               end do
               end if
               ! Loop over cells
               cell_loop: do ind=1,twotondim
                  iskip=ncoarse+(ind-1)*ngridmax
                  ! setup ind_grid
                  do i=1,ngrid
                     ind_cell(i)=ind_grid(i)+iskip
                  end do
                  ncell_loc=0
                  do i=1,ngrid
                     if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
                        cell_level(ind_cell(i))=ilevel

                        ncell=ncell+1
                        ncell_loc=ncell_loc+1
                        flag1(ncell)=ind_cell(i)
                        
                        bisec_ind_cell(ncell)=ind_cell(i)
                     end if
                  end do
               end do cell_loop
            end do grid_loop
         end do cpu_loop
      end do level_loop
      ! End loop over levels
      ! Ok, bisec_ind_cell is good, init the bound arrays for a level-0 histogram
      bisec_hist_bounds=0
      ! only one region (region id=1) at level 0
      bisec_hist_bounds(1)=1; bisec_hist_bounds(2)=ncell+1

   end subroutine



   subroutine build_bisection_histogram(lev,dir)
      ! build the 2**lev histograms of level lev, for a split in direction dir

#ifndef WITHOUTMPI
      include 'mpif.h'
#endif
      integer, intent(in) :: lev, dir
      integer :: ix,iy,iz
      integer :: nc, nx_loc, nxny
      integer :: i, ibicell, icell, igrid, isubcell, cell_cost
      integer, dimension(1:3) :: iarray, icoarse_array
      integer :: ierr

      real(dp) :: subcell_c, cell_coord, dx, scale
      integer  :: cell_slot

      nc=2**lev; nxny=nx*ny
      iarray=0; icoarse_array=(/ icoarse_min, jcoarse_min, kcoarse_min /)

      nx_loc=icoarse_max-icoarse_min+1; scale=boxlen/dble(nx_loc)

      ! reset histograms before adding up loads
      bisec_hist=0

      ! This is the histogram-building loop on bisection cells
      bisection_domains_loop: do ibicell=1,nc
      
         ! loop over all cells of the ibicell domain
         do i=bisec_hist_bounds(ibicell),bisec_hist_bounds(ibicell+1)-1
            ! corresponding cell id
            icell=bisec_ind_cell(i)
            dx = 0.5d0**cell_level(icell)

            if(icell<=ncoarse) then
               ! extract ix,iy,iz for coarse cells
               iz=(icell-1)/nxny
               iy=((icell-1)-iz*nxny)/nx
               ix=((icell-1)-iy*nx-iz*nxny)

               iarray = (/ ix, iy, iz /)
               cell_coord = scale* ( dble(iarray(dir))-0.5d0 ) * dx
            else
               ! refined cell : extract grid id and subcell id
               isubcell = ((icell-ncoarse)/ngridmax)+1
               igrid = icell - ncoarse - ngridmax*(isubcell-1)
               ! compute subcell center shift
               iz=(isubcell-1)/4
               iy=(isubcell-1-4*iz)/2
               ix=(isubcell-1-2*iy-4*iz)
   
               iarray = (/ ix, iy, iz /)
               subcell_c = (dble(iarray(dir))-0.5d0)*dx - dble(icoarse_array(dir)) 
               ! compute cell center
               cell_coord = scale*( xg(igrid,dir) + subcell_c )
            end if

            cell_slot = floor(cell_coord/bisec_res)+1

            ! NOTE : CELL COST FIXED TO 1
            cell_cost = 1
            ! cell_cost = 2**(cell_level(icell)-1)
            ! stack up the cell in the histogram
            bisec_hist(ibicell,cell_slot) = bisec_hist(ibicell,cell_slot)+cell_cost
         end do

         ! cumulate histogram
         do i=2,bisec_nres
            bisec_hist(ibicell,i) = bisec_hist(ibicell,i) + bisec_hist(ibicell,i-1)
         end do
      end do bisection_domains_loop

   end subroutine build_bisection_histogram
end module bisection

