subroutine clump_finder(create_output)
  use amr_commons
  use pm_commons, ONLY:nsink,xsink
  use poisson_commons, ONLY:phi
  use clfind_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::create_output

  !----------------------------------------------------------------------------
  ! Description of clump_finder:
  ! The clumpfinder assigns a test particle to each cell having a density above
  ! a given threshold. These particles are moved to the densest neighbors until
  ! all particles sit in a local density maximum. The particles (now knowing
  ! the peak they belong to) are moved back to their original position and all
  ! the relevant properties are computed. If a so called peak patch is
  ! considered irrelevant, it is merged to the neighbor which it is connected
  ! to through the saddle point with the highest density.
  ! Andreas Bleuler & Davide Martizzi & Romain Teyssier 
  !---------------------------------------------------------------------------- 

  integer::itest,istep,nskip,ilevel,info,icpu,nmove,nmove_all
  integer::i,j,jj,ntest,ntest_all,peak_nr
  integer,dimension(1:ncpu)::ntest_cpu,ntest_cpu_all
  integer,dimension(1:ncpu)::npeaks_per_cpu,npeaks_per_cpu_tot
  logical::all_bound,ok
  real(dp),dimension(1:nvector,1:3)::pos
  integer,dimension(1:nvector)::cell_index,cell_levl,cc
  integer::flag_form,flag_form_tot
  real(dp)::dx,dx_min,dist,scale,tff,acc_r
  real(dp)::fourpi,threepi2
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  if(verbose.and.myid==1)write(*,*)' Entering clump_finder'

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !-------------------------------------------------------------------------------
  ! count the number of test particles to be created, flag the cells, share info
  !------------------------------------------------------------------------------- 
  ntest=0
  do ilevel=levelmin,nlevelmax
     call count_test_particle(ilevel,ntest,0,1) !action 1: count and flag 
  end do
  ntest_cpu=0; ntest_cpu_all=0
  ntest_cpu(myid)=ntest
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntest_cpu(1)=ntest_cpu_all(1)
#endif
  do icpu=2,ncpu
     ntest_cpu(icpu)=ntest_cpu(icpu-1)+ntest_cpu_all(icpu)
  end do
  ntest_all=ntest_cpu(ncpu)
  if(myid==1)then
     if(ntest_all.gt.0)then
        write(*,'(" Total number of test particles=",I10)')ntest_all
     endif
  end if

  !-------------------------------------------------------------------------------
  ! Allocate arrays and create test particles
  !-------------------------------------------------------------------------------
  if (ntest>0) then
     allocate(denp(ntest),levp(ntest),iglobalp(ntest),icellp(ntest))
     denp=0.d0; levp=0; iglobalp=0; icellp=0
  endif
  itest=0
  nskip=ntest_cpu(myid)-ntest
  do ilevel=levelmin,nlevelmax
     call count_test_particle(ilevel,itest,nskip,2) !case2 : create the parts 
  end do
  do ilevel=nlevelmax,levelmin,-1
     call make_virtual_fine_int(flag2(1),ilevel)
  end do

  !-------------------------------------------------------------------------------
  ! Sort particles according to density
  !-------------------------------------------------------------------------------
  if (ntest>0) then
     allocate(testp_sort(ntest)) 
     do i=1,ntest
        denp(i)=-denp(i)
        testp_sort(i)=i
     end do
     call quick_sort_dp(denp(1),testp_sort(1),ntest) 
     deallocate(denp)
  endif

  !-------------------------------------------------------------------------------
  ! Count number of density peaks
  !-------------------------------------------------------------------------------
  npeaks=0; nmove=0
  if(ntest>0)call scan_for_peaks(ntest,npeaks,1) !case 1: count peaks
  npeaks_per_cpu=0
  npeaks_per_cpu(myid)=npeaks
  if(clinfo .and. npeaks>0)write(*,*)'n_peaks on processor number',myid,'= ',npeaks

  !----------------------------------------------------------------------------                       
  ! Share number of peaks per cpu and create a list  
  !----------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeaks_per_cpu,npeaks_per_cpu_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  npeaks_per_cpu_tot=npeaks_per_cpu
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeaks,npeaks_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  npeaks_tot=npeaks
#endif
  if (myid==1.and.npeaks_tot>0)write(*,'(" Total number of density peaks found=",I6)')npeaks_tot

  !----------------------------------------------------------------------------
  ! Determine peak-ids positions for each cpu
  !----------------------------------------------------------------------------
  peak_nr=0
  do icpu=1,myid-1
     peak_nr=peak_nr+npeaks_per_cpu_tot(icpu)
  end do

  !----------------------------------------------------------------------------
  ! flag peaks with global peak id
  !----------------------------------------------------------------------------
  nmove=0
  nskip=peak_nr
  flag2=0
  if(ntest>0)call scan_for_peaks(ntest,nskip,2) !case 2: flag peaks
  do ilevel=nlevelmax,levelmin,-1
     call make_virtual_fine_int(flag2(1),ilevel)
  end do

  !-------------------------------------------------------------------------------               
  ! main step:
  ! - order cells in descending density
  ! - get peak id from densest neighbor 
  ! - nmove is number of peak id's passed along
  ! - done when nmove=0 (for single cores, only one sweep necessary)   
  !-------------------------------------------------------------------------------
  nmove=1
  istep=0
  do while (nmove.gt.0)
     nmove=0
     nskip=peak_nr
     if(ntest>0)call scan_for_peaks(ntest,nmove,3)
     do ilevel=nlevelmax,levelmin,-1
        call make_virtual_fine_int(flag2(1),ilevel)
     end do
     istep=istep+1
#ifndef WITHOUTMPI 
     call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nmove=nmove_all
#endif   
     if(myid==1)write(*,*)"istep=",istep,"nmove=",nmove   
  end do

  ! Allocate peak-patch property arrays
  call allocate_peak_patch_arrays

  ! Compute peak-patch mass etc. and output these properties before merging 
  call compute_clump_properties(ntest) 
  if (clinfo)call write_clump_properties(.false.)


  !-------------------------------------------------------------------------------
  ! Find the saddle point densities and merge irrelevant clumps
  !-------------------------------------------------------------------------------
  if (npeaks_tot > 0)then
     !communicate across boundaries
     do ilevel=nlevelmax,levelmin,-1
        call make_virtual_fine_int(flag2(1),ilevel)
        call make_virtual_fine_dp(phi(1),ilevel)
     end do
     call saddlepoint_search(ntest) 
     call merge_clumps(ntest)


     !if all clumps need to be gravitationally bound to survive - merge again
     if (merge_unbound)then
        do while (.not. all_bound)
           call compute_clump_properties_round2(ntest,all_bound)
           call write_clump_properties(.false.)
           do j=npeaks_tot,1,-1
              if (isodens_check(j)<1.)relevance_tot(j)=1.
           end do
           call merge_clumps(ntest)
        end do
     endif

     !for sink formation, in star formation case, intersection of 4cell ball and clump is considered
     if ((.not. smbh) .and. sink .and. (.not. create_output))then
        if(myid==1) print*,'now trimming clumps'
        call trim_clumps(ntest)
        call compute_clump_properties(ntest)
     end if


     call compute_clump_properties_round2(ntest,all_bound)
     ! write properties to screen 
     call write_clump_properties(.false.)
     ! ..and if wanted to disk
     if (create_output)then
        call write_clump_properties(.true.)
        call write_clump_map(ntest)
     end if

  end if

  !------------------------------------------------------------------------------
  ! if the clumpfinder is used to produce sinks, flag all the cells which contain
  ! a relevant density peak whose peak patch doesn't yet contain a sink.
  !------------------------------------------------------------------------------
  if(sink.and.(.not. create_output))then
     allocate(occupied(1:npeaks_tot),occupied_all(1:npeaks_tot))
     occupied=0; occupied_all=0;
     ! loop over sinks and mark all clumps containing a sink
     pos=0.0
     if(myid==1 .and. clinfo)write(*,*)'looping over ',nsink,' sinks and marking their clumps'
     do j=1,nsink
        pos(1,1:3)=xsink(j,1:3)
        call cmp_cpumap(pos,cc,1)
        if (cc(1) .eq. myid)then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           if (flag2(cell_index(1))>0)then
              occupied(flag2(cell_index(1)))=1
              if(clinfo)write(*,*)'CPU # ',myid,'blocked clump # ',flag2(cell_index(1)),' for sink production because of sink # ',j
           end if
        end if
     end do


#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(occupied,occupied_all,npeaks_tot,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     occupied_all=occupied
#endif

     !------------------------------------------------------------------------------
     ! determine whether a peak patch is eligible to form a new sink.
     ! if a new sink has to be created, flag2 is set to 1 at the peak position
     !------------------------------------------------------------------------------     
     pos=0.0
     flag2=0
     allocate(form(1:npeaks_tot),form_all(1:npeaks_tot))
     form=0; form_all=0; 
     flag_form=0; flag_form_tot=0;  
     call heapsort_index(max_dens_tot,sort_index,npeaks_tot)
     do j=npeaks_tot,1,-1
        jj=sort_index(j)
        if (.not.smbh)then
           ok=.true.
           ok=ok.and.relevance_tot(jj)>0.
           ok=ok.and.occupied_all(jj)==0
           ok=ok.and.max_dens_tot(jj)>(n_sink/scale_nH)
           ok=ok.and.contracting(jj)
           ok=ok.and.Icl_dd_tot(jj)<0.
           if (ok)then
              pos(1,1:3)=peak_pos_tot(jj,1:3)
              call cmp_cpumap(pos,cc,1)
              if (cc(1) .eq. myid)then
                 call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
                 flag2(cell_index(1))=jj
                 write(*,*)'cpu ',myid,' produces a new sink for clump number ',jj
              end if
           end if
        else 
           ok=.true.
           ok=ok.and.relevance_tot(jj)>0.
           ok=ok.and.occupied_all(jj)==0
           ok=ok.and.peak_check(jj)>1.
           !ok=ok.and.ball4_check(jj)>1.
           !ok=ok.and.isodens_check(jj)>1.
           fourpi=4.0d0*ACOS(-1.0d0)
           threepi2=3.0d0*ACOS(-1.0d0)**2
           if(cosmo)fourpi=1.5d0*omega_m*aexp
           tff=sqrt(threepi2/8./fourpi/(max_dens_tot(jj)+1.0d-30))
           acc_r=clump_mass_tot4(jj)*dble(scale_d)*(dble(scale_l)**3.0)*3600.0*24.0*365.0/1.98892d33/tff/dble(scale_t)
           ok=ok.and.acc_r > 30.d0

           if (ok)then
              pos(1,1:3)=peak_pos_tot(jj,1:3)
              call cmp_cpumap(pos,cc,1)
              if (cc(1) .eq. myid)then
                 call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
                 ! Geometrical criterion 
                 if(ivar_refine>0)then
                    if(uold(cell_index(1),ivar_refine)>var_cut_refine)then
                       flag2(cell_index(1))=jj
                       form(jj)=1
                       flag_form=1
                    end if
                 else
                    flag2(cell_index(1))=jj
                    form(jj)=1
                    flag_form=1
                 end if
              end if
           end if
        end if
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(form,form_all,npeaks_tot,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     form_all=form
#endif
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(flag_form,flag_form_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     flag_form_tot=flag_form
#endif
     if(myid == 1)then
        if(flag_form_tot>0)write(*,'(135A)')'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [cm] size_y [cm] size_z [cm] |v|_CM [u.u.] rho- [H/cc] rho+ [H/cc] rho_av [H/cc] M_cl [M_sol] V_cl [AU^3] rel.  peak_check   isodens_check   clump_check '
        do j=npeaks_tot,1,-1
           jj=sort_index(j)
           if(form_all(jj) == 1)write(*,'(I6,X,I10,17(1X,1PE14.7))')jj&
                ,n_cells_tot(jj)&
                ,peak_pos_tot(jj,1),peak_pos_tot(jj,2),peak_pos_tot(jj,3)&
                ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*scale_l &
                ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*scale_l &
                ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l &
                ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+ &
                clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)*scale_l/scale_t&
                ,min_dens_tot(jj)*scale_nH,max_dens_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33&
                ,clump_vol_tot(jj)*(scale_l)**3&
                ,relevance_tot(jj)&
                ,peak_check(jj)&
                !,ball4_check(jj)&
                ,isodens_check(jj)&
                ,clump_check(jj)
        end do
     end if
     deallocate(occupied,occupied_all)
     deallocate(form,form_all)
  endif
  
  ! Deallocate test particle and peak arrays
  if (ntest>0)then
     deallocate(icellp)
     deallocate(levp)
     deallocate(testp_sort)
     deallocate(iglobalp)
  endif
  call deallocate_all
  if(create_output.and.smbh)deallocate(clump_mass_tot4)

end subroutine clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
subroutine count_test_particle(ilevel,ntot,nskip,action)
  use amr_commons
  use hydro_commons, ONLY:uold
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,ntot,nskip,action

  !----------------------------------------------------------------------
  ! Description: This routine loops over all cells above and checks wether
  ! their density lies above the threshold. If so:
  ! case 1: count the new test particles and flag the cell
  ! case 2: create the test particle
  !----------------------------------------------------------------------

  ! local constants
  real(dp)::d0
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! other variables
  integer ::ncache,ngrid
  integer ::igrid,ind,i,iskip
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  logical ,dimension(1:nvector)::ok

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return

  if(verbose .and. myid==1)write(*,*)' Entering count test particle'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Clump density threshold from H/cc to code units
  d0 = density_threshold/scale_nH
  if(cosmo)d0=d0/aexp**3

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        !checks
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0 !check if leaf cell
           ok(i)=ok(i).and.uold(ind_cell(i),1)>d0 !check density
        end do

        select case (action) 
        case (1) !count and flag
           ! Compute test particle map
           do i=1,ngrid
              flag2(ind_cell(i))=0
              if(ok(i))then
                 flag2(ind_cell(i))=1 
                 ntot=ntot+1
              endif
           end do
        case(2) !create 'testparticles'
           do i=1,ngrid
              if (ok(i))then
                 ntot=ntot+1                    ! Local test particle index
                 levp(ntot)=ilevel              ! Level
                 iglobalp(ntot)=ntot+nskip      ! Global test particle index
                 flag2(ind_cell(i))=ntot+nskip  ! Initialize flag2 to GLOBAL test particle index
                 icellp(ntot)=ind_cell(i)       ! Local cell index
                 denp(ntot)=uold(ind_cell(i),1) ! Save density values here!
              end if
           end do
        end select
     end do
  end do

end subroutine count_test_particle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine scan_for_peaks(npartt,n,action)
  use amr_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
#endif
  integer::npartt,n,action

  !----------------------------------------------------------------------
  ! vectorization of the neighborsearch for the action cases
  ! 1: count the peaks (no denser neighbor)
  ! 2: count and flag peaks with global peak index number
  ! 3: get global clump index from densest neighbor
  !----------------------------------------------------------------------

  integer::ilevel,next_level,ipart,ip
  integer,dimension(1:nvector)::ind_cell
  
  ip=0
  do ipart=1,npartt
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0
     if(ipart<npartt)next_level=levp(testp_sort(ipart+1)) !level of next particle
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(ind_cell,ip,n,ilevel,action)
        ip=0
     endif
  end do
  if (ip>0)call neighborsearch(ind_cell,ip,n,ilevel,action)

  if(verbose)write(*,*)'   Exiting scan_for_peaks',n

end subroutine scan_for_peaks
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine neighborsearch(ind_cell,np,count,ilevel,action)
  use amr_commons
  use clfind_commons, ONLY: icellp
  use hydro_commons, ONLY: uold
  implicit none
  integer::np,count,ilevel,action
  integer,dimension(1:nvector)::ind_grid,ind_cell

  !------------------------------------------------------------
  ! This routine constructs all neighboring leaf cells at levels 
  ! ilevel-1, ilevel, ilevel+1.
  ! Depending on the action case value, fuctions performing
  ! further checks for the neighbor cells are called.
  !------------------------------------------------------------

  integer::j,ind,nx_loc,i1,j1,k1,i2,j2,k2,i3,j3,k3,ix,iy,iz
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,dx_loc,scale,vol_loc
  !nvector length arrays
  integer ,dimension(1:nvector)::cell_index,cell_levl,ind_max,clump_nr,indv
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xtest
  real(dp),dimension(1:nvector)::density_max
  real(dp),dimension(1:3)::skip_loc
  logical ,dimension(1:nvector)::okpeak,ok

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

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=0; i3max=0
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=0; j3max=0
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=0; k3max=0
  if(ndim>0)then
     i1max=1; i2max=2; i3max=3
  end if
  if(ndim>1)then
     j1max=1; j2max=2; j3max=3
  end if
  if(ndim>2)then
     k1max=1; k2max=2; k3max=3
  end if

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! some preliminary action...
  do j=1,np
     indv(j)=(ind_cell(j)-ncoarse-1)/ngridmax+1 ! cell position in grid
     ind_grid(j)=ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
     density_max(j)=uold(ind_cell(j),1)*1.0001 !get cell density (1.0001 probably not necessary, just a safety measure)
     ind_max(j)=ind_cell(j) !save cell index   
     if (action.ge.4)clump_nr(j)=flag2(ind_cell(j)) ! save clump number
  end do
  
  !initialze logical array
  okpeak=.true.

  !================================
  ! generate neighbors level ilevel-1
  !================================
  if(ilevel>levelmin)then
     ! Generate 2x2x2 neighboring cells at level ilevel-1
     do k1=k1min,k1max
        do j1=j1min,j1max
           do i1=i1min,i1max
              ok=.false.
              do j=1,np
                 xtest(j,1)=(xg(ind_grid(j),1)+2*xc(indv(j),1)-skip_loc(1))*scale+(2*i1-1)*dx_loc
                 xtest(j,2)=(xg(ind_grid(j),2)+2*xc(indv(j),2)-skip_loc(2))*scale+(2*j1-1)*dx_loc
                 xtest(j,3)=(xg(ind_grid(j),3)+2*xc(indv(j),3)-skip_loc(3))*scale+(2*k1-1)*dx_loc
              end do
              call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
              do j=1,np 
                 !check wether neighbor is in a leaf cell at the right level
                 if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel-1))ok(j)=.true.
              end do     
              !check those neighbors
              if (action<4) call peakcheck(cell_index,okpeak,ok,density_max,ind_max,np)
              if (action==4) call saddlecheck(ind_cell,cell_index,clump_nr,ok,np)
              if (action==5) call phi_ref_check(ind_cell,cell_index,clump_nr,ok,np)
           end do
        end do
     end do
  endif

  !================================
  ! generate neighbors at level ilevel
  !================================
  ! Generate 3x3x3 neighboring cells at level ilevel
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           ok=.false.
           do j=1,np
              xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc
              xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc
              xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc
           end do
           call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
           do j=1,np
              !check wether neighbor is in a leaf cell at the right level
              if(son(cell_index(j))==0.and.cell_levl(j)==ilevel)ok(j)=.true.
           end do
           !check those neighbors
           if (action<4)call peakcheck(cell_index,okpeak,ok,density_max,ind_max,np)
           if (action==4)call saddlecheck(ind_cell,cell_index,clump_nr,ok,np)
           if (action==5) call phi_ref_check(ind_cell,cell_index,clump_nr,ok,np)
        end do
     end do
  end do

  !===================================
  ! generate neighbors at level ilevel+1
  !====================================
  if(ilevel<nlevelmax)then
     ! Generate 4x4x4 neighboring cells at level ilevel+1
     do k3=k3min,k3max
        do j3=j3min,j3max
           do i3=i3min,i3max
              ok=.false.
              do j=1,np
                 xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0
                 xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0
                 xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0
              end do
              call get_cell_index(cell_index,cell_levl,xtest,ilevel+1,np)
              do j=1,np
                 !check wether neighbor is in a leaf cell at the right level
                 if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel+1))ok(j)=.true.
              end do
              !check those neighbors
              if (action<4)call peakcheck(cell_index,okpeak,ok,density_max,ind_max,np)
              if (action==4)call saddlecheck(ind_cell,cell_index,clump_nr,ok,np)
              if (action==5) call phi_ref_check(ind_cell,cell_index,clump_nr,ok,np)
           end do
        end do
     end do
  endif


  !===================================
  ! choose action for different cases
  !====================================
  select case (action)
  case (1)   ! Count peaks  
     do j=1,np
        if(okpeak(j))count=count+1
     end do
  
  case (2)   ! Initialize flag2 to peak global index
     do j=1,np
        if(okpeak(j))then 
           count=count+1
           flag2(ind_cell(j))=count
        end if
     end do
     
  case (3) ! Propagate flag2
     do j=1,np
        if(flag2(ind_cell(j)).ne.flag2(ind_max(j)))count=count+1
        flag2(ind_cell(j))=flag2(ind_max(j))
     end do
  end select

end subroutine neighborsearch
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine peakcheck(cell_index,okpeak,ok,density_max,ind_max,np)
  use amr_commons, ONLY:flag2,nvector,dp
  use hydro_commons, ONLY: uold
  implicit none

  !small routine to check wether neighbor is denser or not

  logical,dimension(1:nvector)::ok,okpeak
  integer,dimension(1:nvector)::cell_index,ind_max
  real(dp),dimension(1:nvector)::density_max
  integer::np,j

  do j=1,np
     !check if neighboring cell is denser
     ok(j)=ok(j).and.uold(cell_index(j),1)>density_max(j)
  end do
  do j=1,np
     if(ok(j))then !so if there is a denser neighbor
        okpeak(j)=.false. !no peak
        density_max(j)=uold(cell_index(j),1) !change densest neighbor dens
        ind_max(j)=cell_index(j) !change densest neighbor index
     endif
  end do

end subroutine peakcheck
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine saddlecheck(ind_cell,cell_index,clump_nr,ok,np)
  use amr_commons, ONLY:flag2,nvector,dp
  use clfind_commons
  use hydro_commons, ONLY: uold
  implicit none

  !small routine to check wether neighbor is connected through new densest saddle

  logical,dimension(1:nvector)::ok
  integer,dimension(1:nvector)::cell_index,clump_nr,ind_cell,neigh_cl
  real(dp),dimension(1:nvector)::av_dens
  integer::np,j
  integer::imat,flaggo
  integer::index,index_prev,index_next,imatp,imatn

  do j=1,np
     neigh_cl(j)=flag2(cell_index(j))!nuber of clump the neighboring cell is in 
  end do
  do j=1,np
     ok(j)=ok(j).and. neigh_cl(j)/=0 !neighboring cell is in a clump
     ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) !neighboring cell is in another clump
     av_dens(j)=(uold(cell_index(j),1)+uold(ind_cell(j),1))*0.5 !average density of cell and neighbor cell
  end do
  do j=1,np
     if(ok(j))then ! if all criteria met, replace saddle density array value
        if(last_imat_free==0)then ! initialize first element of the linked list
           last_imat_free=last_imat_free+1
           saddle_dens(last_imat_free)=max(saddle_dens(last_imat_free),av_dens(j))
           icurrent(last_imat_free)=clump_nr(j)+(neigh_cl(j)-1)*npeaks_tot
           imat_next(last_imat_free)=0
           imat_prev(last_imat_free)=0
        else
           flaggo=0
           index_next=1e7
           index_prev=1
           imatp=1
           imatn=1
           do imat=1,last_imat_free ! check if we need to update the value in the linked list
              if(clump_nr(j)+(neigh_cl(j)-1)*npeaks_tot==icurrent(imat))then ! update 
                 saddle_dens(imat)=max(saddle_dens(imat),av_dens(j))
                 flaggo=1
              else ! compute the closest neighbors in the list
                 index=icurrent(imat)
                 if(index>index_prev.and.index<clump_nr(j)+(neigh_cl(j)-1)*npeaks_tot)then
                    imatp=imat
                    index_prev=max(index_prev,index)
                 end if
                 if(index<index_next.and.index>clump_nr(j)+(neigh_cl(j)-1)*npeaks_tot)then
                    imatn=imat
                    index_next=min(index_next,index)
                 end if
              end if
           end do
           if(flaggo == 0)then ! create new element of the list and link it 
              last_imat_free=last_imat_free+1
              saddle_dens(last_imat_free)=max(saddle_dens(last_imat_free),av_dens(j))
              icurrent(last_imat_free)=clump_nr(j)+(neigh_cl(j)-1)*npeaks_tot
              imat_next(last_imat_free)=imatn
              imat_prev(last_imat_free)=imatp
              imat_next(imatp)=last_imat_free
              imat_prev(imatn)=last_imat_free
           end if
        end if
     end if
  end do

end subroutine saddlecheck
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine phi_ref_check(ind_cell,cell_index,clump_nr,ok,np)
  use amr_commons, ONLY:flag2,nvector,dp
  use clfind_commons, ONLY: phi_ref
  use poisson_commons, ONLY: phi
  implicit none

  !small routine to check wether neighbor is connected through new densest saddle

  logical,dimension(1:nvector)::ok
  integer,dimension(1:nvector)::cell_index,clump_nr,ind_cell,neigh_cl
  real(dp),dimension(1:nvector)::av_phi
  integer::np,j

  do j=1,np
     neigh_cl(j)=flag2(cell_index(j))!nuber of clump the neighboring cell is in 
  end do
  do j=1,np
     ok(j)=ok(j).and. clump_nr(j)>0 !check that cell is not in a clump that has been merged to zero
     ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) !neighboring cell is in another clump (can be zero)
     av_phi(j)=(phi(cell_index(j))+phi(ind_cell(j)))*0.5 !average pot of cell and neighbor cell
  end do
  do j=1,np
     if(ok(j))then ! if criteria met, reference potential for clump
        phi_ref(clump_nr(j))=min(av_phi(j),phi_ref(clump_nr(j)))
     end if
  end do

end subroutine phi_ref_check
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index(cell_index,cell_levl,xpart,ilevel,n)
  use amr_commons
  implicit none

  integer::n,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart

  !----------------------------------------------------------------------------
  ! This routine returns the index and level of the cell, (at maximum level
  ! ilevel), in which the input the position specified by xpart lies
  !----------------------------------------------------------------------------

  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if

  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
     zz = xpart(i,3)/boxlen + (nz-1)/2.0
     igrid=igrid0
     do j=1,ilevel 
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine read_clumpfind_params()
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  
  !--------------------------------------------------                           
  ! Namelist definitions                                                        
  !--------------------------------------------------                           

  namelist/clumpfind_params/relevance_threshold,density_threshold,mass_threshold,merge_unbound,clinfo

  ! Read namelist file 
  rewind(1)
  read(1,NML=clumpfind_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &CLUMPFIND_PARAMS in parameter file'
  call clean_stop
102 rewind(1)

end subroutine read_clumpfind_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine surface_int(ind_cell,np,ilevel)
  use amr_commons
  use clfind_commons, ONLY: icellp,center_of_mass_tot,Psurf,peak_pos_tot
  use hydro_commons, ONLY: uold,gamma
  implicit none
  integer::np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_cell

  !------------------------------------------------------------
  ! This routine constructs all neighboring leaf cells that 
  ! have a common cell surface at levels 
  ! ilevel-1, ilevel, ilevel+1.
  ! Depending on the action case value, fuctions performing
  ! further checks for the neighbor cells are called.
  !------------------------------------------------------------

  integer::j,ind,nx_loc,i2,j2,k2,ix,iy,iz,idim,jdim,i3,j3,k3

  real(dp)::dx,dx_loc,scale,vol_loc
  !nvector length arrays
  integer ,dimension(1:nvector)::cell_index,cell_levl,clump_nr,indv,neigh_cl
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xtest,r
  real(dp),dimension(1:nvector)::ekk_cell,ekk_neigh,P_cell,P_neigh,r_idim,r_dot_n
  real(dp),dimension(1:3)::skip_loc,n
  logical ,dimension(1:nvector)::ok

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

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ekk_cell=0.; P_neigh=0; P_cell=0
  ! some preliminary action...
  do j=1,np
     indv(j)=(ind_cell(j)-ncoarse-1)/ngridmax+1 ! cell position in grid
     ind_grid(j)=ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
     clump_nr(j)=flag2(ind_cell(j)) ! save clump number
     do jdim=1,ndim
        ekk_cell(j)=ekk_cell(j)+0.5*uold(ind_cell(j),jdim+1)**2
     end do
     ekk_cell(j)=ekk_cell(j)/uold(ind_cell(j),1)
     P_cell(j)=(gamma-1.0)*(uold(ind_cell(j),ndim+2)-ekk_cell(j))
  end do




  
  
  !================================
  ! generate neighbors at level ilevel
  !================================
  ! Generate 3x3 neighboring cells at level ilevel
  do k2=0,2
     do j2=0,2
        do i2=0,2
           if((k2-1.)**2+(j2-1.)**2+(i2-1.)**2==1)then !check whether common face exists 
              
              n=0.
              if (k2==0)n(3)=-1.
              if (k2==2)n(3)=1.
              if (j2==0)n(2)=-1.
              if (j2==2)n(2)=1.
              if (i2==0)n(1)=-1.
              if (i2==2)n(1)=1.
              if (n(1)**2+n(2)**2+n(3)**2/=1)print*,'n has wrong lenght'
              
              
              r=0.
              do j=1,np                 
                 xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc
                 xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc
                 xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc

                 if (clump_nr(j)>0)then                    
                    r(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc*0.5&
                         -center_of_mass_tot(clump_nr(j),1)
                    r(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc*0.5&
                         -center_of_mass_tot(clump_nr(j),2)
                    r(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc*0.5&
                         -center_of_mass_tot(clump_nr(j),3)
                 endif                 
              end do
              
              call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
              do j=1,np           
                 ok(j)=(son(cell_index(j))==0)
              end do
              do j=1,np
                 neigh_cl(j)=flag2(cell_index(j))!nuber of clump the neighboring cell is in 
                 ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) !neighboring cell is in another clump
                 ok(j)=ok(j).and. 0/=clump_nr(j) !clump number is not zero
              end do
              
              r_dot_n=0.
              do j=1,np
                 do idim=1,3
                    r_dot_n(j)=r_dot_n(j)+n(idim)*r(j,idim)
                 end do
              end do
              
              ekk_neigh=0.
              do j=1,np
                 if (ok(j))then 
                    do jdim=1,ndim
                       ekk_neigh(j)=ekk_neigh(j)+0.5*uold(cell_index(j),jdim+1)**2
                    end do
                    ekk_neigh(j)=ekk_neigh(j)/uold(cell_index(j),1)
                    P_neigh(j)=(gamma-1.0)*(uold(cell_index(j),ndim+2)-ekk_neigh(j))
                    Psurf(clump_nr(j))=Psurf(clump_nr(j))+r_dot_n(j)*dx_loc**2*0.5*(P_neigh(j)+P_cell(j))
                 endif
              end do
           endif
        end do
     end do
  end do
  

  !===================================
  ! generate neighbors at level ilevel+1
  !====================================  
  if(ilevel<nlevelmax)then  
     ! Generate 4x4x4 neighboring cells at level ilevel+1 
     do k3=0,3
        do j3=0,3
           do i3=0,3
              if((k3-1.5)**2+(j3-1.5)**2+(i3-1.5)**2==2.75)then !check whether common face exists

                 n=0.
                 if (k3==0)n(3)=-1. 
                 if (k3==3)n(3)=1.
                 if (j3==0)n(2)=-1. 
                 if (j3==3)n(2)=1.
                 if (i3==0)n(1)=-1. 
                 if (i3==3)n(1)=1.
                 if (n(1)**2+n(2)**2+n(3)**2/=1)print*,'n has wrong lenght'

                 r=0.
                 do j=1,np 

                    xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0
                    xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0
                    xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0
                    
                    if (clump_nr(j)>0)then                       
                       r(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass_tot(clump_nr(j),1)
                       r(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass_tot(clump_nr(j),2)
                       r(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass_tot(clump_nr(j),3)
                    endif
                 end do
                 call get_cell_index(cell_index,cell_levl,xtest,ilevel+1,np)

                 ok=.false.
                 do j=1,np
                    !check wether neighbor is in a leaf cell at the right level
                    if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel+1))ok(j)=.true.
                 end do                 

                 do j=1,np
                    neigh_cl(j)=flag2(cell_index(j))!nuber of clump the neighboring cell is in 
                    ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) !neighboring cell is in another clump
                    ok(j)=ok(j).and. 0/=clump_nr(j) !clump number is not zero 
                 end do
                 
                 r_dot_n=0.
                 do j=1,np
                    do idim=1,3
                       r_dot_n(j)=r_dot_n(j)+n(idim)*r(j,idim)
                    end do
                 end do

                 do j=1,np
                    if (ok(j))then
                       do jdim=1,ndim
                          ekk_neigh(j)=ekk_neigh(j)+0.5*uold(cell_index(j),jdim+1)**2
                       end do
                       ekk_neigh(j)=ekk_neigh(j)/uold(cell_index(j),1)
                       P_neigh(j)=(gamma-1.0)*(uold(cell_index(j),ndim+2)-ekk_neigh(j))
                       Psurf(clump_nr(j))=Psurf(clump_nr(j))+r_dot_n(j)*0.25*dx_loc**2*0.5*(P_neigh(j)+P_cell(j))
                       if(debug.and.((P_neigh(j)-P_cell(j))/P_cell(j))**2>4.)print*,'caution, very high p contrast',(((P_neigh(j)-P_cell(j))/P_cell(j))**2)**0.5
                    endif
                 end do                 

              endif
           end do
        end do
     end do
  endif
     

end subroutine surface_int
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
