subroutine flag_formation_sites
  use amr_commons
  use pm_commons
  use clfind_commons
  use hydro_commons, only:uold
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !=============================================================================
  ! This routine flags (flag 2 = 1 )  the cells where a sink is going to be formed
  !=============================================================================

  real(dp),dimension(1:nvector,1:3)::pos
  integer,dimension(1:nvector)::cell_index,cell_levl,cc
  integer::j,jj,i,nx_loc
  integer::flag_form,flag_form_tot,info
  logical::ok
  real(dp)::dx,dx_min,dist,scale,tff,acc_r
  real(dp)::fourpi,threepi2
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),dimension(1:npeaks_tot)::peakd
  integer,dimension(1:npeaks_tot)::ind_sort

  !gridspacing and physical scales
  dx=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=dx*scale
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)



  ! loop over sinks and mark all clumps which are already occupied by a sink
  allocate(occupied(1:npeaks_tot),occupied_all(1:npeaks_tot))
  occupied=0; occupied_all=0;
  pos=0.0
  if(myid==1 .and. clinfo)write(*,*)'looping over ',nsink,' sinks and marking their clumps'

  if (smbh)then 
     !block clumps that contain a sink for formation
     do j=1,nsink
        pos(1,1:3)=xsink(j,1:3)
        call cmp_cpumap(pos,cc,1)
        if (cc(1) .eq. myid)then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           if (flag2(cell_index(1))>0)then
              occupied(flag2(cell_index(1)))=1
              if(clinfo)write(*,*)'CPU # ',myid,'blocked clump # ',flag2(cell_index(1)),' for sink production because of sink # ',idsink(j)
           end if
        end if
     end do
  else
     !block peaks that are closer than 4 cells from existing sinks
     do j=1,nsink
        do i=1,npeaks_tot
           dist=(xsink(j,1)-peak_pos_tot(i,1))**2+&
                (xsink(j,2)-peak_pos_tot(i,2))**2+&
                (xsink(j,3)-peak_pos_tot(i,3))**2
           if (dist<(ir_cloud*dx_min)**2)then
              occupied(i)=1
              if(myid==1 .and. clinfo)write(*,*)'blocked clump # ',i,' for sink production because of sink # ',idsink(j)
           end if
        end do
     end do
  end if


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(occupied,occupied_all,npeaks_tot,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  occupied_all=occupied
#endif


  !------------------------------------------------------------------------------
  ! determine whether a peak patch is allowed to form a new sink.
  ! if a new sink has to be created, flag2 is set to the clump number at the peak position
  ! -> criteria to be chosen depend on the physics
  ! -> this routine can be patched
  !------------------------------------------------------------------------------     
  pos=0.0
  flag2=0
  allocate(form(1:npeaks_tot),form_all(1:npeaks_tot))
  form=0; form_all=0;
  flag_form=0
  !sort clumps by peak density in ascending order
  do i=1,npeaks_tot
     peakd(i)=max_dens_tot(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks_tot)
  do j=npeaks_tot,1,-1
     jj=ind_sort(j)
     if (smbh)then
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

        if (ok .eqv. .true.)then
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
     else
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
     end if
  end do

  !for the smbh case, create some output for the new sinks

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
     if(flag_form_tot>0)write(*,'(135A)')'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [cm] size_y [cm] size_z [cm] |v|_CM [u.u.] rho- [H/cc] rho+ [H/cc] rho_av [H/cc] M_cl [M_sol] V_cl [AU^3] rel.  peak_check   ball4_c\heck   isodens_check   clump_check '
     do j=npeaks_tot,1,-1
        jj=ind_sort(j)
        if(form_all(jj) == 1)write(*,'(I6,X,I10,16(1X,1PE14.7))')jj&
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
!             ,ball4_check(jj)&
             ,isodens_check(jj)&
             ,clump_check(jj)
     end do
  end if
  deallocate(occupied,occupied_all)
  deallocate(form,form_all)



end subroutine flag_formation_sites
