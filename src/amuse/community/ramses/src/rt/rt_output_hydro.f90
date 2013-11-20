!RT patch: RT variables are output as they are and not divided by gas  
!          density. Also there are checks on zero division to avoid 
!          floating point exceptions.
!          Also added call to output_rtInfo.
!************************************************************************
SUBROUTINE rt_backup_hydro(filename)

!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  use rt_parameters
  implicit none
  character(LEN=80)::filename,filedir,rt_filename

  integer::i,ivar,idim,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
!------------------------------------------------------------------------
  if(verbose)write(*,*)'Entering backup_rt'

  ilun=ncpu+myid+10
     
  if(myid==1)then
     call title(ifout-1,nchar)
     filedir='output_'//TRIM(nchar)//'/'
     rt_filename=TRIM(filedir)//'info_rt_'//TRIM(nchar)//'.txt'
     call output_rtInfo(rt_filename)
  endif                                                           
  
  if(.not.rt)return

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nrtvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do ivar=1,nGroups
                 ! Store photon density in flux units
                 do i=1,ncache
                    xdp(i)=rt_c*rtuold(ind_grid(i)+iskip,iGroups(ivar))
                 end do
                 write(ilun)xdp
                 do idim=1,ndim
                    ! Store photon flux
                    do i=1,ncache
                       xdp(i)=rtuold(ind_grid(i)+iskip,iGroups(ivar)+idim)
                    end do
                    write(ilun)xdp
                 enddo
              end do
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
     
end subroutine rt_backup_hydro

!************************************************************************
SUBROUTINE output_rtInfo(filename)

! Output rt information into info_rt_XXXXX.txt
!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  use rt_parameters
  use rt_cooling_module
  implicit none
  character(LEN=80)::filename
  integer::ilun
  real(dp)::scale_np,scale_pf
  character(LEN=80)::fileloc
!------------------------------------------------------------------------
  if(verbose)write(*,*)'Entering output_rtInfo'

  ilun=myid+10

  ! Conversion factor from user units to cgs units
  call rt_units(scale_np, scale_pf)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  
  ! Write run parameters
  write(ilun,'("nRTvar      =",I11)')nRTvar
  write(ilun,'("nIons       =",I11)')nIons
  write(ilun,'("nGroups     =",I11)')nGroups
  write(ilun,'("iIons       =",I11)')iIons
  write(ilun,*)

  ! Write cooling parameters
  write(ilun,'("X_fraction  =",E23.15)')X
  write(ilun,'("Y_fraction  =",E23.15)')Y
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("unit_np     =",E23.15)')scale_np
  write(ilun,'("unit_pf     =",E23.15)')scale_pf
  write(ilun,'("rt_c_frac   =",E23.15)')rt_c_fraction
  write(ilun,*)

  ! Write polytropic parameters
  write(ilun,'("n_star      =",E23.15)')n_star
  write(ilun,'("T2_star     =",E23.15)')T2_star
  write(ilun,'("g_star      =",E23.15)')g_star
  write(ilun,*)
  call write_group_props(.false.,ilun)

  close(ilun)

end subroutine output_rtInfo

!************************************************************************
SUBROUTINE write_group_props(update,lun)

! Write photon group properties to file or std output.
! lun => File identifier (use 6 for std. output)
!------------------------------------------------------------------------
  use rt_parameters
  use amr_commons,only:myid
  implicit none
  logical::update
  integer::ip,lun
!------------------------------------------------------------------------
  if(myid .ne. 1) RETURN
  if(.not. update) then
     write(lun,*) 'Photon group properties------------------------------ '
  else
     write(lun,*) 'Photon properties have been changed to----------------- '
  endif
  write(lun,901) groupL0(:)
  write(lun,902) groupL1(:)
  write(lun,903) spec2group(:)
  do ip=1,nGroups
     write(lun,907) ip
     write(lun,904) group_egy(ip)
     write(lun,905) group_csn(ip,:)
     write(lun,906) group_cse(ip,:)
  enddo
  write (lun,*) '-------------------------------------------------------'

901 format ('  groupL0  [eV] =', 20f12.3)
902 format ('  groupL1  [eV] =', 20f12.3)
903 format ('  spec2group    =', 20I12)
904 format ('  egy      [eV] =', 20f12.3)
905 format ('  csn    [cm^2] =', 20(1pe12.3))
906 format ('  cse    [cm^2] =', 20(1pe12.3))
907 format ('  ---Group',I2)

END SUBROUTINE write_group_props

!*************************************************************************
SUBROUTINE output_rt_stats

! Output and reset rt statistics. These are cooling statistics and
! star rt feedback statistics
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  implicit none
  integer*8:: max_all, tot_all, cells_all,loopCodes_tot
  integer*8:: loopCodes_all(4)
  integer::info
  real(dp)::step_nPhot_all, step_nStar_all, step_mStar_all
  real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nh, scale_T2
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! Cooling statistics:
  if(rt_output_coolstats) then
     cells_all=0 ; tot_all=0 ; max_all=0 ; loopCodes_all=0
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(n_cool_cells,         cells_all,     1, &
          MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(tot_cool_loopcnt,     tot_all,       1, &
          MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(max_cool_loopcnt,     max_all,       1, &
          MPI_INTEGER,          MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(loopCodes,            loopCodes_all, 4, &
          MPI_INTEGER,          MPI_MAX, MPI_COMM_WORLD, info)
     n_cool_cells     = cells_all ; tot_cool_loopcnt = tot_all
     max_cool_loopcnt = max_all   ; loopCodes        = loopCodes_all
#endif
     if(myid .eq. 1) then
        if(n_cool_cells .eq. 0) n_cool_cells=1.
        write(*, 111) dble(tot_cool_loopcnt)/n_cool_cells,max_cool_loopcnt,rt_advect
        loopCodes_tot = SUM(loopCodes)
        if(loopCodes_tot .gt. 0) then
           write(*, 112) dble(loopCodes)/dble(loopCodes_tot)
        else
           write(*, 112) dble(loopCodes)
        endif
     endif
     max_cool_loopcnt=0; tot_cool_loopcnt=0; n_cool_cells=0; loopCodes(:)=0
  endif ! output_coolstats
111 format(' Coolstats: Avg. # loops = ', f21.6, ', max. # loops = ', I10, ', rt_adv=',L)
112 format(' Subcycling codes [Np, T, xH, xHe]% = ', 4(f7.3, ''))

  ! Stellar rt feedback statistics:
  if(showSEDstats .and. rt_star) then
     step_nPhot_all=0.d0 ; step_nStar_all=0.d0 ; ; step_mStar_all=0.d0
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(step_nPhot,           step_nPhot_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(step_nStar,           step_nStar_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(step_mStar,           step_mStar_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     step_nPhot  = step_nPhot_all
     step_nStar  = step_nStar_all
     step_mStar  = step_mStar_all
#endif
     tot_nPhot = tot_nPhot + step_nPhot
     if(myid .eq. 1)                                                     &
          write(*, 113) step_nPhot, tot_nPhot, step_nStar/dtnew(levelmin)&
          ,step_mStar/dtnew(levelmin), dtnew(levelmin)*scale_t/(3.15569d7)
     step_nPhot = 0.d0 ; step_nStar = 0.d0 ; step_mStar = 0.d0
  endif
113 format(' SED feedback(phot/step/1d50, phot/tot/1d50, *, */Msun , dt[yr])= '  &
                                                             ,10(1pe9.2))
END SUBROUTINE output_rt_stats






