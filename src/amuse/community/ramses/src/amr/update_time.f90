subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  integer::ilevel

  real(dp)::dt,econs,mcons
  real(kind=8)::ttend
  real(kind=8),save::ttstart=0
  integer::i,itest,info

  ! Local constants
  dt=dtnew(ilevel)
  itest=0

#ifndef WITHOUTMPI
  if(myid==1)then
     if(ttstart.eq.0.0)ttstart=MPI_WTIME(info)
  endif
#endif

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(nstep_coarse .ne. nstep_coarse_old)then

     !--------------------------
     ! Check mass conservation
     !--------------------------
     if(mass_tot_0==0.0D0)then
        mass_tot_0=mass_tot
        mcons=0.0D0
     else
        mcons=(mass_tot-mass_tot_0)/mass_tot_0
     end if

     !----------------------------
     ! Check energy conservation
     !----------------------------
     if(epot_tot_old.ne.0)then
        epot_tot_int=epot_tot_int + &
             & 0.5D0*(epot_tot_old+epot_tot)*log(aexp/aexp_old)
     end if
     epot_tot_old=epot_tot
     aexp_old=aexp
     if(const==0.0D0)then
        const=epot_tot+ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(ekin_tot+epot_tot-epot_tot_int-const) / &
             &(-(epot_tot-epot_tot_int-const)+ekin_tot)
     end if

     if(mod(nstep_coarse,ncontrol)==0.or.output_done)then
        if(myid==1)then
           
           !-------------------------------
           ! Output AMR structure to screen
           !-------------------------------
           write(*,*)'Mesh structure'
           do i=1,nlevelmax
              if(numbtot(1,i)>0)write(*,999)i,numbtot(1:4,i)
           end do
           
           !----------------------------------------------
           ! Output mass and energy conservation to screen
           !----------------------------------------------
           if(scheme.eq.'induction')then
#ifdef SOLVERmhd
              write(*,778)nstep_coarse,econs,epot_tot,ekin_tot,emag_tot
#endif
           else if(cooling.or.pressure_fix)then
              write(*,778)nstep_coarse,econs,epot_tot,ekin_tot,eint_tot
#ifdef SOLVERmhd
              write(*,'(" emag=",ES9.2)') emag_tot
#endif
           else
              write(*,777)nstep_coarse,mcons,econs,epot_tot,ekin_tot
#ifdef SOLVERmhd
              write(*,'(" emag=",ES9.2)') emag_tot
#endif
           end if
           if(pic)then
              write(*,888)nstep,t,dt,aexp,&
                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
                   & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
           else
              write(*,888)nstep,t,dt,aexp,&
                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
           endif
           itest=1
        end if
        output_done=.false.
     end if

     !---------------
     ! Exit program
     !---------------
     if(t>=tout(noutput).or.aexp>=aout(noutput).or. &
          & nstep_coarse>=nstepmax)then
        if(myid==1)then
           write(*,*)'Run completed'
#ifndef WITHOUTMPI
           ttend=MPI_WTIME(info)
           write(*,*)'Total elapsed time:',ttend-ttstart
#endif
        endif
        call clean_stop
     end if

  end if
  nstep_coarse_old=nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(nstep,ncontrol)==0)then
     if(myid==1.and.itest==0)then
        if(pic)then
           write(*,888)nstep,t,dt,aexp,&
                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
                & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
        else
           write(*,888)nstep,t,dt,aexp,&
                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
        endif
     end if
  end if

  !------------------------
  ! Update time variables
  !------------------------
  t=t+dt
  nstep=nstep+1
  if(cosmo)then
     ! Find neighboring times
     i=1
     do while(tau_frw(i)>t.and.i<n_frw)
        i=i+1
     end do
     ! Interpolate expansion factor
     aexp = aexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & aexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     hexp = hexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & hexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     texp =    t_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          &    t_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
  else
     aexp = 1.0
     hexp = 0.0
     texp = t
  end if

777 format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2)
778 format(' Main step=',i6,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
888 format(' Fine step=',i6,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')
 
end subroutine update_time
  
subroutine clean_stop
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::info
#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif
  stop
end subroutine clean_stop

subroutine writemem(usedmem)
  real::usedmem
  integer::getpagesize

#ifdef NOSYSTEM
!  call PXFSYSCONF(_SC_PAGESIZE,ipagesize,ierror)
  ipagesize=4096
#else
!  ipagesize = getpagesize()
  ipagesize=4096
#endif
  usedmem=dble(usedmem)*dble(ipagesize)

  if(usedmem>1024.**3.)then
     write(*,999)usedmem/1024.**3.
  else if (usedmem>1024.**2.) then
     write(*,998)usedmem/1024.**2
  else if (usedmem>1024.) then
     write(*,997)usedmem/1024.
  endif

997 format(' Used memory:',F6.1,' kb')
998 format(' Used memory:',F6.1,' Mb')
999 format(' Used memory:',F6.1,' Gb')

end subroutine writemem

subroutine getmem(outmem)
  real::outmem
  character(len=300) :: dir, dir2,  cmd, file
  integer::read_status
  file='/proc/self/stat'
  open(unit=1,file=file,form='formatted')
  read(1,'(A300)',IOSTAT=read_status)dir
  close(1)
  if (read_status < 0)then
     outmem=dble(0.)
     if (myid==1)write(*,*)'Problem in checking free memory'
  else
     ind=300
     j=0
     do while (j<23)
        ind=index(dir,' ')
        dir2=dir(ind+1:300)
        j=j+1
        dir=dir2
     end do
     ind=index(dir,' ')
     dir2=dir(1:ind)
     read(dir2,'(I12)')nmem
     outmem=dble(nmem)
  end if

end subroutine getmem

subroutine cmpmem(outmem)
  use amr_commons
  use hydro_commons
  implicit none

  real::outmem,outmem_int,outmem_dp,outmem_qdp
  outmem_int=0.0
  outmem_dp=0.0
  outmem_qdp=0.0

  outmem_dp =outmem_dp +ngridmax*ndim      ! xg
  outmem_int=outmem_int+ngridmax*twondim   ! nbor
  outmem_int=outmem_int+ngridmax           ! father
  outmem_int=outmem_int+ngridmax           ! next
  outmem_int=outmem_int+ngridmax           ! prev
  outmem_int=outmem_int+ngridmax*twotondim ! son 
  outmem_int=outmem_int+ngridmax*twotondim ! flag1
  outmem_int=outmem_int+ngridmax*twotondim ! flag2
  outmem_int=outmem_int+ngridmax*twotondim ! cpu_map1
  outmem_int=outmem_int+ngridmax*twotondim ! cpu_map2
  outmem_qdp=outmem_qdp+ngridmax*twotondim ! hilbert_key

  ! Add communicator variable here

  if(hydro)then
     
  outmem_dp =outmem_dp +ngridmax*twotondim*nvar ! uold
  outmem_dp =outmem_dp +ngridmax*twotondim*nvar ! unew

  if(pressure_fix)then

  outmem_dp =outmem_dp +ngridmax*twotondim ! uold
  outmem_dp =outmem_dp +ngridmax*twotondim ! uold

  endif

  endif

  write(*,*)'Estimated memory=',(outmem_dp*8.+outmem_int*4.+outmem_qdp*8.)/1024./1024.


end subroutine cmpmem
!------------------------------------------------------------------------
SUBROUTINE getProperTime(tau,tproper)
! Calculate proper time tproper corresponding to conformal time tau (both
! in code units).
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::tau, tproper
  integer::i
  if(.not. cosmo .or. tau .eq. 0.d0) then ! this might happen quite often
     tproper = tau
     return
  endif
  i = 1
  do while( tau_frw(i) > tau .and. i < n_frw )
     i = i+1
  end do
  tproper = t_frw(i  )*(tau-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & t_frw(i-1)*(tau-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
END SUBROUTINE getProperTime
!------------------------------------------------------------------------
SUBROUTINE getAgeGyr(t_birth_proper, age)
! Calculate proper time passed, in Gyrs, since proper time t_birth_proper 
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  implicit none
  real(dp):: t_birth_proper, age
  real(dp), parameter:: yr = 3.15569d+07
  real(dp),save:: scale_t_Gyr
  logical::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then 
     ! The timescale has not been initialized
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_t_Gyr = (scale_t/aexp**2)/yr/1.e9
     scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_Gyr
END SUBROUTINE getAgeGyr
!------------------------------------------------------------------------
SUBROUTINE getAgeSec(t_birth_proper, age)
! Calculate proper time passed, in sec, since proper time t_birth_proper 
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  implicit none
  real(dp):: t_birth_proper, age
  real(dp),save:: scale_t_sec
  logical::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then 
     ! The timescale has not been initialized
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_t_sec = (scale_t/aexp**2)
     scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_sec
END SUBROUTINE getAgeSec
!------------------------------------------------------------------------






