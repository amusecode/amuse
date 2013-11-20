subroutine backup_part(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ilun,ipart
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer ,allocatable,dimension(:)::ii
  integer ,allocatable,dimension(:)::ll

  if(verbose)write(*,*)'Entering backup_part'
  
  ilun=2*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  rewind(ilun)
  ! Write header
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)npart
  write(ilun)localseed
  write(ilun)nstar_tot   
  write(ilun)mstar_tot   
  write(ilun)mstar_lost
  write(ilun)nsink
  ! Write position
  allocate(xdp(1:npart))
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=xp(i,idim)
         end if
     end do
     write(ilun)xdp
  end do
  ! Write velocity
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=vp(i,idim)
        end if
     end do
     write(ilun)xdp
  end do
  ! Write mass
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        xdp(ipart)=mp(i)
     end if
  end do
  write(ilun)xdp
  deallocate(xdp)
  ! Write identity
  allocate(ii(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ii(ipart)=idp(i)
     end if
  end do
  write(ilun)ii
  deallocate(ii)
  ! Write level
  allocate(ll(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ll(ipart)=levelp(i)
     end if
  end do
  write(ilun)ll
  deallocate(ll)

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Write potential (added by AP)
  allocate(xdp(1:npart))
  ipart=0
  do i=1, npartmax
     if(levelp(i)>0) then
        ipart=ipart+1
        xdp(ipart)=ptcl_phi(i)
     end if
  end do
  write(ilun)xdp
  deallocate(xdp)
#endif

  ! Write birth epoch
  if(star.or.sink)then
     allocate(xdp(1:npart))
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=tp(i)
        end if
     end do
     write(ilun)xdp
     ! Write metallicity
     if(metal)then
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=zp(i)
           end if
        end do
        write(ilun)xdp
     end if
     deallocate(xdp)
  end if
  !======================
  ! Write sink properties
  !======================
  if(sink)then
     if(nsink>0)then
        allocate(xdp(1:nsink))
        do i=1,nsink
           xdp(i)=msink(i)
        end do
        write(ilun)xdp ! Write sink mass
        do i=1,nsink
           xdp(i)=tsink(i)
        end do
        write(ilun)xdp ! Write sink birth epoch
        do idim=1,ndim
           do i=1,nsink
              xdp(i)=xsink(i,idim)
           end do
           write(ilun)xdp ! Write sink position
        enddo
        do idim=1,ndim
           do i=1,nsink
              xdp(i)=vsink(i,idim)
           end do
           write(ilun)xdp ! Write sink velocity
        enddo
        do idim=1,ndim
           do i=1,nsink
              xdp(i)=lsink(i,idim)
           end do
           write(ilun)xdp ! Write sink angular momentum
        enddo
        do i=1,nsink
           xdp(i)=delta_mass(i)
        end do
        write(ilun)xdp ! Write sink accumulated rest mass energy
        do i=1,nsink
           xdp(i)=acc_rate(i)
        end do
        write(ilun)xdp ! Write sink accretion rate
        deallocate(xdp)
        allocate(ii(1:nsink))
        do i=1,nsink
           ii(i)=idsink(i)
        end do
        write(ilun)ii ! Write sink index
        do i=1,nsink
           ii(i)=level_sink(i)
        end do
        write(ilun)ii ! Write sink level
        deallocate(ii)
        write(ilun)ncloud_sink ! Write ncloud
     endif
  endif

  close(ilun)

end subroutine backup_part

subroutine output_sink(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ipart,isink
  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
  real(dp)::scale,l_abs,rot_period,dx_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_sink'

  ilun=myid+10

  ! Conversion factor from user units to cgs units                                                                   
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  
  if(verbose)write(*,*)'Entering output_sink'
  
  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace')
  !======================
  ! Write sink properties
  !======================
  write(ilun,*)'Number of sink = ',nsink

  write(ilun,'(" =============================================================================================================================== ")')
  write(ilun,'(" Id     Mass(Msol)     x           y           z           vx        vy        vz     new  rot_period[y] lx/|l|  ly/|l|  lz/|l| ")')
  write(ilun,'(" =============================================================================================================================== ")')
  
  do isink=1,nsink
     l_abs=max((lsink(isink,1)**2+lsink(isink,2)**2+lsink(isink,3)**2)**0.5,1.d-50)
     rot_period=32*3.1415*msink(isink)*(dx_min)**2/(5*l_abs+tiny(0.d0))
     write(ilun,'(I6,2X,F8.4,3(2X,F10.7),3(2X,F6.3),4X,I1,2X,F13.5,3(2X,F6.3))')idsink(isink),msink(isink)*scale_m/2d33,xsink(isink,1:ndim), &
          vsink(isink,1:ndim),new_born_all(isink),rot_period*scale_t/(3600*24*365),lsink(isink,1)/l_abs,lsink(isink,2)/l_abs,lsink(isink,3)/l_abs
  end do
  write(ilun,'(" =============================================================================================================================== ")') 
  close(ilun)

end subroutine output_sink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_sink_csv(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ipart,isink
  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
  real(dp)::scale
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_sink_csv'

  ilun=myid+10

  ! Conversion factor from user units to cgs units 
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  if(verbose)write(*,*)'Entering output_sink'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace')
  !======================
  ! Write sink properties
  !======================
  do isink=1,nsink
     write(ilun,'(I10,8(A1,ES20.10))')idsink(isink),',',msink(isink)*scale_m/2d33,',',xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),',',vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),',',t-tsink(isink)
  end do

  close(ilun)

end subroutine output_sink_csv
