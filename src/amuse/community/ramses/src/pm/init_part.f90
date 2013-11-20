subroutine init_part
  use amr_commons
  use pm_commons
  use clfind_commons
#ifdef RT
  use rt_parameters,only: convert_birth_times
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from grafic files
  !------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  integer::npart2,ndim2,ncpu2
  integer::ipart,jpart,ipart_old,ilevel,idim
  integer::i,igrid,ncache,ngrid,iskip,isink
  integer::ind,ix,iy,iz,ilun,info,icpu,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,indglob,npart_new
  real(dp)::dx,xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp)::scale,dx_loc,rr,rmax,dx_min
  integer::ncode,bit_length,temp
  real(kind=8)::bscale
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nvector)::ind_grid,ind_cell,cc,ii
  integer ,dimension(1:ncpu)::npart_cpu,npart_all
  real(dp),allocatable,dimension(:)::xdp
  integer ,allocatable,dimension(:)::isp

  real(kind=4),allocatable,dimension(:,:)::init_plane,init_plane_x
  real(dp),allocatable,dimension(:,:,:)::init_array,init_array_x
  real(kind=8),dimension(1:nvector,1:3)::xx,vv,xs
  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer,dimension(1:nvector)::ixx,iyy,izz
  real(qdp),dimension(1:nvector)::order
  real(kind=8),dimension(1:nvector)::mm
  real(kind=8)::dispmax=0.0,dispall
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:3)::centerofmass

  integer::ibuf,tag=101,tagf=102,tagu=102
  integer::countsend,countrecv
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
#endif

  logical::error,keep_part,eof,jumped,ic_sink=.false.,read_pos=.false.,ok
  character(LEN=80)::filename,filename_x
  character(LEN=80)::fileloc
  character(LEN=20)::filetype_loc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_part'

  if(allocated(xp))then
     if(verbose)write(*,*)'Initial conditions already set'
     return
  end if

  ! Allocate particle variables
  allocate(xp    (npartmax,ndim))
  allocate(vp    (npartmax,ndim))
  allocate(mp    (npartmax))
  allocate(nextp (npartmax))
  allocate(prevp (npartmax))
  allocate(levelp(npartmax))
  allocate(idp   (npartmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(ptcl_phi(npartmax))
#endif
  xp=0.0; vp=0.0; mp=0.0; levelp=0; idp=0
  if(star.or.sink)then
     allocate(tp(npartmax))
     tp=0.0
     if(metal)then
        allocate(zp(npartmax))
        zp=0.0
     end if
  end if
  if(sink)then
     ! Important persistent sink variables
     allocate(msink(1:nsinkmax))
     allocate(tsink(1:nsinkmax))
     allocate(idsink(1:nsinkmax))
     idsink=0 ! Important: need to set idsink to zero
     nindsink=MAXVAL(idsink) ! Reset max index
     allocate(xsink(1:nsinkmax,1:ndim))
     allocate(vsink(1:nsinkmax,1:ndim))
     allocate(vsold(1:nsinkmax,1:ndim,levelmin:nlevelmax))
     allocate(vsnew(1:nsinkmax,1:ndim,levelmin:nlevelmax))
     allocate(fsink_partial(1:nsinkmax,1:ndim,levelmin:nlevelmax))
     allocate(fsink(1:nsinkmax,1:ndim))
     allocate(acc_rate(1:nsinkmax))
     acc_rate=0.
     allocate(acc_lum(1:nsinkmax))
     acc_lum=0.
     allocate(lsink(1:nsinkmax,1:3))
     lsink=0.d0
     allocate(level_sink(1:nsinkmax))
     allocate(delta_mass(1:nsinkmax))
     ! Temporary sink variables
     allocate(total_volume(1:nsinkmax))
     allocate(wden(1:nsinkmax))
     allocate(wmom(1:nsinkmax,1:ndim))
     allocate(weth(1:nsinkmax))
     allocate(wvol(1:nsinkmax))
     allocate(wden_new(1:nsinkmax))
     allocate(wmom_new(1:nsinkmax,1:ndim))
     allocate(weth_new(1:nsinkmax))
     allocate(wvol_new(1:nsinkmax))
     allocate(divsink(1:nsinkmax,1:nlevelmax))
     allocate(divsink_new(1:nsinkmax,1:nlevelmax))
     divsink=0.;  divsink_new=0.
     allocate(msink_new(1:nsinkmax))
     allocate(msink_all(1:nsinkmax))
     allocate(tsink_new(1:nsinkmax))
     allocate(tsink_all(1:nsinkmax))
     allocate(idsink_new(1:nsinkmax))
     allocate(idsink_all(1:nsinkmax))
     allocate(idsink_old(1:nsinkmax))
     allocate(vsink_new(1:nsinkmax,1:ndim))
     allocate(vsink_all(1:nsinkmax,1:ndim))
     allocate(fsink_new(1:nsinkmax,1:ndim))
     allocate(fsink_all(1:nsinkmax,1:ndim))
     allocate(lsink_new(1:nsinkmax,1:3))
     allocate(lsink_all(1:nsinkmax,1:3))
     allocate(xsink_new(1:nsinkmax,1:ndim))
     allocate(xsink_all(1:nsinkmax,1:ndim))
     allocate(sink_jump(1:nsinkmax,1:ndim,levelmin:nlevelmax))
     sink_jump=0.d0
     allocate(level_sink_all(1:nsinkmax))
     allocate(level_sink_new(1:nsinkmax))
     allocate(dMsink_overdt(1:nsinkmax))
     allocate(r2sink(1:nsinkmax))
     allocate(r2k(1:nsinkmax))
     allocate(v2sink(1:nsinkmax))
     allocate(c2sink(1:nsinkmax))
     allocate(v2sink_new(1:nsinkmax))
     allocate(c2sink_new(1:nsinkmax))
     allocate(v2sink_all(1:nsinkmax))
     allocate(c2sink_all(1:nsinkmax))
     allocate(weighted_density(1:nsinkmax,1:nlevelmax))
     allocate(weighted_volume(1:nsinkmax,1:nlevelmax))
     allocate(weighted_ethermal(1:nsinkmax,1:nlevelmax))
     allocate(weighted_momentum(1:nsinkmax,1:nlevelmax,1:ndim))
     allocate(oksink_new(1:nsinkmax))
     allocate(oksink_all(1:nsinkmax))
     allocate(idsink_sort(1:nsinkmax))
     allocate(xmsink(1:nsinkmax))
     allocate(delta_mass_new(1:nsinkmax),delta_mass_all(1:nsinkmax))
     allocate(vol_gas_agn(1:nsinkmax),mass_gas_agn(1:nsinkmax))
     allocate(ind_blast_agn(1:nsinkmax),mass_blast_agn(1:nsinkmax),vol_blast_agn(1:nsinkmax))
     allocate(p_agn(1:nsinkmax),vol_gas_agn_all(1:nsinkmax),mass_gas_agn_all(1:nsinkmax))
     allocate(ok_blast_agn(1:nsinkmax),ok_blast_agn_all(1:nsinkmax))
     allocate(new_born(1:nsinkmax),new_born_all(1:nsinkmax))
  endif

  !--------------------
  ! Read part.tmp file
  !--------------------
  if(nrestart>0)then

     ilun=2*ncpu+myid+10
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/part_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     open(unit=ilun,file=fileloc,form='unformatted')
     rewind(ilun)
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)npart2
     read(ilun)localseed
     read(ilun)nstar_tot
     read(ilun)mstar_tot
     read(ilun)mstar_lost
     read(ilun)nsink
     if(ncpu2.ne.ncpu.or.ndim2.ne.ndim.or.npart2.gt.npartmax)then
        write(*,*)'File part.tmp not compatible'
        write(*,*)'Found   =',ncpu2,ndim2,npart2
        write(*,*)'Expected=',ncpu,ndim,npartmax
        call clean_stop
     end if
     ! Read position
     allocate(xdp(1:npart2))
     do idim=1,ndim
        read(ilun)xdp
        xp(1:npart2,idim)=xdp
     end do
     ! Read velocity
     do idim=1,ndim
        read(ilun)xdp
        vp(1:npart2,idim)=xdp
     end do
     ! Read mass
     read(ilun)xdp
     mp(1:npart2)=xdp
     deallocate(xdp)
     ! Read identity
     allocate(isp(1:npart2))
     read(ilun)isp
     idp(1:npart2)=isp
     ! Read level
     read(ilun)isp
     levelp(1:npart2)=isp
     deallocate(isp)
     if(star.or.sink)then
        ! Read birth epoch
        allocate(xdp(1:npart2))
        read(ilun)xdp
        tp(1:npart2)=xdp
#ifdef RT
        if(convert_birth_times) then
           do i = 1, npart2 ! Convert birth time to proper for RT postpr.
              call getProperTime(tp(i),tp(i))
           enddo
        endif
#endif
        if(metal)then
           ! Read metallicity
           read(ilun)xdp
           zp(1:npart2)=xdp
        end if
        deallocate(xdp)
     end if
     ! Read sink properties
     if(sink)then
        if(nsink>0)then
           allocate(xdp(1:nsink))
           read(ilun)xdp ! Read sink mass
           msink(1:nsink)=xdp
           read(ilun)xdp ! Read sink birth epoch
           tsink(1:nsink)=xdp
           do idim=1,ndim
              read(ilun)xdp ! Read sink position
              xsink(1:nsink,idim)=xdp
           end do
           do idim=1,ndim
              read(ilun)xdp ! Read sink velocity
              vsink(1:nsink,idim)=xdp
           end do
           do idim=1,3
              read(ilun)xdp ! Read sink angular momentum
              lsink(1:nsink,idim)=xdp
           end do
           read(ilun)xdp ! Read sink accumulated rest mass energy
           delta_mass(1:nsink)=xdp
           read(ilun)xdp ! Read sink accretion rate
           acc_rate(1:nsink)=xdp
           deallocate(xdp)
           allocate(isp(1:nsink))
           read(ilun)isp ! Read sink index
           idsink(1:nsink)=isp
           read(ilun)isp ! Read sink level
           level_sink(1:nsink)=isp
           nindsink=MAXVAL(idsink) ! Reset max index
           deallocate(isp)
           read(ilun)ncloud_sink
        end if
     endif
     close(ilun)
     if(debug)write(*,*)'part.tmp read for processor ',myid
     npart=npart2

     !allow for ncloud_sink to change at restart
     call compute_ncloud_sink


     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     if(sink .and. ir_feedback)then
        do i=1,nsink
           acc_lum(i)=ir_eff*acc_rate(i)*msink(i)/(5*6.955d10/scale_l)
        end do
     end if
  else     

     filetype_loc=filetype
     if(.not. cosmo)filetype_loc='ascii'

     select case (filetype_loc)

     case ('grafic')

        !----------------------------------------------------
        ! Reading initial conditions GRAFIC2 multigrid arrays  
        !----------------------------------------------------
        ipart=0
        ! Loop over initial condition levels
        do ilevel=levelmin,nlevelmax
           
           if(initfile(ilevel)==' ')cycle
           
           ! Mesh size at level ilevel in coarse cell units
           dx=0.5D0**ilevel
           
           ! Set position of cell centers relative to grid center
           do ind=1,twotondim
              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
              if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
              if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
           end do
           
           !--------------------------------------------------------------
           ! First step: compute level boundaries and particle positions
           !--------------------------------------------------------------
           i1_min=n1(ilevel)+1; i1_max=0
           i2_min=n2(ilevel)+1; i2_max=0
           i3_min=n3(ilevel)+1; i3_max=0
           ipart_old=ipart
           
           ! Loop over grids by vector sweeps
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
                    ind_cell(i)=iskip+ind_grid(i)
                 end do
                 do i=1,ngrid
                    xx1=xg(ind_grid(i),1)+xc(ind,1)
                    xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
                    xx2=xg(ind_grid(i),2)+xc(ind,2)
                    xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
                    xx3=xg(ind_grid(i),3)+xc(ind,3)
                    xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
                    i1_min=MIN(i1_min,int(xx1)+1)
                    i1_max=MAX(i1_max,int(xx1)+1)
                    i2_min=MIN(i2_min,int(xx2)+1)
                    i2_max=MAX(i2_max,int(xx2)+1)
                    i3_min=MIN(i3_min,int(xx3)+1)
                    i3_max=MAX(i3_max,int(xx3)+1)
                    keep_part=son(ind_cell(i))==0
                    if(keep_part)then
                       ipart=ipart+1
                       if(ipart>npartmax)then
                          write(*,*)'Maximum number of particles incorrect'
                          write(*,*)'npartmax should be greater than',ipart
                          call clean_stop
                       endif
                       if(ndim>0)xp(ipart,1)=xg(ind_grid(i),1)+xc(ind,1)
                       if(ndim>1)xp(ipart,2)=xg(ind_grid(i),2)+xc(ind,2)
                       if(ndim>2)xp(ipart,3)=xg(ind_grid(i),3)+xc(ind,3)
                       mp(ipart)=0.5d0**(3*ilevel)*(1.0d0-omega_b/omega_m)
                    end if
                 end do
              end do
              ! End loop over cells
           end do
           ! End loop over grids
           
           ! Check that all grids are within initial condition region
           error=.false.
           if(active(ilevel)%ngrid>0)then
              if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
              if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
              if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
           end if
           if(error) then
              write(*,*)'Some grid are outside initial conditions sub-volume'
              write(*,*)'for ilevel=',ilevel
              write(*,*)i1_min,i1_max
              write(*,*)i2_min,i2_max
              write(*,*)i3_min,i3_max
              write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
              call clean_stop
           end if
           if(debug)then
              write(*,*)myid,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
           endif
           
           !---------------------------------------------------------------------
           ! Second step: read initial condition file and set particle velocities
           !---------------------------------------------------------------------
           ! Allocate initial conditions array
           if(active(ilevel)%ngrid>0)then
              allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
              allocate(init_array_x(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
              init_array=0d0
              init_array_x=0d0
           end if
           allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
           allocate(init_plane_x(1:n1(ilevel),1:n2(ilevel)))
           
           ! Loop over input variables
           do idim=1,ndim
              
              ! Read dark matter initial displacement field
              if(multiple)then
                 call title(myid,nchar)
                 if(idim==1)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
                 if(idim==2)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
                 if(idim==3)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
              else
                 if(idim==1)filename=TRIM(initfile(ilevel))//'/ic_velcx'
                 if(idim==2)filename=TRIM(initfile(ilevel))//'/ic_velcy'
                 if(idim==3)filename=TRIM(initfile(ilevel))//'/ic_velcz'

                 if(idim==1)filename_x=TRIM(initfile(ilevel))//'/ic_poscx'
                 if(idim==2)filename_x=TRIM(initfile(ilevel))//'/ic_poscy'
                 if(idim==3)filename_x=TRIM(initfile(ilevel))//'/ic_poscz'

                 INQUIRE(file=filename_x,exist=ok)
                 if(.not.ok)then
                    read_pos = .false.
                 else
                    read_pos = .true.
                    if(myid==1)write(*,*)'Reading file '//TRIM(filename_x)
                 end if

              endif

              if(myid==1)write(*,*)'Reading file '//TRIM(filename)
                               
              if(multiple)then
                 ilun=myid+10
                 open(ilun,file=filename,form='unformatted')
                 rewind ilun
                 read(ilun) ! skip first line
                 do i3=1,n3(ilevel)
                    read(ilun)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                    if(active(ilevel)%ngrid>0)then
                       if(i3.ge.i3_min.and.i3.le.i3_max)then
                          init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                               & init_plane(i1_min:i1_max,i2_min:i2_max)
                       end if
                    endif
                 end do
                 close(ilun)
              else
                 if(myid==1)then
                    open(10,file=filename,form='unformatted')
                    rewind 10
                    read(10) ! skip first line
                 end if
                 do i3=1,n3(ilevel)
                    if(myid==1)then
                       if(debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                       read(10)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                    else
                       init_plane=0.0
                    endif
                    buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                    call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
                    
                    if(active(ilevel)%ngrid>0)then
                       if(i3.ge.i3_min.and.i3.le.i3_max)then
                          init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                               & init_plane(i1_min:i1_max,i2_min:i2_max)
                       end if
                    endif
                 end do
                 if(myid==1)close(10)

                 if(read_pos) then
                    if(myid==1)then
                       open(10,file=filename_x,form='unformatted')
                       rewind 10
                       read(10) ! skip first line
                    end if
                    do i3=1,n3(ilevel)
                       if(myid==1)then
                          if(debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                          read(10)((init_plane_x(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                       else
                          init_plane_x=0.0
                       endif
                       buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                       call MPI_BCAST(init_plane_x,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
                       if(active(ilevel)%ngrid>0)then
                          if(i3.ge.i3_min.and.i3.le.i3_max)then
                             init_array_x(i1_min:i1_max,i2_min:i2_max,i3) = &
                                  & init_plane_x(i1_min:i1_max,i2_min:i2_max)
                          end if
                       endif
                    end do
                    if(myid==1)close(10)
                 end if

              endif
              
              if(active(ilevel)%ngrid>0)then
                 ! Rescale initial displacement field to code units
                 init_array=dfact(ilevel)*dx/dxini(ilevel)*init_array/vfact(ilevel)
                 if(read_pos)then
                    init_array_x = init_array_x/boxlen_ini
                 endif
                 ! Loop over grids by vector sweeps
                 ipart=ipart_old
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
                          ind_cell(i)=iskip+ind_grid(i)
                       end do
                       do i=1,ngrid
                          xx1=xg(ind_grid(i),1)+xc(ind,1)
                          xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
                          xx2=xg(ind_grid(i),2)+xc(ind,2)
                          xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
                          xx3=xg(ind_grid(i),3)+xc(ind,3)
                          xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
                          i1=int(xx1)+1
                          i1=int(xx1)+1
                          i2=int(xx2)+1
                          i2=int(xx2)+1
                          i3=int(xx3)+1
                          i3=int(xx3)+1
                          keep_part=son(ind_cell(i))==0
                          if(keep_part)then
                             ipart=ipart+1
                             vp(ipart,idim)=init_array(i1,i2,i3)
                             if(.not. read_pos)then
                                dispmax=max(dispmax,abs(init_array(i1,i2,i3)/dx))
                             else
                                xp(ipart,idim)=xg(ind_grid(i),idim)+xc(ind,idim)+init_array_x(i1,i2,i3)
                                dispmax=max(dispmax,abs(init_array_x(i1,i2,i3)/dx))
                             endif
                          end if
                       end do
                    end do
                    ! End loop over cells
                 end do
                 ! End loop over grids
              endif

           end do
           ! End loop over input variables
           
           ! Deallocate initial conditions array
           if(active(ilevel)%ngrid>0)then
              deallocate(init_array,init_array_x)
           end if
           deallocate(init_plane,init_plane_x)
           
           if(debug)write(*,*)'npart=',ipart,'/',npartmax,' for PE=',myid
           
        end do
        ! End loop over levels
        
        ! Initial particle number
        npart=ipart
        
        ! Move particle according to Zeldovich approximation
        if(.not. read_pos)then
           xp(1:npart,1:ndim)=xp(1:npart,1:ndim)+vp(1:npart,1:ndim)
        endif

        ! Scale displacement to velocity
        vp(1:npart,1:ndim)=vfact(1)*vp(1:npart,1:ndim)
        
        ! Periodic box
        do ipart=1,npart
#if NDIM>0
           if(xp(ipart,1)<  0.0d0  )xp(ipart,1)=xp(ipart,1)+dble(nx)
           if(xp(ipart,1)>=dble(nx))xp(ipart,1)=xp(ipart,1)-dble(nx)
#endif
#if NDIM>1
           if(xp(ipart,2)<  0.0d0  )xp(ipart,2)=xp(ipart,2)+dble(ny)
           if(xp(ipart,2)>=dble(ny))xp(ipart,2)=xp(ipart,2)-dble(ny)
#endif
#if NDIM>2
           if(xp(ipart,3)<  0.0d0  )xp(ipart,3)=xp(ipart,3)+dble(nz)
           if(xp(ipart,3)>=dble(nz))xp(ipart,3)=xp(ipart,3)-dble(nz)
#endif
        end do
        
#ifndef WITHOUTMPI        
        ! Compute particle Hilbert ordering
        sendbuf=0
        do ipart=1,npart
           xx(1,1:3)=xp(ipart,1:3)
           xx_dp(1,1:3)=xx(1,1:3)
           call cmp_cpumap(xx_dp,cc,1)
           if(cc(1).ne.myid)sendbuf(cc(1))=sendbuf(cc(1))+1
        end do
           
        ! Allocate communication buffer in emission
        do icpu=1,ncpu
           ncache=sendbuf(icpu)
           if(ncache>0)then
              allocate(emission(icpu,1)%up(1:ncache,1:twondim+1))
           end if
        end do

        ! Fill communicators
        jpart=0
        sendbuf=0
        do ipart=1,npart
           xx(1,1:3)=xp(ipart,1:3)
           xx_dp(1,1:3)=xx(1,1:3)
           call cmp_cpumap(xx_dp,cc,1)
           if(cc(1).ne.myid)then
              icpu=cc(1)
              sendbuf(icpu)=sendbuf(icpu)+1
              ibuf=sendbuf(icpu)
              emission(icpu,1)%up(ibuf,1)=xp(ipart,1)
              emission(icpu,1)%up(ibuf,2)=xp(ipart,2)
              emission(icpu,1)%up(ibuf,3)=xp(ipart,3)
              emission(icpu,1)%up(ibuf,4)=vp(ipart,1)
              emission(icpu,1)%up(ibuf,5)=vp(ipart,2)
              emission(icpu,1)%up(ibuf,6)=vp(ipart,3)
              emission(icpu,1)%up(ibuf,7)=mp(ipart)
           else
              jpart=jpart+1
              xp(jpart,1:3)=xp(ipart,1:3)
              vp(jpart,1:3)=vp(ipart,1:3)
              mp(jpart)    =mp(ipart)
           endif
        end do
        
        ! Communicate virtual particle number to parent cpu
        call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

        ! Compute total number of newly created particles
        npart_new=0
        do icpu=1,ncpu
           npart_new=npart_new+recvbuf(icpu)
        end do

        if(jpart+npart_new.gt.npartmax)then
           write(*,*)'No more free memory for particles'
           write(*,*)'Increase npartmax'
           write(*,*)myid
           write(*,*)jpart,npart_new
           write(*,*)bound_key
           call MPI_ABORT(MPI_COMM_WORLD,1,info)
        end if

        ! Allocate communication buffer in reception
        do icpu=1,ncpu
           ncache=recvbuf(icpu)
           if(ncache>0)then
              allocate(reception(icpu,1)%up(1:ncache,1:twondim+1))
           end if
        end do

        ! Receive particles
        countrecv=0
        do icpu=1,ncpu
           ncache=recvbuf(icpu)
           if(ncache>0)then
              buf_count=ncache*(twondim+1)
              countrecv=countrecv+1
              call MPI_IRECV(reception(icpu,1)%up,buf_count, &
                   & MPI_DOUBLE_PRECISION,icpu-1,&
                   & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
           end if
        end do
        
        ! Send particles
        countsend=0
        do icpu=1,ncpu
           ncache=sendbuf(icpu)
           if(ncache>0)then
              buf_count=ncache*(twondim+1)
              countsend=countsend+1
              call MPI_ISEND(emission(icpu,1)%up,buf_count, &
                   & MPI_DOUBLE_PRECISION,icpu-1,&
                   & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
           end if
        end do
        
        ! Wait for full completion of receives
        call MPI_WAITALL(countrecv,reqrecv,statuses,info)
        
        ! Wait for full completion of sends
        call MPI_WAITALL(countsend,reqsend,statuses,info)

        ! Create new particles
        do icpu=1,ncpu
           do ibuf=1,recvbuf(icpu)
              jpart=jpart+1
              xp(jpart,1)=reception(icpu,1)%up(ibuf,1)
              xp(jpart,2)=reception(icpu,1)%up(ibuf,2)
              xp(jpart,3)=reception(icpu,1)%up(ibuf,3)
              vp(jpart,1)=reception(icpu,1)%up(ibuf,4)
              vp(jpart,2)=reception(icpu,1)%up(ibuf,5)
              vp(jpart,3)=reception(icpu,1)%up(ibuf,6)
              mp(jpart)  =reception(icpu,1)%up(ibuf,7)
           end do
        end do
        
        ! Erase old particles
        do ipart=jpart+1,npart
           xp(ipart,1)=0d0
           xp(ipart,2)=0d0
           xp(ipart,3)=0d0
           vp(ipart,1)=0d0
           vp(ipart,2)=0d0
           vp(ipart,3)=0d0
           mp(ipart)  =0d0
        end do
        npart=jpart

        ! Deallocate communicators
        do icpu=1,ncpu
           if(sendbuf(icpu)>0)deallocate(emission(icpu,1)%up)
           if(recvbuf(icpu)>0)deallocate(reception(icpu,1)%up)
        end do

        write(*,*)'npart=',ipart,'/',npartmax,' for PE=',myid
#endif

        ! Compute particle initial level
        do ipart=1,npart
           levelp(ipart)=levelmin
        end do

        ! Compute particle initial age and metallicity
        if(star.or.sink)then
           do ipart=1,npart
              tp(ipart)=0d0
              if(metal)then
                 zp(ipart)=0d0
              end if
           end do
        end if

        ! Compute particle initial identity
        npart_cpu=0; npart_all=0
        npart_cpu(myid)=npart
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        npart_cpu(1)=npart_all(1)
#endif
        do icpu=2,ncpu
           npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
        end do
        if(myid==1)then
           do ipart=1,npart
              idp(ipart)=ipart
           end do
        else
           do ipart=1,npart
              idp(ipart)=npart_cpu(myid-1)+ipart
           end do
        end if

     case ('ascii')

        ! Local particle count
        ipart=0

        if(TRIM(initfile(levelmin)).NE.' ')then

        filename=TRIM(initfile(levelmin))//'/ic_part'
        if(myid==1)then
           open(10,file=filename,form='formatted')
           indglob=0
        end if
        eof=.false.

        do while (.not.eof)
           xx=0.0
           if(myid==1)then
              jpart=0
              do i=1,nvector
                 read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
                 jpart=jpart+1
                 indglob=indglob+1
                 xx(i,1)=xx1+boxlen/2.0
                 xx(i,2)=xx2+boxlen/2.0
                 xx(i,3)=xx3+boxlen/2.0
                 vv(i,1)=vv1
                 vv(i,2)=vv2
                 vv(i,3)=vv3
                 mm(i  )=mm1
                 ii(i  )=indglob
              end do
100           continue
              if(jpart<nvector)eof=.true.
           endif
           buf_count=nvector*3
#ifndef WITHOUTMPI
           call MPI_BCAST(xx,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(vv,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(mm,nvector  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(ii,nvector  ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(eof,1       ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(jpart,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
           call cmp_cpumap(xx,cc,jpart)
#endif

           do i=1,jpart
#ifndef WITHOUTMPI
              if(cc(i)==myid)then
#endif
                 ipart=ipart+1
                 xp(ipart,1:3)=xx(i,1:3)
                 vp(ipart,1:3)=vv(i,1:3)
                 mp(ipart)    =mm(i)
                 levelp(ipart)=levelmin
                 idp(ipart)   =ii(i)
#ifndef WITHOUTMPI
              endif
#endif
           enddo

        end do
        if(myid==1)close(10)

        end if
        npart=ipart

        ! Compute total number of particle
        npart_cpu=0; npart_all=0
        npart_cpu(myid)=npart
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        npart_cpu(1)=npart_all(1)
#endif
        do icpu=2,ncpu
           npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
        end do
        if(debug)write(*,*)'npart=',npart,'/',npart_cpu(ncpu)

     case ('gadget')
        call load_gadget

     case DEFAULT
        write(*,*) 'Unsupported format file ' // filetype
        call clean_stop

     end select

     ! Read initial sink particles
     ! Sink particles that exist at the beginning of the simu are always (independent of the type of
     ! the other ic files) read from a text file.
     if(sink)then
        call compute_ncloud_sink

        nx_loc=(icoarse_max-icoarse_min+1)
        scale=boxlen/dble(nx_loc)
        dx_min=scale*0.5D0**nlevelmax/aexp

        nsink=0
        if(TRIM(initfile(levelmin)).NE.' ')then
           filename=TRIM(initfile(levelmin))//'/ic_sink'
        else
           filename='ic_sink'
        end if
        INQUIRE(FILE=filename, EXIST=ic_sink)
        if (myid==1)write(*,*),'Looking for file ic_sink: ',filename
        if (ic_sink)then
           open(10,file=filename,form='formatted')
           eof=.false.
           if (myid==1)write(*,*)'Reading_file ',filename
           do
              read(10,*,end=102)mm1,xx1,xx2,xx3,vv1,vv2,vv3,ll1,ll2,ll3
              nsink=nsink+1
              idsink(nsink)=nsink
              msink(nsink)=mm1
              xsink(nsink,1)=xx1
              xsink(nsink,2)=xx2
              xsink(nsink,3)=xx3
              vsink(nsink,1)=vv1
              vsink(nsink,2)=vv2
              vsink(nsink,3)=vv3
              lsink(nsink,1)=ll1
              lsink(nsink,2)=ll2
              lsink(nsink,3)=ll3
              tsink(nsink)=0.
              level_sink(nsink)=levelmin
           end do
102        continue
           close(10)
        end if
        nindsink=MAXVAL(idsink) ! Reset max index
        if (myid==1.and.nsink==0)write(*,*)'File ic_sink not found: starting without sink particles!'
        if (myid==1.and.nsink>0.and.verbose)then
           write(*,*),'sinks read from file ic_sink'
           write(*,*),'   id    m       x       y       z       vx      vy      vz      lx      ly      lz  '
           write(*,*),'====================================================================================='
           do isink=1,nsink
              write(*,'(I6,X,F7.3,3(X,F7.3),3(X,F7.3),3(X,F7.3))'),idsink(isink),msink(isink),xsink(isink,1:ndim),&
                   vsink(isink,1:ndim),lsink(isink,1:ndim)
           end do
        end if

        ! Loop over sinks
        do isink=1,nsink
           xs(1,1:ndim)=xsink(isink,1:ndim)
           call cmp_cpumap(xs,cc,1)

           ! Create central cloud particles (negative index)
           if(cc(1).eq.myid)then
              npart=npart+1
              tp(npart)=1d-10
              mp(npart)=msink(isink)     ! Mass
              levelp(npart)=levelmin
              idp(npart)=-isink          ! Identity
              xp(npart,1)=xsink(isink,1) ! Position
              xp(npart,2)=xsink(isink,2)
              xp(npart,3)=xsink(isink,3)
              vp(npart,1)=vsink(isink,1) ! Velocity
              vp(npart,2)=vsink(isink,2)
              vp(npart,3)=vsink(isink,3)
           endif

        end do

     end if
  end if

end subroutine init_part
#define TIME_START(cs) call SYSTEM_CLOCK(COUNT=cs)
#define TIME_END(ce) call SYSTEM_CLOCK(COUNT=ce)
#define TIME_SPENT(cs,ce,cr) REAL((ce-cs)/cr)
subroutine load_gadget
  use amr_commons
  use pm_commons
  use gadgetreadfilemod
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::ok
  TYPE(gadgetheadertype) :: gadgetheader
  integer::numfiles
  integer::ifile
  real(dp),dimension(1:nvector,1:3)::xx_dp
  real, dimension(:, :), allocatable:: pos, vel
  integer, dimension(:), allocatable:: ids
  integer::nparticles, arraysize
  integer::i, icpu, ipart, info, np, start
  integer ,dimension(1:ncpu)::npart_cpu,npart_all
  character(LEN=256)::filename
  integer ,dimension(1:nvector)::cc
  integer :: clock_start, clock_end, clock_rate
  integer :: mpi_cs, mpi_ce
  real:: gadgetvfact
  ! Local particle count
  ipart=0
  call SYSTEM_CLOCK(COUNT_RATE=clock_rate)

  if(TRIM(initfile(levelmin)).NE.' ')then
     filename=TRIM(initfile(levelmin))
     ! read first header to get information
     call gadgetreadheader(filename, 0, gadgetheader, ok)
     if(.not.ok) call clean_stop
     numfiles = gadgetheader%numfiles
     gadgetvfact = SQRT(aexp) / gadgetheader%boxsize * aexp / 100.
     do ifile=0,numfiles-1
        call gadgetreadheader(filename, ifile, gadgetheader, ok)
        nparticles = gadgetheader%npart(2)
        allocate(pos(3,nparticles))
        allocate(vel(3,nparticles))
        allocate(ids(nparticles))
        TIME_START(clock_start)
        call gadgetreadfile(filename,ifile,gadgetheader, pos, vel, ids)
        TIME_END(clock_end)
        if(debug) write(*,*) myid, ':Read ', nparticles, ' from gadget file ', ifile, ' in ', &
        TIME_SPENT(clock_start, clock_end, clock_rate)
        start = 1
        TIME_START(clock_start)
#ifndef WITHOUTMPI
        do i=1,nparticles
           xx_dp(1,1) = pos(1,i)/gadgetheader%boxsize
           xx_dp(1,2) = pos(2,i)/gadgetheader%boxsize
           xx_dp(1,3) = pos(3,i)/gadgetheader%boxsize
           call cmp_cpumap(xx_dp,cc,1)
           if(cc(1)==myid)then
#endif
              ipart=ipart+1
              if (ipart .ge. size(mp)) then
                 write(*,*) "For ", myid, ipart, " exceeds ", size(mp)
                 call clean_stop
              end if
              xp(ipart,1:3)=xx_dp(1,1:3)
              vp(ipart,1)=vel(1, i) * gadgetvfact
              vp(ipart,2)=vel(2, i) * gadgetvfact
              vp(ipart,3)=vel(3, i) * gadgetvfact
              mp(ipart)    = 1.d0/gadgetheader%nparttotal(2)
              levelp(ipart)=levelmin
              idp(ipart)   =ids(i)
#ifndef WITHOUTMPI
            endif
        enddo
        TIME_END(clock_end)
        if(debug) write(*,*) myid, ':Processed ', nparticles, ' in ',&
             &  TIME_SPENT(clock_start, clock_end, clock_rate), " ipart now ", ipart
#endif
        deallocate(pos,vel,ids)
     end do

  end if
  npart=ipart
  ! Compute total number of particleclock_rate
  npart_cpu=0; npart_all=0
  npart_cpu(myid)=npart
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  npart_cpu(1)=npart_all(1)
#endif
  do icpu=2,ncpu
     npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
  end do
  write(*,*)'npart=',npart,'/',npart_cpu(ncpu)

end subroutine load_gadget



!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_ncloud_sink
  use amr_commons, only:dp
  use pm_commons, only:ir_cloud,ncloud_sink
  real(dp)::xx,yy,zz,rr
  integer::ii,jj,kk
  ! Compute number of cloud particles
  ncloud_sink=0
  do kk=-2*ir_cloud,2*ir_cloud
     zz=dble(kk)/2.0
     do jj=-2*ir_cloud,2*ir_cloud
        yy=dble(jj)/2.0
        do ii=-2*ir_cloud,2*ir_cloud
           xx=dble(ii)/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=dble(ir_cloud))ncloud_sink=ncloud_sink+1
        end do
     end do
  end do

end subroutine compute_ncloud_sink

