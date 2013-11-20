program ramses2tipsy

  use random 
  use io_ramses

  implicit none

  integer::ndim,n,i,j,k,iii,twotondim,ncoarse,indcell
  integer::ivar,nvar,ncpu,lmax=100,levelmin,iskip
  integer::nx,ny,nz
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ipos
  integer::ngridmax,nstep_coarse,icpu
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  real(KIND=8)::gamma,smallr,smallc,gammah
  real(KIND=8)::boxlen,boxlen2
  real(KIND=8)::t,aexp,h0,t2,aexp2,hexp2
  real(KIND=8)::omega_m,omega_l,omega_k,omega_b
  real(kind=8)::omega_m2,omega_l2,omega_k2,omega_b2
  real(kind=8)::scale_l,scale_d,scale_t
  real(kind=8)::metmax=0d0

  integer::nx_sample=0,ny_sample=0,ngridtot
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,lmin,nboundary,ngrid_current,lllmin
  integer::ix,iy,iz,ndom,impi,bit_length,maxdom,ii,jj,kk
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,ddx,dxline,ddy,dex,dey,weight,msph=0.d0
  real(KIND=8)::xmin=0,xmax=-1,ymin=0,ymax=-1,zmin=0,zmax=-1,mdm=0.d0,mres=0.d0

  integer::ilevel,ncpu_read,three
  real(kind=8)deltax
  character(len=3)::typ='all'
  character(LEN=5)::nchar
  character(LEN=80)::ordering
  character(LEN=128)::nomfich
  character(LEN=128)::repository,outfich,filetype='bin'
  character(LEN=13)::string
  logical::ok,ok_part,ok_cell,do_max,do_id=.false.
  real(kind=8),dimension(:),allocatable::bound_key,xdp
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list

  integer::ndummypart,nmin=0,nmax=0,nold,nnold,ndummyold
  integer::partcount,respart,denspartcount
  real(KIND=8)::dummy,partmass,volume,facdens=0.d0,averdens
  integer::delm=0,levelsel

  integer::npart,nstar_tot,nsink,npart_tot,nsink_tot
  integer::npart_actual,ndm_actual,nstar_actual,npartnow
  character(LEN=5)::nn
  logical::hydrok=.false.,partok=.false.,metal=.false.,metgas=.false.
  logical::cosmo=.false.,star=.false.,sink=.false.,mhd=.false.,gas=.true.
  
  integer::n_frw
  real(KIND=8)::time,time_tot,time_simu,time_uni
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw

  integer ,dimension(1:1,1:IRandNumSize)::allseed
  integer ,dimension(1:IRandNumSize)::localseed
  integer::iseed=0,poisson


  call read_params

  !-----------------------------------------------
  ! Reading files in RAMSES format
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  hydrok=ok
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
  endif
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  partok=ok
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)ngrid_current
  read(10)boxlen
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)msph
  close(10)
  twotondim=2**ndim

  ! Rescaling mass_sph
  msph=msph/boxlen**3

  ! Default values for box size
  if(xmax<0)xmax=boxlen
  if(ymax<0)ymax=boxlen
  if(zmax<0)zmax=boxlen

  if(hydrok)then
     ! Read nvarh from the Hydro file
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
     open(unit=12,file=nomfich,status='old',form='unformatted')
     read(12)
     read(12)nvarh
     read(12)
     read(12)
     read(12)
     read(12)gamma
     close(12)
     if(ndim==2)then
        if((.not.mhd).and.(nvarh<=4))then
           metgas=.false.
           metal=.false.
        end if
        if((.not.mhd).and.(nvarh>4))then
           metgas=.true.
        end if
        if((mhd).and.(nvarh<=11))then
           metgas=.false.
           metal=.false.
        end if
        if((mhd).and.(nvarh>11))then
           metgas=.true.
        end if
     end if
     if(ndim==3)then
        if((.not.mhd).and.(nvarh<=5))then
           metgas=.false.
           metal=.false.
        end if
        if((.not.mhd).and.(nvarh>5))then
           metgas=.true.
        end if
        if((mhd).and.(nvarh<=11))then
           metgas=.false.
           metal=.false.
        end if
        if((mhd).and.(nvarh>11))then
           metgas=.true.
        end if
     end if
  end if

  if(partok)then
     npart_tot=0
     nsink_tot=0
     do i=1,ncpu
        ! Read number of particles from the Part file
        write(nn,'(I5.5)')i
        nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//nn
        open(unit=11,file=nomfich,status='old',form='unformatted')
        read(11)
        read(11)
        read(11)npart
        read(11)
        read(11)nstar_tot
        read(11)
        read(11)
        read(11)nsink
        close(11)
        nsink_tot=nsink_tot+nsink
        npart_tot=npart_tot+npart
     enddo
     if(nsink_tot>0)sink=.true.
     if(nstar_tot>0)star=.true.
     if(metgas.and.star)metal=.true.
  else
     sink=.false.
     star=.false.
     metal=.false.
  end if

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,*)
  read(10,*)
  read(10,'(A13,I11)')string,levelmin
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'(A13,E23.15)')string,t
  read(10,'(A13,E23.15)')string,aexp
  read(10,'(A13,E23.15)')string,h0
  read(10,'(A13,E23.15)')string,omega_m
  read(10,'(A13,E23.15)')string,omega_l
  read(10,'(A13,E23.15)')string,omega_k
  read(10,'(A13,E23.15)')string,omega_b
  read(10,'(A13,E23.15)')string,scale_l
  read(10,'(A13,E23.15)')string,scale_d
  read(10,'(A13,E23.15)')string,scale_t
  read(10,*)
  read(10,'("ordering type=",A80)'),ordering
  read(10,*)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  close(10)

  lmax=max(min(lmax,nlevelmax),1)
  levelsel=levelmin+delm

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,lmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lllmin=ilevel
     bit_length=lllmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
     endif

     dkey=(dble(2**(lmax+1)/dble(maxdom)))**ndim

     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax

     do i=1,ndom
        if(bit_length>0)then
           ii=idom(i)
           jj=jdom(i)
           kk=kdom(i)
           call hilbert3d(ii,jj,kk,order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0d0)*dkey
     end do
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (bound_key(impi-1).le.bounding_min(i).and.&
                &bound_key(impi).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (bound_key(impi-1).lt.bounding_max(i).and.&
                &bound_key(impi).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if


  if(cosmo)then
     !-----------------------
     ! Cosmological model
     !-----------------------
     ! Allocate look-up tables
     n_frw=1000
     allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
     allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
     
     ! Compute Friedman model look up table
     write(*,*)'Computing Friedman model'
     call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
          & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)

     ! Find neighboring expansion factors
     i=1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i=i+1
     end do
     ! Interploate time
     time_simu=t_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
          & t_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
     write(*,*)'Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
     time_uni=(time_tot+time_simu)! ADD THIS TO CONVERT TO GYR: /(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
  else
     time_simu=t
     time_uni=t!*scale_t
  endif

  npart_actual=0
  denspartcount=0
  nstar_actual=0

  if(partok) then
     call readpart(ncpu,ncpu_read,cpu_list,ndim,repository,metal,star,sink,&
          & lmin,lmax,xmin,xmax,ymin,ymax,zmin,zmax,nmin,nmax,npart_actual,&
          & ndm_actual,nstar_actual)
     !NOTE: READPART SELECTS ONLY PARTICLES WITH ID>0.

     write(*,*)'Number of DM particles in the selected box: ', npart_actual-nstar_actual
     write(*,*)'Number of star particles in the selected box: ', nstar_actual

     do i=1,npart_actual
        if((ageout(i).ne.0.0d0))then
           if(cosmo)then
              iii=1 ! Compute star formation time 
              do while(tau_frw(iii)>ageout(i).and.iii<n_frw)
                 iii=iii+1
              end do
              time=t_frw(iii)*(ageout(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                   & t_frw(iii-1)*(ageout(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
              time=(time_tot+time)!ADD THIS TO CONVERT TO YR: /(h0*1d5/3.08d24)/(365.*24.*3600.)
              ageout(i)=time ! Replace age with formation time
           else
              ageout(i)=ageout(i)!ADD THIS TO CONVERT TO YR: *scale_t/(365.*24.*3600.)
           end if
        end if
     end do

     mdm=msph ! IN THE OLD VERSION: mdm*omega_b/(omega_m-omega_b)
     if(mres.ne.0.d0)mdm=max(mdm,mres)
     if(mdm.eq.0)then
        write(*,*)'Use option -mre, please! STOP!'
        stop
     end if
  else
     mdm=msph ! IN THE OLD VERSION: 1d0/(2d0**(3d0*dble(lmax)))
     if(mres.ne.0.d0)mdm=max(mdm,mres)
     if(mdm.eq.0)then
        write(*,*)'Use option -mre, please! STOP!'
        stop
     end if
  end if

  if(hydrok.and.gas)then
     call gaspart3(ncpu,ncpu_read,cpu_list,repository,ordering,&
          & ndummypart,facdens,levelsel,lmax,xmin,xmax,ymin,&
          & ymax,zmin,zmax,mdm,partmass,averdens,&
          & denspartcount)
     nmin=1
     nmax=denspartcount
  end if

  !-------------------------------------------------------------
  !  Writing output tipsy file in ascii format
  !-------------------------------------------------------------
  if(filetype .EQ. 'ascii')then

  write(*,*)'Outputing data in tipsy ASCII format' 

  if(gas)then
     
     open(66,file=outfich,status='unknown',form='formatted')
     if(do_id)then
        open(55,file='partID_'//TRIM(nchar),status='unknown',form='formatted')
     endif
     
     three=3
     
     !HEADER
     write(66,*)npart_actual+denspartcount,denspartcount,nstar_actual
     write(66,*)three
     write(66,*)time_uni
     
     if(do_id)then
        write(55,*)npart_actual+denspartcount,denspartcount,nstar_actual
        write(55,*)three
        write(55,*)time_uni
     endif
     
     !MASSES AND IDS
     if(hydrok)then
        do i=1,nmax-nmin+1
           write(66,*)partmass
        end do
     end if
     if(partok)then
        do i=1,npart_actual
           if(.not.star)write(66,*)mout(i)
           if((star.or.sink).and.ageout(i)==0.d0)write(66,*)mout(i)
        enddo
        if(do_id)then
           do i=1,npart_actual
              if(.not.star)write(55,*)idout(i)
              if((star.or.sink).and.ageout(i)==0.d0)write(55,*)idout(i)
           enddo
        endif
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)mout(i)
           enddo
           if(do_id)then
              do i=1,npart_actual
                 if(ageout(i)/=0.d0)write(55,*)idout(i)
              enddo
           endif
        endif
     end if
     
     !X COORDINATE
     if(hydrok)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens)write(66,*)xp(i,1)
        end do
     end if
     if(partok)then
        do i=1,npart_actual
           if(.not.star)write(66,*)xout(i,1)
           if((star.or.sink).and.ageout(i)==0.d0)write(66,*)xout(i,1)
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)xout(i,1)
           enddo
        endif
     endif
     
     !Y COORDINATE
     if(hydrok)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens.and.ndim>=2)write(66,*)xp(i,2)
           if(varp(i,1)>=facdens*averdens.and.ndim<2)write(66,*)(ymin+ymax)/2
        end do
     end if
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=2)write(66,*)xout(i,2)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=2)write(66,*)xout(i,2)
           if((.not.star).and.ndim<2)write(66,*)(ymin+ymax)/2
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<2)write(66,*)(ymin+ymax)/2
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=2)write(66,*)xout(i,2)
              if(ageout(i)/=0.d0.and.ndim<2)write(66,*)(ymin+ymax)/2
           enddo
        endif
     end if
     
     !Z COORDINATE
     if(hydrok)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens.and.ndim>=3)write(66,*)xp(i,3)
           if(varp(i,1)>=facdens*averdens.and.ndim<3)write(66,*)(zmin+zmax)/2
        end do
     end if
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=3)write(66,*)xout(i,3)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=3)write(66,*)xout(i,3)
           if((.not.star).and.ndim<3)write(66,*)(zmin+zmax)/2
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<3)write(66,*)(zmin+zmax)/2
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=3)write(66,*)xout(i,3)
              if(ageout(i)/=0.d0.and.ndim<3)write(66,*)(zmin+zmax)/2
           enddo
        endif
     end if
     
     !VELOCITY 
     
     dummy=0.d0
     
     !V_X
     if(hydrok)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,2)
        end do
     end if
     if(partok)then
        do i=1,npart_actual
           if(.not.star)write(66,*)vout(i,1)
           if((star.or.sink).and.ageout(i)==0.d0)write(66,*)vout(i,1)
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)vout(i,1)
           enddo
        endif
     end if
     
     !V_Y
     if(hydrok)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens.and.ndim>=2)write(66,*)varp(i,3)
           if(varp(i,1)>=facdens*averdens.and.ndim<2)write(66,*)dummy
        end do
     end if
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=2)write(66,*)vout(i,2)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=2)write(66,*)vout(i,2)
           if((.not.star).and.ndim<2)write(66,*)dummy
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<2)write(66,*)dummy
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=2)write(66,*)vout(i,2)
              if(ageout(i)/=0.d0.and.ndim<2)write(66,*)dummy
           enddo
        endif
     endif
     
     !V_Z
     if(hydrok)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens.and.ndim>=3)write(66,*)varp(i,4)
           if(varp(i,1)>=facdens*averdens.and.ndim<3)write(66,*)dummy
        end do
     end if
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=3)write(66,*)vout(i,3)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=3)write(66,*)vout(i,3)
           if((.not.star).and.ndim<3)write(66,*)dummy
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<3)write(66,*)dummy
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=3)write(66,*)vout(i,3)
              if(ageout(i)/=0.d0.and.ndim<3)write(66,*)dummy
           enddo
        endif
     endif
     
     !DUMMY GRAVITATIONAL SOFTENING FOR DARK AND STARS
     if(partok)then
        dummy=boxlen/2**lmax !THIS IS A DUMMY VALUE: IT CORRESPONDS TO THE CELL SIZE AT THE MAXIMUM LEVEL.
        do i=1,ndm_actual
           write(66,*)dummy 
        end do
        if(star.and.nstar_actual>0)then
           do i=1,nstar_actual
              write(66,*)dummy
           end do
        end if
     end if

     !GAS DENSITY, TEMPERATURE, DUMMY SPH SMOOTHING LENGTH & GAS METALLICITY. 
     if(hydrok)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,1)
        end do
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens.and.ndim==3.and.(.not.mhd))write(66,*)varp(i,5)/(gamma-1.d0)/varp(i,1) !THIS IS P/RHO=(k_b*T)/(mu*m_h)
           if(varp(i,1)>=facdens*averdens.and.ndim==2.and.(.not.mhd))write(66,*)varp(i,4)/(gamma-1.d0)/varp(i,1) !THIS IS P/RHO=(k_b*T)/(mu*m_h)
           if(varp(i,1)>=facdens*averdens.and.ndim==1.and.(.not.mhd))write(66,*)varp(i,3)/(gamma-1.d0)/varp(i,1) !THIS IS P/RHO=(k_b*T)/(mu*m_h)
           if(varp(i,1)>=facdens*averdens.and.(mhd))write(66,*)varp(i,11)/(gamma-1.d0)/varp(i,1) !THIS IS P/RHO=(k_b*T)/(mu*m_h)
        end do
        dummy=boxlen/2**lmax !THIS IS A DUMMY VALUE: IT CORRESPONDS TO THE CELL SIZE AT THE MAXIMUM LEVEL.
        do i=1,nmax-nmin+1
           write(66,*)dummy 
        end do
        do i=1,nmax-nmin+1
           if(varp(i,1)>=facdens*averdens.and.metgas.and.ndim==3.and.(.not.mhd))write(66,*)varp(i,6)
           if(varp(i,1)>=facdens*averdens.and.metgas.and.ndim==2.and.(.not.mhd))write(66,*)varp(i,5)
           if(varp(i,1)>=facdens*averdens.and.metgas.and.ndim==1.and.(.not.mhd))write(66,*)varp(i,4)
           if(varp(i,1)>=facdens*averdens.and.metgas.and.(mhd))write(66,*)varp(i,12)
        end do
     end if
     
     if(star.and.nstar_actual>0)then
        if(metal)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)metout(i)
           end do
        else
           dummy=0.d0
           do i=1,nstar_actual
              write(66,*)dummy
           end do
        end if
        do i=1,npart_actual
           if(ageout(i)/=0.d0)write(66,*)ageout(i)
        end do
     end if
     
     dummy=1.d0
     do i=1,npart_actual+denspartcount
        write(66,*)dummy
     end do
     
     close(66)

     if(hydrok)deallocate(xp,varp)
     
     if(partok)then
        deallocate(xout,vout,mout,idout)
        if(star.and.nstar_actual>0)then
           deallocate(ageout)
           if(metal)deallocate(metout)
        end if
     end if

  else

     open(66,file=outfich,status='unknown',form='formatted')
     if(do_id)then
        open(55,file='partID_'//TRIM(nchar),status='unknown',form='formatted')
     endif

     three=3
     
     !HEADER
     write(66,*)npart_actual,0,nstar_actual
     write(66,*)three
     write(66,*)time_uni
     
     if(do_id)then
        write(55,*)npart_actual,0,nstar_actual
        write(55,*)three
        write(55,*)time_uni
     endif
     
     !MASSES AND IDS
     if(partok)then
        do i=1,npart_actual
           if(.not.star)write(66,*)mout(i)
           if((star.or.sink).and.ageout(i)==0.d0)write(66,*)mout(i)
        enddo
        if(do_id)then
           do i=1,npart_actual
              if(.not.star)write(55,*)idout(i)
              if((star.or.sink).and.ageout(i)==0.d0)write(55,*)idout(i)
           enddo
        endif
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)mout(i)
           enddo
           if(do_id)then
              do i=1,npart_actual
                 if(ageout(i)/=0.d0)write(55,*)idout(i)
              enddo
           endif
        endif
     end if
     
     !X COORDINATE
     if(partok)then
        do i=1,npart_actual
           if(.not.star)write(66,*)xout(i,1)
           if((star.or.sink).and.ageout(i)==0.d0)write(66,*)xout(i,1)
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)xout(i,1)
           enddo
        endif
     endif
     
     !Y COORDINATE
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=2)write(66,*)xout(i,2)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=2)write(66,*)xout(i,2)
           if((.not.star).and.ndim<2)write(66,*)(ymin+ymax)/2
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<2)write(66,*)(ymin+ymax)/2
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=2)write(66,*)xout(i,2)
              if(ageout(i)/=0.d0.and.ndim<2)write(66,*)(ymin+ymax)/2
           enddo
        endif
     end if
     
     !Z COORDINATE
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=3)write(66,*)xout(i,3)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=3)write(66,*)xout(i,3)
           if((.not.star).and.ndim<3)write(66,*)(zmin+zmax)/2
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<3)write(66,*)(zmin+zmax)/2
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=3)write(66,*)xout(i,3)
              if(ageout(i)/=0.d0.and.ndim<3)write(66,*)(zmin+zmax)/2
           enddo
        endif
     end if
     
     !VELOCITY 
     
     dummy=0.d0
     
     !V_X
     if(partok)then
        do i=1,npart_actual
           if(.not.star)write(66,*)vout(i,1)
           if((star.or.sink).and.ageout(i)==0.d0)write(66,*)vout(i,1)
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)vout(i,1)
           enddo
        endif
     end if
     
     !V_Y
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=2)write(66,*)vout(i,2)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=2)write(66,*)vout(i,2)
           if((.not.star).and.ndim<2)write(66,*)dummy
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<2)write(66,*)dummy
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=2)write(66,*)vout(i,2)
              if(ageout(i)/=0.d0.and.ndim<2)write(66,*)dummy
           enddo
        endif
     endif
     
     !V_Z
     if(partok)then
        do i=1,npart_actual
           if((.not.star).and.ndim>=3)write(66,*)vout(i,3)
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim>=3)write(66,*)vout(i,3)
           if((.not.star).and.ndim<3)write(66,*)dummy
           if((star.or.sink).and.ageout(i)==0.d0.and.ndim<3)write(66,*)dummy
        enddo
        if(star.and.nstar_actual>0)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0.and.ndim>=3)write(66,*)vout(i,3)
              if(ageout(i)/=0.d0.and.ndim<3)write(66,*)dummy
           enddo
        endif
     endif
     
     !DUMMY GRAVITATIONAL SOFTENING FOR DARK AND STARS
     if(partok)then
        dummy=boxlen/2**lmax !THIS IS A DUMMY VALUE: IT CORRESPONDS TO THE CELL SIZE AT THE MAXIMUM LEVEL.
        do i=1,ndm_actual
           write(66,*)dummy 
        end do
        if(star.and.nstar_actual>0)then
           do i=1,nstar_actual
              write(66,*)dummy
           end do
        end if
     end if
     
     if(star.and.nstar_actual>0)then
        if(metal)then
           do i=1,npart_actual
              if(ageout(i)/=0.d0)write(66,*)metout(i)
           end do
        else
           dummy=0.d0
           do i=1,nstar_actual
              write(66,*)dummy
           end do
        end if
        do i=1,npart_actual
           if(ageout(i)/=0.d0)write(66,*)ageout(i)
        end do
     end if
     
     dummy=1.d0
     do i=1,npart_actual
        write(66,*)dummy
     end do
     
     close(66)

  end if

  !-------------------------------------------------------------
  !  Writing output tipsy file in binary format
  !-------------------------------------------------------------
  else

     write(*,*)'Outputing data in tipsy BINARY format' 
     
     open(66,file=outfich,status='unknown',form='unformatted',access='direct',recl=1)
     
     three=3
     
     !HEADER
     write(66,rec=1)real(time_uni)
     write(66,rec=2)real(time_uni)
     write(66,rec=3)npart_actual+denspartcount
     write(66,rec=4)three
     write(66,rec=5)denspartcount
     write(66,rec=6)ndm_actual
     write(66,rec=7)nstar_actual

     write(*,*)'Header done'
     
     !GAS PARTICLES
     if(hydrok)then
        do i=1,denspartcount
           iskip=8+(i-1)*12
           write(66,rec=iskip+1)real(partmass)
           write(66,rec=iskip+2)real(xp(i,1))
           write(66,rec=iskip+3)real(xp(i,2))
           write(66,rec=iskip+4)real(xp(i,3))
           write(66,rec=iskip+5)real(varp(i,2))
           write(66,rec=iskip+6)real(varp(i,3))
           write(66,rec=iskip+7)real(varp(i,4))
           write(66,rec=iskip+8)real(varp(i,1))
           write(66,rec=iskip+9)real(varp(i,5)/varp(i,1)/(gamma-1.0))
           write(66,rec=iskip+10)real(boxlen/2**lmax)
           if(metal)then
              write(66,rec=iskip+11)real(varp(i,6))
           else
              write(66,rec=iskip+11)real(0.*boxlen)
           endif
           write(66,rec=iskip+12)real(0.*boxlen)
        end do

        write(*,*)'Gas particles done'

     endif

     !DM PARTICLES
     if(partok)then
        if(ndm_actual>0)then
           j=1
           do i=1,npart_actual
              iskip=8+denspartcount*12+(j-1)*9
              if(ageout(i)==0)then
                 write(66,rec=iskip+1)real(mout(i))
                 write(66,rec=iskip+2)real(xout(i,1))
                 write(66,rec=iskip+3)real(xout(i,2))
                 write(66,rec=iskip+4)real(xout(i,3))
                 write(66,rec=iskip+5)real(vout(i,1))
                 write(66,rec=iskip+6)real(vout(i,2))
                 write(66,rec=iskip+7)real(vout(i,3))
                 write(66,rec=iskip+10)real(boxlen/2**lmax)
                 write(66,rec=iskip+11)real(0.*boxlen)
                 j=j+1
              endif
           end do

           write(*,*)'Dark matter particles done'

        endif
     endif

     !STAR PARTICLES
     if(partok)then
        if(nstar_actual>0)then
           j=1
           do i=1,npart_actual
              iskip=8+denspartcount*12+ndm_actual*9+(j-1)*11
              if(ageout(i)/=0)then
                 write(66,rec=iskip+1)real(mout(i))
                 write(66,rec=iskip+2)real(xout(i,1))
                 write(66,rec=iskip+3)real(xout(i,2))
                 write(66,rec=iskip+4)real(xout(i,3))
                 write(66,rec=iskip+5)real(vout(i,1))
                 write(66,rec=iskip+6)real(vout(i,2))
                 write(66,rec=iskip+7)real(vout(i,3))
                 if(metal)then
                    write(66,rec=iskip+11)real(metout(i))
                 else
                    write(66,rec=iskip+11)real(0.*boxlen)
                 endif
                 write(66,rec=iskip+11)real(ageout(i))
                 write(66,rec=iskip+10)real(boxlen/2**lmax)
                 write(66,rec=iskip+11)real(0.*boxlen)
                 j=j+1
              endif
           end do

           write(*,*)'New star particles done'

        endif
     endif

     close(66)

     open(66,file=outfich,status='unknown',form='unformatted',access='direct',recl=2)
     write(66,rec=1)time_uni
     close(66)

  endif

  write(*,*)'File dump completed'

  if(hydrok)deallocate(xp,varp)
  
  if(partok)then
     deallocate(xout,vout,mout,idout)
     if(star.and.nstar_actual>0)then
        deallocate(ageout)
        if(metal)deallocate(metout)
     end if
  end if

contains

  subroutine read_params
    
    implicit none
    
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok
    
    n = iargc()
    if (n < 4) then
       print *, 'usage: ramses2tipsy -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-mre gas particle mass] '
       print *, '                 [-fil filetype=bin or ascii] '
       print *, '                 [-cos cosmo run?] '
       print *, '                 [-pid store particle id?] '
       print *, '                 [-mhd mhd?] '
       print *, '                 [-gas store gas?] '
       print *, 'ex: ramses2tipsy -inp output_00001 -out cube.dat'// &
            &   ' -xmi 0.1 -xma 0.7'
       print *, ' '

       stop
    end if
    
    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case ('-inp')
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
       case ('-xmi')
          read (arg,*) xmin
       case ('-xma')
          read (arg,*) xmax
       case ('-ymi')
          read (arg,*) ymin
       case ('-yma')
          read (arg,*) ymax
       case ('-zmi')
          read (arg,*) zmin
       case ('-zma')
          read (arg,*) zmax
       case ('-cos') 
          read (arg,*) cosmo
       case ('-mhd') 
          read (arg,*) mhd
       case ('-fil') 
          filetype = trim(arg)
       case ('-pid') 
          read (arg,*) do_id
       case ('-gas') 
          read (arg,*) gas
       case ('-mre') 
          read (arg,*) mres
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params

end program ramses2tipsy

!================================================================
!================================================================
!================================================================
!================================================================

subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min,age_tot
  real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

!  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
!     write(*,*)'Error: non-physical cosmological constants'
!     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
!     write(*,*)'The sum must be equal to 1.0, but '
!     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
!     stop
!  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

  age_tot=-t
  write(*,666)-t
  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if
  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
