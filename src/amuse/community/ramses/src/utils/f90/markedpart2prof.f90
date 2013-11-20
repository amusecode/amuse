program part2prof
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite surfacique projetee
  ! des particules de matiere noire d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,nstar,nmark,nmark_inrad
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::nx=0,ny=0,ix,iy,ixp1,iyp1,idim,jdim,ncpu_read,n_frw
  integer::nprof=28,irad,ivar,nrad=100
  integer::marked
  integer(kind=8)::nread
  real(KIND=8)::mtot,ddx,ddy,dex,dey,time,time_tot,time_simu,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,rmax=0.5
  real(KIND=8)::xcen=0.5,ycen=0.5,zcen=0.5
  real(KIND=8)::ucen=0.0,vcen=0.0,wcen=0.0
  real(KIND=8)::xx,yy,zz,uu,vv,ww
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin,npart_actual
  real(KIND=8)::xxmin,xxmax,yymin,yymax,dx,dy,deltax,boxlen
  real(KIND=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0
  real(KIND=8)::unit_l,unit_t,unit_d,unit_m,unit_v
  real(KIND=8)::rad2,vol,rprev
  real(KIND=8)::mcumstar,ucumstar,vcumstar,wcumstar,lxcumstar,lycumstar,lzcumstar
  real(KIND=8)::mcumcdm,ucumcdm,vcumcdm,wcumcdm,lxcumcdm,lycumcdm,lzcumcdm

  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:),allocatable::x,v,prof
  real(KIND=8),dimension(:),allocatable::m,age,r
  integer,dimension(:),allocatable::id
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering,format_grille
  character(LEN=128)::nomfich,repository,outfich,filedens,filetype='bin',markfile
  logical::ok,ok_part,periodic=.false.
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  logical::cosmo=.true.
  integer::idcdm=1,iucdm=2,ivcdm=3,iwcdm=4,ilxcdm=5,ilycdm=6,ilzcdm=7
  integer::imcumcdm=8,iucumcdm=9,ivcumcdm=10,iwcumcdm=11,ilxcumcdm=12,ilycumcdm=13,ilzcumcdm=14
  integer::idstar=15,iustar=16,ivstar=17,iwstar=18,ilxstar=19,ilystar=20,ilzstar=21
  integer::imcumstar=22,iucumstar=23,ivcumstar=24,iwcumstar=25,ilxcumstar=26,ilycumstar=27,ilzcumstar=28

  call read_params

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,'("ncpu        =",I11)')ncpu
  read(10,'("ndim        =",I11)')ndim
  read(10,'("levelmin    =",I11)')levelmin
  read(10,'("levelmax    =",I11)')levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,'("boxlen      =",E23.15)')boxlen
  read(10,'("time        =",E23.15)')t
  read(10,'("aexp        =",E23.15)')aexp
  read(10,'("H0          =",E23.15)')h0
  if(h0.eq.1.0)cosmo=.false.
  read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'("unit_t      =",E23.15)')unit_t
  unit_m=unit_d*unit_l**3
  unit_v=unit_l/unit_t
  read(10,*)

  read(10,'("ordering type=",A80)'),ordering
  write(*,'(" ordering type=",A20)'),TRIM(ordering)
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
     time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
     write(*,*)'Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
  else
     time_simu=t
  endif

  !-----------------------
  ! Profile parameters
  !-----------------------
  xmin=MAX(xcen-rmax,0.0d0)
  xmax=MIN(xcen+rmax,1.0d0)
  ymin=MAX(ycen-rmax,0.0d0)
  ymax=MIN(ycen+rmax,1.0d0)
  zmin=MAX(zcen-rmax,0.0d0)
  zmax=MIN(zcen+rmax,1.0d0)

  write(*,*)'time=',t
  write(*,*)'Working array =',nrad
  allocate(r(1:nrad))
  do i=1,nrad
     r(i)=dble(i)*rmax/dble(nrad)
  end do
  allocate(prof(1:nrad,1:nprof))
  prof=0.0d0

  !corbett, removed reference to hilbert ordering, ask Romain, so this code is slightly reduntant but left in
  ncpu_read=ncpu
  do j=1,ncpu
     cpu_list(j)=j
  end do
  !end mod
  write(*,*) 'ncpu,ncpu_read',ncpu,ncpu_read
  npart=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)nstar
     close(1)
     npart=npart+npart2
  end do
  write(*,*)'Found ',npart,' particles.'

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  npart_actual=0
  mtot=0.0d0
  nmark=0
  nmark_inrad=0
  open(20,file=markfile,access='DIRECT',form='unformatted',recl=1)  !corbett: read in markfile here
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
!     write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     allocate(age(1:npart2))
     allocate(id(1:npart2))
     age=0d0
     id=1
     allocate(x(1:npart2,1:ndim2))
     allocate(v(1:npart2,1:ndim2))
     ! Read position
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m/boxlen
     end do
     ! Skip velocity
     do i=1,ndim
        read(1)m
        v(1:npart2,i)=m
     end do
     ! Read mass
     read(1)m
     !if(nstar>0)then !corbett removed
     read(1)id ! Read identity
     if(nstar>0)then !corbett ,pved
        read(1)   ! Skip level
        read(1)age
     endif
     close(1)
     
     do i=1,npart2
        read(20,REC=id(i)+1) marked !corbett here's where I seek to the marked ID
        if(marked<1) cycle !corbett: if it's not marked, then continue to next particle, skipping the rest here
        nmark=nmark+1 !corbett
        rad2=(x(i,1)-xcen)**2+(x(i,2)-ycen)**2+(x(i,3)-zcen)**2
        ok_part=(rad2<rmax**2)
        if(ok_part)then
           nmark_inrad=nmark_inrad+1
           irad=int(dble(nrad)*sqrt(rad2)/rmax)+1
           xx=x(i,1)-xcen
           yy=x(i,2)-ycen
           zz=x(i,3)-zcen
           uu=v(i,1)-ucen/(unit_v/1d5)
           vv=v(i,2)-vcen/(unit_v/1d5)
           ww=v(i,3)-wcen/(unit_v/1d5)
           if(age(i).ne.0.0d0.and.id(i)>0)then
              if(cosmo)then
                 iii=1
                 do while(tau_frw(iii)>age(i).and.iii<n_frw)
                    iii=iii+1
                 end do
                 ! Interpolate time
                 time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                      & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                 time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
              else
                 time=(time_simu-age(i))*unit_t
              end if
              prof(irad,idstar)=prof(irad,idstar)+m(i)
              prof(irad,iustar)=prof(irad,iustar)+m(i)*v(i,1)
              prof(irad,ivstar)=prof(irad,ivstar)+m(i)*v(i,2)
              prof(irad,iwstar)=prof(irad,iwstar)+m(i)*v(i,3)
              prof(irad,ilxstar)=prof(irad,ilxstar)+m(i)*(yy*ww-zz*vv)
              prof(irad,ilystar)=prof(irad,ilystar)-m(i)*(xx*ww-zz*uu)
              prof(irad,ilzstar)=prof(irad,ilzstar)+m(i)*(xx*vv-yy*uu)
           else if(id(i)>0) then
              prof(irad,idcdm)=prof(irad,idcdm)+m(i)
              prof(irad,iucdm)=prof(irad,iucdm)+m(i)*v(i,1)
              prof(irad,ivcdm)=prof(irad,ivcdm)+m(i)*v(i,2)
              prof(irad,iwcdm)=prof(irad,iwcdm)+m(i)*v(i,3)
              prof(irad,ilxcdm)=prof(irad,ilxcdm)+m(i)*(yy*ww-zz*vv)
              prof(irad,ilycdm)=prof(irad,ilycdm)-m(i)*(xx*ww-zz*uu)
              prof(irad,ilzcdm)=prof(irad,ilzcdm)+m(i)*(xx*vv-yy*uu)
           end if
        end if
     end do
     deallocate(x,m,v,age,id)
  end do
  write(*,*) 'number marked particles',nmark
  write(*,*) 'number marked particles within radius cutoff',nmark_inrad
  ! Compute cumulated profiles
  mcumstar=0d0; mcumcdm=0d0
  ucumstar=0d0; vcumstar=0d0; wcumstar=0d0
  ucumcdm=0d0; vcumcdm=0d0; wcumcdm=0d0
  lxcumstar=0d0; lycumstar=0d0; lzcumstar=0d0
  lxcumcdm=0d0; lycumcdm=0d0; lzcumcdm=0d0
  do irad=1,nrad
     mcumstar=mcumstar+prof(irad,idstar)
     ucumstar=ucumstar+prof(irad,iustar)
     vcumstar=vcumstar+prof(irad,ivstar)
     wcumstar=wcumstar+prof(irad,iwstar)
     lxcumstar=lxcumstar+prof(irad,ilxstar)
     lycumstar=lycumstar+prof(irad,ilystar)
     lzcumstar=lzcumstar+prof(irad,ilzstar)

     mcumcdm=mcumcdm+prof(irad,idcdm)
     ucumcdm=ucumcdm+prof(irad,iucdm)
     vcumcdm=vcumcdm+prof(irad,ivcdm)
     wcumcdm=wcumcdm+prof(irad,iwcdm)
     lxcumcdm=lxcumcdm+prof(irad,ilxcdm)
     lycumcdm=lycumcdm+prof(irad,ilycdm)
     lzcumcdm=lzcumcdm+prof(irad,ilzcdm)

     prof(irad,imcumstar)=mcumstar
     prof(irad,iucumstar)=ucumstar
     prof(irad,ivcumstar)=vcumstar
     prof(irad,iwcumstar)=wcumstar
     prof(irad,ilxcumstar)=lxcumstar
     prof(irad,ilycumstar)=lycumstar
     prof(irad,ilzcumstar)=lzcumstar

     prof(irad,imcumcdm)=mcumcdm
     prof(irad,iucumcdm)=ucumcdm
     prof(irad,ivcumcdm)=vcumcdm
     prof(irad,iwcumcdm)=wcumcdm
     prof(irad,ilxcumcdm)=lxcumcdm
     prof(irad,ilycumcdm)=lycumcdm
     prof(irad,ilzcumcdm)=lzcumcdm

  end do
  
  ! Convert profiles into proper astro units
  rprev=0d0
  do irad=1,nrad
     r(irad)=r(irad)*boxlen
     vol=4./3.*3.1415926*(r(irad)**3-rprev**3)
     ! Stars
     if(prof(irad,idstar)>0.0)then
        prof(irad,iustar)=prof(irad,iustar)/prof(irad,idstar)*unit_v/1d5
        prof(irad,ivstar)=prof(irad,ivstar)/prof(irad,idstar)*unit_v/1d5
        prof(irad,iwstar)=prof(irad,iwstar)/prof(irad,idstar)*unit_v/1d5
        prof(irad,ilxstar)=prof(irad,ilxstar)/prof(irad,idstar)*unit_v/1d5/r(irad)
        prof(irad,ilystar)=prof(irad,ilystar)/prof(irad,idstar)*unit_v/1d5/r(irad)
        prof(irad,ilzstar)=prof(irad,ilzstar)/prof(irad,idstar)*unit_v/1d5/r(irad)
     endif
     prof(irad,idstar)=prof(irad,idstar)/vol*unit_d/1.66d-24
     if(prof(irad,imcumstar)>0.0)then
        prof(irad,iucumstar)=prof(irad,iucumstar)/prof(irad,imcumstar)*unit_v/1d5
        prof(irad,ivcumstar)=prof(irad,ivcumstar)/prof(irad,imcumstar)*unit_v/1d5
        prof(irad,iwcumstar)=prof(irad,iwcumstar)/prof(irad,imcumstar)*unit_v/1d5
        prof(irad,ilxcumstar)=prof(irad,ilxcumstar)/prof(irad,imcumstar)*unit_v/1d5/r(irad)
        prof(irad,ilycumstar)=prof(irad,ilycumstar)/prof(irad,imcumstar)*unit_v/1d5/r(irad)
        prof(irad,ilzcumstar)=prof(irad,ilzcumstar)/prof(irad,imcumstar)*unit_v/1d5/r(irad)
     endif
     prof(irad,imcumstar)=sqrt(6.67e-8*prof(irad,imcumstar)*unit_m/r(irad)/unit_l)/1d5
     ! Dark matter
     if(prof(irad,idcdm)>0.0)then
        prof(irad,iucdm)=prof(irad,iucdm)/prof(irad,idcdm)*unit_v/1d5
        prof(irad,ivcdm)=prof(irad,ivcdm)/prof(irad,idcdm)*unit_v/1d5
        prof(irad,iwcdm)=prof(irad,iwcdm)/prof(irad,idcdm)*unit_v/1d5
        prof(irad,ilxcdm)=prof(irad,ilxcdm)/prof(irad,idcdm)*unit_v/1d5/r(irad)
        prof(irad,ilycdm)=prof(irad,ilycdm)/prof(irad,idcdm)*unit_v/1d5/r(irad)
        prof(irad,ilzcdm)=prof(irad,ilzcdm)/prof(irad,idcdm)*unit_v/1d5/r(irad)
     endif
     prof(irad,idcdm)=prof(irad,idcdm)/vol*unit_d/1.66d-24
     if(prof(irad,imcumcdm)>0.0)then
        prof(irad,iucumcdm)=prof(irad,iucumcdm)/prof(irad,imcumcdm)*unit_v/1d5
        prof(irad,ivcumcdm)=prof(irad,ivcumcdm)/prof(irad,imcumcdm)*unit_v/1d5
        prof(irad,iwcumcdm)=prof(irad,iwcumcdm)/prof(irad,imcumcdm)*unit_v/1d5
        prof(irad,ilxcumcdm)=prof(irad,ilxcumcdm)/prof(irad,imcumcdm)*unit_v/1d5/r(irad)
        prof(irad,ilycumcdm)=prof(irad,ilycumcdm)/prof(irad,imcumcdm)*unit_v/1d5/r(irad)
        prof(irad,ilzcumcdm)=prof(irad,ilzcumcdm)/prof(irad,imcumcdm)*unit_v/1d5/r(irad)
     endif
     prof(irad,imcumcdm)=sqrt(6.67e-8*prof(irad,imcumcdm)*unit_m/r(irad)/unit_l)/1d5
     rprev=r(irad)
  end do

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  open(unit=10,file=TRIM(nomfich)//".dark",form='formatted')
  write(10,'(A130)')" r(kpc)      n_d(H/cc)   vc_d(H/cc)  mc_d(Msol)  cu_d(km/s)  cv_d(km/s)  cw_d(km/s)  cl_d(km/s)  lx_d        ly_d        lz_d      "
  do irad=1,nrad
     write(10,999)r(irad)*unit_l/3.08d21,prof(irad,idcdm),prof(irad,imcumcdm),prof(irad,imcumcdm)**2*r(irad)*unit_l/6.67e-8*1d10/2d33 &
          & ,prof(irad,iucumcdm),prof(irad,ivcumcdm),prof(irad,iwcumcdm),sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2) &
          & ,prof(irad,ilxcumcdm)/sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2+1d-30) &
          & ,prof(irad,ilycumcdm)/sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2+1d-30) &
          & ,prof(irad,ilzcumcdm)/sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2+1d-30)
  end do
  close(10)
  open(unit=10,file=TRIM(nomfich)//".star",form='formatted')
  write(10,'(A130)')" r(kpc)      n_*(H/cc)   vc_*(H/cc)  mc_*(Msol)  cu_*(km/s)  cv_*(km/s)  cw_*(km/s)  cl_*(km/s)  lx_*        ly_*        lz_*      "
  do irad=1,nrad
     write(10,999)r(irad)*unit_l/3.08d21,prof(irad,idstar),prof(irad,imcumstar),prof(irad,imcumstar)**2*r(irad)*unit_l/6.67e-8*1d10/2d33 &
          & ,prof(irad,iucumstar),prof(irad,ivcumstar),prof(irad,iwcumstar),sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2) &
          & ,prof(irad,ilxcumstar)/sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2+1e-30) &
          & ,prof(irad,ilycumstar)/sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2+1e-30) &
          & ,prof(irad,ilzcumstar)/sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2+1e-30)
  end do
  close(10)
999 format(30(1PE10.3,2X))

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
         print *, 'usage: part2prof -inp  input_dir'
         print *, '                 -mar  mark_file'
         print *, '                 -out  output_file'
         print *, '                 [-xce xcen] '
         print *, '                 [-yce ycen] '
         print *, '                 [-zce zcen] '
         print *, '                 [-uce ucen] '
         print *, '                 [-vce vcen] '
         print *, '                 [-wce wcen] '
         print *, '                 [-rma rmax] '
         print *, '                 [-nra nrad] '
         print *, '                 [-per flag] '
         print *, 'ex: part2prof -mar mark_file -inp output_00001 -out map.dat'// &
              &   ' -xce 0.1 -yce 0.2 -zce 0.2 -rma 0.1 -nra 100'
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
         case ('-mar')
            markfile = trim(arg)
         case ('-out')
            outfich = trim(arg)
         case ('-xce')
            read (arg,*) xcen
         case ('-yce')
            read (arg,*) ycen
         case ('-zce')
            read (arg,*) zcen
         case ('-uce')
            read (arg,*) ucen
         case ('-vce')
            read (arg,*) vcen
         case ('-wce')
            read (arg,*) wcen
         case ('-nra')
            read (arg,*) nrad
         case ('-rma')
            read (arg,*) rmax
         case ('-per')
            read (arg,*) periodic
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params

  end program part2prof

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif

end subroutine title

!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
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




