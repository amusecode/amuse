program amr2prof
  use io_ramses
  use random

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  implicit none
  integer::ndim,ncell,n,i,j,k,twotondim,ncoarse,type=0,domax=0
  integer::ivar,nvar,ncpu,ncpuh,lmax=120,nboundary,ngrid_current
  integer::nx=0,ny=0,nz=0,ilevel,idim,jdim,kdim,icell
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer::ngridmax,nstep_coarse,icpu,ncpu_read
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  integer::iseed=0,irad
  real::gamma,smallr,smallc,gammah
  real::boxlen,boxlen2
  real::t,aexp,hexp,t2,aexp2,hexp2
  real::omega_m,omega_l,omega_k,omega_b
  real::scale_l,scale_d,scale_t
  real::omega_m2,omega_l2,omega_k2,omega_b2
  integer ,dimension(1:1,1:IRandNumSize)::allseed
  integer ,dimension(1:IRandNumSize)::localseed
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax,nrad=100
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin,levelmin,levelmax
  integer::ix,iy,iz,ixp1,iyp1,izp1,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,dummy,h0,dv
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,rmin=1e-10,rmax=0.5
  real(KIND=8)::xcen=0.5,ycen=0.5,zcen=0.5
  real(KIND=8)::ucen=0.0,vcen=0.0,wcen=0.0
  real(KIND=8)::xx,yy,zz,uu,vv,ww,rr
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez
  real(KIND=8)::unit_l,unit_t,unit_d,unit_m,unit_v
  real(KIND=8)::rad2,vol,rprev,nH,tt,hotcum
  real(KIND=8)::mcum,ucum,vcum,wcum,mcold
  real(KIND=8)::lxcum,lycum,lzcum,lxccum,lyccum,lzccum
  real(KIND=8),dimension(:),allocatable::x,y,z,r
  real(KIND=8),dimension(:,:),allocatable::var,prof
  integer,dimension(:),allocatable::l
  real(kind=4),dimension(:,:,:),allocatable::toto
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  character(LEN=5)::nchar
  logical::ok,ok_part,ok_cell,cosmo
  integer::id=1,iu=2,iv=3,iw=4,ilx=5,ily=6,ilz=7
  integer::imcum=8,iucum=9,ivcum=10,iwcum=11
  integer::ilxcum=12,ilycum=13,ilzcum=14,ip=15,imet=16,ihot=17
  integer::ilxc=18,ilyc=19,ilzc=20,ilxccum=21,ilyccum=22,ilzccum=23
  integer::nprof=23
  logical::logscale=.false.

  call read_params

  ! Initialize random number generator
  call rans(1,iseed,allseed)
  localseed=allseed(1,1:IRandNumSize)

  !-----------------------------------------------
  ! Lecture du fichier info du format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
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
  close(10)

  !-----------------------
  ! Profile parameters
  !-----------------------
  xmin=MAX(xcen-rmax,0.0d0)
  xmax=MIN(xcen+rmax,1.0d0)
  ymin=MAX(ycen-rmax,0.0d0)
  ymax=MIN(ycen+rmax,1.0d0)
  zmin=MAX(zcen-rmax,0.0d0)
  zmax=MIN(zcen+rmax,1.0d0)

  write(*,*)'Working array =',nrad
  allocate(r(1:nrad))
  do i=1,nrad
     r(i)=dble(i)*rmax/dble(nrad)
  end do
  allocate(prof(1:nrad,1:nprof))
  prof=0.0d0

  write(*,*)'Bounding box is:'
  write(*,*)xmin,xmax
  write(*,*)ymin,ymax
  write(*,*)zmin,zmax

  ncell=4000000
  allocate(x(1:ncell),y(1:ncell),z(1:ncell))
  allocate(l(1:ncell),var(1:ncell,1:6))
  x=0D0; y=0D0; z=0D0; l=0; var=0D0

  write(*,*)'Generating random sampling points'
  write(*,*)'ncell=',ncell
  icell=0
  if(logscale)then
     do while (icell<ncell)
        logrmax=log10(rmax)
        logrmin=log10(rmin)
        call ranf(localseed,xx)
        call ranf(localseed,yy)
        call ranf(localseed,zz)
        logr=logrmin+(logrmax-logrmin)*xx
        costheta=2.*(yy-0.5)
        sintheta=sqrt(1.-costheta**2)
        cosphi=cos(zz*2.*!DPI)
        sinphi=sin(zz*2.*!DPI)
        rr=10.**logr
        icell=icell+1
        x(icell)=xcen+r*cosphi*sintheta
        y(icell)=ycen+r*sinphi*sintheta
        z(icell)=zcen+r*costheta
     end do
  else
     do while (icell<ncell)
        call ranf(localseed,xx)
        call ranf(localseed,yy)
        call ranf(localseed,zz)
        rr=(xx-0.5)**2+(yy-0.5)**2+(zz-0.5)**2
        if(rr<0.25)then
           icell=icell+1
           x(icell)=xcen+(2.*xx-1.)*rmax
           y(icell)=ycen+(2.*yy-1.)*rmax
           z(icell)=zcen+(2.*zz-1.)*rmax
        end if
     end do
  endif

  ! Sampling points volume element
  dv=4./3.*3.1415926*(rmax*boxlen)**3./dble(ncell)

  do ilevel=1,lmax
     dx=boxlen*0.5**ilevel
     if(dx**3<dv)exit
  end do
  write(*,*)'Using max level=',ilevel

  call getcell(x,y,z,var,l,ncell,6,repository,levelmax=ilevel)
  
  do i=1,ncell
     rad2=(x(i)-xcen)**2+(y(i)-ycen)**2+(z(i)-zcen)**2
     irad=int(dble(nrad)*sqrt(rad2)/rmax)+1
     xx=x(i)-xcen
     yy=y(i)-ycen
     zz=z(i)-zcen
     uu=var(i,2)-ucen/(unit_v/1e5)
     vv=var(i,3)-vcen/(unit_v/1e5)
     ww=var(i,4)-wcen/(unit_v/1e5)
     prof(irad,id)=prof(irad,id)+var(i,1)*dv
     prof(irad,iu)=prof(irad,iu)+var(i,1)*dv*var(i,2)
     prof(irad,iv)=prof(irad,iv)+var(i,1)*dv*var(i,3)
     prof(irad,iw)=prof(irad,iw)+var(i,1)*dv*var(i,4)
     prof(irad,ilx)=prof(irad,ilx)+var(i,1)*dv*(yy*ww-zz*vv)
     prof(irad,ily)=prof(irad,ily)-var(i,1)*dv*(xx*ww-zz*uu)
     prof(irad,ilz)=prof(irad,ilz)+var(i,1)*dv*(xx*vv-yy*uu)
     prof(irad,ip)=prof(irad,ip)+var(i,5)*dv
     prof(irad,imet)=prof(irad,imet)+var(i,1)*var(i,6)*dv
     tt=var(i,5)/var(i,1)*unit_v**2*1.66d-24/1.38d-16
     nH=var(i,1)*unit_d/1.66d-24*0.76
     if(tt>1d5.and.nH<0.1)then
        prof(irad,ihot)=prof(irad,ihot)+var(i,1)*dv
     else
        prof(irad,ilxc)=prof(irad,ilxc)+var(i,1)*dv*(yy*ww-zz*vv)
        prof(irad,ilyc)=prof(irad,ilyc)-var(i,1)*dv*(xx*ww-zz*uu)
        prof(irad,ilzc)=prof(irad,ilzc)+var(i,1)*dv*(xx*vv-yy*uu)        
     endif
  end do
  
  ! Compute cumulated profiles
  mcum=0d0
  ucum=0d0; vcum=0d0; wcum=0d0
  lxcum=0d0; lycum=0d0; lzcum=0d0
  lxccum=0d0; lyccum=0d0; lzccum=0d0
  hotcum=0d0
  do irad=1,nrad
     mcum=mcum+prof(irad,id)
     ucum=ucum+prof(irad,iu)
     vcum=vcum+prof(irad,iv)
     wcum=wcum+prof(irad,iw)
     lxcum=lxcum+prof(irad,ilx)
     lycum=lycum+prof(irad,ily)
     lzcum=lzcum+prof(irad,ilz)
     lxccum=lxccum+prof(irad,ilxc)
     lyccum=lyccum+prof(irad,ilyc)
     lzccum=lzccum+prof(irad,ilzc)
     hotcum=hotcum+prof(irad,ihot)

     prof(irad,imcum)=mcum
     prof(irad,iucum)=ucum
     prof(irad,ivcum)=vcum
     prof(irad,iwcum)=wcum
     prof(irad,ilxcum)=lxcum
     prof(irad,ilycum)=lycum
     prof(irad,ilzcum)=lzcum
     prof(irad,ilxccum)=lxccum
     prof(irad,ilyccum)=lyccum
     prof(irad,ilzccum)=lzccum
     prof(irad,ihot)=hotcum
  end do

  ! Convert profiles into proper astro units
  rprev=0d0
  do irad=1,nrad
     r(irad)=r(irad)*boxlen
     vol=4./3.*3.1415926*(r(irad)**3-rprev**3)
     if(prof(irad,id)>0.0)then
        prof(irad,imet)=(prof(irad,imet)/prof(irad,id))/0.02
        prof(irad,ip)=sqrt(prof(irad,ip)/prof(irad,id))*unit_v/1d5
        prof(irad,iu)=prof(irad,iu)/prof(irad,id)*unit_v/1d5
        prof(irad,iv)=prof(irad,iv)/prof(irad,id)*unit_v/1d5
        prof(irad,iw)=prof(irad,iw)/prof(irad,id)*unit_v/1d5
        prof(irad,ilx)=prof(irad,ilx)/prof(irad,id)*unit_v/1d5/r(irad)
        prof(irad,ily)=prof(irad,ily)/prof(irad,id)*unit_v/1d5/r(irad)
        prof(irad,ilz)=prof(irad,ilz)/prof(irad,id)*unit_v/1d5/r(irad)
     endif
     prof(irad,id)=prof(irad,id)/vol*unit_d/1.66d-24
     if(prof(irad,imcum)>0.0)then
        prof(irad,iucum)=prof(irad,iucum)/prof(irad,imcum)*unit_v/1d5
        prof(irad,ivcum)=prof(irad,ivcum)/prof(irad,imcum)*unit_v/1d5
        prof(irad,iwcum)=prof(irad,iwcum)/prof(irad,imcum)*unit_v/1d5
        prof(irad,ilxcum)=prof(irad,ilxcum)/prof(irad,imcum)*unit_v/1d5/r(irad)
        prof(irad,ilycum)=prof(irad,ilycum)/prof(irad,imcum)*unit_v/1d5/r(irad)
        prof(irad,ilzcum)=prof(irad,ilzcum)/prof(irad,imcum)*unit_v/1d5/r(irad)
     endif
     mcold=prof(irad,imcum)-prof(irad,ihot)
     if(mcold>0.0)then
        prof(irad,ilxccum)=prof(irad,ilxccum)/mcold*unit_v/1d5/r(irad)
        prof(irad,ilyccum)=prof(irad,ilyccum)/mcold*unit_v/1d5/r(irad)
        prof(irad,ilzccum)=prof(irad,ilzccum)/mcold*unit_v/1d5/r(irad)
     endif
     prof(irad,imcum)=sqrt(6.67e-8*prof(irad,imcum)*unit_m/r(irad)/unit_l)/1d5
     rprev=r(irad)
  end do

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  open(unit=10,file=TRIM(nomfich)//".gas",form='formatted')
  write(10,'(A150)')" r(kpc)      n_g(H/cc)   vc_g(H/cc)  mc_g(Msol)  cu_g(km/s)  cv_g(km/s)  cw_g(km/s)  cl_g(km/s)  lx_g        ly_g        lz_g        c_g(km/s)              "
  do irad=1,nrad
     write(10,999)r(irad)*unit_l/3.08d21,prof(irad,id),prof(irad,imcum),prof(irad,imcum)**2*r(irad)*unit_l/6.67e-8*1d10/2d33 &
          & ,prof(irad,iucum),prof(irad,ivcum),prof(irad,iwcum),sqrt(prof(irad,ilxcum)**2+prof(irad,ilycum)**2+prof(irad,ilzcum)**2) &
          & ,prof(irad,ilxcum)/sqrt(prof(irad,ilxcum)**2+prof(irad,ilycum)**2+prof(irad,ilzcum)**2+1e-30) &
          & ,prof(irad,ilycum)/sqrt(prof(irad,ilxcum)**2+prof(irad,ilycum)**2+prof(irad,ilzcum)**2+1e-30) &
          & ,prof(irad,ilzcum)/sqrt(prof(irad,ilxcum)**2+prof(irad,ilycum)**2+prof(irad,ilzcum)**2+1e-30),prof(irad,ip) &
          & ,prof(irad,imet),prof(irad,ihot)*unit_m/2d33,sqrt(prof(irad,ilxccum)**2+prof(irad,ilyccum)**2+prof(irad,ilzccum)**2) &
          & ,prof(irad,ilxccum)/sqrt(prof(irad,ilxccum)**2+prof(irad,ilyccum)**2+prof(irad,ilzccum)**2+1e-30) &
          & ,prof(irad,ilyccum)/sqrt(prof(irad,ilxccum)**2+prof(irad,ilyccum)**2+prof(irad,ilzccum)**2+1e-30) &
          & ,prof(irad,ilzccum)/sqrt(prof(irad,ilxccum)**2+prof(irad,ilyccum)**2+prof(irad,ilzccum)**2+1e-30)
  end do
  close(10)

999 format(50(1PE10.3,2X))

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
       print *, 'usage: amr2prof -inp  input_dir'
       print *, '                -out  output_file'
       print *, '                 [-xce xcen] '
       print *, '                 [-yce ycen] '
       print *, '                 [-zce zcen] '
       print *, '                 [-uce ucen] '
       print *, '                 [-vce vcen] '
       print *, '                 [-wce wcen] '
       print *, '                 [-rmi rmin] '
       print *, '                 [-rma rmax] '
       print *, '                 [-nra nrad] '
       print *, '                 [-lma lmax] '
       print *, 'ex: amr2prof -inp output_00001 -out prof.dat'// &
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
       case ('-rmi')
          read (arg,*) rmin
       case ('-rma')
          read (arg,*) rmax
       case ('-lma')
          read (arg,*) lmax
       case ('-log')
          read (arg,*) logscale
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params
  
end program amr2prof
