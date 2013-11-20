program amr2cylprof
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
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax,nrad=100,nprof=15
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin,levelmin,levelmax
  integer::ix,iy,iz,ixp1,iyp1,izp1,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,dummy,h0,dv
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,rmax=0.5,hmax=0.0,rrmax
  real(KIND=8)::xcen=0.5,ycen=0.5,zcen=0.5
  real(KIND=8)::ucen=0.0,vcen=0.0,wcen=0.0
  real(KIND=8)::xx,yy,zz,uu,vv,ww,rr
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez
  real(KIND=8)::unit_l,unit_t,unit_d,unit_m,unit_v
  real(KIND=8)::rad2,vol,surf,rprev
  real(KIND=8)::mcum,ucum,vcum,wcum,lxcum,lycum,lzcum
  real(kind=8)::jxin=0.0,jyin=0.0,jzin=1.0,jx,jy,jz,jt
  real(kind=8)::rx,ry,rz,tx,ty,tz,r_cyl,z_coord,u_r,u_t,u_z
  real(kind=8)::vdotu,xxx,yyy,zzz,xpara,ypara,zpara,xperp,yperp,zperp,xcros,ycros,zcros
  real(kind=8)::cos_alpha,sin_alpha
  real(KIND=8),dimension(:),allocatable::x,y,z,r
  real(KIND=8),dimension(:,:),allocatable::var,prof
  integer,dimension(:),allocatable::l
  real(kind=4),dimension(:,:,:),allocatable::toto
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  character(LEN=5)::nchar
  logical::ok,ok_part,ok_cell,cosmo
  integer::id=1,iu=2,iv=3,iw=4,iu2=5,iv2=6,iw2=7,ip=8

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
  if(hmax==0.0)hmax=rmax
  rrmax=sqrt(rmax**2+hmax**2)
  xmin=MAX(xcen-rrmax,0.0d0)
  xmax=MIN(xcen+rrmax,1.0d0)
  ymin=MAX(ycen-rrmax,0.0d0)
  ymax=MIN(ycen+rrmax,1.0d0)
  zmin=MAX(zcen-rrmax,0.0d0)
  zmax=MIN(zcen+rrmax,1.0d0)

  ! Normalized angular momentum vector                                          
  jx=jxin/sqrt(jxin**2+jyin**2+jzin**2)
  jy=jyin/sqrt(jxin**2+jyin**2+jzin**2)
  jz=jzin/sqrt(jxin**2+jyin**2+jzin**2)

  ! Normalized rotation axis
  jt=jx**2+jy**2
  if(jt>0)then
     rx=-jy/sqrt(jt)
     ry=jx/sqrt(jt)
  else
     rx=1d0
     ry=0d0
  endif
  rz=0d0


  ! Rotation angle
  cos_alpha=jz
  sin_alpha=sqrt(1.0-cos_alpha**2)

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

  ncell=2000000
  allocate(x(1:ncell),y(1:ncell),z(1:ncell))
  allocate(l(1:ncell),var(1:ncell,1:6))
  x=0D0; y=0D0; z=0D0; l=0; var=0D0

  write(*,*)'Generating random sampling points'
  icell=0
  do while (icell<ncell)
     call ranf(localseed,xx)
     call ranf(localseed,yy)
     call ranf(localseed,zz)
     rr=(xx-0.5)**2+(yy-0.5)**2
     if(rr<0.25)then
        icell=icell+1
        xx=(2.*xx-1.)*rmax
        yy=(2.*yy-1.)*rmax
        zz=(2.*zz-1.)*hmax

        vdotu=xx*rx+yy*ry+zz*rz
        xpara=vdotu*rx
        ypara=vdotu*ry
        zpara=vdotu*rz
        xperp=xx-xpara
        yperp=yy-ypara
        zperp=zz-zpara
        xcros=ry*zperp-rz*yperp
        ycros=rz*xperp-rx*zperp
        zcros=rx*yperp-ry*xperp
        xxx=xperp*cos_alpha+xcros*sin_alpha+xpara
        yyy=yperp*cos_alpha+ycros*sin_alpha+ypara
        zzz=zperp*cos_alpha+zcros*sin_alpha+zpara

        x(icell)=xcen+xxx
        y(icell)=ycen+yyy
        z(icell)=zcen+zzz
     end if
  end do

  call getcell(x,y,z,var,l,ncell,5,repository,levelmax=lmax)
  
  do i=1,ncell
     xx=x(i)-xcen
     yy=y(i)-ycen
     zz=z(i)-zcen
     z_coord=xx*jx+yy*jy+zz*jz
     r_cyl=sqrt((xx-z_coord*jx)**2+(yy-z_coord*jy)**2+(zz-z_coord*jz)**2)
     irad=int(dble(nrad)*r_cyl/rmax)+1
     ! Galilean invariant frame                                           
     uu=var(i,2)-ucen/(unit_v/1d5)
     vv=var(i,3)-vcen/(unit_v/1d5)
     ww=var(i,4)-wcen/(unit_v/1d5)
     ! Normalized radial vector                                           
     rx=(xx-z_coord*jx)/r_cyl
     ry=(yy-z_coord*jy)/r_cyl
     rz=(zz-z_coord*jz)/r_cyl
     ! Normalized tangential vector                                       
     tx=jy*rz-jz*ry
     ty=jz*rx-jx*rz
     tz=jx*ry-jy*rx
     ! Compute velocity components                                        
     u_z=uu*jx+vv*jy+ww*jz
     u_r=uu*rx+vv*ry+ww*rz
     u_t=uu*tx+vv*ty+ww*tz

     prof(irad,id)=prof(irad,id)+var(i,1)
     prof(irad,iu)=prof(irad,iu)+var(i,1)*u_r
     prof(irad,iv)=prof(irad,iv)+var(i,1)*u_t
     prof(irad,iw)=prof(irad,iw)+var(i,1)*u_z
     prof(irad,iu2)=prof(irad,iu2)+var(i,1)*u_r**2
     prof(irad,iv2)=prof(irad,iv2)+var(i,1)*u_t**2
     prof(irad,iw2)=prof(irad,iw2)+var(i,1)*u_z**2
     prof(irad,ip)=prof(irad,ip)+var(i,5)
  end do
  

  ! Sampling points volume element
  dv=3.1415926*(rmax*boxlen)**2.*2.0*hmax*boxlen/dble(ncell)

  ! Convert profiles into proper astro units
  rprev=0d0
  do irad=1,nrad
     r(irad)=r(irad)*boxlen
     surf=3.1415926*(r(irad)**2-rprev**2)
     if(prof(irad,id)>0.0)then
        prof(irad,ip)=sqrt(prof(irad,ip)/prof(irad,id))*unit_v/1d5
        prof(irad,iu)=prof(irad,iu)/prof(irad,id)*unit_v/1d5
        prof(irad,iv)=prof(irad,iv)/prof(irad,id)*unit_v/1d5
        prof(irad,iw)=prof(irad,iw)/prof(irad,id)*unit_v/1d5
        prof(irad,iu2)=sqrt(prof(irad,iu2)/prof(irad,id)*(unit_v/1d5)**2-prof(irad,iu)**2)
        prof(irad,iv2)=sqrt(prof(irad,iv2)/prof(irad,id)*(unit_v/1d5)**2-prof(irad,iv)**2)
        prof(irad,iw2)=sqrt(prof(irad,iw2)/prof(irad,id)*(unit_v/1d5)**2-prof(irad,iw)**2)
     endif
     prof(irad,id)=prof(irad,id)*dv*unit_m/(surf*unit_l**2)/(2d33/3.08d18**2)
     rprev=r(irad)
  end do

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  open(unit=10,file=TRIM(nomfich)//".gas",form='formatted')
  write(10,'(A97)')" r(kpc)      S_g(Mpc2)   u_r(km/s)   u_t(km/s)   u_z(km/s)   s_r(km/s)   s_t(km/s)   s_z(km/s)   c_g(km/s)" 
  do i=1,nrad
     write(10,999)r(i)*unit_l/3.08d21,(prof(i,ivar),ivar=1,8)
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
       print *, 'usage: amr2prof -inp  input_dir'
       print *, '                -out  output_file'
       print *, '                 [-xce xcen] '
       print *, '                 [-yce ycen] '
       print *, '                 [-zce zcen] '
       print *, '                 [-uce ucen] '
       print *, '                 [-vce vcen] '
       print *, '                 [-wce wcen] '
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
       case ('-jx')
          read (arg,*) jxin
       case ('-jy')
          read (arg,*) jyin
       case ('-jz')
          read (arg,*) jzin
       case ('-nra')
          read (arg,*) nrad
       case ('-rma')
          read (arg,*) rmax
       case ('-hma')
          read (arg,*) hmax
       case ('-lma')
          read (arg,*) lmax
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params
  
end program amr2cylprof
