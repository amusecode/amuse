program amr2cut
  use io_ramses

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  implicit none
  integer::ndim,ncell,n,i,j,jm,k,kk,kkk,twotondim,ncoarse,type=0,domax=0
  integer::ivar,nvar,ncpu,ncpuh,lmax=120,nboundary,ngrid_current
  integer::nx=0,ny=0,nz=0,ilevel,idim,jdim,kdim,icell
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer::ngridmax,nstep_coarse,icpu,ncpu_read
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  real::gamma,smallr,smallc,gammah
  real::boxlen,boxlen2
  real::t,aexp,hexp,t2,aexp2,hexp2
  real::omega_m,omega_l,omega_k,omega_b
  real::scale_l,scale_d,scale_t
  real::omega_m2,omega_l2,omega_k2,omega_b2
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin,levelmin
  integer::ix,iy,iz,ixp1,iyp1,izp1,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::xcenter,ycenter,zcenter
  real(KIND=8)::jxin=0,jyin=0,jzin=0,jx,jy,jz
  real(KIND=8)::kxin,kyin,kzin,kx,ky,kz
  real(KIND=8)::lxin,lyin,lzin,lx,ly,lz
  real(KIND=8)::dkey,order_min,dmax,dummy
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,dz,dm
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,xx,yy,zz
  real(KIND=8)::dx_loc,dx_amr
  real(KIND=8),dimension(:),allocatable::x,y,z
  real(KIND=8),dimension(:,:),allocatable::var
  integer,dimension(:),allocatable::l
  integer::nskip
  real(kind=4),dimension(:,:),allocatable::map
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  logical::ok,ok_part,ok_cell
  logical::rotation=.false.,sideon=.false.

  call read_params

  !-----------------------
  ! Rotation parameters 
  !-----------------------
  if(abs(jxin)+abs(jyin)+abs(jzin)>0)then
     rotation=.true.
     write(*,*)'Performing rotation'
     write(*,*)jxin,jyin,jzin
  else
     write(*,*)'Invalid rotation axis'
     stop
  endif

  kxin=0.
  kyin=-jzin
  kzin=jyin

  lxin=jyin*kzin-jzin*kyin
  lyin=jzin*kxin-jxin*kzin
  lzin=jxin*kyin-jyin*kxin

  jx=jxin/sqrt(jxin**2+jyin**2+jzin**2)
  jy=jyin/sqrt(jxin**2+jyin**2+jzin**2)
  jz=jzin/sqrt(jxin**2+jyin**2+jzin**2)

  kx=kxin/sqrt(kxin**2+kyin**2+kzin**2)
  ky=kyin/sqrt(kxin**2+kyin**2+kzin**2)
  kz=kzin/sqrt(kxin**2+kyin**2+kzin**2)

  lx=lxin/sqrt(lxin**2+lyin**2+lzin**2)
  ly=lyin/sqrt(lxin**2+lyin**2+lzin**2)
  lz=lzin/sqrt(lxin**2+lyin**2+lzin**2)

  xcenter=0.5*(xmin+xmax)
  ycenter=0.5*(ymin+ymax)
  zcenter=0.5*(zmin+zmax)

  if(sideon)then
     allocate(map(1:nx,1:nz))
  else
     allocate(map(1:nx,1:ny))
  endif

  map=0.0d0
  dx=(xmax-xmin)/dble(nx)
  dy=(ymax-ymin)/dble(ny)
  dz=(zmax-zmin)/dble(nz)

  dx_loc=min(dx,dy,dx)
  
  dx_amr=1.0
  lmax=0
  do while (dx_amr.gt.dx_loc)
     lmax=lmax+1
     dx_amr=dx_amr/2.0
!     write(*,*)lmax,dx_amr,dx_loc
  end do
  lmax=lmax-2
  dx_amr=dx_amr*4.0
  write(*,*)'levmax=',lmax
  write(*,*)'dx_amr=',dx_amr
  write(*,*)'dx_usr=',dx_loc

  write(*,*)'Bounding box is:'
  write(*,*)xmin,xmax
  write(*,*)ymin,ymax
  write(*,*)zmin,zmax

  nskip=16

  do kk=1,nz,nskip

     ncell=nx*ny*min(nz-kk+1,nskip)
     nvar=5
     if(.not. allocated(x))then
        write(*,*)'Work with ncell=',ncell
        allocate(x(1:ncell),y(1:ncell),z(1:ncell))
        allocate(l(1:ncell),var(1:ncell,1:nvar))
     endif
     x=0D0
     y=0D0
     z=0D0
     l=0
     var=0D0

     do kkk=1,min(nz-kk+1,nskip)
        do j=1,ny
           do i=1,nx
              ind=1+(i-1)+(j-1)*nx+(kkk-1)*nx*ny
              k=kk+kkk-1
              x(ind)=xcenter+(dble(i-0.5)/dble(nx)-0.5)*(xmax-xmin)*lx+(dble(j-0.5)/dble(ny)-0.5)*(ymax-ymin)*kx+(dble(k-0.5)/dble(nz)-0.5)*(zmax-zmin)*jx
              y(ind)=ycenter+(dble(i-0.5)/dble(nx)-0.5)*(xmax-xmin)*ly+(dble(j-0.5)/dble(ny)-0.5)*(ymax-ymin)*ky+(dble(k-0.5)/dble(nz)-0.5)*(zmax-zmin)*jy
              z(ind)=zcenter+(dble(i-0.5)/dble(nx)-0.5)*(xmax-xmin)*lz+(dble(j-0.5)/dble(ny)-0.5)*(ymax-ymin)*kz+(dble(k-0.5)/dble(nz)-0.5)*(zmax-zmin)*jz
           end do
        end do
     end do

     write(*,*)kk,'/',nz,min(nz-kk+1,nskip)
     
     call getcell(x,y,z,var,l,ncell,nvar,repository,levelmax=lmax)

     do kkk=1,min(nz-kk+1,nskip)
        do j=1,ny
           do i=1,nx
              ind=1+(i-1)+(j-1)*nx+(kkk-1)*nx*ny
              k=kk+kkk-1
              if(sideon)then
                 jm=k
                 dm=dy
              else
                 jm=j
                 dm=dz
              endif
              select case (type)
              case(0)
                 map(i,jm)=max(map(i,jm),real(l(ind)))
              case(1)
                 map(i,jm)=map(i,jm)+var(ind,1)*dm
              case(5)
                 map(i,jm)=map(i,jm)+var(ind,5)*dm
              end select
           end do
        end do
     end do
  end do
  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)

  if (filetype=='bin')then
     open(unit=20,file=nomfich,form='unformatted')
     if(sideon)then
        write(20)nx,nz
        write(20)map
        write(20)xmin,xmax
        write(20)zmin,zmax
     else
        write(20)nx,ny
        write(20)map
        write(20)xmin,xmax
        write(20)ymin,ymax
     endif
     close(20)
  endif
  if (filetype=='ascii')then
     open(unit=20,file=nomfich,form='formatted')
     if(sideon)then
        do j=1,ny
           do i=1,nx
              xx=xmin+(dble(i)-0.5)/dble(nx)*(xmax-xmin)
              yy=ymin+(dble(j)-0.5)/dble(ny)*(ymax-ymin)
              write(20,'(20(1PE15.7,2X))')xx,yy,map(i,j)
           end do
           write(20,*) " "
        end do
     else
        do j=1,nz
           do i=1,nx
              xx=xmin+(dble(i)-0.5)/dble(nx)*(xmax-xmin)
              yy=zmin+(dble(j)-0.5)/dble(nz)*(zmax-zmin)
              write(20,'(20(1PE15.7,2X))')xx,yy,map(i,j)
           end do
           write(20,*) " "
        end do
     endif
     close(20)
  endif
  
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
       print *, 'usage: amr2cut  -inp  input_dir'
       print *, '                -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-nx  nx] '
       print *, '                 [-ny  ny] '
       print *, '                 [-nz  nz] '
       print *, '                 [-fil filetype] '
       print *, 'ex: amr2cut -inp output_00001 -out cube.dat'// &
            &   ' -xmi 0.1 -xma 0.7 -lma 12'
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
       case ('-fil')
          filetype = trim(arg)
       case ('-typ')
          read (arg,*) type
       case('-jx')
          read (arg,*) jxin
       case('-jy')
          read (arg,*) jyin
       case('-jz')
          read (arg,*) jzin
         case ('-sid')
            read (arg,*) sideon
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
       case ('-lma')
          read (arg,*) lmax
       case ('-nx')
          read (arg,*) nx
       case ('-ny')
          read (arg,*) ny
       case ('-nz')
          read (arg,*) nz
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params
  
end program amr2cut
