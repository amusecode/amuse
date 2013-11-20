program amr2tipsy
  use io_ramses

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  implicit none
  integer::ndim,ncell,n,i,j,k,twotondim,ncoarse,type=0,domax=0
  integer::ivar,nvar,ncpu,ncpuh,lmax=0,nboundary,ngrid_current
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
  real(KIND=8)::dkey,order_min,dmax,dummy
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,xx,yy,zz
  real(KIND=8),dimension(:),allocatable::x,y,z
  real(KIND=8),dimension(:,:),allocatable::var
  integer,dimension(:),allocatable::l
  real(kind=4),dimension(:,:,:),allocatable::toto
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  logical::ok,ok_part,ok_cell

  call read_params

  ncell=nx*ny*nz
  write(*,*)'Work with ncell=',ncell
  allocate(x(1:ncell),y(1:ncell),z(1:ncell),l(1:ncell),var(1:ncell,1:4))
  x=0D0
  y=0D0
  z=0D0
  l=0
  var=0D0

  write(*,*)'Bounding box is:'
  write(*,*)xmin,xmax
  write(*,*)ymin,ymax
  write(*,*)zmin,zmax
  do i=1,nx
  do j=1,ny
  do k=1,nz
     ind=1+(i-1)+(j-1)*nx+(k-1)*nx*ny
     x(ind)=xmin+dble(i-0.5)/dble(nx)*(xmax-xmin)
     y(ind)=ymin+dble(j-0.5)/dble(ny)*(ymax-ymin)
     z(ind)=zmin+dble(k-0.5)/dble(nz)*(zmax-zmin)
  end do
  end do
  end do

  call getcell(x,y,z,var,l,ncell,4,repository,levelmax=lmax)

  allocate(toto(1:nx,1:ny,1:nz))
  do i=1,nx
  do j=1,ny
  do k=1,nz
     ind=1+(i-1)+(j-1)*nx+(k-1)*nx*ny
     toto(i,j,k)=var(ind,type)
  end do
  end do
  end do

  ! Output file
  if(TRIM(filetype).eq.'bin')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='unformatted')
     write(20)nx,ny,nz
     write(20)toto
     close(20)
  endif
  
  if(TRIM(filetype).eq.'grafic')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='unformatted')
     dummy=0.0
     write(20)nx,ny,nz,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
     do iz=1,nz
        write(20)((toto(ix,iy,iz),ix=1,nx),iy=1,ny)
     end do
     close(20)
  endif
  
  if(TRIM(filetype).eq.'vtk')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='formatted')
     write(20,'("# vtk DataFile Version 2.0")')
     write(20,'("RAMSES data using vtk file format")')
     write(20,'("ASCII")')
     write(20,'("DATASET STRUCTURED_POINTS")')
     write(20,'("DIMENSIONS ",3(I3,1x))')nx,ny,nz
     write(20,'("ORIGIN 0.0 0.0 0.0")')
     write(20,'("SPACINGS 1.0 1.0 1.0")')
     write(20,'(" ")')
     write(20,'("POINT_DATA ",I8)')nx*ny*nz
     write(20,'("SCALARS values float")')
     write(20,'("LOOKUP_TABLE default")')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              write(20,'(F9.3)')toto(i,j,k)
           end do
        end do
     end do
     close(20)
  endif
  
  if(TRIM(filetype).eq.'tipsy')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='formatted')
     write(20,*)ncell,0,0
     dummy=0.0
     write(20,*)3
     write(20,*)real(dummy,kind=4)
     do i=1,ncell
        write(20,*)var(i,1)
     end do
     do i=1,ncell
        write(20,*)x(i)
     end do
     do i=1,ncell
        write(20,*)y(i)
     end do
     do i=1,ncell
        write(20,*)z(i)
     end do
     do i=1,ncell
        write(20,*)var(i,2)
     end do
     do i=1,ncell
        write(20,*)var(i,3)
     end do
     do i=1,ncell
        write(20,*)var(i,4)
     end do
     do i=1,ncell
        write(20,*)dummy
     end do
     do i=1,ncell
        write(20,*)dummy
     end do
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
       print *, 'usage: amr2tipsy -inp  input_dir'
       print *, '                 -out  output_file'
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
       print *, 'ex: amr2tipsy -inp output_00001 -out cube.dat'// &
            &   ' -typ 1 -xmi 0.1 -xma 0.7 -lma 12'
       print *, ' '
       print *, ' type :-1 = cpu number'
       print *, ' type : 0 = ref. level (default)'
       print *, ' type : 1-9 = variable number'
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
       case ('-typ')
          read (arg,*) type
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params
  
end program amr2tipsy
