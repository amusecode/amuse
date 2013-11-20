program histo_main
  !--------------------------------------------------------------------------
  ! Ce programme calcule un histogramme 2D dans le plan rho-T
  ! d'une simulation RAMSES.
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  integer::ndim,n,i,j,k,twotondim,ncoarse,type=0,domax=0
  integer::ivar,nvar,ncpu,ncpuh,lmax=0,nboundary,ngrid_current
  integer::nx,ny,nz,ilevel,idim,jdim,kdim,icell
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer::ngridmax,nstep_coarse,icpu,ncpu_read
  integer::nhx=0,nhy=0,ihx,ihy,ivar1,ivar2
  real::gamma,smallr,smallc,gammah
  real::boxlen,boxlen2
  real::t,aexp,hexp,t2,aexp2,hexp2
  real::omega_m,omega_l,omega_k,omega_b
  real::scale_l,scale_d,scale_t
  real::omega_m2,omega_l2,omega_k2,omega_b2
  real::tpoly=0d0, npoly=1.0
  
  integer::nx_sample=0,ny_sample=0,nz_sample=0
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax,gcc=0
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin,levelmin
  integer::ix,iy,iz,ixp1,iyp1,izp1,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::levelmax
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dxmax,dummy
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::dymin,dymax,tymin,tymax
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,xx,yy,zz,dxx,dyy,dd,dt
  real(KIND=8),dimension(:,:),allocatable::x,xg,histo,im
  real(KIND=8)::dmin=0.0,dmax=0.0,tmin=0.0,tmax=0.0
  real(KIND=8)::h0,unit_l,unit_d,unit_t,total_mass,mmm
  real(KIND=8),dimension(:,:,:),allocatable::var
  real(KIND=4),dimension(:,:,:),allocatable::toto
  real(KIND=8),dimension(:)  ,allocatable::rho,pre
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:)  ,allocatable::isp
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  logical::ok,ok_part,ok_cell
  real(KIND=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  character(LEN=1)::proj='z'

  type level
     integer::ilevel
     integer::ngrid
     real(KIND=4),dimension(:,:,:),pointer::cube
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level

  type(level),dimension(1:100)::grid

  call read_params


  ! Histogram parameters
  if(nhx==0)then
     nhx=200
  endif
  if(nhy==0)then
     nhy=200
  endif
  if(dmin==0.0)then
     dmin=1d-8
  end if
  if(dmax==0.0)then
     dmax=1d+1
  end if
  if(tmin==0.0)then
     tmin=1d-4
  end if
  if(tmax==0.0)then
     tmax=1d+9
  end if
  dymin=log10(dmin)
  dymax=log10(dmax)
  tymin=log10(tmin)
  tymax=log10(tmax)

  write(*,*)nhx,dymin,dymax
  write(*,*)nhy,tymin,tymax

  allocate(histo(nhx,nhy))
  histo=0.0d0

  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
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
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  if(ndim==2)then
     write(*,*)'Output file contains 2D data'
     write(*,*)'Aborting'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,'("ncpu        =",I11)')ncpu
  read(10,'("ndim        =",I11)')ndim
  read(10,'("levelmin    =",I11)')levelmin
  read(10,'("levelmax    =",I11)')levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,*)
  read(10,'("time        =",E23.15)')t
  read(10,'("aexp        =",E23.15)')aexp
  read(10,'("H0          =",E23.15)')h0
  read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'("unit_t      =",E23.15)')unit_t
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

  !-----------------------
  ! Map parameters
  !-----------------------
  if(lmax==0)then
     lmax=nlevelmax
  endif
  write(*,*)'time=',t
  write(*,*)'Working resolution =',2**lmax
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

  if(TRIM(ordering).eq.'hilbert')then

     dxmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
     do ilevel=1,lmax
        dx=0.5d0**ilevel
        if(dx.lt.dxmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xxmin*dble(maxdom))
        imax=imin+1
        jmin=int(yymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zzmin*dble(maxdom))
        kmax=kmin+1
     endif

     dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
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
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do
     
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
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

  !-----------------------------------------------
  ! Compute histogram
  !----------------------------------------------

  total_mass=0d0
  ! Loop over processor files
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=10,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     do i=1,21
        read(10)
     end do
     ! Read grid numbers
     read(10)ngridlevel
     ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
     read(10)
     if(nboundary>0)then
        do i=1,2
           read(10)
        end do
        read(10)ngridbound
        ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
     endif
     read(10)
! ROM: comment the single follwing line for old stuff
     read(10)
     if(TRIM(ordering).eq.'bisection')then
        do i=1,5
           read(10)
        end do
     else
        read(10)
     endif
     read(10)
     read(10)
     read(10)

     ! Open HYDRO file and skip header
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh
     read(11)
     read(11)
     read(11)
     read(11)

     ! Loop over levels
     do ilevel=1,lmax

        ! Geometry
        dx=0.5**ilevel
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(xg (1:ngrida,1:ndim))
           allocate(son(1:ngrida,1:twotondim))
           allocate(var(1:ngrida,1:twotondim,1:nvarh))
           allocate(x  (1:ngrida,1:ndim))
           allocate(rho(1:ngrida))
           allocate(pre(1:ngrida))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

           ! Read AMR data
           if(ngridfile(j,ilevel)>0)then
              read(10) ! Skip grid index
              read(10) ! Skip next index
              read(10) ! Skip prev index
              ! Read grid center
              do idim=1,ndim
                 if(j.eq.icpu)then
                    read(10)xg(:,idim)
                 else
                    read(10)
                 endif
              end do
              read(10) ! Skip father index
              do ind=1,2*ndim
                 read(10) ! Skip nbor index
              end do
              ! Read son index
              do ind=1,twotondim
                 if(j.eq.icpu)then
                    read(10)son(:,ind)
                 else
                    read(10)
                 end if
              end do
              ! Skip cpu map
              do ind=1,twotondim
                 read(10)
              end do
              ! Skip refinement map
              do ind=1,twotondim
                 read(10)
              end do
           endif

           ! Read HYDRO data
           read(11)
           read(11)
           if(ngridfile(j,ilevel)>0)then
              ! Read hydro variables
              do ind=1,twotondim
                 do ivar=1,nvarh
                    if(j.eq.icpu)then
                       read(11)var(:,ind,ivar)
                    else
                       read(11)
                    end if
                 end do
              end do
           end if
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

              ! Compute cell center
              do i=1,ngrida
                 x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                 x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                 x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
              end do
              ! Check if cell is refined
              do i=1,ngrida
                 ref(i)=son(i,ind)>0.and.ilevel<lmax
              end do
              ! Extract variable
              rho = var(:,ind,1)*unit_d
              pre = var(:,ind,5)*unit_d*(unit_l/unit_t)**2
              ! Store data cube
              do i=1,ngrida
                 ok_cell= .not.ref(i).and. &
                      & x(i,1)>xxmin.and.x(i,1)<xxmax.and. &
                      & x(i,2)>yymin.and.x(i,2)<yymax.and. &
                      & x(i,3)>zzmin.and.x(i,3)<zzmax
                 if(ok_cell)then
                    if (gcc==0)then
                       xx=log10(rho(i)/1.66d-24*0.76)
                    else
                       xx=log10(rho(i))
                    end if
                    yy=log10(pre(i)/rho(i)/1.38d-16*1.66d-24 &
                         & -tpoly*(rho(i)/1.66d-24*0.76/npoly))
                    dxx=(xx-dymin)/(dymax-dymin)*dble(nhx)
                    dyy=(yy-tymin)/(tymax-tymin)*dble(nhy)
                    ihx=int(dxx)
                    ihy=int(dyy)
                    ihx=MAX(ihx,0)
                    ihx=MIN(ihx,nhx-1)
                    ihy=MAX(ihy,0)
                    ihy=MIN(ihy,nhy-1)
                    mmm=rho(i)*(dx*unit_l)**3/2d33
                    total_mass=total_mass+mmm
                    histo(ihx+1,ihy+1)=histo(ihx+1,ihy+1)+mmm
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,var,ref,rho,pre,x)
        endif

     end do
     ! End loop over levels

     close(10)
     close(11)

  end do
  ! End loop over cpus

  write(*,*)'Total mass=',total_mass

  ! Output file
  if (filetype=='bin')then 
     nomfich=TRIM(outfich)
     write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
     open(unit=10,file=nomfich,form='unformatted')
     write(10)nhx,nhy
     write(10)real(histo/total_mass)
     write(10)dymin,dymax
     write(10)tymin,tymax
     close(10)
  end if

  if (filetype=='ascii')then 
     nomfich=TRIM(outfich)
     write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
     open(unit=10,file=nomfich,form='formatted')
     dd=(dymax-dymin)/nhx
     dt=(tymax-tymin)/nhy
     do i=1,nhx
        do j=1,nhy
           write(10,*)dymin+i*dd,tymin+j*dt,histo(i,j)/total_mass
        end do
        write(10,*) " "
     end do
     close(10)
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
       print *, 'usage: histo -inp  input_dir'
       print *, '             -out  output_file'
       print *, '            [-dmi rho_min] '
       print *, '            [-dma rho_max] '
       print *, '            [-tmi T_min  ] '
       print *, '            [-tma T_max  ] '
       print *, '            [-nx  nx_rho ] '
       print *, '            [-ny  ny_T   ] '
       print *, '            [-fil ascii or bin  ] '
       print *, 'ex: histo -inp output_00001 -out map.dat'
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
       case ('-dmi')
          read (arg,*) dmin
       case ('-dma')
          read (arg,*) dmax
       case ('-tmi')
          read (arg,*) tmin
       case ('-tma')
          read (arg,*) tmax
       case ('-nx')
          read (arg,*) nhx
       case ('-ny')
          read (arg,*) nhy
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
       case ('-tpo')
          read (arg,*) tpoly
       case ('-npo')
          read (arg,*) npoly
       case ('-fil')
          read (arg,*) filetype
       case ('-gcc') ! if set to one, density output in g/cc
          read (arg,*) gcc
       case default
          print '("unknown option ",a2," ignored")', opt

          

       end select
    end do
    
  end subroutine read_params
  
end program histo_main

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
