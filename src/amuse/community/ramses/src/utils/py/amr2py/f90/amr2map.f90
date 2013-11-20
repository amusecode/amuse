program amr2map

!! Convert Ramses output on AMR/hydro into a fits
!! Can also convert hydro data to an ascii file and/or compute the pdf
!! Based on Romain Teyssier's amr2map
!! Florent Renaud - 5 Aug 2011

  implicit none
 
  integer::i, j, k, n, impi, icpu, ilevel, iidim, ivar, ind, ipdf
  character(len=5)::nchar,ncharcpu
  character(len=128)::filename, repository, suffix='', outval, outvalunit
  real(KIND=8)::xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1
  character(len=1)::dir
  integer::lmax=0, typ=1, pdfn=500
  logical::maxval=.false., ascii=.false., pdf=.false., makemap=.true., maxrho=.false.
  
  real(kind=8)::pdfmin=1.0D-3, pdfmax=5.0D5, lpdfampli, lpdfmin, lpdfmax
  real(kind=8),dimension(:),allocatable::pdfhist

  integer::ncpu, ndim, nx, ny, nz, nlevelmax, ngridmax, nboundary, ngrid_current
  integer::twotondim, levelmin, bit_length, maxdom, ndom, ncpu_read, nvarh
  integer::ngrida
  integer::idim, jdim, kdim, imin, imax, jmin, jmax, kmin, kmax
  integer::ix, iy, iz
  integer::nx_full, ny_full, nz_full
  integer::nx_sample=0, ny_sample=0
  real(kind=8)::xxmin, xxmax, yymin, yymax, zzmin, zzmax
  real(kind=8)::dkey, order_min, dx, dmax, dxline, boxlen, t
  real::weight
  character(len=80)::ordering
  character(LEN=80)::GMGM
  logical::ok
  
  integer,dimension(1:8)::idom, jdom, kdom, cpu_min, cpu_max
  real(kind=8),dimension(1:3)::xbound
  real(kind=8),dimension(1:8,1:3)::xc
  integer,dimension(:,:),allocatable::ngridfile, ngridlevel, ngridbound
  integer,dimension(:),allocatable::cpu_list
  real(kind=8),dimension(:),allocatable::bound_key
  real(kind=8),dimension(1:8)::bounding_min, bounding_max
  logical,dimension(:),allocatable::cpu_read

  real(kind=8),dimension(:,:),allocatable::x, xg
  real(kind=8),dimension(:,:,:),allocatable::var
  real(kind=8),dimension(:),allocatable::rho, map
  logical,dimension(:),allocatable::ref
  integer,dimension(:,:),allocatable::son

  type level
    integer::ilevel
    integer::ngrid
    real(KIND=8),dimension(:,:),pointer::map
    real(KIND=8),dimension(:,:),pointer::rho
    integer::imin
    integer::imax
    integer::jmin
    integer::jmax
    integer::kmin
    integer::kmax
  end type level

  type(level),dimension(1:100)::grid

  integer::status, unit, blocksize, bitpix, naxis
  integer,dimension(2)::naxes
  integer::group,fpixel,nelements
  integer,dimension(300,200)::array
!  character(80)::filename
  logical::simple,extend
  character(2)::num


  real(KIND=4),dimension(:,:),allocatable::tmpmap
  integer::nxmap 
  real,dimension(:,:),allocatable::map2
  
  real(kind=8)::scale_nH,scale_vkms,scale_T2,scale_t,scale_l,scale_d
  
!=======================================================================  

  call read_params

  if(ascii) makemap=.false.
  
  if(pdf) then
    lpdfmin = log10(pdfmin)
    lpdfmax = log10(pdfmax)
    lpdfampli = (lpdfmax-lpdfmin) / pdfn
    allocate(pdfhist(1:pdfn))
    pdfhist = 0.0
    typ = 1
    makemap=.false.
  endif


  select case (typ)
    case (-1)
      outval = 'cpu'
      outvalunit = ''
    case (0)
      outval = 'level'
      outvalunit = ''
    case (2) ! x-velocity
      outval = 'v_x'
      outvalunit = 'km/s'
    case (3) ! y-velocity
      outval = 'v_y'
      outvalunit = 'km/s'
    case (4) ! z-velocity
      outval = 'v_z'
      outvalunit = 'km/s'
    case (5) ! Pressure
      outval = 'P'
      outvalunit = '?'
    case (6) ! Passive scalar
      outval = 'metal'
      outvalunit = '?'
    case (7) ! Temperature
      outval = 'T'
      outvalunit = 'K'
    case default ! density
      outval = 'rho'
      outvalunit = 'H/cc'
  end select

  if(makemap)then
    if(maxval) suffix = '_mv'//suffix
    if(maxrho) suffix = '_mr'//suffix
    suffix = '_'//TRIM(outval)//suffix
  endif

  ! Read hydro data
  i=INDEX(repository,'output_')
  nchar=repository(i+7:i+13)

  filename=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=filename, exist=ok)
  if(.not. ok)then
    write(*,*) "Error: ", trim(filename), " not found"
    stop
  endif

  filename=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=filename, exist=ok)
  if(.not. ok)then
    write(*,*) "Error: ", trim(filename), " not found"
    stop
  endif

  filename=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=1, file=filename, status='old', form='unformatted')
  read(1) ncpu
  read(1) ndim
  read(1) nx, ny, nz
  read(1) nlevelmax
  read(1) ngridmax
  read(1) nboundary
  read(1) ngrid_current
  read(1) boxlen
  close(1)

  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary > 0) allocate(ngridbound(1:nboundary,1:nlevelmax))

  filename=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=filename, exist=ok)
  if(.not. ok)then
    write(*,*) "Error: ", trim(filename), " not found"
    stop
  endif

  open(unit=1, file=filename, form='formatted', status='old')
  read(1,*)
  read(1,*)
!  read(1,'("levelmin    =",I11)')levelmin
  read(1,'(A13,I11)')GMGM,levelmin
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
!  read(1,'("time        =",E23.15)')t
  read(1,'(A13,E23.15)')GMGM,t
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
!  read(1,'("unit_l      =",E23.15)')scale_l
!  read(1,'("unit_d      =",E23.15)')scale_d
!  read(1,'("unit_t      =",E23.15)')scale_t
  read(1,'(A13,E23.15)')GMGM,scale_l
  read(1,'(A13,E23.15)')GMGM,scale_d
  read(1,'(A13,E23.15)')GMGM,scale_t
  read(1,*)
!  read(1,'("ordering type=",A80)'),ordering
  read(1,'(A13,A80)')GMGM,ordering
  read(1,*)
  

!  output Ramses Merger: kpc, 1e9 Msun
! conversion to CGS:
!  scale_l  = 3.08567752D21  ! = 1 kpc                   in cgs
!  scale_m  = 1.9889D42      ! = 1e9 Msun                in cgs
!  scale_d  = 6.77025D-23    ! = 1e9 Msun / kpc^3        in cgs
!  scale_t  = 4.7043D14      ! = ramses time unit        in cgs
!  scale_v  = 6.559269D6     ! = ramses velocity unit    in cgs
! conversion to useful units: multiply the Ramses output with the scale_X
!  scale_nH   = 30.996345 ! = ramses density     in H/cc
!  scale_vkms = 65.59269  ! = ramses velocity    in km/s
!  scale_T2   = 5.17302D5 ! = ramses temperature in Kelvin
!  scale_t    = 14.9070   ! = ramses time        in Myr

  scale_T2 = (scale_l / scale_t)**2 * 1.66D-24 / 1.3806200D-16
  scale_vkms = scale_l / scale_t / 1D5
  scale_nH = scale_d / 1.66D-24 * 0.76
  scale_l = scale_l / 3.085677581282D21
  scale_t = scale_t / 3.15576D13
! temperature * scale_T2 = K
! velocity * scale_vkms = km / s
! density * scale_nH = H/cm^3
! length * scale_l = kpc
! time * scale_t = Myr

  t = t * scale_t
    
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
    allocate(bound_key(0:ncpu))
    allocate(cpu_read(1:ncpu))
    cpu_read=.false.
    do impi=1,ncpu
      read(1,'(I8,1X,E23.15,1X,E23.15)') i, bound_key(impi-1), bound_key(impi)
    end do
  endif
  close(1)

  ! Compute map parameters
  if(lmax==0) lmax=nlevelmax
!  write(*,*)'time=',t
!  write(*,*)'Working resolution =',2**lmax
  
  if(ndim>2)then  
    select case (dir)
      case ('x')
        idim=2
        jdim=3
        kdim=1
        xxmin=ymin ; xxmax=ymax
        yymin=zmin ; yymax=zmax
        zzmin=xmin ; zzmax=xmax
      case ('y')
        idim=1
        jdim=3
        kdim=2
        xxmin=xmin ; xxmax=xmax
        yymin=zmin ; yymax=zmax
        zzmin=ymin ; zzmax=ymax
      case default
        idim=1
        jdim=2
        kdim=3
        xxmin=xmin ; xxmax=xmax
        yymin=ymin ; yymax=ymax
        zzmin=zmin ; zzmax=zmax
    end select

  else
    idim=1
    jdim=2
    xxmin=xmin ; xxmax=xmax
    yymin=ymin ; yymax=ymax          
! needed ?
    zzmin=0.0  ; zzmax=1.0
  end if

  if(TRIM(ordering).eq.'hilbert')then
    dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
    do ilevel=1,lmax
      dx=0.5d0**ilevel
      if(dx.lt.dmax) exit
    end do
  
    bit_length=ilevel-1
    maxdom=2**bit_length
    imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
    if(bit_length>0)then
      imin=int(xmin*dble(maxdom))
      imax=imin+1
      jmin=int(ymin*dble(maxdom))
      jmax=jmin+1
      kmin=int(zmin*dble(maxdom))
      kmax=kmin+1
    end if

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

    cpu_min=0
    cpu_max=0
    do impi=1,ncpu
      do i=1,ndom
        if (bound_key(impi-1).le.bounding_min(i).and.bound_key(impi).gt.bounding_min(i)) cpu_min(i)=impi
        if (bound_key(impi-1).lt.bounding_max(i).and.bound_key(impi).ge.bounding_max(i)) cpu_max(i)=impi
      end do
    end do
     
    ncpu_read=0
    do i=1,ndom
      do j=cpu_min(i), cpu_max(i)
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
  ! end if on hilbert

  ! Compute hierarchy
  do ilevel=1,lmax
    nx_full=2**ilevel
    ny_full=2**ilevel
    nz_full=2**ilevel
    imin=int(xxmin*dble(nx_full))+1
    imax=int(xxmax*dble(nx_full))+1
    jmin=int(yymin*dble(ny_full))+1
    jmax=int(yymax*dble(ny_full))+1
    allocate(grid(ilevel)%map(imin:imax,jmin:jmax))
    allocate(grid(ilevel)%rho(imin:imax,jmin:jmax))
    grid(ilevel)%map(:,:)=0.0
    grid(ilevel)%rho(:,:)=0.0
    grid(ilevel)%imin=imin
    grid(ilevel)%imax=imax
    grid(ilevel)%jmin=jmin
    grid(ilevel)%jmax=jmax    
    grid(ilevel)%kmin=int(zzmin*dble(nz_full))+1
    grid(ilevel)%kmax=int(zzmax*dble(nz_full))+1
  end do

  ! Compute projected variables
  
  ! open ascii file for particle output
  if(ascii) open(3, file='gas_part_'//TRIM(nchar)//TRIM(suffix)//'.ascii')

  ! Loop over cpu files
  do k=1,ncpu_read
    icpu=cpu_list(k)
    write(ncharcpu,'(I5.5)') icpu
    
    ! Open AMR file and skip header
    filename=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
    open(unit=1, file=filename, status='old', form='unformatted')
!    write(*,*)'Processing file '//TRIM(filename)
    do i=1,21
      read(1)
    end do
    ! Read grid numbers
    read(1) ngridlevel
    ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
    read(1)
    if(nboundary>0)then
      read(1)
      read(1)
      read(1) ngridbound
      ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
    endif
    read(1)     
    read(1) ! comment this line for old stuff
    if(TRIM(ordering).eq.'bisection')then
      do i=1,5
        read(1)
      end do
    else
      read(1)
    endif
    read(1)
    read(1)
    read(1)

    ! Open HYDRO file and skip header
    open(unit=2, file=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu), status='old', form='unformatted')
    read(2)
    read(2) nvarh
    read(2)
    read(2)
    read(2)
    read(2)

    ! Loop over levels
    do ilevel=1, lmax
      ! Geometry
      dx=0.5**ilevel
      dxline=1
      if(ndim==3) dxline=dx
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
        allocate(xg(1:ngrida,1:ndim))
        allocate(son(1:ngrida,1:twotondim))
        allocate(var(1:ngrida,1:twotondim,1:nvarh))
        allocate(x  (1:ngrida,1:ndim))
        allocate(rho(1:ngrida))
        allocate(map(1:ngrida))
        allocate(ref(1:ngrida))
      endif

      ! Loop over domains
      do j=1,nboundary+ncpu
        ! Read AMR data
        if(ngridfile(j,ilevel)>0)then
          read(1) ! Skip grid index
          read(1) ! Skip next index
          read(1) ! Skip prev index
          ! Read grid center
          do iidim=1,ndim
            if(j.eq.icpu)then
              read(1) xg(:,iidim)
            else
              read(1)
            endif
          end do
          read(1) ! Skip father index
          do ind=1,2*ndim
            read(1) ! Skip nbor index
          end do
          ! Read son index
          do ind=1,twotondim
            if(j.eq.icpu)then
              read(1) son(:,ind)
            else
              read(1)
            end if
          end do
          ! Skip cpu map
          do ind=1,twotondim
            read(1)
          end do
          ! Skip refinement map
          do ind=1,twotondim
            read(1)
          end do
        endif

        ! Read hydro data
        read(2)
        read(2)
        if(ngridfile(j,ilevel)>0)then
          ! Read hydro variables
          do ind=1,twotondim
            do ivar=1,nvarh
              if(j.eq.icpu)then
                read(2) var(:,ind,ivar)
              else
                read(2)
              end if
            end do
          end do
        end if
      end do
      ! end loop over domains
      
      ! Compute map
      if(ngrida>0)then

        ! Loop over cells
        do ind=1,twotondim
          ! Compute cell center
          do i=1,ngrida
            x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
            x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
            if(ndim>2)x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
          end do
          ! Check if cell is refined
          do i=1,ngrida
            ref(i)=son(i,ind)>0.and.ilevel<lmax
          end do
          ! Extract variable

          ! var(i,1): d, var(i,2:ndim+1): u,v,w and var(i,ndim+2): P.
          rho = var(:,ind,1)
          select case (typ)
            case (-1)
              map = icpu
            case (0)
              map = ilevel
            case (2) ! Mass weighted x-velocity
              map = var(:,ind,2)*var(:,ind,1)*scale_vkms
            case (3) ! Mass weighted y-velocity
              map = var(:,ind,3)*var(:,ind,1)*scale_vkms
            case (4) ! Mass weighted z-velocity
              map = var(:,ind,4)*var(:,ind,1)*scale_vkms
            case (5) ! Pressure
              map = var(:,ind,5)
            case (6) ! Passive scalar
              map = var(:,ind,6)*var(:,ind,1)
!              metmax=max(metmax,maxval(var(:,ind,6)))
            case (7) ! Mass weighted temperature
              map = var(:,ind,5)*scale_T2  ! / density * density   (divide to get the temperature from pressure, and multiply to have the mass-weighted temperature.)
!            case (7)
!              map = 0.5*(var(:,ind,1)**2+var(:,ind,2)**2+var(:,ind,3)**2)
            case default ! Mass weighted density
              map = var(:,ind,1)*var(:,ind,1)*scale_nH
          end select

          ! Store data map
          do i=1,ngrida
            if(.not.ref(i))then
              ix=int(x(i,idim)*dble(nx_full))+1
              iy=int(x(i,jdim)*dble(ny_full))+1
              iz=int(x(i,kdim)*dble(nz_full))+1

              ! 2D selection
              if(ix>=grid(ilevel)%imin.and.iy>=grid(ilevel)%jmin.and.ix<=grid(ilevel)%imax.and.iy<=grid(ilevel)%jmax) then  
                ! 3D selection
                if(iz>=grid(ilevel)%kmin.and.iz<=grid(ilevel)%kmax) then
                  ! Compute the PDF
                  if(pdf) then
                    if((var(i,ind,1)*scale_nH > pdfmin).and.(var(i,ind,1)*scale_nH < pdfmax)) then
                      ipdf = int( (log10(var(i,ind,1)*scale_nH)-lpdfmin) / lpdfampli )+1
                      pdfhist(ipdf) = pdfhist(ipdf) + var(i,ind,1) * (boxlen / 2.0**ilevel)**3
                    endif
                  endif

                  ! Output of particle data
                  if(ascii) then
                    ! x, y, z, vx(mass-weighted), vy(mass-weighted), vz(mass-weighted), rho, level
                    write(3, '(7e20.6e3,1I5)') x(i,1)*boxlen,x(i,2)*boxlen,x(i,3)*boxlen,var(i,ind,2)*scale_vkms,var(i,ind,3)*scale_vkms,var(i,ind,4)*scale_vkms,var(i,ind,1)*scale_nH,ilevel               
                  endif
                         
                endif
                ! end of 3D selection
                
                if(makemap)then
                  if(ndim==3)then
                    weight=(min(x(i,kdim)+dx/2.,zzmax)-max(x(i,kdim)-dx/2.,zzmin))/dx
                    weight=min(1.0d0,max(weight,0.0d0))
                  else
                    weight=1.0
                  endif

                  if(maxval)then
                    if(grid(ilevel)%map(ix,iy)<map(i))then
                      grid(ilevel)%map(ix,iy)=map(i) ! update the variable map
                      grid(ilevel)%rho(ix,iy)=rho(i) ! update the weight map with the density at *this* position
                    endif
                  else
                    if(maxrho)then
                      if(grid(ilevel)%rho(ix,iy)<rho(i))then
                        grid(ilevel)%map(ix,iy)=map(i) ! update the variable map
                        grid(ilevel)%rho(ix,iy)=rho(i) ! update the weight map with the density at *this* position
                      endif
                    else ! average
                      grid(ilevel)%map(ix,iy)=grid(ilevel)%map(ix,iy)+map(i)*dxline*weight/(zzmax-zzmin)
                      grid(ilevel)%rho(ix,iy)=grid(ilevel)%rho(ix,iy)+rho(i)*dxline*weight/(zzmax-zzmin)
                    endif
                  endif
                endif

              endif
              ! end of 2D selection
              
            end if
          end do
        end do
        ! End loop over cell
        deallocate(xg, son, var, ref, rho, map, x)
      end if
    end do
    ! End loop over levels

    close(1)
    close(2)

  end do
  ! End loop over cpu

  if(ascii) then
    close(3)
    stop
  endif
  
  if(pdf) then
    open(1, file='pdf_'//TRIM(nchar)//TRIM(suffix)//'.hist')
    do ipdf=1, pdfn
      write(1,*) 10**( (ipdf-1) * lpdfampli + lpdfmin ), pdfhist(ipdf) ! left side of the bin, bin value
    end do
    deallocate(pdfhist)
    close(1)
    stop
  endif

  nx_full=2**lmax
  ny_full=2**lmax
  imin=int(xxmin*dble(nx_full))+1
  imax=int(xxmax*dble(nx_full))
  jmin=int(yymin*dble(ny_full))+1
  jmax=int(yymax*dble(ny_full))

  do ix=imin,imax
    xmin=((ix-0.5)/2**lmax)
    do iy=jmin,jmax
      ymin=((iy-0.5)/2**lmax)
      do ilevel=1,lmax-1
        ndom=2**ilevel
        i=int(xmin*ndom)+1
        j=int(ymin*ndom)+1
        if(maxval) then
          if(grid(lmax)%map(ix,iy)<grid(ilevel)%map(i,j))then
            grid(lmax)%map(ix,iy)=grid(ilevel)%map(i,j) ! update the variable map
            grid(lmax)%rho(ix,iy)=grid(ilevel)%rho(i,j) ! update the weight map with the density at *this* position
          endif
        else
          if(maxrho)then
            if(grid(lmax)%rho(ix,iy)<grid(ilevel)%rho(i,j))then
              grid(lmax)%map(ix,iy)=grid(ilevel)%map(i,j) ! update the variable map
              grid(lmax)%rho(ix,iy)=grid(ilevel)%rho(i,j) ! update the weight map with the density at *this* position
            endif
          else ! average
            grid(lmax)%map(ix,iy)=grid(lmax)%map(ix,iy) + grid(ilevel)%map(i,j)
            grid(lmax)%rho(ix,iy)=grid(lmax)%rho(ix,iy) + grid(ilevel)%rho(i,j)
          endif
        endif
      end do
    end do
  end do

  if(nx_sample==0)then
    allocate(tmpmap(imax-imin+1,jmax-jmin+1))
    nxmap=max(imax-imin+1,jmax-jmin+1) ! projected average density
    if(typ > 0) then
      tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)/grid(lmax)%rho(imin:imax,jmin:jmax)
    else
      tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)
    endif    
  else
    if(ny_sample==0) ny_sample = nx_sample
    allocate(tmpmap(0:nx_sample,0:ny_sample))
    nxmap=max(nx_sample+1,ny_sample+1)
    do i=0,nx_sample
      ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
      ix=min(ix,imax)
      do j=0,ny_sample
        iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
        iy=min(iy,jmax)
        if(typ > 0) then
          tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)/grid(lmax)%rho(imin:imax,jmin:jmax)
        else
          tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)
        endif    
      end do
    end do
  endif

  allocate(map2(nxmap,nxmap))
  map2=tmpmap

  ! write data in a fits
  !status=0
  !filename='map_'//TRIM(nchar)//TRIM(suffix)//'.fits'
  !write(*,*) TRIM(filename)//' has been created.'
  
  !call deletefile(filename,status)
  !call ftgiou(unit,status)
  !blocksize=1
  !call ftinit(unit,filename,blocksize,status)
  !simple=.true.
  !bitpix=-32
  !naxis=2
  !naxes(1)=NXmap
  !naxes(2)=NXmap
  !extend=.true.
  !call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  
  !call ftpkyd(unit,'time',t,6,'time',status)
  !call ftpkyd(unit,'boxlen',boxlen,6,'boxlen',status)
  !call ftpkyd(unit,'xmin',xxmin,6,'xmin',status)
  !call ftpkyd(unit,'xmax',xxmax,6,'xmax',status)
  !call ftpkyd(unit,'ymin',yymin,6,'ymin',status)
  !call ftpkyd(unit,'ymax',yymax,6,'ymax',status)
  !call ftpkyd(unit,'zmin',zzmin,6,'zmin',status)
  !call ftpkyd(unit,'zmax',zzmax,6,'zmax',status)
  !call ftpkyj(unit,'lmax',lmax,'lmax',status)
  !call ftpkys(unit,'outval',outval,'value',status)
  !call ftpkys(unit,'outvalunit',outvalunit,'value unit',status)
  !call ftpkyl(unit,'maxval',maxval,'maxval',status)
  !call ftpkyl(unit,'maxrho',maxrho,'maxrho',status)
  
  !group=1
  !fpixel=1
  !nelements=naxes(1)*naxes(2)
  !call ftppre(unit,group,fpixel,nelements,MAP2,status)
  !call ftclos(unit, status)
  !call ftfiou(unit, status)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains
  subroutine read_params
    implicit none

    integer::iargc
    character(len=8)::opt
    character(len=128)::arg
    namelist /size/ dir, xmin, xmax, ymin, ymax, zmin, zmax
    namelist /amr/ ascii, pdf, lmax, typ, maxval, maxrho, pdfmin, pdfmax, pdfn

    n = iargc()
    if (n < 1) then
      print *, 'usage: amr2map -inp input_dir [-nml namelist] [-out suffix] [-ascii] [-pdf] '
      print *, ''
      print *, '   Set ascii to output the cell data as a particle list'
      print *, '   Set pdf to compute the pdf'
      print *, '   All namelist parameters can be forced if given in the cmd line.'
      print *, ''
      print *, '      -dir     projection (default: z)'
      print *, '      -xmin    selection, also for y and z (default: 0.0)'
      print *, '      -xmax    selection, also for y and z (default: 1.0)'
      print *, '      -lmax    maximum refinement level'
      print *, '      -typ     data to be plotted (see below)'
      print *, '      -maxval  get the maximum value along LOS (default: no)'
      print *, '      -maxrho  get the value where rho is maximum along LOS (default: no)'
      print *, '      -pdfmin  min density for the pdf [H/cc] (if pdf is set)'
      print *, '      -pdfmax  max density for the pdf [H/cc] (if pdf is set)'
      print *, '      -pdfn    number of log bin for the pdf (if pdf is set)'
      print *, ''
      print *, ' typ :-1 = cpu number'
      print *, '       0 = ref. level'
      print *, '       1 = gas density (default)'
      print *, '       2 = X velocity'
      print *, '       3 = Y velocity'
      print *, '       4 = Z velocity'
      print *, '       5 = gas pressure'
      print *, '       6 = gas metallicity'
      print *, '       7 = gas temperature'
      stop
    end if

    i = 1
    do while(i.le.n)
      call getarg(i,opt)
      select case (opt)
        case ('-inp')
          call getarg(i+1,arg)
          repository = trim(arg)        
        case ('-nml')
          call getarg(i+1,arg)
          open(1,file=trim(arg))
          read(1,size)
          read(1,amr)
          close(1)
        case ('-out')
          call getarg(i+1,arg)
          suffix = trim(arg)
          if(len(TRIM(suffix)).ne.0)suffix = '_'//TRIM(suffix)
        case ('-ascii')
          ascii = .true.
          i = i-1
        case ('-pdf')
          pdf = .true.
          i = i-1

        case ('-dir')
          call getarg(i+1,arg)
          dir = trim(arg) 
        case ('-xmin')
          call getarg(i+1,arg)
          read (arg,*) xmin
        case ('-xmax')
          call getarg(i+1,arg)
          read (arg,*) xmax
        case ('-ymin')
          call getarg(i+1,arg)
          read (arg,*) ymin
        case ('-ymax')
          call getarg(i+1,arg)
          read (arg,*) ymax
        case ('-zmin')
          call getarg(i+1,arg)
          read (arg,*) zmin
        case ('-zmax')
          call getarg(i+1,arg)
          read (arg,*) zmax
        case ('-lmax')
          call getarg(i+1,arg)
          read (arg,*) lmax
        case ('-typ')
          call getarg(i+1,arg)
          read (arg,*) typ
        case ('-maxval')
          maxval = .true.
          i = i-1
        case ('-maxrho')
          maxrho = .true.
          i = i-1
        case ('-pdfmin')
          call getarg(i+1,arg)
          read (arg,*) pdfmin
        case ('-pdfmax')
          call getarg(i+1,arg)
          read (arg,*) pdfmax
        case ('-pdfn')
          call getarg(i+1,arg)
          read (arg,*) pdfn
        case default
          print '("unknown option ",a8," ignored")', opt
          i = i-1
      end select
      i = i+2
    end do
    return
    
  end subroutine read_params
  
end program amr2map


!=======================================================================
!=======================================================================
!=======================================================================
!
!subroutine deletefile(filename,status) !  Delete a FITS file
!
!  integer::status,unit,blocksize
!  character(*)::filename
!  
!  if (status .gt. 0) return
!
!  call ftgiou(unit,status) ! Get an unused Logical Unit Number
!  call ftopen(unit,filename,1,blocksize,status) ! Try to open the file
!  
!  if (status .eq. 0)then ! file is opened: delete it 
!    call ftdelt(unit,status)
!  else if (status .eq. 103)then ! file doesn't exist: reset status and clear errors
!    status=0
!    call ftcmsg
!  else ! there was some other error opening the file: delete the file anyway
!    status=0
!    call ftcmsg
!    call ftdelt(unit,status)
!  end if
!  
!  call ftfiou(unit, status) ! Free the unit number
!
!end

!=======================================================================
!=======================================================================
!=======================================================================

subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer,intent(in)::bit_length, npoint
  integer,intent(in),dimension(1:npoint)::x, y, z
  real(kind=8),intent(out),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask, y_bit_mask, z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i, ip, cstate, nstate, b0, b1, b2, sdigit, hdigit

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
      b2=0
      if(i_bit_mask(3*i+2))b2=1
      b1=0
      if(i_bit_mask(3*i+1))b1=1
      b0=0
      if(i_bit_mask(3*i  ))b0=1
      
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
