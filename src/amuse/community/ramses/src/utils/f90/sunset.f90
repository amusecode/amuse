program sunset
  !------------------------------------------------------------------------------------------------------------------
  ! Sunset is a program that computes realistic galaxy images using spectra
  ! of a single stellar population (SSP) of your favorite model (GALAXEV, STARDUST, PEGASE).
  !
  ! The images are computed for various filters (U,V,B,R,I,J...) and take into account,
  ! if desired, dust attenuation. In order to work, you need to provide the following data:
  ! - a particle file containing stars, with position, velocity, birth date and metallicity
  ! - a SSP spectral 3D cube (age, metallicity, lambda)
  ! _ the spectral response datebase for each filer
  ! If you need dust absorption, you also need:
  ! - a dust opacity model (kappa g/cm2 versus lambda) 
  ! - a dust mass density 3D array covering the same bounding box than the image
  ! The output data are pretty standard stuff.
  !
  ! The current version works with:
  ! - any RAMSES particle file
  ! - a simplified STARDUST model for SSP (download at http://www.ucolick.org/~patrik/sunrise/)
  ! - any dust Draine model for the MW, LMC or SMC (download at http://www.astro.princeton.edu/~draine/dust/dustmix.html)
  ! - the PEGASE filter response function database (download at http://www2.iap.fr/pegase/)
  ! - a binary 3D array of metal mass density from RAMSES data.
  !
  ! Don't forget to check your unit system ! Good luck with sunset ! 
  !
  !                                              Zurich, Nov. 1st, 2008
  !                                              Romain Teyssier
  !------------------------------------------------------------------------------------------------------------------
  implicit none
  integer,parameter::NDUSTMAX=10000
  integer,parameter::nmaxlambda=2800
  integer,parameter::NMAXFILTERS=100
  integer,parameter::nmaxlines=100
  integer,parameter::nmaxtimes=100
  
  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,nstar,nfilters
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii,n1,n2,n3
  integer::nx=0,ny=0,nz=0,ix,iy,iz,ixp1,iyp1,idim,jdim,kdim,ncpu_read,n_frw
  real(KIND=8)::mtot,ddx,ddy,ddz,dex,dey,dez,time,age,time_tot,time_simu,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,col,dust_opacity
  real(KIND=8)::max_column,intssp_dust,boxlen,xx,yy
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin
  integer::ndust,nmetal,ntime,nlambda,ivel,nvel,ihist,nhist=10000
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,dz,deltax,vrot_tot,velocity,vel_min,vel_max
  real(KIND=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d,lll,fff,unit_m
  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:),allocatable::map
  real(KIND=8),dimension(:,:),allocatable::x,v
  real(KIND=8),dimension(:)  ,allocatable::m,bdate,z
  integer,dimension(:)  ,allocatable::id
  real(KIND=8),dimension(:)  ,allocatable::lssp,tssp,mssp,sed,vrot,kappa_dust,filter_contrib
  real(KIND=8),dimension(:)  ,allocatable::time_hist,metal_hist
  real(KIND=8),dimension(:,:),allocatable::intssp
  real(KIND=8),dimension(:,:,:),allocatable::ssp,column
  real(KIND=4),dimension(:,:,:),allocatable::rho

  real(KIND=8),dimension(1:NDUSTMAX)::ldust,kdust
  real(KIND=8),dimension(1:NMAXFILTERS,1:1000)::lambdafilter,transfilter
  integer,dimension(1:NMAXFILTERS)::nlambdafilter,typetrans,typecalib

  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering,format_grille,str
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  character(LEN=128)::dustfile='/home/rteyssie/sunset/dust/kext_albedo_WD_MW_3.1_60.txt'
  character(LEN=128)::sspfile='/home/rteyssie/sunset/stardust99/Patrik-imfKroupa-Zmulti-subsampled.dat'
!  character(LEN=128)::sspfile='/home/rteyssie/sunset/bc03/bruzual_charlot_2003.dat'
  character(LEN=128)::filterfile='/home/rteyssie/sunset/pegase/filters.dat'
  character(LEN=128)::rhofile=''
  character(LEN=10),dimension(1:NMAXFILTERS)::fname
  character(LEN=10)::filtername
  logical::ok,ok_part,star=.false.,ageweight=.false.,do_dust=.false.,err
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,dum,kappa,lambda,normfilter,nu
  real(KIND=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  integer::ilambda_min,ilambda_max,ifilter=6,itime,imetal,ilambda
  logical::cosmo=.true.,metal=.true.

  call read_params

  ! Define velocity space
  nvel=100
  vel_min=-500.
  vel_max=+500.
  allocate(vrot(1:nvel))
  vrot=0d0

  !------------------------------------------
  ! Read filter transmission from PEGASE file
  !------------------------------------------
  write(*,*)'Read filters for the following bands'
  open(30,status='old',file=filterfile)
  read(30,*)nfilters
  ifilter=1
  do i=1,nfilters
     read(30,*) nlambdafilter(i),typetrans(i),typecalib(i),fname(i)
     write(*,*)i,fname(i)
     if(TRIM(fname(i))==TRIM(filtername))ifilter=i
     do j=1,nlambdafilter(i)
        read(30,*)lambdafilter(i,j),transfilter(i,j)
     end do
  end do  
  close(30)
  write(*,*)'Working in band',ifilter,TRIM(filtername)

  !-----------------------------------
  ! Read dust opacity file from Draine
  ! Units are micron and cm2/g
  !-----------------------------------
  if(do_dust)then
     open(1,file=dustfile,form='formatted',status='old')
     err=.true.
     do while(err)
        read(1,'(A80)',END=81)str
        if(TRIM(str)==' lambda   albedo    g     C_ext/H    K_abs')err=.false.
     end do
!     write(*,*)'  i    '//TRIM(str)
     read(1,*)
     read(1,*)
     err=.true.
     n=0
     do while(err)
        read(1,*,END=81)lambda,dum,dum,dum,kappa
!        write(*,'(I4," ",5(1PE10.3," "))')n,lambda,kappa
        n=n+1
        ldust(n)=lambda
        kdust(n)=kappa
     end do
81   continue
     close(1)
     ndust=n
     ! Convert microns to Angstroms
     ldust(1:ndust)=ldust(1:ndust)*1d4
  end if

  !----------------------------------------------
  ! Read SSP file (STARDUST 99 format)
  ! The mass of the SSP is 10^6 Msol.
  ! Fluxes are in erg/s/A.
  ! lambda are in Angstroms, 
  ! ages in years, 
  ! metallicity in solar units
  !----------------------------------------------
  write(*,*)'Reading spectra from file '//TRIM(sspfile)
  open(1,file=sspfile,form='unformatted',status='old')
  read(1)nlambda,ntime,nmetal
  allocate(lssp(1:nlambda))
  allocate(sed (1:nlambda))
  allocate(filter_contrib(1:nlambda))
  sed=0d0
  allocate(tssp(1:ntime))
  allocate(time_hist(0:nhist))
  allocate(metal_hist(0:nhist))
  time_hist=0d0
  metal_hist=0d0
  allocate(mssp(1:nmetal))
  allocate(ssp(1:ntime,1:nmetal,1:nlambda))
  read(1)lssp
  read(1)tssp
  read(1)mssp
  read(1)ssp
  close(1)

  !-------------------------------
  ! Read RAMSES info (header) file
  !-------------------------------
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
  if(h0==1)cosmo=.false.
  read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'("unit_t      =",E23.15)')unit_t
  unit_m=unit_d*unit_l**3
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

  !----------------------------------------------
  ! Convert log lum into flux per code mass unit 
  !----------------------------------------------
  do itime=1,ntime
  do imetal=1,nmetal
  do ilambda=1,nlambda
     ssp(itime,imetal,ilambda)=unit_d*unit_l**3/1d6/2d33*10d0**ssp(itime,imetal,ilambda) &
          & /(4.*3.1415926*(10.*3.08d18)**2)
  end do
  end do
  end do

  !----------------------------------------------
  ! Compute dust opacity at the sed wavelength
  !----------------------------------------------
  if(do_dust)then
     allocate(kappa_dust(1:nlambda))
     ilambda_min=1
     do ilambda=1,nlambda
        j=ilambda_min
        do while(ldust(j)<lssp(ilambda))
           j=j+1
        end do
        lll=lssp(ilambda)
        if(j==1)then
           kappa_dust(ilambda)=0d0
        else
           kappa_dust(ilambda)=kdust(j-1)*(lll-ldust(j))/(ldust(j-1)-ldust(j))+ &
                & kdust(j)*(lll-ldust(j-1))/(ldust(j)-ldust(j-1))
        endif
        ilambda_min=j
        !     write(77,*)lll,kappa_dust(ilambda)
     end do
  endif

  !----------------------------------------------
  ! Compute ssp emission in the chosen filter
  !----------------------------------------------
  allocate(intssp(1:ntime,1:nmetal))
  intssp=0d0
  ilambda_min=1
  do while(lssp(ilambda_min)<lambdafilter(ifilter,1))
     ilambda_min=ilambda_min+1
  end do
  ilambda_min=ilambda_min-1
  ilambda_max=nlambda
  do while(lssp(ilambda_max)>lambdafilter(ifilter,nlambdafilter(ifilter)))
     ilambda_max=ilambda_max-1
  end do
  ilambda_max=ilambda_max+1

  ! Normalize filter
  do ilambda=ilambda_min,ilambda_max
     lll=lssp(ilambda)
     iii=int(dble(nlambdafilter(ifilter)-1)*(lll-lambdafilter(ifilter,1))/ &
          & (lambdafilter(ifilter,nlambdafilter(ifilter))-lambdafilter(ifilter,1)))
     iii=iii+1
     iii=max(iii,1)
     iii=min(iii,nlambdafilter(ifilter)-1)
     fff=transfilter(ifilter,iii)*(lll-lambdafilter(ifilter,iii+1)) &
          & /(lambdafilter(ifilter,iii)-lambdafilter(ifilter,iii+1))+ &
          & transfilter(ifilter,iii+1)*(lll-lambdafilter(ifilter,iii)) &
          & /(lambdafilter(ifilter,iii+1)-lambdafilter(ifilter,iii))
     fff=max(fff,0d0)
     normfilter=normfilter+fff*(lssp(ilambda+1)-lssp(ilambda))*2.0/(lssp(ilambda+1)+lssp(ilambda))
  end do

  do itime=1,ntime
     do imetal=1,nmetal
        do ilambda=ilambda_min,ilambda_max
           lll=lssp(ilambda)
           iii=int(dble(nlambdafilter(ifilter)-1)*(lll-lambdafilter(ifilter,1))/ &
                & (lambdafilter(ifilter,nlambdafilter(ifilter))-lambdafilter(ifilter,1)))
           iii=iii+1
           iii=max(iii,1)
           iii=min(iii,nlambdafilter(ifilter)-1)
           fff=transfilter(ifilter,iii)*(lll-lambdafilter(ifilter,iii+1)) &
                & /(lambdafilter(ifilter,iii)-lambdafilter(ifilter,iii+1))+ &
                & transfilter(ifilter,iii+1)*(lll-lambdafilter(ifilter,iii)) &
                & /(lambdafilter(ifilter,iii+1)-lambdafilter(ifilter,iii))
           fff=max(fff,0d0)
           nu=3d10/1d-8*2.0/(lssp(ilambda+1)+lssp(ilambda))
           intssp(itime,imetal)=intssp(itime,imetal)+fff*ssp(itime,imetal,ilambda) &
                & *(lssp(ilambda+1)-lssp(ilambda))/normfilter/nu
           filter_contrib(ilambda)=fff*(lssp(ilambda+1)-lssp(ilambda))/normfilter/nu
        end do
     end do
  end do

  !-----------------------------
  ! Compute cosmological model
  !-----------------------------
  if(cosmo)then
     ! Allocate look-up tables
     n_frw=1000
     allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
     allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
     
     ! Compute Friedman model look up table
     write(*,*)'Computing Friedman model'
     call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
          & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
     ! Find neighboring conformal time
     i=1
     do while(tau_frw(i)>t.and.i<n_frw)
        i=i+1
     end do
     ! Interploate time
     time_simu=t_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
          & t_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
     write(*,*)'Time simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
     write(*,*)'Hubble time=',(time_tot)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
  else
     write(*,*)'Time simu=',t*unit_t/(365.*24.*3600.*1d9)
     time_simu=t
  endif
  !-----------------------
  ! Set map parameters
  !-----------------------
  if(nx==0)then
     nx=2**levelmin
  endif
  if(ny==0)then
     ny=nx
  end if
  if(nz==0)then
     nz=ny
  end if
  write(*,*)'Working cube =',nx,ny,nz
  allocate(map(0:nx,0:ny))
  allocate(column(0:nx-1,0:ny-1,0:nz))
  column=0d0
  map=0.0d0
  if (proj=='x')then
     idim=2
     jdim=3
     kdim=1
     xxmin=ymin ; xxmax=ymax
     yymin=zmin ; yymax=zmax
     zzmin=xmin ; zzmax=xmax
  else if (proj=='y') then
     idim=1
     jdim=3
     kdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=zmin ; yymax=zmax
     zzmin=ymin ; zzmax=ymax
  else
     idim=1
     jdim=2
     kdim=3
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
     zzmin=zmin ; zzmax=zmax
  end if
  dx=(xxmax-xxmin)/dble(nx)
  dy=(yymax-yymin)/dble(ny)
  dz=(zzmax-zzmin)/dble(nz)

  !-----------------------------
  ! Read dust mass density array
  !-----------------------------
  if(do_dust)then
     write(*,*)'Reading dust density in file '//TRIM(rhofile) 
     open(unit=20,file=rhofile,form='unformatted')
     read(20)n1,n2,n3
     if (proj=='x') then
        if(n1.ne.ny.or.n2.ne.nz.or.n3.ne.nx)then
           write(*,*)'Dimension of density file not compatible'
           write(*,*)'stop'
           stop
        endif
     else if (proj=='y') then
        if(n1.ne.nx.or.n2.ne.nz.or.n3.ne.ny)then
           write(*,*)'Dimension of density file not compatible'
           write(*,*)'stop'
           stop
        endif
     else if (proj=='z') then
        if(n1.ne.nx.or.n2.ne.ny.or.n3.ne.nz)then
           write(*,*)'Dimension of density file not compatible'
           write(*,*)'stop'
           stop
        endif
     end if
     allocate(rho(1:n1,1:n2,1:n3))
     read(20)rho
     close(20)
     do i=1,n1
     do j=1,n2
     do k=1,n3
        if (proj=='x') then
           ix=j-1
           iy=k-1
           iz=i
        else if (proj=='y') then
           ix=i-1
           iy=k-1
           iz=j
        else if (proj=='z') then
           ix=i-1
           iy=j-1
           iz=k
        end if
        column(ix,iy,iz)=unit_d*rho(i,j,k)
     end do
     end do
     end do
     deallocate(rho)
     write(*,*)'Computing column density'
     max_column=0d0
     do ix=0,nx-1
     do iy=0,ny-1
     do iz=1,nz
        column(ix,iy,iz)=column(ix,iy,iz-1)+column(ix,iy,iz)*dz*unit_l
        max_column=max(max_column,column(ix,iy,iz))
     end do
     end do
     end do
     write(*,*)'Maximum dust column density ',max_column/1.66d-24
  end if

  !-----------------------------------
  ! Compute Hilbert key of the domain
  !-----------------------------------
  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,levelmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
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
     
     dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
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

  !------------------------------------------
  ! Read total number of particle in files
  !------------------------------------------
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
  if(nstar==0)then
     write(*,*)'Found no star particle'
     stop
  endif

  !-----------------------------------------------
  ! Compute image in the chosen filter
  ! While reading particle in chunks
  !----------------------------------------------
  mtot=0.0d0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     allocate(id(1:npart2))
     allocate(bdate(1:npart2))
     allocate(z(1:npart2))
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
     ! Read id
     read(1)id 
     read(1) ! Skip level
     ! Read birth date
     read(1)bdate
     ! Read metallicity
     if(metal)then
        read(1)z
     else
        z=0.02
     endif
     close(1)

     dust_opacity=1d0
     do i=1,npart2
        ok_part=(x(i,1)>=xmin.and.x(i,1)<xmax.and. &
             &   x(i,2)>=ymin.and.x(i,2)<ymax.and. &
             &   x(i,3)>=zmin.and.x(i,3)<zmax)
        
        ok_part=ok_part.and.(bdate(i).ne.0.0d0)!.and.id(i).gt.0

        if(ok_part)then

           ! Compute star age in years
           if(cosmo)then
              iii=1
              do while(tau_frw(iii)>bdate(i).and.iii<n_frw)
                 iii=iii+1
              end do
              time=t_frw(iii)*(bdate(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                   & t_frw(iii-1)*(bdate(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
              age=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.)        
              time=(time_tot+time)/(h0*1d5/3.08d24)/(365.*24.*3600.)        
           else
              time=(bdate(i))*unit_t/(365.*24.*3600.)
              age=(time_simu-bdate(i))*unit_t/(365.*24.*3600.)
           endif

!           if(age<15d6)m(i)=0d0

           itime=1
           do while(tssp(itime)<age.and.itime<ntime)
              itime=itime+1
           end do

           ihist=int(time/15d9*dble(nhist))
           time_hist(ihist)=time_hist(ihist)+m(i)*unit_m/2d33

           ! Compute star metallicity in solar units
           imetal=1
           do while(mssp(imetal)<z(i).and.imetal<nmetal)
              imetal=imetal+1
           end do

           ihist=int(z(i)/0.1*dble(nhist))
           metal_hist(ihist)=metal_hist(ihist)+m(i)*unit_m/2d33

           ddx=(x(i,idim)-xxmin)/dx
           ddy=(x(i,jdim)-yymin)/dy
           ddz=(x(i,kdim)-zzmin)/dz
           ix=ddx
           iy=ddy
           iz=ddz
           dex=ddx-ix
           dey=ddy-iy
           dez=ddz-iz

           if(ix>=0.and.ix<nx.and.iy>=0.and.iy<ny.and.iz>=0.and.iz<nz)then
              if(do_dust)then
                 col=column(ix,iy,iz+1)!(1d0-dez)*column(ix,iy,iz)+dez*column(ix,iy,iz+1)
              else
                 map(ix,iy)=map(ix,iy)+m(i)*intssp(itime,imetal)
                 mtot=mtot+m(i)*intssp(itime,imetal)
              endif
              do ilambda=1,nlambda
                 dust_opacity=1d0
                 if(do_dust)then
                    dust_opacity=exp(-kappa_dust(ilambda)*col)
                 endif
                 sed(ilambda)=sed(ilambda)+m(i)*ssp(itime,imetal,ilambda)*dust_opacity
              end do

              if(do_dust)then
                 intssp_dust=0d0
                 do ilambda=ilambda_min,ilambda_max
                    dust_opacity=exp(-kappa_dust(ilambda)*col)
                    intssp_dust=intssp_dust+filter_contrib(ilambda)*ssp(itime,imetal,ilambda)*dust_opacity
                 end do
                 map(ix,iy)=map(ix,iy)+m(i)*intssp_dust
                 mtot=mtot+m(i)*intssp_dust
              end if

              velocity=v(i,kdim)*unit_l/unit_t/1d5
              ivel=int((velocity-vel_min)/(vel_max-vel_min)*dble(nvel))
              ivel=min(max(ivel,0),nvel-1)+1
              vrot(ivel)=vrot(ivel)+m(i)*intssp(itime,imetal)
              vrot_tot=vrot_tot+m(i)*intssp(itime,imetal)
           endif

        end if

     end do
     deallocate(x,m,v,id)
     deallocate(bdate,z)
  end do

  !-------------------------------------
  ! Output total magnitude in the band
  !-------------------------------------
  write(*,*)'Total luminosity=',mtot
  write(*,*)'Magnitude=',-2.5*log10(mtot)-48.6

  !----------------------------------------------
  ! Output SED to file in units of erg/sec/cm2/A
  !----------------------------------------------
  open(10,file='sed.dat',form='formatted')
  do ilambda=1,nlambda
     write(10,*)lssp(ilambda),log10(sed(ilambda))
  end do
  close(10)

  !----------------------------------------------
  ! Output mass weighted time histogram
  !----------------------------------------------
  open(10,file='time_histo.dat',form='formatted')
  do ihist=0,nhist
     write(10,*)dble(ihist)/dble(nhist)*15d9,time_hist(ihist)
  end do
  close(10)

  !----------------------------------------------
  ! Output mass weighted metallicity histogram
  !----------------------------------------------
  if(metal)then 
     open(10,file='metal_histo.dat',form='formatted')
     do ihist=0,nhist
        write(10,*)dble(ihist)/dble(nhist)*0.1,metal_hist(ihist)
     end do
     close(10)
  endif

  !----------------------------------------------
  ! Output luminosity weighted velocity histogram
  !----------------------------------------------
  open(10,file='vrot.dat',form='formatted')
  do ivel=1,nvel
     write(10,*)vel_min+(dble(ivel)-0.5)*(vel_max-vel_min)/dble(nvel),vrot(ivel)/vrot_tot
  end do
  close(10)

  !---------------------
  ! Output map to file
  !---------------------
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  if(TRIM(filetype).eq.'bin')then
     open(unit=10,file=nomfich,form='unformatted')
     write(10)nx,ny
     allocate(toto(nx,ny))
     toto=log10(map(0:nx-1,0:ny-1)+1d-20)
     write(10)toto
     close(10)
  endif
  if(TRIM(filetype).eq.'asc'.or.TRIM(filetype).eq.'ascii')then
     open(unit=10,file=nomfich,form='formatted')
     allocate(toto(nx,ny))
     toto=log10(map(0:nx-1,0:ny-1)+1d-20)
     do j=0,ny
        do i=0,nx
           xx=xxmin+dble(i)/dble(nx)*(xxmax-xxmin)
           yy=yymin+dble(j)/dble(ny)*(yymax-yymin)
           write(10,*)xx,yy,toto(i,j)
        end do
        write(10,*) " "
        end do
     close(10)
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
       print *, 'usage: sunset    -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-dir axis] '
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-nx  nx  ] '
       print *, '                 [-ny  ny  ] '
       print *, '                 [-fil filetype] '
       print *, '                 [-bnd filtername] '
       print *, '                 [-dst rhodustfilename] '
       print *, 'ex: sunset -inp output_00001 -out map.dat'// &
            &   ' -dir z -xmi 0.1 -xma 0.7'
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
       case ('-dir')
          proj = trim(arg) 
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
       case ('-nx')
          read (arg,*) nx
       case ('-ny')
          read (arg,*) ny
       case ('-fil')
          filetype = trim(arg)
       case ('-bnd')
          filtername = trim(arg)
       case ('-dst')
          do_dust=.true.
          rhofile = trim(arg)
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params
  
end program sunset

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




