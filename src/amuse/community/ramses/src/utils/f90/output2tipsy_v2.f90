program output2tipsy_v2
  use random 
  use io_ramses

  implicit none

  integer::ndim,n,i,j,k,twotondim,ncoarse,indcell
  integer::ivar,nvar,ncpu,lmax=20,levelmin
  integer::nx,ny,nz
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ipos
  integer::ngridmax,nstep_coarse,icpu
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  real(KIND=8)::gamma,smallr,smallc,gammah
  real(KIND=8)::boxlen,boxlen2
  real(KIND=8)::t,aexp,hexp,t2,aexp2,hexp2
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
  real(KIND=8)::dkey,order_min,dmax,ddx,dxline,ddy,dex,dey,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,mdm

  integer::ilevel,ncpu_read
  real(kind=8)deltax
  character(len=3)::typ='all'
  character(LEN=5)::nchar
  character(LEN=80)::ordering
  character(LEN=128)::nomfich
  character(LEN=128)::repository,outfich,filetype='bin'
  character(LEN=13)::string
  logical::ok,ok_part,ok_cell,do_max
  real(kind=8),dimension(:),allocatable::bound_key,xdp
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list

  integer::ndummy=64,ndummypart,nmin,nmax,nold,nnold,ndummyold
  real(KIND=8),dimension(:,:),allocatable::xp
  real(KIND=8),dimension(:,:),allocatable::varp
  integer::partcount,respart,denspartcount
  real(KIND=8)::dummy,partmass,volume,facdens=0.d0,averdens
  integer::delm=0,levelsel

  integer::npart,nstar_tot,nsink,npart_tot,nsink_tot
  integer::npart_actual,ndm_actual,nstar_actual,npartnow
  integer,dimension(:),allocatable::idpart
  real(KIND=8),dimension(:,:),allocatable::xpart,vpart
  real(KIND=8),dimension(:),allocatable::mpart,age,met
  character(LEN=5)::nn
  logical::metal=.true.,star=.true.,sink=.true.
  
  integer ,dimension(1:1,1:IRandNumSize)::allseed
  integer ,dimension(1:IRandNumSize)::localseed
  integer::iseed=0,poisson

  ndummypart=ndummy**3
  ndummyold=ndummy
  nold=1
  nnold=ndummypart
  nmin=nold
  nmax=nnold

  call read_params

  ndummypart=ndummy**3
  if(nmin==nold.and.nmax==nnold.and.ndummy/=ndummyold)then
     nmin=1
     nmax=ndummypart
  endif

  !-----------------------------------------------
  ! Reading files in RAMSES format
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     typ='gas'
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
  read(10,'(A13,E23.15)')string,hexp
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
  ndummypart=ndummy**3

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,lmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lllmin=ilevel!max(ilevel,levelsel)
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

!MASS IN CGS UNITS: write(*,*)partmass*scale_d*(scale_l**3)

  if(typ=='gas'.and.ndim==3) then
     call gaspart(ncpu,ncpu_read,cpu_list,repository,ordering,ndummypart,facdens,&
          levelsel,lmax,xmin,xmax,ymin,ymax,zmin,zmax,nmin,nmax,&
          partmass,averdens,xp,varp,denspartcount)
     open(66,file=outfich,status='unknown',form='formatted')
     !  write(66,*)ndummypart,ndummypart,0
     write(66,*)denspartcount,denspartcount,0
     write(66,*)ndim
     write(66,*)t
     do i=1,nmax-nmin+1
        write(66,*)partmass
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)xp(i,1)
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)xp(i,2)
        !     write(66,*)yp(i)
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)xp(i,3)
        !     write(66,*)zp(i)
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)varp(i,2)
        !     write(66,*)varp(i,2)
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)varp(i,3)
        !     write(66,*)varp(i,3)
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)varp(i,4)
        !     write(66,*)varp(i,4)
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)varp(i,1)
        !     write(66,*)varp(i,1)
     end do
     do i=1,nmax-nmin+1
        if(varp(i,1)>facdens*averdens)write(66,*)varp(i,5)/(gamma-1.d0)/varp(i,1)
        !     write(*,*)varp(i,5)/(dble(gamma)-1.d0)/varp(i,1)
     end do
     if(nvarh>=6)then
        do i=1,nmax-nmin+1
           if(varp(i,1)>facdens*averdens)write(66,*)0.d0
        end do
        do i=1,nmax-nmin+1
           if(varp(i,1)>facdens*averdens)write(66,*)varp(i,6)
        end do
     endif
     close(66)
  end if

  if(typ=='all'.and.ndim==3) then
     call readpart(ncpu,ncpu_read,cpu_list,ndim,repository,metal,star,sink,&
          lmin,lmax,xmin,xmax,ymin,ymax,zmin,zmax,nmin,nmax,npart_actual,&
          ndm_actual,nstar_actual,xpart,vpart,mpart,idpart,age,met)

     mdm=1d30
     do i=1,npart_actual
        if(idpart(i)>=0.and.age(i)==0.d0)mdm=min(mdm,mpart(i))
     end do
     mdm=mdm*omega_b/(omega_m-omega_b)

     call gaspart3(ncpu,ncpu_read,cpu_list,repository,ordering,ndummypart,facdens,&
          levelsel,lmax,xmin,xmax,ymin,ymax,zmin,zmax,mdm,&
          partmass,averdens,xp,varp,denspartcount)

     nmin=1
     nmax=denspartcount

     !-------------------------------------------------------------
     !  Writing output tipsy file 
     !-------------------------------------------------------------
  
123  open(66,file=outfich,status='unknown',form='formatted')
     open(55,file='partID_'//TRIM(nchar),status='unknown',form='formatted')
     !  write(66,*)ndummypart,ndummypart,0
     write(66,*)npart_actual+denspartcount,denspartcount,nstar_actual
     write(66,*)ndim
     write(66,*)t

     write(55,*)npart_actual+denspartcount,denspartcount,nstar_actual
     write(55,*)ndim
     write(55,*)t

     do i=1,nmax-nmin+1
        write(66,*)partmass
     end do
     do i=1,npart_actual
        if((.not.star).and.(.not.sink))write(66,*)mpart(i)
        if((.not.star).and.(sink).and.idpart(i)>0)write(66,*)mpart(i)
        if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(66,*)mpart(i)
        if((.not.star).and.(.not.sink))write(55,*)idpart(i)
        if((.not.star).and.(sink).and.idpart(i)>0)write(55,*)idpart(i)
        if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(55,*)idpart(i)
     enddo
     if(star.and.nstar_actual>0)then
        do i=1,npart_actual
           if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)mpart(i)
           if(age(i)/=0.d0.and.idpart(i)>0)write(55,*)idpart(i)
        enddo
     endif
     do i=1,nmax-nmin+1
        if(varp(i,1)>=facdens*averdens)write(66,*)xp(i,1)
     end do
     do i=1,npart_actual
        if((.not.star).and.(.not.sink))write(66,*)xpart(i,1)
        if((.not.star).and.(sink).and.idpart(i)>0)write(66,*)xpart(i,1)
        if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(66,*)xpart(i,1)
     enddo
     if(star.and.nstar_actual>0)then
        do i=1,npart_actual
           if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)xpart(i,1)
        enddo
     endif
     do i=1,nmax-nmin+1
        if(varp(i,1)>=facdens*averdens)write(66,*)xp(i,2)
     end do
     do i=1,npart_actual
        if((.not.star).and.(.not.sink))write(66,*)xpart(i,2)
        if((.not.star).and.(sink).and.idpart(i)>0)write(66,*)xpart(i,2)
        if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(66,*)xpart(i,2)
     enddo
     if(star.and.nstar_actual>0)then
        do i=1,npart_actual
           if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)xpart(i,2)
        enddo
     endif

  do i=1,nmax-nmin+1
     if(varp(i,1)>=facdens*averdens)write(66,*)xp(i,3)
  end do
  do i=1,npart_actual
     if((.not.star).and.(.not.sink))write(66,*)xpart(i,3)
     if((.not.star).and.(sink).and.idpart(i)>0)write(66,*)xpart(i,3)
     if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(66,*)xpart(i,3)
  enddo
  if(star.and.nstar_actual>0)then
     do i=1,npart_actual
        if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)xpart(i,3)
     enddo
  endif

  do i=1,nmax-nmin+1
     if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,2)
  end do
  do i=1,npart_actual
     if((.not.star).and.(.not.sink))write(66,*)vpart(i,1)
     if((.not.star).and.(sink).and.idpart(i)>0)write(66,*)vpart(i,1)
     if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(66,*)vpart(i,1)
  enddo
  if(star.and.nstar_actual>0)then
     do i=1,npart_actual
        if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)vpart(i,1)
     enddo
  endif

  do i=1,nmax-nmin+1
     if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,3)
  end do
  do i=1,npart_actual
     if((.not.star).and.(.not.sink))write(66,*)vpart(i,2)
     if((.not.star).and.(sink).and.idpart(i)>0)write(66,*)vpart(i,2)
     if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(66,*)vpart(i,2)
  enddo
  if(star.and.nstar_actual>0)then
     do i=1,npart_actual
        if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)vpart(i,2)
     enddo
  endif

  do i=1,nmax-nmin+1
     if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,4)
  end do
  do i=1,npart_actual
     if((.not.star).and.(.not.sink))write(66,*)vpart(i,3)
     if((.not.star).and.(sink).and.idpart(i)>0)write(66,*)vpart(i,3)
     if((star.or.sink).and.age(i)==0.d0.and.idpart(i)>0)write(66,*)vpart(i,3)
  enddo
  if(star.and.nstar_actual>0)then
     do i=1,npart_actual
        if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)vpart(i,3)
     enddo
  endif

  dummy=0.d0
  do i=1,ndm_actual
     write(66,*)dummy
  end do
  do i=1,nstar_actual
     write(66,*)dummy
  end do

  do i=1,nmax-nmin+1
     if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,1)
  end do
  do i=1,nmax-nmin+1
     if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,5)/(gamma-1.d0)/varp(i,1)
  end do

  do i=1,nmax-nmin+1
     write(66,*)dummy
  end do

  do i=1,nmax-nmin+1
     if(varp(i,1)>=facdens*averdens)write(66,*)varp(i,6)
  end do

  if(star.and.nstar_actual>0)then
     if(metal)then
        do i=1,npart_actual
           if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)met(i)
        end do
     else
        do i=1,nstar_actual
           write(66,*)dummy
        end do
     end if
     do i=1,npart_actual
        if(age(i)/=0.d0.and.idpart(i)>0)write(66,*)age(i)
     end do
  end if
  
  close(66)
  endif


  deallocate(xp)
  deallocate(varp)

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
       print *, 'usage: output2tipsy -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-nmi nmin] '
       print *, '                 [-nma nmax] '
       print *, '                 [-str stars] '
       print *, '                 [-snk sink] '
       print *, '                 [-met metals in stars] '
       print *, '                 [-dum ndummy (only if typ=gas)] '
       print *, '                 [-fde facdens] '
       print *, '                 [-dlm delm] '
       print *, '                 [-typ output type] '
       print *, 'ex: output2tipsy -inp output_00001 -out cube.dat'// &
            &   ' -typ "gas" -xmi 0.1 -xma 0.7 -lma 12'
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
       case ('-lma')
          read (arg,*) lmax
       case ('-nmi')
          read (arg,*) nmin
       case ('-nma')
          read (arg,*) nmax
       case ('-str') 
          read (arg,*) star
       case ('-snk') 
          read (arg,*) sink
       case ('-met') 
          read (arg,*) metal
       case ('-dum')
          read (arg,*) ndummy
       case ('-fde')
          read (arg,*) facdens
       case ('-dlm')
          read (arg,*) delm
       case ('-typ')
          read (arg,*) typ
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params

end program output2tipsy_v2