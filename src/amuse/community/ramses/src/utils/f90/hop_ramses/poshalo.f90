program poshalo
  implicit none
  integer::ncpu,ndim,npart,i,j,icpu,ipos,nstar,ngroup
  integer::ncpu2,npart2,ndim2,id_group
  integer,dimension(:),allocatable::group_id,n_group
  real(kind=8),dimension(:),allocatable::xref_group,yref_group,zref_group
  logical,dimension(:),allocatable::flag_group
  real(KIND=8)::mtot,mcut=1000000,contamine
  real(KIND=8)::period=1,xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,rvir
  real(kind=8),dimension(:),allocatable::m_group,x_group,y_group,z_group,u_group,v_group,w_group,mpure_group
  real(KIND=8),dimension(:,:),allocatable::x,v
  real(KIND=8),dimension(:),allocatable::m,age
  integer,dimension(:),allocatable::id
  character(LEN=5)::nstring,ncharcpu
  character(LEN=80)::directory,file_groupe_in
  character(LEN=128)::nomfich,repository
  logical::ok

  call read_params

  ! Reading groups
  file_groupe_in=TRIM(directory)//'.tag'
  Open(15,file=trim(file_groupe_in),FORM='UNFORMATTED')
  Read(15)npart, ngroup
  Write(*,*)'npart=',npart
  Write(*,*)'ngroup=',ngroup
  allocate(group_id(1:npart))
  allocate(n_group(1:ngroup),m_group(1:ngroup),mpure_group(1:ngroup))
  allocate(flag_group(1:ngroup))
  allocate(xref_group(1:ngroup),yref_group(1:ngroup),zref_group(1:ngroup))
  allocate(x_group(1:ngroup),y_group(1:ngroup),z_group(1:ngroup))
  allocate(u_group(1:ngroup),v_group(1:ngroup),w_group(1:ngroup))
  Read(15)(group_id(j),j=1,npart)
  Close(15)
  x_group=0.0
  y_group=0.0
  z_group=0.0
  u_group=0.0
  v_group=0.0
  w_group=0.0
  m_group=0.0
  mpure_group=0.0
  n_group=0.0
  flag_group=.true.  ! True if group reference particle is UNSET
  xref_group=0.0
  yref_group=0.0
  zref_group=0.0
  group_id=group_id+1  ! WARNING: convert from 0 to ng-1 to 1 to ng

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nstring=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part_'//TRIM(nstring)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nstring)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,'(13X,I11)')ncpu
  read(10,'(13X,I11)')ndim
  close(10)

  npart=0
  do icpu=1,ncpu
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nstring)//'.out'//TRIM(ncharcpu)
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
  if(nstar>0)then
     write(*,*)'Discard star particles.'
  endif

  mtot=0.0d0
  npart=0
  do icpu=1,ncpu
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nstring)//'.out'//TRIM(ncharcpu)
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
     if(nstar>0)then
        allocate(age(1:npart2))
        allocate(id(1:npart2))
     endif
     allocate(x(1:npart2,1:ndim2))
     allocate(v(1:npart2,1:ndim2))
     ! Read position
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m
     end do
     ! Read velocity
     do i=1,ndim
        read(1)m
        v(1:npart2,i)=m
     end do
     ! Read mass
     read(1)m
     if(nstar>0)then
        read(1)id
        read(1) ! Skip level
        read(1)age
     endif
     close(1)

     do i=1,npart2
        if(nstar>0) then
           if(age(i)/=0.or.id(i)<1) cycle
        end if
        npart=npart+1
        id_group=group_id(npart)
        if(id_group>0.and.id_group.LE.ngroup)then
           if(flag_group(id_group)) then
              flag_group(id_group)=.false.
              xref_group(id_group)=x(i,1)
              yref_group(id_group)=x(i,2)
              zref_group(id_group)=x(i,3)
           end if
           n_group(id_group)=n_group(id_group)+1
           m_group(id_group)=m_group(id_group)+m(i)
           if(m(i)<=mcut)mpure_group(id_group)=mpure_group(id_group)+m(i)
           x_group(id_group)=x_group(id_group)+m(i)*shift1d(x(i,1),xref_group(id_group),period)
           y_group(id_group)=y_group(id_group)+m(i)*shift1d(x(i,2),yref_group(id_group),period)
           z_group(id_group)=z_group(id_group)+m(i)*shift1d(x(i,3),zref_group(id_group),period)
           u_group(id_group)=u_group(id_group)+m(i)*v(i,1)
           v_group(id_group)=v_group(id_group)+m(i)*v(i,2)
           w_group(id_group)=w_group(id_group)+m(i)*v(i,3)
        endif
     enddo
     deallocate(x,m,v)
     if(nstar>0)deallocate(age,id)
  end do

  open(18,file=TRIM(directory)//'.pos')
  write(*,'(A100)')'   #   npart       mass  cont.frac         xc         yc         zc         uc         vc         wc'  
  write(18,'(A100)')'   #   npart       mass  cont.frac         xc         yc         zc         uc         vc         wc'  
  do i=1,ngroup
     x_group(i)=shift1d(x_group(i)/m_group(i),0.5d0,period) ! Ensure group center is within bounds
     y_group(i)=shift1d(y_group(i)/m_group(i),0.5d0,period)
     z_group(i)=shift1d(z_group(i)/m_group(i),0.5d0,period)
     u_group(i)=u_group(i)/m_group(i)
     v_group(i)=v_group(i)/m_group(i)
     w_group(i)=w_group(i)/m_group(i)
     if(         x_group(i).gt.xmin.and.x_group(i).lt.xmax &
          & .and.y_group(i).gt.ymin.and.y_group(i).lt.ymax &
          & .and.z_group(i).gt.zmin.and.z_group(i).lt.zmax)then
        rvir=(m_group(i)/200./(4./3.*3.1415926))**(1./3.)
        contamine=(m_group(i)-mpure_group(i))/m_group(i)
        if(i<=10)write(*,'(I5,A,I7,A,1PE10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3)')&
             & i,' ',n_group(i),' ',m_group(i),' ',contamine,' ',x_group(i),' ',y_group(i),' ',z_group(i),' ', &
             & u_group(i),' ',v_group(i),' ',w_group(i)
                write(18,'(I5,A,I7,A,1PE10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3)') &
             & i,' ',n_group(i),' ',m_group(i),' ',contamine,' ',x_group(i),' ',y_group(i),' ',z_group(i),' ', &
             & u_group(i),' ',v_group(i),' ',w_group(i)
     endif
  end do
  close(18)

contains

  real(kind=8) function shift1d(x, xref, per)
      !> Shifts x by a multiple of per such that abs(x-xref) is minimal
      real(kind=8) :: x, xref, per
      integer :: ishift
      if(per>0d0) then
         ishift = floor((xref-x)/per+0.5)
         shift1d = x+ishift*per
      else
         shift1d = x
      end if
  end function shift1d

  subroutine read_params

      implicit none

      integer       :: i,n
      integer       :: iargc
      character(len=4)   :: opt
      character(len=128) :: arg
      
      n = iargc()
      if (n < 4) then
         print *, 'usage: poshalo -inp ramses_input -pre hop_prefix'
         print *, 'ex: poshalo -inp output_00001 -pre zregroup'
         stop
      end if

      do i = 1,n,2
         call getarg(i,opt)
         if (i == n) then
            print '("option ",a4," has no argument")', opt
            stop 2
         end if
         call getarg(i+1,arg)
         select case (opt)
         case ('-pre')
            directory = trim(arg)
         case ('-inp')
            repository = trim(arg)
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
         case ('-cut')
            read (arg,*) mcut
         case ('-zma')
            read (arg,*) zmax
         case ('-per')
            read (arg,*) period
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params
end program poshalo
  

subroutine title(n,nstring)
  implicit none
  integer::n
  character*5::nstring

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nstring = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nstring = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nstring = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nstring = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nstring = '0000'//nchar1
  endif
end subroutine title

