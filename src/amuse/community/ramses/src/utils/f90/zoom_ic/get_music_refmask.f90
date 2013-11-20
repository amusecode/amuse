program get_music_refmask
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite surfacique projetee
  ! des particules de matiere noire d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,nstar,nstart,inull,ico
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,ndark,ismooth
  integer::nx=0,ny=0,ix,iy,iz,ixp1,iyp1,idim,jdim,ncpu_read
  real(KIND=8)::mtot,ddx,ddy,dex,dey,t,soft,poty,mass,btime,unit_l,aexp,unit_t
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,r,xc=0.5,yc=0.5,zc=0.5,rad=-1
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin,ipart
  real(KIND=8)::xxmin,xxmax,yymin,yymax,dy,deltax,fakeage
  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=8),dimension(:,:),allocatable::map
  real(KIND=8),dimension(:)  ,allocatable::x
  real(KIND=8),dimension(:)  ,allocatable::y
  real(KIND=8),dimension(:)  ,allocatable::z
  real(KIND=8),dimension(:)  ,allocatable::vx
  real(KIND=8),dimension(:)  ,allocatable::vy
  real(KIND=8),dimension(:)  ,allocatable::vz
  real(KIND=8),dimension(:)  ,allocatable::m
  real(KIND=8),dimension(:)  ,allocatable::temp
  real(KIND=8),dimension(:)  ,allocatable::bt
  real(KIND=8),dimension(:)  ,allocatable::tempx,tempy,tempz,tempvx,tempvy,tempvz,tempm,tempbt,tempr
  integer ,allocatable,dimension(:)::tempid,temp2,indtempid
  integer ,allocatable,dimension(:)::id
  integer ,allocatable,dimension(:)::idpart,indidpart
  integer::outputmode
  real ,allocatable,dimension(:,:,:)::imark
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering,format_grille
  character(LEN=80)::GMGM
  character(LEN=128)::nomfich,repository,filetype='bin',grafic
  character(LEN=128)::nomfich2,repository2, outputname
  logical::ok,ok_part,periodic=.false.,star=.false.,okerode=.false.
  integer::impi,ndom,bit_length,maxdom,maxid,idd
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,vfact
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  integer(kind=4)::np1,np2,np3
  real::dx,dx2,x1o,x2o,x3o,astart,omegam,omegav,h0,x1or,x2or,x3or,dxor,omegak

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
  write(*,*)ncpu,ndim,levelmin,levelmax

  read(10,*)
  read(10,'("time        =",E23.15)')t
  read(10,'("aexp        =",E23.15)')aexp
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'("unit_l      =",E23.15)')unit_l
  read(10,*)
  read(10,'("unit_t      =",E23.15)')unit_t

  read(10,*)
  read(10,'("ordering type=",A80)'),ordering
  read(10,*)
  write(*,'(" ordering type=",A20)'),TRIM(ordering)
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

  if(rad>0) then
     xmin=xc-rad
     xmax=xc+rad
     ymin=yc-rad
     ymax=yc+rad
     zmin=zc-rad
     zmax=zc+rad
  endif

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

  npart=0
  do k=1,ncpu_read
     write(*,*) 'CPU=',k
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
  write(*,*) npart,' particles in the region'
  allocate(m(1:npart))
  allocate(x(1:npart))
  allocate(y(1:npart))
  allocate(z(1:npart))
  allocate(vx(1:npart))
  allocate(vy(1:npart))
  allocate(vz(1:npart))
  allocate(id(1:npart))
  if(nstar>0) then
     allocate(bt(1:npart))
  endif

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot=0.0d0
  nstart=1
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
     allocate(temp(1:npart2))
     allocate(temp2(1:npart2))
     ! Read positions
     read(1)temp
     x(nstart:nstart+npart2-1)=temp
     read(1)temp
     y(nstart:nstart+npart2-1)=temp
     read(1)temp
     z(nstart:nstart+npart2-1)=temp
     ! Read velocity
     read(1)temp
     vx(nstart:nstart+npart2-1)=temp
     read(1)temp
     vy(nstart:nstart+npart2-1)=temp
     read(1)temp
     vz(nstart:nstart+npart2-1)=temp
!     Read mass
     read(1)temp
     m(nstart:nstart+npart2-1)=temp 
     !Read identity
     read(1)temp2
     id(nstart:nstart+npart2-1)=temp2
     !Read level
     read(1)temp2   
     if(nstar>0) then
        ! Read BT
        read(1)temp
        bt(nstart:nstart+npart2-1)=temp 
     endif
! ----------------------------
     nstart=nstart+npart2  !Fill up the next set
     deallocate(temp)
     deallocate(temp2)
  enddo

40 format(3e16.8)
50 format(2I16)
  !Outputs IDs of selected particles
  ipart=0
  mass=0.0
  write(*,*) 'Getting IDs...'
  open(18,file='partID.dat',form='formatted')
  maxid=0
  do i=1,npart !To get maximum identity of the particle        
     if(nstar.eq.0) then  !Only DM particles
        btime=0
     else
        btime=bt(i)
     endif
     if(btime.eq.0) then
        ok_part=(x(i)>=xmin.and.x(i)<=xmax.and. &
             &   y(i)>=ymin.and.y(i)<=ymax.and. &
             &   z(i)>=zmin.and.z(i)<=zmax)
        if(rad>0) then
           r=(x(i)-xc)**2+(y(i)-yc)**2+(z(i)-zc)**2
           ok_part=(sqrt(r)<=rad)
        endif
        if(ok_part) then
           maxid=max(maxid,id(i))
           ipart=ipart+1
           mass=mass+m(i)
        endif
     endif
  enddo
  write(*,*) 'We have',ipart,' particles in selected region'
  write(*,*) 'Total mass =', mass
  
30 format(i16)
  write(18,50) ipart,npart,maxid

  allocate(idpart(1:ipart))
  j=1
  do i=1,npart  !Start finding the IDs
     if(nstar.eq.0) then  !Only DM particles
        btime=0
     else
        btime=bt(i)
     endif
     if(btime.eq.0) then
        ok_part=(x(i)>=xmin.and.x(i)<=xmax.and. &
             &   y(i)>=ymin.and.y(i)<=ymax.and. &
             &   z(i)>=zmin.and.z(i)<=zmax)
        if(rad>0) then
           r=(x(i)-xc)**2+(y(i)-yc)**2+(z(i)-zc)**2
           ok_part=sqrt(r)<=rad
        endif
        if(ok_part) then
           write(18,30) id(i)   !Write IDs
           idpart(j) = id(i)
           j = j+1
        endif
     endif
  enddo

  npart = ipart
  allocate(indidpart(1:npart))
  CALL quick_sort(idpart, indidpart, npart)

  !------- read IC data and match -----------

  deallocate(m)
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(id)
  if(nstar>0) then
     deallocate(bt)
  endif

  
  allocate(m(1:npart))
  allocate(x(1:npart))
  allocate(y(1:npart))
  allocate(z(1:npart))
  allocate(vx(1:npart))
  allocate(vy(1:npart))
  allocate(vz(1:npart))
  allocate(id(1:npart))
  allocate(bt(1:npart))
  

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot=0.0d0
  !nstart=1
  ipart=0

  ipos=INDEX(repository2,'output_')
  nchar=repository2(ipos+7:ipos+13)
  
  do icpu=1,ncpu
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository2)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
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
     allocate(tempx(1:npart2))
     allocate(tempy(1:npart2))
     allocate(tempz(1:npart2))
     allocate(tempvx(1:npart2))
     allocate(tempvy(1:npart2))
     allocate(tempvz(1:npart2))
     allocate(tempm(1:npart2))
     allocate(tempid(1:npart2))
     allocate(indtempid(1:npart2))
     allocate(temp2(1:npart2))
     allocate(tempbt(1:npart2))
     ! Read positions
     read(1)tempx
     read(1)tempy
     read(1)tempz
     ! Read velocity
     read(1)tempvx
     read(1)tempvy
     read(1)tempvz
     ! Read mass
     read(1)tempm
     ! Read identity
     read(1)tempid
     ! Read level
     read(1)temp2   
     if(nstar.gt.0) then
        ! Read BT
        read(1)tempbt
     else
        tempbt=0
     endif

     close(1)

     call quick_sort(tempid, indtempid, npart2)
     
     ico=1
     i=1
     do while (i.le.npart2.and.ico.le.npart)
        if(tempid(i).lt.idpart(ico))then 
           i=i+1
        else
           if(tempid(i).gt.idpart(ico))then
              ico=ico+1
           else
              if(tempid(i).eq.idpart(ico).and.tempbt(i).ne.0)then
                 i=i+1
              else
                 if(tempid(i).eq.idpart(ico).and.tempbt(i).eq.0)then
                    ipart=ipart+1
                    m(ipart)=tempm(indtempid(i))
                    x(ipart)=tempx(indtempid(i))
                    y(ipart)=tempy(indtempid(i))
                    z(ipart)=tempz(indtempid(i))
                    vx(ipart)=tempvx(indtempid(i))
                    vy(ipart)=tempvy(indtempid(i))
                    vz(ipart)=tempvz(indtempid(i))
                    id(ipart)=tempid(i)
                    bt(ipart)=tempbt(indtempid(i))
                    i=i+1
                    ico=ico+1
                 end if
              end if
           end if
        end if
     end do

     deallocate(tempx)
     deallocate(tempy)
     deallocate(tempz)
     deallocate(tempvx)
     deallocate(tempvy)
     deallocate(tempvz)
     deallocate(tempm)
     deallocate(tempid)
     deallocate(indtempid)
     deallocate(temp2)
     deallocate(tempbt)
     
  end do

  open(20,file=TRIM(outputname),form='formatted')

  if( outputmode .eq. 1 ) then
     do i=1,ipart
        write(20,1001) x(i), y(i), z(i), vx(i), vy(i), vz(i)
     end do
  else
      do i=1,ipart
        write(20,1002) x(i), y(i), z(i)
     end do
  end if

  close(20)

  write(*,*)'Wrote data to file ',trim(outputname)

1001 format (f16.8,f16.8,f16.8,e16.7,e16.7,e16.7)
1002 format (f16.8,f16.8,f16.8)

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
         print *, 'usage: geticref  -inf  input_dir_final_snapshot'
         print *, '                 -ini  input_dir_first_snapshot'
         print *, '                 [-out output_name] '
         print *, '                 [-xc xc] '
         print *, '                 [-yc yc] '
         print *, '                 [-zc zc] '
         print *, '                 [-rad rad] '
         print *, '                 [-vel 0|1] output also velocity data'
         print *, 'ex: geticref -inf output_00010 -ini output_00001 -xc 0.5 -yc 0.5 -zc 0.5 -rad 0.1'
         stop
      end if

      outputname = 'music_region_file.txt'
      outputmode = 0

      do i = 1,n,2
         call getarg(i,opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
         end if
         call getarg(i+1,arg)
         select case (opt)
         case ('-inf')
            repository = trim(arg)
         case ('-ini')
            repository2 = trim(arg)
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
         case ('-xc')
            read (arg,*) xc
         case ('-yc')
            read (arg,*) yc
         case ('-zc')
            read (arg,*) zc
         case ('-rad')
            read (arg,*) rad
         case ('-out')
            outputname = trim(arg)
         case ('-vel')
            read (arg,*) outputmode
         case ('-per')
            read (arg,*) periodic 
         case ('-gfc')
            grafic = trim(arg)
         case ('-fil')
            filetype = trim(arg)
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params
  
  end program get_music_refmask

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

SUBROUTINE quick_sort(list,order,n)
  
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  
  IMPLICIT NONE
  INTEGER :: n
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
  
  ! Local variable
  INTEGER :: i
  
  DO i = 1, n
     order(i) = i
  END DO
  
  CALL quick_sort_1(1, n)
  
CONTAINS
  
  RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER              :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6
    
    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort(left_end, right_end)
       
    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1
       
       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO
          
          
          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       
       IF (left_end < j) CALL quick_sort_1(left_end, j)
       IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
  END SUBROUTINE quick_sort_1
  
  
  SUBROUTINE interchange_sort(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER           :: temp
    
    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO
    
  END SUBROUTINE interchange_sort
  
END SUBROUTINE quick_sort

