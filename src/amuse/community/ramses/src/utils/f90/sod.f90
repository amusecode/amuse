program sphericaloverdensity
  !---------------------------------------------------------------------
  ! Ce programme calcule la fonction de masse avec
  ! l'algorithme SPHERICAL OVERDENSITY.
  ! Ce programme doit lire en input les fichiers suivants:
  !          - un fichier resultat RAMSES:    part_00600.out
  ! Il genere en output les fichiers suivants:
  !          - un fichier fonction de masse:  mf_sod_00600.dat
  !          - un fichier structure:          st_sod_00600.dat
  ! E. Audit, R. Teyssier
  ! Meudon, le 30/05/96.
  !---------------------------------------------------------------------
  ! Version F90 par R. Teyssier le 28/02/00   
  ! f90 sod.f90 -o sod
  ! Memory: 44*npart in bytes
  !         npart=256**3   707 Mo
  !         npart=512**3   5.6 Go
  !---------------------------------------------------------------------
  implicit none
  integer::ndim,npart,ngrid,n,i,icpu,ipos
  integer::ny,nz,ncpu,npart_new
  integer::ncpu2,npart2,ndim2
  real(kind=8)::r,scale,Mmin,mtot,rhomoyen
  real,dimension(:,:),allocatable::x
  real,dimension(:),allocatable::mp,tmp
  integer,dimension(:),allocatable::isort
  real(kind=8),dimension(:),allocatable::xsort

  real(kind=8)::seuil=600.
  integer::nmin=100,nx=0,ncut=4

  character*5::nchar,ncharcpu
  character*50::nomfich,repository
  logical::ok

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

  open(unit=1,file=nomfich,status='old',form='unformatted')
  read(1)ncpu
  read(1)ndim
  close(1)

  npart=0
  do icpu=1,ncpu
  call title(icpu,ncharcpu)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
  write(*,*)'Reading file '//TRIM(nomfich)
  open(unit=1,file=nomfich,status='old',form='unformatted')
  read(1)ncpu2
  read(1)ndim2
  read(1)npart2
  close(1)
  npart=npart+npart2
  end do

  write(*,*)'Found ',npart,' particles'
  write(*,*)'Reading positions and masses...'
  allocate(x(1:npart,1:ndim))
  allocate(mp(1:npart),xsort(1:npart))
  allocate(isort(1:npart))

  npart=0
  do icpu=1,ncpu
  call title(icpu,ncharcpu)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
  open(unit=1,file=nomfich,status='old',form='unformatted')
  read(1)ncpu2
  read(1)ndim2
  read(1)npart2
  allocate(tmp(1:npart2))
  do i=1,ndim
     read(1)tmp
     x(npart+1:npart+npart2,i)=tmp
  end do
  do i=1,ndim
     read(1)tmp
  end do
  read(1)tmp
  mp(npart+1:npart+npart2)=tmp
  close(1)
  npart=npart+npart2
  deallocate(tmp)
  end do

  mtot=0.0d0
  do i=1,npart
     mtot=mtot+mp(i)
  enddo
  if(nx==0)then
     nx=int(dble(npart)**(1./3.))
  endif

  ny=nx
  nz=nx
  rhomoyen=mtot
  write(*,*)'Working mesh=',nx,ny,nz
  write(*,*)'Found rho_bar=',rhomoyen,' -> rescaling to 1.0'
  mp=mp/rhomoyen

  scale=nx
  x=x*scale
  mp=mp*scale**3
  Mmin=dble(nmin)/dble(npart)*scale**3

  ngrid=nx*ny*nz

  !------------------------------------------
  ! Tri des particules selon leur densite NGP
  !------------------------------------------
  write(*,*) 'Sorting particles by NGP density'
  call tri_ngp(x,mp,xsort,isort,npart,nx,ny,nz,npart_new)

  !---------------------------------------------------------
  ! Calcul de la fonction de masse par SPHERICAL OVERDENSITY
  !---------------------------------------------------------
  write(*,*) 'Computing halo mass function'
  write(*,*) 'using a spherical overdensity delta_bar=',seuil
  call spherover(x,mp,xsort,isort,nchar,npart,nx,ny,nz,npart_new,ncut,seuil,Mmin,scale)
  write(*,*) 'done'


  stop

contains
  subroutine read_params

      implicit none

      integer       :: i,n
      integer       :: iargc
      character(len=4)   :: opt
      character(len=128) :: arg
      LOGICAL       :: bad, ok
      
      n = iargc()
      if (n < 2) then
         print *, 'usage: sod  [-inp input_dir]'
         print *, '            [-nx  nx_grid]  (optional)'
         print *, '            [-min np_min]   (optional)'
         print *, '            [-sel ncell]    (optional)'
         print *, '            [-del del_bar]  (optional)'
         print *, '            [-hlp]          (optional)'
         stop
      end if

      do i = 1,n,2
         call getarg(i,opt)
         if (opt == '-hlp') then
            print 1, repository,nx,seuil,mmin
            stop
         end if
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
         end if
         call getarg(i+1,arg)
         select case (opt)
         case ('-inp')
            repository = trim(arg)
         case ('-nx')
            read (arg,*) nx
         case('-min')
            read (arg,*) nmin
         case('-sel')
            read (arg,*) ncut
         case ('-del') 
            read (arg,*) seuil
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

1     format(/,&
           & " -inp [] ",A,/,&
           & " -nx  [] ",I6,/,&
           & " -min [] ",I6,/,&
           & " -sel [] ",I6,/,&
           & " -del [] ",e10.2,/)

      return

    end subroutine read_params

end program sphericaloverdensity

!------------------------------------------------------------
! TRI DES PARTICULES SELON LA DENSITE NGP SUR UNE GRILLE FINE
!------------------------------------------------------------
subroutine tri_ngp(x,mp,pdens,isort,npart,nx,ny,nz,npart_new)
  implicit none
  integer::npart,nx,ny,nz,npart_new
  integer,dimension(1:npart)::isort
  real,dimension(1:npart,1:3)::x
  real,dimension(1:npart)::mp
  real(kind=8),dimension(1:npart)::pdens

  integer(kind=8),dimension(1:npart)::indx
  integer::ngxngy,ngx,ngy,ngz
  integer(kind=8)::nmax,i1,i2,i3,ind
  integer::i,j,ifail
  integer::imin,imax,ntot
  real(kind=8)::amass
  
  ngx=4*nx
  ngy=4*ny
  ngz=4*nz

  ngxngy=ngx*ngy
   
  ! Remplissage du tableau des indices
  !-----------------------------------
  do i = 1,npart
     i1 = int(x(i,1)*dble(ngx)/dble(nx))
     i2 = int(x(i,2)*dble(ngy)/dble(ny))
     i3 = int(x(i,3)*dble(ngz)/dble(nz))
     indx (i) = 1+i1+i2*ngx+i3*ngxngy
  end do
  
  ! Tri des particules selon les indices
  !-------------------------------------
  ifail=0
  pdens=indx
  call quick_sort(pdens,isort,npart)
  indx=pdens

  ! Calcul de l'inverse de la densite NGP pour chaque particule
  !------------------------------------------------------------
  imin = 1
  do while(imin .le. npart)
     ind  = indx(imin)
     i    = imin
     do while(indx(i).eq.ind.and.i.lt.npart)
        i = i + 1
     end do
     if(i.eq.npart)i=npart+1
     imax = i - 1
     ntot = i - imin
     amass= 0.0
     do i = imin ,imax
       amass=amass+mp(isort(i))
     end do
     do i = imin ,imax
        pdens(isort(i)) = 1.0d0/amass
     end do
     imin = imax + 1
  end do
  print *, 1./minval(pdens)
  print *, 1./maxval(pdens)

  ! Tri des particules selon la densite
  !------------------------------------
  call quick_sort(pdens,isort,npart)

  i=1
  do while (pdens(i) < 0.5)
     i=i+1
  end do
  npart_new=i
  write(*,*)'npart_active=',npart_new
  return

end subroutine tri_ngp

!---------------------------------------------------------------------------
subroutine spherover(x,mp,xsort,isort,nchar,npart,nx,ny,nz,npart_new,ncut,seuil,Mmin,scale)
!---------------------------------------------------------------------------
  implicit none

  character*5::nchar
  integer::npart,nx,ny,nz,npart_new
  real(kind=8)::seuil,Mmin,scale
  integer::ncut          ! Amplitude de selection des particules
  real,dimension(1:npart,1:3)::x
  real,dimension(1:npart)::mp
  real(kind=8),dimension(1:npart)::xsort
  integer,dimension(1:npart)::isort
  
  integer::np
  real(kind=8)::twopi=6.283185307179586
  real(kind=8)::dtmax=0.01 ! dr/rayon pour la convergence du barycentre
  integer::Nsel=4000000    ! Nombre maximum de particules selectionnees

!  integer,dimension(1:npart)::structure
  integer,dimension(1:npart)::indx_ngp,isort_ngp
  integer,dimension(:),allocatable::first_ngp,num_ngp

  integer,dimension(:),allocatable::indx_sel,isort_sel,irank_sel
  real,dimension(:),allocatable::m_sel
  real,dimension(:,:),allocatable::x_sel
  real(kind=8),dimension(:),allocatable::distance

  integer::i,j,k,npart_sel,ipart  
  integer::is,iamas,masse,ip,ifail
  integer::ic,jc,kc,ii,jj,kk,ind
  integer::i1,i2,i3,imin,imax,ntot
  integer::niter,itest
  real(kind=8)::xc,yc,zc,xx,yy,zz,xc0,yc0,zc0,dc0
  real(kind=8)::xb,yb,zb,dt
  real(kind=8)::taillex,tailley,taillez,amass
  real(kind=8)::dx,dy,dz,d,d2,rayon,rayon2,volume,overdens
  real(kind=8)::dmin,x_tmp,y_tmp,z_tmp,dis_tmp,qqq
  real(kind=8)::Mtotsel

  character*80::nomf

!  structure = 0
  taillex   = dble(nx)
  tailley   = dble(ny)
  taillez   = dble(nz)
  nomf      = 'sod_'//TRIM(nchar)//'.dat'

  allocate(num_ngp  (1:nx*ny*nz))
  allocate(first_ngp(1:nx*ny*nz))

  allocate(indx_sel (Nsel))
  allocate(isort_sel(Nsel))
  allocate(distance (Nsel))
  allocate(x_sel    (Nsel,1:3))
  allocate(m_sel    (Nsel))

  ! Calcul du pointeur cellules vers particules
  !--------------------------------------------
  write(*,*)'Building particle linked list...'
  do i  = 1,npart
     i1 = int(x(i,1))
     i2 = int(x(i,2))
     i3 = int(x(i,3))
     indx_ngp(i) = 1+i1+i2*nx+i3*nx*ny
  end do

  xsort=indx_ngp
  call quick_sort(xsort,isort_ngp,npart)
  indx_ngp=xsort

  first_ngp(:)=0
  num_ngp(:)=0
  imin = 1
  do while(imin.le.npart)
     ind  = indx_ngp(imin)
     i    = imin
     first_ngp(ind) = imin
     do while(indx_ngp(i).eq.ind.and.i.lt.npart)
        i = i + 1
     end do
     if(i.eq.npart)i=npart+1
     imax = i - 1
     ntot = i - imin
     num_ngp(ind) = ntot
     imin = imax + 1
  end do

  write(*,*) '  halo#    npart  mass       x          y          z         it  radius     del'
  open(10,file=nomf,form='formatted')
  iamas = 0
  xsort = 1.0
  do ipart = 1,npart_new

     ! On choisit la particule la plus dense de la liste
     !-------------------------------------------------
     if(xsort(isort(ipart)).gt.0.)then

        ! On initialise le barycentre de l'amas avec 
        ! la position de la particule courante
        xc = x(isort(ipart),1)
        yc = x(isort(ipart),2)
        zc = x(isort(ipart),3)

        xc0=xc; yc0=yc; zc0=zc

        ! Preselection des particules les plus proches du barycentre
        !-----------------------------------------------------------
        is = 0
        ic = int(xc)
        jc = int(yc)
        kc = int(zc)
        mtotsel=0.
        do i=-ncut,ncut
           ii=i+ic
           if(ii.ge.nx)ii=ii-nx
           if(ii.lt.0) ii=ii+nx
           do j=-ncut,ncut
              jj=j+jc
              if(jj.ge.ny)jj=jj-ny
              if(jj.lt.0) jj=jj+ny
              do k=-ncut,ncut
                 kk=k+kc
                 if(kk.ge.nz)kk=kk-nz
                 if(kk.lt.0) kk=kk+nz
                 ind = 1+ii+jj*nx+kk*nx*ny
                 imin = first_ngp(ind)
                 imax = imin+num_ngp(ind)-1
                 do ip = imin,imax
                    if(xsort(isort_ngp(ip)).gt.0.)then
                       is = is + 1
                       if(is > Nsel)then
                          write(*,*)'Increase Nsel in main program'
                          stop
                       end if
                       indx_sel(is) =  isort_ngp(ip)
                       x_sel(is,1) = x(isort_ngp(ip),1)
                       x_sel(is,2) = x(isort_ngp(ip),2)
                       x_sel(is,3) = x(isort_ngp(ip),3)
                       m_sel(is) = mp (isort_ngp(ip))
                       mtotsel = mtotsel + m_sel(is)
                    end if
                 end do
              end do
           end do
        end do
        npart_sel = is
!        write(*,*)'npart_sel=',npart_sel
!        write(*,*)'Mtotsel=',Mtotsel,mmin

        if(Mtotsel .ge. Mmin) then

        dt    = 1.
        niter = 0

        do while(dt.gt.dtmax.and.niter.le.100)

           dx = abs(xc-xc0)
           dx = min(dx,taillex-dx)
           dy = abs(yc-yc0)
           dy = min(dy,tailley-dy)
           dz = abs(zc-zc0)
           dz = min(dz,taillez-dz)
           dc0 = dx*dx+dy*dy+dz*dz

           ! Calcul des distances des particules au barycentre
           !--------------------------------------------------
           do i = 1,npart_sel
              xx = x_sel(i,1)
              yy = x_sel(i,2)
              zz = x_sel(i,3)
              dx = abs(xx-xc)
              dx = min(dx,taillex-dx)
              dy = abs(yy-yc)
              dy = min(dy,tailley-dy)
              dz = abs(zz-zc)
              dz = min(dz,taillez-dz)
              distance(i) = dx*dx+dy*dy+dz*dz
           end do

           ! Tri des distances
           !---------------------------------------------------
           call quick_sort(distance,isort_sel,npart_sel)

           ! Calcul du rayon de Viriel et de la masse de l'amas
           !---------------------------------------------------
           i        = 1
           overdens = 2.0*seuil
           amass    = m_sel(isort_sel(i))
           rayon    = 0.0
           do while(i.lt.10.or.(overdens.gt.seuil.and.i.lt.npart_sel.and.rayon.lt.(real(ncut)-dc0)))
              i        = i+1
              rayon    = SQRT(distance(i))
              volume   = 2./3.*twopi*rayon**3
              amass    = amass+m_sel(isort_sel(i))
              overdens = amass/volume
!              write(*,*)i,rayon,amass,overdens
           end do
           if(overdens.gt.seuil)then
              write(*,*)'Increase selection radius (-sel option in command line)'
              write(*,*)'mass     =',i
              write(*,*)'nsel     =',npart_sel
              write(*,*)'overdens =',overdens
              write(*,*)'radius   =',rayon+dc0
              write(*,*)'nsel     =',ncut
              write(*,*)'niter    =',niter
              stop
           else
              rayon =SQRT(distance(i-1))
              volume=2./3.*twopi*rayon**3
              amass =amass-m_sel(isort_sel(i))
              if(volume > 0.)then
                 overdens=amass/volume
              else
                 overdens=0.
              end if
              masse=i-1
           end if
!           write(*,*)'iter=',niter,'mass=',masse

           ! Calcul du nouveau barycentre:
           !------------------------------
           xb=0.
           yb=0.
           zb=0.
           do i=1,masse
              xx = x_sel(isort_sel(i),1)
              yy = x_sel(isort_sel(i),2)
              zz = x_sel(isort_sel(i),3)
              dx = (xx-xc)
              if(dx.gt.  taillex/2. )dx = dx-taillex
              if(dx.lt.(-taillex/2.))dx = dx+taillex
              dy = (yy-yc)
              if(dy.gt.  tailley/2. )dy = dy-tailley
              if(dy.lt.(-tailley/2.))dy = dy+tailley
              dz = (zz-zc)
              if(dz.gt.  taillez/2. )dz = dz-taillez
              if(dz.lt.(-taillez/2.))dz = dz+taillez
              xb = xb+dx*m_sel(isort_sel(i))
              yb = yb+dy*m_sel(isort_sel(i))
              zb = zb+dz*m_sel(isort_sel(i))
           enddo
           xb = xb/amass
           yb = yb/amass
           zb = zb/amass
           if(rayon > 0)then
              dt = sqrt(xb*xb+yb*yb+zb*zb)/rayon
           else
              dt = 0.
           end if
           xc = xc+xb
           yc = yc+yb
           zc = zc+zb

           niter = niter + 1
		 
           ! On sort de la boucle si la position du barycentre converge
           !-----------------------------------------------------------
        enddo

        if(amass.ge.Mmin)then
           iamas = iamas + 1
           write(*,997)iamas,masse,amass/scale**3,xc/scale,yc/scale,zc/scale &
                & ,niter,rayon/scale,overdens
        end if

        ! On retire de la liste des particules celles de l'amas
        !------------------------------------------------------
        do i = 1,masse
!           structure(indx_sel(isort_sel(i))) = iamas
           xsort(indx_sel(isort_sel(i))) = 0.
        enddo
	      
        ! Ecriture des caracteristiques de l'amas sur fichier
        !----------------------------------------------------
        if(amass.ge.Mmin)then
           write(10,999)iamas,masse,amass/scale**3,xc/scale,yc/scale,zc/scale
        end if

     end if

     end if
	   
     if(mod(ipart,100000)==0)write(*,'(0PF5.1,"% complete")')100.*dble(ipart)/dble(npart_new)

  end do
  
  close(10)

!  nomf = 'st_sod_'//TRIM(nchar)//'.dat'
!  open (10,file=nomf,form='unformatted')
!  write(10)structure
!  close(10)

997 format (I8,1x,I8,4(1x,1pe10.3),1x,I2,1x,5(E10.3,1x))
999 format (I8,1x,I8,4(1x,1pe10.3))
	
  return
end subroutine spherover

subroutine title(n,nchar)
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

SUBROUTINE quick_sort(list, order, n)
  
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  
  IMPLICIT NONE
  INTEGER :: n
  REAL*8, DIMENSION (1:n), INTENT(INOUT)  :: list
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
    REAL*8              :: reference, temp
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
    REAL*8              :: temp
    
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
