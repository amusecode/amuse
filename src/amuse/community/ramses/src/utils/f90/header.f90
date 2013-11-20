program header
  !---------------------------------------------------------------------
  ! Ce programme lit le header d'un fichier amr.out d'une simulation RAMSES
  ! Version F90 par R. Teyssier le 28/02/00
  ! Version parallele par R. Teyssier le 13/06/03
  ! f90 head_amr.f90 -o ~/bin/header
  !---------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,n
  integer::nx,ny,nz
  integer::nlevelmax
  integer::ngridmax,nstep_coarse
  real::t,aexp,hexp,boxlen
  real::omega_m,omega_l,omega_k,omega_b
  real::scale_l,scale_d,scale_t
  logical::ok
  character*5::nchar
  character*50::nomfich,nomdir
  integer::iargc,ipos
  character(len=128)::arg

  n = iargc()
  if (n.NE.1) then
     print *, 'usage: header output_dir'
     stop
  end if
  call getarg(1,arg)
  nomdir=trim(arg)//'/'
  !-----------------------------------------------
  ! Lecture des fichiers outputs au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(nomdir,'output_')
  nchar=nomdir(ipos+7:ipos+13)
  ! Lecture du fichier AMR
  nomfich=TRIM(nomdir)//'amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(nomfich)//' not found'
  else
     open(unit=10,file=nomfich,status='old',form='unformatted')
     read(10)ncpu
     read(10)ndim
     read(10)nx,ny,nz
     read(10)nlevelmax
     read(10)ngridmax
     read(10)nstep_coarse    

     read(10)boxlen
     read(10)t,aexp,hexp
     read(10)omega_m,omega_l,omega_k,omega_b
     read(10)scale_l,scale_d,scale_t
     close(10)

     write(*,991)ncpu,ndim,nstep_coarse
     write(*,993)nlevelmax,ngridmax
     write(*,997)boxlen
     write(*,994)t,aexp,hexp
     write(*,995)omega_m,omega_l,omega_k,omega_b
     write(*,996)scale_l,scale_d,scale_t
  end if
     
990 format(' Enter output number:',i6)
991 format(' ncpu=',i6,' ndim=',i1,' nstep=',i6)
992 format(' nx=',i3,' ny=',i3,' nz=',i3)
993 format(' nlevelmax=',i3,' ngridmax=',i8)
994 format(' t=',1pe10.3,' aexp=',1pe10.3,' hexp=',1pe10.3)
995 format(' omega_m=',F6.3,' omega_l=',F6.3,' omega_k=',F6.3,' omega_b=',F6.3)
996 format(' scale_l=',E9.3,' scale_d=',E9.3,' scale_t=',E9.3)
997 format(' boxlen=',1pe10.3)

end program header

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

