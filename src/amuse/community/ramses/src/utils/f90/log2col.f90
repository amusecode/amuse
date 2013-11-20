program log2col
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite projetee pour les
  ! variables hydro d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  integer::nstep,ipos,jpos
  logical::ok
  character(LEN=4)::char4
  character(LEN=5)::nchar,ncharcpu
  character(LEN=6)::char6
  character(LEN=9)::char9
  character(LEN=11)::char11
  character(LEN=12)::char12
  character(LEN=128)::nomfich,repository,outfich
  character(LEN=128)::full_line,char
  real::ekin,t,a,dt,epot,econs,mem_grid,mem_part

  call read_params

  !-----------------------------------------------
  ! Lecture du fichier log au format RAMSES
  !-----------------------------------------------  
  nomfich=TRIM(repository)
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  print *,'Reading file '//TRIM(nomfich)
  open(unit=10,file=nomfich,status='old',form='formatted')
  
  !-----------------------------------------------
  ! Ecriture du fichier log en colonnes
  !-----------------------------------------------  
  open(unit=11,file=TRIM(outfich),form='formatted')
  write(11,'("   nstep            t          dt        aexp        ekin        epot       econs  memgrid(%)  mempart(%) ")')

  do

     read(10,'(a)',END=999)full_line
     ipos=index(full_line,"Main step")
     if(ipos.gt.0)then
        char6=full_line(ipos+11:ipos+16)
        read(char6,'(I6)')nstep
        jpos=index(full_line,"econs=")
        char9=full_line(jpos+6:jpos+15)
        read(char9,'(1pe9.2)')econs
        jpos=index(full_line,"epot=")
        char9=full_line(jpos+5:jpos+14)
        read(char9,'(1pe9.2)')epot
        jpos=index(full_line,"ekin=")
        char11=full_line(jpos+5:jpos+14)
        read(char11,'(1pe9.2)')ekin
        read(10,'(a)',END=999)full_line
        jpos=index(full_line,"t=")
        char12=full_line(jpos+2:jpos+14)
        read(char12,'(1pe12.5)')t
        jpos=index(full_line,"dt=")
        char11=full_line(jpos+3:jpos+14)
        read(char11,'(1pe11.4)')dt
        jpos=index(full_line,"a=")
        char11=full_line(jpos+2:jpos+13)
        read(char11,'(1pe11.4)')a
        jpos=index(full_line,"mem=")
        char4=full_line(jpos+4:jpos+8)
        read(char4,'(F4.1)')mem_grid
        char4=full_line(jpos+10:jpos+14)
        read(char4,'(F4.1)')mem_part

        write(11,'(I8,1X,1pe12.5,7(1X,1pe11.4))')nstep,t,dt,a,ekin,epot,econs,mem_grid,mem_part

     endif

  enddo
999  print *,'Complete file '//TRIM(outfich)

  close(10)
  close(11)
  
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
       print *, 'usage: log2col -inp amrlog_file'
       print *, '               -out column_file'
       print *, 'ex: log2col -inp run.log -out run.col'
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
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    return
    
  end subroutine read_params
  
end program log2col

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

