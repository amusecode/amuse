program extract_grafic
  !---------------------------------------------------------------------
  ! Ce programme extrait un sous-cube des fichiers ic_... generes
  !  par GRAFIC.
  ! Ce programme doit lire en input les fichiers suivants:
  !          - un fichier deltab: input/ic_deltab
  !          - un fichier velbx:  input/ic_velbx
  !          - un fichier velby:  input/ic_velby
  !          - un fichier velbz:  input/ic_velbz
  !          - un fichier velbx:  input/ic_velcx
  !          - un fichier velby:  input/ic_velcy
  !          - un fichier velbz:  input/ic_velcz
  ! Il genere en output les fichiers suivants:
  !          - un fichier deltab: output/ic_deltab
  !          - un fichier velbx:  output/ic_velbx
  !          - un fichier velby:  output/ic_velby
  !          - un fichier velbz:  output/ic_velbz
  !          - un fichier velcx:  output/ic_velcx
  !          - un fichier velcy:  output/ic_velcy
  !          - un fichier velcz:  output/ic_velcz
  !                    
  !         
  ! M. Gonzalez
  ! Saclay, le 31/08/01.
  !---------------------------------------------------------------------
  !  f90 extract_grafic.f90 -o ~/bin/extract_grafic
  !---------------------------------------------------------------------
  implicit none
  integer::xc1,xc2,xc3,i1,i2,i3,np1,np2,np3,iargc,narg,nfiles
  integer::np1_cube,np2_cube,np3_cube
  integer::min_x,max_x,min_y,max_y,min_z,max_z
  integer::i,j,k,i_file
  real::x1o,x2o,x3o,x1o_cube,x2o_cube,x3o_cube,dx,astart,omegam,omegav,h0
  real,dimension(:,:),allocatable::f,f_cube
  character*80::input,output
  character*80,dimension(14)::filename 
  logical::cosmo_ics=.false.

  narg = iargc()
  IF(narg .NE. 2)THEN
     write(*,*)'You should type: a.out input output'
     write(*,*)'where directory input should contain GRAFIC files'
     write(*,*)'and directory output should be empty'
     STOP
  END IF

  CALL getarg(1,input)
  CALL getarg(2,output)

  !  SAFETY CONDITION
  if (input == output) then 
     write(*,*)'If input and output directories are the same'
     write(*,*)'input files will be erased by output ones'
     write(*,*)'so type DIFFERENT directories !!!!'
     stop
  endif


  !check wether its cosmoligical ics or not 
  INQUIRE(FILE=TRIM(input)//'/ic_deltab', EXIST=cosmo_ics)
  if(cosmo_ics)then
     write(*,*)'looking for ic_deltab, ic_velcx,...'
     ! READING INPUT FILES   
     filename(1) =TRIM(input)//'/ic_deltab'
     filename(2) =TRIM(input)//'/ic_velcx'
     filename(3) =TRIM(input)//'/ic_velcy'
     filename(4) =TRIM(input)//'/ic_velcz'
     filename(5) =TRIM(input)//'/ic_velbx'
     filename(6) =TRIM(input)//'/ic_velby'
     filename(7) =TRIM(input)//'/ic_velbz'
     filename(8) =TRIM(output)//'/ic_deltab'
     filename(9) =TRIM(output)//'/ic_velcx'
     filename(10)=TRIM(output)//'/ic_velcy'
     filename(11)=TRIM(output)//'/ic_velcz'
     filename(12)=TRIM(output)//'/ic_velbx'
     filename(13)=TRIM(output)//'/ic_velby'
     filename(14)=TRIM(output)//'/ic_velbz'
  else
     write(*,*)'looking for ic_d, ic_u,...'
  ! READING INPUT FILES 
     filename(1) =TRIM(input)//'/ic_d'
     filename(2) =TRIM(input)//'/ic_u'
     filename(3) =TRIM(input)//'/ic_v'
     filename(4) =TRIM(input)//'/ic_w'
     filename(5) =TRIM(input)//'/ic_p'
     filename(6) =''
     filename(7) =''
     filename(8) =TRIM(output)//'/ic_d'
     filename(9) =TRIM(output)//'/ic_u'
     filename(10)=TRIM(output)//'/ic_v'
     filename(11)=TRIM(output)//'/ic_w'
     filename(12)=TRIM(output)//'/ic_p'
     filename(13)=''
     filename(14)=''
  end if

  open(11,file=filename(1),form='unformatted')
  read(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  close(11)
  write(*,*)'Input array size is :',np1,np2,np3
  write(*,*)'Enter center of cube xc1,xc2,xc3 (input mesh units)'
  read(*,*) xc1,xc2,xc3
  write(*,*)'Enter length of cube nc1,nc2,nc3 (input mesh units)'
  read(*,*) np1_cube,np2_cube,np3_cube
  write(*,*) np1_cube,np2_cube,np3_cube
  
  min_x=max(xc1-np1_cube/2,0)
  max_x=min(xc1+np1_cube/2,np1)
  min_y=max(xc2-np2_cube/2,0)
  max_y=min(xc2+np2_cube/2,np2)
  min_z=max(xc3-np3_cube/2,0)
  max_z=min(xc3+np3_cube/2,np3)
  np1_cube=max_x-min_x
  np2_cube=max_y-min_y
  np3_cube=max_z-min_z
  write(*,*) np1_cube,np2_cube,np3_cube

  !  COMPUTING NEW OFFSETS
  x1o_cube=x1o+min_x*dx
  x2o_cube=x2o+min_y*dx
  x3o_cube=x3o+min_z*dx

  allocate(f(np1,np2))
  allocate(f_cube(np1_cube,np2_cube))
  
  nfiles=5
  if(cosmo_ics)nfiles=4
  do i_file=1,nfiles

     write(*,*)'Reading input file '//TRIM(filename(i_file))
     open(11,file=filename(i_file),form='unformatted')
     read(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
     
     write(*,*)'Writing ouput file '//TRIM(filename(7+i_file))
     open(12,file=filename(7+i_file),form='unformatted')
     write(12) np1_cube,np2_cube,np3_cube,dx,x1o_cube,x2o_cube,x3o_cube,astart,omegam,omegav,h0

     do i3=1,min_z
        read(11)
     end do
     do i3=min_z+1,max_z
        read (11)((f(i1,i2),i1=1,np1),i2=1,np2)
        do i1=min_x+1,max_x
           do i2=min_y+1,max_y
              f_cube(i1-min_x,i2-min_y)=f(i1,i2)
           enddo
        enddo
        write(12)((f_cube(i1,i2),i1=1,np1_cube),i2=1,np2_cube)
     end do
     do i3=max_z+1,np3
        read(11)
     end do

     close(11)
     close(12)
     
  enddo

  deallocate(f,f_cube)

end program extract_grafic
