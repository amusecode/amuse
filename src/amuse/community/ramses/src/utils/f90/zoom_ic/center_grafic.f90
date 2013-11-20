program center_grafic
  !---------------------------------------------------------------------
  ! Ce programme recentre autour d'un pixel donne les fichiers ic_... 
  ! generes par GRAFIC.
  ! Ce programme doit lire en input les fichiers suivants:
  !          - un fichier deltab: input/ic_deltab
  !          - un fichier velbx:  input/ic_velbx
  !          - un fichier velby:  input/ic_velby
  !          - un fichier velbz:  input/ic_velbz
  !          - un fichier velbx:  input/ic_velcx
  !          - un fichier velby:  input/ic_velcy
  !          - un fichier velbz:  input/ic_velcz
  !          - un fichier refmap: input/ic_refmap
  !          - un fichier metal:  input/ic_pvar_00001
  ! Il genere en output les fichiers suivants:
  !          - un fichier deltab: output/ic_deltab
  !          - un fichier velbx:  output/ic_velbx
  !          - un fichier velby:  output/ic_velby
  !          - un fichier velbz:  output/ic_velbz
  !          - un fichier velcx:  output/ic_velcx
  !          - un fichier velcy:  output/ic_velcy
  !          - un fichier velcz:  output/ic_velcz
  !          - un fichier refmap: output/ic_refmap
  !          - un fichier metal:  output/ic_pvar_00001
  !                    
  !         
  ! M. Gonzalez
  ! Saclay, le 31/08/01.
  !---------------------------------------------------------------------
  !  f90 center_grafic.f90 -o ~/bin/center_grafic
  !---------------------------------------------------------------------
  implicit none
  integer::xc1,xc2,xc3,i1,i2,i3,np1,np2,np3,skip1,skip2
  integer::min_x,max_x,min_y,max_y,min_z,max_z
  integer::i,j,k,i_file,narg,iargc
  real::x1o,x2o,x3o,dx,astart,omegam,omegav,h0
  real,dimension(:,:),allocatable::f1,f2
  character*80::input,output
  character*80,dimension(18)::filename
  logical::ok,found

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

  ! COMPUTE FILES TO OPEN AND TO WRITE 
  filename(1) =TRIM(input)//'/ic_deltab'
  filename(2) =TRIM(input)//'/ic_velcx'
  filename(3) =TRIM(input)//'/ic_velcy'
  filename(4) =TRIM(input)//'/ic_velcz'
  filename(5) =TRIM(input)//'/ic_velbx'
  filename(6) =TRIM(input)//'/ic_velby'
  filename(7) =TRIM(input)//'/ic_velbz'
  filename(8) =TRIM(input)//'/ic_refmap'
  filename(9) =TRIM(input)//'/ic_pvar_00001'

  filename(10) =TRIM(output)//'/ic_deltab'
  filename(11) =TRIM(output)//'/ic_velcx'
  filename(12)=TRIM(output)//'/ic_velcy'
  filename(13)=TRIM(output)//'/ic_velcz'
  filename(14)=TRIM(output)//'/ic_velbx'
  filename(15)=TRIM(output)//'/ic_velby'
  filename(16)=TRIM(output)//'/ic_velbz'
  filename(17)=TRIM(output)//'/ic_refmap'
  filename(18)=TRIM(output)//'/ic_pvar_00001'

  ! GET INPUT FILE PARAMETERS
  found=.false.
  i_file=0
  do while(.not.found)
     i_file=i_file+1

     INQUIRE(file=filename(i_file),exist=ok)
     if(ok)then
        print*,'CENT GRAF opening file:',filename(i_file)
        open(11,file=filename(i_file),form='unformatted')
        read(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
        close(11)
        found=.true.
     endif
  enddo
  write(*,*)'Input array size is          :',np1,np2,np3
  write(*,*)'Old center coordinates are   :',np1/2,np2/2,np3/2
  write(*,*)'Enter new center coordinates : (i,j,k)'
  read(*,*) xc1,xc2,xc3

  min_x=xc1-np1/2
  max_x=xc1+np1/2
  min_y=xc2-np2/2
  max_y=xc2+np2/2
  min_z=xc3-np3/2
  max_z=xc3+np3/2

  skip1=max(0,min_z)
  skip2=min(max_z,np3)

  allocate(f1(np1,np2))
  allocate(f2(np1,np2))
     
  do i_file=1,9

     INQUIRE(file=filename(i_file),exist=ok)

     if(ok) then

        ! READING INPUT FILES
        write(*,*)'Reading input file '//TRIM(filename(i_file))
        open(11,file=TRIM(filename(i_file)),form='unformatted')
        read(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
        write(*,*)'Writing output file '//TRIM(filename(9+i_file))
        open(12,file=TRIM(filename(9+i_file)),form='unformatted')
        write(12) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
        
        do i3=1,skip2
           read(11)
        end do
        do i3=skip2+1,np3
           read (11)((f1(i1,i2),i1=1,np1),i2=1,np2)
           do i=min_x+1,max_x
              i1=i
              if(i1 < 1) i1=i1+np1
              if(i1 > np1) i1=i1-np1
              do j=min_y+1,max_y
                 i2=j
                 if(i2 < 1) i2=i2+np2
                 if(i2 > np2) i2=i2-np2
                 f2(i-min_x,j-min_y)=f1(i1,i2)
              enddo
           enddo
           write(12)((f2(i1,i2),i1=1,np1),i2=1,np2)
        end do
        close(11)
        
        open(11,file=TRIM(filename(i_file)),form='unformatted')
        read(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
        do i3=1,skip1
           read(11)
        end do
        do i3=skip1+1,skip2
           read (11)((f1(i1,i2),i1=1,np1),i2=1,np2)
           do i=min_x+1,max_x
              i1=i
              if(i1 < 1) i1=i1+np1
              if(i1 > np1) i1=i1-np1
              do j=min_y+1,max_y
                 i2=j
                 if(i2 < 1) i2=i2+np2
                 if(i2 > np2) i2=i2-np2
                 f2(i-min_x,j-min_y)=f1(i1,i2)
              enddo
           enddo
           write(12)((f2(i1,i2),i1=1,np1),i2=1,np2)
        end do
        close(11)
        
        open(11,file=TRIM(filename(i_file)),form='unformatted')
        read(11) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
        do i3=1,skip1
           read (11)((f1(i1,i2),i1=1,np1),i2=1,np2)
           do i=min_x+1,max_x
              i1=i
              if(i1 < 1) i1=i1+np1
              if(i1 > np1) i1=i1-np1
              do j=min_y+1,max_y
                 i2=j
                 if(i2 < 1) i2=i2+np2
                 if(i2 > np2) i2=i2-np2
                 f2(i-min_x,j-min_y)=f1(i1,i2)
              enddo
           enddo
           write(12)((f2(i1,i2),i1=1,np1),i2=1,np2)
        end do
        do i3=skip1+1,np3
           read(11)
        end do
        
        close(11)

        write(*,*)'done'
        close(12)

     endif

  enddo

  deallocate(f1,f2)
     

end program center_grafic
