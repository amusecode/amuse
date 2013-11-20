program center_whitenoise
  !---------------------------------------------------------------------
  ! Ce programme recentre autour d'un pixel donne les fichiers ic_... 
  ! generes par GRAFIC.
  ! Ce programme doit lire en input les fichiers suivants:
  !          - un fichier white noise
  ! Il genere en output les fichiers suivants:
  !          - un fichier white noise
  !             
  ! M. Gonzalez
  ! Saclay, le 31/08/01.
  !---------------------------------------------------------------------
  !  f90 center_grafic.f90 -o ~/bin/center_grafic
  !---------------------------------------------------------------------
  implicit none
  integer::xc1,xc2,xc3,i1,i2,i3,np1,np2,np3,skip1,skip2
  integer::min_x,max_x,min_y,max_y,min_z,max_z
  integer::i,j,k,i_file,narg,iargc,iseed
  real::x1o,x2o,x3o,dx,astart,omegam,omegav,h0
  real,dimension(:,:),allocatable::f1,f2
  character*80::input,output

  narg = iargc()
  IF(narg .NE. 2)THEN
     write(*,*)'You should type: a.out input output'
     write(*,*)'where directory input should be a grafics white noise file'
     STOP
  END IF

  CALL getarg(1,input)
  CALL getarg(2,output)

  !  SAFETY CONDITION
  if (input == output) then 
     write(*,*)'If input and output files are the same'
     write(*,*)'input will be erased'
     write(*,*)'so type DIFFERENT directories !!!!'
     stop
  endif

  ! GET INPUT FILE PARAMETERS
  open(11,file=trim(input),form='unformatted')
  read(11) np1,np2,np3,iseed
  close(11)

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
     
     ! READING INPUT FILE
  write(*,*)'Reading input file '//TRIM(input)
  open(11,file=input,form='unformatted')
  read(11)  np1,np2,np3,iseed
  write(*,*)'Writing output file '//TRIM(output)
  open(12,file=output,form='unformatted')
  write(12) np1,np2,np3,iseed

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
  
  open(11,file=input,form='unformatted')
  read(11)  np1,np2,np3,iseed
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
  
  open(11,file=input,form='unformatted')
  read(11) np1,np2,np3,iseed
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
    
  deallocate(f1,f2)
     

end program center_whitenoise
