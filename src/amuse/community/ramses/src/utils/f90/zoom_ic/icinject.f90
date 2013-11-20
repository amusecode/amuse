program icdegrade
  implicit none

  integer::i1,i2,i3,i,j,k,narg,iargc,i_file,ii,jj,kk,xstart,ystart,zstart
  integer(kind=4)::np1,np2,np3
  integer(kind=4)::hnp1,hnp2,hnp3
  real::dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  real::hdx,hx1o,hx2o,hx3o
  real,dimension(:,:,:),allocatable::f,f2
  character*80::input,output
  character*80,dimension(3)::filename
  logical::ok

  narg = iargc()
  IF(narg .NE. 2)THEN
     write(*,*)'You should type: a.out input1 input2'
     write(*,*)'where directory input1 should contain the ic_refmap file to be injected'
     write(*,*)'and directory input2 should contain an ic_deltab file of the output size'
     write(*,*)'and dimensions. The resulting ic_refmap will be put in the input2 directory.'
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

  ! READING INPUT FILES
  filename(1) =TRIM(input)//'/ic_refmap'
  filename(2) =TRIM(output)//'/ic_deltab'
  filename(3) =TRIM(output)//'/ic_refmap'
  INQUIRE(file=filename(1),exist=ok)
  if(ok)then
  
     open(10,file=filename(1),form='unformatted')
     read (10)np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
     write(*,*)'ic_refmap array size is:',np1,np2,np3
     allocate(f(np1,np2,np3))
     do k=1,np3
        read(10) ((f(i,j,k),i=1,np1),j=1,np2)
     enddo
     close(10)
     
     open(11,file=filename(2),form='unformatted')
     read (11)hnp1,hnp2,hnp3,hdx,hx1o,hx2o,hx3o,astart,omegam,omegav,h0
     write(*,*)'Host array size is:',hnp1,hnp2,hnp3
     allocate(f2(hnp1,hnp2,hnp3))
     close(11)
     f2=0
     write(*,*) 'ic_refmap box size=', np1*dx
     write(*,*) 'Host box size=', hnp1*hdx
     write(*,*)'Injecting ic_refmap into host array...'
     if(dx.ne.hdx) then
        write(*,*) 'Error! Cell sizes are different!'
        write(*,*) 'ic_refmap dx=',dx
        write(*,*) 'Host dx=',hdx
        stop
     endif
     xstart=(x1o-hx1o)/dx
     ystart=(x2o-hx2o)/dx
     zstart=(x3o-hx3o)/dx
     
     write(*,*) 'Start cells are {i,j,k}=',xstart,ystart,zstart
     
     ! Loop over planes
     kk=1
     do k=zstart+1,zstart+np3
        jj=1
        do j=ystart+1,ystart+np2
           ii=1
           do i=xstart+1,xstart+np1
              f2(i,j,k)=f(ii,jj,kk)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
     
     write(*,*) 'Outputting Grafics file ic_refmap'
     open(33,file=filename(3),form='unformatted')
     write(33) hnp1,hnp2,hnp3,hdx,hx1o,hx2o,hx3o,astart,omegam,omegav,h0 
     do k=1,hnp3
        write(33) ((f2(i,j,k),i=1,hnp1),j=1,hnp2)
     enddo
     close(33)
     deallocate(f,f2)

  endif

  filename(1) =TRIM(input)//'/ic_pvar_00001'
  filename(2) =TRIM(output)//'/ic_deltab'
  filename(3) =TRIM(output)//'/ic_pvar_00001'
  INQUIRE(file=filename(1),exist=ok)
  if(ok)then
  
     open(10,file=filename(1),form='unformatted')
     read (10)np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
     write(*,*)'ic_pvar_00001 array size is:',np1,np2,np3
     allocate(f(np1,np2,np3))
     do k=1,np3
        read(10) ((f(i,j,k),i=1,np1),j=1,np2)
     enddo
     close(10)
     
     open(11,file=filename(2),form='unformatted')
     read (11)hnp1,hnp2,hnp3,hdx,hx1o,hx2o,hx3o,astart,omegam,omegav,h0
     write(*,*)'Host array size is:',hnp1,hnp2,hnp3
     allocate(f2(hnp1,hnp2,hnp3))
     close(11)
     f2=0
     write(*,*) 'ic_pvar_00001 box size=', np1*dx
     write(*,*) 'Host box size=', hnp1*hdx
     write(*,*)'Injecting ic_pvar_00001 into host array...'
     if(dx.ne.hdx) then
        write(*,*) 'Error! Cell sizes are different!'
        write(*,*) 'ic_refmap dx=',dx
        write(*,*) 'Host dx=',hdx
        stop
     endif
     xstart=(x1o-hx1o)/dx
     ystart=(x2o-hx2o)/dx
     zstart=(x3o-hx3o)/dx
     
     write(*,*) 'Start cells are {i,j,k}=',xstart,ystart,zstart
     
     ! Loop over planes
     kk=1
     do k=zstart+1,zstart+np3
        jj=1
        do j=ystart+1,ystart+np2
           ii=1
           do i=xstart+1,xstart+np1
              f2(i,j,k)=f(ii,jj,kk)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
     
     write(*,*) 'Outputting Grafics file ic_pvar_00001'
     open(33,file=filename(3),form='unformatted')
     write(33) hnp1,hnp2,hnp3,hdx,hx1o,hx2o,hx3o,astart,omegam,omegav,h0 
     do k=1,hnp3
        write(33) ((f2(i,j,k),i=1,hnp1),j=1,hnp2)
     enddo
     close(33)
     
  endif

end program icdegrade


       
