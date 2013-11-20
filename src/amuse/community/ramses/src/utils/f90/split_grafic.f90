!================================================================
!================================================================
program spg
  implicit none
  !===============================================================
  ! This code splits grafic initial conditions files in ncpu files
  ! following a 3D space filling curve.
  ! 
  ! Syntax : split_grafic input output ncpu
  ! 
  ! - input is a directory which contains the 7 grafic files
  ! - output will contain the splitted files
  ! - ncpu must be a power of 2
  ! - a verbose mode may be activated by adding 1 at the end of the 
  ! command line
  !
  ! v 0.1 (Dominique Aubert/SAP) initial version
  !===============================================================
  integer::np1,np2,np3
  real::dx,x1off,x2off,x3off
  real::astart,omegam,omegav,h0
  real::zstart
  character(len=80)::filegraf,fname,fdir,input,output,ncpu_string,commande,fcurr
  character(len=5)::extnum
  character(len=80)::debug_string
  character(len=80),dimension(14)::filename
  logical::verbose=.false.,yorick=.false.
  integer,dimension(:),allocatable::open_flag,open_line
  integer::ncoarse,ncpu,debug

  real(kind=8),dimension(0:2048)::bound_key
  real(kind=8)::dx_loc,x_plan,y_plan,z_plan
  real(kind=8),dimension(:,:),allocatable::plan
  real(kind=8),dimension(1,3)::x

  real(kind=8)::boxlen,dl
  integer::nlevelmax,icoarse_min,icoarse_max,nx_loc
  integer::jcoarse_min,jcoarse_max,kcoarse_min,kcoarse_max,ind
  real(kind=8),dimension(:),allocatable::order,order_min,order_max,order_plan,o
  real(kind=8)::order_all_min,order_all_max
  
  integer,dimension(:),allocatable::cpu_map,cpu_plan
  integer,dimension(:,:),allocatable::cpu_plan_min,cpu_plan_max

  integer::i1,i2,i3,il,ix,iy,iz,i,minunit,ifile,c

  real(kind=4)::x1o_loc,x2o_loc,x3o_loc
  integer::np1_loc,np2_loc,np3_loc

  integer::nplan,info
  integer::narg
  integer::iargc
  real,dimension(:,:),allocatable::ics
  integer::ncode,bit_length

  common /graficparam/np1,np2,np3,dx,x1off,x2off,x3off,astart,omegam,omegav,h0
  common /ioparam/filegraf
  common /simuparam/boxlen,nlevelmax,icoarse_min,icoarse_max
  common /cpumap/ncpu
  common /ordcom/order_all_min,order_all_max
  common /bound/bound_key
  common /hilbertparam1/nx_loc,ncode,bit_length

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------

  !! Command line analysis

  narg=iargc()
  IF(narg .LT. 3)THEN
     write(*,*)'You should type: split_grafic input output ncpu'
     write(*,*)'where directory input should contain GRAFIC files'
     write(*,*)'and directory output should be empty'
     write(*,*)'ncpu MUST BE a power of 2'
     STOP
  END IF

  CALL getarg(1,input)
  CALL getarg(2,output)
  CALL getarg(3,ncpu_string)

  if(narg.eq.4) then
     CALL getarg(4,debug_string)
     read(debug_string,*) debug
     if(debug.eq.1) then 
        verbose=.true.
     elseif(debug.eq.2) then
        verbose=.true.
        yorick=.true.
     end if
  end if

  !!Getting the number of cpus  
  read(ncpu_string,*) ncpu
  dl=log(real(ncpu))/log(2.)
  if(ncpu.eq.1.or.(dl-int(dl)).ne.0) then
     write(*,*) 'ERROR: ncpu should be a power of 2'
     stop
  end if


  !! Reading the Grafic Header  
  filegraf=trim(input)//'/ic_deltab'
  call read_grafic_header

  
  !! Defining some constants
  zstart=1./astart-1.
  boxlen=np1*dx

  ncoarse=np1*np2*np3

  icoarse_min=1
  icoarse_max=np1
  
  jcoarse_min=icoarse_min
  jcoarse_max=icoarse_max

  kcoarse_min=icoarse_min
  kcoarse_max=icoarse_max

  ! Hilbert curve domain decomposition
  nx_loc=icoarse_max-icoarse_min+1
  ncode=2*nx_loc
  bit_length=int(log(dble(ncode))/log(2.))+1

  if(verbose) then
     write(*,*) '*************************'
     write(*,*) 'Cosmology found in grafic files'
     write(*,*) 'astart=',astart
     write(*,*) 'zstart=',zstart
     write(*,*) 'omegam=',omegam
     write(*,*) 'omegav=',omegav
     write(*,*) 'h0=',h0
     write(*,*) 'dx (Mpc)=',dx
     write(*,*) 'Boxlen (Mpc)=',boxlen
     write(*,*) '*************************'
  end if
  
  !! allocation of CPU boundaries
  allocate(cpu_plan(1:ncpu))
  allocate(cpu_plan_min(2,1:ncpu))
  allocate(cpu_plan_max(2,1:ncpu))
  allocate(open_flag(ncpu))
  allocate(open_line(ncpu))
  allocate(ics(np1,np2))

  !! Get Min and Max order
  dx_loc=boxlen
  x(1,:)=0.5*boxlen
  order_all_min=dble(0d0)
  order_all_max=dble(2d0*np1)**3
  if(verbose) write(*,*) 'Minimal Hilbert order =',order_all_min
  if(verbose) write(*,*) 'Maximal Hilbert order =',order_all_max

  !!Bound_keys
  do i=0,ncpu-1
     bound_key(i)=order_all_min+dble(i)/dble(ncpu)*(order_all_max-order_all_min)
  end do
  bound_key(ncpu)=order_all_max

  !! **********
  !! Z Climbing
  !! **********

  ! init files and directories names
  filename(1) = '/ic_deltab'
  filename(2) = '/ic_velcx'
  filename(3) = '/ic_velcy'
  filename(4) = '/ic_velcz'
  filename(5) = '/ic_velbx'
  filename(6) = '/ic_velby'
  filename(7) = '/ic_velbz'

  filename(8) = '/dir_deltab'
  filename(9) = '/dir_velcx'
  filename(10)= '/dir_velcy'
  filename(11)= '/dir_velcz'
  filename(12)= '/dir_velbx'
  filename(13)= '/dir_velby'
  filename(14)= '/dir_velbz'

  if(verbose) write(*,*) 'Z Splitting over ',ncpu,' cpus : Start'

  ! Main loop over grafic files
  do ifile=1,4

     write(*,*) 'Splitting '//filename(ifile)

     !! Creating directory structure
     fdir=trim(output)//filename(ifile+7)
     commande='mkdir '//fdir
!     call system(commande)
     call PXFMKDIR(TRIM(fdir),LEN(TRIM(fdir)),O'755',info)
     
     !opening ICs files
     open(unit=10,file=trim(input)//filename(ifile),form='unformatted')
     read(10) np1,np2,np3,dx,x1off,x2off,x3off,astart,omegam,omegav,h0
     
     ! some inits
     minunit=20
     open_flag(:)=-1
     
     ! main Loop over Grafic planes
     do iz=kcoarse_min,kcoarse_max
        
        if(verbose) write(*,*) 'Parsing plane #',iz
        z_plan=dble(iz+0.5-kcoarse_min)*2d0
        
        ! some inits in the plane
        cpu_plan(:)=-1
        cpu_plan_min(:,:)=-1
        cpu_plan_max(:,:)=-1

        ! looking for cpus in current plane
        do iy=jcoarse_min,jcoarse_max
        y_plan=dble(iy+0.5-jcoarse_min)*2d0
           do ix=icoarse_min,icoarse_max
           x_plan=dble(ix+0.5-icoarse_min)*2d0
           
           call cmp_cpumap(x_plan,y_plan,z_plan,c)
           
           ! The cpu c has been found
           cpu_plan(c)=1 
           
           ! Looking for Rectangles corners
           if(cpu_plan_min(1,c).eq.-1) then
              cpu_plan_min(1,c)=ix
              cpu_plan_min(2,c)=iy
           end if
           
           if(ix.gt.cpu_plan_max(1,c))cpu_plan_max(1,c)=ix
           if(iy.gt.cpu_plan_max(2,c))cpu_plan_max(2,c)=iy
           
           end do
        end do

        ! reading data
        read(10) ((ics(i1,i2),i1=1,np1),i2=1,np2)

        !loop over cpus
        do i=1,ncpu
           
           ! have the cpu been found in this plane ?
           if(cpu_plan(i).ne.-1) then
              
              ! is the file already open ? if no we create it
              if(open_flag(i).eq.-1) then
                 
                 open_flag(i)=minunit
                 call title(i,extnum)

                 fcurr=trim(fdir)//filename(ifile)
                 fcurr=trim(fcurr)//'.'
                 fname=trim(fcurr)//trim(extnum)
                 
                 if(verbose) write(*,*) '************  opening file '//fname
                 open(unit=minunit,file=fname,form='unformatted')
                 minunit=minunit+1
                 
                 ! we should write the header
                 np1_loc=cpu_plan_max(1,i)-cpu_plan_min(1,i)+1
                 np2_loc=cpu_plan_max(2,i)-cpu_plan_min(2,i)+1
                 np3_loc=ncoarse/ncpu/np1_loc/np2_loc
                 
                 x1o_loc=real((cpu_plan_min(1,i)-icoarse_min)*dx)
                 x2o_loc=real((cpu_plan_min(2,i)-jcoarse_min)*dx)
                 x3o_loc=real((iz               -kcoarse_min)*dx)
                 
                 write(open_flag(i)) np1_loc,np2_loc,np3_loc,dx,x1o_loc,x2o_loc,x3o_loc,astart,omegam,omegav,h0
                 
                 ! we write the first line of data
                 write(open_flag(i)) ((ics(i1,i2),i1=cpu_plan_min(1,i),cpu_plan_max(1,i)),i2=cpu_plan_min(2,i),cpu_plan_max(2,i))
                 
                 open_line(i)=np3_loc-1
                                  
              else
                 
                 ! we write the current line of data
                 write(open_flag(i)) ((ics(i1,i2),i1=cpu_plan_min(1,i),cpu_plan_max(1,i)),i2=cpu_plan_min(2,i),cpu_plan_max(2,i))
                 open_line(i)=open_line(i)-1
                 
                 !we check if the current file is finished and should be closed
                 if(open_line(i).eq.0) then
                    if(verbose) write(*,*) '*********** closing cpu file #',i
                    close(open_flag(i))
                 end if
                 
              end if
           end if
        end do
        !Loop end over the z-planes
     end do
  
     close(10)
     ! Loop end over the current grafic file
  end do

  if(verbose) write(*,*) 'Z Climbing : Done'
  
end program spg

!================================================================
!================================================================
!===========SUBROUTINES SECTION ================================
!================================================================
!================================================================
subroutine cmp_cpumap(x,y,z,c)
  implicit none
  integer::c
  real(kind=8)::x,y,z
  integer::ncpu,icpu
  real(kind=8)::order
  real(kind=8),dimension(0:2048)::bound_key
  common /cpumap/ncpu
  common /bound/bound_key 
  call cmp_ordering(x,y,z,order)
  c=ncpu ! default value
  do icpu=1,ncpu
     if(    order.ge.bound_key(icpu-1).and. &
          & order.lt.bound_key(icpu  ))then
        c=icpu
     endif
  end do
end subroutine cmp_cpumap
!================================================================
!================================================================
subroutine cmp_ordering(x,y,z,order)
  implicit none
  real(kind=8)::x,y,z,order
  !--------------------------------------------------------
  ! This routine computes the index key of the input cell
  ! according to its position in space and for the chosen
  ! ordering. Position x are in user units.
  !-----------------------------------------------------
  integer::i,ncode,bit_length,nx_loc
  integer,dimension(1)::ix,iy,iz
  real(kind=8),dimension(1)::xorder
  common /hilbertparam1/nx_loc,ncode,bit_length
  ix(1)=int(x)
  iy(1)=int(y)
  iz(1)=int(z)
  call hilbert3d(ix,iy,iz,xorder,bit_length,1)
  order=xorder(1)
end subroutine cmp_ordering
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

!================================================================
!================================================================
subroutine read_grafic_header
  implicit none
  
  integer::np1,np2,np3
  real::dx,x1off,x2off,x3off
  real::astart,omegam,omegav,h0
  character(len=80)::filegraf

  common /graficparam/np1,np2,np3,dx,x1off,x2off,x3off,astart,omegam,omegav,h0
  common /ioparam/filegraf
  

  open(unit=10,file=trim(filegraf),form='unformatted')
  read(10) np1,np2,np3,dx,x1off,x2off,x3off,astart,omegam,omegav,h0
  close(10)

end subroutine read_grafic_header
!================================================================
!================================================================
!================================================================
!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character(LEN=5)::nchar

  character(LEN=1)::nchar1
  character(LEN=2)::nchar2
  character(LEN=3)::nchar3
  character(LEN=4)::nchar4
  character(LEN=5)::nchar5

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

