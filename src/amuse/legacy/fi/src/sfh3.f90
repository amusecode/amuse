program sfrinfo
 real  :: tmin,tmax,dt,rmax,dr
 integer :: ntime,nr
 integer :: ig,is,i,j
 real :: pos(3),vel(3),rho,temp,xe,fuv,h2,smass,gmass,time
 character :: filename*80
 real, allocatable :: sfh(:,:),taxis(:),raxis(:)
 real,parameter :: pi=3.141592654
 integer :: naxes(3)
 
 print*,' sf history file?'
 read*, filename
 print*,' tmin,tmax,rmax?'
 read*, tmin,tmax,rmax
 print*,'# time,r bins?'
 read*, ntime,nr
 dt=(tmax-tmin)/ntime
 dr=2*rmax/nr
 allocate(sfh(ntime,nr),taxis(ntime),raxis(nr))
 
 do i=1,ntime
  taxis(i)=tmin+dt*(i-.5)
 enddo
 
 do i=1,nr
  raxis(i)=-rmax+dr*(i-.5)
 enddo
  sfh=0.
 open(unit=20,file=trim(filename)//'.sfh',status='OLD',IOSTAT=ioerror)
 if(ioerror.NE.0) then
  print*,' file error'
  stop
 endif
 
10 read(20,*,IOSTAT=ioerror) ig,is
   if(ioerror.EQ.0) then
    read(20,*) time
    read(20,*) pos
    read(20,*) vel
    read(20,*) rho
    read(20,*) temp
    read(20,*) xe
    read(20,*) fuv
    read(20,*) h2
    read(20,*) smass,gmass
    
    i=(time-tmin)/dt
    j=(pos(1)+rmax)/dr    
    if(i.ge.1.and.i.le.ntime.and.j.ge.1.and.j.le.nr) sfh(i,j)=sfh(i,j)+smass
    
    goto 10
   endif
  close(20)

  do i=1,nr
   do j=1,ntime
    sfh(j,i)=sfh(j,i)/dt
   enddo
  enddo

  filename=trim(filename)//'_sfh.fits'
  
  naxes(1)=ntime
  naxes(2)=nr
  call writefits(filename,2,naxes,sfh)
  
  deallocate(sfh,taxis)

end program

 subroutine writefits(filename,naxis,naxes,pic)
      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer i,j,k,group,fpixel,nelements
      character filename*80
      logical simple,extend
      real pic(*)

      status=0
      print*,filename
      call deletefile(filename,status)
      call ftgiou(unit,status)

      blocksize=1
      call ftinit(unit,filename,blocksize,status)

      simple=.true.
      bitpix=-64
      extend=.true.

      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      group=1
      fpixel=1
      nelements=1
      do i=1,naxis
       nelements=naxes(i)*nelements
      enddo
      call ftpprd(unit,group,fpixel,nelements,pic,status)

      call ftclos(unit, status)
      call ftfiou(unit, status)
      
      if (status .gt. 0)call printerror(status)
      end subroutine
          
      subroutine printerror(status)


      integer status
      character errtext*30,errmessage*80

      if (status .le. 0)return
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      enddo
      end subroutine
      
      
      subroutine deletefile(filename,status)
      integer status,unit,blocksize
      character*(*) filename

      if (status .gt. 0)return

      call ftgiou(unit,status)

      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
          call ftdelt(unit,status)
      else if (status .eq. 103)then
          status=0
          call ftcmsg
      else
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

      call ftfiou(unit, status)
      end subroutine

