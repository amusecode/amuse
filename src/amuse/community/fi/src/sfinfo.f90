program sfrinfo
 integer :: ig,is,i,ioerror
 real :: pos(3),vel(3),rho,temp,xe,fuv,h2,smass,gmass,time
 character :: filename*80
 
 print*,' sf history file?'
 read*, filename
 
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
    
    write(*,5),time,rho,h2,xe
    
    goto 10
   endif
  close(20)
5 format( e12.4,' ',e12.4,' ',e12.4,' ',e12.4)


end program
