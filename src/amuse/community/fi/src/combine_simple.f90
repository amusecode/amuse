program combsimple
 integer :: i,nbod1,nbod2,nsph1,nsph2,nstar1,nstar2,ioer
 real :: pmass,x(3),v(3)
 character*200 :: filename1,filename2,outfile
 
print*,'** combine to <simple> snapshot files **'
print*,'file 1?'
read*, filename1
print*,'file 2?'
read*,filename2
print*,'output file?'
read*,outfile

 open (1, file=filename1, form='unformatted')
 read (1,iostat=ioer) nbod1,nsph1,nstar1
   if(ioer) then 
    rewind(1)
    read (1,iostat=ioer) nbod1
    nsph1=0
    nstar1=0
   endif
 print*, 'file1 has:',nbod1,nsph1,nstar1

 open (2, file=filename2, form='unformatted')
 read (2,iostat=ioer) nbod2,nsph2,nstar2
   if(ioer) then 
    rewind(2)
    read (2,iostat=ioer) nbod2
    nsph2=0
    nstar2=0
   endif
 print*, 'file2 has:',nbod2,nsph2,nstar2

 open (3, file=outfile, form='unformatted')
 write(3), nbod1+nbod2,nsph1+nsph2,nstar1+nstar2

 do i=1,nsph1
  read(1) pmass,x,v
  write(3) pmass,x,v
 enddo
 do i=1,nsph2
  read(2) pmass,x,v
  write(3) pmass,x,v
 enddo
 do i=1,nbod1-nsph1-nstar1
  read(1) pmass,x,v
  write(3) pmass,x,v
 enddo
 do i=1,nbod2-nsph2-nstar2
  read(2) pmass,x,v
  write(3) pmass,x,v
 enddo
 do i=1,nstar1
  read(1) pmass,x,v
  write(3) pmass,x,v
 enddo
 do i=1,nstar2
  read(2) pmass,x,v
  write(3) pmass,x,v
 enddo

 close(1)
 close(2)
 close(3)

end program
