subroutine basicinput(filenaam)
  include 'globals.h'
  character*200 filenaam
  integer i,mapnr,nfeedback
  real,allocatable :: weights(:)
 
!  call set_parameters(0)

  print*,'filenaam?'
  read(*,'(a)') filenaam
 
  call readbods(filenaam)
  
!  call heattabel
         
!  call initpars

!  if(periodic) call initbc

!  call initnewstar

!  if(usepm) then
!    if(verbosity.GT.0) print*,' ...initPM...'
!    call initpm
!  endif

!  call postprocessread

end subroutine

program snapreader
  include 'globals.h'
  integer p
  character*200 filenaam
 
  call initmem(nbodsmax,nsphmax,ncells)
  call basicinput(filenaam)

! insert code after this

  open(unit=1,file=trim(filenaam)//'.simple',form='unformatted')
  write(1) nbodies,nsph,nstar
  do p=1,nbodies
    write(1), mass(p),pos(p,1:3),vel(p,1:3)
  enddo
  close(1) 
end program
