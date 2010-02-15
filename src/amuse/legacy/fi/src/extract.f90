module extractMod

 implicit none
 private
 public :: initextract,extractdensity,endextract

 integer, parameter :: dim=3
 real :: vertex(dim+1,dim)
 integer :: n(0:dim)
 real, allocatable, save :: pic(:,:,:)

 contains

subroutine initextract
 integer :: i,j

 n(0:dim)=0
 vertex(1:dim+1,1:dim)=0
 do i=1,dim
  print*,'# punten dim ',i; read*,n(i)
  if(n(i).LE.1) exit
 enddo  

 do i=1,dim
  n(i)=max(n(i),1)
 enddo 
 
 allocate(pic(n(1),n(2),n(3)))
 
 if(n(1).EQ.0) then
   print*,' extractdens heeft niks te doen!'
   return
 endif  
 
 do i=1,dim+1
  if(n(i-1).ge.2.OR.i.EQ.1) then 
    print*,'punt',i,'coordinaten?'
    read*,(vertex(i,j),j=1,dim)
    if(i.gt.1) vertex(i,1:dim)=vertex(i,1:dim)-vertex(1,1:dim)
  endif
 enddo
 
end subroutine

subroutine extractdensity(weights)
 
 real :: x(dim),t(dim),hsearch,dens,ddensdh,hn,weights(*),meana
 integer i,j,k,tot,nn,index(dim)
 character,dimension(dim) :: axes='xyz'

 hsearch=0.
 t(1:dim)=0  
 tot=1;   
 do i=1,dim
 tot=tot*n(i)
 enddo 
 
  
 do i=0,tot-1
  k=i
  x(1:dim)=vertex(1,1:dim)
  do j=1,dim
   if(n(j).EQ.1) then
    index(j:dim)=1
    exit
   endif
   t(j)=MOD(k,n(j))/REAL(n(j)-1)
   index(j)=MOD(k,n(j))+1
   x(1:dim)=x(1:dim)+t(j)*vertex(j+1,1:dim)
   k=k/n(j)   
  enddo
  nn=0
  hn=0
  call hsmdenspos(x(1:3),hsearch,dens,ddensdh,nn,hn)
  call gatter3(nn,x(1:3),hsearch,dens,ddensdh,meana,weights)
  pic(index(1),index(2),index(3))=meana
!  print*,'extractdensity:',(x(j),j=1,dim),gatscatdensity(x(1:3))
 enddo

end subroutine extractdensity

subroutine endextract(filename)
  character, optional :: filename*80
  character :: file*80
  integer :: naxes(3)
 
  if(present(filename)) then
   file=filename
  else
   print*,' uit filenaam?'
   read*, file
  endif
 
  naxes(1)=n(1)
  naxes(2)=n(2)
  naxes(3)=n(3)
  call writefits(file,3,naxes)
 
  deallocate(pic)
  
end subroutine 

 subroutine writefits(filename,naxis,naxes)
      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer i,j,k,group,fpixel,nelements
      character filename*80
      logical simple,extend

      status=0

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



end module

program extract
 use extractMod
 include 'globals.h'
 character*24 filenaam
 integer i,mapnr,nfeedback
 real,allocatable :: weights(:)
 call initmem(nbodsmax,nsphmax,ncells)
 
 call set_parameters(0)

 print*,'filenaam?'
 read*,filenaam

 call readbods(filenaam)

 call heattabel
        
 call initpars

 if(periodic) call initbc

 call initnewstar

 if(usepm) then
  if(verbosity.GT.0) print*,' ...initPM...'
  call initpm
 endif

 call postprocessread

10  print*,' map of? (0=rho, 1=HI, 2=H_2, 3=xe, 4= A, 5= A+Asn)'
    read*,mapnr

 nfeedback=0
 if(mapnr.EQ.5) then
  do i=nbodies-nstar+1,nbodies
   if(snentropy(i).GT.0) then
    nfeedback=nfeedback+1
    bodlist(nfeedback)=i
   endif
  enddo 
 endif     
 
 allocate(weights(nsph+nfeedback))

 weights(1:nsph)=mass(1:nsph)

 if(mapnr.EQ.1) then
  do i=1,nsph
   weights(i)=mass(i)*MAX(0.,(1-h2frac(i)-elecfrac(i)))
  enddo
 endif
 
 if(mapnr.EQ.2) then
  do i=1,nsph
   weights(i)=mass(i)*h2frac(i)
  enddo
 endif
 
 if(mapnr.EQ.3) then
  do i=1,nsph
   weights(i)=mass(i)*elecfrac(i)
  enddo
 endif

 if(mapnr.EQ.4.OR.mapnr.EQ.5) then
  do i=1,nsph
   weights(i)=mass(i)*entropy(i)
  enddo
  if(mapnr.EQ.5) then
   do i=1,nfeedback
    weights(i+nsph)=snentropy(bodlist(i))
    mass(i+nsph)=0.
    pos(i+nsph,1)=pos(bodlist(i),1)
    pos(i+nsph,2)=pos(bodlist(i),2)
    pos(i+nsph,3)=pos(bodlist(i),3)
   enddo
   nsph=nsph+nfeedback
  endif
 endif
 
 call makesphtree

 call initextract

 call extractdensity(weights)

 call endextract

 deallocate(weights)

 if(mapnr.EQ.5) nsph=nsph-nfeedback

 print*,'again? (0=yes)'
 read*, mapnr
 if(mapnr.EQ.0) goto 10

end program

