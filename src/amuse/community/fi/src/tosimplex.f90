subroutine basicinput
include 'globals.h'
character*200 filenaam
integer i,mapnr,nfeedback
real,allocatable :: weights(:)
 
 call set_parameters(0)

 print*,'filenaam?'
 read(*,'(a)') filenaam
 

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

end subroutine

program snapreader
 use StarsMod
 use ionMod
 use ElementsMod
 include 'globals.h'
 integer i,p,nvertex,nn,nlum,hy
 real timefac,starage, h,d,dd,meant,meanxe,ttemp,temp,xe,ndens

 real, allocatable :: x(:),y(:),z(:),dens(:),xf(:),vtemp(:),lum(:)
 real, allocatable :: weights1(:),weights2(:)
 
 call initmem(nbodsmax,nsphmax,ncells)
 call basicinput

 print*,densconst,crionrate,mumhkg1

 timefac=timescale/year

! insert code after this
 call makesphtree

 nvertex=nsph
 do p=nbodies-nstar+1,nbodies 
  starage=(tnow-tform(p))*timefac
  if(starage.LT.3.e7) nvertex=nvertex+1
 enddo
 
 allocate( x(nvertex),y(nvertex),z(nvertex), &
      lum(nvertex),dens(nvertex),vtemp(nvertex),xf(nvertex))
 allocate(weights1(nsph),weights2(nsph))

 x(1:nsph)=pos(1:nsph,1)
 y(1:nsph)=pos(1:nsph,2)
 z(1:nsph)=pos(1:nsph,3)
 lum(1:nsph)=0.
 dens(1:nsph)=rho(1:nsph)*densconst*meanmwt*amu
 vtemp(1:nsph)=temperat(1:nsph)
 xf(1:nsph)=elecfrac(1:nsph)
 
 weights1(1:nsph)=temperat(1:nsph)*mass(1:nsph)
 weights2(1:nsph)=elecfrac(1:nsph)*mass(1:nsph)

 nlum=0
 do i=1,nstar-nbh
  p=nbodies-nstar+i 
  starage=(tnow-tform(p))*timefac
  if(starage.LT.3.e7) then
   nlum=nlum+1
   x(nsph+nlum)=pos(p,1)
   y(nsph+nlum)=pos(p,2)
   z(nsph+nlum)=pos(p,3)
   lum(nsph+nlum)=Nlya(starage)*unitm_in_msun*mass(p)
   call hsmdenspos2(pos(p,1:3),h,d,dd,nn)
   bodlist(1:nn)=srlist(1:nn) ! dirty fix
   call gatter3(nn,pos(p,1:3),h,d,dd,meant,weights1)
   call gatter3(nn,pos(p,1:3),h,d,dd,meanxe,weights2)
   meant=meant/d
   meanxe=meanxe/d
   dens(nsph+nlum)=d*densconst*meanmwt*amu
   vtemp(nsph+nlum)=meant
   xf(nsph+nlum)=meanxe
  endif
 enddo

 hy=elementnrQ("H   ")
 
 
 open(unit=1,file='vertices.txt',status='unknown')
 write(1, '(a, g14.6)') "meanmwt: ", meanmwt
 write(1, '(a, g14.6)') "f_hydrogen: ", fractionQ(hy)
 write(1, *) "nvertex: ", nvertex
 write(1,'(7a14)') "x(kpc)","y(kpc)","z(kpc)","Nlya(#/s)","rho(g/cm^3)","T(K)","x_e"
 do i=1,nvertex
   write(1,'(4X,7g14.6)') x(i),y(i),z(i),lum(i),dens(i),vtemp(i),xf(i)
 enddo 
 close(1)
 print*,nlum

 
end program
