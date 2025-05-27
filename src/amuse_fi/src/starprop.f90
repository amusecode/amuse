module StarsMod
 implicit none
 private
 public:: InitStars,FUV,dFUV,nbandsQ,mbands,Ha,axavQ,Nlya

! call subroutine InitStars before use
! FUV(t): far UV in erg/sec/Msolar of initial mass for t in years
! dFUV(t) far UV corrected for absorption by parental cloud
! nbandsQ(): Queries number of colorbands read in
! mbands(): gives luminosities in UBV.. in 10^magnitude/Msolar ? 
! Ha(): halpha flux in erg/sec/Msolar/ Lya->Ha factor
! Nlya(): Ly a continuum in # photons/sec (yr?)/Msol
 character(len=30) :: Stellarpropertiesfile
 character(len=30) :: modelprefix='ssp_salp_' ! single stellar burst
                                              ! salpeter IMF

 integer, save :: nbands,ntimetable
 real, allocatable,save :: axav(:)
 real, allocatable,save :: timetable(:),bands(:,:)
 real, allocatable,save :: bolometric(:),farUV(:), &
                           mass(:),snr(:),halpha(:),lya(:) 
 character(len=200),save :: datadir

 real, parameter :: tdark=4.e6
 real, parameter :: obscuration=.75

 logical, parameter :: debug=.FALSE.

 contains

subroutine InitStars(indatadir,Z)
 character :: dummy
 character(len=200) :: indatadir
 real, optional :: Z ! Z=Zsolar=0.02
 integer :: i,j

 datadir=indatadir
 if(allocated(axav)) then
   print*, '** warning Starsmod already initialized **'
   return
 endif
 Stellarpropertiesfile=trim(modelprefix)//'z02'//'.starmodel'
 if(present(Z)) then 
  if(Z.LE.0.0002)  &
   Stellarpropertiesfile=trim(modelprefix)//'z0001'//'.starmodel'
  if(Z.GT.0.0002.AND.Z.LE.0.001) &
   Stellarpropertiesfile=trim(modelprefix)//'z0004'//'.starmodel'
  if(Z.GT.0.001.AND.Z.LE.0.006) &
   Stellarpropertiesfile=trim(modelprefix)//'z004'//'.starmodel'
  if(Z.GT.0.006.AND.Z.LE.0.015) &
   Stellarpropertiesfile=trim(modelprefix)//'z008'//'.starmodel'
  if(Z.GT.0.015.AND.Z.LE.0.03) &
   Stellarpropertiesfile=trim(modelprefix)//'z02'//'.starmodel'
  if(Z.GT.0.03.AND.Z.LE.0.075) &
   Stellarpropertiesfile=trim(modelprefix)//'z05'//'.starmodel'
  if(Z.GT.0.075) &
   Stellarpropertiesfile=trim(modelprefix)//'z10'//'.starmodel'
 endif ! these are somewhat subjective limits!

 write(*,*) ' > starcluster model: ',Stellarpropertiesfile
 open(unit=19,file=trim(datadir)//Stellarpropertiesfile,status='old')
 
 read(19,*) nbands
 allocate(axav(nbands))
 read(19,*) (axav(i),i=1,nbands)
 read(19,*) ntimetable
 allocate(bands(ntimetable,nbands),bolometric(ntimetable),     &
          farUV(ntimetable),mass(ntimetable),snr(ntimetable),  &
          timetable(ntimetable),halpha(ntimetable),lya(ntimetable))
 read(19,*) dummy
 do i=1,ntimetable
 read(19,*) timetable(i), (bands(i,j),j=1,nbands),bolometric(i), &
            farUV(i),mass(i),snr(i),lya(i)
 enddo
 close(19)
 farUV=10**farUV
 bands=10**(-.4*bands)
 halpha=10**(lya-11.52)
 lya=10**lya
 
end subroutine

function mbands(t) result(x)
 real :: t,x(nbands),dt
 integer :: ntime
 call findtbin(t,dt,ntime)
  x=bands(ntime,1:nbands)+dt*(bands(ntime+1,1:nbands)-bands(ntime,1:nbands))
end function 

function dFUV(t) result(x)
 real :: t,x,obsc
 obsc=0.
 if(t.lt.tdark) obsc=obscuration*(1-t/tdark)
 if(t.lt.0) obsc=obscuration
 x=FUV(t)*(1.-obsc)
end function 

function FUV(t) result(x)
 real :: t,x,dt
 integer :: ntime
 call findtbin(t,dt,ntime)
 x=farUV(ntime)+dt*(farUV(ntime+1)-farUV(ntime))
end function

function Ha(t) result(x)
 real :: t,x,dt
 integer :: ntime
 call findtbin(t,dt,ntime)
 x=halpha(ntime)+dt*(halpha(ntime+1)-halpha(ntime))
end function

function Nlya(t) result(x)
 real :: t,x,dt
 integer :: ntime
 call findtbin(t,dt,ntime)
 x=lya(ntime)+dt*(lya(ntime+1)-lya(ntime))
end function

function nbandsQ() result(x)
 integer :: x
 x=nbands
end function

function axavQ(i) result(x)
 integer :: i
 real :: x
 x=axav(i)
end function


subroutine findtbin(t,dt,ntime)
 real,intent(inout) :: t
 real,intent(out) :: dt
 integer, intent(out) :: ntime
 integer :: ntime1,itime
 
 if(debug) then
  if(t.lt.timetable(1).OR.t.gt.timetable(ntimetable)) print*,' starprop:',t
 endif 
 
 t=max(timetable(1),t)
 t=min(timetable(ntimetable),t)
 
 ntime=1
 ntime1=ntimetable 
 do while((ntime1-ntime).GT.1)
  itime=(ntime+ntime1)/2
  if(t.GT.timetable(itime)) then
   ntime=itime  
  else
   ntime1=itime
  endif
 enddo
 dt=(t-timetable(ntime))/(timetable(ntime1)-timetable(ntime))
end subroutine


end module

!program test
! use StarsMod
! real t,ff
! character(len=200) :: indatadir

! indatadir="/home/pelupes/fish/stars/data/"
! call InitStars(indatadir)
 
!10 read*,t
! if(t.eq.0) stop
! ff=FUV(t)
! print*,t,ff
! goto 10 
!end program
