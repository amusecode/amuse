! declaration of large arrays
module memMod

public
private :: ndim,nsubcell

  integer, parameter :: ndim=3, nsubcell=2**ndim

! srlist: threadprivate array for search results
! bodlist: threadprivate integer array of treewalk results
! tempvect: threadprivate real array, scratch space
! templist: temporary integer array
!
! rest is cell or body data

  integer, allocatable :: bodlist(:),cellstart(:),  &
    itimestp(:),nbexist(:),npercell(:),order_bodlist(:),  &
    otimestp(:),pactive(:),srlist(:),pshape(:),subp(:,:),templist(:)

  real, allocatable :: acc(:,:),bottom(:,:),cellsize(:),csound(:),  &
    derad(:),drhodh(:),elecfrac(:),epsgrav(:),  &
    esnthdt(:),fuvheat(:),h2frac(:),hsmcurlv(:),  &
    hsmdivv(:),hsmooth(:),mass(:),mumaxdvh(:),oldderad(:),phi(:),  &
    phiext(:),pos(:,:),quad(:,:),rho(:),snentropy(:),starfuv(:),  &
    tcollaps(:),temperat(:),tempvect(:),tfeedb(:),tform(:),tvel(:),  &
    vdisp(:),vel(:,:),veltpos(:,:)
    
  real, allocatable, target :: ethermal(:),ethold(:),dethdt(:),dethold(:)
  real, pointer :: entropy(:),entold(:),dentdt(:),dentold(:)

!$omp threadprivate(bodlist,srlist,tempvect)

  integer :: ppropcount, pordercount, treestatecount, sphtreecount
      ! ppropcount -> +1 if particle prop changed
      ! (mass,pos,epsgrav,starfuv,hsmooth)
      ! pordercount -> +1 particle order changed


 contains

 subroutine initmem(nbodsmax,nsphmax,ncells)
  integer nbodsmax,nsphmax,ncells
  integer nbods1,nbodcell

  nbods1=nbodsmax+1
  nbodcell=nbodsmax+ncells


 allocate( &
   cellstart(nbods1:nbodcell),  &
   itimestp(1:nbodsmax),  &
   nbexist(1:nbodsmax),  &
   npercell(nbods1:nbodcell),  &
   order_bodlist(1:nbodsmax),  &
   otimestp(1:nbodsmax),  &
   pactive(1:nbodsmax),  &
   pshape(nbods1:nbodcell),  &
   subp(nbods1:nbodcell,1:nsubcell),  &
   templist(1:nbodsmax) &
 )
 allocate( &
   acc(1:nbodsmax,1:ndim+1),  &
   bottom(nbods1:nbodcell,ndim),  &
   cellsize(nbods1:nbodcell),  &
   csound(1:nsphmax),  &
   derad(1:nsphmax),  &
   dethdt(1:nsphmax),  &
   dethold(1:nsphmax), &
   drhodh(1:nsphmax),  &
   elecfrac(1:nsphmax),  &
   epsgrav(1:nbodcell),  &
   esnthdt(1:nsphmax),  &
   ethermal(1:nsphmax),  &
   ethold(1:nsphmax),  &
   fuvheat(1:nbodsmax),  &
   h2frac(1:nsphmax),  &
   hsmcurlv(1:nsphmax),  &
   hsmdivv(1:nsphmax),  &
   hsmooth(1:nbodcell),  &
   mass(1:nbodcell),  &
   mumaxdvh(1:nsphmax),  &
   oldderad(1:nsphmax),  &
   phi(1:nbodsmax),  &
   phiext(1:nbodsmax),  &
   pos(1:nbodcell,1:ndim),  &
   quad(nbods1:nbodcell,1:2*ndim-1),  &
   rho(1:nsphmax),  &
   snentropy(1:nbodsmax),  &
   starfuv(1:nbodcell),  &
   tcollaps(1:nsphmax),  &
   temperat(1:nsphmax),  &
   tfeedb(1:nbodsmax),  &
   tform(1:nbodsmax),  &
   tvel(1:nbodsmax),  &
   vdisp(1:nsphmax),  &
   vel(1:nbodsmax,1:ndim),  &
   veltpos(1:nsphmax,1:ndim) &
 )

!$ call omp_set_dynamic(.FALSE.)

!$omp parallel
 allocate( &
   tempvect(1:nbodsmax+nsphmax),  &
   bodlist(1:nbodsmax),  &
   srlist(1:nbodcell) &
 )
!$omp end parallel

 entropy=>ethermal
 entold=>ethold
 dentdt=>dethdt
 dentold=>dethold

 ppropcount=0
 pordercount=0
 treestatecount=-1
 sphtreecount=-1
 
 end subroutine
 
end module 
