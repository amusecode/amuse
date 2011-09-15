!> \file update_particles.F90

!> \brief matches particles from one snapshot to the next
!<

module update_particles_mod
use myf03_mod
use gadget_public_input_mod
use gadget_cosmoBH_input_mod
use gadget_owls_input_mod
use gadget_vbromm_input_mod
use gadget_public_input_hdf5_mod
use global_mod, only: psys, GV
implicit none
private

public :: update_particles

contains

!> Reads in an update snapshot.  The idea is to keep the ionization fractions
!! and temperature from the already loaded snapshot while updating the 
!! positions, velocities, smoothing lengths ... from the snapshot on disk.
!! The particles are reordered during the raytracing so when a new snapshot
!! is readin they must be matched
!===========================================================================
subroutine update_particles()

  character(clen), parameter :: myname="update_particles" 
  logical, parameter :: crash=.true.

  integer(i4b), allocatable :: idold(:)
  real(r4b), allocatable :: Told(:)
  real(r4b), allocatable :: yeold(:)
  real(r4b), allocatable :: xHIold(:)
  real(r4b), allocatable :: xHIIold(:)

  real(r4b), allocatable :: xHeIold(:)
  real(r4b), allocatable :: xHeIIold(:)
  real(r4b), allocatable :: xHeIIIold(:)

  integer(i8b), allocatable :: lasthitold(:)
  integer(i8b), allocatable :: orderold(:)  

  integer(i8b) :: minIDold, minIDnew
  integer(i8b) :: maxIDold, maxIDnew

  integer(i8b) :: i, idnew, indx, err




  ! store carry over variables
  !==============================

  allocate( idold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate IDold",myname,crash)     
  idold = psys%par%id
  minIDold = minval(idold)
  maxIDold = maxval(idold)

  allocate( Told(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate Told",myname,crash)     
  Told = psys%par%T

  allocate( yeold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate yeold",myname,crash)     
  yeold = psys%par%ye

  allocate( xHIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHIold",myname,crash)     
  xHIold = psys%par%xHI

  allocate( xHIIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHIIold",myname,crash)     
  xHIIold = psys%par%xHII

#ifdef incHe
  allocate( xHeIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHeIold",myname,crash)     
  xHeIold = psys%par%xHeI

  allocate( xHeIIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHeIIold",myname,crash)     
  xHeIIold = psys%par%xHeII

  allocate( xHeIIIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHeIIIold",myname,crash)     
  xHeIIIold = psys%par%xHeIII
#else
  allocate( xHeIold(1), xHeIIold(1), xHeIIIold(1) )
#endif



  allocate( lasthitold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate lasthitold",myname,crash)     
  lasthitold = psys%par%lasthit


  ! deallocate the current particle array and read new particle data
  !==================================================================
  deallocate( psys%par )
  if (GV%InputType == 1) then
     call read_Gpublic_particles()
  else if (GV%InputType == 2) then
     call read_GcosmoBH_particles()
  else if (GV%InputType == 3) then
     call read_Gowls_particles()
  else if (GV%InputType == 4) then
     call read_Gvbromm_particles()
  else if (GV%InputType == 5) then
     call read_Gpubhdf5_particles()
  end if
  minIDnew = minval(psys%par%id)
  maxIDnew = maxval(psys%par%id)


  ! match id's from new snapshot data to the old IDs
  !===================================================


  ! create an indexing array for the old particle array ID's
  ! ID's that are no longer present will be indexed with minIDold - 1
  !--------------------------------------------------------------------
  write(*,*) "creating indexing array for current particle system..."
  write(*,*) "min/max ID (old snap) = ", minIDold, maxIDold

  allocate( orderold(minIDold:maxIDold), stat=err)
  if (err/=0) call myerr("allocating orderold array",myname,crash)

  orderold = minIDold - 1
  do i = 1,size(idold)
     orderold(idold(i)) = i
  end do

  ! transfer all carryover variables
  !-----------------------------------
  do i = 1,size(psys%par)

     idnew = psys%par(i)%id   
     indx  = orderold(idnew)

     if ( idnew /= idold(indx) ) then
        write(*,*) "i,idnew,idold", i, idnew, idold(indx)
        if (err/=0) call myerr("reordering id error",myname,crash)
     end if

     if (.not. GV%FixSnapTemp) then
        psys%par(i)%T = Told(indx)
     end if

     psys%par(i)%ye   = yeold(indx) 
     psys%par(i)%xHI  = xHIold(indx)
     psys%par(i)%xHII = xHIIold(indx)

#ifdef incHe
     psys%par(i)%xHeI   = xHeIold(indx)
     psys%par(i)%xHeII  = xHeIIold(indx)
     psys%par(i)%xHeIII = xHeIIIold(indx)
#endif

     psys%par(i)%lasthit = lasthitold(indx)

  end do

  deallocate (idold, Told, yeold, xHIold, xHIIold, lasthitold, orderold)

  deallocate (xHeIold, xHeIIold, xHeIIIold)



end subroutine update_particles


end module update_particles_mod
