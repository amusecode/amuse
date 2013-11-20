!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters,only:boxlen
  use poisson_parameters  
  use pm_commons,only:msink, xsink

  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp):: rr,rx,ry,rz

  ! Constant vector
  if(gravity_type==1)then 
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then 
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xsink(1,1)
        ry=x(i,2)-xsink(1,2)
        rz=x(i,3)-xsink(1,3)
        rr=sqrt(rx**2+ry**2+rz**2+gravity_params(1)**2) ! Softened gravity
        f(i,1)=-msink(1)*rx/rr**3
        f(i,2)=-msink(1)*ry/rr**3
        f(i,3)=-msink(1)*rz/rr**3
     end do
  end if

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana
!#########################################################
!#########################################################
!#########################################################
