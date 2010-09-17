module opacity_co
   ! -----------------------------------------------------------------
   ! This module contains all variables relating to the CO opacity
   ! tables. Adapted from Eldridge implemented by SdM
   ! -----------------------------------------------------------------
   ! CO_MT:          opacity table size in log10 T direction
   ! CO_MR:          "" in log R = rho^3/T6 direction
   ! CO_MH:          "" in Hyrdrogen abundance direction
   ! CO_CO:          "" in Carbon and oxygen direction 
   ! SPLINE_OPAC_CO: Array for spline coefficients
   ! CBASE, OBASE:   C and O abundance because of metallicity only
   ! opT, opR:       log T and log R on which opacity table is defined
   ! opH, COcompos:  Compositions on which opacity table is defined
   ! -----------------------------------------------------------------
   
   use real_kind
   
   implicit none
   integer, parameter :: co_mh = 5  
   integer, parameter :: co_mt = 141
   integer, parameter :: co_mr = 31
   integer, parameter :: co_co = 8
   
   real(double) :: cbase,obase
   
   real(double) :: opt(co_mt), opr(co_mr)
   real(double), parameter :: oph(co_mh) = &
        (/0., 0.03, 0.1, 0.35, 0.7/)
   real(double), parameter :: cocompos(co_co) =  &
        (/0.0d0, 0.01d0, 0.03d0, 0.1d0, 0.2d0, 0.4d0, 0.6d0, 1.0d0/)
   
   real(double), save, allocatable :: spline_opac_co(:,:,:,:,:)
end module opacity_co



! As OPSPLN for the CO-enhanced tables SdM
! Calculate a bicubic spline interpolation fit for the temperature
! and density opacity fit.  This is based on routines from Princeton
! Splines, D. McCune found by L. Dray.
subroutine opspln_co
   use real_kind
   use settings
   use opacity_co
   
   implicit none
   real(double) ::  mat(4,co_mt)
   integer :: i, j, k, ic, ic1, ic2
   
   ! define T and R intervals in opac tables
   ! NB: R = rho/T_6^4 and not rho as in the ld tables
   do i=1,co_mt
      opt(i) = 3.0d0 + 5.0d-2*dble(i-1)
   end do
   do j=1,co_mr
      opr(j) = -8.0d0 + 5.0d-1*dble(j-1)
   end do
   ! Setup CO spline tables
   do k=1,305
      do j=1,co_mr
         do i=1,co_mt
            mat(1,i) = spline_opac_co(1,1,i,j,k)
         end do
         ! Construct splines in the T direction
         call spline(co_mt, opt, mat)
         do i=1,co_mt-1
            do ic1=2,4
               spline_opac_co(ic1,1,i,j,k) = mat(ic1,i)
            end do
         end do
      end do
      
      ! Construct splines in the R (= rho/T_6^3) direction
      do i=1,co_mt-1
         do ic = 1,4
            do j=1,co_mr
               mat(1,j) = spline_opac_co(ic,1,i,j,k)
            end do
            mat(2,1)  = 0.0d0
            mat(3,1)  = 0.0d0
            mat(2,co_mr) = 0.0d0
            mat(3,co_mr) = 0.0d0   
            
            call spline(co_mr, opr, mat)            
            do j=1,co_mr-1
               do ic2 = 2,4                       
                  spline_opac_co(ic,ic2,i,j,k) = mat(ic2,j)
               end do         !IC2                     
            end do            !J           
         end do               !IC 
      end do                  !I
   end do                     !K
end subroutine opspln_co



! Read CO enhanced opacity tables
subroutine load_opacity_co(co_file)
   use real_kind
   use settings
   use constants
   use opacity_co
   
   implicit none
   integer :: i,j,k
   real(double) :: temp
   integer, intent(in) :: co_file
   
   if (.not. allocated(spline_opac_co)) allocate(spline_opac_co(4,4,co_mt,co_mr,305))
   
   read (co_file,'(e10.2)') czs
   
   do k=1,305
      read (co_file,*)   
      do i=1,co_mt
         read (co_file,'(f5.2, 31f7.3)') temp, (spline_opac_co(1,1,i,j,k), j=1,co_mr)
      end do
   end do
   rewind (co_file)
   
   ! Calculate initial H and basic CO abundance as scaled from solar
   ! Don't override CH value from init.dat (CH>0)!
   if (ch<0.0d0) ch = 0.76d0 - 3.0d0*czs
   cbase = czs*cc         
   obase = czs*co 
   
   ! log (Z/Z_sun), used now and then, eg. Vink mass loss rate
   clogz = log10(max(czs/czsn, 1.0d-40))
   
   ! Set up coefficients for cubic spline interpolation in opacity
   call opspln_co
   
end subroutine load_opacity_co



!> -----------------------------------------------------------------------
!! Computes the opacity for a given composition by interpolating in
!! Temperature and density/R recently changed by SdM
!! -----------------------------------------------------------------------
!! INPUT variables - JX: index referring to a certain composition
!!                 - FT: log10 T 
!!                 - FR: log10 R = log (T6/rho^3)
!!
!! OUTPUT Variables
!!   - FKL: log opacity for nearest lower composition value in table
!!   - FKH: log opacity for nearest higher composition value in table
!! The actual opacity is found by doing a linear interpolation between
!! FKL and FKH
!! 
!! Uses:
!!  + MODULE OPACITY_CO
!!    - SPLINE_OPAC_CO: spline coefficients for CO enhanced opacity tables
!!    - opT(1:CO_MT): log10 T at which CO opacities are defined
!!    - opR(1:CO_MR): log10 R at which CO opacities are defined 
!!  + Former common block STAT3
!!
!!  \todo FIXME: CB no longer exists, variable names changed
!!
!!    - F: spline coefficients for Pols opacity tables
!!    - FTM(127): log10 T at which Pols opacities are defined
!!    - FRM(90): log10 Rho at which Pols opacities are defined
!! -----------------------------------------------------------------------
!<
subroutine opacty_co ( jjx, ft, fr, fkl, fkh )
   use real_kind
   use opacity_co
   
   implicit none
   integer, intent(in) :: jjx
   real(double), intent(in) :: ft, fr
   real(double), intent(out) :: fkl, fkh
   
   real(double) :: xft, xfr, dt, dr
   integer :: jx, next_jx, i1, i2
   
   ! Calculate a bicubic spline interpolation fit for the temperature
   ! and density opacity fit.  Do not stop if the input lies outside the
   ! array but rather use the value at the nearest edge-point.
   xft = dmax1(dmin1(ft, opt(co_mt)), opt(1))
   xfr = dmax1(dmin1(fr, opr(co_mr)), opr(1))
   ! Find interval in which target point lies. Index I1 is index of
   ! nearest lower(!) entry in opT/opR
   
   i1 = 1 + int((co_mt - 1)*(xft - opt(1))/(opt(co_mt) - opt(1)))
   i2 = 1 + int((co_mr - 1)*(xfr - opr(1))/(opr(co_mr) - opr(1)))
   
   dt = xft - opt(i1)
   dr = xfr - opr(i2)         
   
   ! Evaluate the splines.      
   next_jx = 61
   jx = jjx
   if (jx+next_jx>305) then
      jx = jx-next_jx
   end if
   
   fkl = spline_opac_co(1, 1, i1, i2, jx) + dr*(spline_opac_co(1, 2, i1, i2, jx)&
        + dr*(spline_opac_co(1, 3, i1, i2, jx) + dr* spline_opac_co(1, 4, i1, i2, jx)))&
        + dt*(spline_opac_co(2, 1, i1, i2, jx) + dr*(spline_opac_co(2, 2, i1, i2, jx)&
        + dr*(spline_opac_co(2, 3, i1, i2, jx) + dr* spline_opac_co(2, 4, i1, i2, jx)))&
        + dt*(spline_opac_co(3, 1, i1, i2, jx) + dr*(spline_opac_co(3, 2, i1, i2, jx)&
        + dr*(spline_opac_co(3, 3, i1, i2, jx) + dr* spline_opac_co(3, 4, i1, i2, jx)))&
        + dt*(spline_opac_co(4, 1, i1, i2, jx) + dr*(spline_opac_co(4, 2, i1, i2, jx)&
        + dr*(spline_opac_co(4, 3, i1, i2, jx) + dr* spline_opac_co(4, 4, i1, i2, jx))))))
   
   fkh = spline_opac_co(1, 1, i1, i2, jx + next_jx) + dr*(spline_opac_co(1, 2, i1, i2, jx + next_jx)&
        + dr*(spline_opac_co(1, 3, i1, i2, jx + next_jx) + dr* spline_opac_co(1, 4, i1, i2, jx + next_jx)))&
        + dt*(spline_opac_co(2, 1, i1, i2, jx + next_jx) + dr*(spline_opac_co(2, 2, i1, i2, jx + next_jx)&
        + dr*(spline_opac_co(2, 3, i1, i2, jx + next_jx) + dr* spline_opac_co(2, 4, i1, i2, jx + next_jx)))&
        + dt*(spline_opac_co(3, 1, i1, i2, jx + next_jx) + dr*(spline_opac_co(3, 2, i1, i2, jx + next_jx)&
        + dr*(spline_opac_co(3, 3, i1, i2, jx + next_jx) + dr* spline_opac_co(3, 4, i1, i2, jx + next_jx)))&
        + dt*(spline_opac_co(4, 1, i1, i2, jx + next_jx) + dr*(spline_opac_co(4, 2, i1, i2, jx + next_jx)&
        + dr*(spline_opac_co(4, 3, i1, i2, jx + next_jx) + dr* spline_opac_co(4, 4, i1, i2, jx + next_jx))))))
   
   return      ! We normally don't want any output here...
   if ( (xft /= ft) .or. (xfr /= fr) ) write (10,100) ft, fr
   
100 format ('opacity out of range',' ft fr ', 2f9.4)
end subroutine opacty_co


