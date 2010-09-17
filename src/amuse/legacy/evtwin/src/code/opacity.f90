module opacity_data
   use real_kind
   
   implicit none
   integer, parameter :: kt=127  ! Number of temperature poinst
   integer, parameter :: kr=90   ! Number of density points

   real(double), save :: fspl(4,4,kt,kr,10)  ! spline interpolation tables
   real(double), save :: tfm(kt)          ! temperatures
   real(double), save :: frm(kr)          ! densities
   real(double), save :: csx(10)          ! Composition values for opacity tables
   real(double), save :: cs(kr,kt,10)     ! Opacity tables
   integer, save :: kcsx            ! Number of compositions
end module opacity_data

! Calculate a bicubic spline interpolation fit for the temperature
! and density opacity fit.  This is based on routines from Princeton
! Splines, D. McCune found by L. Dray.
subroutine opspln
   use real_kind
   use opacity_data
   
   implicit none
   integer :: jt,jr,jx,jq,ir,it,ic,iq
   real(double) :: mat(4,kt)

   do jt = 1, kt
      tfm(jt) = 2.95d0 + 0.05d0*dble(jt)
   end do
   do jr = 1, kr
      frm(jr) = 0.25d0*dble(jr) - 12.25d0
   end do
   do jx = 1, 10
      do jq = 1, kt
         do iq = 1, kr
            fspl(1, 1, jq, iq, jx) = cs(iq, jq, jx)
         end do
      end do
      ! Construct splines in the T direction.
      do ir = 1, kr
         do it = 1, kt
            mat(1, it) = fspl(1, 1, it, ir, jx)
         end do
         call spline ( kt, tfm, mat )
         do it = 1, kt - 1
            fspl(2, 1, it, ir, jx) = mat(2, it)
            fspl(3, 1, it, ir, jx) = mat(3, it)
            fspl(4, 1, it, ir, jx) = mat(4, it)
         end do
      end do
      ! Construct splines in the rho direction.
      do it = 1, kt - 1
         ! Construct splines for each T coeff
         do ic = 1, 4
            do ir = 1, kr
               mat(1, ir) = fspl(ic, 1, it, ir, jx)
            end do
            mat(2, 1) = 0.0d0
            mat(3, 1) = 0.0d0
            mat(2, kr) = 0.0d0
            mat(3, kr) = 0.0d0
            call spline ( kr, frm, mat )
            do ir = 1, kr - 1
               fspl(ic, 2, it, ir, jx) = mat(2, ir)
               fspl(ic, 3, it, ir, jx) = mat(3, ir)
               fspl(ic, 4, it, ir, jx) = mat(4, ir)
            end do
         end do
      end do
   end do

end subroutine opspln



! Calculate the coefficients of a 1-D cubic spline:
! Forsythe, Malcolm, Moler, Computer Methods for Mathematical
! Computations, Prentice-Hall, 1977, p.76
subroutine spline ( k, x, f )
   use real_kind

   implicit none
   integer :: k
   real(double) :: x(*),f(4,*)

   integer :: i,ib
   real(double) :: t

   f(2:4,k) = 0.0d0
   ! Set up a tridiagonal system for A*y=B where y(i) are the second
   ! derivatives at the knots.
   ! f(2,i) are the diagonal elements of A
   ! f(4,i) are the off-diagonal elements of A
   ! f(3,i) are the B elements/3, and will become c/3 upon solution
   f(4,1) = x(2)-x(1)
   f(3,2) = (f(1,2) - f(1,1))/f(4,1)
   do i = 2, k - 1
      f(4,i) = x(i+1) - x(i)
      f(2,i) = 2.0d0*(f(4,i-1) + f(4,i))
      f(3,i+1) = (f(1,i+1) - f(1,i))/f(4,i)
      f(3,i) = f(3,i+1) - f(3,i)
   end do
   ! Boundaries.
   f(2,2) = f(4,1) + 2.0d0*f(4,2)
   f(3,2) = f(3,2)*f(4,2)/(f(4,1) + f(4,2))
   f(2,k-1) = 2.0d0*f(4,k-2) + f(4,k-1)
   f(3,k-1) = f(3,k-1)*f(4,k-2)/(f(4,k-1) + f(4,k-2))
   ! Forward elimination.
   t = f(4,2)/f(2,2)
   f(2,3) = f(2,3) - t*(f(4,2) - f(4,1))
   f(3,3) = f(3,3) - t*f(3,2)
   do i = 4, k - 2
      t = f(4,i-1)/f(2,i-1)
      f(2,i) = f(2,i)-t*f(4,i-1)
      f(3,i) = f(3,i)-t*f(3,i-1)
   end do
   t = (f(4,k-2) - f(4,k-1))/f(2,k-2)
   f(2,k-1) = f(2,k-1) - t*f(4,k-2)
   f(3,k-1) = f(3,k-1) - t*f(3,k-2)
   ! Back substitution.
   f(3,k-1) = f(3,k-1)/f(2,k-1)
   do ib = 1, k - 4
      i = k - 1 - ib
      f(3,i) = (f(3,i) - f(4,i)*f(3,i+1))/f(2,i)
   end do
   f(3,2) = (f(3,2) - (f(4,2) - f(4,1))*f(3,3))/f(2,2)
   ! Reset d array to step size.
   f(4,1) = x(2) - x(1)
   f(4,k-1) = x(k) - x(k-1)
   ! Set f(3,1) for not-a-knot.
   f(3,1) = (f(3,2)*(f(4,1) + f(4,2)) - f(3,3)*f(4,1))/f(4,2)
   f(3,k) = f(3,k-1) + (f(3,k-1) - f(3,k-2))*f(4,k-1)/f(4,k-2)
   ! Compute the polynomial coefficients.
   do i = 1, k - 1
      f(2,i) = (f(1,i+1) - f(1,i))/f(4,i) - f(4,i)*(f(3,i+1) &
           + 2.0d0*f(3,i))
      f(4,i) = (f(3,i+1) - f(3,i))/f(4,i)
      f(3,i) = 3.0d0*f(3,i)
      f(4,i) = f(4,i)
   end do

end subroutine spline



! Read opacity tables
subroutine load_opacity(file)
   use real_kind
   use settings
   use constants
   use opacity_data
   
   implicit none
   integer, intent(in) :: file

   integer :: i,jr,jt
   real(double) :: ch_opac

990 format (1x, 10f7.3)
   ! Read opacity, nuclear reaction and neutrino loss rate data
   read (file,*) kcsx, czs, ch_opac
   ! Don't override CH value from init.dat (CH>0)!
   if (ch<0.0d0) ch = ch_opac
   read (file,990) csx
   do i = 1, kcsx
      read (file,990) ((cs(jr, jt, i), jr = 1, 90), jt = 1, 127)
   end do
   ! log (Z/Z_sun), used now and then, eg. Vink mass loss rate
   clogz = log10(max(czs/czsn, 1.0d-40))
   rewind (file)
   ! Set up coefficients for cubic spline interpolation in opacity
   call opspln
end subroutine load_opacity



! Calculate a bicubic spline interpolation fit for the temperature
! and density opacity fit.  Do not stop if the input lies outside the
! array but rather use the value at the nearest edge-point.
! Uses the older opacity tables
subroutine opacty ( jx, tf, frho, fkl, fkh )
   use real_kind
   use opacity_data
   
   implicit none
   integer :: jx
   real(double) :: tf, frho, fkl, fkh

   integer :: i1,i2
   real(double) :: xtf,xfr,dt,dr

   xtf = dmax1(dmin1(tf, tfm(kt)), tfm(1))
   xfr = dmax1(dmin1(frho, frm(kr)), frm(1))
   ! Find interval in which target point lies.
   i1 = 1 + (kt - 1)*(xtf - tfm(1))/(tfm(kt) - tfm(1))
   i2 = 1 + (kr - 1)*(xfr - frm(1))/(frm(kr) - frm(1))
   dt = tf - tfm(i1)
   dr = frho - frm(i2)
   ! Evaluate the splines.
   fkl = fspl(1, 1, i1, i2, jx) + dr*(fspl(1, 2, i1, i2, jx)               &
        + dr*(fspl(1, 3, i1, i2, jx) + dr*fspl(1, 4, i1, i2, jx)))         &
        + dt*(fspl(2, 1, i1, i2, jx) + dr*(fspl(2, 2, i1, i2, jx)          &
        + dr*(fspl(2, 3, i1, i2, jx) + dr*fspl(2, 4, i1, i2, jx)))         &
        + dt*(fspl(3, 1, i1, i2, jx) + dr*(fspl(3, 2, i1, i2, jx)          &
        + dr*(fspl(3, 3, i1, i2, jx) + dr*fspl(3, 4, i1, i2, jx)))         &
        + dt*(fspl(4, 1, i1, i2, jx) + dr*(fspl(4, 2, i1, i2, jx)          &
        + dr*(fspl(4, 3, i1, i2, jx) + dr*fspl(4, 4, i1, i2, jx))))))
   fkh = fspl(1, 1, i1, i2, jx + 1) + dr*(fspl(1, 2, i1, i2, jx + 1)       &
        + dr*(fspl(1, 3, i1, i2, jx + 1) + dr*fspl(1, 4, i1, i2, jx + 1))) &
        + dt*(fspl(2, 1, i1, i2, jx + 1) + dr*(fspl(2, 2, i1, i2, jx + 1)  &
        + dr*(fspl(2, 3, i1, i2, jx + 1) + dr*fspl(2, 4, i1, i2, jx + 1))) &
        + dt*(fspl(3, 1, i1, i2, jx + 1) + dr*(fspl(3, 2, i1, i2, jx + 1)  &
        + dr*(fspl(3, 3, i1, i2, jx + 1) + dr*fspl(3, 4, i1, i2, jx + 1))) &
        + dt*(fspl(4, 1, i1, i2, jx + 1) + dr*(fspl(4, 2, i1, i2, jx + 1)  &
        + dr*(fspl(4, 3, i1, i2, jx + 1) + dr*fspl(4, 4, i1, i2, jx + 1))))))
   return      ! We normally don't want any output here...
   if ( (xtf /= tf) .or. (xfr /= frho) ) write (10,100) tf, frho
   if(fkl >  10. .or. fkh >  10.) write(10,*) "opacity too high", fkl, fkh
   if(fkl < -15. .or. fkh < -15.) write(10,*) "opacity too low", fkl, fkh
100 format ('opacity out of range',' tf frho ', 2f9.4)

end subroutine opacty



! ------------------------------------------------------------------
! GET_OPACITY
! Calculates interpolated opacities as function of density and
! temperature for the given composition.
! ------------------------------------------------------------------
! INPUT:
!    FRHO  - log10 rho, density
!    FT    - log10 T, temperature
!    ABUND - abundance struct, for the requested abundance
! RETURNS:
!    The total opacity, in cm^2/g
!
! Control variable:
! KOP = 1: for opacities as implemented by Pols & al.
!       2: wrong! opacities assuming a never changing composition
!       3: CO enhanced opacities, after Eldridge & Tout
!       4: as 3, while every element heavier than oxygen contributes
!          to the opacity as oxygen does
! ------------------------------------------------------------------
function get_opacity(frho, ft, abund)
   use real_kind
   use constants
   use settings
   use opacity_data
   use opacity_co
   use atomic_data
   use eostate_types
   
   implicit none
   real(double), intent(in)  :: frho, ft ! density, temperature
   type(abundance), intent(in) :: abund

   integer :: nc, no, iic, iio, jx=1
   integer :: k2,jjx,jk
   real(double) :: tkappa(2,2)
   real(double) :: xh,xhe,xf,xt,xu,fkl,fkh,fk,fr,xxc,xxo,fc,fo
   real(double) :: fkl_co,fkh_co,fk_co,ec,eo

   real(double) :: get_opacity

   !Common blocks:
   real(double) :: na(9), avm

   na(:) = abund%na(:)
   avm = abund%avm

   xh = na(1)*can(1)/avm      ! mass fraction of Hydrogen
   xhe = na(2)*can(2)/avm     ! mass fraction of Helium
   fk = 0.0d0

   if (kop == 1) then
      ! Opacity tables from Alexander & Ferguson (1994; molecular), Itoh (1983;
      ! electron conduction) and Iglesias & Rogers (1992; the rest)
      ! Find XF interval in opacity table.
      xf  = xh + xh + xhe + 1.0d0
      !> \todo FIXME: should use a more efficient binary search here...
      do jx=1,kcsx-1
         if( csx(jx+1)<xf ) exit
      end do
      xt =  (xf-csx(jx))/(csx(jx+1)-csx(jx))
      xu = 1.0d0 - xt
      ! XT = MIN(1.0D0,(MAX(0.0D0,XF))) ! make sure XT between 0 and 1.
      ! XU = MIN(1.0D0,(MAX(0.0D0,XU))) ! make sure XU between 0 and 1.

      ! bicubic spline interpolation
      call opacty ( jx, ft, frho, fkl, fkh)
      fk = xt*10.0d0**min(fkh, 3.0d2) + xu*10.0d0**min(fkl, 3.0d2)
   else if (kop > 1 .or. kop < 0) then     ! KOP<0 for debug purposes
      ! ------------------------------------------------------------------------
      ! Use new opacity tables as implemented after ELdridge and Tout MNRAS 2003
      ! ------------------------------------------------------------------------
      ! OPAL tables (Iglesias & Rogers 1996) for different H, C, and O
      ! fractions extended to lower temperatures with Alexander & Ferguson
      ! 1994 (! or more up to date?) (only for different H fractions) and
      ! extended to higher temperatures with Buchler & Yueh 1976. Finally the
      ! opacity owing to electron conduction is added according to Hubbard &
      ! Lampe (1969) and Canuto 1970

      ! The new opacity tables use the variable FR rather than the density
      fr = frho - 3.0d0*(ft -6.0d0) ! FR = log10 (Density/T_6^3)
      ! Determine in which interval of XH we are and prepare for linear interpolation
      !> \todo FIXME: no good reason not to use a binary search here...
      do jx=1,co_mh-1
         if ( oph(jx+1)>xh ) exit
      end do
      ! JX will be the index of the element in opH with nearest lower H
      ! abundance on which the opacity table is defined. Note that opH has a
      ! descending H abundance
      jx = min(jx, co_mh-1)
      xt =  (xh-oph(jx))/(oph(jx+1)-oph(jx))
      xt = min(1.0d0,(max(0.0d0,xt))) ! make sure XT between 0 and 1.
      ! ?? Could we extrapolate for X > 0.7, at low Z
      xu = 1.0d0 - xt
      xu = min(1.0d0,(max(0.0d0,xu))) ! make sure XU between 0 and 1.

      ! Composition parameters used for interpolation: (FC=0: Carbon
      ! abundance as solar but scaled with metallicity, FC=1 max enhancement)
      xxc = na(3)*can(3)/avm  ! mass fraction of carbon
      xxo = na(5)*can(5)/avm  ! mass fraction of oxygen
      fc = (xxc - cbase)/(1.0d0-czs-xh)
      fo = (xxo - obase)/(1.0d0-czs-xh)
      fc = min(1.0d0, max(0.0d0, fc))
      fo = min(1.0d0, max(0.0d0, fo))

      ! Include every element heavier than oxygen in that abundance
      ! slightly different from JJE's code, but it seems right to me=SdM

      if (kop >= 4.and.xh<1.0d-10) then
         xxo = 1.0d0 - xh- xhe -xxc
         fo = (xxo - obase)/(1.0d0-czs-xh)
      end if

      ! For KOP=2, or if the enhancement of C and O is very small the opacity
      ! is taken for a mixture where C and O are never enhanced and only CBASE
      ! and OBASE scale with Z
      if(kop == 2 .or. (fc<1.0d-6 .and. fo<1.0d-6)) then
         call opacty_co(1+(jx-1)*61, ft, fr, fkl_co, fkh_co)
         fk_co = xt*10.0d0**min(fkh_co,3.0d2) + &
              xu*10.0d0**min(fkl_co,3.0d2)
         fk = fk_co
      else                   !Do full 5d interpolation for KOP>=3
         ! for very small  C and O, neglect them
         if (fc < 1.0d-20) fc = 0.0d0
         if (fo < 1.0d-20) fo = 0.0d0
         ! Find relevant CO fraction interval in opacity table
         iic = 1
         iio = 1
         do k2=1,7
            if(cocompos(k2) < fc)  iic = k2
            if(cocompos(k2) < fo)  iio = k2
         end do

         iic = max(min(7,iic), 1)
         iio = max(min(7,iio),1)
         ! For each FC there are tables with 8 different FO's except for
         ! FC = 0.6 (7 tables) and FC = 1.0 (6 tables)
         ! For construction of the table see explanation JJE, 10 October 2003
         ! "implementation of Enhanced CO Opacity Table", website
         ! IIO and IIC are swapped compared to JJE's implementation
         if(iic<=7) then
            jjx = (jx-1)*61 + (iic-1)*8 + iio
         else
            jjx = (jx-1)*61 +       6*8 + 7 + iio
         end if

         ! Calculate opacity at all four lower and upper levels for C and O abundance
         ! Varying O direction +/- 1, Varying in C direction +/- 8 (or 7)
         ! IIO and IIC are swapped compared to JJE's implementation
         do nc=0,1
            do no = 0,1
               if(iic <= 6) then
                  jk = jjx + 8*nc + no
               else
                  jk = jjx + 7*nc + no
               end if

               call opacty_co(jk, ft,fr,fkl_co,fkh_co)
               fk_co = xt*10.0d0**min(fkh_co,3.0d2) + &
                    xu*10.0d0**min(fkl_co,3.0d2)
               ! Store opacities of nearest point in the table for later interpolation
               tkappa(nc+1,no+1) = fk_co
            end do
         end do
         ! Final 2D lineair interpolation for C and O
         ec = ( fc-cocompos(iic) )/( cocompos(iic+1)-cocompos(iic) )
         eo = ( fo-cocompos(iio) )/( cocompos(iio+1)-cocompos(iio) )

         ec = min(max(ec, 1.0d-50), 1.0d0)
         eo = min(max(eo, 1.0d-50), 1.0d0)

         fk_co = (1.0d0-ec) * (1.0d0-eo) * tkappa(1,1) +&
              ec  * (1.0d0-eo) * tkappa(2,1) +&
              (1.0d0-ec) *        eo  * tkappa(1,2) +&
              ec  *        eo  * tkappa(2,2)

         fk = fk_co
      end if
   end if
   get_opacity = fk

end function get_opacity
