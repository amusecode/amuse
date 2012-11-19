! ------------------------------------------------------------------------------
! UPDATE
!  Update stellar properties after a succesul iteration of the solver.
!  Store current model in case we need to step back and reduce the timestep.
! ------------------------------------------------------------------------------
! Input:
!  DTY - The timestep, in years
! Output:
!  COMMON H(:,:) - the stellar model after the last timestep
! ------------------------------------------------------------------------------
subroutine update ( dty )
   use real_kind
   use mesh
   use mesh_enc
   use test_variables
   use current_model_properties
   use binary_history, only: hpr
   
   implicit none
   real(double) :: dty

   
   ! Store certain current and previous values, for possible emergency backup:
   age = age + dty
   jmod = jmod + 1
   pprev(1:81) = prev(1:81)
   prev(2:81) = (/hspn, rlf, zet, xit, age, bm, mc, om, bper, sm, enc,      &       ! zq(1:17)
        tc, tfr, t0, m0, mta, om0, omta, a0, ata, e0, eta, cdd,             &       ! zq(18:30)
        bp, horb, ro, ra2, rs, secc, tn, wmh, wmhe, mh, mhe, mco,           &       ! zq(31:43)
        vmg, be, lh, lhe, lc, lnu, lth, mcb, msb, rcb, tct/)                        ! zq(44:81)

   jm2 = jm1
   jm1 = jmod

   ! Update the stored models: previous and anteprevious:
   hpr(:, 1:kh) = h(:, 1:kh)
   h(:, 1:kh) = h(:, 1:kh) + dh(:, 1:kh)

   ! Update nucleosynthesis:
   jhold = jhold + 1
end subroutine update



! ------------------------------------------------------------------------------
! BACKUP
!  When convergence fails, back track to the previous model and try again
!  with a reduced timestep.
! ------------------------------------------------------------------------------
! Input:
!  DTY - The timestep, in years
! Output:
!  DTY - The new timestep, in years (normally CT3 times previous timestep)
!  JO  - Return code, 2 indicates that the timestep is reduced below UC(12),
!        or cannot be reduced (either CT3=1, or we're at the first model).
!  COMMON H(:,:) - the stellar model from the previous timestep
!  COMMON DH(:,:) - Forced to 0 (changes in variables not known)
! ------------------------------------------------------------------------------
subroutine backup ( dty, jo )
   use mesh
   use mesh_enc
   use constants
   use settings
   use nucleosynthesis
   use test_variables
   use current_model_properties
   use binary_history, only: hpr
   use real_kind
   
   implicit none
   integer :: jo
   real(double) :: dty
   
   
   hspn = pprev(2:3)
   rlf = pprev(4:5)
   zet = pprev(6:7)
   xit = pprev(8:9)
   age = pprev(10)
   bm = pprev(11)
   mc = pprev(12:13)
   om = pprev(14)
   bper = pprev(15)
   sm = pprev(16)
   enc = pprev(17)
   tc = pprev(18:19)
   tfr = pprev(20)
   t0 = pprev(21)
   m0 = pprev(22)
   mta = pprev(23)
   om0 = pprev(24)
   omta = pprev(25)
   a0 = pprev(26)
   ata = pprev(27)
   e0 = pprev(28)
   eta = pprev(29)
   cdd = pprev(30)
   bp = pprev(31)
   horb = pprev(32)
   ro = pprev(33)
   ra2 = pprev(34)
   rs = pprev(35)
   secc = pprev(36)
   tn = pprev(37:38)
   wmh = pprev(39)
   wmhe = pprev(40)
   mh = pprev(41)
   mhe = pprev(42)
   mco = pprev(43)
   vmg = pprev(44)
   be = pprev(45:46)
   lh = pprev(47)
   lhe = pprev(48)
   lc = pprev(49)
   lnu = pprev(50)
   lth = pprev(51)
   mcb = pprev(52:59)
   msb = pprev(60:65)
   rcb = pprev(66:73)
   tct = pprev(74:81)

   ! *don't* backup DT!!
   jmod = jm2
   dt = ct3*dt
   dty = dt/csy
   if ( dt < uc(12) .or. jmod.le.2 .or. ct3 > 0.9999d0 ) jo = 2
   dh(:, 1:kh) = 0.0d0
   h(:, 1:kh) = hpr(:, 1:kh)
   call backup2
   jhold = -1

end subroutine backup



!> ------------------------------------------------------------------------------
!! NEXTDT
!!  Decide how large the next timestep will be, taking into account various
!!  factors, such as the values of CT1, CT2, CT3, the timestep control parameter
!!  CDD, the age limit UC(2) and the state of the binary (close to RLOF or
!!  not, whether star 2 is about to overtake star 1, etc.)
!! ------------------------------------------------------------------------------
!!
!! Input:
!!  DTY - The old timestep, in years
!!
!!       \todo FIXME: there are better ways to do the following:
!!  IT  - The FORTRAN unit number that init.dat was read from. Only used to
!!        decide whether we're in the middle of ZAHB construction (IT=24)
!!
!!  COMMON DH(:,:) - changes in independent variables; the ones we care
!!         about are listed in KJN(:), which is read from init.dat
!!
!! Output:
!!  DTY - The new timestep, in years
!!  JO  - Return code, 5 indicates that the age is now equal to or exceeds UC(2)
!! ------------------------------------------------------------------------------
!<

subroutine nextdt ( dty, jo, it )
   use mesh
   use compare_floats;
   use constants
   use settings
   use test_variables
   use current_model_properties
   use binary_history
   use real_kind
   
   implicit none
   integer :: jo,it
   real(double) :: dty

   integer :: ik,ij,kb,ka,jd,ips
   real(double) :: sum,fac,dty_min,dty_max,fdt,tdf
   
   
   ! Find change from last timestep:
   sum = 0.0d0
   do ik = 1, kh
      do ij = 1, kn
         sum = sum + dabs(dh(kjn(ij), ik))
      end do
   end do

   if ( sum == 0.0d0 ) then
      fac = 1.0d0
   else
      fac = sum/(cdd*kh)
      ! Specify next DTY; simple, for *1
      if ( jhold > 3 ) dty = dmax1( ct1, dmin1(ct2, 1.0d0/fac))*dty
      
      ! Limit DTY in such a way that we don't exceed the age specified in the
      ! input file (ignore for ZAHB construction).
      ! In this way, we get a model of *exactly* a specific age.
      ! Control the timestep so that the approach to the final age is a bit smooth
      if (it /= 24) then
         dty_min = uc(12)/csy
         dty_max = uc(2) - age
         if ( age+2*dty < uc(2) .and. age+3*dty >= uc(2)) then
            ! We expect three more timesteps, start constraining the timestep
            dty = 0.4*max(dty_max, dty_min )
         end if
         if ( age+dty < uc(2) .and. age+2*dty >= uc(2)) then
            ! We expect maybe two more timesteps, constrain
            dty = 0.6*max(dty_max, dty_min )
         end if
         if ( age+dty >= uc(2)) then             ! .and. age>uc(12)/csy) then
            ! This is our final timestep
            dty = min(dty_max, dty )
         end if
         if ( dty_max <= dty_min ) jo = 5
      end if
   end if

   if ( jb /= 1 ) then
      kb = 1
      if ( age > sect(kb) ) then
         ! For *2, find KB such that SECT(KB) < AGE <= SECT(KB + 1)
         do while ( .not. (age > sect(kb) .and. age.le.sect(kb + 1)) )
            kb = kb + 1
         end do
      end if

      ka = kb + 1
      if ( age + dty > sect(kb) ) then
      ! For *2, find KA such that SECT(KA - 1) < AGE + DTY <= SECT(KA)
         do
            if ( age + dty > sect(ka - 1) .and. age + dty.le.sect(ka) ) exit
            ka = ka + 1
         end do
      end if
      jd = 2.0d0 + 1.0d0/(fac + 0.1d0)
      ! If KA >> KB, pull back a bit
      if ( ka > kb + jd ) then
         ka = kb + jd
         dty = 0.5d0*(sect(ka) + sect(ka - 1)) - age
         dh(:, 1:kh) = 0.0d0
      end if
      if ( sect(ka) > 1.0d19 ) jo = 3
      write (2, 962) kb, ka, jd, age, dty, age + dty, sect(kb), sect(ka)
 962 format (3i5, 1p, 10d16.8)
   end if
   dt = csy*dty
   
   if ( jb == 1 .or. jo == 3 ) return

   ! For *2, some factors relating to accretion from *1. Ignored if this is *1
   if ( ka == 1 ) ka = 2
   ips = ka - 1
   t0 = csy*sect(ips)
   m0 = cmsn*ms(ips)
   fdt = sect(ka) - sect(ips)
   if ( fdt <= 0.0d0 ) jo = 11
   tdf = 0.0d0
   if ( fdt > 0.0d0 ) tdf = 1.0d0/(fdt*csy)
   mta = cmsn*(ms(ka) - ms(ips))*tdf
   om0 = cmsn*scm(ips)
   omta = cmsn*(scm(ka) - scm(ips))*tdf
   a0 = sang(ips)
   ata = (sang(ka) - sang(ips))*tdf
   e0 = se(ips)
   eta = (se(ka) - se(ips))*tdf
end subroutine nextdt

