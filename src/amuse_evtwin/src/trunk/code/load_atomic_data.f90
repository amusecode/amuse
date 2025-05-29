! Read nuclear reaction (QRT) and neutrino (QNT) Q values, in MeV; constants
! for electron screening (CZA, CZB, CZC, CZD, VZ); atomic parameters (CBN, KZN),
! with masses (CAN) consistent with Q-values; ionization potentials (CHI) and
! statistical weights (COM); molecular hydrogen parameters (CH2)
subroutine load_atomic_data(file)
   use real_kind
   use reaction_rate_data
   use atomic_data
   
   implicit none
   integer, intent(in) :: file
   integer z1(92), z2(92), j
   real(double) :: cxa,cxb,cxc,cxd

   data z1 /1,2,2,4,4,6,7,8,4,6,7,8,10,6,8,8,10,12,0,0,1,1,1,0,1,1,1,  &
        2,2,2,2,1,1,2,1,1,2,2,1,1,2,1,2,2,1,2,2,0,0,1,1,1,1,1,4,1,1,1,  &
        4,4,4,1,1,1,1,4,4,4,0,0,0,1,0,0,0,1,0,0,0,1,1,1,4,4,4,4,1,1,1,  &
        1,1,1/
   data z2 /1,2,2,0,1,1,1,1,2,2,2,2, 2,6,6,8, 0, 0,0,0,1,3,4,4,6,7,8,  &
        6,7,8,3,5,6,6,8,8,8,8,9,9,9,10,10,10,10,10,10,11,11,11,11,11,  &
        11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,  &
        13,13,13,13,13,13,13,13,13,13,13,11,11,11,14,14,14,7,6,10/

992 format (1p, 12(10d12.4,/), 32(9d12.4,/), 4d12.4,/, 9i12)


   read (file,992)   &
        qrt, qnt, cza(1:20), czb(1:20), czc(1:20), czd(1:20),  &
        vz, cbn, can, com, chi, ch2, kzn
   close (file)

   !     Compute screening factors: mainly because not all of those that are
   !     needed for the nucleosynthesis code are in the data file.
   cxa = 5.0/3.0
   cxc = cxa - 1.0
   cxb = 2.0*cxc
   cxd = 1.86
   do j = 1, 92
      cza(j) = (z1(j)+z2(j))**cxa - z1(j)**cxa - z2(j)**cxa
      czb(j) = (z1(j)+z2(j))**cxb - z1(j)**cxb - z2(j)**cxb
      czc(j) = -((z1(j)+z2(j))**cxc - z1(j)**cxc - z2(j)**cxc)
      czd(j) = (z1(j)+z2(j))**cxd - z1(j)**cxd - z2(j)**cxd
   end do

   ! We need the charges KZN in floating point operations a few times
   ! Do the int->float conversion once on startup to safe a few cycles
   dkzn(:) = kzn(:)

   ! Log of A**1.5, used several times in the EoS because it is a factor
   ! that appears in the free energy
   lcan = 1.5d0*log(can)
end subroutine load_atomic_data

