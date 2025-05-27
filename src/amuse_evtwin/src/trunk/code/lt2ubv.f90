module ltm2ubv

   use real_kind
   
   implicit none
   real(double) :: tgr(34), ggr(13), tab(34,13,5)

contains

subroutine load_colour_conversion_table(fu)
   use real_kind
   
   implicit none
   integer, intent(in) :: fu
   integer :: i,j, k

991 format (3(10f10.5,/), 4f10.5,/, 10f10.5,/, 3f10.5,/, 442(5f8.3,/))
   read (fu,991) tgr, ggr, (((tab(i,j,k), k=1,5), j=1,13), i=1,34)
   close (fu)
end subroutine


subroutine lt2ubv ( logl, logt, mass, mv, bminv, uminb )
   ! yields values of MV, B-V and U-B for given log L, log T, and mass
   
   use real_kind
   
   implicit none
   real(double) :: logl, logt, mass, mv, bminv, uminb
   
   integer :: i,ing1,ing2,int1,int2,k,j,k0,k1,k2
   real(double) :: logm,logg,logg1,logg2,logt1,logt2
   real(double) :: dg1,dg2,bc1,bc2,ub1,ub2,bv1,bv2,dt1,dt2,ubx,bvx,bcx,mbol
   real(double) :: bc(4), ub(4), bv(4), vr(4), ri(4)
   
   logm = log10(mass)
   logg = logm + 4.0d0*logt - logl - 10.6071d0

   ! determine values of log g to interpolate between
   ing1 = 1
   ing2 = 13
1  continue                         ! FIXME - replace with do while()
   if ( ing2 - ing1.gt.1 ) then     ! FIXME - replace with do while()
      i = (ing1 + ing2)/2
      if ( ggr(i).gt.logg ) then
         ing2 = i
      else
         ing1 = i
      end if
      goto 1                        ! FIXME - replace with end do
   end if                           ! FIXME - replace with end do
   logg1 = ggr(ing1)
   logg2 = ggr(ing2)
   
   ! determine values of log T to interpolate between
   int1 = 1
   int2 = 34
2  continue                         ! FIXME - replace with do while()
   if ( int2 - int1.gt.1 ) then     ! FIXME - replace with do while()
      i = (int1 + int2)/2
      if ( tgr(i).gt.logt ) then
         int2 = i
      else
         int1 = i
      end if
      goto 2                        ! FIXME - replace with end do
   end if                           ! FIXME - replace with end do
   
   logt1 = tgr(int1)
   logt2 = tgr(int2)
   
   do k = 1, 2
      do j = 1, 2
         k0 = (k - 1)*2 + j
         k1 = int1 - 1 + k
         k2 = ing1 - 1 + j
         bc(k0) = tab(k1, k2, 1)
         ub(k0) = tab(k1, k2, 2)
         bv(k0) = tab(k1, k2, 3)
         vr(k0) = tab(k1, k2, 4)
         ri(k0) = tab(k1, k2, 5)
      end do
   end do
   
   dg1 = (logg - logg1)/(logg2 - logg1)
   dg2 = 1.0d0 - dg1
   bc1 = bc(2)*dg1 + bc(1)*dg2
   ub1 = ub(2)*dg1 + ub(1)*dg2
   bv1 = bv(2)*dg1 + bv(1)*dg2
   bc2 = bc(4)*dg1 + bc(3)*dg2
   ub2 = ub(4)*dg1 + ub(3)*dg2
   bv2 = bv(4)*dg1 + bv(3)*dg2
   dt1 = (logt - logt1)/(logt2 - logt1)
   dt2 = 1.0d0 - dt1
   bcx = bc2*dt1 + bc1*dt2
   ubx = ub2*dt1 + ub1*dt2
   bvx = bv2*dt1 + bv1*dt2
   mbol = 4.75d0 - 2.5d0*logl
   mv = mbol - bcx
   bminv = bvx
   uminb = ubx
   
end subroutine lt2ubv
      
end module
