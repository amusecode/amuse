
!> Load model number IM from the file IP
subroutine load_star_model(ip,im, h, dh, hnuc, sm,dty,age,per,bms,ecc,p,enc,kh,kp,jmod,jb,jn,jf)
   use real_kind
   use mesh, only:nm, nvar
   use indices
   
   implicit none
   integer, intent(in) :: ip, im
   integer, intent(out) :: kh,jn,jf, jmod,kp,jb
   real(double), intent(out) :: sm,dty,age,per,bms,ecc,p,enc
   real(double), intent(out) :: h(nvar,nm), dh(nvar,nm), hnuc(50, nm)
   real(double) :: var(nvar)
   character(len=500) :: ip_name
   integer :: ioerror
   integer :: i, j, k
   
   do i = 1, im
      ! Skip all models upto IM
      read (ip, *, iostat=ioerror) sm, dty, age, per, bms, ecc, p, enc, kh, kp, jmod, jb, jn, jf
      if ( ioerror /= 0 ) then
         inquire(unit=ip, name=ip_name)
         write(0, *) 'error reading model', im, 'from "', trim(ip_name),'"'
         stop
      end if
      do k = 1, kh
         read (ip, *) var(1:jn)
         do j = 1, jn
            h(var_input_order(j), k) = var(j)
         end do
      end do
      if (iand(jf, 4) == 4) then
         do k = 1, kh
            read (ip, *) var(1:jn)
            do j = 1, jn
               dh(var_input_order(j), k) = var(j)
            end do
         end do
      end if
      if (iand(jf, 8) == 8) then
         do k = 1, kh
            read (ip, *) hnuc(1:50, k)
         end do
      end if
   end do
   
   if (.not. (iand(jf, 4) == 4)) dh(:,:) = 0.0d0
   if (.not. (iand(jf, 8) == 8)) hnuc(:,:) = 0.0d0
   
   rewind (ip) 
end subroutine load_star_model

