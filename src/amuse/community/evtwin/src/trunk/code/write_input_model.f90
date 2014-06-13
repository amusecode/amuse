! ------------------------------------------------------------------------------
! OUTPUT
!  Write a stellar model in a format that is suitable for restarting the
!  evolution code.
!  Also stores some information at the beginning of the ZAHB construction
! (JO=8) that is restored when ZAHB consstruction is done (JO=13).
! During ZAHB construction no input models are output.
! ------------------------------------------------------------------------------
!  Input:
!     KP    - The total number of timesteps for the current evolution run.
!             Written to the input file, but is that really useful (this
!             parameter is normally read from init.run).
!     JOP   - The FORTRAN file unit where the output is written to. Unless
!             JOP=13 or 14, in which case it's written to 14 or 13, unless
!             JO==13, in which case it really is written to 13 or 14. ??!!?!??
!     JO    - Used to decide whether the code is starting (8) or stopping (13)
!             ZAHB construction. Properties are then stored for JO=8 and
!             written out for JO=13.
!     JF    - A flag field, passed directly to the output file.
! ------------------------------------------------------------------------------

subroutine output ( kp, jop, jo, jf )
   ! Write an intermediate or final model to disc
   use real_kind
   use mesh
   use mesh_enc
   use constants
   use test_variables
   use current_model_properties
   use nucleosynthesis, only: Hnuc
   use indices
   
   implicit none
   integer :: kp,jop,jo,jf

   integer :: jp,ii,ij,ik
   real(double) :: smout,aper
   real(double) :: var(nvar)

   ! Save some parameters so that they can be restored after the He flash
   real(double), save :: aj
   integer, save :: jmd, kp_bck

   integer, parameter :: i1(2) = (/ 0, 24/)  ! Offset for *1, *2
   integer, parameter :: i2(2) = (/24, 16/)  ! Number of variables for *1, *2


92 format (1x, 1p, 40d23.15, 0p)
95 format (1x, 1p, 8d23.15, 0p, 6i6)


   if ( jo == 8 ) then        ! He flash
      
      ! Backup age, model number and number of models, to be restored after the flash
      aj = age
      jmd = jmod
      kp_bck = kp
      return   ! Comment this line out to save post-He-flash model.
   else if ( jo == 13 ) then
      
      ! Restore age and model number, increase the timestep
      age = aj
      jmod = jmd
      dt = 5.0d4*csy

      ! The number of remaining models to calculate. Do at least a few.
      kp = max(kp_bck - jmod, 10)
      if (kp_bck < 0) kp = kp_bck
   end if

   if ( mutate .and. jmod>0 ) dt = min(dt, 1.0d3*csy)
   jp = jop
   if ( ( jop == 13.or.jop == 14 ) .and. jo /= 13 ) jp = 27 - jop
   do ii = 1, ktw
      smout = h(idx_for_star(VAR_MASS, ii), 1)/cmsn
      aper = 2.0*cpi/(abs(h(idx_for_star(VAR_OMEGA, ii), 1)) * csday)
      write (jp, 95) smout, dt/csy, age, bper, bm/cmsn, secc, aper, enc, kh, kp, jmod, jb, i2(ii), jf
      do ik = 1, kh        ! Model data
         aper = 2.0*cpi/(abs(h(idx_for_star(VAR_OMEGA, ii), ik)) * csday)
         do ij = 1, i2(ii)
            var(ij) = h(idx_for_star(var_input_order(ij), ii), ik)
         end do
         var(13) = aper
         write (jp, 92) var(1:i2(ii))
      end do
      if (iand(jf, 4) == 4) then
         do ik = 1, kh     ! Derivative data
            do ij = 1, i2(ii)
               var(ij) = dh(idx_for_star(var_input_order(ij), ii), ik)
            end do
            write (jp, 92) var(1:i2(ii))
         end do
      end if
      if (iand(jf, 8) == 8) then
         do ik = 1, kh     ! Nucleosynthesis data
            write (jp, '(1x, 1p, 50d23.15, 0p)') (Hnuc(ii, ij, ik), ij = 1, 50)
         end do
      end if
   end do

   flush(jp)

end subroutine output

