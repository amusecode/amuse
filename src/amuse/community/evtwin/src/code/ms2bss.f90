!! ------------------------------------------------------------------------------
!  MS2BSS
!  Convert a stellar model for which we have entropy (actually, lnf/lnT) and
!  and composition profiles but not necessarily with a mesh spacing that we
!  can use directly and/or boundary conditions that match with this code
!  to a model that we can use, and that has the same total mass, entropy profile
!  and composition profiles.
!  In practice there will be small deviations from the input model, of course.
!  Originally intended to convert a Make Me A Star (MMAS) collision product
!  into a STARS input model for the study of blue straggler stars.
! ------------------------------------------------------------------------------
!  Input:
!   JOP  - The FORTRAN unit for the output file, for STAR12
!  Output:
!   JO3  - The return code, typically -1 if the model converged
! ------------------------------------------------------------------------------
subroutine ms2bss ( jop, jo3 )
   ! Convert a normal main sequence star to a different type of star, say
   ! a collision remnant.
   use real_kind
   use mesh
   use mesh_enc
   use control
   use constants
   use test_variables
   use current_model_properties
   use starting_values0
   use starting_values1
   use save_init
   use nucleosynthesis
   
   implicit none
   integer :: jop, jo3
   integer :: jc1, ksv, jip, flag
   integer :: k, j
   
   ! First stage: Evolve a normal main sequence star to the right mass and
   !  core hydrogen content. Read settings from the normal init.dat
   call read_target_model()
   uc(15) = get_iph(0.0d0, 5)
   rewind (jop)
   jip = jop
   call star12 ( jo3, jc1, jop, jip, ksv, 22 )
   rewind (22)
   rewind (jop)
   uc(15) = -1
   if ( jo3 == 51 ) write (jb, *) 'begin mutation'

   ! Termination code should be 51 (end of MS)
   if (jo3 /= 51) return

   ! Second stage: mutate the model by tweaking the composition profiles and
   !  puting an artificial contribution to luminosity in the energy equation
   ! Stop when we have matched the target radius

   ! Enable some flags that control various parts of the mutation process
   adj_mea = .true.
   adj_comp = .true.
   usemenc = .true.
   curr_diffsqr = 1.0d3;
   jip = 27 - jip
   jo3 = -1
   best_diffsqr = curr_diffsqr
   mutant_h(1:24, 1:kh) = h(1:24, 1:kh)
   rewind (jip)
   call star12 ( jo3, jc1, jop, jip, ksv, mutate_dat )
   rewind (mutate_dat)
   rewind (jop)
   rewind (jip)

   if (jo3 == 0 .or. jo3 == 2 .or. jo3 == 5 .or. jo3 == 12 .or. jo3 == 15) jo3 = 53

   if (jo3 == 53) then
      h(1:24, 1:kh) = mutant_h(1:24, 1:kh)
      dt = 1.0d1*csy
   end if

   age = 0
   flag = 0

   ! Set nucleosynthesis abundances, if we have them
   if (interpolation_has_nucabund()) then
      flag = 8
      if (allocated(Hnuc)) deallocate(Hnuc)
      allocate(Hnuc(2, nvar_nuc, kh))
      do k = 1, kh
         do j = 1, nvar_nuc
            Hnuc(1, j, k) = get_ipnucabund(H(4, k), j)
         end do
      end do
   end if

   ! Rotational period; set it as when it was read in from the target model
   ! Optionally replace with the value from init.run
   if (t0p >= 0.0) then
      h(13, 1:kh) = t0p
   else
      h(13, 1:kh) = t1p
   end if
   call output ( s_kpt, 14, 0, 0 )
   call output ( s_kpt, 13, 0, 0 )
   rewind (14)
   rewind (13)
   call output ( s_kpt, 65, 0, flag )
   close (65)
   if (interpolation_has_nucabund()) deallocate(Hnuc)

   ! Termination code should be 53 (Found convergence minumum)
   if (jo3 /= 53) return

   write (jb, *) 'mutation done', best_mod, best_diffsqr
   ! Third stage: evolve the remnant
   ! Let the normal code handle this
   jo3 = -1
   mutate = .false.
   adj_mea = .false.
   adj_comp = .false.
   usemenc = .false.

end subroutine ms2bss


