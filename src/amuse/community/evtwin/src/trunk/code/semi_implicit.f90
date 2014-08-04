! ------------------------------------------------------------------------------
!  COMP_SEMI_IMPLICIT_QUANTITIES
!   computes the values of quantities that are updated in between
!   iterations of the solver (hence semi-implicit).
! ------------------------------------------------------------------------------
!  Input:
!     //          - H and DH
!  Output:
!     ddXn        - Second order corrections to diffusion equation of species n
!  Returns:
!     true if the implicit quantities were updated and the iteration can be
!     considered "converged", or "false" if more iterations are needed (for
!     instance because the change in implicit functions is still large)
! ------------------------------------------------------------------------------
logical function comp_semi_implicit_quantities()
   use real_kind
   use mesh
   use semi_implicit_variables
   use control
   use interpolate
   use distortion
   use constants
   use settings
   use indices
   use roche
   use structure_functions
   use eostate_types
   implicit none
   ! Local variables:
   type(eostate) :: eos
   real(double) :: var(Nvar)
   real(double) :: rho(KH, 2), p(KH, 2), rr(KH, 2), m(KH, 2), cs(KH, 2)
   real(double) :: r0(KH), sigma_rot(KH), mu_rot(KH), sigma(KH,2:4), mu(KH,2:4)
   real(double) :: avgg(KH), invg(KH), spsi(KH)
   real(double) :: fu2(2,KH), fv(2, KH), omega(KH, 2)
   real(double) :: fn(NFUNC)
   real(double) :: a_orb, oa, mp(2), phi_l1, phi_l2, RL, RL2, phi_scale, po, q, rossby, thetaf
   real(double) :: delta_phi(2)
   integer :: k, kk, Jstar

   real(double) :: old_fp_star(NM, 2)
   real(double) :: old_ft_star(NM, 2)
   real(double) :: old_Vmc2_star(NM, 2)
   real(double) :: old_diff_omega(NM, 2)
   real(double) :: max_Vmc2, max_domega, epsilon
   logical, parameter :: monitor_convergence = .false.

   real(double) :: phi(NM, 2), ent(KH, 2), hp(KH, 2)
   real(double) :: mdot(2), old_mdot(2)
   real(double) :: ent_l1(2), cs_l1(2), vm, xx, fac, fac1, fac2, rrho, oenth, ell
   real(double) :: m1, m2
   type(interpolate_t) :: enthalpy(2), sound_speed, pressure(2), density(2)

   comp_semi_implicit_quantities = .true.

   var(:) = H(:, 1) + dh(:, 1)
   oa = var(VAR_HORB)
   mp(1) = var(VAR_PMASS)
   mp(2) = var(VAR_BMASS) - mp(1)
   a_orb = oa*oa*(mp(1) + mp(2))/(mp(1)*mp(2)*cg1)**2
   ! Potential at L1, L2
   m1 = max(mp(1), mp(2))
   m2 = min(mp(1), mp(2))
   phi_l1 = calc_scaled_potential(m1, m2, calc_scaled_xl1(m1, m2), 0.0d0, 0.0d0)
   phi_l2 = calc_scaled_potential(m1, m2, calc_scaled_xl2(m1, m2), 0.0d0, 0.0d0)
   !if (ktw == 1) mp(2) = 0.0d0

   !print *, calc_scaled_xl1(mp(1), mp(2)), calc_scaled_xl1(mp(2), mp(1))
   !print *, calc_scaled_potential(mp(1), mp(2), calc_scaled_xl1(mp(1), mp(2)), 0.0d0, 0.0d0)
   !print *, calc_scaled_potential(mp(2), mp(1), calc_scaled_xl1(mp(2), mp(1)), 0.0d0, 0.0d0)
   !print *, phi_l1
   !stop

   old_mdot = mdot_rlof
   if (monitor_convergence) then
      old_fp_star = fp_star
      old_ft_star = ft_star
      old_Vmc2_star = Vmc2_star
      old_diff_omega = diff_omega

      max_Vmc2 = 1.0d0
      max_domega = 1.0d0
   end if

   ! Break out now if we don't need to calculate detailed properties at each grid point
   if (.not. use_contact_flow .and. .not. use_clairaut_distortion) return

   ! Stellar properties at each gridpoint
   do k=1, kh
      var(:) = H(:, k) + dh(:, k)
      call funcs1(k, 0, var(:), DH(:, k), fn, eos)
      !print *, k, (fn(FN_PHI_ROCHE) - phi_l1) * CG*(mp(1)+mp(2))/a_orb*1.0d22, fn(FN_ENT)
      do Jstar = 1, ktw
         phi(k, Jstar) = fn(fn_idx_for_star(FN_PHI_ROCHE, Jstar))
         ent(k, Jstar) = fn(fn_idx_for_star(FN_ENT, Jstar))
         fv(Jstar, k)  = fn(fn_idx_for_star(FN_FV, Jstar))
         fu2(Jstar, k) = fn(fn_idx_for_star(FN_FU2K, Jstar))
         rho(k, Jstar) = fn(fn_idx_for_star(FN_RHO, Jstar))
         p(k, Jstar)   = fn(fn_idx_for_star(FN_P, Jstar))
         cs(k, Jstar)  = eos%gamma1 * eos%p / eos%rho
         omega(k, Jstar) = H(idx_for_star(VAR_OMEGA, Jstar), k)
         m(k, Jstar)  = H(idx_for_star(VAR_MASS, Jstar), k)
         rr(k, Jstar) = sqrt(max(exp(2.0d0 * H(idx_for_star(VAR_LNR, Jstar), k)) - CT(8), 0.0d0))
         hp(k, Jstar) = fn(fn_idx_for_star(FN_HP, Jstar))
      end do
   end do

   ! Horizontal component of meridional circulation
   do k=2, kh-1
      do Jstar=1, ktw
         Vmc2_star(k, Jstar) = fv(Jstar,k) * (fu2(Jstar, k) - fu2(Jstar, k+1))
         diff_omega(k, Jstar) = 0.5d0 * (omega(k+1, Jstar) - omega(k-1, Jstar))
      end do
   end do
   Vmc2_star(1, :) = 0.0d0
   Vmc2_star(kh, :) = 0.0d0
   diff_omega(1, Jstar) = (omega(2, Jstar) - omega(1, Jstar))
   diff_omega(kh, Jstar) = (omega(kh, Jstar) - omega(kh-1, Jstar))

   ! Calculate the shape of the stars, from the multipole expansion of the gravitational potential.
   if (use_clairaut_distortion) then
      ft_star(1:kh,1:KTW) = 1.0
      fp_star(1:kh,1:KTW) = 1.0
      do Jstar=1, ktw
         m(kh, Jstar) = m(kh-1, Jstar)      ! Avoid singularity at m = 0
         rr(kh, Jstar) = rr(kh-1, Jstar)    ! Avoid singularity at r = 0
         call calculate_distortion_coefficients(kh, rho(:,Jstar), rr(:,Jstar), m(:,Jstar), omega(:,Jstar), mp(3-Jstar), a_orb,&
                                                   r0, a_rot2(1:kh,Jstar), a_tide(1:kh,2:4,Jstar),                   &
                                                   sigma_rot, mu_rot, sigma, mu)
         call calculate_effective_gravity(kh, rr(:,Jstar), m(:,Jstar), omega(:,Jstar), mp(3-Jstar), a_orb,           &
                                             r0, a_rot2(1:kh,Jstar), a_tide(1:kh,2:4,Jstar),                         &
                                             sigma_rot, mu_rot, sigma, mu,                                           &
                                             avgg, invg, spsi, fp_star(1:kh,Jstar), ft_star(1:kh,Jstar))

         r0_over_rv(1:KH-1, Jstar) = r0(1:KH-1) / rr(1:KH-1, Jstar)
         r0_over_rv(KH, Jstar) = 1.0d0
      end do
   else
      r0_over_rv = 1.0d0
      a_rot2 = 0.0d0
      a_tide = 0.0d0
   end if

   ! Calculate mass transfer, both semi-detached and in contact
   ! We obviously need to be doing a binary in TWIN mode for this to make sense
   mdot_rlof = 0.0d0
   if (ktw == 2 .and. use_contact_flow .and. (phi(1,1) <= phi_l1 .or. phi(1,2) <= phi_l1)) then
   call make_interpolation_table(kh, phi(:, 1), ent(:, 1), enthalpy(1))
   call make_interpolation_table(kh, phi(:, 2), ent(:, 2), enthalpy(2))
   call make_interpolation_table(kh, phi(:, 1), cs(:, 1), sound_speed)
   call make_interpolation_table(kh, phi(:, 1), p(:, 1), pressure(1))
   call make_interpolation_table(kh, phi(:, 2), p(:, 2), pressure(2))
   call make_interpolation_table(kh, phi(:, 1), rho(:, 1), density(1))
   call make_interpolation_table(kh, phi(:, 2), rho(:, 2), density(2))
   ent_l1 = 0.0d0
   cs_l1 = 0.0d0

   ! Calculate sound speed and enthalpy in L1
   if (phi(1, 1) <= phi_l1) then
      ent_l1(1) = evaluate_interpolation_table(phi_l1, enthalpy(1))
      cs_l1(1)  = evaluate_interpolation_table(phi_l1, sound_speed)
   end if
   if (phi(1, 2) <= phi_l1) then
      call update_interpolation_table(kh, phi(:, 2), cs(:, 2), sound_speed)
      ent_l1(2) = evaluate_interpolation_table(phi_l1, enthalpy(2))
      cs_l1(2)  = evaluate_interpolation_table(phi_l1, sound_speed)
   end if

   ! The scale factor for the potential function, give results in CGS units (erg/g)
   phi_scale = CG*(mp(1)+mp(2))/a_orb*1.0d22

   ! Mass ratio, in the sense of (least massive) / (most massive)
   q = mp(2) / mp(1)
   if (q > 1.0d0) q = mp(1) / mp(2)

   ! Calculate the mass and energy flows in both components
   mdot = 0.0d0
   do Jstar=1, ktw
      do k=1, kh
         ! Test whether this grid point overflows the critical lobe
         if (phi(k, Jstar) <= phi_l1) then

            ! Test whether we are within the shared envelope or not. This changes the calculation of the flow velocity:
            ! if we are not in the shared envelope, then we have a free expansion through the surface through L1 and L1 is a sonic
            ! point in the flow. If we are in the shared envelope the velocity field is determined by the continuity equation.
            ! The direction of the flow is determined by the horizontal pressure gradient.
            ! In either case we use Bernoulli's equation to calculate the velocity at the current point.
            po = 0.0d0
            vm = 0.0d0
            if (phi(k, Jstar) >= phi(1, 3-Jstar)) then   ! Contact
               po    = evaluate_interpolation_table(phi(k,Jstar), pressure(3-Jstar))
               rrho  = evaluate_interpolation_table(phi(k,Jstar), density(3-Jstar))
               oenth = evaluate_interpolation_table(phi(k,Jstar), enthalpy(3-Jstar))
               vm = rrho**2 * 2.d0 * abs((ent(k, Jstar) - oenth) / (rrho**2 - rho(k, Jstar)**2))
            else                                         ! No contact
               vm = cs(k, Jstar)
            end if

            ! Because of the distortion near L1, the density at the isobars is different (much lower) than the density far from
            ! L1. The following is an approximation to the expansion factor based on fits to Roche geometry; this depends on the
            ! mass ratio and the (dimensionless) value of the potential excess.
            ! Because we use different fits for the primary (most massive) and secondary (least massive) component that are not
            ! exactly equal at q=1 we do a smooth interpolation between the two near q=1.
            xx = (phi(k,Jstar) - phi_l1) / (phi_l2 - phi_l1)      ! Dimensionless coordinate of potential
            if (mp(Jstar) > mp(3-Jstar)) then   ! Currently most massive component
               if (q > 0.95) then               ! Nearly equal masses, do a smooth transition between the two interpolations
                  fac1 = stream_cross_section_primary(q, xx)
                  fac2 = 0.5d0 * (fac1 + stream_cross_section_secondary(q, xx))
                  fac  = fac1 + (q - 0.95) * (fac2 - fac1) / 0.05
                  !fac = 0.5 * ((2.0d0 - q) * stream_cross_section_primary(q, xx) + q * stream_cross_section_secondary(q, xx))
               else
                  fac = stream_cross_section_primary(q, xx)
               end if
            else
               if (q > 0.95) then               ! Nearly equal masses, do a smooth transition between the two interpolations
                  fac1 = stream_cross_section_secondary(q, xx)
                  fac2 = 0.5d0 * (stream_cross_section_primary(q, xx) + fac1)
                  fac  = fac1 + (q - 0.95) * (fac2 - fac1) / 0.05
                  !fac = 0.5 * (q * stream_cross_section_primary(q, xx) + (2.0d0 - q) * stream_cross_section_secondary(q, xx))
               else
                  fac = stream_cross_section_secondary(q, xx)
               end if
            end if

            ! Mass flow follows the pressure gradient
            Rossby = 0.0d0
            thetaf = 0.0d0
            if (vm > 0.0d0 .and. p(k,Jstar) > po) then
               vm = sqrt(vm)
               mdot(Jstar) = mdot(Jstar) + fac*rho(k, Jstar) * vm * a_orb**2 * (phi(k, Jstar) - phi(k+1, Jstar)) * 1.0d-11

               ! Calculate the Rosby number associated with the flow; this determines the opening angle of the tidal stream
               ! The scale-height is set to the pressure scale hight here. An alternative is the size of the star, which is
               ! typically much larger and so produces a much smaller Rossby number.
               ell = hp(k, Jstar)
               !ell = 1.0d11 * rr(k, Jstar)
               Rossby = vm / (2.d0 * omega(k, Jstar) * ell)

               ! Width (in radians) of the geostrophic flow in a contact system
               thetaf = acos(Rossby / (Rossby + 1.0d0))

               ! Calculate coefficients in the energy transport equation
            end if

            !write (*, '(1x,1p,I3,I3,10E18.10)') Jstar, k, xx, sqrt(vm), sqrt(cs(k, Jstar)), sqrt(vm/cs(k, Jstar)), p(k, Jstar), &
            !   po, rho(k, Jstar), mdot(Jstar) * CSY/CMSN, Rossby, thetaf
         end if
      end do
   end do

   ! Translate the energy flow in a sequence of source/sink terms

   ! Calculate the with and temperature of the tidal stream

   !RL = a_orb * rlobe((mp(1)/mp(2))**(1./3.))
   !RL2 = a_orb * rlobe((mp(2)/mp(1))**(1./3.))
   !print *, sqrt(cs_l1(1)), sqrt(cs_l1(2))
   !print *, mp(1)/ mp(2), RL, rr(1,1)
   !print *, mp(1)/ mp(2), RL2, rr(1,2)
   !RL = (max(rr(1,1)-RL, 0.0d0)/RL)**3
   !print *, mdot(1) * CSY / CMSN,&
   !   RL * rho(1,1) / (5./3. * p(1, 1) / rho(1, 1))**1.5*(CG*1.0d33*(mp(1)+mp(2)))**2 / (CMSN*1.0d33)*CSY,&
   !   RL*mp(1)/sqrt(CG * rho(1,1)) * CSY / CMSN*1.0d-11, RL
   !print *, 1, phi_l1, phi(1, 1)!, phi_l2
   !print *, 2, phi_l1, phi(1, 2)!, phi_l2
   !if (phi(1, 1) <= phi_l1 .and. phi(1, 2) <= phi_l1) print *, 'Contact'
   !print *, ''

   if (monitor_convergence) then
      do k=1, kh
         do Jstar=1, ktw
            if (abs(diff_omega(k, Jstar)) > max_domega) max_domega = abs(diff_omega(k, Jstar))
            if (abs(Vmc2_star(k, Jstar)) > max_Vmc2) max_Vmc2 = abs(Vmc2_star(k, Jstar))
         end do
      end do

      do k=1, kh
         do Jstar=1, ktw
            epsilon = epsilon + abs(diff_omega(k, Jstar) - old_diff_omega(k, Jstar)) / max_domega
            epsilon = epsilon + abs(Vmc2_star(k, Jstar) - Vmc2_star(k, Jstar)) / max_Vmc2
            epsilon = epsilon + abs(fp_star(k, Jstar) - old_fp_star(k, Jstar))
            epsilon = epsilon + abs(ft_star(k, Jstar) - old_ft_star(k, Jstar))
         end do
      end do
      epsilon = epsilon / (ktw * kh * 4)

      write (1, *) 'Explicit epsilon = ', log10(epsilon + 1.0d-16)
   end if 

   ! Store actual mass loss rate
   mdot_rlof0 = mdot

   ! Extract coefficient for scaling relation
   delta_phi(1) = max(phi_l1 - phi(1, 1), 0.d0)
   delta_phi(2) = max(phi_l1 - phi(1, 2), 0.d0)
   where (delta_phi > 0.0d0) mdot = mdot / delta_phi**2.5

   epsilon = abs(mdot(1) - old_mdot(1))/(abs(old_mdot(1))+1.0d-32) + abs(mdot(2) - old_mdot(2))/(abs(old_mdot(2))+1.0d-32)
   write (1, *) 'Explicit epsilon = ', log10(epsilon + 1.0d-16)
   write (1, '(1x,1p,4e12.4)') mdot(1) * CSY / CMSN, mdot(2) * CSY / CMSN, old_mdot(1) * CSY / CMSN, old_mdot(2) * CSY / CMSN
   !if (abs(mdot(1) - mdot_rlof(1)) > 0.05d0 * mdot_rlof(1)) then
   if (epsilon > 1.0d-2 .and. abs(mdot(1)) > abs(old_mdot(1))) then
      if (abs(old_mdot(1)) > 0.0d0) mdot = 0.5d0 * mdot + 0.5 * old_mdot
      comp_semi_implicit_quantities = .false.
   end if
   if (epsilon < 1.0d-4) mdot = old_mdot ! Hack: don't update if accurate enough
   write (1, *) comp_semi_implicit_quantities

   call destroy_interpolation_table(sound_speed)
   call destroy_interpolation_table(enthalpy(1))
   call destroy_interpolation_table(enthalpy(2))
   call destroy_interpolation_table(density(1))
   call destroy_interpolation_table(density(2))
   call destroy_interpolation_table(pressure(1))
   call destroy_interpolation_table(pressure(2))
   mdot_rlof = mdot
   end if

contains

   ! ------------------------------------------------------------------------------
   ! FRACSTEP
   ! ------------------------------------------------------------------------------
   ! Input:
   !  F - Numerator of fraction to be caculated
   !  G - Denominator of fraction to be calculated
   ! Returns:
   !  F/G if F>0, 0 otherwise
   ! ------------------------------------------------------------------------------
   pure function fracstep (f, g)
      use real_kind

      implicit none
      real(double), intent(in) :: f, g
      real(double) :: fracstep

      fracstep = 0.0d0
      if (f > 0.0d0) fracstep = f / g
   end function fracstep

end function comp_semi_implicit_quantities



