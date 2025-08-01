! -*- Mode: Fortran90; tab-width: 3 -*-
#include "assert.h"

subroutine solver ( iter, ig, kt5, jo )
   use real_kind
   use mesh
   use settings
   use solver_global_variables
   use nucleosynthesis
   use control
   use current_model_properties
   use indices
   use control

   implicit none

   ! Arguments:
   integer, intent(in) :: iter, ig(130), kt5
   integer, intent(inout) :: jo

   ! External functions:
   real(double), external :: linesearch
   logical, external :: comp_semi_implicit_quantities

   ! Local variables:
   integer :: i, ij, ik, ikm(NVAR), ikd(NVAR), k
   real(double) :: err, errpr, ert(NVAR), fac, frr, frrpr, vx, vz
   real(double) :: max_errpr

   logical :: converged, converged_ok

   logical :: updated_explicit
   logical :: updated_once

   real(double) :: outres, besterr
   real(double) :: besth(NVAR,kh), bestdh(NVAR,kh)

   character :: varnames(1:NVAR)*(7),nucname(47)*(5)

   data nucname/'  g  ','  n  ','  D  ',' He3 ',' Li7 ',' Be7 ',' B11 ',&
       ' C13 ',' C14 ',' N15 ',' O17 ',' O18 ',' F19 ',' Ne21',&
       ' Ne22',' Na22',' Na23',' Mg24',' Mg25',' Mg26','Al26M',&
       'Al26G',' Al27',' Si28',' Si29',' Si30',' P31 ',' S32 ',&
       ' S33 ',' S34 ',' Fe56',' Fe57',' Fe58',' Fe59',' Fe60',' Co59',&
       ' Ni58',' Ni59',' Ni60',' Ni61', ' H1  ',' He4 ',' C12 ',&
       ' N14 ',' O16 ','Ne20 ',' Ca40'/

   ! Set the names of independent variables dynamically.
   do i=1, nvar
      write(varnames(i), '(1x,1p,I4)') i
   end do

   varnames(idx_for_star(VAR_LNF,   1)) = '   ln f'
   varnames(idx_for_star(VAR_LNT,   1)) = '   ln T'
   varnames(idx_for_star(VAR_O16,   1)) = '    X16'
   varnames(idx_for_star(VAR_MASS,  1)) = '      M'
   varnames(idx_for_star(VAR_H1,    1)) = '     X1'
   varnames(idx_for_star(VAR_QK,    1)) = '      C'
   varnames(idx_for_star(VAR_LNR,   1)) = '   ln r'
   varnames(idx_for_star(VAR_LUM,   1)) = '      L'
   varnames(idx_for_star(VAR_HE4,   1)) = '     X4'
   varnames(idx_for_star(VAR_C12,   1)) = '    X12'
   varnames(idx_for_star(VAR_NE20,  1)) = '    X20'
   varnames(idx_for_star(VAR_INERT, 1)) = '      I'
   varnames(idx_for_star(VAR_OMEGA, 1)) = '   Prot'
   varnames(idx_for_star(VAR_PHI,   1)) = '    phi'
   varnames(idx_for_star(VAR_PHIS,  1)) = '  phi_s'
   varnames(idx_for_star(VAR_N14,   1)) = '    X14'
   varnames(idx_for_star(VAR_MG24,  1)) = '    X24'
   varnames(idx_for_star(VAR_SI28,  1)) = '    X28'
   varnames(idx_for_star(VAR_FE56,  1)) = '    X56'
   varnames(idx_for_star(VAR_TAM,   1)) = '  Hspin'

   varnames(idx_for_star(VAR_LNF,   2)) = '*2:ln f'
   varnames(idx_for_star(VAR_LNT,   2)) = '*2:ln T'
   varnames(idx_for_star(VAR_O16,   2)) = '*2: X16'
   varnames(idx_for_star(VAR_MASS,  2)) = '*2:   M'
   varnames(idx_for_star(VAR_H1,    2)) = '*2:  X2'
   varnames(idx_for_star(VAR_QK,    2)) = '*2:   C'
   varnames(idx_for_star(VAR_LNR,   2)) = '*2:ln r'
   varnames(idx_for_star(VAR_LUM,   2)) = '*2:   L'
   varnames(idx_for_star(VAR_HE4,   2)) = '*2:  X4'
   varnames(idx_for_star(VAR_C12,   2)) = '*2: X22'
   varnames(idx_for_star(VAR_NE20,  2)) = '*2: X20'
   varnames(idx_for_star(VAR_INERT, 2)) = '*2:   I'
   varnames(idx_for_star(VAR_OMEGA, 2)) = '*2:Prot'
   varnames(idx_for_star(VAR_PHI,   2)) = '*2: phi'
   varnames(idx_for_star(VAR_PHIS,  2)) = '*2:phis'
   varnames(idx_for_star(VAR_N14,   2)) = '*2: X24'
   varnames(idx_for_star(VAR_MG24,  2)) = '*2: X24'
   varnames(idx_for_star(VAR_SI28,  2)) = '*2: X28'
   varnames(idx_for_star(VAR_FE56,  2)) = '*2: X56'
   varnames(idx_for_star(VAR_TAM,   2)) = '*2:Hspn'

   varnames(idx_for_star(VAR_HORB,  1)) = '   Horb'
   varnames(idx_for_star(VAR_ECC,   1)) = '      e'
   varnames(idx_for_star(VAR_XI,    1)) = '     xi'
   varnames(idx_for_star(VAR_BMASS, 1)) = '   Mbin'

   !-----------------------------------------------------------------------
   ! explicitly initialize local variables (Steve, 5/08).

   ikm(:) = 0
   errpr = 0
   vx = 0
   vz = 0

   frr = 0.0_dbl
   !-----------------------------------------------------------------------

   ! Set output unit for star 1 or star 2, it it is still at its default value
   if (solver_output_unit(1) == 1 .or. solver_output_unit(1) == 2) solver_output_unit(1) = jb

   ! Set number of equations/variables.
   ! Used to be a block copy, but does not work reliably when using
   ! aggressive optimisations in the compiler.
   ke1 = ig(1)
   ke2 = ig(2)
   ke3 = ig(3)                   ! Third order equations, always 0
   kbc = ig(4)
   kev = ig(5)
   kl  = ig(7)
   jh1 = ig(8)
   jh2 = ig(9)
   jh3 = ig(10)
   kd( 1:40)  = ig(11:50)              ! Permutations of variables/equations - variables
   kd(nkd+1:nkd+40)  = ig(51:90)       ! Permutations of variables/equations - equations
   kd(2*nkd+1:2*nkd+40) = ig(91:130)   ! Permutations of variables/equations - boundary conditions
   if (joc > 1) then ! Nucleosynthesis post-processing: no permutations
      do ij = 1, nkd
         kd(ij) = ij
         kd(nkd+ij) = ij
         kd(2*nkd+ij) = ij
      end do
   end if
   kq = 1 - 2*kl
   keq = ke1 + ke2               ! Total number of equations in block matrix
   if ( keq == 0 ) return        ! No equns to solve
   kvb = keq + kev               ! Total number of variables (incl. EVs)
   kvc = kvb + kvb               ! ?
   ke4 = ke3 + ke2               ! Total number of equations
   ! Locate vertical (kj) and horizontal (ki) partitions in matrix
   kj1 = kbc + 1
   kj2 = keq + 1
   kj3 = keq + kbc
   kj4 = kj3 + 1
   kj5 = keq + keq
   kj6 = kj5 + 1
   kj7 = kj5 + kbc
   kj8 = kj7 + 1
   kj9 = kj5 + keq
   kj10 = kj9 + 1
   kj11 = kj9 + kev
   kj12 = kj10 + kev
   ki1 = keq - kbc
   ki2 = ki1 + 1
   ki3 = ki1 + kev
   ki4 = ki3 + 1

   ! Construct reverse lookup table for the permutation of independent variables
   ikd = 0
   forall (i = 1:kvb) ikd(kd(i)) = i

   ! If the number of variables/meshpoints is larger than the size of the
   ! arrays (eg, we loaded a new model), we need to reallocate them.
   if ( allocated(S) .and. (kh>block_nmesh .or. kvb>block_vars) ) then
      deallocate( S )
      deallocate( C )
      deallocate( xd )
      deallocate( ddh )
      deallocate( er )
      deallocate( eqn_scale )
      deallocate( sbc_scale )
      deallocate( func )
      deallocate( dfunc )
      deallocate( diffx )
      !     deallocate( GRADF )
   end if

   ! Allocate memory for storage of the Jacobian
   if (.not. allocated(S)) then
      block_vars = kvb
      block_nmesh = kh
      allocate( S(3*kvb+1, kvb) )
      allocate( C(kh+1, kvb, kvb+1) )
      allocate( xd(kh) )
      allocate( ddh(kvb, kh) )
      allocate( er(NVAR) )
      allocate( eqn_scale(NVAR) )
      allocate( sbc_scale(NVAR) )
      allocate( func(NFUNC, kh+1) )
      allocate( dfunc(NFUNC, kvb, kh+1) )
      allocate( diffx(kvb, 0:kh+1) )
      !allocate( gradf(keq*kh + kev) )
   end if

   if (.not. allocated(default_er)) allocate(default_er(NVAR))
   default_er = 1.0_dbl

   ! First meshpoint is ik_first, last is ik_last
   ik_first = 1 + kl*(kh - 1)
   ik_last = kh + 1 - ik_first

   ! Determine 'typical', usually maximal, values for each variable
   if (joc == 1) then
      do ij = 1, kvb
         i = kd(ij)
         er(i) = default_er(i)
         do k = 1, kh
            er(i) = max(er(i), abs(H(i, k)))
         end do
         if ( er(i) == 0.0d0 ) er(i) = 1.0d0
         if ( jnn > 1 .and. es(i) > 0.0d0 .and. es(i) /= 1.0d0 ) er(i) = 0.5d0*(er(i) + es(i))
      end do
   else
      do ij = 1, kvb
         i = kd(ij)
         er(i) = 1.0d-8
         do k = 1, kh
            er(i) = max(er(i), abs(hnuc(joc-1,i, k)))
         end do
      end do
   end if

   ! Now determine typical values for each equation, based on current
   ! values in H(,) and the current timestep.
   ! Unfortunately we need to go through FUNCS1() to do this (perhaps this
   ! could better be collected in printb and then passed in here).
   call scale_eqn ( eqn_scale, sbc_scale )

   ! KEE is the left-most KEQxKVB sub-block we need to solve. if there are
   ! no higher order equations to solve, the left-most block (1, derivatives
   ! wrt variables on the previous meshpoint) is always filled with 0s, so
   ! we skip it in that case.
   kee = 1
   err = 1.0d0
   max_err = 1.0d0
   frrpr = 0.0d0
   if (ke4 == 0) kee = 2

   ! Begin iterative loop
   fac = 1.0d0
   besterr = 1.0d0
   converged = .false.
   C(:,:,:) = 0.0d0
   S(:, :) = 0.0d0
   residue = 0.0d0
   prev_residue = 0.0d0

   updated_once = .false.
   jo = 0
   do jter = 1, iter
      errpr = err
      err = 0.0d0
      jkh = 0
      !if ( JNN == JH1 .and. JTER == JH2 ) JKH = 1
      ! Calculate new values of quantities that are evaluated between iterations
      updated_explicit = .false.
      if (joc == 1) then
         if (jter == 1 .or. errpr < 8.*eps .or. updated_once) then
            if (comp_semi_implicit_quantities()) updated_explicit = .true.
            if (jter > 1 .and. updated_explicit) updated_once = .true.
         end if
      else
         updated_explicit = .true.
      end if

      ! Solve system of linear equations by Gaussian elimination
      call invert_matrix( jo )
      if (joc == 1 .and. jo == 1) exit

      ! Make sure changes (especially abundances) are "reasonable"
      converged_ok = .true.
      do ij = 1, kvb
         if (is_abundance(kd(ij))) then
            do ik = 1, kh
               if (abs(C(ik, ij, 1)) > 1.0d-2) then
                  !print *, ik, H(kd(ij), ik), C(ik, ij, 1)
                  C(ik, ij, 1) = sign(1.0d-2*H(kd(ij), ik), C(ik, ij, 1))
                  converged_ok = .false.
               end if
            end do
         end if
      end do

      ! Force to 0 changes that will be discarded anyway
      C(kh, ikd(idx_for_star(VAR_LNR, 1)), 1) = 0.0d0
      C(kh, ikd(idx_for_star(VAR_MASS, 1)), 1) = 0.0d0
      C(kh, ikd(idx_for_star(VAR_LUM, 1)), 1) = 0.0d0

      ! Map corrections to the range [-1,1].
      ! Corrections that are small are hardly changed at all (the mapping is
      ! linear for small C), while corrections that are (very) large are reduced
      ! significantly (to 1.0), but smoothly.
      ! It is an interesting idea to pass the result through another function
      ! that has the effect of reducing the limiting value from +/-1, but this
      ! seems to not improve convergence (or make it worse) so this is commented
      ! out for the time being. It's still possible that there is a way to make
      ! that work.
      C(1:kh, 1:kvb, 1) = tanh(C(1:kh, 1:kvb, 1))

      if (maxval(abs(C(1:kh, 1:kvb, 1))) > climit) then
         C = 0.5d0 * C
         converged_ok = .false.
      end if

      ! Estimate accuracy of iteration
      var_max_err = 1
      max_errpr = max_err
      max_err = 0.0d0
      do ij = 1, kvb
         vx = 0.0d0
         do ik = 1, kh
            ddh(ij, ik) = C(ik, ij, 1) * er(kd(ij))
            vz = abs(C(ik, ij, 1))
            if ( vz >= vx ) then
               vx = vz
               ikm(ij) = ik
            end if
            err = err + vz
         end do
         ert(ij) = vx
         if (vx > max_err) then
            point_max_err = ikm(ij)
            var_max_err = kd(ij)
            max_err = vx
         end if
      end do
      err = err/(kvb*kh)
      if (err /= 0.0d0 .and. jo == 1) exit
      if (err == 0.0d0) err = wanted_eps**2

      outres = residue
      if (residue == 0.0d0) outres = 1.0d-16

      ! Alternative to using a line search algorithm: decrease FAC, more or less
      ! ad-hoc
      if (jter > 1 .and. (err > errpr .or. max_err > max_errpr)) then
         fac = fac_down * fac
      else
         fac = min(del/max(err, max_err * eps/max_eps, del), fac_up*fac)
      end if
      if (err < (10.d0*eps)) fac = max(fac, 1.0d0 - err / (10.d0*eps))
      if (err < eps) fac = 1.0d0
      fac = max(fac, fac_min)
      if (joc > 1) then
         fac = 0.9
         if (err < 1.0d-5) fac = min(1.0d0, fac*1.05)
      end if
      if (use_linesearch) then
         fac = 1.0d0
         if (jter > 1) fac = linesearch( jo )
      end if
      frrpr = frr
      frr = log10(err)
      ert(1:kvb) = log10(ert(1:kvb) + 1.3d-10)

      if (abs(frr-frrpr) < 1.0d-2 .and. err > eps) fac = fac * fac_down
      if (joc == 1 .and. abs(frr-frrpr) < 1.0d-3 .and. fac < 1.0d-2) fac = 1.0

      ! write convergence information to summary file .out1/2
      if (joc == 1) then
         if ( jter == kt5+1 ) write (solver_output_unit(1), '(a4,2x,a3,1x,a4,2x,a3,1x, 17(1x,a7),/,3(20x,17(1x,a7),/))') &
             'iter', 'err', 'ferr', 'fac', (varnames(kd(ij)), ij=1,kvb)
         if ( jter > kt5 ) write (solver_output_unit(1), 992) jter, frr, log10(outres), fac, (ikm(ij), ert(ij), ij = 1, kvb)
         flush ( solver_output_unit(1) )
      else
         if ( jter == kt5+1 ) write (solver_output_unit(joc), "('*', I1)") joc-1
         if ( jter == kt5+1 ) write (solver_output_unit(joc), '(a4,2x,a3,1x,a4,2x,a3,1x, 17(1x,a7),/,3(20x,17(1x,a7),/))') &
             'iter', 'err', 'ferr', 'fac', (nucname(kd(ij)), ij=1,kvb)
         if ( jter > kt5 ) write (solver_output_unit(joc), 992) jter, frr, log10(outres), fac, (ikm(ij), ert(ij), ij = 1, kvb)
         flush ( solver_output_unit(joc) )
      end if
992  format (i3, 3f6.2, 17(i4, f4.1),/, 3(21x, 17(i4, f4.1),/))

      ! Apply corrections, scaled down if too large
      ! Perform sanity checks on abundances and boundary conditions
      if (joc == 1) then
         forall (i = 1:kvb) DH(kd(i), 1:kh) = DH(kd(i), 1:kh) - fac*ddh(i, 1:kh)
         call checks
      else
         forall (i = 1:kvb) DHnuc(joc-1, kd(i), 1:kh) = DHnuc(joc-1, kd(i), 1:kh) - fac*ddh(i, 1:kh)
         call checks2
      end if

      if ( err < besterr .and. updated_explicit .and. converged_ok) then
         if (joc == 1) then
            besth(1:NVAR,1:kh) = H(1:NVAR,1:kh)
            bestdh(1:NVAR,1:kh) = DH(1:NVAR,1:kh)
         else
            besth(1:50,1:kh) = Hnuc(joc-1, 1:50, 1:kh)
            bestdh(1:50,1:kh) = DHnuc(joc-1, 1:50, 1:kh)
         end if
         besterr = err
      end if

      ! For the nucleosynthesis, don't break the run if corrections are huge
      ! since this could be fine.
      if (joc == 1) then
         if (jnn > 10 .and. jter > 1) then
            if ( err > 1.0d1*del ) jo = 1
            !    Check all variables: if one of them has a very big error, then abort.
            !    In practice we will not converge in this situation anymore.
            !    FIXME: if the "best" solution we have found would be acceptable, then accept it
            do i = 1, kvb
               if (ert(i) > max_erk) jo = 1
            end do
         end if
      else
         if (err > 5.0d1) jo = 1
      end if

      ! Check whether the solution has converged sufficiently accurately that
      ! we can terminate the iteration
      if (err < wanted_eps .and. max_err < max_eps) converged = .true.
      if (residue > 0.0d0 .and. residue < accept_res) converged = .true.
      converged = converged .and. updated_explicit .and. converged_ok

      ! Nucleosynthesis should be a little bit more accurate
      if (joc > 1) then
         converged = .false.
         if (err <= 0.01*wanted_eps) converged = .true.
      end if

      ! Break out if not making progress
      if (joc > 1 .and. abs(frr-frrpr) < 1.0d-3) exit

      ! RJS sort out neutrons once convergence is complete
      ! (could/should be done earlier?)
      if ( converged .and. joc > 1 ) then
         !call checks3
         call neutron(joc-1)
      end if

      if ( jo == 1 ) exit                 ! Convergence failure
      if ( converged ) then
         if (joc == 1) es(:) = er(:)      ! Store typical size of variables
         return
      end if

      ! Continue iterating if not yet accurate enough
   end do


   ! Aledgedly not converged (to accuracy WANTED_EPS), but we may have
   ! converged to accuracy EPS!
   if ( (joc == 1 .and. besterr < eps) .or. (joc > 1 .and. besterr < 1.0d-1*eps) ) then
      frr = log10(besterr)
      if (kt5 < iter) write (solver_output_unit(1), "('Accepting solution with err = ', f6.2)") frr

      if (joc == 1) then
         H(1:NVAR,1:kh) = besth(1:NVAR,1:kh)
         DH(1:NVAR,1:kh) = bestdh(1:NVAR,1:kh)
         es(:) = er(:)
      else
         Hnuc(joc-1, 1:50,1:kh) = besth(1:50,1:kh)
         DHnuc(joc-1, 1:50,1:kh) = bestdh(1:50,1:kh)
         call neutron(joc-1)
         !call checks3
      end if
      jo = 0
      return
   end if
   jo = 1

end subroutine solver




subroutine invert_matrix(jo)
   use real_kind
   use mesh
   use solver_global_variables

   implicit none
   integer, intent(inout) :: jo

   ! Local variables
   integer :: jk, kk, jklast, ii, ij

   prev_residue = residue
   residue = 0.0d0
   residue_deriv = 0.0d0
   !  GRADF(:) = 0.0D0
   xd(:) = 0.0d0
   do jk = ik_first + kl, ik_last+kl+kq, kq
      call full_differences(jk)
   end do

   ! Evaluate functions at first meshpoint:
   jk = ik_first + kl
   call differences ( jk, ki2, keq, 2*nkd + kev )
   call divide_rows ( ki2, kj6, keq, jk, jo )
   if ( jo == 1 ) write(solver_output_unit(1), *) 'Singular matrix', jk
   if ( jo == 1 ) return
   jk = jk + kq
   if ( ik_first == ik_last ) return

   ! Ditto second, eliminating some unknowns:
   call differences ( jk, 1, keq, nkd )
   call eliminate   ( 1, kj2, keq, kj3, ki2, kj4, kj5, jk - kq)   !,   1 )
   call divide_rows ( 1, kj4, keq, jk, jo)
   if ( jo == 1 ) write(solver_output_unit(1), *) 'Singular matrix', jk
   if ( jo == 1 ) return

   ! Ditto remaining meshpoints except last:
   do jk = ik_first + kl + 2*kq, ik_last + kl, kq
      call differences ( jk, 1, keq, nkd )
      call eliminate   ( 1,   1, ke4, kbc, ki2, kj1, keq, jk - 2*kq)  !  , 1 )
      call eliminate   ( 1, kj1, ke4, keq,   1, kj4, kj5, jk - kq)    !,   2 )
      call eliminate   ( 1, kj2, keq, kj3, ki2, kj4, kj5, jk - kq)    !,   3 )
      call divide_rows ( 1, kj4, keq, jk, jo)
      if ( jo == 1 ) write(solver_output_unit(1), *) 'Singular matrix', jk
      if ( jo == 1 ) return
   end do
   jk = ik_last + kl + kq

   ! Ditto last meshpoint:
   call differences ( jk, 1, ki3, 2*nkd )
   call eliminate   ( 1, kj2, ke4, kj3, ki2, kj4, kj5, jk - 2*kq)   !, 1 )
   call eliminate   ( 1, kj4, ke4, kj5,   1, kj8, kj9, jk - kq)     !,   2 )
   call eliminate   ( 1, kj6, ki3, kj7, ki2, kj8, kj9, jk - kq)     !,   3 )

   ! Solve for corrections at last meshpoint
   call divide_rows ( 1, kj8, ki3, jk, jo)
   if ( jo == 1 ) write(solver_output_unit(1), *) 'Singular matrix', jk
   if ( jo == 1 ) return

   ! By back-substitution find corrections throughout
   ii = 1
   jklast = ik_first + kl
   if (kbc == 0) jklast = jklast-kq
   do jk = ik_last + kl, jklast, -kq
      if ( jk == ik_first + kl ) ii = ki2
      kk = jk + kq
      do ij = 1, ki2-1
         C(jk, ii:keq, ki4) = C(jk, ii:keq, ki4) - C(jk, ii:keq, ij)*C(kk, ij, ki4)
      end do
      kk = ik_last + 1 - kl
      do ij = ki2, ki3
         C(jk, ii:keq, ki4) = C(jk, ii:keq, ki4) - C(jk, ii:keq, ij)*C(kk, ij, ki4)
      end do
      !call PRINTC ( JK, JKH, II, KEQ, 1, KI4 )
   end do
   if ( kev /= 0 ) then
      forall (jk = 1:kh) C(jk, kj2:kvb, 1) = C(ik_last + 1 - kl, kj2 - kbc:kvb - kbc, ki4)
   end if
   if ( kbc /= 0 ) then
      C(1:kh, 1:kbc, 1) = C(kl+1:kh+kl, 1 + ki1:kbc + ki1, ki4)
   end if
   C(1:kh, 1 + kbc:ki1 + kbc, 1) = C(2 - kl:kh + 1 - kl, 1:ki1, ki4)
   !call PRINTC ( JK, JKH, 1, KVB, 1, 1 )
end subroutine invert_matrix



!> \brief Line search algoritm
!!
!! Line search algoritm, as described in Numerical Recipes, 2nd edition, Sect.9.7
!! Slightly different from the implemantation there.
!! This function returns the value "lambda" (or "FAC") by which the current
!! set of corrections are to be multiplied.
!!
!! \todo FIXME: in order for this to work properly, equations need to be scaled
!! with their typical values, otherwise checking whether the residuals are
!! "large" doesn't make sense. Needs to be fixed for TWIN mode.
function linesearch(jo)
   use real_kind
   use mesh
   use solver_global_variables

   implicit none
   integer, intent(inout) :: jo
   real(double) :: linesearch

   ! Local variables
   real(double), parameter :: alpha = 1.0d-6
   real(double) :: a, b, disc, rhs1, rhs2, slope
   real(double) :: lambda, lambda_min, new_lambda, prev_lambda
   real(double) :: init_residue
   real(double) :: dh_old(kvb, kh)
   integer :: ik, ij, i

   ! Backup DH
   ij = 0  ! Always initialised
   forall (ij = 1:kvb) dh_old(ij, 1:kh) = DH(kd(ij), 1:kh)

   ! Determine downward slope, grad f.dx
   ! slope = 0.0d0
   ! Contribution from normal variables
   !  do IK=1, KH
   !     do IJ = 1, KEQ
   !        SLOPE = SLOPE + GRADF((IK-1)*KEQ+IJ) * DDH(IJ, IK)
   !     end do
   !  end do
   ! Contribution from eigenvalues
   !  do IJ = KEQ+1, KVB
   !     SLOPE = SLOPE + GRADF((KH-1)*KEQ+IJ) * DDH(IJ, KH)
   !  end do

   ! Find allowed range of values for "lambda"
   lambda_min = 0.0d0
   do ik = 1, kh
      do ij = 1, kvb
         i = kd(ij)
         lambda_min = max( lambda_min, abs(ddh(ij, ik))/er(i) )
      end do
   end do
   lambda_min = 1.0d-4 / lambda_min
   lambda = 1.0d0
   linesearch = 1.0d0
   slope = -2*residue
   prev_lambda = lambda

   ! Store residual at original point
   init_residue = residue

   do
      !    New value of DH
      forall (ij = 1:kvb) DH(kd(ij), 1:kh) = dh_old(ij, 1:kh) + lambda*ddh(ij, 1:kh)

      !    Evaluate functions for new parameters
      !print *, 'Line search:', LAMBDA, RESIDUE
      call invert_matrix( jo )
      if ( jo == 1 ) exit

      !    Check for convergence on lambda (parameters)
      if (lambda < lambda_min) exit

      !    Check for convergence on function value
      if (residue <= init_residue + alpha*lambda*slope) exit

      !    Determine new value of lambda
      if (lambda == 1.0d0) then     !> \todo FIXME: float comparison
         new_lambda = -slope/(2.0d0*(residue - init_residue - slope))
      else
         rhs1 = residue - init_residue - lambda*slope
         rhs2 = prev_residue - init_residue - prev_lambda*slope
         a = (rhs1/lambda**2 - rhs2/prev_lambda**2) / (lambda - prev_lambda)
         b = -(prev_lambda*rhs1/lambda**2 - lambda*rhs2/prev_lambda**2) /&
             (lambda - prev_lambda)
         if (a == 0.0d0) then    ! This will NEVER happen in practice...
            new_lambda = -slope/(2.0d0*b)
         else
            disc = b**2 - 3.0d0*a*slope
            if (disc < 0.0d0) then
               new_lambda = 0.5_dbl*lambda
            else if (b<= 0.0d0) then
               new_lambda = (-b + sqrt(disc))/(3.0d0*a)
            else
               new_lambda = -slope/(b+sqrt(disc))
            end if
            new_lambda = min(0.5d0*lambda, new_lambda)
         end if
      end if
      prev_lambda = lambda
      lambda = max(new_lambda, 0.1d0 * lambda)
   end do

   ! Restore DH and return
   forall (ij = 1:kvb) DH(kd(ij), 1:kh) = dh_old(ij, 1:kh)
   linesearch = lambda
end function linesearch



subroutine full_differences(k)
   use real_kind
   use mesh
   use solver_global_variables
   use structure_functions
   use current_model_properties
   use settings, only: dh0
   use indices

   implicit none
   ! parameters:
   integer, intent(in) :: k
   real(double) :: var(NVAR), dvar(NVAR)
   integer :: i
   integer :: ji
   real(double) :: dx, dvx, vx

   if ( (kl == 0.and.k <= block_nmesh) .or. (kl == 1.and.k >= 2) ) then
      ! Evaluate argument list of variables and increments at current meshpoint.
      ! Next evaluate and store the required funcs of the independent variables
      ! varying independent variables in turn to evaluate numeric derivatives of
      ! funcs (if using numerical derivatives, joc == 1) or store results
      ! of analytical derivatives (joc > 1, nucleosynthesis code).
      ! These values can be recomputed (typically) or stored in a cache
      ! (not implemented).
      if ( joc == 1 ) then
         dvar(1:NVAR) = DH(1:NVAR, k - kl)
         var(1:NVAR) =  H(1:NVAR, k - kl) + dvar(1:NVAR)
         ! Evaluate and store the required functions of the independent variables
         call funcs1 ( k - kl, 0, var(:), dvar(:), func(:, k - kl) )

         ! Varying independent variables in turn, evaluate numeric derivatives of functions
         do i = 1, kvb
            ji = kd(i)
            vx = var(ji)
            dvx = dvar(ji)
            !dx = sign(dh0*max(abs(vx), 1.0d0), dvx)
            dx = -dh0*max(abs(vx), 1.0d0)
            if (is_abundance(ji)) dx = -dx
            var(ji) = vx + dx
            dvar(ji) = dvx + dx
            call funcs1 ( k - kl, ji, var(:), dvar(:), dfunc(:,i, k - kl) )

            diffx(i, k - kl) = dx
            dvar(ji) = dvx
            var(ji) = vx
         end do
      end if  ! joc /= 1
   end if

end subroutine full_differences



!> \brief Compute numerical or analytic derivatives
! JZ is the `start' index in the list of permutations
subroutine differences ( k, jx, jy, jz )
   use real_kind
   use mesh
   use solver_global_variables
   use nucleosynthesis
   use structure_functions
   use current_model_properties
   use settings, only: dh0
   use indices

   implicit none
   ! parameters:
   integer, intent(in) :: k, jx, jy, jz

   ! Local variables:
   integer :: i, iee, ieq, ii, ivb, j, jeq, jj, jvb3, jvb4
   integer :: ji, Jstar
   real(double) :: dx, dvx, vx
   real(double) :: ds

   ! common variables:
   real(double) :: var(NVAR), dvar(NVAR)
   real(double) :: equ(NEQ), dequ(NVAR,3,NEQ)
   real(double) :: varnuc(50), dvarnuc(50), fn1nuc(NFUNC), dfn1nuc(50, NFUNC)
   real(double) :: eqs(NEQ)

   !> \todo FIXME: these need to be persistent across function calls, but using
   !! the save attribute is probably not the best way to do this...
   real(double), save :: fn1(NFUNC), dfn1(NVAR,NFUNC)
   real(double), save :: fn2(3,NFUNC), dfn2(3,NVAR,NFUNC)


   !-----------------------------------------------------------------------
   ! explicitly initialize all local variables (Steve, 5/08).

   i = 0
   iee = 0
   ieq = 0
   ii = 0
   ivb = 0
   j = 0
   jeq = 0
   jj = 0
   jvb3 = 0
   jvb4 = 0
   ds = 0
   Jstar = joc - 1
   !-----------------------------------------------------------------------

   if ( jkh >= 1 .and. ( k < 4.or.k > block_nmesh-2 )) then
      jkh = jh3 + 1
   else
      jkh = 1
   end if
   if ( k == ik_first + kl) then
      fn2 = 0.0d0
      dfn2 = 0.0d0
   end if
   ! Redefine previous current meshpoint as current previous meshpoint
   if ( (kl == 0.and.k <= block_nmesh) .or. (kl == 1.and.k >= 2) ) then
      fn2(kee:2, 1:NFUNC)        =  fn2(kee+1:3, 1:NFUNC)
      dfn2(kee:2, 1:kvb, 1:NFUNC) = dfn2(kee+1:3, 1:kvb, 1:NFUNC)

      ! Evaluate argument list of variables and increments at current meshpoint.
      ! Next evaluate and store the required funcs of the independent variables
      ! varying independent variables in turn to evaluate numeric derivatives of
      ! funcs (if using numerical derivatives, joc == 1) or store results
      ! of analytical derivatives (joc > 1, nucleosynthesis code).
      ! These values can be recomputed (typically) or stored in a cache
      ! (not implemented).
      dvar(1:NVAR) = DH(1:NVAR, k - kl)
      var(1:NVAR) =  H(1:NVAR, k - kl) + dvar(1:NVAR)
      if ( joc == 1 ) then
         xd(1:kvc) = xd(1 + kvb:kvc + kvb)

         ! Evaluate and store the required functions of the independent variables
         fn1(:) = func(:, k-kl)

         ! Varying independent variables in turn, evaluate numeric derivatives of functions
         do i = 1, kvb
            ji = kd(i)
            vx = var(ji)
            dvx = dvar(ji)
            dx = diffx(i, k-kl)
            xd(i + kvc) = er(ji)/dx
            dfn1(i,:) = dfunc(:,i,k-kl)
         end do
      else  ! joc /= 1
         ! Alternative to be used if derivs computed analytically
         dvarnuc(1:50) = dhnuc(joc-1, 1:50,k-kl)
         varnuc(1:50) = hnuc(joc-1, 1:50,k-kl) + dvarnuc(1:50)
         fn1nuc = 0.0_dbl
         dfn1nuc = 0.0_dbl
         call funcs2 ( joc-1, k - kl, 1, var, varnuc, dvarnuc, fn1nuc, dfn1nuc )
      end if  ! joc /= 1
   end if

   if ( joc == 1 ) then  ! Numerical derivatives
      fn2(3, 1:NFUNC) = fn1(1:NFUNC)
      dfn2(3, 1:kvb, 1:NFUNC) = dfn1(1:kvb, 1:NFUNC)

      ! Evaluate and store the difference equns which are to be satisfied
      call equns1 ( k, kl, kq, fn2(:,:), equ(:) )
      forall (ieq = jx:jy)
         eqs(ieq) = eqn_scale(kd(ieq + jz))
         S(kj12, ieq) = equ(kd(ieq + jz))
      end forall
      if (k == 1) eqs(jx:jy) = 1.0d0
      residue = residue + 0.5d0*sum( (S(kj12, jx:jy)/eqs(jx:jy))**2 )

      ! Varying independent variables in turn, evaluate numeric derivatives of equns
      jvb3 = kee*kvb - kvb                   ! First column in this block
      S(kj5+keq+1:kj5+kvb, jx:jy) = 0.0d0    ! Clear derivatives to EVs
      do iee = kee, 3
         ! Loop over all normal variables
         fn1(1:NFUNC) = fn2(iee, 1:NFUNC)
         do ivb = 1, keq
            jvb4 = keq*(iee - 1) + ivb
            fn2(iee, 1:NFUNC) = dfn2(iee, ivb, 1:NFUNC)
            call equns1 ( k, kl, kq, fn2(:,:), equ(:) )
            forall (ieq = jx:jy) S(jvb4, ieq) = (equ(kd(ieq + jz)) - S(kj12, ieq))*xd(jvb3 + ivb)

            ! Force to 0 the derivative of the spin surface boundary condition with respect to the moment of inertia.
            ! We want the moment of inertia to be calculated directly, with the rotation rate derived from the boundary
            ! condition enforcing conservation of angular momentum. Keeping the moment of inertia in the boundary condition
            ! has a too strong effect on the moment of inertia itself and turns out to cause numerical instability
            if (kd(ivb) == VAR_INERT) then
               do ieq = jx, jy
                  if (kd(ieq + jz) == SBC_SSPIN) S(jvb4, ieq) = 0.0d0
               end do
            end if

            ! Store EQU.J vector
            ! I = BLOCK_NMESH*KEQ - (K + (2 - IEE))*KEQ + ivb
            ! GRADF(I) = GRADF(I) + S(KJ12,IEQ)*DS
            residue_deriv = residue_deriv - sum(S(kj12, jx:jy)**2)
         end do
         jvb3 = jvb3 + keq

         ! Loop over all eigenvalues; these derivatives have their own
         ! column off to the right of the matrix.
         do ivb = keq+1, kvb
            jvb3 = jvb3 + 1
            jvb4 = kj5 + ivb                 ! Eigenvalues in column
            fn2(iee, 1:NFUNC) = dfn2(iee, ivb, 1:NFUNC)
            call equns1 ( k, kl, kq, fn2(:,:), equ(:) )
            do ieq = jx, jy
               jeq = kd(ieq + jz)
               S(jvb4, ieq) = S(jvb4, ieq) + (equ(jeq) - S(kj12, ieq))*xd(jvb3)
               residue_deriv = residue_deriv - S(kj12, ieq)**2/block_nmesh
               ! Store EQU.J vector (eigenvalues go at the end)
               !              I = BLOCK_NMESH*KEQ + ivb - KEQ
               !              GRADF(I) = GRADF(I) + S(KJ12, IEQ)*DS
            end do
         end do
         fn2(iee, 1:NFUNC) = fn1(1:NFUNC)
      end do

   else                      ! joc /= 0: analytic derivatives

      dfn2(3, 1:kvb, 1:NFUNC) = dfn1nuc(1:kvb, 1:NFUNC)
      fn2(3, 1:NFUNC) = fn1nuc(1:NFUNC)
      ! Alternative to be used if derivs computed analytically
      call equns2 ( k, kl, kq, keq, Jstar, fn2, dfn2, equ, dequ )
      do j = jx, jy
         jj = kd(j + jz)
         forall (i = 1:keq)
            S(i, j) = dequ(i, 1, jj) * er(kd(i))
            S(i + keq, j) = dequ(i, 2, jj) * er(kd(i))
            S(i + kj5, j) = dequ(i, 3, jj) * er(kd(i))
         end forall
         if ( kev /= 0 ) then
            forall (i = kj2:kvb)
               S(i + kj5, j) = (dequ(i, 1, jj) + dequ(i, 2, jj) + dequ(i, 3, jj)) * er(kd(i))
            end forall
         end if
         S(kj12, j) = equ(jj)
      end do
   end if
end subroutine differences



subroutine divide_rows ( id1, jd1, id2, k, jo )
   use real_kind
   use mesh
   use solver_global_variables
   use current_model_properties

   implicit none
   ! Row flags
   integer, parameter :: forced_pivot = 1
   integer, parameter :: current_pivot = 2
   integer, parameter :: later_pivot = 3
   !integer, parameter :: used_pivot = 4  ! Unused
   integer, parameter :: possible_pivot = 5
   integer, parameter :: singular = 6

   ! Value to promote a "future" pivot to a current candidate pivot
   integer, parameter :: next_pivot = 2

   ! parameters:
   integer, intent(in) :: id1, jd1, id2, k
   integer, intent(out) :: jo

   ! Local variables:
   integer :: i, ii, im, j, jc, jd2, jd3, jj, jl, jmi
   real(double) :: vm, vs, vt, vx
   logical :: all_done

   integer :: itype(NVAR), jmax(NVAR+1)
   logical :: idone(NVAR)

   !-----------------------------------------------------------------------
   ! explicitly initialize all local variables (Steve, 5/08).

   i = 0
   ii = 0
   im = 0
   j = 0
   jc = 0
   jd2 = 0
   jd3 = 0
   jj = 0
   jl = 0
   jmi = 0
   vm = 0
   vs = 0
   vt = 0
   vx = 0
   itype = 0

   !-----------------------------------------------------------------------

   jc = 1
   jd2 = jd1 + id2 - id1
   if ( id1 > id2 .or. jd2 + 1 > kj12 ) return
   !if ( JKH >= 2 ) call PRINTS ( K, JKH, ID1, ID2, 0 )
   idone(id1:id2) = .false.
   itype(id1:id2) = possible_pivot
   all_done = .false.

   do while (jc < 5000 .and. .not. all_done)
      jc = jc + 1
      do i = id1, id2
         if ( itype(i) < later_pivot ) itype(i) = itype(i) + next_pivot
         if ( itype(i) < possible_pivot ) idone(jmax(i)) = .true.
      end do
      vt = 0.0

      ! Locate the most significant (remaining) row. `Significance' of a row is
      ! the ratio of largest |element| = VM to sum of remaining |elements| = VS
      do i = id1, id2
         if ( itype(i) >= possible_pivot ) then
            vm = 0.0
            jl = 0
            do j = jd1, jd2
               jj = j - jd1 + id1
               if ( .not. idone(jj) ) then
                  vx = abs(S(j, i))
                  if ( vx >= vm ) then
                     vm = vx
                     jl = jj
                  end if
               end if
            end do

            if ( jl < 1 .or. jl > block_vars+1 ) then
               jo = 1
               return
            end if

            jmax(i) = jl
            vs = 0.0d0
            do j = jd1, jd2
               if ( j - jd1 + id1 /= jl ) vs = vs + abs(S(j, i))
            end do
            if ( vm ==  0.0d0 ) then
               itype(i) = singular
               vx = 0.0d0
            end if
            if ( vs == 0.0d0 ) then    ! Only one non-zero element in row
               if ( vm > 0.0d0 ) then
                  itype(i) = forced_pivot
                  idone(jl) = .true.
                  vx = 2.0d0
               end if
            else
               vx = vm/(vm + vs)
            end if
            if ( vx >= vt ) then
               vt = vx
               im = i
            end if
         end if
      end do

      if ( im < 1 .or. im > NVAR )  then
         jo = 1
         return
      end if

      if ( itype(im) == possible_pivot ) itype(im) = current_pivot

      ! Largest element moduluswise of most significant row is the leading
      ! pivot; eliminate elements above and below it
      do i = id1, id2
         im = itype(i)
         if ( im < later_pivot ) then
            jmi = jmax(i) + jd1 - id1
            jd3 = jd1
            if ( im == forced_pivot ) jd3 = jd2 + 1
            vx = 1.0d0/S(jmi, i)
            S(jd3:kj12, i) = vx*S(jd3:kj12, i)
            S(jmi, i) = 1.0d0
            do ii = id1, id2
               if ( itype(ii) > later_pivot ) then
                  S(jd3:kj12, ii) = S(jd3:kj12, ii) - S(jmi, ii)*S(jd3:kj12, i)
                  S(jmi, ii) = 0.0d0
               end if
            end do
         end if
         idone(i) = .false.
      end do

      ! Are we done now?
      all_done = .true.
      do i = id1, id2
         if ( itype(i) == possible_pivot .or. itype(i) <= current_pivot ) then
            all_done = .false.
            exit
         end if
      end do

   end do

   forall(i=id1:id2, j=jd2+1:kj12)  C(k, jmax(i), j - kj12 + ki4) = S(j, i)

   if ( jkh >= 2 ) return
   if ( jc < 5000 ) return

   jo = 1

   ! Some emergency printout:
   !write (10, '(7I5, /, 40I3, /, 41I3, /, 40I3)') K, JNN, JTER, JKH, JL, IM, JC, ITYPE, jmax, idone
   !call PRINTC ( K, JKH, ID1, ID2, JD2 + 1 - KJ12 + KI4, KI4 )

end subroutine divide_rows



subroutine eliminate ( is1, js1, is2, js2, is3, js4, js5, k )
   use real_kind
   use mesh
   use solver_global_variables
   use current_model_properties

   implicit none
   ! Parameters:
   integer :: is1, js1, is2, js2, is3, js4, js5, k

   ! Local variables:
   integer :: i, jj, jjj
   real(double) :: vx

   !-----------------------------------------------------------------------
   ! explicitly initialize all local variables (Steve, 5/08).

   i = 0
   jj = 0
   jjj = 0
   vx = 0
   !-----------------------------------------------------------------------

   !if ( JKH >= 2 ) call PRINTS ( K, JKH, IS1, IS2, JCL )
   if ( js1 > js2 .or. 1 > is2 ) return
   if ( js4 <= js5 ) then
      do jj = js1, js2
         jjj = jj - js1 + is3
         forall (i=is1:is2)  S(js4:js5, i) = S(js4:js5, i) - S(jj, i)*C(k, jjj, 1:js5 - js4 + 1)
      end do
   end if
   do jj = js1, js2
      jjj = jj - js1 + is3
      forall (i=is1:is2)  S(kj10:kj12, i) = S(kj10:kj12, i) - S(jj, i)*C(k, jjj, kj10 - kj12 + ki4:ki4)
   end do
end subroutine eliminate




subroutine prints ( jk, jkh, ips1, ips2, jcl )
   use real_kind
   use mesh
   use solver_global_variables
   implicit none

   !Common block:
   integer :: jk, jkh, ips1, ips2, jcl
   integer :: i,j

   write (10,*) 'printS', jcl, jk, jkh, kee, keq, kj2, kj5, kj6, kj12
   if ( kee /= 2 ) then
      do i = ips1, ips2
         write (10, 991) jk, (S(j, i), j = 1, keq)
      end do
      write (10, 991) jk
   end if
   do i = ips1, ips2
      write (10, 991) jk, (S(j, i), j = kj2, kj5)
   end do
   write (10, 991) jk
   do i = ips1, ips2
      write (10, 991) jk, (S(j, i), j = kj6, kj12)
   end do
   write (10, 991) jk
991 format (i5, 1p, 17d9.1, /, 17d9.1)

end subroutine prints



subroutine printc ( jk, jkh, ipc1, ipc2, jpc1, jpc2 )
   use real_kind
   use mesh
   use solver_global_variables
   implicit none

   integer :: jk, jkh, ipc1, ipc2, jpc1, jpc2
   integer :: i,j

   if ( jkh < 2 .or. (jk >= 4 .and. jk <= kh - 2) ) return

   write (10,*) 'printC', jk, jkh, jpc1, jpc2
   do i = ipc1, ipc2
      write (10, '(1p, 15d10.2)') (C(jk, i, j), j = jpc1, jpc2)
   end do
end subroutine printc


