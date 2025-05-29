!!> \file stars_structure_functions.f90   Contains funcs1()
#include "assert.h"

module structure_functions

contains

!> ------------------------------------------------------------------------------
!!  FUNCS1:
!!   Evaluate the functions of independent variables needed to set up the
!!   difference equations.
!! ------------------------------------------------------------------------------
!!  INPUT:
!!     JK - the number of the meshpoint to calculate values for
!!     JI - determines from where the function is called, or which
!!          variable was changed:
!!           >0: The number of the variable that has been varied for this call
!!            0: No variable changed
!!           -1: PRINTB
!!           -2: REMESH
!!           The last two options affect what is calculated by the code,
!!           the others are used for caching the output from the
!!           equation-of-state.
!!     VAR(:)  - Array with values for independent variables
!!     DVAR(:) - Array with changes in independent variables
!!  OUTPUT:
!!     FN1(:) - Ouput quantities ("functions")
!!     EOSOUT - (optional) struct to pass the output from statef() to the caller
!!     PX(:) - Stellar structure is stored in the PX array (JI<=0 only)
!!
!! \todo FIXME: because of the presence of an optional argument, this function
!! needs an "explicit interface" (in C terms, a declaration) so that the
!! argument list is known when the function is called and the compiler can
!! do the right thing when the optional argument is not present (the code
!! segfaults otherwise).
!! Ideally we would place funcs1 in a module, which means the compiler will
!! generate the interface for us, but there is a circular dependance of
!! funcs1 on modules defined in output_properties.f90. That needs to be fixed
!! first. For now, we write the interface ourselves (in structure_functions.f90)
!! ------------------------------------------------------------------------------

subroutine funcs1 ( jk, ji, var, dvar, fn1, eosout, abundout, px )
   use real_kind
   use cuberoot
   use mesh
   use mesh_enc
   use control
   use constants
   use settings
   use massloss
   use plotvariables
   use explicit_functions
   use nucleosynthesis
   use semi_implicit_variables
   use atomic_data
   use test_variables
   use eostate_types
   use current_model_properties
   use accretion_abundances
   use radiative_acceleration
   use step_functions
   use roche
   use indices
   use stopping_conditions

   implicit none
   ! Declaring subroutine variables
   integer, intent(in) :: jk, ji
   real(double), intent(in) :: var(NVAR), dvar(NVAR)
   real(double), intent(out) :: fn1(NFUNC)
   type(eostate), optional, intent(out) :: eosout
   type(abundance), optional, intent(out) :: abundout
   real(double), optional, intent(out) :: px(npx)

   ! Declaring common blcok variables explicitly

   real(double) :: af, at, vx16, m, vx1, vqk, ar, L, vx4, vx12, vx20, ai, vomega, phi, phis, vx14, vx24, vx28, vx56   ! stellar vars
   real(double) :: oa, ecc, xi, mb, pmass                                                                             !  binary vars
   real(double) :: daf, dat, dx16, dm, dx1, dvqk, dar, dL, dx4, dx12, dx20, dai, domega, dphi, dphis          ! stellar var. changes
   real(double) :: dx14, dx24, dx28, dx56                                                                     ! stellar var. changes
   real(double) :: doa, decc, dxi, dmb                                                                        !  binary var. changes
   real(double) :: tam, dtam

   real(double) :: vp, vpk, vr, vrk, vt, vtk, vl, lk, lq, mt, vm, vmk, sg, sgth
   real(double) :: x1, x1t, x4, x4t, x12, x12t, x14, x14t, x16, x16t, x20, x20t, x24, x24t, x28, x28t, x56, x56t
   real(double) :: avmu, avmu_new
   real(double) :: vik, omega, omegat
   real(double) :: bcp, bct, bcm, phik, bcf, bcs, bcph
   real(double) :: bca, bce, xim, xik, dlrk, bcmb, bcmp

   real(double) :: xfs(nfunc)                           ! Extra output variables for star (S)

   real(double) :: xh, xhe, xc, xn, xo, xne, xmg, xsi, xfe, na(9), neo, nio, nzz, avma, ne, xai(9, 26)

   ! Various things, needed for gravitational settling:
   integer :: nref                        ! Reference abundance (H)
   integer :: i
   real(double) :: gi, apr, atr
   real(double) :: wpi(9), ddiff(9), xa(9)
   real(double) :: wpp, kf_num, kf_denom
   real(double) :: avgz(9), nn(9), nne, diffp(9),difft(9),diffg(9),diffc(9,9)
   real(double) :: logkappa, Frad

   real(double) :: ap, arho, u, uint,urec,uass, p, rho, fk, t, sf, st, zt, grada, scp, rf, rt, xhi, s, pr, pg, pf, pt, en
   real(double) :: rpp, r33, r34,  rbe, rbp, rpc, rpna,  rpo,  r3a,  rac
   real(double) :: ran, rao, rane, rcca, rco, roo, rgne, rgmg, ramg, rgsi, rccg, rpng
   real(double) :: ex, enx, wmu, delta, phie, ext, fkt, fkr, prandtl
   type(eostate) :: eos
   type(abundance) :: abund

   ! Declaring local variables
   integer :: Jstar, j, Jstar_accretor

   real(double) :: r, r2, r3, wth, target_daf, target_dat, ds, mor3, r0, lor3,  r2mu, fp, ft
   real(double) :: gor, apmu, gradr, s2, wc1, mu, lom, wc2, wc3, wc4, gta, atmu, vpp, vtt, xmc, vmf, vmmu, vrr, muk, mk, b
   real(double) :: wos, con, fac, vml, df0, sgf, deg, ug, fac1, fac2, v2, sgth_numer, sgth_denom, dth, sgth2
   real(double) :: app, mmk, absl, avm, xh_surface, xhe_surface, z_eff, rht, rad, w2, tet, rot, mdtddw, alfac, mdtl
   real(double) :: mdtsw, mdot_jvink, mdot_wr, mdot_cmi, mdot, s_x, s_he, s_c, s_o, f_photon_tiring, fdt, e2, ef, we25, rl
   real(double) :: rgr, etgr, oatgr, rlfac, xi1, ww, gam, we5, we65, we2, we3, wea, web, dmt0, cif, alpha, tau, mm

   real(double) :: rpole         ! Polar radius
   real(double) :: roche_phi(2)  ! Scaled Roche potential
   real(double) :: phi_l1        ! Potential at L1
   real(double) :: delta_phi(2)  ! Scaled Roche potential difference wrt L1
   real(double) :: tm(2)         ! Total mass for each component
   real(double) :: m1, m2        ! Masses for calculation of potential, such that m1>m2
   real(double) :: xx, yy        ! Coordinates on the surface of the star for which to calculate the Roche potential

   real(double) :: sam           ! Specific Angular Momentum (j)
   real(double) :: dsam          ! Variation of Specific angular momentum (dj)
   real(double) :: samt          ! Time variation of Specific angular momentum (dj/dt)
   real(double) :: sgam          ! Diffusion coefficient for angular momentum
   real(double) :: si            ! Specific moment of inertia
   real(double) :: dsi           ! Variation of specific moment of inertia
   real(double) :: meff          ! "Effective" mass; reduced for rotation
   real(double) :: Vmc2          ! P2(cos theta) coeff. of horizontal component of meridional circulation
   real(double) :: Umc2          ! P2(cos theta) coeff. of vertical component of meridional circulation
   real(double) :: Umc4          ! P4(cos theta) coeff. of vertical component of meridional circulation
   real(double) :: adam          ! Advection flux of angular momentum

   real(double) :: grad          ! grad T (used to be called gradt in some routines, change back?)
   real(double) :: dg            ! Convective grad_r - grad_a  for Schwarzschild criterion
   real(double) :: wcv           ! Convective velocity
   real(double) :: wl            ! Convective velocity x mixing length
   real(double) :: hp            ! Pressure scale height
   real(double) :: htemp         ! Temperature scale height
   real(double) :: NT2, Nmu2     ! Brunt-Vaisala  frequency
   real(double) :: Dsh, Dsv      ! Horizontal and vertical components of the shear
   real(double) :: Deff          ! Effective diffusion coefficient for chemical mixing due to shear
   real(double) :: zeta          ! Shear rate, dlog omega / dlog r
   real(double) :: alpha_rot     ! Gradient of specific angular momentum, 0.5 dlog r**2 omega / dlog r = 1 + 0.5 zeta

   real(double) :: qq            ! Determines mesh-point metric: mesh-point interval
   real(double) :: qm            ! Determines mesh-point metric: derivative of mesh spacing function wrt mass
   real(double) :: phim          ! Derivative of gravitational potential with respect to m**(2/3)
   real(double) :: gmr           ! Effective gravity at the surface(?)
   real(double) :: m3            ! m^(1/3), m is the mass coordinate

   real(double) :: gradl, dg_ledoux, dsc, gradmu, gradtsc, a, lom_enc, ds_enc, ex_enc, wth_enc
   real(double) :: aper, dap
   real(double) :: eg            ! Convective grad_r - grad_a mixing due to for overshooting

   ! For determining mass transfer rates
   real(double) :: rr, rpr, gmrr, gmrp(2), mtcv(2), q, logq, gq

   real(double) :: ra2s(2), ais(2), oatmb(2)
   real(double) :: et(2), mp(2), pa(2), br(2), zep(2), w1p(2), w2p(2), w3p(2), w4p(2), spp(2), dsp(2), sup(2), dmp(2), brp(2)
   real(double) :: enp(2), mkp(2), phip(2), phisp(2), pap(2), sdp(2), accl(2)
   real(double) :: mdotedd(2), ledd(2), fmdotlost(2), phiratio(2)
   real(double) :: LoLedd

   real(double) :: egr           ! Binding energy density (per unit mass)
   real(double) :: egr0          ! Binding energy density (per unit mass): gravitational energy
   real(double) :: egr1          ! Binding energy density (per unit mass): internal energy
   real(double) :: egr2          ! Binding energy density (per unit mass): recombination energy
   real(double) :: egr3          ! Binding energy density (per unit mass): H2 association energy
   real(double) :: sep           ! Orbital separation

   ! Related to spinup:
   real(double) :: rmin_stream, r_spinup
   real(double) :: wowcrit, wowcritp(1:2)


   ! Initialise local variables:
   deg = 0.0d0
   bca = 0.0d0
   bcp = 0.0d0
   bct = 0.0d0
   bcm = 0.0d0
   bcmb = 0.0d0
   bce = 0.0d0
   bcf = 0.0d0
   bcs = 0.0d0
   bcph = 0.0d0
   m2 = 0.0d0
   xi1 = 0.0d0
   apr = 0.0d0
   atr = 0.0d0

   accl = 0.0d0

   ! Initialise output functions to 0, they may be used in equns1.
   fn1 = 0.0d0



   do Jstar = 1, ktw
   ! Input variables and their changes, star 1 or 2
   af     = var(idx_for_star(VAR_LNF,   Jstar));   daf    = dvar(idx_for_star(VAR_LNF, Jstar))
   at     = var(idx_for_star(VAR_LNT,   Jstar));   dat    = dvar(idx_for_star(VAR_LNT, Jstar))
   vx16   = var(idx_for_star(VAR_O16,   Jstar));   dx16   = dvar(idx_for_star(VAR_O16, Jstar))
   m      = var(idx_for_star(VAR_MASS,  Jstar));   dm     = dvar(idx_for_star(VAR_MASS, Jstar))
   vx1    = var(idx_for_star(VAR_H1,    Jstar));   dx1    = dvar(idx_for_star(VAR_H1, Jstar))
   vqk    = var(idx_for_star(VAR_QK,    Jstar));   dvqk   = dvar(idx_for_star(VAR_QK, Jstar))
   ar     = var(idx_for_star(VAR_LNR,   Jstar));   dar    = dvar(idx_for_star(VAR_LNR, Jstar))
   L      = var(idx_for_star(VAR_LUM,   Jstar));   dL     = dvar(idx_for_star(VAR_LUM,  Jstar))
   vx4    = var(idx_for_star(VAR_HE4,   Jstar));   dx4    = dvar(idx_for_star(VAR_HE4, Jstar))
   vx12   = var(idx_for_star(VAR_C12,   Jstar));   dx12   = dvar(idx_for_star(VAR_C12, Jstar))
   vx20   = var(idx_for_star(VAR_NE20,  Jstar));   dx20   = dvar(idx_for_star(VAR_NE20, Jstar))
   ai     = var(idx_for_star(VAR_INERT, Jstar));   dai    = dvar(idx_for_star(VAR_INERT, Jstar))
   vomega = var(idx_for_star(VAR_OMEGA, Jstar));   domega = dvar(idx_for_star(VAR_OMEGA, Jstar))
   phi    = var(idx_for_star(VAR_PHI,   Jstar));   dphi   = dvar(idx_for_star(VAR_PHI, Jstar))
   phis   = var(idx_for_star(VAR_PHIS,  Jstar));   dphis  = dvar(idx_for_star(VAR_PHIS, Jstar))
   vx14   = var(idx_for_star(VAR_N14,   Jstar));   dx14   = dvar(idx_for_star(VAR_N14, Jstar))
   vx24   = var(idx_for_star(VAR_MG24,  Jstar));   dx24   = dvar(idx_for_star(VAR_MG24, Jstar))
   vx28   = var(idx_for_star(VAR_SI28,  Jstar));   dx28   = dvar(idx_for_star(VAR_SI28, Jstar))
   vx56   = var(idx_for_star(VAR_FE56,  Jstar));   dx56   = dvar(idx_for_star(VAR_FE56, Jstar))
   tam    = var(idx_for_star(VAR_TAM,   Jstar));   dtam   = dvar(idx_for_star(VAR_TAM, Jstar))
   Vmc2   = Vmc2_star(jk, Jstar)

   ! Input variables and their changes, binary orbit
   oa     = var(idx_for_star(VAR_HORB,  Jstar));   doa    = dvar(idx_for_star(VAR_HORB, Jstar))
   ecc    = var(idx_for_star(VAR_ECC,   Jstar));   decc   = dvar(idx_for_star(VAR_ECC, Jstar))
   xi     = var(idx_for_star(VAR_XI,    Jstar));   dxi    = dvar(idx_for_star(VAR_XI, Jstar))
   mb     = var(idx_for_star(VAR_BMASS, Jstar));   dmb    = dvar(idx_for_star(VAR_BMASS, Jstar))
   pmass  = var(idx_for_star(VAR_PMASS, Jstar))

   ! Sanity checks
   assert(vx1 <= 1.0d0)
   assert(vx4 <= 1.0d0)
   assert(vx12 <= 1.0d0)
   assert(vx14 <= 1.0d0)
   assert(vx16 <= 1.0d0)
   assert(vx20 <= 1.0d0)
   assert(vx24 <= 1.0d0)
   assert(vx28 <= 1.0d0)
   assert(vx56 <= 1.0d0)

   ! M is the mass variable
   mt = dm/dt                                                  ! 1e33 g/s
   dmp(Jstar) = dm

   ! Radius and higher powers of the radius
   r2  = exp(2.0d0*ar) - ct(8)
   if (jk == kh) then
      r2 = 0.0d0
      dar = 0.0d0
   end if
   r   = sqrt(r2)
   r3  = r*r2

   ! Compute the rotational period from the rotational frequency
   omega = vomega
   aper  = 2.0*cpi/(abs(omega) * csday)                        ! days
   dap   = -aper/abs(omega) * domega

   ! Calculate the mean-molecular weight (for thermohaline mixing)
   ! Use the values from the previous timestep instead of the current one; this
   ! seems to be more numerically stable
   !> \todo FIXME: now that we have the value of the molecular weight gradient (for the
   !! purpose of semiconvection) we should rework the TH mixing prescription.
   !<
   xh   = vx1  - dx1
   xhe  = vx4  - dx4
   xc   = vx12 - dx12
   xn   = vx14 - dx14
   xo   = vx16 - dx16
   xne  = vx20 - dx20
   xsi  = vx28 - dx28
   xfe  = vx56 - dx56
   xmg  = max(0.0d0, 1.0d0 - xh - xhe - xc - xn - xo - xne - xsi - xfe)
   avmu = 1.0/(0.5 + 1.5*xh + (42.0*xhe + 14.0*xc + 12.0*xn + 10.5*xo + 8.4*xne + 7.0*xmg + 6.0*xsi - 3.0*xfe)/168.0)

   ! set up the composition variables
   xh  = vx1
   xhe = vx4
   xc  = vx12
   xn  = vx14
   xo  = vx16
   xne = vx20
   xmg = vx24
   xsi = vx28
   xfe = vx56
   avmu_new = 1.0/(0.5 + 1.5*xh + (42.0*xhe + 14.0*xc + 12.0*xn + 10.5*xo + 8.4*xne + 7.0*xmg + 6.0*xsi - 3.0*xfe)/168.0)

   ! If magnesium is not solved for by an equation, it is 1-rest
   ! Also, we don't do any alpha captures on magnesium
   if (.not. use_mg24_eqn) then
     xmg = 1.0d0 - xh - xhe - xc - xn - xo - xne - xsi - xfe
     if ( xmg < 0.0d0 ) xmg = 0.0d0
      ! Force O16 to burn to Mg24 directly rather than Si28
      ! Do not consider photodisintegration of Si28 (Si28(g,a)M24)
      eos%ramg = -eos%roo
   end if

   ! Store composition variables in an array as well, for convenience
   xa(1:9) = (/ xh, xhe, xc, xn, xo, xne, xmg, xsi, xfe /)

   ! Equation of state
   call statel ( jk, ji, af, at, xa, Jstar, abund, eos )

   ! Copy EOS output to local variables (for convenience)
   ap    = eos%ap;    arho    = eos%arho; u     = eos%u;     p    = eos%p
   uint  = eos%u1;    urec    = eos%u2;   uass  = eos%u3;
   rho   = eos%rho;   fk      = eos%fk;   t     = eos%t;     sf   = eos%sf
   st    = eos%st;    zt      = eos%zt;   grada = eos%grada; scp  = eos%scp
   rf    = eos%rf;    rt      = eos%rt;   xhi   = eos%xhi;   s    = eos%s
   pr    = eos%pr;    pg      = eos%pg;   pf    = eos%pf;    pt   = eos%pt
   en    = eos%en;    rpp     = eos%rpp;  r33   = eos%r33;   r34  = eos%r34
   rbe   = eos%rbe;   rbp     = eos%rbp;  rpc   = eos%rpc;   rpna = eos%rpna
   rpo   = eos%rpo;   r3a     = eos%r3a;  rac   = eos%rac;   ran  = eos%ran
   rao   = eos%rao;   rane    = eos%rane; rcca  = eos%rcca;  rco  = eos%rco
   roo   = eos%roo;   rgne    = eos%rgne; rgmg  = eos%rgmg;  rccg = eos%rccg
   rpng  = eos%rpng;  ex      = eos%ex;   enx   = eos%enx;   wmu  = eos%wmu
   ramg  = eos%ramg;  rgsi    = 0.0d0
   delta = eos%delta; phie    = eos%phi;  ext   = eos%ext;   fkt  = eos%fkt
   fkr   = eos%fkr;   prandtl = eos%prandtl

   na    = abund%na;    neo   = abund%neo;   nio = abund%nio; nzz = abund%nzz
   ne    = abund%ne;    xai   = abund%xai
   avma  = abund%avm;

   ! Force O16 to burn to Mg24 directly rather than Si28
   ! Do not consider photodisintegration of Si28 (Si28(g,a)M24)
   if (.not. use_mg24_eqn) eos%ramg = -eos%roo

   enp(Jstar) = u + p/rho           ! Enthalpy for this star

   wth = 2.0d0*kth/(kth**4 + 1.0d0) * luminosity_fudge
   ! Thermal energy generation rate
   ds = sf*daf + st*dat
   if (usemenc) then
      if (m<maximum_match_mass .and. m > th(VAR_MASS,ip_mesh_size)) then
         target_daf = get_iph(m, VAR_LNF)-af
         target_dat = get_iph(m, VAR_LNT)-at
         ds_enc = (1.0d0 - impose_entropy_factor)*ds - impose_entropy_factor * (sf*target_daf + st*target_dat)
         wth_enc = entropy_force * impose_entropy_factor*luminosity_fudge
         ex_enc = ex*(1.1d0 - impose_entropy_factor*impose_entropy_factor)/1.1d0
      else
         wth = 1.0
      end if
      enc = 0.0
   end if

   ! Gradient of luminosity, LOM = dL/dM
   lom = ex + en + enx + enc + menc(Jstar, jk) - wth*t*ds/dt
   lom_enc = ex_enc + en + enx - wth_enc*t*ds_enc/dt 
   if (usemenc .and. m<maximum_match_mass .and. m > th(VAR_MASS,ip_mesh_size)) lom=lom_enc

   if ( jk == kh ) then
     ! Approximations for the centre (M=0); the normal expressions become singular
     mor3 = cpi4*rho*c3rd
     lor3 = lom*mor3
     r0   = 1.0d-10
     m3   = cbrt(mor3 * r0**3)
     mp(Jstar) = m3**3
     meff = mor3 * r0**3 * (1.0d0 - omega**2 / (2.*CPI*CG*rho))
   else
     ! Values away from the centre
     mor3 = m/r3
     lor3 = L/r3
     r0   = r
     m3   = cbrt(m)
     mp(Jstar) = m
     meff = m * (1.0d0 - omega**2 / (2.*CPI*CG*rho))
   end if
   r2mu = 3.0d0*mor3**c3rd/(cpi4*rho)                          ! CGS *and* code units

   ! Correction factors FP and FT to the pressure and temperature gradients for
   ! rotating stars.
   ! Limit FP to 0.1, smaller values seem to give convergence problems -SdM
   if (jmod == 1 .or. .not. use_clairaut_distortion) then
      fp = max( 0.1d0, 1.0d0 - 2.0d0*c3rd/( mor3*(cg2*aper)**2 ))
      ft = 1.0d0
      rpole = r
   else
      fp = fp_star(jk, Jstar)
      ft = ft_star(jk, Jstar)

      ! Polar radius, at theta = phi = 0, so that P2(0) = -0.5, P3(0) = 0, P4(0) = 0.375
      rpole = (1.0d0 - 0.5d0 * (a_rot2(jk,Jstar) + a_tide(jk,2,Jstar)) + 0.375d0 * a_tide(jk,4,Jstar)) * &
              r0_over_rv(jk, Jstar) * r 
   end if

   ! Calculate value of the scaled Roche potential at the current grid point.
   ! This is used to describe contact binary components.
   ! For this we need the total mass of the current component, the total mass of
   ! the binary and the orbital separation.
   ! We calculate the potential at the pole of the star, where the expansion of
   ! the potential agrees fairly well with Roche geometry, even for Roche
   ! lobe filling or overfilling stars.
   tm(1) = pmass                                ! Primary mass
   tm(2) = bm - pmass                           ! Secondary mass
   app = oa*oa*bm/(tm(1)*tm(2)*cg1)**2          ! Orbital separation
   xx = calc_scaled_x(tm(1), tm(2), Jstar)      ! Centre-of-mass coordinate of star
   yy = rpole / app                             ! y-coordinate of north pole
   roche_phi(Jstar) = calc_scaled_potential(tm(1), tm(2), xx, yy, 0.0d0) ! Scaled Roche potential
   phi_l1 = calc_scaled_potential(tm(1), tm(2), calc_scaled_xl1(tm(1), tm(2)), 0.0d0, 0.0d0)

   ! pressure gradient equation; gravity corrected for rotation.
   gor  = cg*mor3 * fp                                         ! CGS
   phim = 5.0d21*gor*r2mu                                      ! Code units
   apmu = - rho*phim/p                                         ! Code units
   hp   = min(p/(gor*1.0d11*r0*rho), sqrt(p/(cg*rho*rho)))     ! CGS

   ! temperature gradient equation
   gradr  = 0.25d0*fk*p*lor3/(cpi4*cl*gor*pr) * ft             ! Radiative temperature gradient
   gradmu = expl_var(jk, explv_gradmu, Jstar)                  ! Molecular weight gradient
   gradl  = grada + convection_ledoux * phie*gradmu/delta      ! Critical gradient for Ledoux convection
   dg     = gradr - grada
   dg_ledoux = gradr - gradl

   if (usemenc .and. m<maximum_match_mass .and. m > th(VAR_MASS,ip_mesh_size)) lom=lom_enc

   ! mixing-length theory; the cubic equation for wl/chi has one real solution
   s2  = grada*scp*t
   wc1 = min(0.5d0*s2*(calp*calp*hp/(9.0d0*xhi))**2, 1.0d302)
   wc2 = min(546.75d0*wc1*max(0.0d0, dg_ledoux) + 73.0d0, 1.0d302)
   wc3 = max(cbrt(wc2 + sqrt(wc2*wc2 + 12167.0d0)), 1.0d-200)
   wc4 = max(wc3 - 23.0d0/wc3 - 2.0d0, 0.0d0)
   wcv = wc4*xhi/(3.0d0*calp*hp)
   wl  = calp*hp*wcv

   ! if GRADR <= GRADA, i.e. DG <= 0, we get WCV = 0 and GRAD = GRADR
   grad = gradr - 4.0d0*hp*wcv**3/(calp*s2*xhi)

   ! Heat transport due to semi-convection, after Langer, Sugimoto & Fricke (1983)
   if (csmce > 0.0d0 .and. dg >= 0.0d0 .and. dg_ledoux < 0.0d0) then
     b = pg/p
     a = b*(8. - 3.*b) / (32. - 24.*b - 3.*b**2)
     gradtsc = ( (0.5d0*csmc*a*gradmu + dg_ledoux)**2 + 2.0d0 * csmc*(dg**2 - dg_ledoux*dg -a*gradmu*dg) )/(csmc-2.0d0)**2
     gradtsc = sqrt(gradtsc)
     gradtsc = -gradtsc + gradr + 0.5d0*(2.0d0*dg_ledoux - 2.0*csmc*dg + csmc*a*gradmu)/(csmc-2.0d0)
     grad    = gradtsc
   end if

   ! Heat transport due to thermohaline convection, based on the same idea
   if (dg < 0.0d0 .and. gradmu < 0 .and. cthe > 0.0d0) then
     b = pg/p
     a = b*(8. - 3.*b) / (32. - 24.*b - 3.*b**2)
     grad = grada + 0.5*dg + 6.*cth*gradmu - 0.5*sqrt( (12.*cth*gradmu + dg)**2 - 48.*cth*a*gradmu**2 )
   end if
   gta  = grad - grada

   ! Temperature gradient, dlogT/d(m**2/3)
   atmu = grad*apmu

   ! Temperature scale-height
   htemp = hp / grad

   ! Components of the Brunt-Vaisala frequency
   NT2  = -gor * delta * r / hp * gta
   Nmu2 =  gor * phie  * r / hp * gradmu

   ! mesh spacing equation, with modified pressure, temperature gradient equations
   vp  = ct(4)*ap + ct(5)*log(p + ct(9))
   vpp = ct(4) + ct(5)*p/(p + ct(9))
   vt  = ct(7)*log(t/(t + ct(10)))
   vtt = ct(7)*ct(10)/(t + ct(10))

   ! Extra mesh spacing terms for the AGB
   vp  = vp  + ct(12)*log(p + ct(11)) + ct(14)*log(p + ct(13)) + ct(16)*log(p + ct(15))
   vpp = vpp + ct(12)*p/(p + ct(11))  + ct(14)*p/(p + ct(13))  + ct(16)*p/(p + ct(15))

   ! VR and VM must go to zero at centre like r**2, m**(2/3)
   xmc = ct(6)*cbrt(mc(Jstar)*mc(Jstar))
   vmf = xmc + m3*m3
   if ( jk == kh ) vmf = xmc
   vm   = log(abs(xmc/vmf))
   vmmu = - 1.0d0/vmf
   vr   = - ct(3)*log(r2/ct(8) + 1.0d0)
   vrr  = - ct(3)/(r2 + ct(8))

   ! QQ is the quantity that the meshpoints should be at equal intervals of
   qq = vp + vt + vm + vr
   qm = vpp*apmu + vtt*atmu + vmmu + vrr*r2mu

   ! L dependence for mesh function
   qm  = qm + ct(2)*ex*1.5d-8*m3
   muk = vqk/qm                                                  ! Code units
   mk  = 1.5d0*m3*muk
   mkp(Jstar) = mk

   ! Convert d/dm to d/dk = (d/dm)(dm/dk)
   vmk = vmmu*muk
   vpk = vpp*apmu*muk
   vtk = vtt*atmu*muk
   vrk = vrr*r2mu*muk

   ! potential equation
   phik = phim*muk
   phip(Jstar)  = phi
   phisp(Jstar) = phis

   ! Moment of inertia and specific angular momentum
   ! FIXME: this assumes we are dealing with thin spherical shells, but if the configuration is
   ! distorted we should calculate the moment of inertia differently.
   si  = r2/1.5d0                                              ! 1e22 cm^2
   dsi = 2.d0 * (r2 + ct(8)) * dar / 1.5d0                     ! 1e22 cm^2
   sam = si * omega                                            ! 1e22 cm^2 / s
   dsam = si * domega + dsi * omega                            ! 1e22 cm^2 / s
   samt = dsam / dt                                            ! 1e22 cm^2 / s^2
   vik = si*mk                                                 ! Code units, 1e22 cm^2
   omegat = ( domega + omega* ( dsi - c3rd*r2mu/m3*dm/1.5d0 ) )/dt * mk

   ! Differential rotation
   ! zeta = dlogw/dlogr = dlogw/dk dk/dm dm/dr r = (dlogw/dk) / mk * 4pi r**3 rho
   zeta = (diff_omega(jk, Jstar)/omega)/mk * cpi4*r**3*rho
   alpha_rot = 1.0d0 + 0.5d0 * zeta

   ! Vertical component of the meridional circulation
   Umc2 = 0.0
   Umc4 = 0.0
   if (dg < 0.0d0) Umc2 = -L / (meff * gor*r0*1e11) * P / (rho * scp*T) / dg *&
                              (1.0 - lom * mor3/Lor3 - omega**2 / (2.*CPI*CG*rho)) *  4.*omega**2 / (3. * CG * mor3)

   ! Rotational shear
   ! Based on Mathis, Palacios & Zahn, 2004, A&A, 425, 243 equations (19) and (7)
   ! NB: omega may become negative during solver iterations.
   Dsh = Chs * sqrt( (1.5d-5 / 10.d0) * 1.0d33*r3*abs(omega) * abs(2.d0*Vmc2 - alpha_rot*Umc2) )  ! cm^2 / s
   Dsv = 0.0d0
   if (dg < 0.0d0 .and. Dsh /= 0.0d0) Dsv = Cvs * Dsh * omega**2 / (NT2 * Dsh / (Dsh + xhi) + Nmu2) * zeta**2 / 6.0

   ! Effective diffusion coefficient, from Chaboyer & Zahn 1992, A&A, 253, 178; Mathis & Zahn 2004, A&A, 425, 229
   Deff = 0.0d0
   if (Dsh > 0.0d0) Deff = Cshear * r2 / Dsh * (Umc2 / 30.d0 + Umc4 / 84.d0)

   ! Composition equations:
   ! hydrogen equation with ZAMS baryon correction when KCN = 1
   x1  = vx1
   x1t = (2.0*kx*((1-kcn)*rpp + rpc + rpng + (1-2*kcn)*(rpna + rpo)) + dx1/dt)*mk

   ! helium equation:
   x4  = vx4
   x4t = (4.0*(-kx*(0.5*rpp + rpna + rpo)*(1-kcn)              &
       + ky*(3.0*r3a + rac + 1.5*ran + rao + rane)             &
       - kz*(rcca + rco + roo + rgne + rgmg + rgsi - ramg)) + dx4/dt)*mk

   ! carbon equation:
   x12  = vx12
   x12t = (12.0*(kx*(rpc - rpna) - ky*(r3a - rac)              &
        + kz*(2.0*(rcca + rccg) + rco)) + dx12/dt)*mk

   ! nitrogen equation:
   x14  = vx14
   x14t = (14.0*(kx*(rpna + rpng - rpc - rpo) + ky*ran) + dx14/dt)*mk

   ! oxygen equation:
   x16  = vx16
   x16t = (16.0*(kx*(rpo - rpng) - ky*(rac - rao)              &
        + kz*(rco + 2.0*roo - rgne)) + dx16/dt)*mk

   ! neon equation:
   x20  = vx20
   x20t = (20.0*(ky*(rane - ran - rao) + kz*(rgne - rgmg - rcca)) + dx20/dt)*mk

   ! Magnesium equation:
   x24  = vx24
   x24t = (24.0*(kz*(ramg + rgmg - rgsi - rane - rccg - rco)) + dx24/dt)*mk

   ! TODO: Re-implement the following reaction:
   !   O16(O16,ag)Si28(g,a)Mg24  [ROO]
   ! Replace with
   !   O16(O16,ag)Si28           [ROO]
   !   Mg24(a,g)Si28             [RAMG24]
   ! The rate Mg24(a,g)Si28 is actually available in the nucleosynthesis code
   ! Do we really need the inverse reaction as well? It is not actually difficult to calculate...
   ! (Yes, we need both for temperatures above 1e9; below Mg24(a,g) is (vastly) dominant)
   ! However, Mg24(a,g) is not in Nacre...

   ! Silicon equation:
   x28  = vx28
   x28t = (28.0*kz*(rgsi - roo - ramg) + dx28/dt)*mk

   ! Iron equation:
   x56  = vx56
   x56t = (dx56/dt)*mk

   ! Artificial composition adjustment:
   if (adj_comp) then
     cif = mixing_fudge*impose_composition_factor
     cif = cif**1.5
     tau = min(dt, 1.0d6*csy)
     mm  = min(m, th(4, 1))
     x1t  = (dx1  + cif*(x1  - get_iph(mm, VAR_H1  ))) / tau*mk
     x4t  = (dx4  + cif*(x4  - get_iph(mm, VAR_HE4 ))) / tau*mk
     x12t = (dx12 + cif*(x12 - get_iph(mm, VAR_C12 ))) / tau*mk
     x14t = (dx14 + cif*(x14 - get_iph(mm, VAR_N14 ))) / tau*mk
     x16t = (dx16 + cif*(x16 - get_iph(mm, VAR_O16 ))) / tau*mk
     x20t = (dx20 + cif*(x20 - get_iph(mm, VAR_NE20))) / tau*mk
     x24t = (dx24 + cif*(x24 - get_iph(mm, VAR_MG24))) / tau*mk
     x28t = (dx28 + cif*(x28 - get_iph(mm, VAR_SI28))) / tau*mk
     x56t = (dx56 + cif*(x56 - get_iph(mm, VAR_FE56))) / tau*mk
   end if

   ! Convective mixing coefficient, energy equation, mass-transfer equation, central BCs.
   b = pr/pg
   if ( jk == kh ) then
     wos = 1.0d10
     con = 0.0d0
     vl  = L
     lk  = lom*muk*sqrt(xmc)
     lq  = 0.0d0
     fac = 0.0d0
     xik = 0.0d0
     vm  = m
     vr  = r2
   else
     ! For convective mixing, we don't use the diffusion coefficient calculated above, but rather
     ! use the value based on the assumption that mixing is efficient.
     wos  =(2.5d0+b*(2.0d1+1.6d1*b))*(1.5d0*cu*m3/abs(apmu*m)+1.0d0)
     con  = 6.0d-22*r2*m3/(r2mu*r2mu*muk)*cbrt(wc1)*xhi
     vml  = sqrt(vmf)/m3
     vl   = L*vml
     lk   = vml*mk*(lom - c3rd*xmc/vmf*lor3/mor3)
     lq   = vml*wth*scp*t*apmu*muk*gta
     df0  = cdf*1.0d22*cg*m/r
     fac  = step(phi, df0)
     xik  = cmt*sqrt(pstv(2.0d0*phis, df0))/r*fac*mk
   end if
   if (usemenc .and. m<th(4,1)) lq = 0.0d0

   ! Conversion factor for diffusion coefficients (converting to code units)
   sgf = (cpi4*r2*rho)**2/mk           ! 1e11 g / cm**2

   ! For angular momentum transport
   sgam = (Dsv + wl) * sgf * 1.0e-22   ! Diffusion coefficient, should be in 1e33 g/s for AM transport
   adam = .2d-11*cpi4*r2*rho*Umc2*sam  ! Advection flux, should be in 1e55 g cm**2/s**2
   adam = adam * Cadam * adam_smooth
   if (jmod < 2) adam = 0.0d0

   ! Enforce rigid rotation for the outer few grid points
   if (jk <= 3) sgam = min(sgam, con)
   if (jk >= KH-10) sgam = min(sgam, con)

   ! Fudged convective diffusion coefficient: COS > 0 gives overshooting
   !> \todo FIXME: when Ledoux convection is used, currently does overshooting from
   !! the Ledoux boundary, which doesn't make much sense - it should do
   !! overshooting from the Schwarzschild boundary.
   !<
   deg = cos/wos
   if ( xh < 1.0d-7 ) deg = cps/wos
   eg = dg_ledoux + deg
   ug = pstv(eg, 0.0d0)

   ! The convective mixing coefficient should have cbrt(ug), but ug**2 gives a smoother
   ! transition close to the convective boundary.
   sg = crd*con*ug*ug

   ! Semi-convection, after Langer, Sugimoto & Fricke 1983 (A&A)
   ! CSMC = 0.04 is the effciency of semiconvection, alpha in (10) of LSF83.
   ! Stability condition: DG > 0 (Schwarzschild unstable), DG_LEDOUX<0
   ! (Ledoux unstable)
   Dsc = 0.0d0
   if (dg >= 0.0d0 .and. dg_ledoux < 0.0d0) then
     Dsc = -csmc/6.0d0 * xhi * dg / dg_ledoux
   end if
   ! Convert to code units and add to the overall diffusion coefficient
   ! SGF in 1e11 g/cm**2, DSC in g/cm**2, so SGF*DSC is in 1e-22 [1e33 g/s]
   sg = sg + 1.0d-22*Dsc*sgf

   ! Artificial extra mixing: ARTMIX is the diffusion coefficient in cm^2/s
   ! 0 by default, but can be set from init.dat.
   sg = sg + artmix*sgf

   ! Mixing due to rotational shear
   ! SGF in 1e11 g/cm**2, Deff in g/cm**2, so SGF*Deff is in 1e-22 [1e33 g/s]
   sg = sg + 1.0d-22*Deff*sgf

   ! If we're matching the composition profile to a reference model, do not include convective mixing.
   if (adj_comp) sg = 0.0d0

   ! For ZAMS runs (KCN = 1) always weakly mix inner and outer meshpoints
   if (kcn == 1) then
     if ( jk < 0.075d0*kh .or. jk > 0.8d0*kh ) sg = min(sg, 1.0d-4*crd*con)
   end if
   if ( jk >= (kh - inner_points_mixed) ) sg = min(sg, 1.0d-4*crd*con)
   if ( jk <=       outer_points_mixed  ) sg = min(sg, 1.0d-4*crd*con)

   ! Always weakly mix the next-to-surface meshpoints when gravitational
   ! settling is enabled. Without this the surface abundances are
   ! very unstable and do not properly represent the true surface
   ! values.
   ! Also mix the inner few mesh points, which otherwise show the same type
   ! of problem in some test runs.
   if ( cgrs>0.0d0 ) then
      if ( jk <  1+0.01d0*(kh+1) ) sg = min(sg, 1.0d-2*crd*con)
      if ( jk > kh-0.01d0*(kh+1) ) sg = min(sg, 1.0d-2*crd*con)
   end if

   ! The overall mixing coefficient can be reduced by the solver when convergence is difficult
   sg = sg * mixing_fudge
   wl = sg

   ! Mixing coefficient for thermohaline mixing (except for mu gradient)
   sgth = 0.0d0
   if (cth > 0.0d0) then
     sgth_numer = 16.0 * cpi * r2 * hp * ca * cl * t**3
     sgth_denom = -gta * scp * fk * rho * avmu * abs(mk)
     dth = 0.0d0
     if (sgth_denom /= 0.0d0) dth = sgth_numer/sgth_denom

     ! correct for radiation pressure, see Kippenhahn, Ruschenplatt &
     ! Thomas, A&A 91, 1980 equations (34) and (A6).
     ! Correction factor for ideal gas with radiation is
     ! phi/delta = beta / (4-3beta) = 1/(1+4*B).
     ! NB: B = Prad/Pgas, not beta = Pgas/P
     ! We use the actual phi/delta from the EoS.
     dth = dth * phie / delta

     ! Convert to code mass units (1e33g)
     dth = 1.0d-33 * dth

     ! SGF is in 1e11 g/cm**2, so SGTH is in [1e33 g]/s
     sgth2 = dth*sgf

     ! >0 is 'unphysical'
     sgth = min(cth * sgth2, 0.0d0) * mixing_fudge
   end if

   ! Gravitational settling and radiative levitation
   Frad = 0.0
   if (.not. grs_burgers) then
      ! Use (3.19) from Montmerle & Michaud, 1976 ApJS 31, ignoring the atomic
      ! diffusion term (which should be negligible) and including only the
      ! gravitational settling term.
      ! Reference species is either H or He, depending on which is more abundant
      ! (so this will not work if neither is present anymore).
      nref = 1
      if (na(2)>na(1)) nref = 2
      if (cgrs > 0.0d0 .and. jmod > 2 .and. xa(nref)>0.0d0) then
         ! Compute atomic diffusion coefficients relative to the reference
         ! species for all elements [cm^2/s].
         ! Paquette & al. (1986) equation (40)
         call diffusion_coefficients(rho, t, abund, ddiff, nref)

         ! Compute advection velocities [1e11 cm/s] in the restframe of the
         ! reference particles Montmerle & Michaud (1976) equation (3.19)
         apr = 2 * r * 1.0d-11 * apmu/r2mu     ! dlogP/dr, cgs
         do j=1,9
            gi = na(j)/(na(nref)+na(j))
            kf_num = (1.0_dbl+can(j)*kzn(j)*gi) * (1.0_dbl+0.5_dbl*(kzn(j)+1)*gi)
            kf_denom = (1.0_dbl+can(j)*gi) * (1.0_dbl+0.5_dbl*(kzn(j)+1)*kzn(j)*gi)
            wpi(j) = ddiff(j)*(1.0_dbl+gi)*(cbn(j)-1.0_dbl + kf_num/kf_denom * (can(j)-kzn(j))) / (1.0_dbl+can(j)*gi) * apr
            ! Convert to code units
            wpi(j) = 1.0d-11 * wpi(j)
         end do

         ! Radiative levitation
         Frad = 1e11*L * gradr/grad / (cpi4 * r2)     ! Radiative flux
         if (CRLEV > 0.0d0 .and. sg == 0.0d0) then
            nn(1:9) = na(1:9) / nio                   ! Number fractions
            if (ji == rlev_update) then               ! Only recalculate for "current" values; effectively uses semi-implicit values
               call get_radiative_accelerations(at/cln, arho/cln, dot_product(nn(1:9), can(1:9)), Frad, 9, nn(1:9), kzn(1:9),&
                                                                                                radacc(1:9,jk, Jstar), logkappa)
               radacc(1:9,jk,Jstar) = radacc(1:9,jk,Jstar) / can(1:9)   ! Divide by atomic mass to get acceleration
            end if
            wpi(1:9) = wpi(1:9) + 1.0d-11 * ddiff(1:9) * radacc(1:9,jk,Jstar) * can(1:9) / (CR*T)
         end if


        ! The velocity of the reference particles in the rest frame of the
        ! star is found from the requirement that the total mass flux must
        ! vanish.
        ! NB: Fe and Si are inert, so don't contribute to the mass flux.
        ! They *do* contribute to the sum of all mass fractions and hence
        ! the normalisation of the sum from which the hydrogen velocity is
        ! calculated.
        wpi(nref) = 0
        wpp = 0
        do j=1,9
           wpp = wpp - wpi(j) * xa(j)
        end do
        ! Correct the normalisation for the presence of Fe and Si.
        !wpp = wpp / (1.0-(xfe+xsi))
        ! Now transform all velocities to the rest frame of the star.
        wpi(1:9) = wpi(1:9) + wpp

        ! Export fluxes to equns1
        do j=1,7
           xfs(fn_fh+j-1) = mixing_fudge * cgrs*wpi(j) * rho * xa(j) * cpi4*r2
        end do
     else
        xfs(fn_fh:fn_ffe) = 0.0d0
     end if
   else
      ! Use the treatment of Thoul, Bacall & Loeb, 1994 ApJ 421 to solve the
      ! full set of Burgers equations.
      ! See also Hu & al., 2010 A&A.
      apr = 2.0d0 * r * 1.0d-11 * apmu/r2mu     ! dlogP/dr, cgs
      atr = 2.0d0 * r * 1.0d-11 * atmu/r2mu     ! dlogT/dr, cgs
      Frad = 1e11*L * gradr/grad / (cpi4 * r0**2)        ! Radiative flux
      if (cgrs > 0.0_dbl .and. .false.) then
         ! Using the partial ionisation information available from the EoS,
         ! calculate the actual number of electrons in the plasma and the average
         ! degree of ionisation of each species.
         ! Assume metals are at least singly ionised, because concentrations are
         ! defined wrt electrons and because the Paquette coefficients (Coulomb
         ! scattering) are not valid for a neutral gas.
         avgz(1:9) = dble(kzn(1:9))
         do i = 3, 9
            if ( xa(i) .gt. 1.0d-10 ) then
               avgz(i) = 0.0d0
               do j = 1, kzn(i)
                  avgz(i) = avgz(i) + j*xai(i, j)
               end do
               if (i>2) avgz(i) = max(1.0d0, avgz(i)/na(i))
            end if
         end do
         ! Number of electrons (per baryon)
         ! Can't use output from EoS because that wasn't calculated assuming "proper"
         ! detailed ionisation.
         nne = dot_product(avgz(1:9), na(1:9))
         ! Convert particle number fractions to number densities
         nn(1:9) = na(1:9) * rho/(avma*amu)
         nne = nne * rho/(avma*amu)

         ! Find radiative accelerations
         if (CRLEV > 0.0d0 .and. .false.) then ! .and. sg == 0.0d0) then
            nn(1:9) = na(1:9) / nio                      ! Number fractions
            if (ji == rlev_update) then                  ! Only recalculate for "current" values; effectively semi-implicit values
               call get_radiative_accelerations(at/cln, arho/cln, dot_product(nn(1:9), can(1:9)), Frad, 9, nn(1:9), kzn(1:9),&
                                                radacc(1:9,jk, Jstar), logkappa)
               radacc(1:9,jk,Jstar) = radacc(1:9,jk,Jstar) / can(1:9)   ! Divide by atomic mass to get acceleration
            end if
         else
            radacc(1:9,jk, Jstar) = 0.0d0
         end if

         ! Find terms for the diffusion equations, 9 in total
         call diffusion_terms_burgers(rho,t,9,nn,nne, radacc(1:9,jk, Jstar), avgz,can, diffp,difft,diffg,diffc)
         ! Compute diffusion fluxes and export to equns1
         ! The factor P/Pg corrects for radiation pressure.
         ! TODO: concentration gradient
         wpi(1:9) = 1.0d-11 * (diffp(1:9)*apr + difft(1:9)*atr + crlev*diffg(1:9)) * p/pg
         if (jmod > 2) then
            forall (j=1:9) xfs(fn_fh+j-1) = cgrs*mixing_fudge * rho*xa(j)*cpi4*r2 * wpi(j)
         end if
      else
         xfs(fn_fh:fn_ffe) = 0.0d0
      end if
   end if



   ! Mass transfer terms
   xim  = xi
   dlrk = 0.0d0
   fac1 = 0.0d0
   fac2 = 0.0d0
   v2   = 0.0d0

   if ( Jstar == 2 ) then
     ! Put here the *internal* quantities that depend on *TWIN* variables Quantities
     ! XX needed from both stars should have been stored earlier as XXP(JSTAR) = XX.
     ! Approx. to binary separation, and potential difference between L1, L2
     app = oa*oa*(mp(1) + mp(2))/(mp(1)*mp(2)*cg1)**2
     call potent ( mp(1)/mp(2), dphi )
     df0 = cdf*1.0d22*cg*(mp(1) + mp(2))*dphi/app
     ! Differential mass flux XI, semidetached or in contact: XIK is dXI/dK
     ! Constants CMT, CLT are read in as data. SMF is a smoothing
     ! function to prevent a jump discontinuity at L1 surface.
     fac1 = step(phip(1), df0)
     fac2 = step(phip(2), df0)
     mmk  = 0.5d0*(mkp(1) + mkp(2))
     v2   = (pstv(phisp(1), df0) - pstv(phisp(2), df0))*(fac1 + fac2)
     ! Mass transfer: if XIK>0, transfer from *2->*1, otherwise from *1->*2
     xik  = cmt*v2/(sqrt(abs(v2) + 1.0d11)*app)*mmk * mdot_smooth
     ! Heat transfer due to differential rotation, CLT rad/sec
     dlrk = - clt*(enp(1) - enp(2))*fac1*fac2*mmk

     ! Accretion luminosity
     accl(1) = clac * ( enp(2)+phip(2) - enp(1)-phip(1) ) * max(xit(2) - xit(1), 0.d0)
     accl(2) = clac * ( enp(1)+phip(1) - enp(2)-phip(2) ) * max(xit(1) - xit(2), 0.d0)
   end if

   if ( jk == 1 ) then     ! Surface boundary conditions
      ! pressure, temperature surface boundary conditions
      absl = abs(L)
      bcp = log(fk/(r*1.0d11*gor)*(1.5d0*pg + 0.75d0*pr))
      bct = log(1.0d11*absl/(0.75d0*cpi4*cl*r2*pr))

      fn1(fn_bct) = fn1(fn_bct) + log(1.d0 + accl(1) / absl)
      bct         = bct         + log(1.d0 + accl(2) / absl)

      ! Eddington factor, used in (some) mass loss prescriptions below
      LoLedd = fk*lor3/(cpi4*cl*cg*mor3)

      ! FIXME!  LoLedd can become >1, and we're taking sqrt(1-LoLedd) etc. below
      ! This might happen in brief occasions only, and hence this 'fix' may be good enough
      LoLedd = max(min(LoLedd,1.0_dbl),0.0_dbl)
      ledd(Jstar) = L / LoLedd

      ! Eddington accretion rate, optionally used to limit accretion rate
      mdotedd(Jstar) = cpi4*r*cl/fk * 1.0d-22   ! in 1e33g/s

      ! Store accreted abundances for companion star
      ! For a single star with accretion composition switched on (CCAC=1.0), the
      ! accreted abundance is read from init.dat
      Jstar_accretor = 3-Jstar
      xac(1:9, Jstar_accretor) = xa(1:9)

      ! ----------------------------------------------!
      ! Determine mass loss rate due to stellar winds !
      ! ----------------------------------------------!
      ! Metalicity dependent mass loss scaling
      ! get mass-fraction from baryon-fraction (X1, X4) !SdM
      !     AVM = sum ( mass in amu * baryonfract / baryon nr )
      avm = can(1)*x1/cbn(1)  +&
            can(2)*x4/cbn(2)  +&
            can(3)*x12/cbn(3) +&
            can(4)*x14/cbn(4) +&
            can(5)*x16/cbn(5) +&
            can(6)*x20/cbn(6) +&
            can(7)*xmg/cbn(7) +&
            can(8)*xsi/cbn(8) +&
            can(9)*xfe/cbn(9)
      ! Determine effective surface metallicity i.o.t. scale Mdot with Z_EFF !SdM
      xh_surface  = x1/cbn(1)*can(1)/avm  ! massfraction H surface
      xhe_surface = x4/cbn(2)*can(2)/avm  ! massfraction He surface
      z_eff = 1.d0 - xh_surface - xhe_surface
      z_eff = min(1.d0, max(z_eff, 0.d0))

      ! mass loss rate for dynamo-driven wind (MDTDDW); Alfven radius squared (RA2)
      rht = (0.755d0*absl**0.47d0 + 0.05d0*L**0.8d0)/m**0.31d0
      rad = r/rht
      w2 = r*absl/m
      tet = 120.0d0*r/w2**c3rd*rad**2.7d0
      ro = 1.67d0*abs(aper/tet)
      rot = 1.0d0/(ro**2 + 1.0d0)
      mdtddw = 1.54d-17*w2*rad**2*rot**3.67d0
      bp = 5.4d-3*sqrt(mor3)*w2**c3rd*rad**3.4d0*rot**1.21d0
      if (bp<1.0d-99) bp = 0.0d0
      alfac = 6.2d-08*(r2/mor3*(bp*bp/mdtddw)**2)**c3rd
      ra2 = 0.0d0
      if (chl>0.0d0) ra2 = chl*(r*alfac)**2

      ! Eggleton's Reimers-like wind
      ! cool superwind (Mdot prop. to Lum/(Env. binding energy))
      mdtsw = 0.0d0
      if ( (be(Jstar) /= 0.0d0) .and. (tn(Jstar) /= 0.0d0) ) &
          mdtsw = m*min(1.3d-5*L/abs(be(Jstar)), 1.0d1/(tn(Jstar)*csy))

      ! Use smart mass loss routine that determines the mass loss rate based
      ! on recipes appropriate for the stellar parameters, falling back to the
      ! de Jager rate when no other applicable rate can be found.
      mdot = 0.0
      if (smart_mass_loss>0.0) then
         s_x = can(1)*x1/cbn(1)
         s_he = can(2)*x4/cbn(2)
         s_c = can(3)*x12/cbn(3)
         s_o = can(5)*x16/cbn(5)
         mdot = calculate_mdot(t, L/clsn, m/cmsn, r/crsn, s_x, s_he, s_o, s_c, czs)
         mdot = mdot * cmsn/csy * smart_mass_loss
      else
         ! Mass-loss rate for luminous stars from de Jager et al (1988)
         mdtl = 0.d0
         if (cmj > 0.d0) then
            mdtl = calc_mdot_dejager( log10(T), log10(L/CLSN) )
            if (mdot_errno /= 0) mdtl = 0.0d0
         end if

         ! Mass loss rate for massive stars, Vink et al. (1999, 2000, 2001)
         mdot_jvink = 0.0d0
         if (cmv > 0.0d0) mdot_jvink = cmv*calc_mdot_vink(m/cmsn, L/clsn, t, clogz)

         ! Mass loss rate for Wolf-Rayet stars, Nugis&Lamers 2000
         mdot_wr = 0.0d0
         if (cmnl > 0.0d0) then
            s_x = can(1)*x1/cbn(1)
            s_he = can(2)*x4/cbn(2)
            s_c = can(3)*x12/cbn(3)
            s_o = can(5)*x16/cbn(5)
            mdot_wr = cmnl*calc_mdot_wr_nl(log10(L/clsn),log10(t),s_x,s_he,s_c,s_o,clogz)
         end if

         ! Scale with the effective metallicity: (as Vink et al)
         ! The Vink rate and the Nugis & Lamers rate already depend on the
         ! metallicity
         if (zscaling_mdot /= 0.d0) then
            mdtl = mdtl*(z_eff/czsn)**zscaling_mdot
         end if

         ! The expected mass loss rate is the maximum of the different recipes
         ! Convert from Msun/yr to code units (1e33g/s)
         mdot = max(cmj*mdtl, mdot_jvink, mdot_wr) * cmsn/csy
     end if

      ! Total wind mass loss rate, in 1e33 g/s
      zep(Jstar) = cml*mdtddw + max(mdot, cmr*mdtsw)
      if ( mutate .or. jmod == 0 ) zep(Jstar) = 0.0d0
      zep(Jstar) = zep(Jstar) * mdot_smooth * cmdot_wind

      !     Rotation rate over the critical rotation rate, simply defined by
      !     setting gravitational acceleration equal to the centrifugal one,
      !     taking into account the Eddington factor (e.g. Heger, Langer,
      !     Woosley 2000, eq 50)
      wowcrit = omega/sqrt(cg*mor3*(1.0d0-LoLedd))
      wowcritp(Jstar) = wowcrit

      !     Rotational enhanced mass loss, according to options:
      !      1. Heger, Langer & Langer (2000) [they cite Friend & Abbott 1986, but
      !         the expression does not appear in that paper]
      !      2. Maeder & Meynet (2000) (Paper VI from their series)
      if (cmdotrot_hlw>0.0d0) then
         zep(Jstar) = zep(Jstar) / (1.0d0 - wowcrit)**0.43_dbl
      end if
      !     Rotational enhanced mass loss from Maeder & Meynet (2000)
      !     Empirical force multipliers, from Lamers et al. 1995
      if (cmdotrot_mm>0.0d0) then
         alpha = 0.52
         if (at < 4.30*cln) then
            alpha = -0.71 + 0.22*at/cln
         else if (at > 4.30*cln .and. at < 4.35*cln) then
            alpha = 0.236+0.284*(at/cln - 4.30)/0.05
         else if (at > 5.0*cln) then
            alpha = 1.0
         end if
         zep(Jstar) = zep(Jstar) * ( (1.0_dbl-LoLedd) / abs(1.0_dbl - (2*wowcrit*c3rd)**2 - LoLedd) )**(1.0d0/alpha-1.0d0)
      end if

      ! Artificial mass loss/gain - rate depends on configuration
      mdot_cmi = 0.0d0
      if ( jmod /= 0 .and. cmi /= 0.0d0 ) then
         select case (cmi_mode)
            case (1) ! Exponential
               mdot_cmi = mdot_smooth*cmi*mp(1)
               if ((m+mdot_cmi*dt) > uc(13)*cmsn .and. cmi > 0.0d0) mdot_cmi = max(1.0d-9*cmsn/csy, (m-uc(13)*cmsn)/dt)
               if (m > uc(13)*cmsn .and. cmi > 0.0d0) mdot_cmi = 0.0d0
            case (2) ! Linear
               mdot_cmi = mdot_smooth*cmi*cmsn
               if ((m+mdot_cmi*dt) > uc(13)*cmsn .and. cmi > 0.0d0) mdot_cmi = max(0.0d0, (m-uc(13)*cmsn)/dt)
         end select
         if (jmod < 10) mdot_cmi = mdot_cmi * jmod / 10.d0
      end if

      ! Adjust luminosity for the kinitic energy of the wind (photon tiring)
      ! FIXME: this should be part of the energy equation (LOM) above
      f_photon_tiring = cphotontire*min(1.0d22*gor/lor3 * zep(Jstar)/r, 1.0d0)
      L = (1.0d0 - f_photon_tiring) * L
      vl = L

      ! Transport wind mass loss in Mo/yr to printb
      wml = -zep(Jstar)*csy/cmsn

      ! Determination of mass of current component (M), and other cpt (OM).
      if ( ktw == 1 ) then
         if ( jb == 1 ) om = mb - m               ! M2; M = M1
         if ( jb == 2 ) then
            zep(Jstar) = 0.0d0
            ra2 = 0.0d0
            bp = 0.0d0
            fdt = age*csy + dt - t0
            if ( fdt.le.0.0d0 ) fdt = dt
            m2 = om0 + fdt*omta
            oa = a0 + fdt*ata
            ecc = e0 + fdt*eta
            om = m0 + fdt*mta                     ! M1; M = M2
            mb = m2 + om
         end if
      end if
      if ( ktw == 2 ) om = mb - m
      ! end of M, OM determination

      bm = mb
      horb = oa
      secc = ecc
      ! quantities depending on the orbit
      mu = m*om/mb
      e2 = ecc*ecc
      ef = 1.0d0/sqrt(1.0d0 - e2)
      we25 = ef**5
      sep = (oa*ef/(cg1*mu))**2/mb
      bper = sqrt(sep**3/mb)/cg2
      rl = sep*rlobe(cbrt(m/om))
      ! quantities for spin up/down due to mass transfer -SdM
      ! impact parameter for a mass transfer stream: following Ulrichj&Bulger76
      rmin_stream = sep * 4.25d-2*( (m/om)*(1.0d0+(m/om)) )**0.25 !if current star is accretor

      if (rmin_stream>r) then
         !     Stream misses star -> accretion disk formed Assume AM is accreted
         !     rom the inner kepler orbit of the disk, which has radius equal to
         !     stellar radiusr
         r_spinup = r
      else
         !     Stream hits the star, accretion by impact.  According to U&B76 a
         !     disk would have formed with radius of 1.7*rmin_stream, if the
         !     accreting star was not in its way.
         r_spinup = 1.7*rmin_stream
      end if

      ! Wind accretion, Bipolar reemission,  AM transfer
      ! -------------------------------------------------
      ! IMPROVE THIS
      if ( ktw == 1) then
         pap(2) = cpa*(rlobe(cbrt(om/m)))**2
         brp(2) = cbr
         sup(2) = csu*oa
         sdp(1) = csd*oa
      else ! If doing TWIN
         ! PAP: Partial Accretion of the wind of the companion:
         ! [TODO We should replace this by Bondi-Hoyle accretion or something like that]
         pap(Jstar) = cpa*(r/sep)**2

         ! BRP: Bipolar Reemision: i.e. material accreted through RL overflow can
         ! be reemitted again along the poles. This is the way in TWIN to do
         ! non-conservative mass transfer. It assumed (later on in the code) that
         ! the material lost has the specific angular momentum of the orbit of
         ! the accreting star. For comparison: beta (as for example in de mink
         ! etal 2008) = Maccr / Mtransfer = (1-BRP).
         ! [Todo: make this a parameter in init.dat]
         brp(Jstar) = 0.0d0

         ! SUP: Spin Up: specific angular momentum of the material accreted,
         ! which is added to the spin of the accreting star.  Warning:
         ! Switching this off (CSU=0) can give strange effects in wider
         ! binaries: accreting mass but conserving spin leads to slower
         ! rotation, opposite to what is realistic. Only current solution is
         ! setting the efficency of spin orbit coupling very high to prevent
         ! this. Switching it on also has unwanted (but physical) effects: the accreting
         ! star will reach critical rotation, especially in wider system where the tides are ineffective
         sup(Jstar) = csu*cg1*sqrt(m*r_spinup)

         ! SDP: Spin Down: specific angular momentum lost from the donor star.
         ! The specific angular momentum lost due to RLOF is that of a disk with the
         ! radius of the star, *not* the specific angular momentum of a sphere with the
         ! radius of the star.
         sdp(Jstar) = csd*(cg1*r**2)/(cg2*aper)
      end if

      ! Radius for determination of RLOF: average the star's current radius and
      ! its radius on the previous timestep to smooth out kinks in the mass
      ! transfer rate.
      rpr = sqrt(exp(2.0d0*h(VAR_LNR, 1)) - ct(8))
      rr  = 0.5d0 * (r + rpr)

      ! effective gravity (centrifugal and gravity) at the surface - efective
      ! gravity at the Roche-lobe surface (assuming synchronous rotation of
      ! stars with the rotation of orbit)
      gmr  = 1.0d22*cg*(mb*c3rd/sep**3*(r2 - rl*rl) + m*(rl - r)/(r*rl))
      gmrr = 1.0d22*cg*(mb*c3rd/sep**3*(rr*rr - rl*rl) + m*(rl - rr)/(rr*rl))


      ! GR terms for circularisation and ang. mom. loss
      rgr   = cgrt*mu/(mb*bper)*(sep/bper)**5*we25
      etgr  = rgr/9.6d1*(3.04d2 + 1.21d2*e2)*ecc
      oatgr = cgw*rgr*(1.0d0 + 0.875d0*e2)*oa


      ! Magnetic braking/1d50 according to Rappaport, Verbunt & Joss, 1983, using multiplier CMB:
      oatmb(Jstar) = 0.0_dbl
      if(cmb > 0.0_dbl) then
         oatmb(Jstar) = cmb * 3.8d-3 * m * r**4 * (cpi/(aper*4.32d4))**3  ! Exponents in constants: -30+33+4*11-50 = -3

         ! Change MB strength depending on mass ratio of the convective envelope:
         if(qcnv > 0.d0 .and. qcnv < 0.02d0)  oatmb(Jstar) = oatmb(Jstar)*exp(1.d0 - 2.d-2/qcnv)
         if(qcnv < 1.d-9) oatmb(Jstar) = 0.d0       ! No MB when no convective envelope
         if(qcnv > 1.d0-1.d-9) oatmb(Jstar) = 0.d0  ! No MB when star fully convective
      end if


      ! BC for surface potential eigenvalue
      ! BCPH = PHIS + GMRR
      bcph = phis + gmr

      ! BC for the potential on L1 surface
      bcf = phi + gmr

      ! Primary mass
      bcmp = pmass - mp(1)

      ! Ratio of difference in potential with RLOF surface and potential at surface
      gmrp(Jstar) = gmr
      phiratio(Jstar) = 1.0d-22*gmr / (cg*mor3*r2 + cg*mb*r2/sep**3)

      ! The cheaper model for RLOF, used if CMS > 0
      rlfac = ar - log(rl)
      if ( cms > 0.0d0 ) xi = cms*(pstv(rlfac, 0.0d0))**3
      if ( cms > 0.0d0 .and. Jstar == 1 ) xi1 = xi
      if ( cms > 0.0d0 .and. Jstar == 2 ) xi = xi1

      ! Mass transfer rate for cataclismic binaries, after Ritter 1988
      mtcv(Jstar) = 0.0d0
      if (r < rl .and. cmt_cv > 0.0d0) then
         q = tm(3-Jstar) / tm(Jstar)
         logq = log(q)
         if (q < 1.0) then
            gq = 0.954 + 0.025*logq - 0.038*logq**2
         else
            gq = 0.954 + 0.039*logq + 0.114*logq**2
         end if
         ! Based on eq. (9) in Ritter 1988
         mtcv(Jstar) = cmt_cv * rho*sqrt(CR*T / avmu) * gor*r/omega**2 * (hp/rl)**3 * exp(-(rl - r)*1.0d11*gq / hp - 0.5d0)
      end if

      ! Model for mass transfer/RLOF in contact
      delta_phi(1) = max(phi_l1 - roche_phi(1), 0.d0)
      delta_phi(2) = max(phi_l1 - roche_phi(2), 0.d0)
      if (use_contact_flow) xi = 1.0d0*(mdot_rlof(2)*delta_phi(2)**2.5 - mdot_rlof(1)*delta_phi(1)**2.5) * mdot_smooth

      mtr = -xi*csy/cmsn  ! Transport to printb

      ! equilibrium-tide model for tidal friction
      ww = 0.4d0*r2*m/(abs(ai) + 1.0d-10)
      gam = ctf*(2.0d0/(ww + ww**3.2d0) + rad**8)

      ! Use CSO as spin-orbit coupling (tidal friction) switch:
      tfr = cso*9.0d0*gam/(3.0d22*m*r2/L)**c3rd*(r/sep)**8*om*mb/(m*m)
      we5 = we25*we25
      we65 = ef**1.3d1
      we2 = (1.0d0 + e2*(7.5d0 + e2*(5.625d0 + e2*0.3125d0)))*we65
      we3 = (1.0d0 + e2*(3.0d0 + e2*0.375d0))*we5
      wea = (5.5d0 + e2*(8.25d0 + e2*0.6875d0))*we5
      web = (9.0d0 + e2*(33.75d0 + e2*(16.875d0 + e2*0.703d0)))*we65
      w1p(Jstar) = om/(mb*m)
      w2p(Jstar) = tfr*we2
      w3p(Jstar) = mu*sep**2/(ai*ef)*tfr*we3
      w4p(Jstar) = (ra2 + 2.0d0*c3rd*r2)/ai

      ! Total spin angular momentum and change in spin angular momentum.
      ! These expressions are different depending on whether differential rotation
      ! is taken into account or not.
      if (rigid_rotation) then
         ! log(Iw) = dI/I - dP/P, angular momentum change for solid body rotation.
         dsp(Jstar) = dai/ai - dap/aper
         ! CG1/CG2 = 2pi/day * 1e5, 2pi/(APER*day) = omega, in s^-1.
         spp(Jstar) = cg1*ai/(cg2*aper)
      else
         !> \todo FIXME: We're not actually interested in the total angular momentum here,
         !! what we're interested in is that the angular-momentum transport through
         !! the surface works out. In other words, the rate of change of the surface
         !! angular momentum should enter the BC. Angular-momentum transport in the
         !! interior then has to make sure that angular momentum is transported back
         !! to the surface.
         !! Global conservation of angular momentum is built into the expressions by
         !! construction (also in the case of advection?) and so *should* be ok.
         !! Caveat: the code should not crash in the situation where the spin
         !! direction of the surface is reversed (APER becomes negative). Note the
         !! way in which APER becomes negative: first APER->+inf, then APER->-inf and
         !! finally -inf<APER<0. For this reason, it may be better to switch to
         !! 1/APER as an independent variable.
         !<
         dsp(Jstar) = dtam / ( tam+1.0d-16)    ! For differential rotation
         spp(Jstar) = cg1 * dtam / cg2
      end if
      et(Jstar) = ecc*tfr*(web - bper/aper*wea )
      if ( ji == 0 ) then
         rlf(Jstar) = rlfac
         hspn(Jstar) = spp(Jstar)
      end if

      if ( ktw == 1 ) then ! Non-twin case
         ! BCs for masses, spins, orbital AM and eccentricity in non-TWIN case
         if ( jb == 1 ) then
            !     ZET(1): wind mass loss that is not accreted by companion
            zet(1) = (1.0d0 - pap(2))*zep(1)
            !     XIT(1): mass transferred to companion (*2) by RLOF + Wind accretion ?
            xit(1) = xi + pap(2)*zep(1)
            !     ZET(2): mass emitted from pole of companion (after accretion)
            zet(2) = brp(2)*xit(1)
            !     BCM: boundary condition for total mass of *1
            bcm = dmp(1) + (zet(1) + xit(1) - mdot_cmi) * dt
            ! The accreted mass flux for the star, used to accrete
            ! material of a given composition.
            ! CCAC is a switch that determines whether this composition
            ! change is taken into account or neglected.
            if (jmod>0) fn1(fn_macc) = ccac*max(mdot_cmi - zet(1) - xit(1), 0.0d0)
            bcs = spp(1)*(dsp(1)/dt + w3p(1) + w4p(1)*zet(1)) - w2p(1)*oa + oatmb(1)
            ! For differential rotation
            !BCS = ( 0.0*DXVS(VAR_TAM) / DT - SI/APER**2*DAP/DT + SI/APER*ZET(1) )
            if (.not. rigid_rotation) bcs = samt*mk - sam*zet(1)
            bca = doa/dt + (w1p(1)*zet(1) + w2p(1))*oa - w3p(1)*spp(1)&
                 + oatgr + m/(om*mb)*zet(2)*oa
            bce = decc/dt + et(1) + etgr
            bcmb = dmb + (zet(1) + zet(2))*dt
            ! Transport different dHorb/dt's to printb - case not TWIN
            ! SP:spin (wind+mag.br.), MB:mag.br., SO:spin-orbit, ML:system mass loss, GW:grav.waves
            dhsp(1) = w4p(1)*spp(1)*zep(1) ! + OATMB(1) ! ?? Onno
            dhmb(1) = ra2*spp(1)*zep(1)/ai + oatmb(1)
            !           DHMT(1) = M/(OM*MB)*ZET(2)*OA
            dhso(1) = (w2p(1)*oa - w3p(1)*spp(1))
            dhml    = w1p(1)*zep(1)*oa
            dhgw    = oatgr
         else
            bcm = m - m2
         end if
      end if

      ! Total orbital angular momentum loss
      dhdt = doa/dt

      ! Need RA2 and AI for both stars to calculate mag.br.
      ra2s(Jstar) = ra2
      ais(Jstar)  = ai

      !     Boundary conditions depending on both components
      if ( Jstar == 2 ) then
         ! Put here the BCs that depend on *both* components (TWIN case, KTW = 2).
         ! For PASW and BPRE, determine ZET, XIT from ZEP, XI, PAP, BRP
         dmt0 = 0.001d0*mp(1)/(tn(1) + 1.0d0) + mp(2)/(tn(2) + 1.0d0)
         do j = 1, 2
            !: PA: Partial accretion of the wind of the companion (if that wind is
            !stronger than its own wind)
            pa(j) = pap(j)*step(zep(j) - zep(3 - j), 0.0d0)
         end do
         !: XI: mass accreted from RLOF + partial wind accretion
         do j = 1, 2
            xit(j) = pstv((3 - 2*j)*xi, 0.0d0) + pa(3 - j)*zep(j) + mtcv(j) - mtcv(3-j)
         end do
         ! Optionally limit the mass transfer rate
         do j = 1, 2
            xit(j) = sign(min(abs(xit(j)),mtr_limit*cmsn/csy), xit(j))
         end do
         ! ----------------------------------------------------------------------------
         !        Non conservative mass transfer
         ! ----------------------------------------------------------------------------
         !        There are several ways to treat non-conservative mass transfer:
         !
         !         1. The mass accretion rate can be limited by the Eddington rate
         !            of the accretor (caveat: it is not clear whether this applies
         !            when accreting from a disk rather than through direct impact of
         !            the tidal stream; naively not, but it depends on whether the
         !            accretion disk becomes thick or not).
         !         2. Angular momentum limit: quench mass accretion rate if the star
         !            is close to critical rotation. This really only works if the star
         !            is spun up by the accretion.
         !         3. Radiative feedback from the hot-spot, in case of a direct impact.
         !         4. Just set a constant fraction of material to be accreted
         !
         !        For now, we implement (1), (2) and (4)
         ! ----------------------------------------------------------------------------

         !        Eddington-limited accretion (based on Beer, Dray, King & Wynn 2007)
         !        Fraction of transferred material lost (confusingly, this is sometimes
         !        called "beta" and sometimes "1-beta", although Beer & al. call it psi).
         do j = 1, 2
            fmdotlost(j) = 0.0d0
            if (xit(3-j) >= mdotedd(j) .and. cmtel > 0.0d0) then
               !        Cannot accrete all, lose a bit
               !                  FMDOTLOST(J) = 1.0d0 - MDOTEDD(J)/XIT(3-J)
               !        This is expression (7) in Beer & al, but with factors rearranged and
               !        using the potential at the Roche lobe surface in lieu of the potential
               !        in L1 (these are similar and ought to be the same anyway).
               fmdotlost(j) = (1.0d0-ledd(j)/(abs(gmrp(j)*xit(3-j))))
               fmdotlost(j) = fmdotlost(j) * phiratio(j)
               fmdotlost(j) = min(1.0d0, max(fmdotlost(j), 0.0d0))
            end if
         end do

         !        Combine all contributions to non-conservative mass transfer to get the
         !        total fraction of material ejected.
         !        Get total fraction of mass accreted by multiplying all terms that may
         !        contribute to non-conservative mass loss (alternative: use the most
         !        restrictive).
         !        CBR stands for "Bipolar reemission" and is a constant in init.dat
         do j = 1, 2
            brp(j) = (1.0d0 - cbr)
            brp(j) = brp(j)*(1.0d0 - cmtel * fmdotlost(j))
            brp(j) = brp(j)*(1.0d0 - cmtwl * wowcritp(Jstar))
            brp(j) = 1.0d0 - brp(j)
         end do

         ! Calculate fraction of material ejected from the system
         forall(j = 1:2) br(j) = brp(j)*step(xit(3 - j) - xit(j), 0.0d0)

         ! Mass ejected from the system: the fraction of the wind that is not accreted
         ! by the companion and the contribution from non-conservative Roche lobe
         ! overflow
         forall(j = 1:2) zet(j) = (1.0d0 - pa(3 - j))*zep(j) + br(j)*xit(3 - j)

         ! The net mass accretion flux for each of the two components
         if (jmod>0) then
            fn1(fn_idx_for_star(fn_macc, 1)) = ccac*max(mdot_cmi - zet(1) + xit(2) - xit(1), 0.0d0)
            fn1(fn_idx_for_star(fn_macc, 2)) = ccac*max(mdot_cmi - zet(2) - xit(2) + xit(1), 0.0d0)
         end if

         ! Mass boundary condition for each of the two components...
         fn1(fn_bcm) = dmp(1) + (zet(1) + xit(1) - xit(2))*dt ! *1
         bcm         = dmp(2) + (zet(2) - xit(1) + xit(2))*dt ! *2
         ! ... and for the binary
         bcmb  = dmb + (zet(1) + zet(2))*dt ! binary

         ! AM boundary condition for each of the two components ...
         fn1(fn_bcs) = spp(1)*(dsp(1)/dt + w3p(1) + w4p(1)*zet(1))-w2p(1)*oa - sup(1)*xit(2) + sdp(1)*xit(1) + oatmb(1)    ! *1
         bcs         = spp(2)*(dsp(2)/dt + w3p(2) + w4p(2)*zet(2))-w2p(2)*oa - sup(2)*xit(1) + sdp(2)*xit(2) + oatmb(2)    ! *2
         ! ... and for the binary
         bca   = doa/dt + (w1p(1)*zet(1) + w2p(1))*oa - w3p(1)*spp(1) + (w1p(2)*zet(2) + w2p(2))*oa - w3p(2)*spp(2) + oatgr &
               + (sup(1) - sdp(2))*xit(2) + (sup(2) - sdp(1))*xit(1) ! orbit

         ! Eccentricity boudary condition
         bce   = decc/dt + et(1) + et(2) + etgr

         ! Export different contributions to dHorb/dt's to printb:
         ! SP:spin (wind+mag.br.), MB:mag.br., SO:spin-orbit, ML:system mass loss, GR:grav.waves
         dhsp(1) = w4p(1)*spp(1)*zet(1)
         dhsp(2) = w4p(2)*spp(2)*zet(2)
         dhmb(1) = ra2s(1)*spp(1)*zet(1)/ais(2) + oatmb(1)
         dhmb(2) = ra2s(2)*spp(2)*zet(2)/ais(2) + oatmb(2)
         dhso(1) = (w2p(1)*oa - w3p(1)*spp(1))
         dhso(2) = (w2p(2)*oa - w3p(2)*spp(2))
         dhml = (w1p(1)*zet(1) + w1p(2)*zet(2))*oa
         dhgw = oatgr
      end if      ! if ( Jstar == 2 )
   end if      ! End of surface boundary condition block:  if ( jk == 1 ), ~500 lines up - should this be a subroutine?

   ! Return output
   ! Stellar terms
   fn1(fn_idx_for_star(fn_bcp, Jstar)) = bcp
   fn1(fn_idx_for_star(fn_bct, Jstar)) = bct
   fn1(fn_idx_for_star(fn_vp, Jstar)) = vp
   fn1(fn_idx_for_star(fn_vpk, Jstar)) = vpk
   fn1(fn_idx_for_star(fn_vr, Jstar)) = vr
   fn1(fn_idx_for_star(fn_vrk, Jstar)) = vrk
   fn1(fn_idx_for_star(fn_vt, Jstar)) = vt
   fn1(fn_idx_for_star(fn_vtk, Jstar)) = vtk
   fn1(fn_idx_for_star(fn_vl, Jstar)) = vl
   fn1(fn_idx_for_star(fn_lk, Jstar)) = lk
   fn1(fn_idx_for_star(fn_lq, Jstar)) = lq
   fn1(fn_idx_for_star(fn_mt, Jstar)) = mt
   fn1(fn_idx_for_star(fn_vm, Jstar)) = vm
   fn1(fn_idx_for_star(fn_vmk, Jstar)) = vmk
   fn1(fn_idx_for_star(fn_sg, Jstar)) = sg
   fn1(fn_idx_for_star(fn_x1, Jstar)) = x1
   fn1(fn_idx_for_star(fn_x1t, Jstar)) = x1t
   fn1(fn_idx_for_star(fn_x16, Jstar)) = x16
   fn1(fn_idx_for_star(fn_x16t, Jstar)) = x16t
   fn1(fn_idx_for_star(fn_x4, Jstar)) = x4
   fn1(fn_idx_for_star(fn_x4t, Jstar)) = x4t
   fn1(fn_idx_for_star(fn_x12, Jstar)) = x12
   fn1(fn_idx_for_star(fn_x12t, Jstar)) = x12t
   fn1(fn_idx_for_star(fn_x20, Jstar)) = x20
   fn1(fn_idx_for_star(fn_x20t, Jstar)) = x20t
   fn1(fn_idx_for_star(fn_bcm, Jstar)) = bcm
   fn1(fn_idx_for_star(fn_vi, Jstar)) = ai
   fn1(fn_idx_for_star(fn_vik, Jstar)) = vik
   fn1(fn_idx_for_star(fn_phi, Jstar)) = phi
   fn1(fn_idx_for_star(fn_phik, Jstar)) = phik
   fn1(fn_idx_for_star(fn_bcf, Jstar)) = bcf
   fn1(fn_idx_for_star(fn_bcs, Jstar)) = bcs
   fn1(fn_idx_for_star(fn_bcph, Jstar)) = bcph
   fn1(fn_idx_for_star(fn_x14, Jstar)) = x14
   fn1(fn_idx_for_star(fn_x14t, Jstar)) = x14t
   fn1(fn_idx_for_star(fn_avmu, Jstar)) = avmu
   fn1(fn_idx_for_star(fn_sgth, Jstar)) = sgth
   fn1(fn_idx_for_star(fn_omega, Jstar)) = omega
   fn1(fn_idx_for_star(fn_omegat, Jstar)) = omegat
   fn1(fn_idx_for_star(fn_sam, Jstar)) = sam
   fn1(fn_idx_for_star(fn_samt, Jstar)) = samt * mk
   fn1(fn_idx_for_star(fn_adam, Jstar)) = adam
   fn1(fn_idx_for_star(fn_sgam, Jstar)) = sgam
   fn1(fn_idx_for_star(fn_si, Jstar)) = si
   fn1(fn_idx_for_star(fn_am, Jstar)) = tam        ! Code units
   fn1(fn_idx_for_star(fn_amk, Jstar)) = vik/aper  ! Code units

   fn1(fn_idx_for_star(fn_fU2k, Jstar)) = Umc2 * rho * r**2       ! 1e22 g/s
   fn1(fn_idx_for_star(fn_fv, Jstar)) = cpi4 * r / (6.d0 * mk)    ! 1e-22 cm/g
   fn1(fn_idx_for_star(fn_Vmc2, Jstar)) = Vmc2                    ! cm/s
   fn1(fn_idx_for_star(fn_Umc2, Jstar)) = Umc2                    ! cm/s

   fn1(fn_idx_for_star(fn_rho, Jstar)) = rho                      ! g/cm^3
   fn1(fn_idx_for_star(fn_p, Jstar)) = p                          ! dyne/cm^2; erg/cm^3
   fn1(fn_idx_for_star(fn_hp, Jstar)) = hp                        ! cm

   fn1(fn_idx_for_star(FN_X24, Jstar)) = x24
   fn1(fn_idx_for_star(FN_X24t, Jstar)) = x24t
   fn1(fn_idx_for_star(FN_x28, Jstar)) = x28
   fn1(fn_idx_for_star(FN_X28t, Jstar)) = x28t
   fn1(fn_idx_for_star(FN_X56, Jstar)) = x56
   fn1(fn_idx_for_star(FN_X56t, Jstar)) = x56t

   fn1(fn_idx_for_star(FN_PHI_ROCHE, Jstar)) = roche_phi(Jstar)
   fn1(fn_idx_for_star(FN_ENT, Jstar)) = enp(Jstar)

   ! Diffusion fluxes
   forall (i=fn_fh:fn_ffe)
      fn1(fn_idx_for_star(i, Jstar)) = xfs(i)
   end forall

   ! Binary functions
   fn1(fn_idx_for_star(fn_bca, Jstar)) = bca
   fn1(fn_idx_for_star(fn_bce, Jstar)) = bce
   fn1(fn_idx_for_star(fn_xim, Jstar)) = xim
   fn1(fn_idx_for_star(fn_xik, Jstar)) = xik
   fn1(fn_idx_for_star(fn_dlrk, Jstar)) = dlrk
   fn1(fn_idx_for_star(fn_bcmb, Jstar)) = bcmb
   fn1(fn_idx_for_star(fn_pmass, Jstar)) = bcmp

   ! Pass through the output from the equation of state calculation
   if (present(eosout)) eosout = eos
   if (present(abundout)) abundout = abund

   ! Store variables for explicit calculations
   if (ji <= 0) then
     expl_var(jk, explv_avmu, Jstar) = avmu_new
     expl_var(jk, explv_logp, Jstar) = ap
   end if

   ! ------------------------------------------------------------------------------
   ! Store computed values if called from PRINTB or REMESH (JI<0) or when
   ! simply calling function values for the current timestep (JI=0)
   ! ------------------------------------------------------------------------------
   if ( ji < 0 .and. present(px)) then
      if ( Jstar == ktw ) then
        ! sundry numbers for PRINTB, FUNCS2; SX(40 to 45,JK) for TWIN variables
        px(48) = sep
        px(45) = fac1
        px(44) = fac2
        px(43) = v2
        px(42) = xik
        px(41) = enp(1) - enp(2)
        px(40) = - dlrk
      end if

      ! do. for non-TWIN variables
      
      egr0 = 1.0d22*gor*r2  ! Gravitational binding energy
      egr1 = -uint          ! Internal binding energy
      egr2 = -urec          ! Recombination binding energy
      egr3 = -uass          ! H2 association binding energy
      egr  = egr0 - u       ! Total binding energy (u=u1+u2+u3)

      ! Thermal energy generation rate
      px(19) = -t*(sf*daf + st*dat)/dt + scp*t*apmu*gta*mt/(1.5d0*m3)

      ! Eddington rate
      px(29) = LoLedd

      px(32) = fp
      px( 7) = grad
      px( 8) = dg
      px(30) = wl
      px(35) = wcv
      px(39) = eg

      ! Pass binding energies to output
      px(80) = egr0
      px(81) = egr1
      px(82) = egr2
      px(83) = egr3
      px(84) = egr

      ! Output needed for remesh
      px(85) = qq
      px(86) = qm
      px(87) = phim
      px(88) = gmr
      px(89) = m3

      !> \todo FIXME: this does NOT work properly in TWIN mode. May not be a problem because this is only used by printb...
      ! Reaction rates
      px(50) = rpp
      px(51) = rpc
      px(52) = rpng
      px(53) = rpna
      px(54) = rpo
      px(55) = ran

      ! Cp ds/dlog p
      px(56) = scp * st/pt

      ! Values of LK and LQ
      px(57) = lk
      px(58) = lq

      ! Thermohaline mixing and convection stuff
      px(31) = avmu
      px(23) = sgth
      px(30) = sg
      px(7) = grad
      px(17) = r

      ! Rotation rate and mixing coefficients
      px(59) = omega
      px(65) = sgf
   end if

   if ( ji <= 0 ) then
      ! values needed for nucleosynthesis calculations in FUNCS2, EQUNS2
      ht(Jstar, 1, jk) = rho
      ht(Jstar, 2, jk) = sg
      ht(Jstar, 3, jk) = zt
      ht(Jstar, 4, jk) = mk
      ht(Jstar, 5, jk) = neo
      ht(Jstar, 6, jk) = nio
      ht(Jstar, 7, jk) = nzz
      ht(Jstar, 8, jk) = avma
      ht(Jstar, 9, jk) = ne
      ht(Jstar, 10:18, jk) = na(1:9)
      ht(Jstar, 19, jk) = mt
      ht(Jstar, 20, jk) = apr
      ht(Jstar, 21, jk) = atr
      ht(Jstar, 22, jk) = sgth * avmu  ! Themohaline mixing coefficient, apart from molecular weight
      ht(Jstar, 23, jk) = Frad
      ht(Jstar, 24, jk) = r
      ht(Jstar, 25, jk) = eos%dv

      ! Always weakly mix inner and outer grid points in nucleosynthesis calculation
      if ( jk < 1.d0+0.01d0*(kh+1) .or. jk > kh-0.05d0*(kh+1) ) ht(Jstar, 2, jk) = min(ht(Jstar, 2, jk), 1.0d-2*crd*con)
   end if


   if ( ji < 0 ) return    ! Only do one star
   end do                  ! Continue with *2, if ktw = 2 (TWIN mode)

end subroutine funcs1


end module structure_functions


!> \brief Solve for dimensionless L1, L2 potentials
pure subroutine potent ( qq, dphi )
   use real_kind
   use constants
   implicit none
   real(double), intent(in) :: qq
   real(double), intent(out) :: dphi
   integer :: ij
   real(double) :: q,q1,q2,q3,q4,q5,q6,xl1,xl2

   q = qq
   if ( q < 1.0d0 ) q = 1.0d0/qq
   q5 = 1.0d0 + q
   xl1 = (q/(3.0d0*q5))**C3RD
   xl2 = 2.0d0 - q*xl1/q5
   q3 = 3.0d0 + q
   q2 = q + q
   q4 = 3.0d0 + q2
   q1 = 2.0d0 + q
   q6 = 1.0d0/q5
   do ij = 1, 4
     ! Newton-Raphson iteration for L1, L2 points
     xl1 = xl1 + &
          (q - xl1*(q2 - xl1*(q - xl1*(q3 - xl1*(q4 - xl1*q5)))))&
          /(q2 - xl1*(q2 - xl1*(3.0d0*q3 - xl1*(4.0d0*q4 - xl1*5.0d0*q5))))
     xl2 = xl2 + &
          (q - xl2*(q2 - xl2*(q1 - xl2*(q3 - xl2*(q4 - xl2*q5)))))&
          /(q2 - xl2*(2.0d0*q1 - xl2*(3.0d0*q3 - xl2*(4.0d0*q4 - xl2*5.0d0*q5))))
   end do
   dphi = q*q6/xl1 + q6/(1.0d0 - xl1) + 0.5d0*(xl1 - q6)**2&
       - (q*q6/xl2 + q6/(xl2 - 1.0d0) + 0.5d0*(xl2 - q6)**2)
end subroutine potent

