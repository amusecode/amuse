!> \file stars_structure.f90   Contains funcs1() and equns1()

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
!!     structure_variables - Stellar structure is stored in the SX array (JI<=0 only)
!!
!! \todo FIXME: because of the presence of an optional argument, this function
!! needs an "explicit interface" (in C terms, a declaration) so that the
!! argument list is known when the function is called and the compiler can
!! do the right thing when the optional argument is not present (the code
!! segfaults otherwise).
!! Ideally we would place funcs1 in a module, which means the compiler will
!! generate the interface for us, but there is a circular dependance of
!! funcs1 on modules defined in output_properties.f90. That needs to be fixed
!! first. For now, we write the interface ourselves (in funcs1_interface.f90)
!! ------------------------------------------------------------------------------
!<
subroutine funcs1 ( jk, ji, var, dvar, fn1, eosout, abundout )
   use real_kind
   use mesh
   use mesh_enc
   use extra_elements
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
   use structure_variables
   use accretion_abundances
   use radiative_acceleration
   
   implicit none
   ! Declaring subroutine variables
   integer, intent(in) :: jk, ji
   real(double), intent(in) :: var(NVAR), dvar(NVAR)
   real(double), intent(out) :: fn1(NFUNC)
   type(eostate), optional, intent(out) :: eosout
   type(abundance), optional, intent(out) :: abundout
   
   ! Declaring common blcok variables explicitly
   
   real(double) :: af, at, vx16, m, vx1, vqk, ar, L, vx4, vx12, vx20, ai, vomega, phi, phis, vx14, vx24, vx28, vx56   ! stellar vars
   real(double) :: oa, ecc, xi, mb                                                                                    !  binary vars
   real(double) :: daf, dat, dx16, dm, dx1, dvqk, dar, dL, dx4, dx12, dx20, dai, domega, dphi, dphis          ! stellar var. changes
   real(double) :: dx14, dx24, dx28, dx56                                                                     ! stellar var. changes
   real(double) :: doa, decc, dxi, dmb                                                                        !  binary var. changes
   
   real(double) :: bcp, bct, vp, vpk, vr, vrk, vt, vtk, vl, lk, lq, mt, vm, vmk, sg, wt1, x1, x1t, x16, x16t, x4, x4t, x12
   real(double) :: x12t, x20, x20t, bcm, vi, vik, vphi, phik, bcf, bcs, bcph, x14, x14t, avmu, sgth, omega, omegat, sgam, si
   real(double) :: bca, bce, xim, xik, dlrk, bcmb, x24, x24t, x28, x28t, x56, x56t
   
   real(double) ::  xvs(nsvar + 1             : nsvar +   nxvstar)                           ! Variables for star (S)
   real(double) ::  xvb(nsvar + 2*nxvstar + 1 : nsvar + 2*nxvstar + nxvbin)                  ! Variables for binary (B)
   real(double) :: dxvs(nsvar + 1             : nsvar +   nxvstar)                           ! Changes in variables for star (S)  
   real(double) :: dxvb(nsvar + 2*nxvstar + 1 : nsvar + 2*nxvstar + nxvbin)                  ! Changes in variables for binary (B)
   

   real(double) :: xfs(nsfunc + 1             : nsfunc +  nxfstar)                           ! Extra output variables for star (S)
   real(double) :: xfb(nsfunc + 2*nxfstar + 1 : nsfunc +2*nxfstar + nxfbin)                  ! Extra output variables for binary (B)

   real(double) :: xh, xhe, xc, xn, xo, xne, xmg, xsi, xfe, na(9), neo, nio, nzz, avma, ne, xai(9, 26)
   
   ! Various things, needed for gravitational settling:
   integer :: nref = 1                    ! Reference abundance (H)
   integer :: i
   real(double) :: gi, apr, atr
   real(double) :: wpi(9), ddiff(9), xa(9)
   real(double) :: wpp, kf_num, kf_denom
   real(double) :: avgz(9), nn(9), nne, diffp(9),difft(9),diffc(9,9)
   real(double) :: logkappa, Frad
 
   real(double) :: ap, arho, u, p, rho, fk, t, sf, st, zt, grada, scp, rf, rt, xhi, s, pr, pg, pf, pt, en
   real(double) :: rpp, r33, r34,  rbe, rbp, rpc, rpna,  rpo,  r3a,  rac
   real(double) :: ran, rao, rane, rcca, rco, roo, rgne, rgmg, rccg, rpng
   real(double) :: ex, enx, wmu, delta, phie, ext, fkt, fkr, prandtl
   type(eostate) :: eos
   type(abundance) :: abund

   ! Declaring local functions
   real(double) :: cbrt, rlobe
   
   ! Declaring external functions
   real(double) :: step, pstv

   ! Declaring local variables
   integer :: Jstar, ikk, j, Jstar_accretor

   real(double) :: avmu_prev, avmu_new, r2, r3, wth, denc, target_daf, target_dat, ds, mor3, r0, lor3,  r2mu, fp
   real(double) :: ft, gor, apmu, gradr, s2, wc1, mu, lom, wc2, wc3, wc4, gta, atmu, vpp, vtt, xmc, vmf, vmmu, vrr, muk, mk, b
   real(double) :: wos, con, fac, f_wl, vml, df0, sgf, deg, ug, w, beta, fac1, fac2, v2, sgth_numer, sgth_denom, dth, sgth2, rec
   real(double) :: ves, vmu, app, mmk, absl, avm, xh_surface, xhe_surface, z_eff, rht, rad, w2, tet, rot, mdtddw, alfac, mdtl
   real(double) :: mdtsw, mdot_jvink, mdot_wr, mdot_cmi, mdot, s_x, s_he, s_c, s_o, f_photon_tiring, fdt, m2, e2, ef, we25, rl
   real(double) :: rgr, etgr, oatgr, rlfac, xi1, ww, gam, we5, we65, we2, we3, wea, web, dmt0, cif, alpha

   real(double) :: gradl, dg_ledoux, dsc, gradmu, tkh, nu, gradtsc, a
   real(double) :: aper, dap

   ! For determining mass transfer rates
   real(double) :: rr, rpr, gmrr, gmrp(2)

   real(double) :: ra2s(2), ais(2), oatmb(2)
   real(double) :: et(2), mp(2), pa(2), br(2), zep(2), w1p(2), w2p(2), w3p(2), w4p(2), spp(2), dsp(2), sup(2), dmp(2), brp(2)
   real(double) :: enp(2), mkp(2), phip(2), phisp(2), pap(2), sdp(2)
   real(double) :: mdotedd(2), ledd(2), fmdotlost(2), phiratio(2)

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

   ! These (definitely xfs, do xfb as well) may be uninitialised and are used in equns1:
   xfs = 0.0_dbl
   xfb = 0.0_dbl




   do Jstar = 1, ktw
   ! Input variables and their changes, star 1 or 2
   af     = var(24*(Jstar - 1) + 1);     daf    = dvar(24*(Jstar - 1) + 1)
   at     = var(24*(Jstar - 1) + 2);     dat    = dvar(24*(Jstar - 1) + 2)
   vx16   = var(24*(Jstar - 1) + 3);     dx16   = dvar(24*(Jstar - 1) + 3)
   m      = var(24*(Jstar - 1) + 4);     dm     = dvar(24*(Jstar - 1) + 4)
   vx1    = var(24*(Jstar - 1) + 5);     dx1    = dvar(24*(Jstar - 1) + 5)
   vqk    = var(24*(Jstar - 1) + 6);     dvqk   = dvar(24*(Jstar - 1) + 6)
   ar     = var(24*(Jstar - 1) + 7);     dar    = dvar(24*(Jstar - 1) + 7)
   L      = var(24*(Jstar - 1) + 8);     dL     = dvar(24*(Jstar - 1) + 8)
   vx4    = var(24*(Jstar - 1) + 9);     dx4    = dvar(24*(Jstar - 1) + 9)
   vx12   = var(24*(Jstar - 1) + 10);    dx12   = dvar(24*(Jstar - 1) + 10)
   vx20   = var(24*(Jstar - 1) + 11);    dx20   = dvar(24*(Jstar - 1) + 11)
   ai     = var(24*(Jstar - 1) + 12);    dai    = dvar(24*(Jstar - 1) + 12)
   vomega = var(24*(Jstar - 1) + 13);    domega = dvar(24*(Jstar - 1) + 13)
   phi    = var(24*(Jstar - 1) + 14);    dphi   = dvar(24*(Jstar - 1) + 14)
   phis   = var(24*(Jstar - 1) + 15);    dphis  = dvar(24*(Jstar - 1) + 15)
   vx14   = var(24*(Jstar - 1) + 16);    dx14   = dvar(24*(Jstar - 1) + 16)
   vx24   = var(NXFstar*(Jstar - 1) + NMg24); dx24   = dvar(NXFstar*(Jstar - 1) + NMg24)
   vx28   = var(NXFstar*(Jstar - 1) + NSi28); dx28   = dvar(NXFstar*(Jstar - 1) + NSi28)
   vx56   = var(NXFstar*(Jstar - 1) + NFe56); dx56   = dvar(NXFstar*(Jstar - 1) + NFe56)

   ! Input variables and their changes, binary orbit
   oa       = var(17); doa      = dvar(17)
   ecc      = var(18); decc     = dvar(18)
   xi       = var(19); dxi      = dvar(19)
   mb       = var(20); dmb      = dvar(20)

   ! Extra variables, star 1 or 2
   i = nsvar+(Jstar-1)*nxvstar
   xvs(nsvar+1:nsvar+nxvstar) = var(i+1:i+nxvstar); dxvs(nsvar+1:nsvar+nxvstar) = dvar(i+1:i+nxvstar)

   ! Extra variables, binary orbit
   i = nsvar+2*nxvstar
   xvb(i+1:i+nxvbin) = var(i+1:i+nxvbin); dxvb(i+1:i+nxvbin) = dvar(i+1:i+nxvbin)

   ! M is the mass variable
   mt = dm/dt                                                  ! 1e33 g/s
   dmp(Jstar) = dm

   ! Compute the rotational period from the rotational frequency
   omega = vomega
   aper = 2.0*cpi/(dabs(omega) * csday)                        ! CGS
   dap = -aper/dabs(omega) * domega

   ! Calculate the mean-molecular weight (for thermohaline mixing)
   ! Use the values from the previous timestep instead of the current one; this
   ! seems to be more numerically stable
   !> \todo FIXME: now that we have the value of the molecular weight gradient (for the
   !! purpose of semiconvection) we should rework the TH mixing prescription.
   !<
   xh =  h(5 + (Jstar-1)*24, jk)
   xhe = h(9 + (Jstar-1)*24, jk)
   xc =  h(10 + (Jstar-1)*24, jk)
   xn =  h(16 + (Jstar-1)*24, jk)
   xo =  h(3 + (Jstar-1)*24, jk)
   xne = h(11 + (Jstar-1)*24, jk)
   xsi = h(NSi28 + (Jstar-1)*NXFSTAR, jk)
   xfe = h(NFe56 + (Jstar-1)*NXFSTAR, jk)
   xmg = max(0.0d0, 1.0d0 - xh - xhe - xc - xn - xo - xne - xsi - xfe)
   avmu_prev = 1.0/(0.5 + 1.5*xh + (42.0*xhe + 14.0*xc + 12.0*xn + 10.5*xo +&
       8.4*xne + 7.0*xmg + 6.0*xsi - 3.0*xfe)/168.0)
   ! set up the composition variables (sanely)
   xh =  max(0.0d0, min(vx1,  1.0d0))
   xhe = max(0.0d0, min(vx4,  1.0d0))
   xc =  max(0.0d0, min(vx12, 1.0d0))
   xn =  max(0.0d0, min(vx14, 1.0d0))
   xo =  max(0.0d0, min(vx16, 1.0d0))
   xne = max(0.0d0, min(vx20, 1.0d0))
   xsi = max(0.0d0, min(vx28, 1.0d0))
   xfe = max(0.0d0, min(vx56, 1.0d0))

   if (use_mg24_eqn) then
     xmg = max(0.0d0, min(vx24, 1.0d0))
   else
     xmg = 1.0d0 - xh - xhe - xc - xn - xo - xne - xsi - xfe
     if ( xmg < 0.0d0 ) xmg = 0.0d0
   end if

   avmu_new = 1.0/(0.5 + 1.5*xh + (42.0*xhe + 14.0*xc + 12.0*xn + 10.5*xo +&
       8.4*xne + 7.0*xmg + 6.0*xsi - 3.0*xfe)/168.0)

   avmu = avmu_smooth*avmu_new + (1.0-avmu_smooth)*avmu_prev

   ! Store composition variables in an array as well, for convenience
   xa(1:9) = (/ xh, xhe, xc, xn, xo, xne, xmg, xsi, xfe /)

   ! Equation of state
   call statel ( jk, ji, af, at, xa, Jstar, abund, eos )
   
   ! Copy EOS output to local variables (for convenience)
   ap    = eos%ap;    arho    = eos%arho; u     = eos%u;     p    = eos%p
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
   delta = eos%delta; phie    = eos%phi;  ext   = eos%ext;   fkt  = eos%fkt
   fkr   = eos%fkr;   prandtl = eos%prandtl

   na    = abund%na;    neo   = abund%neo;   nio = abund%nio; nzz = abund%nzz
   ne    = abund%ne;    xai = abund%xai
   avma  = abund%avm;


   enp(Jstar) = u + p/rho
   r2 = dexp(2.0d0*ar) - ct(8)
   r = dsqrt(dabs(r2))
   r3 = r*abs(r2)
   wth = 2.0d0*kth/(kth**4 + 1.0d0) * luminosity_fudge
   if ( jmod == 0 ) wth = 0.0d0
   denc = 0.0d0
   if (usemenc .and. m<th(4,1)) then
     !LOM = EX + EN + ENX + MENC*MEA - WTH*T*(SF*DAF + ST*DAT)/DT
     target_daf = get_iph(m, 1)-af
     target_dat = get_iph(m, 2)-at
     ds = sf*target_daf + st*target_dat
     denc = impose_entropy_factor*(t*ds/dt)
     lom = ex + en + enx + denc * luminosity_fudge !LOM = dL/dM
   else
     ds = sf*daf + st*dat
     if (interpolate_thermal_energy) ds = (s - get_sm(m))
     lom = ex + en + enx + enc - wth*t*ds/dt
   end if
   lom = lom + enc_parachute

   if ( jk == kh ) then
     ! Approximations for the centre (M=0); the normal expressions become singular
     mor3 = cpi4*rho*c3rd
     lor3 = lom*mor3
     r0 = 1.0d-10
     m3 = cbrt(mor3 * r0**3)
     mp(Jstar) = m3**3
   else
     ! Values away from the centre
     mor3 = abs(m)/r3
     lor3 = L/r3 !note the difference with LOM=dL/dm and not L/M
     m3 = cbrt(abs(m))
     r0 = r
     mp(Jstar) = abs(m)
   end if
   r2mu = 3.0d0*mor3**c3rd/(cpi4*rho)                          ! CGS *and* code units
   
   ! Correction factors FP and FT to the pressure and temperature gradients for
   ! rotating stars. Could (should?) use the more sophisticated scheme from
   ! Endal&Sofia 1976 for this. See code by Alexander Heger. This can be adapted
   ! for doing implicit calculations easily, but will require an extra equation
   ! to solve for R0...
   ! Limit FP to 0.1, smaller values seem to give convergence problems -SdM
   fp = max( 0.1d0, 1.0d0 - 2.0d0*c3rd/( mor3*(cg2*aper)**2 ))
   ft = 1.0d0
   
   ! pressure gradient equation; gravity corrected for rotation.
   gor = cg*mor3 * fp                                          ! CGS
   phim = 5.0d21*gor*r2mu                                      ! Code units
   apmu = - rho*phim/p                                         ! Code units
   hp = dmin1(p/(gor*1.0d11*r0*rho), dsqrt(p/(cg*rho*rho)))    ! CGS
   
   ! temperature gradient equation
   gradr = 0.25d0*fk*p*lor3/(cpi4*cl*gor*pr) * ft
   gradmu = expl_var(jk, explv_gradmu, Jstar)
   gradl = grada + convection_ledoux * phie*gradmu/delta
   dg = gradr - grada
   dg_ledoux = gradr - gradl
   
   ! mixing-length theory; the cubic equation for wl/chi has one real solution
   s2 = grada*scp*t
   wc1 = min(0.5d0*s2*(calp*calp*hp/(9.0d0*xhi))**2, 1.0d302)
   wc2 = min(546.75d0*wc1*dmax1(0.0d0, dg_ledoux) + 73.0d0, 1.0d302)
   wc3 = max(1.0d-200, cbrt(wc2 + dsqrt(wc2*wc2 + 12167.0d0))) ! Overflow guard
   wc4 = dmax1(wc3 - 23.0d0/wc3 - 2.0d0,0.0d0)
   wcv = wc4*xhi/(3.0d0*calp*hp)
   wl = calp*hp*wcv
   
   ! if GRADR <= GRADA, i.e. DG <= 0, we get WCV = 0 and GRAD = GRADR
   grad = gradr - 4.0d0*hp*wcv**3/(calp*s2*xhi)
   
   ! Heat transport due to semi-convection, after Langer, Sugimoto & Fricke (1983)
   if (dg >= 0.0d0 .and. dg_ledoux < 0.0d0 .and. csmce > 0.0d0) then
     b = pg/p
     a = b*(8. - 3.*b) / (32. - 24.*b - 3.*b**2)
     gradtsc = ( (0.5d0*csmc*a*gradmu + dg_ledoux)**2 +&
          2.0d0 * csmc*(dg**2 - dg_ledoux*dg -a*gradmu*dg) )/(csmc-2.0d0)**2
     gradtsc = sqrt(gradtsc)
     gradtsc = -gradtsc + gradr +&
          0.5d0*(2.0d0*dg_ledoux - 2.0*csmc*dg + csmc*a*gradmu)/(csmc-2.0d0)
     grad = gradtsc
   end if
   
   ! Heat transport due to thermohaline convection, based on the same idea
   if (dg < 0.0d0 .and. gradmu < 0 .and. cthe > 0.0d0) then
     b = pg/p
     a = b*(8. - 3.*b) / (32. - 24.*b - 3.*b**2)
     grad = grada + 0.5*dg + 6.*cth*gradmu -&
          0.5*sqrt( (12.*cth*gradmu + dg)**2 - 48.*cth*a*gradmu**2 )
   end if
   gta = grad - grada
   atmu = grad*apmu
   
   ! mesh spacing equation, with modified pressure, temperature gradient equations
   vp = ct(4)*ap + ct(5)*dlog(p + ct(9))
   vpp = ct(4) + ct(5)*p/(p + ct(9))
   vt = ct(7)*dlog(t/(t + ct(10)))
   vtt = ct(7)*ct(10)/(t + ct(10))
   
   ! Extra mesh spacing terms for the AGB
   vp = vp + ct(12)*dlog(p + ct(11)) - ct(14) * dlog(p + ct(13))
   vp = vp + ct(16)*dlog(p + ct(15))
   vpp = vpp + ct(12)*p/(p + ct(11)) - ct(14)*p/(p + ct(13))
   vpp = vpp + ct(16)*p/(p + ct(15))
   !if (ji == 0) then
   !   print *, 'AGB MSF terms:', JK
   !   print *, 'Pres = ', CT(11), CT(13), CT(15)
   !   print *, '       ', P/(P + CT(9))
   !   print *, '       ', P/(P + CT(11)), P/(P + CT(13)), P/(P + CT(15))
   !   print *, ''
   !end if
   
   ! VR and VM must go to zero at centre like r**2, m**(2/3)
   xmc = ct(6)*cbrt(mc(Jstar)*mc(Jstar))
   vmf = xmc + m3*m3
   if ( jk == kh ) vmf = xmc
   vm = dlog(dabs(xmc/vmf))
   vmmu = - 1.0d0/vmf
   vr = - ct(3)*dlog(dabs(r2/ct(8) + 1.0d0))
   vrr = - ct(3)/(r2 + ct(8))
   
   ! QQ is the quantity that the meshpoints should be at equal intervals of
   qq = vp + vt + vm + vr
   qm = vpp*apmu + vtt*atmu + vmmu + vrr*r2mu
   
   ! L dependence for mesh function
   qm = qm + ct(2)*ex*1.5d-8*m3
   muk = vqk/qm                                                  ! Code units
   mk = 1.5d0*m3*muk
   mkp(Jstar) = mk
   
   ! Convert d/dm to d/dk = (d/dm)(dm/dk)
   vmk = vmmu*muk
   vpk = vpp*apmu*muk
   vtk = vtt*atmu*muk
   vrk = vrr*r2mu*muk
   
   ! Calculate ratio of timestep to local radiative diffusion timescale.
   ! This is used to change the differencing scheme in the energy and diffusion
   ! equations. See Sugimoto 1970.
   ! WT1 is the ratio between the radiative energy diffusion time and the
   ! current timestep. If this becomes big, the differencing scheme should
   ! be adjusted. Just what "big" means is set by the input option
   ! OFF_CENTRE_WEIGHT: the larger this number is, the harder it is to
   ! satisfy the criterion.
   
   !> \todo CHECK: the variable wt was called tw in the CB vbles until svn rev.1170;
   !! it was renamed to wt because it had that name in most other routines.  Was the 
   !! name in other routines a mistake, or was it used as a dummy variable there?
   !! The local variable now called wt1 here was originally called wt
   !<
   
   wt1 = 0.0d0
   wt = 0.0d0
   if (gta < 0.0d0 .and. jk < kh) then   ! Radiative
     wt1 = -1.0d22 * r2mu*muk / (xhi*dt) * mor3/(cpi4 * rho)
   end if
   wt1 = wt1/off_centre_weight
   if (wt1 /= 0.0d0) wt = 1.0d0/wt1
   
   ! potential equation, moment of inertia, angular frequency
   vphi = phi
   phip(Jstar) = phi
   phisp(Jstar) = phis
   phik = phim*muk
   vi = ai
   vik = r2*mk/1.5d0                                           ! Code units
   si = vik/mk                                                 ! Code units
   
   ! Total angular momentum, code units
   xfs(fx_am) = xvs(ntam)                                      ! Code units
   xfs(fx_amk) = vik/aper

   ! Composition equations:
   ! hydrogen equation with ZAMS baryon correction when KCN = 1
   x1 = vx1
   x1t = (2.0*kx*((1-kcn)*rpp + rpc + rpng + (1-2*kcn)*(rpna + rpo)) + dx1/dt)*mk
   
   ! helium equation:
   x4 = vx4
   x4t = (4.0*(-kx*(0.5*rpp + rpna + rpo)*(1-kcn)&
       + ky*(3.0*r3a + rac + 1.5*ran + rao + rane)&
       - kz*(rcca + rco + 2.0*roo + rgne + rgmg)) + dx4/dt)*mk
   
   ! carbon equation:
   x12 = vx12
   x12t = (12.0*(kx*(rpc - rpna) - ky*(r3a - rac)&
       + kz*(2.0*(rcca + rccg) + rco)) + dx12/dt)*mk
   
   ! nitrogen equation:
   x14 = vx14
   x14t = (14.0*(kx*(rpna + rpng - rpc - rpo) + ky*ran) + dx14/dt)*mk
   
   ! oxygen equation:
   x16 = vx16
   x16t = (16.0*(kx*(rpo - rpng) - ky*(rac - rao)&
       + kz*(rco + 2.0*roo - rgne)) + dx16/dt)*mk
   
   ! neon equation:
   x20 = vx20
   x20t = (20.0*(ky*(rane - ran - rao) + kz*(rgne - rgmg - rcca))&
       + dx20/dt)*mk
   
   ! Magnesium equation:
   x24 = vx24
   x24t = (24.0*(kz*(rgmg - rane - rccg - rco - roo)) + dx24/dt)*mk

   ! Silicon equation:
   x28 = vx28
   x28t = (dx28/dt)*mk

   ! Iron equation:
   x56 = vx56
   x56t = (dx56/dt)*mk

   ! Artificial composition adjustment:
   if(adj_comp) then
     cif = mixing_fudge*impose_composition_factor
     x1t  = x1t  + cif*(x1  - get_iph(m,  5)) / dt*mk;
     x4t  = x4t  + cif*(x4  - get_iph(m,  9)) / dt*mk;
     x12t = x12t + cif*(x12 - get_iph(m, 10)) / dt*mk;
     x14t = x14t + cif*(x14 - get_iph(m, 16)) / dt*mk;
     x16t = x16t + cif*(x16 - get_iph(m,  3)) / dt*mk;
     x20t = x20t + cif*(x20 - get_iph(m, 11)) / dt*mk;
   end if

   ! Mixing coefficient, energy equation, mass-transfer equation, central BCs.
   b = pr/pg
   if ( jk == kh ) then
     wos = 1.0d10
     con = 0.0d0
     vl = L
     lk = lom*muk*dsqrt(xmc)
     lq = 0.0d0
     fac = 0.0d0
     xik = 0.0d0
     vm = m
     vr = r2
   else
     wos =(2.5d0+b*(2.0d1+1.6d1*b))*(1.5d0*cu*m3/dabs(apmu*m)+1.0d0)
     !CON = 6.0D-22*R2*M3/(R2MU*R2MU*MUK)*WC1**C3RD*XHI
     ! The mixing coefficient is calculated approximately and possibly reduced
     !  by underrelaxation (using the cube root of WCV is convenient because it
     !  has about the right behavior near the edge of the convective zone; physically
     !  it is nonsense of course).
     f_wl = mixing_boost*wc1 + (1.0d0 - mixing_boost)*wcv
     con = 6.0d-22*r2*m3/(r2mu*r2mu*muk)*cbrt(f_wl)*xhi
     vml = dsqrt(vmf)/m3
     vl = L*vml
     lk = vml*mk*(lom - c3rd*xmc/vmf*lor3/mor3)
     lq = vml*wth*scp*t*apmu*muk*gta
     df0 = cdf*1.0d22*cg*m/r
     ! FAC = SMF(1, PHI, DF0)
     fac = step(phi, df0)
     ! XIK = CMT*DSQRT(SMF(2, 2.0D0*PHIS, DF0))/R*FAC*MK
     xik = cmt*sqrt(pstv(2.0d0*phis, df0))/r*fac*mk
   end if
   !   If the thermal energy term is interpolated there is no advection term
   if (interpolate_thermal_energy) lq = 0.0d0
   if (usemenc .and. m<th(4,1)) lq = 0.0d0
   
   
   ! Conversion factor for diffusion coefficients (converting to code units)
   sgf = (cpi4*r2*rho)**2/mk
   
   ! Fudged convective diffusion coefficient: COS > 0 gives overshooting
   !> \todo FIXME: when Ledoux convection is used, currently does overshooting from
   !! the Ledoux boundary, which doesn't make muxh sense - it should do
   !! overshooting from the Schwarzschild boundary.
   !<
   if ( xh >= 1.0d-7 ) deg = cos/wos
   if ( xh < 1.0d-7 ) deg = cps/wos
   eg = dg_ledoux + deg
   ug = pstv(eg, 0.0d0)       ! used to be  UG = PS(EG)

   select case (convection_scheme)
   case (1)       ! Eggleton 1972
     ! This is where the convective mixing coefficient is fudged by
     !  taking the square rather than the cube root, as per the
     !  mixing-length model.
     sg = crd*con*ug*ug

   case (2)       ! Pols&Tout 2001, Stancliffe, Tout & Pols 2004
     w = dg / gradr
     ! Fudge factor beta - beta = inf means no fudge.
     ! Use different values of beta in the exhausted core and in the
     ! envelope; should be set through init.dat, but hardcoded for
     ! now.
     beta = 1.0d0                     ! Envelope
     if (xh < 1.0d-7) beta = 5.0d-5   ! Core; Stancliffe&al 2004
     sg = con*cbrt(ug) * beta*w / ( 1.0d0 - (1.0d0-beta)*w )
   case default   ! Die
     sg = 0
   end select

   ! Semi-convection, after Langer, Sugimoto & Fricke 1983 (A&A)
   ! CSMC = 0.04 is the effciency of semiconvection, alpha in (10) of LSF83.
   ! Stability condition: DG > 0 (Schwarzschild unstable), DG_LEDOUX<0
   ! (Ledoux unstable)
   dsc = 0.0d0
   if (dg >= 0.0d0 .and. dg_ledoux < 0.0d0) then
     dsc = -csmc/6.0d0 * xhi * dg / dg_ledoux
   end if
   ! Convert to code units and add to the overall diffusion coefficient
   ! SGF in 1e11 g/cm**2, DSC in g/cm**2, so SGF*DSC is in 1e-22 [1e33 g/s]
   sg = sg + 1.0d-22*dsc*sgf

   ! Artificial extra mixing: ARTMIX is the diffusion coefficient in cm^2/s
   ! 0 by default, but can be set from init.dat.
   sg = sg + artmix*sgf

   ! For ZAMS runs (KCN = 1) always weakly mix inner and outer meshpoints
   if (kcn == 1) then
     if ( jk < 0.075d0*kh.or.jk > 0.8d0*kh ) sg = min(sg, 1.0d-4*crd*con)
   end if
   !IF ( JK < 0.075D0*KH ) SG = MIN(SG, 1.0D-2*CRD*CON)

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
   
   sg = sg * mixing_fudge

   wl = sg
   xim = xi
   dlrk = 0.0d0
   fac1 = 0.0d0
   fac2 = 0.0d0
   v2 = 0.0d0
   ! Mixing coefficient for thermohaline mixing (except for mu gradient)
   if (cth > 0.0d0) then
     ! radius is ok ~ 10^10 cm
     ! hp is >~ radius (seems ok)
     ! ca is ok (7.566e-15 in cgs)
     ! cl is ok (3e10 in cgs)
     ! T ~ 10^7 K near the centre
     ! gta ~ 10^-3 (small)
     ! scp ~ 10^8
     ! fk ~ 1
     ! rho ~ 10-100
     ! avmu ~ 0.588
     ! mk ~ -1.6e-3
     sgth_numer = 16.0 * cpi * r2 * hp * ca * cl * t**3
     sgth_denom = -gta * scp * fk * rho * avmu * abs(mk)
     ! DTH is in 1e-11 cm^2/s
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
   else
     sgth = 0.0d0
   end if

   ! Angular momentum transport, specific angular momentum
   !  ( (R2+CT(8))*DAR - C3RD*R2MU/M3*DM )/R0**2 can be CGS or code units
   !  Final units: 1/s^2 * [SI]*[M] = 10^55 g cm^2/s^2
   omegat = ( domega + 2*omega/r0**2*&
       ( (r2+ct(8))*dar - c3rd*r2mu/m3*dm ) )/dt * si * mk
   ! Innermost meshpoint is always at constant M&R (or should be!)
   if (jk == kh) then
     omegat = domega/dt * si * mk
   end if

   ! The transport coefficient for angular momentum transport is the sum of
   !  different coefficients, but most of these depend on the local shear,
   !  which has to be found from equns1.
   ! In terms of units, the dimension of the mixing coefficients should be
   !  in code units.
   sgam = con*ug**2 * mixing_fudge

   ! When relaxing the AM diffusion term, always have some (weak!) mixing,
   !  this helps convergence.
   sgam = sgam - (1.0-amadv_smooth)*1.0d-16
   !SGAM = SGAM - 1.0D-16

   ! Enforce rigid rotation on the first model
   !IF (JMOD==1) SGAM = -1.0D21
   if (jmod==1) sgam = con

   ! Always make outer few meshpoints rotate as a solid body
   ! This is designed to aid numerical stability
   !IF ( JK < 0.075D0*KH ) SGAM = CON
   !IF ( JK > 0.95D0*KH ) SGAM = CON

   ! Dynamical shear instability. For this, we need both the Richardson number
   !  and the size of the unstable region. For the former we need the gradient
   !  of omega, for the latter we need the size of the unstable region, here
   !  estimated to be one pressure scaleheight (compare notes from Yoon).
   ! We omit the composition gradient from the Richardson number for simplicity.

   ! All quantities in the Richardson number should be in CGS units.
   !  RHO/P * PSTV(-GTA, 0.0D0)*GOR**2 is in CGS
   !  R2MU is in CGS (or code units)
   !  MUK/R0 is in 10^22 g^(2/3)/10^11 cm = 10^11 g^(2/3)/cm
   !  So R2MU*MUK/R0 is in 10^11 cm
   xfs(fx_rich) = rho/p * pstv(-gta, 0.0d0)*gor**2 * (1.0d11*r2mu*muk/r0)**2

   ! The mixing coefficient should be in [1e33g]/s.
   !  HP**2/t_dyn is in cm**2/s
   !  SGF is in 1e11 g/cm**2
   ! So SGF*HP**2/t_dyn is in 1e11 g/s = 1e-22 [1e33g]/s
   xfs(fx_ddsi) = amadv_smooth*cdsi * 1.0d-22 * hp**2/sqrt(cg*mor3) * sgf * mixing_fudge

   ! Secular shear instability, apart from domega/dk
   xfs(fx_dssi) = 0.0d0
   if (gta<0.0d0) then
     rec = 2500.0d0    ! Critical Reynolds number
     xfs(fx_dssi) = cssi*c3rd * xhi * 1.0d-33*r**3 * hp / ( gor*pstv(-gta, 0.0d0)*(muk*r2mu)**2 )
     xfs(fx_dssi) = xfs(fx_dssi) * sgf * mixing_fudge
     xfs(fx_sssi) = prandtl*rec * 1.0d-18*rho/p * pstv(-gta, 0.0d0) * (1.0d11*gor*r2mu*muk/r0)**2/8.0d0
   end if

   ! Eddington-Sweet circulation; size of the unstable region should be
   !  limited to a velocity scale-height; the velocity itself should be
   !  reduced by an adverse mu gradient.
   ! To estimate the size of the unstable region, we use the pressure scale-height
   !  instead.
   ves = 0.0d0
   vmu = 0.0d0
   if (gta<0.0d0) then
     ! Circulation velocity VES is in 1e-11 cm/s
     ! We use the Zahn (1992) expression (3.40), which is equivalent to
     ! Heger (2000) except in the last term (after some thermodynamics).
     ! Following Maeder&Zahn (1998) we include the gravitational
     ! term in the energy generation rate. This means we can use LOM
     ! instead of ENUC+EN+ENX. NB: LOM != L/M.
     ! This expression now differs from Maeder&Zahn in not including
     ! the stabilising effect of the mu gradient (which is
     ! incorporated here by using "mu currents", wrong according to
     ! Maeder&Zahn) and different factors that arise due to
     ! differential rotation.
     ves = -grada/gta * lor3 * omega**2/(cg*mor3 * gor) * &
          (2.0*lom/lor3 - 2.0/mor3 + omega**2/(cg*mor3*cpi*rho))/r0**2
     ! Stabilising current driven by molecular weight gradient.
     ! Based on Kippenhahn & al. (1980).
     !> \todo FIXME: add the "phi/delta grad_mu" term to -GTA in both VES and VMU
     ! First, calculate the local thermal timescale. We need the
     ! typical size of a fluid blob, expression (29) of Kippenhahn & al.
     ! We omit a factor Hp because this gets divided out in the
     ! expression for v_mu anyway.
     ! NU is the dimensionless plasma viscosity. All quantities here
     ! are in physical units, giving VMU in cm/s = 1e11 (1e-11 cm)/s.
     nu = xhi*prandtl * 2.5d0 *cl*fk*rho**2 / pr
     nu = max(nu,0.0_dbl)   ! FIXME?  nu can become negative!  This suppresses the symptoms, but doesn't fix the problem!
     tkh = 24*sqrt(-0.6_dbl*nu*grada/gta) * (b/(1.0_dbl+b)) * &
          scp / (16*cl*ca*t**3)
     vmu = 1.0d11*phi*gradmu/(delta * gta * tkh)

     !         ! Get this from D_{th} = d * v_{mu} ~ r * v_{mu}, see Heger PhD thesis
     !         ! DTH is now in 1e22 cm^2/s, R0 is in 1e11 cm, so DTH/R0 is in
     !         !  1e11 cm^2/s = 1e22 (1e-11 cm^2/s)
     !         VMU = 1.0D22*DTH/R0
   end if
   ! Diffusion coefficient, according to Heger (originally due to Mestel)
   ! Export variables to EQUNS1, to complete calculation there.
   ! VES * HP * SGF should be in 1e33 g/s
   !  VES (and VMU) should be in 1e-11 cm/s
   !  HP is in cm
   !  SGF is in 1e11 g/cm**2
   ! -> 1e-33 * VES*HP*SGF is in 1e33 g/s
   xfs(fx_des) = 1.0d-33*cesc*max(0.0d0, dabs(ves) - dabs(vmu))*hp*sgf
   !XFS(FX_VES) = AMADV_SMOOTH * VES * MIXING_FUDGE
   !XFS(FX_VMU) = VMU
   !XFS(FX_HP) = HP
   !XFS(FX_SGF) = SGF
   ! Diffusion coefficient, according to Maeder & Zahn (1998)
   ! Since this uses R0 rather than HP, the conversion factor is 1e-22
   ! This doesn't include the stabilising effect of the mu gradient, which
   ! needs to be added in equns1
   !XFS(FX_DES) = AMADV_SMOOTH*CESC * VES*HP*1.0D-11 * MIXING_FUDGE * SGF
   xfs(fx_des) = amadv_smooth*mixing_fudge * cesc * 1d-22*dabs(ves*r0) * sgf
   xfs(fx_dgmu) = 1.0/(gta * delta * muk * apmu * avmu)

   ! Gravitational settling
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
           xfs(fx_fh+j-1) = mixing_fudge * cgrs*wpi(j) * rho * xa(j) * cpi4*r2
        end do
     else
        xfs(fx_fh:fx_ffe) = 0.0d0
     end if
   else
     ! Use the treatment of Thoul, Bacall & Loeb, 1994 ApJ 421 to solve the
     ! full set of Burgers equations.
     ! See also Hu & al., 2010 A&A.
     if (cgrs > 0.0_dbl) then
        apr = 2.0d0 * r * 1.0d-11 * apmu/r2mu     ! dlogP/dr, cgs
        atr = 2.0d0 * r * 1.0d-11 * atmu/r2mu     ! dlogT/dr, cgs
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
        ! Find terms for the diffusion equations, 9 in total
        call diffusion_terms_burgers(rho,t,9,nn,nne, avgz,can, diffp,difft,diffc)
        ! Compute diffusion fluxes and export to equns1
        ! The factor P/Pg corrects for radiation pressure.
        ! TODO: concentration gradient
        wpi(1:9) = 1.0d-11 * (diffp(1:9)*apr + difft(1:9)*atr) * p/pg
        if (jmod > 2) then
           forall (j=1:9) xfs(fx_fh+j-1) = cgrs*mixing_fudge * rho*xa(j)*cpi4*r2 * wpi(j)
        end if
     else
        xfs(fx_fh:fx_ffe) = 0.0d0
     end if
   end if

   ! Radiative levitation
   ! Will only be used if gravitational settling is enabled as well (because we need the diffusion coefficients and it doesn't
   ! really make much sense to switch it on otherwise anyway). Currently relies on the diffusion coefficients from the Burgers
   ! formalism, so we depend on those too.
   Frad = 1e11*L * gradr/grad / (cpi4 * r2)        ! Radiative flux
   if (CRLEV > 0.0d0 .and. cgrs > 0.0_dbl .and. grs_burgers .and. sg == 0.0d0) then
      nn(1:9) = na(1:9) / nio                      ! Number fractions
      if (ji == rlev_update) then                  ! Only recalculate for "current" values; effectively uses semi-implicit values
         call get_radiative_accelerations(at/cln, arho/cln, dot_product(nn(1:9), can(1:9)), Frad, 9, nn(1:9), kzn(1:9),           &
                                                                                                   radacc(1:9,jk, Jstar), logkappa)
         radacc(1:9,jk,Jstar) = radacc(1:9,jk,Jstar) / can(1:9)   ! Divide by atomic mass to get acceleration
      end if
      wpi(1:9) = 1.0d-11 * diffp(1:9) * radacc(1:9,jk,Jstar) * can(1:9) / (CR*T)
      ! There should be no net mass flux; make sure the net flow velocity is 0 (it is in practice, this just compensates for
      ! numerical noise)
      wpp = -dot_product(wpi(1:9), xa(1:9))
      wpi(1:9) = wpi(1:9) + wpp
      if (jmod > 2) forall (j = 1:9) xfs(fx_fh+j-1) = xfs(fx_fh+j-1) + crlev * mixing_fudge * rho*xa(j)*cpi4*r2 * wpi(j)
   end if

   if ( Jstar /= 1 ) then
     ! Put here the *internal* quantities that depend on *TWIN* variables Quantities
     ! XX needed from both stars should have been stored earlier as XXP(JSTAR) = XX.
     ! Approx. to binary separation, and potential difference between L1, L2
     app = oa*oa*(mp(1) + mp(2))/(mp(1)*mp(2)*cg1)**2
     call potent ( mp(1)/mp(2), dphi )
     df0 = cdf*1.0d22*cg*(mp(1) + mp(2))*dphi/app
     ! Differential mass flux XI, semidetached or in contact: XIK is dXI/dK
     ! Constants CMT, CLT are read in as data (fort.22). SMF is a smoothing
     ! function to prevent a jump discontinuity at L1 surface.
     fac1 = step(phip(1), df0)
     fac2 = step(phip(2), df0)
     mmk = 0.5d0*(mkp(1) + mkp(2))
     v2 = (pstv(phisp(1), df0) - pstv(phisp(2), df0))*(fac1 + fac2)
     ! Mass transfer: if XIK>0, transfer from *2->*1, otherwise from *1->*2
     xik = cmt*v2/(dsqrt(dabs(v2) + 1.0d11)*app)*mmk * mdot_smooth
     ! Heat transfer due to differential rotation, CLT rad/sec
     dlrk = - clt*(enp(1) - enp(2))*fac1*fac2*mmk
   end if

   ! ------------------------------------------------------------------------------
   ! Store computed values if called from PRINTB or REMESH (JI<0) or when
   ! simply calling function values for the current timestep (JI=0)
   ! ------------------------------------------------------------------------------
   if ( ji > 0 ) goto 4
   ikk = kh + 2 - jk
   if ( Jstar == ktw ) then
     ! sundry numbers for PRINTB, FUNCS2; SX(40 to 45,JK) for TWIN variables
     sx(45, ikk) = fac1
     sx(44, ikk) = fac2
     sx(43, ikk) = v2
     sx(42, ikk) = xik
     sx(41, ikk) = enp(1) - enp(2)
     sx(40, ikk) = - dlrk
   end if
   
   ! do. for non-TWIN variables
   eth = - t*(sf*daf + st*dat)/dt + scp*t*apmu*gta*mt/(1.5d0*m3)
   egr = 1.0d22*gor*r2 - u

   !> \todo FIXME: this does NOT work properly in TWIN mode. May not be a problem because this is only used by printb...
   ! Reaction rates
   sx(50, ikk) = rpp
   sx(51, ikk) = rpc
   sx(52, ikk) = rpng
   sx(53, ikk) = rpna
   sx(54, ikk) = rpo
   sx(55, ikk) = ran

   ! Cp ds/dlog p
   sx(56, ikk) = scp * st/pt

   ! Values of LK and LQ
   sx(57, ikk) = lk
   sx(58, ikk) = lq

   ! Thermohaline mixing and convection stuff
   sx(31, ikk) = avmu
   sx(23, ikk) = sgth
   sx(30, ikk) = sg
   sx(7, ikk) = grad

   ! Rotation rate and mixing coefficients
   sx(59, ikk) = omega
   sx(60, ikk) = xfs(fx_rich)
   sx(61, ikk) = xfs(fx_ddsi)
   sx(62, ikk) = xfs(fx_dssi)
   sx(63, ikk) = xfs(fx_ves)
   sx(64, ikk) = xfs(fx_vmu)
   sx(65, ikk) = sgf

   sx(66, ikk) = xfs(fx_sssi) ! needed in printb SdM

   sx(76, ikk) = xfs(fx_des)
   sx(78, ikk) = xfs(fx_dgmu)
   
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

4  continue
   if ( jk == 1 ) then
   ! pressure, temperature surface boundary conditions
   absl = abs(L)
   bcp = log(fk/(r*1.0d11*gor)*(1.5d0*pg + 0.75d0*pr))
   bct = log(1.0d11*absl/(0.75d0*cpi4*cl*r2*pr))

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
   avm = can(1)*x1/cbn(1)  +&!independent
         can(2)*x4/cbn(2)  +&
         can(3)*x12/cbn(3) +&
         can(4)*x14/cbn(4) +&
         can(5)*x16/cbn(5) +&
         can(6)*x20/cbn(6) +&
         can(7)*xmg/cbn(7) +&!sometimes independent
         can(8)*xsi/cbn(8) +&!sometimes independent
         can(9)*xfe/cbn(9)   !sometimes independent
   ! Determine effective surface metallicity i.o.t. scale Mdot with Z_EFF !SdM
   xh_surface = x1/cbn(1)*can(1)/avm   ! massfraction H surface
   xhe_surface = x4/cbn(2)*can(2)/avm  ! massfraction He surface
   z_eff = 1.d0 - xh_surface - xhe_surface
   z_eff = min(1.d0, max(z_eff, 0.d0))
   ! mass loss rate for dynamo-driven wind (MDTDDW); Alfven radius squared (RA2)
   rht = (0.755d0*absl**0.47d0 + 0.05d0*L**0.8d0)/m**0.31d0
   rad = r/rht
   w2 = r*absl/m
   tet = 120.0d0*r/w2**c3rd*rad**2.7d0
   ro = 1.67d0*dabs(aper/tet)
   rot = 1.0d0/(ro**2 + 1.0d0)
   mdtddw = 1.54d-17*w2*rad**2*rot**3.67d0
   bp = 5.4d-3*dsqrt(mor3)*w2**c3rd*rad**3.4d0*rot**1.21d0
   if (bp<1.0d-99) bp = 0.0d0
   alfac = 6.2d-08*(r2/mor3*(bp*bp/mdtddw)**2)**c3rd
   ra2 = 0.0d0
   if (chl>0.0d0) ra2 = chl*(r*alfac)**2

   ! Eggleton's Reimers-like wind
   ! cool superwind (Mdot prop. to Lum/(Env. binding energy))
   mdtsw = 0.0d0
   if ( (be(Jstar) /= 0.0d0) .and. (tn(Jstar) /= 0.0d0) ) &
       mdtsw = m*dmin1(1.3d-5*L/dabs(be(Jstar)), 1.0d1/(tn(Jstar)*csy))

   ! Use smart mass loss routine that determines the mass loss rate based
   ! on recipes appropriate for the stellar parameters, falling back to the
   ! de Jager rate when no other applicable rate can be found.
   mdot = 0.0
   if (smart_mass_loss>0.0) then
     s_x = can(1)*x1/cbn(1)
     s_he = can(2)*x4/cbn(2)
     s_c = can(3)*x12/cbn(3)
     s_o = can(5)*x16/cbn(5)
     mdot = calculate_mdot( &
          !t, L/clsn, L/LoLedd/clsn, m/cmsn, r/crsn, s_x, s_he, s_o, s_c, czs &
          t, L/clsn, m/cmsn, r/crsn, s_x, s_he, s_o, s_c, czs &
          )
     mdot = mdot * cmsn/csy * smart_mass_loss * mdot_smooth
   else
     ! Mass-loss rate for luminous stars from de Jager et al (1988)
     mdtl = 0.d0
     if (cmj > 0.d0) then
        mdtl = calc_mdot_dejager( dlog10(T), dlog10(L/CLSN) )
        if (mdot_errno /= 0) mdtl = 0.0d0
     end if
     ! Mass loss rate for massive stars, Vink et al. (1999, 2000, 2001)
     mdot_jvink = 0.0d0
     if (cmv > 0.0d0) mdot_jvink =&
          cmv*calc_mdot_vink(m/cmsn, L/clsn, t, clogz)
     ! Mass loss rate for Wolf-Rayet stars, Nugis&Lamers 2000
     mdot_wr = 0.0d0
     if (cmnl > 0.0d0) then
        s_x = can(1)*x1/cbn(1)
        s_he = can(2)*x4/cbn(2)
        s_c = can(3)*x12/cbn(3)
        s_o = can(5)*x16/cbn(5)
        mdot_wr = cmnl*&
             calc_mdot_wr_nl(log10(L/clsn),log10(t),s_x,s_he,s_c,s_o,clogz)&
             * mdot_smooth
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
   zep(Jstar) = cml*mdtddw + dmax1(mdot, cmr*mdtsw)
   if ( mutate .or. jmod == 0 ) zep(Jstar) = 0.0d0
   zep(Jstar) = zep(Jstar) * mdot_smooth

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
   if ( jmod /= 0 ) then
     select case (cmi_mode)
     case (1) ! Exponential
        mdot_cmi = mdot_smooth*cmi*mp(1)
     case (2) ! Linear
        mdot_cmi = mdot_smooth*cmi*cmsn
     end select
     if ((m+mdot_cmi*dt) > uc(13)*cmsn) mdot_cmi = (m-uc(13)*cmsn)/dt
   end if

   ! Adjust luminosity for the kinitic energy of the wind (photon tiring)
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
   ef = 1.0d0/dsqrt(1.0d0 - e2)
   we25 = ef**5
   sep = (oa*ef/(cg1*mu))**2/mb
   bper = dsqrt(sep**3/mb)/cg2
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
     ! Switching this of CSU=0 can give strange effects in wider
     ! binaries: accreting mass but conserving spin leads to slower
     ! rotation, opposite to what is realistic. Only current solution is
     ! setting the efficency of spin orbit coupling very high to prevent
     ! this. Switching it on also has unwanted (but physical) effects: the accreting
     ! star will reach critical rotation, especially in wider system where the tides are ineffective
     sup(Jstar) = csu*cg1*dsqrt(m*r_spinup)
     ! SDP: Spin Down: specific angular momentum lost from the donor star.
     ! The specific angular momentum lost due to RLOF is that of a disk with the
     ! radius of the star, *not* the specific angular momentum of a sphere with the
     ! radius of the star.
     sdp(Jstar) = csd*(cg1*r**2)/(cg2*aper)
   end if

   ! Radius for determination of RLOF: average the star's current radius and
   ! its radius on the previous timestep to smooth out kinks in the mass
   ! transfer rate.
   rpr = dsqrt(dexp(2.0d0*h(7, 1)) - ct(8))
   rr = 0.5d0 * (r + rpr)
   
   ! effective gravity (centrifugal and gravity) at the surface - efective
   ! gravity at the Roche-lobe surface (assuming synchronous rotation of
   ! stars with the rotation of orbit)
   gmr  = 1.0d22*cg*(mb*c3rd/sep**3*(r2 - rl*rl) + m*(rl - r)/(r*rl))
   gmrr = 1.0d22*cg*(mb*c3rd/sep**3*(rr*rr - rl*rl) + m*(rl - rr)/(rr*rl))
   
   
   ! GR terms for circularisation and ang. mom. loss
   rgr = cgrt*mu/(mb*bper)*(sep/bper)**5*we25
   etgr = rgr/9.6d1*(3.04d2 + 1.21d2*e2)*ecc
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
   ! Ratio of difference in potential with RLOF surface and potential at surface
   gmrp(Jstar) = gmr
   phiratio(Jstar) = 1.0d-22*gmr / (cg*mor3*r2 + cg*mb*r2/sep**3)
   ! The cheaper model for RLOF, used if CMS > 0
   rlfac = ar - dlog(rl)
   if ( cms > 0.0d0 ) xi = cms*(pstv(rlfac, 0.0d0))**3
   if ( cms > 0.0d0 .and. Jstar == 1 ) xi1 = xi
   if ( cms > 0.0d0 .and. Jstar == 2 ) xi = xi1
   mtr = -xi*csy/cmsn  ! Transport to printb
   ! equilibrium-tide model for tidal friction
   ww = 0.4d0*r2*m/(dabs(ai) + 1.0d-10)
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
     ! dlog(Iw) = dI/I - dP/P, angular momentum change for solid body rotation.
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
     dsp(Jstar) = dxvs(ntam) / ( xvs(ntam)+1.0d-16)    ! For differential rotation
     spp(Jstar) = cg1 * dxvs(ntam) / cg2
   end if
   et(Jstar) = ecc*tfr*(web - bper/aper*wea )
   if ( ji == 0 ) then
     rlf(Jstar) = rlfac
     hspn(Jstar) = spp(Jstar)
   end if
   xfs(fx_macc) = 0.0d0

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
        bcm = dmp(1)/dt + zet(1) + xit(1) - mdot_cmi
        ! The accreted mass flux for the star, used to accrete
        ! material of a given composition.
        ! CCAC is a switch that determines whether this composition
        ! change is taken into account or neglected.
        if (jmod>0) xfs(fx_macc) = ccac*max(mdot_cmi - zet(1) - xit(1), 0.0d0)
        bcs = spp(1)*(dsp(1)/dt + w3p(1) + w4p(1)*zet(1)) - w2p(1)*oa + oatmb(1)
        ! For differential rotation
        !BCS = ( 0.0*DXVS(NTAM) / DT - SI/APER**2*DAP/DT + SI/APER*ZET(1) )
        bca = doa/dt + (w1p(1)*zet(1) + w2p(1))*oa - w3p(1)*spp(1)&
             + oatgr + m/(om*mb)*zet(2)*oa
        bce = decc/dt + et(1) + etgr
        bcmb = dmb/dt + zet(1) + zet(2)
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
     if(compute_twinbc_theoldway) then
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
           xit(j) = pstv((3 - 2*j)*xi, 0.0d0) + pa(3 - j)*zep(j)
        end do
        ! Optionally limit the mass transfer rate
        do j = 1, 2
           xit(j) = sign(min(dabs(xit(j)),mtr_limit*cmsn/csy), xit(j))
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
           !               if (ji == 0) print '(1X,1P,I3,3E24.16)', J,
           !     &            XIT(3-J)*CSY/CMSN, MDOTEDD(J)*CSY/CMSN,
           !     &            (LEDD(3-J)/(1.0d22*PHISP(J)*DABS(XIT(3-J)))+1.0d0)*PHIRATIO(J)
           if (xit(3-j) >= mdotedd(j) .and. cmtel > 0.0d0) then
              !        Cannot accrete all, lose a bit
              !                  FMDOTLOST(J) = 1.0d0 - MDOTEDD(J)/XIT(3-J)
              !        This is expression (7) in Beer & al, but with factors rearranged and
              !        using the potential at the Roche lobe surface in lieu of the potential
              !        in L1 (these are similar and ought to be the same anyway).
              fmdotlost(j) = (1.0d0-ledd(j)/(dabs(gmrp(j)*xit(3-j))))
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

        !: Calculate fraction of material ejected from the system
        do j = 1, 2
           br(j) = brp(j)*step(xit(3 - j) - xit(j), 0.0d0)
        end do
        !: Mass ejected from the system: the fraction of the wind that is not accreted
        !  by the companion and the contribution from non-conservative Roche lobe
        !  overflow
        do j = 1, 2
           zet(j) = (1.0d0 - pa(3 - j))*zep(j) + br(j)*xit(3 - j)
        end do
        ! The net mass accretion flux for each of the two components
        if (jmod>0) then
           fn1(fx_macc)   = ccac*max(mdot_cmi - zet(1) + xit(2) - xit(1), 0.0d0)
           xfs(fx_macc) = ccac*max(mdot_cmi - zet(2) - xit(2) + xit(1), 0.0d0)
        end if
        ! Mass boundary condition for each of the two components...
        fn1(27) = dmp(1)/dt + zet(1) + xit(1) - xit(2) ! *1
        bcm   = dmp(2)/dt + zet(2) - xit(1) + xit(2)   ! *2
        ! ... and for the binary
        bcmb  = dmb/dt + zet(1) + zet(2) ! binary
        ! AM boundary condition for each of the two components ...
        fn1(33) = spp(1)*(dsp(1)/dt + w3p(1) + w4p(1)*zet(1)) - w2p(1)*oa - sup(1)*xit(2) + sdp(1)*xit(1) + oatmb(1) ! *1
        bcs   = spp(2)*(dsp(2)/dt + w3p(2) + w4p(2)*zet(2))-w2p(2)*oa - sup(2)*xit(1) + sdp(2)*xit(2) + oatmb(2)     ! *2
        ! ... and for the binary
        bca   = doa/dt + (w1p(1)*zet(1) + w2p(1))*oa - w3p(1)*spp(1) + (w1p(2)*zet(2) + w2p(2))*oa - w3p(2)*spp(2) + oatgr &
             + (sup(1) - sdp(2))*xit(2) + (sup(2) - sdp(1))*xit(1) ! orbit
        !Eccentricity boudary condition
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

     else !New twin boundary equations (needed for spinup)

        ! Reformulation of boundary conditions will be implemented here, .....-SdM

     end if

   end if      ! if ( Jstar == 2 )
   end if      ! End of surface boundary condition block:  if ( jk == 1 ), ~500 lines up - should this be a subroutine?
   ! Copy output to fn1 in INF

   ! Stellar terms
   fn1(58*(Jstar-1) + 1) = bcp
   fn1(58*(Jstar-1) + 2) = bct
   fn1(58*(Jstar-1) + 3) = vp
   fn1(58*(Jstar-1) + 4) = vpk
   fn1(58*(Jstar-1) + 5) = vr
   fn1(58*(Jstar-1) + 6) = vrk
   fn1(58*(Jstar-1) + 7) = vt
   fn1(58*(Jstar-1) + 8) = vtk
   fn1(58*(Jstar-1) + 9) = vl
   fn1(58*(Jstar-1) + 10) = lk
   fn1(58*(Jstar-1) + 11) = lq
   fn1(58*(Jstar-1) + 12) = mt
   fn1(58*(Jstar-1) + 13) = vm
   fn1(58*(Jstar-1) + 14) = vmk
   fn1(58*(Jstar-1) + 15) = sg
   fn1(58*(Jstar-1) + 16) = wt1
   fn1(58*(Jstar-1) + 17) = x1
   fn1(58*(Jstar-1) + 18) = x1t
   fn1(58*(Jstar-1) + 19) = x16
   fn1(58*(Jstar-1) + 20) = x16t
   fn1(58*(Jstar-1) + 21) = x4
   fn1(58*(Jstar-1) + 22) = x4t
   fn1(58*(Jstar-1) + 23) = x12
   fn1(58*(Jstar-1) + 24) = x12t
   fn1(58*(Jstar-1) + 25) = x20
   fn1(58*(Jstar-1) + 26) = x20t
   fn1(58*(Jstar-1) + 27) = bcm
   fn1(58*(Jstar-1) + 28) = vi
   fn1(58*(Jstar-1) + 29) = vik
   fn1(58*(Jstar-1) + 30) = vphi
   fn1(58*(Jstar-1) + 31) = phik
   fn1(58*(Jstar-1) + 32) = bcf
   fn1(58*(Jstar-1) + 33) = bcs
   fn1(58*(Jstar-1) + 34) = bcph
   fn1(58*(Jstar-1) + 35) = x14
   fn1(58*(Jstar-1) + 36) = x14t
   fn1(58*(Jstar-1) + 37) = avmu
   fn1(58*(Jstar-1) + 38) = sgth
   fn1(58*(Jstar-1) + 39) = omega
   fn1(58*(Jstar-1) + 40) = omegat
   fn1(58*(Jstar-1) + 41) = sgam
   fn1(58*(Jstar-1) + 42) = si

   !> \todo FIXME: magnesium values are stored in the same location for both stars
   !! (in the binary part), this needs to be fixed!
   !<
   fn1(49) = x24
   fn1(50) = x24t
   xfs(FX_X24) = x24
   xfs(FX_X24t) = x24t
   xfs(FX_X28) = x28
   xfs(FX_X28t) = x28t
   xfs(FX_X56) = x56
   xfs(FX_X56t) = x56t

   ! Binary functions
   fn1(43) = bca
   fn1(44) = bce
   fn1(45) = xim
   fn1(46) = xik
   fn1(47) = dlrk
   fn1(48) = bcmb

   ! Extra functions
   i = nsfunc+(Jstar-1)*nxfstar
   fn1(i+1:i+nxfstar) = xfs(nsfunc+1:nsfunc+nxfstar)
   i = nsfunc+2*nxfstar
   fn1(i+1:i+nxfbin) = xfb(i+1:i+nxfbin)

   ! Pass through the output from the equation of state calculation
   if (present(eosout)) eosout = eos
   if (present(abundout)) abundout = abund

   ! Store variables for explicit calculations
   if (ji <= 0) then
     expl_var(jk, explv_avmu, Jstar) = avmu_new
     expl_var(jk, explv_logp, Jstar) = ap
   end if

   if ( ji < 0 ) return    ! Only do one star
   end do                  ! Continue with *2, if ktw = 2 (TWIN mode)
end subroutine funcs1




subroutine equns1 ( jk, kl, kq, fn2, equv )
   use real_kind
   use mesh
   use mesh_enc
   use extra_elements
   use constants
   use settings
   use control
   use semi_implicit_variables
   use current_model_properties
   use accretion_abundances
   
   implicit none
   integer, intent(in) :: jk, kl, kq
   real(double), intent(in) :: fn2(3, NFUNC)
   real(double), intent(out) :: equv(NEQ)

   integer :: Jstar,ii
   real(double) :: dmu12,dmu23,dw12,dw23,rich12,rich23,ddsi12,ddsi23,ton
   real(double) :: ris12,ris23,dssi12,dssi23,des12,des23,mt_smooth,pstv
   real(double) :: s12,s23,vgsf12,vgsf23,ves12,ves23,hp12,hp23
   real(double) :: wta,wtb,wtc,wtd
   
   ! Difference-equation terms:
   real(double) :: bcp(3), bct(3), vp(3), vpk(3), vr(3), vrk(3)
   real(double) :: vt(3), vtk(3), L(3), lk(3), lq(3), mt(3), vm(3), vmk(3), sg(3)
   real(double) :: wt(3), x1(3), x1t(3), x16(3), x16t(3), x4(3), x4t(3), x12(3)
   real(double) :: x12t(3), x20(3), x20t(3), bcm(3), vi(3), vik(3), phi(3), phik(3)
   real(double) :: bcf(3), bcs(3), bcph(3), x14(3), x14t(3), avmu(3), sgth(3), omega(3)
   real(double) :: omegat(3), sgam(3), si(3)
   real(double) :: bca(3), bce(3), xi(3), xik(3), dlrk(3), bcmb(3), x24(3), x24t(3)
   real(double) :: x28(3), x28t(3), x56(3), x56t(3)
   
   !Common blocks:
   ! Extra variables
   real(double) :: xvars(3, nsfunc+1:nsfunc+nxfstar), xvarb(3, nsfunc+2*nxfstar+1:nsfunc+2*nxfstar+nxfbin)

   ! equations and extra equations
   real(double) :: equ(neq)

   ! VAR(3),(2),(1): values at current, previous and anteprevious meshpoints.
   ! (3) is nearest centre if KL=0, KQ=1; nearest surface if KL=1, KQ=-1
   
   ! Initialise local variables
   ves12 = 0.0d0
   hp12 = 0.0d0
   ves23 = 0.0d0
   hp23 = 0.0d0

   do Jstar = 1, ktw
   ! Copy input variables, star 1 or 2
   bcp(:)   = fn2(:, 58*(Jstar-1) + 1)
   bct(:)   = fn2(:, 58*(Jstar-1) + 2)
   vp(:)    = fn2(:, 58*(Jstar-1) + 3)
   vpk(:)   = fn2(:, 58*(Jstar-1) + 4)
   vr(:)    = fn2(:, 58*(Jstar-1) + 5)
   vrk(:)   = fn2(:, 58*(Jstar-1) + 6)
   vt(:)    = fn2(:, 58*(Jstar-1) + 7)
   vtk(:)   = fn2(:, 58*(Jstar-1) + 8)
   L(:)     = fn2(:, 58*(Jstar-1) + 9)
   lk(:)    = fn2(:, 58*(Jstar-1) + 10)
   lq(:)    = fn2(:, 58*(Jstar-1) + 11)
   mt(:)    = fn2(:, 58*(Jstar-1) + 12)
   vm(:)    = fn2(:, 58*(Jstar-1) + 13)
   vmk(:)   = fn2(:, 58*(Jstar-1) + 14)
   sg(:)    = fn2(:, 58*(Jstar-1) + 15)
   wt(:)    = fn2(:, 58*(Jstar-1) + 16)
   x1(:)    = fn2(:, 58*(Jstar-1) + 17)
   x1t(:)   = fn2(:, 58*(Jstar-1) + 18)
   x16(:)   = fn2(:, 58*(Jstar-1) + 19)
   x16t(:)  = fn2(:, 58*(Jstar-1) + 20)
   x4(:)    = fn2(:, 58*(Jstar-1) + 21)
   x4t(:)   = fn2(:, 58*(Jstar-1) + 22)
   x12(:)   = fn2(:, 58*(Jstar-1) + 23)
   x12t(:)  = fn2(:, 58*(Jstar-1) + 24)
   x20(:)   = fn2(:, 58*(Jstar-1) + 25)
   x20t(:)  = fn2(:, 58*(Jstar-1) + 26)
   bcm(:)   = fn2(:, 58*(Jstar-1) + 27)
   vi(:)    = fn2(:, 58*(Jstar-1) + 28)
   vik(:)   = fn2(:, 58*(Jstar-1) + 29)
   phi(:)   = fn2(:, 58*(Jstar-1) + 30)
   phik(:)  = fn2(:, 58*(Jstar-1) + 31)
   bcf(:)   = fn2(:, 58*(Jstar-1) + 32)
   bcs(:)   = fn2(:, 58*(Jstar-1) + 33)
   bcph(:)  = fn2(:, 58*(Jstar-1) + 34)
   x14(:)   = fn2(:, 58*(Jstar-1) + 35)
   x14t(:)  = fn2(:, 58*(Jstar-1) + 36)
   avmu(:)  = fn2(:, 58*(Jstar-1) + 37)
   sgth(:)  = fn2(:, 58*(Jstar-1) + 38)
   omega(:) = fn2(:, 58*(Jstar-1) + 39)
   omegat(:)= fn2(:, 58*(Jstar-1) + 40)
   sgam(:)  = fn2(:, 58*(Jstar-1) + 41)
   si(:)    = fn2(:, 58*(Jstar-1) + 42)

   !> \todo FIXME: magnesium values are stored in the same location for both stars
   !! (in the binary part), this needs to be fixed!
   !<
   x24(:)   = fn2(:, 49)
   x24t(:)  = fn2(:, 50)

   ! Input variables, binary orbit
   bca(:)   = fn2(:, 43)
   bce(:)   = fn2(:, 44)
   xi(:)    = fn2(:, 45)
   xik(:)   = fn2(:, 46)
   dlrk(:)  = fn2(:, 47)
   bcmb(:)  = fn2(:, 48)

   ! Extra variables
   ii = nsfunc + (Jstar-1)*nxfstar
   xvars(:, nsfunc+1:nsfunc+nxfstar) = fn2(:, ii+1:ii+nxfstar)
   ii = nsfunc + 2*nxfstar
   xvarb(:, ii+1:ii+nxfbin) = fn2(:, ii+1:ii+nxfbin)

   x28(:) = xvars(:, fx_x28)
   x28t(:) = xvars(:, fx_x28t)
   x56(:) = xvars(:, fx_x56)
   x56t(:) = xvars(:, fx_x56t)

   ! second-order difference equations at interior points
   if ( 3 <= jk + kl .and. jk + kl <= kh ) then

     ! Molecular weight gradient, for thermohaline mixing and angular momentum transport
     dmu12 = kq * (avmu(1)-avmu(2))
     dmu23 = kq * (avmu(2)-avmu(3))

     ! Gradient in rotational velocity, for rotational shear
     dw12 = (omega(2) - omega(1))**2
     dw23 = (omega(3) - omega(2))**2

     ! Richardson number (minus omega gradient), for shear instability
     rich12 = 0.5*(xvars(1, fx_rich)+xvars(2, fx_rich))
     rich23 = 0.5*(xvars(2, fx_rich)+xvars(3, fx_rich))

     ! Mixing coefficient for dynamical shear instability
     ddsi12 = 0.0d0
     ddsi23 = 0.0d0
     if ( dw12 > 0.0d0 .and. rich12>0.0d0 .and. rich12/dw12 < 1.0d0 ) then
        ton = (1.0d0 - rich12/dw12)**2   ! Turn on factor
        ddsi12 = 0.5d0*( xvars(1, fx_ddsi)+xvars(2, fx_ddsi) )*ton
     end if
     if ( dw23 > 0.0d0 .and. rich23>0.0d0 .and. rich23/dw23 < 1.0d0 ) then
        ton = (1.0d0 - rich23/dw23)**2   ! Turn on factor
        ddsi23 = 0.5d0*( xvars(3, fx_ddsi)+xvars(2, fx_ddsi) )*ton
     end if

     ! Mixing coefficients secular shear instability
     ris12 = 0.5d0*(xvars(1, fx_sssi)+xvars(2, fx_sssi))
     ris23 = 0.5d0*(xvars(2, fx_sssi)+xvars(3, fx_sssi))
     dssi12 = 0.0d0
     dssi23 = 0.0d0
     if ( dw12 > 0.0d0 .and. ris12>0.0d0 .and. ris12/dw12 < 1.0d0 ) then
        dssi12 = 0.5d0*(xvars(1, fx_dssi)+xvars(2, fx_dssi))*dw12
     end if
     if ( dw23 > 0.0d0 .and. ris23>0.0d0 .and. ris23/dw23 < 1.0d0 ) then
        dssi23 = 0.5d0*(xvars(2, fx_dssi)+xvars(3, fx_dssi))*dw23
     end if

     !HP12 = 0.5*(XVARS(1, FX_HP)+XVARS(2, FX_HP))
     !VES12 = 0.5*(XVARS(1, FX_VES)+XVARS(2, FX_VES))
     !VMU12 = 0.5*(XVARS(1, FX_VMU)+XVARS(2, FX_VMU)) * CFMU*DMU12
     !VES12 = PSTV(DABS(VES12) - DABS(VMU12), 0.0D0)
     !print *, JK, 1.0D-11*HP12*VES12, 1.0d22*DDSI12/XVARS(1, FX_SGF)

     !         ! Eddington-Sweet circulation, compute net circulation velocity
     !         HP12 = 0.5*(XVARS(1, FX_HP)*XVARS(1, FX_SGF)+XVARS(2, FX_HP)*XVARS(2, FX_SGF))
     !         HP23 = 0.5*(XVARS(2, FX_HP)*XVARS(2, FX_SGF)+XVARS(3, FX_HP)*XVARS(3, FX_SGF))
     !         VMU12 = 0.5*(XVARS(1, FX_VMU)+XVARS(2, FX_VMU)) * CFMU*DMU12
     !         VMU23 = 0.5*(XVARS(2, FX_VMU)+XVARS(3, FX_VMU)) * CFMU*DMU23
     !         VES12 = 0.5*(XVARS(1, FX_VES)+XVARS(2, FX_VES))
     !         VES23 = 0.5*(XVARS(2, FX_VES)+XVARS(3, FX_VES))
     !         VES12 = PSTV(DABS(VES12) - DABS(VMU12), 0.0D0)
     !         VES23 = PSTV(DABS(VES23) - DABS(VMU23), 0.0D0)

     ! Eddington-Sweet circulation, diffusion approximation following Heger.
     des12 = 0.5*(xvars(1, fx_des) + xvars(2, fx_des))
     des23 = 0.5*(xvars(2, fx_des) + xvars(3, fx_des))

     ! Alternatively, the expression from Maeder & Zahn (1998)
     ! Factor by which to multiply VES to include stabilising effect of
     ! mu gradient
     !         VESMU121 = 1.0d0/max(1.0d0, 1.0d0-XVARS(1, FX_DGMU)*DMU12)
     !         VESMU122 = 1.0d0/max(1.0d0, 1.0d0-XVARS(2, FX_DGMU)*DMU12)
     !         VESMU231 = 1.0d0/max(1.0d0, 1.0d0-XVARS(3, FX_DGMU)*DMU23)
     !         VESMU232 = 1.0d0/max(1.0d0, 1.0d0-XVARS(2, FX_DGMU)*DMU23)
     !         VESMU12 = 0.5d0*(VESMU121+VESMU122)
     !         VESMU23 = 0.5d0*(VESMU231+VESMU232)
     !
     !         DES12 = 0.5*(VESMU121*XVARS(1, FX_DES) + VESMU122*XVARS(2, FX_DES))
     !         DES23 = 0.5*(VESMU231*XVARS(3, FX_DES) + VESMU232*XVARS(2, FX_DES))

     ! Combined diffusion coefficients for chemical mixing
     ! Convection and thermohaline mixing
     mt_smooth = 0.0d0
     s12 = 0.5d0*(sg(1) + sg(2)) - mixing_fudge*pstv(kq*mt(2), mt_smooth)&
          + 0.5d0*(sgth(1)+sgth(2))*pstv(dmu12, 0.0d0)
     s23 = 0.5d0*(sg(2) + sg(3)) - mixing_fudge*pstv(-kq*mt(3), mt_smooth)&
          + 0.5d0*(sgth(2)+sgth(3))*pstv(dmu23, 0.0d0)

     ! Dynamical shear instability
     s12 = s12 + cfc*ddsi12
     s23 = s23 + cfc*ddsi23

     ! Secular shear instability
     s12 = s12 + cfc*dssi12
     s23 = s23 + cfc*dssi23

     ! Eddington-Sweet circulation
     s12 = s12 + cfc*des12
     s23 = s23 + cfc*des23

     equ(1)     = s23*(x1(3)  - x1(2))  - s12*(x1(2) - x1(1))  - x1t(2)
     equ(2)     = s23*(x16(3) - x16(2)) - s12*(x16(2) -x16(1)) - x16t(2)
     equ(3)     = s23*(x4(3)  - x4(2))  - s12*(x4(2) - x4(1))  - x4t(2)
     equ(4)     = s23*(x12(3) - x12(2)) - s12*(x12(2) -x12(1)) - x12t(2)
     equ(5)     = s23*(x20(3) - x20(2)) - s12*(x20(2) -x20(1)) - x20t(2)
     equ(en14)  = s23*(x14(3) - x14(2)) - s12*(x14(2) -x14(1)) - x14t(2)
     equ(emg24) = s23*(x24(3) - x24(2)) - s12*(x24(2) -x24(1)) - x24t(2)
     equ(esi28) = s23*(x28(3) - x28(2)) - s12*(x28(2) -x28(1)) - x28t(2)
     equ(efe56) = s23*(x56(3) - x56(2)) - s12*(x56(2) -x56(1)) - x56t(2)
     !print *, jk, x28(3), x28(2), x28t(2)

     if (apply_second_order_corrections) then
        ! Second order corrections to hydrogen equation
        equ(1) = equ(1) - ddx1(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(1) = equ(1) - ddx1(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(1) = equ(1) + ddx1(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(1) = equ(1) + ddx1(jk)    * pstv(-kq*mt(2), 0.0d0)
        ! Second order corrections to oxygen equation
        equ(2) = equ(2) - ddx16(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(2) = equ(2) - ddx16(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(2) = equ(2) + ddx16(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(2) = equ(2) + ddx16(jk)    * pstv(-kq*mt(2), 0.0d0)
        ! Second order corrections to helium equation
        equ(3) = equ(3) - ddx4(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(3) = equ(3) - ddx4(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(3) = equ(3) + ddx4(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(3) = equ(3) + ddx4(jk)    * pstv(-kq*mt(2), 0.0d0)
        ! Second order corrections to carbon equation
        equ(4) = equ(4) - ddx12(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(4) = equ(4) - ddx12(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(4) = equ(4) + ddx12(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(4) = equ(4) + ddx12(jk)    * pstv(-kq*mt(2), 0.0d0)
        ! Second order corrections to neon equation
        equ(5) = equ(5) - ddx20(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(5) = equ(5) - ddx20(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(5) = equ(5) + ddx20(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(5) = equ(5) + ddx20(jk)    * pstv(-kq*mt(2), 0.0d0)
        ! Second order corrections to nitrogen equation
        equ(en14) = equ(en14) - ddx14(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(en14) = equ(en14) - ddx14(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(en14) = equ(en14) + ddx14(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(en14) = equ(en14) + ddx14(jk)    * pstv(-kq*mt(2), 0.0d0)
        ! Second order corrections to magnesium equation
        equ(emg24) = equ(emg24) - ddx24(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(emg24) = equ(emg24) - ddx24(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(emg24) = equ(emg24) + ddx24(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(emg24) = equ(emg24) + ddx24(jk)    * pstv(-kq*mt(2), 0.0d0)
     end if

     ! Add advection terms (contributions from gravitational settling)
     if (cgrs > 0.0d0) then
        equ(1)     = equ(1)    - kq*(xvars(3, fx_fh)  - xvars(2, fx_fh))
        equ(2)     = equ(2)    - kq*(xvars(3, fx_fo) - xvars(2, fx_fo))
        equ(3)     = equ(3)    - kq*(xvars(3, fx_fhe)  - xvars(2, fx_fhe))
        equ(4)     = equ(4)    - kq*(xvars(3, fx_fc) - xvars(2, fx_fc))
        equ(5)     = equ(5)    - kq*(xvars(3, fx_fne) - xvars(2, fx_fne))
        equ(en14)  = equ(en14) - kq*(xvars(3, fx_fn) - xvars(2, fx_fn))
        equ(emg24)  = equ(emg24) - kq*(xvars(3, fx_fmg) - xvars(2, fx_fmg))
        equ(esi28)  = equ(esi28) - kq*(xvars(3, fx_fsi) - xvars(2, fx_fsi))
        equ(efe56)  = equ(efe56) - kq*(xvars(3, fx_ffe) - xvars(2, fx_ffe))
     end if
     equ(esumx) = equ(1)+equ(2)+equ(3)+equ(4)+equ(5)+equ(en14)+equ(emg24)+equ(esi28)+equ(efe56)

     ! Angular momentum transport
     ! The mixing coefficients here need to be weighed with the specific
     !  moment of inertia, so we need to recalculate some quantities
     mt_smooth = 1.0d-22
     mt_smooth = 1.0d-3*dabs(mt(2))
     mt_smooth = 0.0d0
     s12 = 0.5d0*(sgam(1)*si(1) + sgam(2)*si(2)) - 1.0d0*pstv(kq*mt(2), mt_smooth)*si(2)
     s23 = 0.5d0*(sgam(2)*si(2) + sgam(3)*si(3)) - 1.0d0*pstv(-kq*mt(3), mt_smooth)*si(2)
     !S12 = 0.5D0*(SGAM(1)*SI(1) + SGAM(2)*SI(2)) - MT(2)*SI(2)
     !S23 = 0.5D0*(SGAM(2)*SI(2) + SGAM(3)*SI(3)) - MT(3)*SI(3)
     !S12 = ( 0.5D0*(SGAM(1) + SGAM(2)) - PSTV(KQ*MT(2), 0.0D0) )*SI(2)
     !S23 = ( 0.5D0*(SGAM(2) + SGAM(3)) - PSTV(-KQ*MT(3), 0.0D0) )*SI(2)
     ! Dynamical shear instability
     if ( dw12 > 0.0d0 .and. rich12>0.0d0 .and. rich12/dw12 < 1.0d0 ) then
        ton = (1.0d0 - rich12/dw12)**2   ! Turn on factor
        ddsi12 = 0.5d0*( xvars(1, fx_ddsi)*si(1)+xvars(2, fx_ddsi)*si(2) )*ton
        !DDSI12 = 0.5D0*( XVARS(1, FX_DDSI)+XVARS(2, FX_DDSI) )*SI(2)*TON
     end if
     if ( dw23 > 0.0d0 .and. rich23>0.0d0 .and. rich23/dw23 < 1.0d0 ) then
        ton = (1.0d0 - rich23/dw23)**2   ! Turn on factor
        ddsi23 = 0.5d0*( xvars(3, fx_ddsi)*si(3)+xvars(2, fx_ddsi)*si(2) )*ton
        !DDSI23 = 0.5D0*( XVARS(3, FX_DDSI)+XVARS(2, FX_DDSI) )*SI(2)*TON
     end if
     s12 = s12 + ddsi12
     s23 = s23 + ddsi23

     ! Secular shear instability
     dssi12 = 0.0d0
     dssi23 = 0.0d0
     if ( dw12 > 0.0d0 .and. ris12>0.0d0 .and. ris12/dw12 < 1.0d0 ) then
        dssi12 = 0.5d0*(xvars(1, fx_dssi)*si(1)+xvars(2, fx_dssi)*si(2))*dw12
     end if
     if ( dw23 > 0.0d0 .and. ris23>0.0d0 .and. ris23/dw23 < 1.0d0 ) then
        dssi23 = 0.5d0*(xvars(2, fx_dssi)*si(2)+xvars(3, fx_dssi)*si(3))*dw23
     end if
     !print *, JK, DSSI12, DSSI23
     s12 = s12 + dssi12
     s23 = s23 + dssi23

     ! Eddington-Sweet circulation:
     !> \todo FIXME: use the diffusion coefficient exported from funcs1, as
     !! for chemical transport - or use the "proper" advection treatment.
     !<
     ! HP12 = 0.5*(XVARS(1, FX_HP)*XVARS(1, FX_SGF)*SI(1)+XVARS(2, FX_HP)*XVARS(2, FX_SGF)*SI(2))
     ! HP23 = 0.5*(XVARS(2, FX_HP)*XVARS(2, FX_SGF)*SI(2)+XVARS(3, FX_HP)*XVARS(3, FX_SGF)*SI(3))
     ! S12 = S12 + 1.d-33*CESC*VES12 * HP12
     ! S23 = S23 + 1.d-33*CESC*VES23 * HP23

     ! Eddington-Sweet circulation, diffusion approximation following Heger.
     des12 = 0.5*(xvars(1, fx_des)*si(1) + xvars(2, fx_des)*si(2))
     des23 = 0.5*(xvars(2, fx_des)*si(2) + xvars(3, fx_des)*si(3))
     s12 = s12 + des12
     s23 = s23 + des23


     ! Goldreich-Schubert-Fricke instability
     !> \todo FIXME
     !! The mixing coefficient is similar to that of the Eddington-Sweet
     !!  circulation, and has a similar order of magnitude. For now, we take
     !!  the easy way out and just apply the Eddington-Sweet term again
     !<
     ! Stability condition: j=i omega increases outward
     vgsf12 = ves12 * (omega(2) - omega(1))/(vr(2) - vr(1))
     vgsf23 = ves23 * (omega(3) - omega(2))/(vr(3) - vr(2))
     if ( (si(2+kl)*omega(2+kl) - si(3-kl)*omega(3-kl)) < 0) then
        s23 = s23 + 1.d-33*cgsf*vgsf23 * hp23
     end if
     if ( (si(1+kl)*omega(1+kl) - si(2-kl)*omega(2-kl)) < 0) then
        s12 = s12 + 1.d-33*cgsf*vgsf12 * hp12
     end if

     ! Diffusion equation for angular momentum
     equ(eamt) = s23*(omega(3) - omega(2)) - s12*(omega(2) - omega(1)) - omegat(2)
   end if
   ! next-to-surface boundary conditions for second-order equations
   if ( jk + kl == 2 ) then
     ! Modified equations for accreting of matter with different composition

     ! Thermohaline mixing of surface layer
     ! Accretion abundance set in FUNCS1
     !X1AC =    0.598490994028748
     !X4AC =    0.381392784203137
     !X16AC =   1.009965586641970E-002
     !X12AC =   3.539808668161780E-003
     !X20AC =   1.851094969072195E-003
     !X14AC =   1.045628497798272E-003
     !X24AC =   6.838266745312392E-004

     if (kl == 0) then    ! (3) is nearest to centre, so we have the grid points | 3 | 2 | AC |
        s12 = 0.5d0*(sg(1) + sg(2)) - mixing_fudge*pstv(kq*mt(2), 0.0d0)
        s23 = 0.5d0*(sg(2) + sg(3)) - mixing_fudge*pstv(-kq*mt(3), 0.0d0)
        equ(1)     = s23*(x1(3)  - x1(2))  - s12*(x1(2) - xac(1,Jstar))   - x1t(2)
        equ(2)     = s23*(x16(3) - x16(2)) - s12*(x16(2) -xac(5,Jstar))  - x16t(2)
        equ(3)     = s23*(x4(3)  - x4(2))  - s12*(x4(2) - xac(2,Jstar))   - x4t(2)
        equ(4)     = s23*(x12(3) - x12(2)) - s12*(x12(2) -xac(3,Jstar))  - x12t(2)
        equ(5)     = s23*(x20(3) - x20(2)) - s12*(x20(2) -xac(6,Jstar))  - x20t(2)
        equ(en14)  = s23*(x14(3) - x14(2)) - s12*(x14(2) -xac(4,Jstar))  - x14t(2)
        equ(emg24) = s23*(x24(3) - x24(2)) - s12*(x24(2) -xac(7,Jstar))  - x24t(2)
        equ(esi28) = s23*(x28(3) - x28(2)) - s12*(x28(2) -xac(8,Jstar)) - x28t(2)
        equ(efe56) = s23*(x56(3) - x56(2)) - s12*(x56(2) -xac(9,Jstar)) - x56t(2)
     else                 ! (3) is nearest to surface, we have the grid points | 1 | 2 | 3 | AC |
        s12 = 0.5d0*(sg(2) + sg(3)) - pstv(kq*mt(3), 0.0d0)&
             + 0.5d0*(sgth(2)+sgth(3))*pstv(kq * (avmu(2)-avmu(3)), 0.0d0)
        s23 = -xvars(3, fx_macc)

        equ(1)     = s23*(xac(1,Jstar)   - x1(3))  - s12*(x1(3) - x1(2))  - x1t(3)
        equ(2)     = s23*(xac(5,Jstar)  - x16(3)) - s12*(x16(3) -x16(2)) - x16t(3)
        equ(3)     = s23*(xac(2,Jstar)   - x4(3))  - s12*(x4(3) - x4(2))  - x4t(3)
        equ(4)     = s23*(xac(3,Jstar)  - x12(3)) - s12*(x12(3) -x12(2)) - x12t(3)
        equ(5)     = s23*(xac(6,Jstar)  - x20(3)) - s12*(x20(3) -x20(2)) - x20t(3)
        equ(en14)  = s23*(xac(4,Jstar)  - x14(3)) - s12*(x14(3) -x14(2)) - x14t(3)
        equ(emg24) = s23*(xac(7,Jstar)  - x24(3)) - s12*(x24(3) -x24(2)) - x24t(3)
        equ(esi28) = s23*(xac(8,Jstar)  - x28(3)) - s12*(x28(3) -x28(2)) - x28t(3)
        equ(efe56) = s23*(xac(9,Jstar)  - x56(3)) - s12*(x56(3) -x56(2)) - x56t(3)
     end if

     ! Advection terms (from gravitational settling)
     if (cgrs > 0.0d0) then
        equ(1)     = equ(1)    - kq*(xvars(3, fx_fh)  - xvars(2, fx_fh))
        equ(2)     = equ(2)    - kq*(xvars(3, fx_fo) - xvars(2, fx_fo))
        equ(3)     = equ(3)    - kq*(xvars(3, fx_fhe)  - xvars(2, fx_fhe))
        equ(4)     = equ(4)    - kq*(xvars(3, fx_fc) - xvars(2, fx_fc))
        equ(5)     = equ(5)    - kq*(xvars(3, fx_fne) - xvars(2, fx_fne))
        equ(en14)  = equ(en14) - kq*(xvars(3, fx_fn) - xvars(2, fx_fn))
        equ(emg24)  = equ(emg24) - kq*(xvars(3, fx_fmg) - xvars(2, fx_fmg))
        equ(esi28)  = equ(esi28) - kq*(xvars(3, fx_fsi) - xvars(2, fx_fsi))
        equ(efe56)  = equ(efe56) - kq*(xvars(3, fx_ffe) - xvars(2, fx_ffe))
     end if

     equ(esumx) = equ(1)+equ(2)+equ(3)+equ(4)+equ(5)+equ(en14)+equ(emg24)+equ(esi28)+equ(efe56)

     ! Angular momentum transport boundary condition (use normal BC ?)
     s12 = 0.5d0*(sgam(1)*si(1) + sgam(2)*si(2)) - pstv(kq*mt(2), 0.0d0)*si(2)
     s23 = 0.5d0*(sgam(2)*si(2) + sgam(3)*si(3)) - pstv(-kq*mt(3), 0.0d0)*si(2)
     ! Dynamical shear instability
     !> \todo FIXME: add dynamical shear instability for surface layer
     ! Always need at least a little bit of mixing to correlate meshpoints
     equ(eamt) = s23*(omega(3) - omega(2)) - s12*(omega(2) - omega(1)) - omegat(2)
     !equ(EAMT) = -S23*(OMEGA(3) - OMEGA(2)) - OMEGAT(3)
     !equ(EAMT) = OMEGAT(2)
     !print *, equ(EAMT), S12, S23
     !print *, OMEGAT(3), OMEGA(3), OMEGA(2)
     !print *, OMEGAT(3), SI(3)*OMEGA(3)
     !print *, MT(1), MT(2), MT(3)
     equ(eamt) = omega(3) - omega(2)
   end if
   ! next-to-central boundary conditions for second-order equations
   if ( jk + kl == kh + 1 ) then
     s23 = kq*0.5d0*(sg(2)+sg(3)) &
          + kq*0.5d0*(sgth(2)+sgth(3))*pstv(kq * (avmu(2)-avmu(3)), 0.0d0)
     equ(1)     = s23*(x1(3)  - x1(2))  + x1t(3 - kl)
     equ(2)     = s23*(x16(3) - x16(2)) + x16t(3 - kl)
     equ(3)     = s23*(x4(3)  - x4(2))  + x4t(3 - kl)
     equ(4)     = s23*(x12(3) - x12(2)) + x12t(3 - kl)
     equ(5)     = s23*(x20(3) - x20(2)) + x20t(3 - kl)
     equ(en14)  = s23*(x14(3) - x14(2)) + x14t(3 - kl)
     equ(emg24) = s23*(x24(3) - x24(2)) + x24t(3 - kl)
     equ(esi28) = s23*(x28(3) - x28(2)) + x28t(3 - kl)
     equ(efe56) = s23*(x56(3) - x56(2)) + x56t(3 - kl)
     equ(esumx) = equ(1)+equ(2)+equ(3)+equ(4)+equ(5)+equ(en14)+equ(emg24)+equ(esi28)+equ(efe56)

     ! Angular momentum transport
     s23 = kq*0.5d0*(sgam(2)*si(2) + sgam(3)*si(3))
     !> \todo FIXME: add other mixing coefficients for central boundary condition
     equ(eamt) = s23*(omega(3) - omega(2)) + omegat(3 - kl)
   end if

   ! First-order difference equations at interior points.
   ! KL=0, KQ=1 means that XX(3) is deeper in star than XX(2); KL=1, KQ=-1
   ! means nearer surface. XXK is sort-of dXX/dK, except that its sign is
   ! *independent* of KL, KQ, so a factor of KQ is put in.
   if ( 2.le.jk .and. jk.le.kh ) then
     wta = 0.5d0*kq
     wtb = 0.5d0*(wt(2) + wt(3))
     wtc = wta / (1.0 + wtb)
     if (dabs(wtc) < 1.0d-3) wtc = 0.0d0 ! Snap when value becomes small
     wtd = kq - wtc
     equ(6) = vp(3) - vp(2) - wtc*vpk(3 - kl) - wtd*vpk(2 + kl)
     equ(7) = vr(3) - vr(2) - wtd*vrk(3 - kl) - wtc*vrk(2 + kl)
     equ(8) = vt(3) - vt(2) - wtc*vtk(3 - kl) - wtd*vtk(2 + kl)
     ! Luminosity equn: MT gives advective term, DLRK heat transfer in contact
     ii = 2*Jstar - 3
     equ(9) = L(3) - L(2) - wtd*lk(3 - kl) - wtc*lk(2 + kl) &
          - lq(2)*pstv(kq*mt(2),0.0d0) + lq(3)*pstv(-kq*mt(3),0.0d0)&
          + ii*kq*dlrk(2 + kl)
     equ(10) = vm(3) - vm(2) - wta*(vmk(3) + vmk(2))
     equ(11) = vi(3) - vi(2) - wta*(vik(3) + vik(2))
     equ(12) = phi(3) - phi(2) - wta*(phik(3) + phik(2))
     equ(19) = xi(3) - xi(2) - wta*(xik(3) + xik(2))
     equ(etam) = xvars(3, fx_am) - xvars(2, fx_am) - wta*(xvars(3, fx_amk) + xvars(2, fx_amk))

     if (apply_second_order_corrections) then
        ! Second order corrections to luminosity equation
        equ(9) = equ(9) - ddl(jk-kq) * pstv( kq*mt(2), 0.0d0)
        equ(9) = equ(9) - ddl(jk+kq) * pstv(-kq*mt(3), 0.0d0)
        equ(9) = equ(9) + ddl(jk)    * pstv( kq*mt(3), 0.0d0)
        equ(9) = equ(9) + ddl(jk)    * pstv(-kq*mt(2), 0.0d0)
     end if
   end if
   ! surface boundary conditions for first-order equations and `eigenvalues'
   if ( jk == 1 ) then
     equ(6) = bcm(3)
     equ(7) = bcp(3)
     equ(8) = bct(3)
     equ(9) = bcf(3)
     equ(10) = bcs(3)
     equ(11) = bcph(3)
     equ(17) = bca(3)
     equ(18) = bce(3)
     equ(20) = bcmb(3)
   end if
   ! central boundary conditions for first-order equations
   if ( jk == kh + 1 ) then
     equ(6) = vm(3)
     equ(7) = L(3)
     equ(8) = vr(3)
     equ(9) = vi(3)
     equ(19) = xi(3)
     equ(etam) = xvars(3, fx_am)
   end if

   ! Copy output back to INE

   ! Equations, star 1 or 2
   ii = 24*(Jstar-1)
   equv(ii+1:ii+16) = equ(1:16)

   ! Equations, binary orbit
   equv(17:24) = equ(17:24)

   ! Extra equations
   ii = nseq+(Jstar-1)*nxestar
   equv(ii+1:ii+nxestar) = equ(ii+1:ii+nxestar)
   ii = nseq+2*nxestar
   equv(ii+1:ii+nxebin) = equ(ii+1:ii+nxebin)

   end do
end subroutine equns1



!> Solve for dimensionless L1, L2 potentials
subroutine potent ( qq, dphi )
   use real_kind
   implicit none
   real(double) :: qq,dphi
   integer :: ij
   real(double) :: q,q1,q2,q3,q4,q5,q6,xl1,xl2

   q = qq
   if ( q < 1.0d0 ) q = 1.0d0/qq
   q5 = 1.0d0 + q
   xl1 = (q/(3.0d0*q5))**0.333333333d0
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
   end do
   xl2 = xl2 + (q - xl2*(q2 - xl2*(q1 - xl2*(q3 - xl2*(q4 - xl2*q5)))&
       ))/(q2 - xl2*(2.0d0*q1 - xl2*(3.0d0*q3 - xl2*(4.0d0*q4 &
       - xl2*5.0d0*q5))))
   dphi = q*q6/xl1 + q6/(1.0d0 - xl1) + 0.5d0*(xl1 - q6)**2&
       - (q*q6/xl2 + q6/(xl2 - 1.0d0) + 0.5d0*(xl2 - q6)**2)
end subroutine potent


!> Compute the cubic root 
function cbrt(x)
   use real_kind
   use constants
   
   implicit none
   real(double) :: cbrt,x
   
   cbrt = abs(x)**c3rd
end function cbrt


!> Compute the Roche-lobe radius
function rlobe(x)
   use real_kind
   use constants
   
   implicit none
   real(double) :: rlobe,x,x2
   
   x2 = x**2
   rlobe = 0.49_dbl*x2 / (0.6_dbl*x2 + log(1.0_dbl + x))
end function rlobe

