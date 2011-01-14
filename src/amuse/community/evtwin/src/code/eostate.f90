! STATEF: calculate equation of state variables
! AF is the logarithm of the variable f that parametrizes the el. degeneracy
! AT is the logarithm of the temperature
! Composition variables are passed using a common block.
subroutine statef ( af, at, xa, abund, eos )
   use real_kind
   use constants
   use settings
   use opacity_co
   use neutrinos
   use atomic_data
   use eostate_types
   
   implicit none
   real(double), intent(in) :: af, at
   real(double), intent(in) :: xa(9)
   type(eostate), intent(out) :: eos
   type(abundance), intent(out) :: abund
   real(double), parameter :: eps = 1.0d-30

   integer :: i,j
   real(double) :: f,wf,psi,g,se,sef,set,ue,sp,up,spt,spf
   real(double) :: xip,paf,puf,pwf,ff,ti,de,dc,dvt,dvf,dpa,dpat,dpaf
   real(double) :: dsat,dsaf,dua,dv,nef,si,sif,sit,ui,vm,sha
   real(double) :: dsa,net,vx,sjha,scha,nx,nxf,nxt,dh2ti,expdh2ti
   real(double) :: d1ti,d2ti,d3ti,zet,dzh2t,dzh2tt,zh2,zh2t,zh2tt,zh2s
   real(double) :: h2a,h2bt,h2at,qa,qb,qc,hg,hi,np,h2,ni,qd,qaf,qat,qf,qbf,qbt
   real(double) :: hgf,hgt,hif,hit,h2f,h2t,db,dl,dpb,dpbt,dpbf,dsb,dsbt,dsbf,dub
   real(double) :: tcr,p0,pi,t2,t4,b,dsf,dst,q,scv,gamma,tst,frho,tf
   real(double) :: llambda,mu_plasma,mu_rad,nu
   real(double) :: find_positron_f,dvt_over_dvf,scha_over_sha,sjha_over_sha,get_opacity
   real(double) :: inv_ne
   real(double) :: ha(26), va(26)

   real(double) :: ap, arho, u, p, rho, fk, T, sf, st, zt, grada
   real(double) :: scp, rf, rt, xhi, s, pr, pg, pf, pt, uf, ut, en
   real(double) :: wmu, delta, phi, fkt, fkr, prandtl

   real(double) :: fdi(16)
   real(double) :: re, pe, qe, ret, pet, qet, ref, pef, qef, xtt, xft, xff, xttt, xftt, xfft, xfff
   real(double) :: rp, pp, qp, rpt, ppt, qpt, rpf, ppf, qpf

   real(double) :: na(9), neo, nio, nzz, avm, ne, xai(9, 26)
   real(double) :: lna(9)           ! log(na), used several times
   real(double) :: mu(9)            ! Chemical potentials (in units of kT)
   real(double) :: zeta(9)          ! Average charge for atomic species i
   real(double) :: dlrdni(9)        ! Derivative of density, -(dlnrho/dNi)_P,T
   real(double) :: dsdnip(9)        ! Derivative of entropy, (dS/dNi)_p
   real(double) :: dsdnif(9)        ! Derivative of entropy, (dS/dNi)_f
   real(double) :: nif(9)           ! dNi/dlnf, for atomic species i
   real(double) :: nit(9)           ! dNi/dlnT, for atomic species i

   ! Calculate F and G parameters, needed to calculate Fermi-Dirac integrals
   f = dexp(af)
   T = dexp(at)
   !T = MAX(1000.0, T)         ! Don't trust T < 1000K
   uf = f/(1.0d0 + f)
   wf = dsqrt(1.0d0 + f)
   psi = af + 2.0d0*(wf - dlog(1.0d0 + wf))
   g = cte*T*wf

   ! Evaluate Fermi-Dirac integrals according to Eggleton, Faulkner &
   ! Flannery (1973): re, pe, se and ue correspond to rho*, p*, s* and u*.
   ! psi is the usual electron-degeneracy parameter
   call fdirac ( f, g, fdi )
   re   = fdi(1); pe   = fdi(2); qe   = fdi(3)
   ret  = fdi(4); pet  = fdi(5); qet  = fdi(6)
   ref  = fdi(7); pef  = fdi(8); qef  = fdi(9)
   xtt  = fdi(10); xft  = fdi(11); xff  = fdi(12)
   xttt = fdi(13); xftt = fdi(14); xfft = fdi(15)
   xfff = fdi(16)

   pe = g*pe
   pet = pet + 1.0d0
   pef = pef + 0.5d0*uf
   qe = qe/(re*wf)
   se = qe + 2.0d0*wf - psi
   sef = qe*(qef - ref - 0.5d0*uf) - 1.0d0/wf
   set = qe*(qet - ret)
   ue = se + psi - pe/(re*cte*T)

   ! Contributions from positrons
   rp = 0.0d0
   pp = 0.0d0
   sp = 0.0d0
   up = 0.0d0
   rpt = 0.0d0
   ppt = 0.0d0
   spt = 0.0d0
   rpf = 0.0d0
   ppf = 0.0d0
   spf = 0.0d0
   xip = -psi - 2.0 / (cte*T)
   if (eos_include_pairproduction .and. (xip > -15.0)) then
      !        Calculate effective degeneracy parameter for positrons
      !        See also Cox & Giuli ¤24.9
      !        NB: PF and PG are later returned as dP/dlogf and Pgas
      paf = find_positron_f(xip)
      pf = dexp(paf)
      puf = pf/(1.0d0 + pf)
      pwf = dsqrt(1.0d0 + pf)
      pg = cte*T*pwf

      call positron_fdirac(pf, pg, fdi)
      rp   = fdi(1); pp   = fdi(2); qp   = fdi(3)
      rpt  = fdi(4); ppt  = fdi(5); qpt  = fdi(6)
      rpf  = fdi(7); ppf  = fdi(8); qpf  = fdi(9)

      !        Derivatives are now known with respect to the *positron* degeneracy
      !        parameter pf, but we need to know them with respect to the *electron*
      !        degeneracy parameter f. The conversion factor is easy, however,
      !        FF = dlog pf / dlog f = - sqrt( (1+f)/(1+pf) )
      !        (because dxi/dxip = -1)
      if (rp > 0.0d0) then
         ff = -wf/pwf

         rpf = rpf * ff
         ppf = ppf * ff
         qpf = qpf * ff

         !        Derivatives with respect to temperature get an extra term because
         !        the effective positron degeneracy parameter depends on the temperature
         !        as well as on the electron degeneracy parameter.
         rpt = rpt - 2.0d0 * rpf / ( cte*T*g )
         ppt = ppt - 2.0d0 * ppf / ( cte*T*g )
         qpt = qpt - 2.0d0 * qpf / ( cte*T*g )

         pp = pg*pp
         ppt = ppt + 1.0d0  + puf / pg
         ppf = ppf + 0.5d0 * puf * ff
         qp = qp/(rp*pwf)
         sp = qp + 2.0d0*pwf - xip
         spf = qp*(qpf - rpf - 0.5d0*puf) + 1.0d0/wf
         spt = qp*(qpt - rpt)
         up = sp + xip - pp/(rp*cte*T)
      end if
   end if

   ! Evaluate some quantities that do not depend on the state of ionization:
   ! the na are the element number densities (per baryon); avm is the average
   ! mass per baryon (in amu); neo and nio are the numbers of electrons and
   ! ions assuming complete ionization.
   na(1:9) = xa(1:9)/cbn(1:9)
   nio = sum(na(1:9))
   neo = dot_product(dkzn(1:9), na(1:9))
   nzz = dot_product(dkzn(1:9)**2, na(1:9))
   avm = dot_product(can(1:9), na(1:9))

   ! Pressi gives a crude model for pressure ionization and
   ! a model for coulomb interactions, returning corrections to the electron
   ! chemical potential, pressure, entropy and internal energy.
   ! ti is 1ev/kt, de is 1 amu * number of electrons/cm3
   ti = cevb/T
   de = re*cd
   call pressi ( 1, ti, psi, neo, nio, nzz, de, ref, ret, xtt, xft, xff, xftt, xfft, xfff, f, &
                  dc, dvt, dvf, dpa, dpat, dpaf, dsa, dsat, dsaf, dua )
   de = max(re-rp, 1.0d-50)*cd
   dv = dc - psi
   dvf = dvf - wf

   lna = log(na + eps)
   nef = 0.0d0
   net = 0.0d0
   si = 0.0d0
   sif = 0.0d0
   sit = 0.0d0
   ui = 0.0d0
   vm = 0.0d0
   mu = 0.0d0
   nif = 0.0d0
   nit = 0.0d0
   dsdnip = 0.0d0
   dsdnif = 0.0d0
   ! Contributions of the completely ionized species
   ne = dot_product(dkzn(kion+1:9), na(kion+1:9))
   si = -dot_product(na(kion+1:9), lna(kion+1:9) - lcan(kion+1:9))
   dsdnif(kion+1:9) = -(lna(kion+1:9) - lcan(kion+1:9))
   mu(kion+1:9) = -dsdnif(kion+1:9) + chi(1, kion+1:9)*ti
   zeta(kion+1:9) = dkzn(kion+1:9)

   ! We don't really care about the ionisation state of all other elements,
   ! but we need it for some calculations (gravitational settling,
   ! radiative levitation), so we need to calculate what the ionisation
   ! fractions would be anyway. We use a flag to optionally skip this part
   ! of the calculation.
   if (.false.) then
      do i = kion+1, 9
         va(1) = -chi(1,i)*ti + dv
         ha(1) = fxp(va(1))*com(kzn(i))/com(kzn(i)+1)
         sha = 1.0d0 + ha(1)
         do j = 2, kzn(i)
            va(j) = -chi(j,i)*ti + j*dv
            ha(j) = ha(j-1)*fxp(va(j) - va(j-1)) * com(kzn(i)+1-j)/com(kzn(i)+2-j)
            sha = sha + ha(j)
         end do
         vx = na(i)/sha          ! Number density of ground state per baryon
         do j = 1, kzn(i)
            xai(i,j) = ha(j)*vx  ! Number density of state J per baryon
         end do
      end do
   else
      forall (i = kion+1:9)
         xai(i, kzn(i)-1) = 0.0d0
         xai(i, kzn(i)) = na(i)
      end forall
   end if
   ! Calculate ionization of the first kion elements, in reverse order so
   ! hydrogen is last (the H calculation is continued after the loop ends)
   dvt_over_dvf = dvt / dvf
   do i = kion, 1, -1
      ! Compute potentials va and number ratios ha of ionization state j
      ! relative to the ground state 1
      zeta(i) = 0.0d0
      va(1) = -chi(1,i)*ti + dv
      ha(1) = fxp(va(1))*com(kzn(i))/com(kzn(i)+1)
      sha = 1.0d0 + ha(1)
      sjha = ha(1)
      scha = chi(1,i)*ti*ha(1)
      do j = 2, kzn(i)
         va(j) = -chi(j,i)*ti + j*dv
         ha(j) = ha(j-1)*fxp(va(j) - va(j-1)) * com(kzn(i)+1-j) / com(kzn(i)+2-j)
         sha = sha + ha(j)
         sjha = sjha + j*ha(j)
         scha = scha + chi(j,i)*ti*ha(j)
         zeta(i) = zeta(i) + j * ha(j)
      end do
      scha_over_sha = scha / sha
      sjha_over_sha = sjha / sha
      si = si + na(i)*lcan(i)
      dsdnif(i) = dsdnif(i) + lcan(i)
      mu(i) = mu(i) - lcan(i)
      zeta(i) = zeta(i) / ( sha + eps )
      if ( i > 1 ) then
         ! Contributions to electron number density ne, entropy si and
         ! internal energy ui for helium and heavier
         vx = na(i)/sha
         si = si - vx*dlog(vx/com(kzn(i)+1) + eps)

         do j = 1, kzn(i)
            nx = ha(j)*vx        ! Number in state J per baryon
            xai(i,j) = nx        ! Store
            nxf = nx*dvf*(j - sjha_over_sha)
            nxt = nxf*dvt_over_dvf + nx*(chi(j,i)*ti - scha_over_sha)
            ne = ne + j*nx
            nef = nef + j*nxf
            net = net + j*nxt
            sif = sif - va(j)*nxf
            sit = sit - va(j)*nxt
            si = si - nx*dlog(max(nx/com(kzn(i) + 1 - j), eps))
            dsdnif(i) = dsdnif(i) - ha(j)/sha * dlog(max(nx/com(kzn(i) + 1 - j), eps))
            mu(i) = mu(i) + ha(j)/sha * ( log(ha(j)/sha + eps) )
            ui = ui + chi(j,i)*nx
            nif(i) = nxf
            nit(i) = nxt
         end do
      end if
   end do

   ! Set derivatives to 0 for species that are not present.
   ! Not entirely correct, but the easiest way to avoid numerical problems
   where (xa <= 0.0d0)
      dsdnif = 0.0d0
      dsdnip = 0.0d0
   end where

   ! Ionization and molecular dissciation of hydrogen.
   ! Partition function for H2 from Vardya (1960), Webbink (1975)
   vm = can(1)*dsqrt(can(1))
   dh2ti = ch2(1)*ti
   expdh2ti = dexp(-min(dh2ti, 100.0d0))
   d1ti = ch2(2)*ti
   d2ti = (ch2(3)*ti)**2
   d3ti = (ch2(4)*ti)**3
   zet = 1.0d0 - (1.0d0 + dh2ti)*expdh2ti
   dzh2t = -dh2ti**2*expdh2ti/zet
   dzh2tt = (dh2ti - 2.0d0 - dzh2t)*dzh2t
   zh2 = 6608.8d0*zet*dh2ti**(-2.5d0)*dexp(-d1ti - d2ti - d3ti)
   zh2 = min(zh2, 1.0d100)       ! Avoid overflows later
   zh2t = 2.5d0 + d1ti + 2.0d0*d2ti + 3.0d0*d3ti + dzh2t
   zh2tt = -d1ti - 4.0d0*d2ti - 9.0d0*d3ti + dzh2tt
   zh2s = zh2*dsqrt(8.0d0)/vm
   h2a = cen*(zh2s/4.0d0)*de/(T*dsqrt(T))/expdh2ti
   h2bt = dh2ti + 1.5d0 - zh2t
   h2at = ret - h2bt

   ! Solve for densities of H+, H, and H2
   qa = 2.0d0*h2a + ha(1)*(1.0d0 + ha(1))
   qb = ne + ha(1)*(ne - na(1))
   qc = na(1)*ne
   hg = 2.0d0*qc/(dsqrt(qb*qb + 4.0d0*qa*qc) + qb)
   hi = ha(1)*hg + tiny(0.0_dbl)    ! Make sure it's never 0
   xai(1,1) = hi                    ! Number density of protons (ionised H)
   ne = ne + hi                     ! Ionisation electrons / baryon
   np = rp / dabs(re - rp) * ne     ! Positrons / baryon
   inv_ne = 1.0d0/ne
   h2 = h2a*hg*hg*inv_ne
   ni = nio - h2
   mu(1) = mu(1) + hi/(hi+hg) * ( log(hi/com(1) + eps) - chi(1, 1)*ti )
   mu(1) = mu(1) + hg/(hi+hg) * ( log(hg/com(2) + eps) )
   dsdnif(1) = dsdnif(1) - hi/(hi+hg) * log(hi/com(1) + eps)
   dsdnif(1) = dsdnif(1) - hg/(hi+hg) * log(hg/com(2) + eps)
   zeta(1) = hi / (hi + hg + h2 + eps)

   ! Derivatives w.r.T. f and T
   qa = ne + 4.0d0*hg*h2a
   qb = ha(1)*(ne - 2.0d0*h2)
   qd = 1.0d0/(qa + qb)
   qc = 2.0d0*h2*qd
   qaf = (nef - ne*ref )*qc
   qat = (net - ne*h2at)*qc
   qf = hg*qd
   qbf = dvf*qf
   qbt = (chi(1,1)*ti + dvt)*qf
   hgf = qaf - qb*qbf
   hgt = qat - qb*qbt
   hif = ha(1)*(qaf + qa*qbf)
   hit = ha(1)*(qat + qa*qbt)
   nef = nef + hif
   net = net + hit
   h2f = h2*ref  + inv_ne*(2.0d0*h2a*hg*hgf - h2*nef)
   h2t = h2*h2at + inv_ne*(2.0d0*h2a*hg*hgt - h2*net)

   ! Hydrogen contribution to entropy, internal energy
   sif = sif - va(1)*hif - h2bt*h2f
   sit = sit - va(1)*hit - h2bt*h2t + h2*(zh2t + zh2tt)

   ! Avoid overflow problems when HI, HG, H2 -> 0
   si = si - hi*dlog(dabs(hi/com(1)) + eps)              &
           - hg*dlog(dabs(hg/com(2)) + eps)              &
           - h2*(dlog(dabs(h2/zh2s) + eps) - zh2t)
   ui = ui + chi(1,1)*hi + 0.5*ch2(1)*(hi + hg)

   ! db is 1 amu * number of baryons/cm3; rho is the mass density in g/cm3
   db = inv_ne*de
   dl = dlog(db)
   rho = db*avm
   arho = dlog(rho)
   rt = (re*ret + rp*rpt)/(re + rp) - inv_ne*net
   rf = (re*ref + rp*rpf)/(re + rp) - inv_ne*nef

   ! Second call to pressi compensates for spurious pressure and entropy terms
   de = db*neo
   !wmu = 1.0d0/(neo+ni+h2)
   wmu = 1.0d0/(ne+nio)
   call pressi ( 0, ti, psi, neo, nio, nzz, de, rf, rt, xtt, xft, xff, xftt, xfft, xfff, f, &
                  dc, dvt, dvf, dpb, dpbt, dpbf, dsb, dsbt, dsbf, dub )

   ! Pressure terms
   pe = cb*pe              ! Partial pressure due to electrons
   pp = cb*pp              ! Partial pressure due to positrons
   tcr = T*cr
   p0 = tcr*db
   pi = ni*p0
   t2 = T*T
   t4 = t2*t2
   pr = ca*t4/3.0d0
   b = 4.0d0*pr/p0
   ! When positrons are present, ignore pressure ionisation - everything is
   ! fully ionised anyway.
   if (pp > 0.0d0) then
      dpa = 0.0
      dpb = 0.0
      dpaf = 0.0
      dpbf = 0.0
      dpat = 0.0
      dpbt = 0.0
      dsaf = 0.0
      dsbf = 0.0
      dsat = 0.0
      dsbt = 0.0
      dsa = 0.0
      dsb = 0.0
   end if

   ! Finalise calculation of chemical potentials
   mu = mu + log(cen) + dl - 1.5d0*at

   ! Derivative of the density with respect to the composition
   vm     =  pi + (pe*pef + tcr*(dpaf-dpbf)) / ref
   dlrdni = (p0 + (pe*pef + tcr*(dpaf-dpbf)) / ref * zeta / ne) / vm

   ! Finalise calculation of entropy derivatives
   ! Common terms are collected in dS/dNi_f, then we add
   ! the terms proportional to dNe/dNi to dS/dNi_f and terms
   ! proportional to dlnrho/dNi_T,p and dNe/dNi to dS/dNi_p
   dsdnif = CR*(dsdnif + (1.5d0*at - dl + 1.5d0 - dlog(cen))) / avm * wmu
   dsdnip = dsdnif
   dsdnif = dsdnif + CR*((b + na + ne*se + ne*dsa - neo*dsb)*zeta/ne) / (cbn*avm)
   dsdnip = dsdnip + (CR*(b + na)/(ne*cbn) - (ne*sef + ne*dsaf-neo*dsbf)/ref) * dlrdni/avm &
                   + (CR*ne*se/avm + (ne*sef + ne*dsaf-neo*dsbf)/ref * (pi/vm - 1.0d0) ) * zeta/ne/avm

   ! Convert from derivative with respect to number fraction per unit mass to derivative of mass fraction
   dsdnif = dsdnif*cbn*wmu

   ! Pressure, in dyne/cm^2
   pg = pe + pp + pi + tcr*(dpa - dpb)
   p = max(pg + pr, 1.0d-200)
   pf = (pe*pef + pp*ppf      + pi*rf - h2f*p0 + tcr*(dpaf - dpbf))/p
   pt = (pe*pet + pp*ppt + pi + pi*rt - h2t*p0 + tcr*(dpat - dpbt)+pr*4.d0)/p
   ap = dlog(p)

   ! Entropy, in erg/g/K:
   dsf = nef*dsa + ne*dsaf - neo*dsbf - rf*b
   dst = net*dsa + ne*dsat - neo*dsbt - (rt - 3.0d0)*b
   sf = cr*(-ni*rf           + nef*se + ne*sef + np*spf + sif + dsf)/avm
   st = cr*( ni*(1.5d0 - rt) + net*se + ne*set + np*spt + sit + dst)/avm
   s = cr*(se*ne + sp*np + dsa*ne - dsb*neo + b + (1.5d0*at - dl + 2.5d0  &
        - dlog(cen))*ni + si)/avm

   ! Internal energy, in erg/g:
   u = tcr*(ue*ne + up*np + dua*ne - dub*neo + 1.5d0*ni + zh2t*h2   &
        + 0.75d0*b + ti*ui)/avm
   uf = T*sf + p/rho * rf + CR*T*( dot_product(mu, nif) - dv*nef )
   ut = T*st + p/rho * rt + CR*T*( dot_product(mu, nit) - dv*net )

   ! Other thermodynamic quantities:
   q = min(pt*sf - pf*st, 1.0d300)
   scp = -q/pf             ! Specific heat at constant pressure
   scv = st - sf*rt/rf     ! Specific heat at constant volume
   grada = sf/q
   gamma = q/(rt*sf-rf*st)
   zt = dsqrt(dabs((ne*ref/wf + nzz)/ni))

   ! Coefficients of thermal expansion (K&W (6.6))
   delta = rf*pt/pf - rt
   phi   = -dot_product(xa*dlrdni, -(ne+nio)/(1.0d0+zeta))

   ! tst ought to be zero, if all the above programming is correct
   tst = sf/cr - p*(rt*pf-rf*pt)/(tcr*rho)

   ! Copy abundance results to output struct
   abund%xa = xa
   abund%na = na
   abund%neo = neo
   abund%nio = nio
   abund%nzz = nzz
   abund%avm = avm
   abund%ne = ne
   abund%xai = xai

   !*** End of thermodynamic calculation. beginning of table-based calculations
   frho = arho/cln
   tf = at/cln

   ! Cap values to allowed range for opacity tables
   frho = min(10.0d0, max(frho,-12.0d0))
   tf = min(9.3d0, max(tf,  3.0d0))
   fk = get_opacity(frho, tf, abund)
   xhi = 4.0d0*cl*pr/(fk*rho*rho*scp*T)   ! cm**2/s
   ! Plasma and photon viscosity
   ! The plasma viscosity is from Spitzer (1962)
   ! Photon viscosity from Thomas (1930) and Weinberg (1971) with the
   ! correction factor (10/9) of van den Horn & van Weert (1981).
   llambda = log( 2.0d0/(3.0d0*echar**3) * sqrt(AMU*wmu*(boltzm*T)**3/(cpi*rho*zt**5)) )
   mu_plasma = 0.406*sqrt( AMU*wmu*(boltzm*T)**5 ) / ( (zt*echar)**4*llambda )
   mu_rad = 8.0d0*pr/(3.0d0*cl*fk*rho)    ! gram/cm s
   nu = ( mu_plasma + mu_rad ) / rho
   ! Prandtl number
   prandtl = nu/xhi
   
   ! Derivatives of opacity: (dlogk/dlogT)_rho and (dlogk/dlogrho)_T
   !> \todo FIXME: these are the results for Kramer's opacity, should take
   !! derivative of spline instead to get value from tables.
   !<
   fkr =  1.0d0
   fkt = -3.5d0

   ! Neutrino-loss rates from Itoh et al (1983-1992)
   call get_neutrino_rate(tf, frho, abund, en)
   ! Neutrino rate, from Itoh et al. (1996)
   !      CALL GET_NEUTRINO_RATE(T, RHO, NIO, NEO, EN)

   ! Copy all EoS results to output struct
   eos%at = at
   eos%af = af
   eos%ap    = ap;    eos%arho    = arho; eos%u     = u;     eos%p    = p
   eos%rho   = rho;   eos%fk      = fk;   eos%T     = T;     eos%sf   = sf
   eos%st    = st;    eos%zt      = zt;   eos%grada = grada; eos%scp  = scp
   eos%rf    = rf;    eos%rt      = rt;   eos%xhi   = xhi;   eos%s    = s
   eos%pr    = pr;    eos%pg      = pg;   eos%pf    = pf;    eos%pt   = pt
   eos%en    = en;    eos%wmu     = wmu;  eos%gamma1 = gamma
   eos%delta = delta; eos%phi     = phi;  eos%fkt   = fkt
   eos%fkr   = fkr;   eos%nu      = nu;   eos%prandtl = prandtl
   eos%dv    = dv
   return

contains

   ! FXP: Fudged exponential function, used to avoid too large or too small numbers
   ! in the ionisation state of the Saha equation. Usesprecomputed values of
   ! the limiting exponentials.
   function fxp(x)
      use real_kind
      
      implicit none
      real(double) :: fxp
      real(double), intent(in) :: x
      real(double), parameter :: fxp_low = -50.0d0
      real(double), parameter :: fxp_high = 50.0d0
      real(double), parameter :: low = 1.928749847963917820563444d-22
      real(double), parameter :: high = 5.184705528587072045056000d+21
      if (x>fxp_high) then
         fxp = high
         return
      else if (x>fxp_low) then
         fxp = exp(x)
         return
      end if
      fxp = low
   end function fxp
end subroutine statef



! ------------------------------------------------------------------------------
! FIDRAC
! Calculate the Fermi-Dirac integrals, as parametrised by functions f and g.
! ------------------------------------------------------------------------------
!  Returns:
!   FDI(9) - array of Fermi-Dirac integrals.
! ------------------------------------------------------------------------------
subroutine fdirac ( f, g, fdi )
   use real_kind
   
   implicit none
   real(double), intent(in) :: f,g
   real(double), intent(out) :: fdi(16)

   integer :: i,ij,ik,il,im,ind
   real(double) :: ff(4), gg(4), vw(4,4), vx(4,4)
   real(double) :: vf,vg,uf,ug,fdf,wv,ww

   real(double) :: d(3,3)
   real(double) :: xtt, xft, xff, xttt, xftt, xfft, xfff

   real(double), parameter :: c(4,4,3) = reshape(                     &
        (/ 2.315472d0,  7.128660d0,  7.504998d0,  2.665350d0,   &
           7.837752d0, 23.507934d0, 23.311317d0,  7.987465d0,   &
           9.215560d0, 26.834068d0, 25.082745d0,  8.020509d0,   &
           3.693280d0, 10.333176d0,  9.168960d0,  2.668248d0,   &
           2.315472d0,  6.748104d0,  6.564912d0,  2.132280d0,   &
           7.837752d0, 21.439740d0, 19.080088d0,  5.478100d0,   &
           9.215560d0, 23.551504d0, 19.015888d0,  4.679944d0,   &
           3.693280d0,  8.859868d0,  6.500712d0,  1.334124d0,   &
           1.157736d0,  3.770676d0,  4.015224d0,  1.402284d0,   &
           8.283420d0, 26.184486d0, 28.211372d0, 10.310306d0,   &
           14.755480d0, 45.031658d0, 46.909420d0, 16.633242d0,  &
           7.386560d0, 22.159680d0, 22.438048d0,  7.664928d0/), &
           (/ 4, 4, 3 /) )


   ! Evaluate Fermi-Dirac integrals (Eggleton, Faulkner & Flannery 1973).
   ! Matrix D contains rho*, P* and Q* and 1st logarithmic derivatives w.r.t.
   ! T, f. XTT etc are 2nd and 3rd log derivs of rho*
   vf = 1.0d0/(1.0d0 + f)
   vg = 1.0d0/(1.0d0 + g)
   uf = f*vf
   ug = g*vg
   fdf = g*g + g
   fdf = uf*fdf*dsqrt(fdf)
   ff(1) = 1.0d0
   gg(1) = 1.0d0
   xtt = 0.0d0
   xft = 0.0d0
   xff = 0.0d0
   xttt = 0.0d0
   xftt = 0.0d0
   xfft = 0.0d0
   xfff = 0.0d0
   do i = 2, 4
      ff(i) = f*ff(i - 1)
      gg(i) = g*gg(i - 1)
      fdf = fdf*vf*vg
   end do
   ind = 4
   do i = 1, 3
      vx(1:ind, 1:ind) = 0.0d0
      do ij = 1, 4
         do ik = 1, 4
            vw(1, 1) = c(ik, ij, i)*gg(ij)*ff(ik)
            do il = 1, ind - 1
               vw(il + 1, 1) = (ij - 1)*vw(il, 1)
               do im = 1, ind - il
                  vw(il, im + 1) = (ik - 1)*vw(il, im)
               end do
            end do
            vx(1:ind, 1:ind) = vx(1:ind, 1:ind) + vw(1:ind, 1:ind)
         end do
      end do

      wv = 1.0d0/vx(1, 1)
      vx(1, 2:ind) = vx(1, 2:ind)*wv
      vx(2:ind, 1:ind-1) = vx(2:ind, 1:ind-1)*wv
      d(i, 1) = fdf*vx(1, 1)
      d(i, 2) = vx(2, 1) + 1.5d0 - 1.5d0*ug
      ww = 0.5d0*d(i, 2) - 4.0d0
      d(i, 3) = vx(1, 2) + 1.0d0 + ww*uf

      if ( i == 1 ) then
         ! Second- and third-order derivatives of density, needed in pressi:
         xtt = vx(3, 1) - vx(2, 1)**2 - 1.5d0*ug*vg
         xft = vx(2, 2) - vx(2, 1)*vx(1, 2) + 0.5d0*uf*xtt
         xff = vx(1, 3) - vx(1, 2)**2 + uf*(xft + ww*vf - 0.25d0*xtt*uf)
         xttt = vx(4, 1) + vx(2, 1)*(2.0d0*vx(2, 1)**2 - 3.0d0*vx(3, 1))   &
              - 1.5d0*(1.0-g)*ug*vg**2
         xftt = vx(3, 2) - 2.0d0*vx(2, 1)*vx(2, 2) - vx(1, 2)*(vx(3, 1)    &
              - 2.0d0*vx(2, 1)**2) + 0.5d0*uf*xttt
         xfft = vx(2, 3) - 2.0d0*vx(2, 2)*vx(1, 2) - vx(2, 1)*(vx(1, 3)    &
              - 2.0d0*vx(1,2)**2) + uf*(xftt + 0.5d0*vf*xtt - 0.25d0*uf*xttt)
         xfff = vx(1, 4) + vx(1, 2)*(2.0d0*vx(1, 2)**2 - 3.0d0*vx(1, 3))   &
              + uf*(1.5d0*(xfft + vf*xft) - uf*(0.75d0*(xftt + vf*xtt) -   &
              0.125d0*uf*xttt) + ww*(1.0d0 - f)*vf**2)
      end if
      ind = 2
   end do

   ! Copy output
   fdi = (/ d(1, 1), d(2, 1), d(3, 1), d(1, 2), d(2, 2), d(3, 2),          &
            d(1, 3), d(2, 3), d(3, 3), xtt, xft, xff, xttt, xftt, xfft, xfff /)
end subroutine fdirac



! --------------------------------------------------------------------------
! PRESSI
! --------------------------------------------------------------------------
! Non-ideal corrections to the equation-of-state: pressure ionisation and
! Coulomb interactions.
! As explained in Pols&al. (1995), the effect of pressure ionisation is
! calculated twice: once for the actual number of electrons Ne and once for
! the total number of electrons Ne0, to make sure the corrections to
! pressure and entropy vanish in the limit of complete ionisation.
!
! Input parameters:
!  IPR      - 1 or 0, if 1 also calculates Coloumb interactions
!  TI       - 1eV/kT, reciprocal temperature
!  PSI      - electron degeneracy parameter
!  NEO      - electron number fraction
!  NIO      - ion number fraction
!  NZZ      - average charge squared, ~electron density * average ion charge
!  XI       - electron number density in funny units, 1 amu * electrons/cm3
!  XF       - logarithmic derivative of electron density to f, dln rhoe/dln f
!  XT       - logarithmic derivative of electron density to T, dln rhoe/dln T
!  xtt      - higher order logarithmic derivative of electron density
!  xft      - higher order logarithmic derivative of electron density
!  xff      - higher order logarithmic derivative of electron density
! (xttt     - higher order logarithmic derivative of electron density)
!  xftt     - higher order logarithmic derivative of electron density
!  xfft     - higher order logarithmic derivative of electron density
!  xfff     - higher order logarithmic derivative of electron density
!  F        - Eggleton F parameter (alternative degeneracy parameter)
! Output parameters:
!  DC       - Correction for electron chemical potential (sign flipped)
!  DCT      - Derivative of above with respect to log T
!  DCF      - Derivative of above with respect to log f
!  DP       - Correction for the pressure
!  DPT      - Derivative of above with respect to log T
!  DPF      - Derivative of above with respect to log f
!  DS       - Correction for the entropy
!  DST      - Derivative of above with respect to log T
!  DSF      - Derivative of above with respect to log f
!  DU       - Correction for the internal energy
! The output parameters are calculated by taking derivatives of the free
! energy, which requires high order derivatives of the electron density.
! --------------------------------------------------------------------------
subroutine pressi(ipr, ti, psi, neo, nio, nzz, xi, xf, xt, xtt, xft, xff, xftt, xfft, xfff, f, &
     dc, dct, dcf, dp, dpt, dpf, ds, dst, dsf, du)
   use real_kind
   use constants
   ! Computes effect (`pressure ionisation') on EoS of a free energy contribution
   ! delta(f) = -R.T.ne.g(x,y), where x = xi = Ne/v = ne, and y = yi = xih/(R.T).
   ! also does Coulomb interaction.
   ! delta(el. chem. pot.) = -dc = 1/(R.T) df/dNe
   ! delta(p) = R.T.dp = -dF/dv
   ! delta(s) = R.ne.ds = -dF/dT.
   ! Works firstly with independent variables x, y (called xi, yi), which are
   ! effectively Ne, v and T, but then has to get derivatives of dc, dp, ds
   ! with respect to independent variables f and T.

   implicit none
   integer, intent(in) :: ipr
   real(double), intent(in) :: ti, psi, neo, nio, nzz
   real(double), intent(in) :: xi, xf, xt, f, xtt, xft, xff, xftt, xfft, xfff
   real(double), intent(out) :: dc, dct, dcf, dp, dpt, dpf, ds, dst, dsf, du

   real(double) :: yi,ee,wx,cc,bb,gi,ff1,aa,th,ff2,af,thf,tht,rxf,thx,thy,psx,psy
   real(double) :: w2,psxx,psxy,psyy,eex,bx,by,bxx,bxy,byy,cx,cxx
   real(double) :: dgdx,dgdy,dgdxx,dgdxy,dgdyy,aff,thff,thft,thtt,tof
   real(double) :: wxx,wxy,thxx,thxy,thyy,thc,ww,ge,eeg,beg,begg,geg,dgdg
   real(double) :: gam,rbe,dgdgg,rthc,wy,gamx,gamy,gamxx,gamxy,gamyy,wt

   real(double), parameter :: ca1 = 0.89752d0
   real(double), parameter :: ca2 = 0.768d0
   real(double), parameter :: ca3 = 0.208d0
   real(double), parameter :: cp1 = 3.0d0
   real(double), parameter :: cp2 = 0.35d0
   real(double), parameter :: cp3 = 2.0d0
   real(double), parameter :: cp4 = 3.0d-2

   real(double) :: cbrt
   

   ! Pressure ionization:
   yi = 13.595d0*ti
   ee = 1.0d0/(1.0d0 + xi/cp4)
   wx = min((cp1/xi)**cp2, 300.0d0)
   cc = dexp(-wx)
   bb = yi + psi - cp3*dlog(ee)

   ! extra free energy is -R.T.Ne.GI(X,Y)
   ! GI = exp{-WX(X)}{Y + psi(X,Y) + 2ln(1+X/const)}
   ! It OUGHT to have 3psi/5 for psi, but works worse!
   gi = cc*bb
   ff1 = 1.0d0/(f + 1.0d0)
   ! AA = dlnf/dpsi; THeta = dlnne/dpsi; also needed in Coulomb corrections
   ! f, T derivs, then X, Y.
   aa = dsqrt(ff1)
   th = aa*xf
   ff2 = -0.5d0*f*ff1
   af = aa*ff2
   thf = af*xf + aa*xff
   tht = aa*xft
   rxf = 1.0d0/xf
   thx = thf*rxf
   thy = thx*xt - tht
   ! first and second derivatives dPSI/dlnX ... d2PSI/dlnY2
   psx = 1.0d0/th
   psy = xt*psx
   w2 = psx*psx
   psxx = -thx*w2
   psxy = -thy*w2
   psyy = psxy*xt + (xft*xt*rxf - xtt)*psx
   ! derivative -dlnEE/dlnX; -d2lnEE/dlnX2 = -EE*dlnEE/dlnX
   eex = cp3*(xi/cp4)*ee
   ! derivatives of BB
   bx = psx + eex
   by = yi + psy
   bxx = psxx + ee*eex
   bxy = psxy
   byy = yi + psyy
   ! derivatives of CC
   cx = cp2*wx*cc
   cxx = cp2*cx*(wx - 1.0d0)
   ! derivatives of GI
   dgdx = cc*bx + cx*bb
   dgdy = cc*by
   dgdxx = cc*bxx + 2.0d0*cx*bx + cxx*bb
   dgdxy = cc*bxy + cx*by
   dgdyy = cc*byy
   if ( ipr == 1 ) then
      ! Coulomb interaction.
      ! further derivatives of AA, THeta
      aff = 3.0d0*af*ff2 + af
      thff = aff*xf + 2.0d0*af*xff + aa*xfff
      thft = af*xft + aa*xfft
      thtt = aa*xftt
      ! d2TH/dlnX2, etc
      tof = xt*rxf
      wxx = thff - thx*xff
      wxy = thx*xft - thft
      thxx = wxx*rxf*rxf
      thxy = wxy*rxf + thxx*xt
      thyy = tof*(tof*wxx + 2.0d0*wxy) + thtt - thx*xtt
      ! gam is the plasma interaction parameter. note that thc = zt**2*nio/neo
      thc = th + dabs(nzz/neo)
      ww = 1.0d0/ca2
      gam = cbrt(xi*(cpl*neo/nio)**2*c3rd)*ti/cevb*thc
      ! new BB and EE, and their derivates
      bb = (ca1*dsqrt(3.0d0/gam))**ww
      ee = (gam/(gam + ca3))**ww
      rbe = 1.0d0/(ee + bb)
      ! further addition GE to free en; adds to previous GI.
      ! GE = GE(GAMma), GAM = const. * X**0.33 * Y * (THeta(X,Y) + const)
      ge = (nio/neo)*ca1*gam*rbe**ca2
      eeg = ca3/(gam + ca3)
      beg = (eeg*ee - 0.5d0*bb)*rbe
      begg = (eeg*eeg*(1.0d0 - gam*ca2/ca3)*ee + 0.25d0*bb)*rbe
      geg = 1.0d0 - beg
      dgdg = ge*geg
      dgdgg = ge*(geg*geg + (beg*beg - begg)*ww)
      rthc = 1.0d0/thc
      wx = thx*rthc
      wy = thy*rthc
      ! dlnGAM/dlnX, etc
      gamx = c3rd + wx
      gamy = 1.0d0 + wy
      gamxx = thxx*rthc - wx*wx
      gamxy = thxy*rthc - wx*wy
      gamyy = thyy*rthc - wy*wy
      gi = gi + ge
      ! derivatives w.r.t. X, Y; in effect w.r.t. Ne, V, and T, since X = Ne/V
      dgdx = dgdx + dgdg*gamx
      dgdy = dgdy + dgdg*gamy
      dgdxx = dgdxx + dgdgg*gamx**2 + dgdg*gamxx
      dgdxy = dgdxy + dgdgg*gamx*gamy + dgdg*gamxy
      dgdyy = dgdyy + dgdgg*gamy**2 + dgdg*gamyy
   end if
   ! Evaluate changes to el. chem. potential (-dc), pressure (R.T.dp), and
   ! entropy (R.ne.ds), and their derivatives w.r.t. log(f) and log(T)
   dc = dgdx + gi
   ww = dgdxx + dgdx
   wt = ww*xt
   dcf = ww*xf
   dct = wt - dgdxy - dgdy
   dp = -xi*dgdx
   dpf = -xi*dcf
   dpt = xi*(dgdxy - wt) + dp
   ds = gi - dgdy
   dst = dgdx - dgdxy
   dsf = dst*xf
   dst = dst*xt - dgdy + dgdyy
   du = -dgdy
end subroutine pressi



! Compute rates of (at present) 20 nuclear rections, and the corresponding
! energy and neutrino release
subroutine nucrat ( tl, abund, eos )
   use real_kind
   use constants
   use reaction_rate_data
   use eostate_types
   
   implicit none
   real(double), intent(in) :: tl
   type(abundance), intent(in) :: abund
   type(eostate), intent(inout) :: eos
   integer :: it
   real(double) :: rhb,wc,wb,wa,xb,vl,za,zb,zc,zd,tf,tt,tu
   real(double) :: scrn(20),strn(20),dstr(20)
   real(double) :: fpng,rpn,fccg,rcc,f34,pp2,pp3,qpp,qnpp

   real(double) :: rtt(21)
   real(double), parameter :: csa = 0.624
   real(double), parameter :: csb = 0.316
   real(double), parameter :: csc = 0.460
   real(double), parameter :: csd = 0.38
   real(double), parameter :: cxd = 0.86

   real(double) :: rrt(21)
   real(double) :: rpp, r33, r34,  rbe, rbp, rpc, rpna,  rpo,  r3a,  rac
   real(double) :: ran, rao, rane, rcca, rco, roo, rgne, rgmg, rccg, rpng

   real(double) :: n1, n4, n12, n14, n16, n20, n24, n28, n56
   real(double) :: ne, ni, avm

   real(double) :: cbrt
   
   ! Not always defined:
   dstr = 0.0_dbl
   strn = 0.0_dbl

   ni = abund%nio
   ne = abund%neo
   n1 = abund%na(1)
   n4 = abund%na(2)
   n12 = abund%na(3)
   n14 = abund%na(4)
   n16 = abund%na(5)
   n20 = abund%na(6)
   n24 = abund%na(7)
   n28 = abund%na(8)
   n56 = abund%na(9)
   avm = abund%avm

   ! rhb is 'baryon density': 1 amu * number of baryons per cm3
   rhb = eos%rho/avm
   ! Electron-screening theory from Graboske, Dewitt, Grossman & Cooper (1973),
   ! for strong (za, zb, zc) are intermediate screening (zd). The reaction
   ! dependent charge parameters are stored in cza ... czd.
   wc = dot_product((/n1, n4, n12, n14, n16, n20, n24, n28, n56/), vz(1:9))
   wc = wc/ni
   wb = ne/ni
   wa = eos%zt*eos%zt/(wb*wb)
   xb = cbrt(dabs(wb))
   vl = cpl*dsqrt(dabs(ni)*rhb*dexp(-3.0d0*tl))
   za = csa*xb*vl**(2.0/3.0)
   zb = csb*xb*za
   zc = csc/(xb*xb)
   zd = csd*wc*wa*(vl/(wa*eos%zt))**cxd

   ! Reaction rates interpolated in T, mostly from Caughlan & Fowler (1988):
   tf = tl/cln
   rrt(:) = 0.0d0
   rtt(:) = 0.0d0
   tt = 50.0d0*(tf - 6.0d0) + 1.0d0
   if ( tt >= 1.0d0 ) then
      it = max0(1, min0(199, int(tt)))
      tt = tt - it
      tu = 1.0d0 - tt
      rrt(2:20) = tu*crt(it, 1:19) + tt*crt(it + 1, 1:19)
      scrn(1:19) = zd*czd(1:19)
      strn(1:19) = za*cza(1:19) + zb*czb(1:19)
      dstr(1:19) = zc*czc(1:19)
      where (dstr < 0.29d0*strn) scrn = min(scrn, strn - dstr)
      rrt(2:20) = exp(cln*rrt(2:20) + scrn(1:19))

      ! Derivative (dR/dlogT)_rho, by differentiation of the interpolation function:
      rtt(2:20) = 50.0/cln * (crt(it + 1, 1:19) - crt(it, 1:19))
   end if

   ! Copy rates to sensible variable names
   !> \todo FIXME: we really shouldn't need this, do we?
   rpp = rrt(2)
   r33 = rrt(3)
   r34 = rrt(4)
   rBe = rrt(5)
   rBp = rrt(6)
   rpc = rrt(7)
   rpN = rrt(8)
   rpO = rrt(9)
   r3a = rrt(10)
   raC = rrt(11)
   raN = rrt(12)
   raO = rrt(13)
   raNe = rrt(14)
   rCC = rrt(15)
   rCO = rrt(16)
   rOO = rrt(17)
   rgNe = rrt(18)
   rgMg = rrt(19)
   rCCg = rrt(20)
   
   ! Multiply with density and abundances to get rates per baryon per second,
   ! note that abundances of He3 and Be7 are not needed in equilibrium
   rpp = rhb*n1*n1*rpp/2.0d0
   r33 = rhb*r33/2.0d0
   r34 = rhb*n4*r34
   rbe = rhb*ne*rbe
   rbp = rhb*n1*rbp
   rpc = rhb*n1*n12*rpc
   rpn = rhb*n1*n14*rpn
   rpo = rhb*n1*n16*rpo
   r3a = rhb*rhb*n4*n4*n4*r3a/6.0d0
   rac = rhb*n4*n12*rac
   ran = rhb*n4*n14*ran
   rao = rhb*n4*n16*rao
   rane = rhb*n4*n20*rane
   rcc = rhb*n12*n12*rcc/2.0d0
   rco = rhb*n12*n16*rco
   roo = rhb*n16*n16*roo/2.0d0
   rgne = n20*rgne
   rgmg = n24*rgmg
   ! Branching of pn and cc reactions
   fpng = 8.0d-4
   rpna = (1.0d0 - fpng)*rpn
   rpng = fpng*rpn
   rpn = rpna
   fccg = rccg
   rcca = (1.0d0 - fccg)*rcc
   rccg = fccg*rcc
   rcc = rcca
   ! pp chain in equilibrium, rpp becomes effective rate of 2 h1 -> 0.5 he4
   f34 = 0.0d0
   rpp = min(rpp, 1.0d303)
   r33 = min(r33, 1.0d303)
   r34 = min(r34, 1.0d303)
   if (r34 > 1.0d-20)  f34 = 2.0d0/(1.0d0 + dsqrt(1.0d0 + 8.0d0*rpp*r33/(r34*r34)))
   rpp = rpp*(1.0d0 + f34)
   pp2 = 1.0d0
   if (rbe + rbp > 1.0d-20) pp2 = rbe/(rbe + rbp)
   pp3 = 1.0d0 - pp2
   qpp = qrt(1) + 0.5d0*qrt(2)
   qnpp = (qnt(1) + f34*(qnt(4)*pp2 + qnt(5)*pp3))/(1.0d0 + f34)

   ! Copy rates back to rrt array for loop below
   rrt(2:21) = (/ rpp, r33, r34,  rBe, rBp, rpc, rpN,  rpO,  r3a,  raC,&
                  raN, raO, raNe, rCC, rCO, rOO, rgNe, rgMg, rCCg, rpng /)

   ! Calculate energy release and neutrino loss, in erg/gram/sec
   eos%ex = qpp*rpp
   eos%enx = qnpp*rpp
   eos%ext = qpp*rpp*rtt(2)
   eos%ex = eos%ex + dot_product(qrt(6:20), rrt(7:21))
   eos%enx = eos%enx + dot_product(qnt(6:20), rrt(7:21))
   eos%ext = eos%ext + dot_product(qrt(6:20), rrt(7:21)**2)
   eos%ext = eos%ext / max(eos%ex, 1.0d-50)
   ! Guard against overflows
   eos%ex = min(cme*eos%ex/avm, 1.0d303)
   eos%enx = -min(cme*eos%enx/avm, 1.0d303)

   ! Copy reaction rates to output struct
   eos%rpp  = rpp;  eos%r33  = r33;
   eos%r34  = r34;  eos%rbe  = rbe;
   eos%rbp  = rbp;  eos%rpc  = rpc;
   eos%rpna = rpna; eos%rpo  = rpo;
   eos%r3a  = r3a;  eos%rac  = rac;
   eos%ran  = ran;  eos%rao  = rao;
   eos%rane = rane; eos%rcca = rcca;
   eos%rco  = rco;  eos%roo  = roo;
   eos%rgne = rgne; eos%rgmg = rgmg;
   eos%rccg = rccg; eos%rpng = rpng;
end subroutine nucrat



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



subroutine load_reaction_neutrino_rates(file)
   use real_kind
   use reaction_rate_data
   use neutrinos
   
   implicit none
   integer, intent(in) :: file

   ! Read nuclear reaction and neutrino loss rate data
   read (file, '(1x, 10f7.3)') cnu     ! Neutrino rates
   read (file, '(1x, 10f7.3)') crt     ! Nuclear reaction rates
   rewind (file)
end subroutine load_reaction_neutrino_rates


