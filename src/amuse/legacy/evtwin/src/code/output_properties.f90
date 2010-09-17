! For some variables it is inconvenient to calculate them implicitly, so
! we calculate them explicitly, between timesteps/iterations
module explicit_functions
   use real_kind
   use mesh
   
   implicit none
   integer, parameter :: explv_gradmu = 1
   integer, parameter :: explv_avmu   = 2
   integer, parameter :: explv_logp   = 3
   integer, parameter :: num_explv = 4

   real(double) :: expl_var(NM, num_explv, 2)
   real(double), save :: radacc(9, nm, 2)
end module explicit_functions



module printb_global_variables
   ! Module for exchanging variables between subroutines within PRINTB only
   use real_kind
   
   implicit none
   real(double) ::  dmtr, perc, rcz, drcz, dmt, tet
   real(double) :: dmsw, oms, fr, fl, sdc, sdm, sds, stc, stm, sts, psurf, pcntr, vk2
   real(double) :: raf, dphi, ephi, f1, df, htot, aj, dl, tkh, mex(12), ecc

   ! Meshpoint were maximum temperature is reached
   integer :: im_tmax

   ! Pressures at the edge of the H and He exhausted cores
   real(double) :: ph
   real(double) :: phe
end module printb_global_variables





subroutine printb ( jo, jcm, rz, Jstar, ift )
   use real_kind
   use mesh
   use mesh_enc
   use control
   use constants
   use settings
   use test_variables
   use current_model_properties
   
   implicit none
   ! Subroutine Arguments
   real(double), intent(out) :: rz
   integer, intent (in) :: Jstar, ift
   integer, intent (inout) :: jo, jcm

   integer :: ig, jmad
   integer, save :: lmod(2)

   ! ----------------------------------------------------------
   ! Compute quantities
   ! ----------------------------------------------------------
   call compute_output_quantities ( Jstar )
   rz = rlf(Jstar)

   ! -------------------------------------------------------------------
   ! Update quantities that have to be calculated explicitly
   ! -------------------------------------------------------------------
   call update_explicit_quantities( jo, Jstar )

   ! -------------------------------------------------------------------
   ! Update the control parameters for the next timestep
   ! -------------------------------------------------------------------
   call update_timestep_parameters( Jstar )

   ! -------------------------------------------------------------------
   ! Adjust the weights of the mesh spacing function
   ! -------------------------------------------------------------------
   call adjust_mesh_spacing

   ! -------------------------------------------------------------------
   ! Check termination conditions
   ! -------------------------------------------------------------------
   call check_stop_conditions ( Jstar, jo, jcm, ift )

   ! -------------------------------------------------------------------
   ! Output
   ! -------------------------------------------------------------------
   ig = max(Jstar, jb)
   if ( ift == 24 ) ig = 25

   jmad = jmod + 1           ! Model number for output, count from 1
   if ( jnn == 0 ) jmad = -jmod ! Backup, print negative model number

   ! Decide whether to print the interior or not, 1,2: file.out1,2
   if ( mod(jmad, kt1) == 0 .or. jmod == 0 ) call write_internal_details(ig)

   ! Print the short summary, for every KT4'th model, 1,2: file.out1,2
   if ( kt4>0 .and. mod(jmod, kt4) == 0 ) call write_summary ( Jstar, jmad, ig )

   ! Same, for nucleosynthesis
   if ( kt4>0 .and. mod(jmod, kt4) == 0 ) call write_nucsyn_summary ( Jstar, jmad, 34 + ig )

   ! Funny composition profile - flag this as a problem
   !if ( mh < mhe .or. mhe < mco ) write (ig,*) 'ccc', mh, mhe, mco
   call flush ( ig )

   ! Write some quantities of each model for plotting purposes, 31,32: file.plt1,2
   ! Print a line for every KT4'th model
   ! Don't output models during ZAHB construction
   if (ig <= 2 .and. mod(jmad, kt4) == 0 .and. jmad > 0) then
      ! Overwrite last line if model didn't converge:
      if (kt4>0) then
         do while (jmad <= lmod(ig))
            backspace(30+ig)
            lmod(ig) = lmod(ig) - kt4
         end do
      end if
      lmod(ig) = jmad
      call write_plt_file(Jstar, jmad, 30+ig)
      call write_nucsyn_plt_file(Jstar, jmad, 36+ig)
   end if

   ! Write detailed model data for plotting purposes, 33,34: file.mdl1,2
   ! Write one model every KT1'th model
   ! Don't output models during ZAHB construction
   !> \todo FIXME: in case of convergence failure, go back and delete previous
   !! model, as done for the PLT file above
   !<
   if (ig <= 2 .and. mod(jmad, kt1) == 0) call write_mdl_file(Jstar, jmad, 32+ig)
   if (ig <= 2 .and. mod(jmad, kt1) == 0) call write_nucsyn_mdl_file(Jstar, jmad, 38+ig)
   return
end subroutine printb



subroutine compute_output_quantities ( Jstar )
   use real_kind
   use mesh
   use mesh_enc
   use control
   use constants
   use settings
   use plotvariables
   use printb_global_variables
   use solver_global_variables, only: default_er
   use reaction_rate_data, only: qrt
   use atomic_data
   use test_variables
   use eostate_types
   use funcs1_interface
   use current_model_properties
   use structure_variables
   
   implicit none
   ! Subroutine Arguments
   integer, intent (in) :: Jstar
   
   ! /abund /
   real(double) :: ya(9), wne0, avm, wne
   
   real(double) ::       ap, arho, u, p, rho, fk, t, sf, st, zt, grada,    &
        scp, rf, rt, xhi, s, pr, pg, pf, pt, en,                     &
        rpp, r33, r34,  rbe, rbp, rpc, rpna,  rpo,  r3a,  rac,       &
        ran, rao, rane, rcca, rco, roo, rgne, rgmg, rccg, rpng,      &
        ex, enx, wmu, delta, phi, ext, fkt, fkr, prandtl
   type(eostate) :: eos
   type(abundance) :: abund

   integer :: icnv, icnv1

   integer :: ic, is, io, je, i, ikk, ik, imax

   real(double), external :: pstv

   real(double) :: mcnv, mcnv1
   !real(double) :: ftout(NM), fpout(NM), xr(NM), xm(NM)
   real(double) :: exlim, tmax, hpc, wf, pdg, peg, rich, domega2
   real(double) :: cgs_sgf, dmuleft, dmuright, dmu, dwleft
   real(double) :: dwright, dw, dcgs_con, dcgs_thl, dcgs_dsi, dcgs_ssi, dcgs_shi, dcgs_gsf, dcgs_esc, ton, ris
   !real(double) :: es_hp, es_vmu, es_ves
   real(double) :: drmax, ddr, dty
   real(double) :: var(nvar), dvar(nvar), fn1(nfunc)
   real(double) :: qpp, q33, q34, qbe, qbp, qpc, qpna, qpo, q3a, qac,qan, qao, qane, qcca, qco, qoo, qgne, qgmg, qccg, qpng


   ic = 1
   is = 1
   io = 1
   im_tmax = 1
   ! Set lower limit for nuclear burning to 10*L/M, for printing purposes only (Onno)
   exlim = 10.0*h(8 + 24*(Jstar-1),1)/h(4 + 24*(Jstar-1),1)
   je = 1
   mex(1:12) = 0.0d0
   tmax = 0.0d0
   px = 0.0d0
   wmh= 0.0d0
   wmhe= 0.0d0
   mh= 0.0d0
   mhe= 0.0d0
   mco= 0.0d0
   vmg= 0.0d0
   be(1:2) = 0.0d0
   lh = 0.0d0
   lhe = 0.0d0
   lc = 0.0d0
   lnu = 0.0d0
   lth = 0.0d0
   mcb(1:8) = 0.0d0
   msb(1:6) = 0.0d0
   rcb(1:8) = 0.0d0
   tct(1:8) = 0.0d0

   ! Simple names for reaction rate Q values
   qpp = qrt(1)
   q33 = qrt(2)
   q34 = qrt(3)
   qbe = qrt(4)
   qbp = qrt(5)
   qpc = qrt(6)
   qpna = qrt(7)
   qpo = qrt(8)
   q3a = qrt(9)
   qac = qrt(10)
   qan = qrt(11)
   qao = qrt(12)
   qane = qrt(13)
   qcca = qrt(14)
   qco = qrt(15)
   qoo = qrt(16)
   qgne = qrt(17)
   qgmg = qrt(18)
   qccg = qrt(19)
   qpng = qrt(20)

   dl = 0.0d0
   mcnv = 0.d0
   mcnv1 = 0.d0
   icnv = 0
   icnv1 = 0

   ! Rotational periods for the centre and surface
   psurf = 0.0d0
   pcntr = 0.0d0
   dty = dt/csy
   do ikk = 1, kh
      ik = kh + 1 - ikk

      do i = 1, 16
         dvar(i) = dh(i + 24*(Jstar-1), ik)
         var(i) = h(i + 24*(Jstar-1), ik) + dvar(i)
      end do
      do i = 1, nxfstar
         dvar(40+i) = dh(i + 40+24*(Jstar-1), ik)
         var(40+i) = h(i + 40+24*(Jstar-1), ik) + dvar(40+i)
      end do
      dvar(17:24) = dh(17:24, ik)
      var(17:24) = h(17:24, ik) + dvar(17:24)
      call funcs1 ( ik, -1, var(:), dvar(:), fn1(:), eos, abund )
      
      ! Assign EoS output to local variables (for convenience)
      !> \todo FIXME: below, should use the output from the struct directly to
      !! simpify the code
      !<
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
      delta = eos%delta; phi     = eos%phi;  ext   = eos%ext;   fkt  = eos%fkt
      fkr   = eos%fkr;   prandtl = eos%prandtl

      ya = abund%na
      wne0 = abund%neo
      avm = abund%avm
      wne = abund%ne

      if ( ikk == 1 ) then
         hpc = dsqrt(p/(cg*rho*rho))
         mc(Jstar) = 3.5d-33*rho*hpc**3
         !PCNTR = exp(var(13))
         !PCNTR = var(13)
         pcntr = 2.0*cpi/(dabs(var(13)) * csday)
      end if
      if ( ikk == kh ) then
         !PSURF = exp(var(13))
         !PSURF = var(13)
         psurf = 2.0*cpi/(dabs(var(13)) * csday)
         !print *, JMOD, var(42), var(4)
         !print *, PCNTR, PSURF
      end if
      
      ! Evaluate the functions to be printed
      wf = dsqrt(1.0d0 + dexp(var(1)))
      px(1) = dmin1(99.0d0, 2.0d0*(wf - dlog(wf + 1.0d0)) + var(1))
      px(2) = p
      px(3) = rho
      px(4) = t
      tmax = max(tmax, t)
      if ( t == tmax ) im_tmax = ikk + 1
      px(5) = fk
      px(6) = grada
      px(7) = grad
      px(8) = dg
      pdg = px(8)
      if ( ikk == kh .and. pdg > 0.0d0 ) pdg = 0.0d0
      px(39) = eg
      peg = px(39)
      if ( ikk == kh .and. peg > 0.0d0 ) peg = 0.0d0
      px(9) = var(4)/cmsn
      
      ! Correct abundances for actual atomic masses
      do i = 1, 7
         px(i + 9) = can(i)*ya(i)/avm
         if (px(i + 9)<1.0d-99) px(i + 9) = 0.0d0
      end do
      px(46) = can(8) * ya(8) / avm
      px(47) = can(9) * ya(9) / avm

      px(17) = r/crsn
      px(18) = var(8)/clsn
      px(19) = eth
      px(20) = ex
      px(21) = -en
      px(27) = u
      px(28) = s
      px(29) = LoLedd
      px(30) = wl
      px(31) = cr*rho*t/pg       ! Mean molecular weight
      px(32) = wt
      px(33) = wne
      px(34) = dvar(4)/(cmsn*dty)  ! dM/dt
      !PX(34) = WNE0
      px(35) = wcv
      px(36) = var(12)
      px(37) = var(14)
      px(38) = var(19)

      ! Copy TWIN variables
      px(40:45) = sx(40:45, ikk + 1)
      
      ! Reaction rates:
      px(50:54) = sx(50:54, ikk + 1)
      
      ! Thermohaline mixing:
      px(23) = sx(23, ikk + 1)
      px(31) = sx(31, ikk + 1)
      px(30) = sx(30, ikk + 1)
      px(7) = sx(7, ikk + 1)
      
      ! Differential rotation:
      px(59) = sx(59, ikk + 1)
      px(60) = sx(60, ikk + 1)

      dl = dl + 0.5d0*(px(40) + sx(40, ikk))/clsn
      
      ! Locate convective/radiative radius boundaries (DG = 0):
      if ( ikk > 1 ) then
         
         ! FIXME: we start the iterations with ic=1 and fill up tct(2-9) and rbc(2-8) here!
         if ( ic < 8 ) then                         ! ic = 1-7
            if ( sx(8, ikk)*pdg <= 0.0d0 ) then
               ic = ic + 1                          ! Start a new convective layer; ic is 2-8 here (never 1), so ic+1 below can be 9
               rcb(ic) = (sx(17, ikk)*pdg-px(17)*sx(8, ikk))/(pdg-sx(8, ikk))
               if ( pdg > 0.0d0 .and. ic < 8) tct(ic + 1) = (px(17) - rcb(ic))/wcv                 ! otherwise ic+1 may be 9 - FIXME
               if ( sx(8, ikk) > 0.0d0 ) tct(ic) = tct(ic) + (rcb(ic) - sx(17, ikk))/sx(35, ikk)
            else
               ! Integrate dr/(conv. vel.), for convective turnover time:
               if ( wcv > 0.0d0 ) tct(ic + 1) = tct(ic + 1) + &                                                      ! ic+1 may be 8
                    (px(17) - sx(17,ikk)) * (wcv + sx(35,ikk)) / ( wcv*wcv + sx(35,ikk)*(sx(35,ikk) + wcv) )
            end if
         end if
         
         ! locate convective/semiconvective mass boundaries (EG = CGR = 0.01):
         if ( (px(39)-cgr)*(sx(39,ikk)-cgr) <= 0.0d0 .and. is <= 6 ) then
            msb(is) = (sx(9,ikk)*(px(39) - cgr) - px(9)*(sx(39,ikk) - cgr))/abs(px(39) - sx(39,ikk))
            is = is + 1
         end if
         
         ! Switch the counting of MCNV on and off:
         if(peg > 0.d0.and.px(4) < 1.d7) mcnv1 = mcnv1 + px(9) - sx(9,ikk)
         if(peg > 0.d0.and.px(4) < 1.d7 .and. sx(39,ikk)*peg < 0.0d0) then
            mcnv = mcnv + sx(9,ikk) - (sx(9,ikk)*peg - px(9)*sx(39,ikk)) / (peg - sx(39,ikk))
            icnv = 1
            icnv1 = 1
         end if
         if(peg < 0.d0.and.px(4) < 1.d7 .and. sx(39,ikk)*peg < 0.0d0) then
            mcnv = mcnv + (sx(9,ikk)*peg - px(9)*sx(39,ikk)) / (peg - sx(39,ikk)) - sx(9,ikk)
            if(icnv == 0.and.px(4) < 1.d5) mcnv = mcnv + sx(9,ikk)
            icnv = 0
         end if
         if (icnv == 1) mcnv = mcnv + px(9) - sx(9,ikk)
         if (ikk == kh.and.icnv1 == 0.and.mcnv == 0.d0.and.mcnv1 > 0.d0) mcnv = mcnv1 + sx(9,2)
         
         ! locate overshoot mass boundaries (EG = 0):
         if ( io < 8 .and. sx(39, ikk)*peg <= 0.0d0 ) then
            mcb(io) = (sx(9,ikk)*peg-px(9)*sx(39,ikk))/abs(peg-sx(39,ikk))
            io = io + 1
         end if
         
         ! locate burning shell boundaries (XH = CXB, XHE = CXB, XC = 0.05):
         if ( px(10) > cxb .and. sx(10,ikk) < cxb ) then
            mh = (px(9)*(cxb-sx(10,ikk)) + sx(9,ikk)*(px(10)-cxb)) / (px(10)-sx(10,ikk))
            ph = (px(2)*(cxb-sx(10,ikk)) + sx(2,ikk)*(px(10)-cxb)) / (px(10)-sx(10,ikk))
         end if
         if ( px(11) > cxb .and. sx(11,ikk) < cxb ) then
            mhe = (px(9)*(cxb-sx(11,ikk)) + sx(9,ikk)*(px(11)-cxb)) / (px(11)-sx(11,ikk))
            phe = (px(2)*(cxb-sx(11,ikk)) + sx(2,ikk)*(px(11)-cxb)) / (px(11)-sx(11,ikk))
         end if
         if ( px(12) > 0.05d0 .and. sx(12,ikk) < 0.05d0 ) then
            mco = (px(9)*(0.05d0-sx(12, ikk)) + sx(9, ikk)*(px(12)-0.05d0)) / (px(12) - sx(12, ikk))
         end if
         mh = max(mh, 0.0d0)
         mhe = max(mhe, 0.0d0)
         mco = max(mco, 0.0d0)
         
         ! Locate boundaries of nuclear-burning regions, EX > EXLIM (Onno)
         if((px(20)-exlim)*(sx(20,ikk)-exlim) <= 0.0d0.and.je < 12.and.&
              ikk >= 2) then
            mex(je)= (sx(9,ikk)*(px(20)-exlim)-px(9)*(sx(20,ikk)-exlim))&
                 /(px(20)-sx(20,ikk))
            je = je  + 1
         end if
         ! some homology invariants
         wf = 1.0d0/dlog10(px(2)/sx(2,ikk))
         px(24) = wf*dlog10(px(3)/sx(3,ikk))
         px(25) = -wf*dlog10(px(17)/(dabs(sx(17,ikk))+1.0d-10))
         px(26) = -wf*dlog10(px(9)/(dabs(sx(9,ikk))+1.0d-10))
         ! Some integrated quantities, to be printed in the short summary
      else     ! IKK>1
         ph = px(2)
         phe = max(1.0d19, px(2))
      end if   ! IKK>1
      !print *, 'printb:', ph, phe, PX(2)
      px(22) = px(9) - sx(9,ikk)
      if ( ikk == 1 ) px(22) = px(9)
      wmh = wmh + px(22)*px(10)
      wmhe = wmhe + px(22)*px(11)
      if ( .not. (ikk == 1 .or. mh == 0.0d0))  then
         if ( mh > sx(9, ikk) ) then
            be(Jstar) = (px(9) - mh)*egr
         else
            be(Jstar) = be(Jstar) + px(22)*egr
         end if
      end if
      wf = px(22)*cmsn/clsn

      lh = lh + wf*cme*((qpp + 0.5*q33)*rpp + qpc*rpc + qpo*rpo&
           + qpna*rpna + qpng*rpng)
      lhe = lhe + wf*cme*(q3a*r3a + qac*rac + qan*ran&
           + qao*rao + qane*rane)
      lc = lc + wf*cme*(qcca*rcca + qccg*rccg + qco*rco +qoo*roo&
           + qgne*rgne + qgmg*rgmg)
      lnu = lnu - wf*en
      lth = lth + wf*eth
      ! save places 40 - 45 for TWIN variables
      sx(1:39, ikk + 1) = px(1:39)
      sx(46:NPX, ikk + 1) = px(46:NPX)
   end do
   ! set the luminosity scale for the solver
   default_er(8+(Jstar-1)*24) = max(dabs(lnu), dabs(lth), lh+lhe+lc)*clsn
   ! Compute FT and FP correction factors
   !xr(1:kh) = sx(17, 2:kh+1)
   !xm(1:kh) = sx(9, 2:kh+1)
   !call potential(kh, xm, xr, sx(3, 2:nm+1), sx(59, 2:nm+1), ftout, fpout)
   !do ik=1, kh
   !   print *, ik, fpout(ik), ftout(ik)/fpout(ik)
   !end do
   ! Richardson number, angular momentum diffusion coefficients
   do ikk = 1, kh-1
      rich = sx(60, ikk+1)
      domega2 = (sx(59, ikk+1) - sx(59, ikk+2))**2
      if (rich>0.0d0 .and. domega2>0.0d0 .and. rich/domega2 < 1.0d0) then
         ! Unstable
         sx(60, ikk+1) = rich/(1.0d-32 + domega2);
         sx(61, ikk+1) = sx(61, ikk+1)/sx(65, ikk+1)
      else
         ! Stable
         if (rich<0.0d0 .or. domega2<0.0d0) then
            sx(60, ikk+1) = 1.0d32
         else
            sx(60, ikk+1) = rich/(1.0d-32 + domega2);
         end if
         sx(61, ikk+1) = 0
      end if
      !  SX(59, IKK) = OMEGA
   end do


   !------------------------------------------------------------------
   ! DIFFUSION COEFFICIENTS [cm^2/s]
   ! The diffussion coefficients are partially calculated in FUNCS1
   ! (available here through SX) except for the gradients which are
   ! calculated in EQUNS1. In this section the gradients are calculated
   ! again for the inner meshpoints and stored in SX to able to print
   ! them in the mdl file for plotting purposes. -SdM
   !------------------------------------------------------------------
   ! Memory help: Variables available from calculations in FUNCS1
   !  SX(23, IKK) = SGTH
   !  SX(59, IKK) = OMEGA
   !  SX(60, IKK) = XFS(FX_RICH)
   !  SX(61, IKK) = XFS(FX_DDSI)
   !  SX(62, IKK) = XFS(FX_DSSI)
   !  SX(63, IKK) = XFS(FX_VES)
   !  SX(64, IKK) = XFS(FX_VMU)
   !  SX(65, IKK) =  XFS(FX_SGF)
   !  SX(66, IKK) =  XFS(FX_SSSI)
   !------------------------------------------------------------------
   ! We loop from the center to surface when IKK runs from 2 to KH-2
   do ikk = 3, kh-1
      ! Conversion factor: SGF = (4 pi r^2 rho)^2/(dm/dk) [1e11 g/cm**2]
      cgs_sgf =  sx(65, ikk)*1e22 + tiny(0.0_dbl)  ! Make sure it's never 0
      ! Molecular weight gradient (DMU)
      dmuleft  =  ( sx(31, ikk) - sx(31, ikk-1)  )
      dmuright =  ( sx(31,ikk+1 )  - sx(31, ikk) )
      dmu      = 0.5*( dmuleft  +dmuright )
      ! Omega gradient**2 (DW)
      dwleft  = ( sx(59, ikk-1) - sx(59, ikk)   )
      dwright = ( sx(59,  ikk)  - sx(59, ikk+1) )
      dw      = 0.5*( dwleft  +dwright )
      dw      = dw**2
      !------------------------------------------------------------------
      ! [Dcgs_CON] Convective mixing
      ! Either eggletons 72 scheme or pols&Tout01 set by convection_scheme
      ! parameter + possibly artificial extra mixing + weakly mixing of
      ! inner and outer boundary points
      ! = SG/SGF,
      dcgs_con = sx(30, ikk) / cgs_sgf
      ! =====================================
      ! [Dcgs_THL] Thermohaline mixing
      ! = 0.0 or SGTH * DMU / SGF
      dcgs_thl = sx(23, ikk) * pstv(dmu, 0.0d0) / cgs_sgf
      ! =====================================
      ! [Dcgs_DSI] Dynamical Shear instability
      !  =  0.0 or  XFS(FX_DDSI)* (1.0D0 - RICH/DW)**2 /cgs_SGF
      rich = sx(60, ikk)
      dcgs_dsi = 0.0d0
      if( dw > 0.0d0 .and. rich > 0.0d0 .and. rich/dw < 1.0d0 ) then
         ton = (1.0d0 - rich/dw)**2 ! Turn on factor
         dcgs_dsi= ton*sx(61, ikk)/cgs_sgf
      end if
      ! =====================================
      ! [Dcgs_SSI] Secular Shear Instability
      ! = 0.0 or XFS(FX_DSSI)*DW
      ris = sx(66, ikk)
      dcgs_ssi = 0.0d0
      if ( dw > 0.0d0 .and. ris>0.0d0 .and. ris/dw < 1.0d0 ) then
         dcgs_ssi =  sx(62, ikk)*dw / cgs_sgf
      end if
      ! =====================================
      ! [Dcgs_ESC] Eddington-Sweet circulation,
      ! compute net circulation velocity
      !     HP = SX(67, IKK)
      !     ES_HP  = HP * SX(65, IKK) !!!should this be the non converted ////sgf????
      !     ES_VMU = SX(64, IKK) * CFMU*DMU
      !     ES_VES = SX(63, IKK)
      !     ES_VES = PSTV( DABS(ES_VES) - DABS(ES_VMU), 0.0D0)
      !     Dcgs_ESC = 1d-33*ES_VES * ES_HP /  cgs_SGF
      ! SX(76, ) is the diffusion coefficient from Maeder&Zahn (1998),
      ! apart from a factor to account for the molecular weight gradient.
      dcgs_esc = sx(76, ikk)/(max(1.0d0, 1.0d0-sx(78, ikk)*dmu) * cgs_sgf)
      ! =====================================
      ! [Dcgs_GSF] Goldreich-Schubert-Fricke instability
      ! Not yet effective for chemical transport
      dcgs_gsf = 0.0d0 ! not yet properly implemented
      ! =====================================
      ! [Dcgs_SHI] Solberg-Hoiland Instability
      dcgs_shi = 0.0d0 !not yet implemented
      !------------------------------------------------------------------
      ! store gradients in SX:
      sx(69, ikk) = dmu
      sx(70, ikk) = dw
      !    store diffusion coefficients in SX:
      sx(71, ikk) =   dcgs_con ! Convection + artificial mixing
      sx(72, ikk) =   dcgs_thl ! Thermohaline
      sx(73, ikk) =   dcgs_shi ! Solberg-Hoiland
      sx(74, ikk) =   dcgs_dsi ! Dynamical Shear
      sx(75, ikk) =   dcgs_ssi ! Secular Shear
      sx(76, ikk) =   dcgs_esc ! Eddington-Sweet
      sx(77, ikk) =   dcgs_gsf ! Goldberg-Schubert-Fricke
   end do
   !------------------------------------------------------------------
   ! set (next-to) boundary values to zero
   do ikk = 1,2
      sx(69, ikk) =   0.0d0
      sx(70, ikk) =   0.0d0
      sx(71, ikk) =   0.0d0 ! Convection + artificial mixing
      sx(72, ikk) =   0.0d0 ! Thermohaline
      sx(73, ikk) =   0.0d0 ! Solberg Hoiland
      sx(74, ikk) =   0.0d0 ! Dynamical Shear
      sx(75, ikk) =   0.0d0 ! Secular Shear
      sx(76, ikk) =   0.0d0 ! Eddington Sweet
      sx(77, ikk) =   0.0d0 ! Goldberg Schubert Frick
   end do
   do ikk = kh, kh+1
      sx(69, ikk) =   0.0d0
      sx(70, ikk) =   0.0d0
      sx(71, ikk) =   0.0d0 ! Convection + artificial mixing
      sx(72, ikk) =   0.0d0 ! Thermohaline
      sx(73, ikk) =   0.0d0 ! Solberg Hoiland
      sx(74, ikk) =   0.0d0 ! Dynamical Shear
      sx(75, ikk) =   0.0d0 ! Secular Shear
      sx(76, ikk) =   0.0d0 ! Eddington Sweet
      sx(77, ikk) =   0.0d0 ! Goldberg Schubert Frick
   end do

   !------------------------------------------ END DIFFUSION COEFFICIENTS

   ! SX(60, KH+1) = 1.0D0;!SX(60, KH+1)/MAX((SX(59, KH-1) - SX(59, KH))**2, 1.0d-64)
   ! variables needed for MB in FUNCS1
   sm = sx(9, kh+1)
   Teff = px(4)
   qcnv = mcnv/sm

   ! convective envelope turnover time
   drmax = 0.0d0
   imax = 1
   do i = 1, ic
      if ( i == 1 ) ddr = rcb(1)
      if ( i > 1 ) ddr = rcb(i) - rcb(i - 1)
      if ( ddr >= drmax .and. tct(i) /= 0.0d0)  imax = i
      if ( ddr >= drmax .and. tct(i) /= 0.0d0)  drmax = ddr
   end do
   tet = 0.4d0*tct(imax)*crsn/5.76d-7+1.0d-2
   rcz = rcb(imax)/px(17)
   drcz = drmax/px(17)

   tkh = 1.0d22*cg*var(4)*var(4)/(r*var(8)*csy)
   tn(Jstar) = 4.0d10*sm/var(8)
   dmt = sx(34, kh+1)  ! dM/dt
   !  dmt = dvar(4)/(cmsn*dty)
   dmtr = xit(Jstar)*csy/cmsn  ! dMtrans/dt
   dmsw = zet(Jstar)*csy/cmsn  ! dMwind/dt
   oms = om/cmsn
   fr = dlog10(sx(17, kh+1))
   fl = dlog10(sx(18, kh+1))
   sdc = dlog10(sx(3,2))
   sdm = dlog10(sx(3,im_tmax))
   sds = dlog10(sx(3, kh+1))
   stc = dlog10(sx(4,2))
   stm = dlog10(sx(4,im_tmax))
   sts = dlog10(sx(4, kh+1))
   perc = 0.5d0*dlog10(px(17)**3/(8.157d0*sm))
   ! IF M.I. isn't being solved for, VK2 = gyr. rad. squared is nonsense.
   vk2 = dmin1(sx(36,kh)/(var(4)*r*r), 1.0d0)
   raf = dsqrt(dabs(ra2))/r
   call potent ( var(4)/om, dphi )
   ephi = 1.0d22*cg*(var(4) + om)*dphi/sep
   f1 = var(14)/ephi
   df = 0.0
   if ( ktw == 2 ) df = (h(62 - 24*Jstar, 1) + &
        dh(62 - 24*Jstar, 1))/ephi - f1
   horb = var(17)            ! This seems to not have been updated before?
   ecc  = var(18)
   htot = horb + hspn(1) + hspn(2)

   aj = age + dty
   if ( jnn == 0 .and. ktw == Jstar ) age = aj
end subroutine compute_output_quantities





subroutine update_explicit_quantities( jo, Jstar )
   use real_kind
   use mesh
   use mesh_enc
   use control
   use constants
   use settings
   use explicit_functions
   use fgb2hb_composition
   use printb_global_variables
   use test_variables
   use current_model_properties
   
   implicit none
   integer, intent(inout) :: jo
   integer, intent(in) :: Jstar

   real(double) :: w1, dlogMU, dlogP
   integer :: ik

   ! Update the composition of accreted material, for instance during ZAHB
   ! construction. This does not matter for the core, but it may matter for
   ! the envelope or the hydrogen burning shell (according to Haili Hu).
   call update_accretion_abundance

   ! update constant extra energy source
   w1 = dt*cet/csy
   enc = enc*(1.0d0 + w1)*cea
   if (cea /= 0.0d0) enc = enc/(cea + w1*enc)

   if (jnn > 1) enc_parachute = enc_parachute*0.1
   if (enc_parachute < 1.0d-4) enc_parachute = 0.0d0

   if (adj_mea) then
      ! Monitor convergence to target model.
      call check_conv_to_target_structure()
      if ( curr_diffsqr < best_diffsqr ) then
         best_diffsqr = curr_diffsqr
         best_mod = jmod+1
         mutant_h(1:24, 1:kh) = h(1:24, 1:kh)
      end if
      if ( best_diffsqr<1.0d-4 ) jo = 53
   end if

   if ( eps>1.0d-4 ) then
      eps = max(eps*0.8, 1.0d-4)
      print *, 'Set EPS to ', eps
   end if

   ! Adjust the timescale for the artificial energy term. This causes it to be
   ! turned on smoothly.
   if (usemenc .and. impose_entropy_factor<1.0d0) then
      impose_entropy_factor = min(2.0d0*impose_entropy_factor, 1.0d0)
      !PRINT *, 'Set eart factor to ', IMPOSE_ENTROPY_FACTOR
   end if

   ! Adjust the fudge-factor for the artificial composition adjustment. This causes
   ! the composition profile to be changed more gradually by turning the factor
   ! on smoothly.
   if (adj_comp .and. impose_composition_factor<1.0d0) then
      impose_composition_factor = min(3.0d0*impose_composition_factor, 1.0d0);
      !print *, 'Set composition factor to ', IMPOSE_COMPOSITION_FACTOR
   end if

   ! Explicit functions otherwise used in funcs1
   ! Calculate molecular weight gradient.
   ! Don't bother with boundary points (can they even matter?)
   expl_var(1, explv_gradmu, Jstar) = 0.0d0
   expl_var(kh, explv_gradmu, Jstar) = 0.0d0
   do ik=2, kh-1
      dlogMU = 0.5d0 *&
           (expl_var(ik-1,explv_avmu,Jstar) - expl_var(ik+1,explv_avmu,Jstar))
      dlogMU = dlogMU / expl_var(ik, explv_avmu, Jstar)
      dlogP = 0.5d0 *&
           (expl_var(ik-1,explv_logp,Jstar) - expl_var(ik+1,explv_logp,Jstar))
      if (dlogP /= 0.0d0) expl_var(ik, explv_gradmu, Jstar) = dlogMU / dlogP
   end do
end subroutine update_explicit_quantities



! Adjust the timestep control parameter according to rules-of-thumb
subroutine update_timestep_parameters( Jstar )
   use real_kind
   use mesh
   use settings
   use constants
   use plotvariables
   use test_variables
   use structure_variables
   
   implicit none
   integer, intent(in) :: Jstar
   
   real(double), save :: rlf_prev(2) = (/0.0, 0.0/)
   real(double), save :: qcnv_prev(2) = (/0.0, 0.0/)
   real(double), save :: lnuc_prev(2) = (/0.0, 0.0/)
   real(double), save :: lhe_prev(2) = (/0.0, 0.0/)
   real(double), save :: lh_prev(2) = (/0.0, 0.0/)
   real(double) :: lnuc, tnuc

   if ( sx(10, kh+1) < cxb ) mh = sx(9, kh+1)
   if ( sx(11, kh+1) < cxb ) mhe = sx(9, kh+1)
   lnuc = lh + lhe + lc
   cdd = cdc(1)
   ! Allow larger or smaller time-increments, according to rules-of-thumb
   if ( mhe <= mh ) then
      ! Pre-main sequence
      !IF (SX(4,2) < 1.0d7) CDD = CDC(1)*2.0
      ! End of main sequence, reduce timestep to resolve hook
      if (sx(10,2) < 5.0d-2 .and. sx(10,2) > 1.0d-5) cdd = cdc(1)*cdc(6)
      ! Hertzsprung gap/end of core hydrogen burning
      if (sx(10,2) < 1.0d-5 .and. sx(11,2) > 1.0d-5&
           .and. mh < 0.22 .and. sx(20,2) < 1.0) cdd = cdc(1)*cdc(7)
      ! First dredge up; can probably be extended t any dredgeup episode
      if (mh > 0.12 .and. mh < 0.4 .and. qcnv>0.01 .and. qcnv > qcnv_prev(Jstar)) cdd = cdc(1)*cdc(8)
      ! Core helium burning
      if (sx(10,2) < 1.0d-5.and.sx(11,2) < 0.95d0) cdd = cdc(1)*cdc(2)
      ! End of core helium burning
      if ( sx(11,2) < 0.1d0 ) cdd = cdc(1)*cdc(3)
      ! Double shell burning/shells become close (TP-AGB)
      if ( mhe > 0.75d0*mh ) cdd = cdc(1)*cdc(4)
      ! Themal pulse
      if ( mhe>0.75d0*mh .and. lhe>lhe_prev(Jstar) .and. lh<lh_prev(Jstar))&
           cdd = 0.5d0*cdc(1)*cdc(4)
      ! Reduce timestep when approaching Roche lobe overflow -SdM
      if (rlf(Jstar) > -2.0d-2 .and. rlf(Jstar) > rlf_prev(Jstar))&
           cdd = min(cdd, cdc(1)*cdc(9))
      ! Keep timestep small if star is close to RLOF but shrinking
      if (rlf(Jstar) > -2.0d-2 .and. rlf(Jstar) < rlf_prev(Jstar))&
           cdd = min(cdd, cdc(1)*cdc(10))
      ! Switch off Reimers-like wind for He stars
      if ( mhe > 0.96d0*sx(9, kh+1) ) cmr = 0.0d0
   end if
   ! Check for thermo-nuclear "flashes" and cut back the timestep as they are
   ! detected. Based on ideas by Attay Kovetz
   if (lnuc_prev(Jstar) > 0.0d0 .and. lnuc>1.1d0*sx(18, NM+1)) then
      tnuc = dt
      ! Adjust timestep such that we will only allow a timestep that changes the
      ! luminosity by at most 10%.
      if (lnuc>1.1d0*lnuc_prev(Jstar)) then
         tnuc = dlog(1.1d0)*dt / dlog(lnuc / lnuc_prev(Jstar))
         !print *, LNUC/LNUC_PREV(Jstar), TNUC/CSY, DT/CSY
      end if
      cdd = cdd * min(tnuc, dt) / dt
   end if
   rlf_prev(Jstar) = rlf(Jstar)
   qcnv_prev(Jstar) = qcnv
   lnuc_prev(Jstar) = lnuc
   lhe_prev(Jstar) = lhe
   lh_prev(Jstar) = lh

   ! Dirty hacks to help get through He-flash; not all work well and none
   ! work well enough
   if ( lhe > 1.0d-2 .and. dlog10(sx(3, 2)) > 5.3 ) cdd = cdd*0.1
   if ( lhe > 1.0 .and. dlog10(sx(3, 2)) > 5.3 ) cdd = cdd*0.5
   !IF ( LHE > 1.0D5 ) CDD = CDD*0.5
   !IF ( LHE > 1.0D7 ) EPS = 1.0D-4
   !IF ( LHE < 1.0D7 .AND. EPS>1.0D-6) EPS = 1.0D-6
   !IF ( LHE > 1.0D6 ) EPS = 1.0D-2
   !IF ( LHE < 1.0D6 .AND. EPS>1.0D-4) EPS = 1.0D-4
   !IF ( LHE > 1.0D0 ) WANTED_EPS = 1.0D-12
end subroutine update_timestep_parameters




! Adjust the weights in the MSF according to rules-of-thumb
!> \todo FIXME: two stars in a binary may want different MSFs, this is not
!! currently taken into consideration.
!<
subroutine adjust_mesh_spacing(  )
   use real_kind
   use mesh
   use settings
   use printb_global_variables
   use test_variables
   use current_model_properties
   
   implicit none
   
   ! Set defaults
   !CT(:) = INITIAL_CT(:)
   ! Set AGB pressure terms in MSF
   !IF (JMOD < 10) RETURN
   ct(11) = 0.5d0*(0.1d0 * ph + ct(11))
   ct(13) = 0.5d0*(3.0d0 * phe + ct(13))
   ct(15) = 0.5d0*(0.3d0 * phe + ct(15))
   ! PPN phase, increase resolution in the outer regions to reduce numerical
   ! diffusion
   !IF (MHE > 0.99*MH) CT(3) = CT(3) * (1.0 + 50.0*(MHE/MH - 0.99))
   !IF (MH > 0.95*SM .and. CT(3) < 4.0*INITIAL_CT(3)) CT(3) = CT(3) * 1.05
end subroutine adjust_mesh_spacing





subroutine check_stop_conditions ( Jstar, jo, jcm, ift )
   use real_kind
   use mesh
   use settings
   use constants
   use test_variables
   use printb_global_variables
   use current_model_properties
   use structure_variables
   use accretion_abundances
   
   implicit none
   integer, intent(in) :: Jstar
   integer, intent(inout) :: jo, jcm
   integer, intent(in) :: ift
   
   real(double) :: m
   
   
   ! Conditions for terminating *1 or *2; the UC are consts from fort.23
   if ( jb == 1 .and. rlf(1) > uc(1) ) jo = 4 !0.1D0
   if ( age >= uc(2) ) jo = 5 !2.0D10
   if ( lc > uc(3) ) jo = 6  !1.0D2
   if ( jb == 2 .and. rlf(Jstar) > uc(4) ) jo = 7 !0.0D0
   if ( lhe > uc(5) .and. sdc > uc(6) .and. mhe == 0.0d0 ) jo = 8
   if ( mhe > uc(7) .and. sdc > uc(8) ) jo = 9
   if ( dabs(dmt) > uc(9)*sm/tkh ) jo = 10 !30.0D0
   jcm = 1
   if ( sx(10, 2) == 0.0d0 ) jcm = 2
   if ( sx(11, 2) == 0.0d0 ) jcm = 3
   ! Reduce accuracy for He depleted cores (don't want this/not needed)
   !IF ( SX(11,2) < UC(10) ) EPS = UC(11)     !0.15D0, 1.0D-4 !! fudge
   ! restore old value when XHe < 1e-6
   !IF ( SX(11,2) < 1.0D-6 .AND. EPS > 1.0D-6 ) EPS = 1.0D-6
   if ( sx(10,2) < uc(15) ) jo = 51                      ! Stop at TAMS
   if ( px(17) > uc(16) .and. uc(16) > 0.0d0) jo = 52    ! Radius > limit
   ! ZAHB construction
   if (ift == 24) then
      m = sm + 1.0d-6                              ! Target mass + eps
      ! Switch on composition adjustment when the mass exceeds the target core
      ! mass; if we adjust the composition before this point the code runs
      ! very unstable (we're not reconstructing the H burning shell properly).
      if (m > mh) then
         ccac = 1.0d0
      end if
      ! Check for various intermediate stages and goals
      if ( m >= uc(13) ) cmi = 0.0d0   ! Target mass reached, switch off CMI
      if ( mh >= uc(14) ) kx = 0       ! Core mass reached, stop burning
      if ( m>=uc(13) .and. mh>=uc(14) .and. dt>1.0d6*csy) jo = 13 !  Done
      if ( jo == 5) jo = 0    ! Don't look at age during ZAHB construction
   end if
   !if ( mh < mhe .or. mhe < mco ) jo = 14
end subroutine check_stop_conditions







!> FIND_STELLAR_TYPE:
!!
!! Determine the stellar type, more or less as in HPT00 (some of these
!! will never be returned obviously).
!! Return value:
!!  0 - Main sequence star, Qenv > 0.1
!!  1 - Main sequence star, other
!!  2 - Hertzsprung gap
!!  3 - RGB
!!  4 - Core helium burning/horizontal branch
!!  5 - Early AGB
!!  6 - Thermal pulsing AGB
!!  7 - Helium main sequence star / WR star
!!  8 - Helium Hertzsprung gap star
!!  9 - Helium giant
!! 10 - Helium white dwarf
!! 11 - CO white dwarf
!! 12 - O/Ne white dwarf
!! 13 - neutron star
!! 14 - black hole
!! 15 - massless remnant
!< 

function find_stellar_type ( )
   use real_kind
   use mesh
   use settings
   use plotvariables
   use test_variables
   use structure_variables
   
   implicit none
   integer :: find_stellar_type
   
   
   ! Main-sequence stars: hydrogen in the core
   if ( sx(10, 2) > 1.0d-5 ) then
      find_stellar_type = 0
      if (qcnv < 0.1d0) then
         find_stellar_type = 1
      end if
      return
   end if

   ! Hertzsprung gap/end of core hydrogen burning
   if (sx(10,2) < 1.0d-5 .and. sx(11,2) > 1.0d-5&
        .and. mh < 0.22 .and. sx(20,2) < 1.0) then
      find_stellar_type = 2
      return
   end if

   ! Check for evolved stars without nuclear burning
   if (lhe < 1.0 .and. lh < 1.0 .and. lc < 1.0) then
      if (sx(11,2) > 1.0d-5) then      ! He white dwarf
         find_stellar_type = 10
         return
      end if
      if (sx(12,2) > 1.0d-5) then      ! C/O white dwarf
         find_stellar_type = 11
         return
      end if
      ! If none of the above, call it an O/Ne white dwarf
      ! It's not likely to be a neutron star
      find_stellar_type = 12
      return
   end if

   ! Red giant branch: helium core, but no helium burning
   if ( lhe < 3.0 .and. sx(11, 2) > 1.0d-5 ) then
      find_stellar_type = 3
      return
   end if

   ! AGB: inert C/O core
   if ( sx(11, 2) < 1.0d-5 .and. lc < 1.0 ) then
      ! Early AGB: shell sources are well seperated
      if ( mhe < 0.75*mh ) then     ! 0.75? Should be closer?
         find_stellar_type = 5
      else                          ! T-P AGB
         find_stellar_type = 6
      end if
      return
   end if

   ! Some kind of helium star (hopefully, otherwise it's weird)
   if ( sx(11, 2) > 1.0d-5 .and. sx(10, kh+1) < 0.1d0 ) then
      find_stellar_type = 7
      return
   end if

   ! Some kind of helium star (hopefully, otherwise it's weird)
   if ( sx(11, 2) < 1.0d-5 .and. sx(10, kh+1) < 0.1d0 ) then
      if (lc < 1.0) then
         find_stellar_type = 8
      else
         find_stellar_type = 9
      end if
      return
   end if

   ! We should probably ever get here... what ever this is, call it a
   ! Hertzsprung gap star? Presumably it's out of thermal equilibrium...
   !> \todo  Give it a code that says "don't know - probably HG"?
   find_stellar_type = 2
end function find_stellar_type






subroutine write_internal_details ( ig )
   use real_kind
   use settings
   use constants
   use mesh_enc
   use plotvariables
   use printb_global_variables
   use test_variables
   use current_model_properties
   use structure_variables
   
   implicit none
   integer, intent(in) :: ig
   
   integer :: i, j, ij, ikk, ip, ik
   character(len=5) :: CHAR(0:NPX)
   DATA CHAR/'     ',' psi ','  P  ',' rho ','  T  ','kappa','grada',' grad',&
        'gr-ga','  m  ','  H1 ',' He4 ',' C12 ',' N14 ',' O16 ',' Ne20',&
        ' Mg24','  r  ','  L  ',' Eth ',' Enuc',' Eneu','  dm ','SGTH ',&
        'n/n+1','lr/lp','lm/lp',' U   ',' S   ','L/Edd',' w.l ','  mu ',&
        '  wt ',' Ne  ','dm/dt','  wcv',' M.I.','phi  ','  xi ',' DGOS',&
        'dDLdk','Denth',' xik ','v**2 ',' F2  ',' F1  ','  Si ','  Fe ',&
        '  48 ','  49 ',' RPP ',' RPC ',' RPNG',' RPN ',' RPO ',' RAN ',&
        'dS/dP','  LK ','  LQ ','Omega','RichN','DDSI ','  62 ','  63 ',&
        '  64 ','  65 ','  66 ','  67 ','  68 ','  69 ','  70 ', ' 71 ',&
        '  72 ', ' 73 ','  74 ','  75 ','  76 ','  77 ','  78 ','  79 ',&
        '  80 ' /


99001 format (/, '  K', 15(4x, a5, 1x),/)
99002 format (i4, 1p, 15d10.3)
99005 format (/, ' K ', 12(5x, a5, 2x),/)
99006 format (i4, 1p, 12d12.5)
   if(kt3 < 1) return

   write (ig, 99001) (char(ksx(i)), i = 1, 15)
   ! Print the interior details on first `page', if required
   do ikk = 1, kh
      ik = kh + 1 - ikk
      if ( mod(ik-1,kt2) == 0 ) write(ig,99002) ik,&
           (sx(ksx(j), ikk + 1), j=1,15)
      if ( mod(ik,10*kt2) == 0 ) write(ig, 99001) &
           (char(ksx(j)), j = 1, 15)
   end do
   ! Write further `pages' for each detailed model, if required
   do i = 2, kt3
      ij = 15*i - 15
      write(ig, 99005) (char(ksx(ip + ij)), ip=1,12)
      do ikk = 1, kh
         ik = kh + 1 - ikk
         ! Have to choose the SX's that are wanted for printing
         if ( mod(ik-1, kt2) == 0 ) write(ig, 99006) ik,&
              (sx(ksx(ip + ij), ikk + 1), ip = 1, 12)
         if ( mod(ik, 10*kt2) == 0 ) write(ig, 99005) &
              (char(ksx(ip + ij)), ip = 1, 12)
      end do
   end do
end subroutine write_internal_details



subroutine write_summary ( Jstar, jmad, ig )
   use real_kind
   use constants
   use mesh_enc
   use plotvariables
   use printb_global_variables
   use test_variables
   use current_model_properties
   use structure_variables
   
   implicit none
   integer, intent(in) :: Jstar, jmad, ig
   
11993 format (&
           I6, ' M/dty/    Po/P*/e/   xi/zet/   tN/tKH/  LH/LHE/LC Pcr/Rcz ',&
           '  MH/Me/MC/  XH     XHe     XC      XN      XO      XNe     XMg',&
           '      psi   logrho  logT ', /,                                   &
           '  age/cm/Horb     rf1/Fn   mt/rf2/dF tET/DL/BE /Lnu/Lth  ',      &
           'DRcz/RA/Bp  MH/Me convective/semiconvective boundaries         ',&
           '     k**2           logR    logL ',  /,                          &
           3(1p, d16.9, 5d10.3, 0p, f8.3, 7f8.5, 2f8.4, f9.5, 2x, a4, /),   &
           (1p, d16.9, 5d10.3, 0p, f8.3, 7f8.3, 2f8.4, f9.5, 2x, a4, /),   &
           1p, d16.9, 5d10.3, 0p, 8f8.3, 24x, i6, i2, /)

   integer :: jj, i
   real(double) :: dty
   jj = max0(jb, Jstar)
   dty = dt/csy
   ! KTW=2: suffix Jstar for things set in PRINTB, 1 for things set in FUNCS1
   write (ig, 11993) jmad,                                                                                           &
        sm, bper,  dmtr, tn(Jstar), lh,  perc, mh,   (sx(i, 2),       i = 10, 16), sx(1,2),        sdc, stc, 'cntr', &
       dty, psurf, dmsw, tkh,       lhe, rcz,  mhe,  (sx(i, im_tmax), i = 10, 16), sx(1, im_tmax), sdm, stm, 'Tmax', &
        aj, secc,  dmt,  tet,       lc, drcz,  mco,  (sx(i, kh),      i = 10, 16), px(1),          sds, sts, 'srfc', &
       oms, rlf(1:2),    dl,        lnu, raf,  wmh,  mcb,                                           fr, fl,  '    ', &
      htot, f1,    df,   be(Jstar), lth, bp,   wmhe, msb,                     vk2, jmad, jj
   if (adj_mea) then
      write (ig, '(1p, 3d16.9, i6)') curr_maxdiffsqr, curr_diffsqr, best_diffsqr, best_mod
   end if
   call flush ( ig )
end subroutine write_summary



subroutine write_nucsyn_summary ( Jstar, jmad, ig )
   use real_kind
   use constants
   use mesh_enc
   use plotvariables
   use printb_global_variables
   use constants
   use nucleosynthesis
   use test_variables
   implicit none
   integer, intent(in) :: Jstar, jmad, ig
   real(double) :: avghnuc(nvar_nuc)
   real(double) :: m, dm
   real(double) :: dty
   integer :: k, kc, ks, km, i
   ! RJS Column headings for nucleosynthesis output - plus more precise output format
   ! One day I'll make these selectable...
   character(len=5) :: nucname(47)= (/                                                            &
   '  g  ', '  n  ', '  D  ', ' He3 ', ' Li7 ', ' Be7 ', ' B11 ', ' C13 ', ' C14 ', ' N15 ', &
   ' O17 ', ' O18 ', ' F19 ', ' Ne21', ' Ne22', ' Na22', ' Na23', ' Mg24', ' Mg25', ' Mg26', &
   'Al26M', 'Al26G', ' Al27', ' Si28', ' Si29', ' Si30', ' P31 ', ' S32 ', ' S33 ', ' S34 ', &
   ' Fe56', ' Fe57', ' Fe58', ' Fe59', ' Fe60', ' Co59', ' Ni58', ' Ni59', ' Ni60', ' Ni61', &
   ' H1  ', ' He4 ', ' C12 ', ' N14 ', ' O16 ', 'Ne20 ', ' Ca40'/)

   integer :: nucvar(45) = (/&
        41,3,  4,42,  5,  6,  7,  43,8,9, 44,10, 45,11,12, 13, 46,14,15,&
        16,17, 18,19,20, 22,23, 24,25,26, 27, 28,29,30, 47, 31,32,33,34,35,&
        36, 37,38,39,40,1&
        /)

   dty = dt/csy

   ! Nothing to do if we're not computing nucleosynthesis output
   if ( .not. allocated(hnuc) ) return

   ! Calculate mass-averaged abundances
   avghnuc(:) = 0.0d0
   m = 0.0d0
   do k=1, kh_nuc-1
      dm = dabs(ht(Jstar, 4, k))
      m = m + dm
      avghnuc(:) = avghnuc(:) + hnuc(Jstar,:,k)*dm
   end do

   avghnuc(:) = avghnuc(:)/m
   kc = kh_nuc-1
   ks = 1
   km = kc - (im_tmax-1) + 1

   write (ig, 11994) jmad,&
        (nucname(nucvar(i)), i = 1, 14),&
        sm, (hnuc(Jstar, nucvar(i), kc), i=1, 14),      'cntr',&
        dty, (hnuc(Jstar, nucvar(i), km), i = 1, 14),    'Tmax',&
        aj, (hnuc(Jstar, nucvar(i), ks), i = 1, 14),    'surf',&
        (avghnuc(nucvar(i)), i = 1, 14),     'Mavg',&
        (nucname(nucvar(i)), i = 15, 28),&
        stc, (hnuc(Jstar, nucvar(i), kc), i=15, 28),     'cntr',&
        stm, (hnuc(Jstar, nucvar(i), km), i = 15, 28),   'Tmax',&
        sts, (hnuc(Jstar, nucvar(i), ks), i = 15, 28),   'surf',&
        fl, (avghnuc(nucvar(i)), i = 15, 28),    'Mavg',&
        (nucname(nucvar(i)), i = 29, 42),&
        sdc, (hnuc(Jstar, nucvar(i), kc), i=29, 42),     'cntr', &
        sdm, (hnuc(Jstar, nucvar(i), km), i = 29, 42),   'Tmax', &
        sds, (hnuc(Jstar, nucvar(i), ks), i = 29, 42),   'surf',&
        fr, (avghnuc(nucvar(i)), i = 29, 42),    'Mavg'
   call flush ( ig )

11994 format (/,&
           I6,' M/dty/age ', 14('  ',A5,'   '),/,&
           3(1P, D16.9, 14D10.3, 2X, A4, /),&
           (1P, 16X,   14D10.3, 2X, A4, /),&
           6X,' logT/logL ', 14('  ',A5,'   '),/,&
           4(1P, D16.9, 14D10.3, 2X, A4, /),&
           6X,'logrho/logR', 14('  ',A5,'   '),/,&
           4(1P, D16.9, 14D10.3, 2X, A4, /))
end subroutine write_nucsyn_summary





subroutine write_plt_file ( Jstar, jmad, ig )
   use real_kind
   use constants
   use mesh_enc
   use plotvariables
   use printb_global_variables
   use test_variables
   use current_model_properties
   use structure_variables
   
   implicit none
   integer, intent(in) :: Jstar, jmad, ig
   
   logical :: first_time(2) = (/.true., .true./)
   integer :: i, j
   real(double) :: dty
   
   
   dty = dt/csy

   ! Write number of cols to the first line of the star.plt[12] files
   if(first_time(Jstar)) then
      rewind(ig)
      write(ig,'(i4)') 83
      first_time(Jstar) = .false.
   end if
   write (ig,12100) jmad,aj,dty,sm,mh,mhe,mco,&
        log10(px(17)),log10(px(18)),log10(px(4)),log10(sx(4,2)),stm,&
        log10(sx(3,2)),sdm,be(Jstar),lh,lhe,lc,lnu,lth,&
        psurf,vk2,rcz,drcz,tet,raf,bp,&
        bper,rlf(Jstar),f1,dmt,dmsw,dmtr,horb,dhdt,dhgw,dhml,dhso(Jstar),&
        dhmt(Jstar),oms,ecc,&
        (px(j),j=10,16),(sx(i,im_tmax),i=10,16),(sx(j,2),j=10,16),&
        (mcb(j),j=1,6),(msb(j),j=1,6),(mex(j),j=1,6),qcnv,sx(2,2), pcntr
   call flush(ig)
12100 format (i6,es17.9,es14.6,12es13.5,7es12.4,3es13.5,17es12.4,&
           39es13.5,es14.6,es13.5,es14.6)
end subroutine write_plt_file



subroutine write_nucsyn_plt_file ( Jstar, jmad, ig )
   use real_kind
   use constants
   use mesh_enc
   use plotvariables
   use printb_global_variables
   use nucleosynthesis
   use test_variables
   use current_model_properties
   use structure_variables
   
   implicit none
   integer, intent(in) :: Jstar, jmad, ig
   character(len=5) :: nucname(47)= (/                                                            &
   '  g  ', '  n  ', '  D  ', ' He3 ', ' Li7 ', ' Be7 ', ' B11 ', ' C13 ', ' C14 ', ' N15 ', &
   ' O17 ', ' O18 ', ' F19 ', ' Ne21', ' Ne22', ' Na22', ' Na23', ' Mg24', ' Mg25', ' Mg26', &
   'Al26M', 'Al26G', ' Al27', ' Si28', ' Si29', ' Si30', ' P31 ', ' S32 ', ' S33 ', ' S34 ', &
   ' Fe56', ' Fe57', ' Fe58', ' Fe59', ' Fe60', ' Co59', ' Ni58', ' Ni59', ' Ni60', ' Ni61', &
   ' H1  ', ' He4 ', ' C12 ', ' N14 ', ' O16 ', 'Ne20 ', ' Ca40'/)

   integer :: nucvar(45) = (/&
        41,3,  4,42,  5,  6,  7,  43,8,9, 44,10, 45,11,12, 13, 46,14,15,&
        16,17, 18,19,20, 22,23, 24,25,26, 27, 28,29,30, 47, 31,32,33,34,35,&
        36, 37,38,39,40,1&
        /)
   
   logical :: first_time(2) = (/.true., .true./)
   integer :: i
   real(double) :: dty

   dty = dt/csy

   if ( .not. allocated(hnuc) ) return

   ! Write number of cols to the first line of the star.plt[12] files
   if(first_time(Jstar)) then
      rewind(ig)
      write(ig, '(1p,1x,46a6)') '#', (nucname(nucvar(i)), i=1,45)
      write(ig,'(i4)') 14+45
      first_time(Jstar) = .false.
   end if
   
   ! The first 14 output quantities are the same as for the normal .plt file
   write (ig,12100) jmad,aj,dty,sm,mh,mhe,mco, &
        log10(px(17)),log10(px(18)),log10(px(4)),log10(sx(4,2)),stm, &
        log10(sx(3,2)),sdm,(hnuc(Jstar, nucvar(i), 1), i=1,45)
   
   call flush(ig)
12100 format (i6,es17.9,es14.6,12es13.5,7es12.4,3es13.5,17es12.4, &
           39es13.5,es14.6,es13.5,es14.6)
end subroutine write_nucsyn_plt_file



! Write structure output (the values of some physical parameters on
! each meshpoint) to a file usually called star.mdl1 or star.mdl2
!
! INPUT VARIABLES
! IG = indicates which starr: 1 for primary, 2 for secondary component
! JMAD = model number
subroutine write_mdl_file ( Jstar, jmad, ig )
   use real_kind
   use mesh
   use settings
   use plotvariables
   use printb_global_variables
   use structure_variables
   
   implicit none
   integer, intent(in) :: Jstar, ig, jmad
   
   logical :: first_time(2) =(/ .true., .true./)
   integer :: kk, j
   integer, parameter :: nwrt5 = 42
   integer :: ipx(nwrt5)
   ! IPX: List of Indices of the entries in the array SX containing
   ! physcial pars as calculated during latest call cfuncs. This
   ! determines the order in which these pars are written to the output
   ! file .mdl[12].
   data ipx/9, 17,  2,  3,  4,    5,  6 , 8, 10, 11,   &
        12, 13, 14, 15, 16,   18, 19, 20, 21, 28,   &
        27, 50, 51, 52, 53,   54, 55, 56, 31, 23,&
        30,  7, 59, 69, 70,   71, 72, 73, 74, 75,&
        76, 77/! 67, 68, 60, 60, 61,   66, 63, 64, 65/!, 24, 25, 26/

   ! If this function is called for the first time, write a header
   ! (expected by onno's yorick routines)
   if ( first_time(Jstar)) then
      rewind(ig)
      write (ig,'(2i6,f7.3)') kh, nwrt5, cos
      first_time(Jstar) = .false.
   end if
   ! Every time this fuction is called write a block header, and a datablock
   write (ig,12201) jmad, aj
   do kk = 1, kh
      ! Limit range of output exponent so it fits in the output format
      write (ig,12202) &
           (dsign(min(dabs(sx(ipx(j),kk+1)),1.0d99),&
           sx(ipx(j),kk+1) ), j=1, nwrt5)
   end do
   call flush(ig)
12201 format (i6,1p,e17.9,0p)
12202 format (1p,e13.6,4e11.4,80e11.3,0p)
end subroutine write_mdl_file


! Write structure output for nucleosynthesis (the values of some
! physical parameters on each meshpoint) to a file usually called
! star.nucmdl1 or star.nucmdl2
!
! INPUT VARIABLES
! IG = indicates which star: 1 for primary, 2 for secondary component
! JMAD = model number
subroutine write_nucsyn_mdl_file ( Jstar, jmad, ig )
   use real_kind
   use mesh
   use settings
   use plotvariables
   use printb_global_variables
   use nucleosynthesis
   use structure_variables
   
   implicit none
   integer, intent(in) :: Jstar, ig, jmad
   
   logical :: first_time(2) =(/ .true., .true./)
   integer :: kk, j
   integer, parameter :: nwrt5 = 8     ! Parameters to write from SX
   integer :: nucvar(45) = (/         &! Nucleosynthesis abund. perm.
        41,3,  4,42,  5,  6,  7,  43,8,9, 44,10, 45,11,12, 13, 46,14,15,&
        16,17, 18,19,20, 22, 23, 24,25,26, 27, 28,29,30, 47, 31,32,33,34,35,&
        36, 37,38,39,40,1&
        /)
   integer :: ipx(nwrt5)
   ! IPX: List of Indices of the entries in the array SX containing
   ! physcial pars as calculated during latest call cfuncs. This
   ! determines the order in which these pars are written to the output
   ! file .mdl[12].
   data ipx/9, 17,  2,  3,  4,    5,  6 , 8/

   if ( .not. allocated(hnuc) ) return

   ! If this function is called for the first time, write a header
   ! (expected by onno's yorick routines)
   if ( first_time(Jstar)) then
      rewind(ig)
      write (ig,'(2i6,f7.3)') kh, nwrt5, cos
      first_time(Jstar) = .false.
   end if
   ! Every time this fuction is called write a block header, and a datablock
   ! The first 8 variables are the same as in the structure output for the
   ! .mdl1 file
   write (ig,12201) jmad, aj
   do kk = 1, kh
      ! Limit range of output exponent so it fits in the output format
      write (ig,12202) &
           (dsign(min(dabs(sx(ipx(j),kk+1)),1.0d99),&
           sx(ipx(j),kk+1) ), j=1, nwrt5),&
           ( hnuc(Jstar, nucvar(j), kh+1 - kk), j=1,45 )
   end do
   call flush(ig)
12201 format (i6,1p,e17.9,0p)
12202 format (1p,e13.6,4e11.4,80e11.3,0p)
end subroutine write_nucsyn_mdl_file


