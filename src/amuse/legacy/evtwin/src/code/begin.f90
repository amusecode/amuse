!> Read first model from JIP and write to JOP
!! Read init.dat from IT
!<
subroutine star12 ( jo, jcm, jop, jip, ksv, it )
   use real_kind
   use mesh
   use mesh_enc
   use control
   use settings
   use extrapolate_dh
   use explicit_functions
   use nucleosynthesis
   use test_variables
   use current_model_properties
   use constants
   
   implicit none
   integer :: jop,jip,kp,kr,kr1,kr2,ksv,kt5,it,jo,jf,jcm
   integer :: iter,ik,Jstar,Jstop,itry
   real(double) :: dty,rz,prz
   real(double) :: cdd_1, cdd_2
   real(double) :: dx12, dx23
   
   if(debug.ge.2) write(6,'(/,A20,6I12)')'star12: ',jo, jcm, jop, jip, ksv, it

   
   
   ! Read an initial model:
   jnn = 0
   expl_var(:, :, :) = 0.0d0
   call beginn ( jip, dty, kp, kr1, kr2, ksv, kt5, it, jo, jf )
   if ( jo == -2 ) goto 3
   call printb ( jo, jcm, rz, 1, it)
   if ( jo.eq.1 ) goto 3
   if ( ktw.eq.2 ) then
     !average change wanted by *1
      cdd_1 = cdd  !zq(30)
      call printb ( jo, jcm, rz, 2, it )
     !average change wanted by *2
      cdd_2 = cdd  !zq(30)
     !set average change to minimum of the two components
      cdd = min(cdd_1, cdd_2)  !zq(30)
   end if
   if ( kp.eq.0 ) goto 3
   call nextdt ( dty, jo, it )
   if ( jo.eq.3 ) goto 3
   iter = kr1
  
   ! Begin evolutionary loop of KP time steps:
   jnn = 1
   jter = 0
   prz = 0
   do_kp: do
      jo = 0
      ! Solve for structure, mesh, and major composition variables
      joc = 1

      call smart_solver ( iter, id, kt5, jo )
      call relax_timestep ( iter, id, kt5, jo, dty )

      ! If no convergence, restart from 2 steps back, DT decreased substantially
5     continue
      if (jo /= 0) then
         call backup ( dty, jo )
         if (use_quadratic_predictions) then
            call restore_previous_model()
            call predict_dh(dty, dh)
         end if
         if ( jo.eq.2 ) exit do_kp    ! Abort if timestep below limit
         goto 4
      end if
     
      ! If required, solve for minor composition variables; mesh, structure fixed.
      if (ktw > 1) call update_accretion_abundances2()
      !call smooth_spikes()
      call update_concentration_gradients()
      call update_radative_accelerations2()
      if (allocated(dhnuc)) dhnuc(1:ktw, :, 1:kh_nuc) = 0.0d0
      do_jstar: do Jstar = 1, ktw
         joc = Jstar + 1
         call solver ( kr_nucsyn, ie, kt5, jo )

         ! Set small abundances to 0 if no convergence
         if (jo /= 0) then
            jo = 0
            dhnuc(Jstar, :, 1:kh_nuc) = 0.0d0
            where (hnuc < 1.0d-20) hnuc = 0.0d0
            do ik = 1, kh_nuc
               if (h(5 + (jstar-1)*24, ik) == 0.0d0) hnuc(jstar, 41, ik) = 0.0d0
            end do
            call solver ( kr_nucsyn, ie, kt5, jo )
         end if

         if (jo /= 0) then
            ! Smooth "spikes", if they're likely to cause the problem
            !call smooth_spikes()
            !jo = 0
            !dhnuc(Jstar, :, 1:kh_nuc) = 0.0d0
            !call solver ( kr_nucsyn, ie, kt5, jo )
            do_itry: do itry = 1, 5
               if (point_max_err > 1 .and. point_max_err < kh) then
                  dx12 = hnuc(jstar, var_max_err, point_max_err-1) - hnuc(jstar, var_max_err, point_max_err  )
                  dx23 = hnuc(jstar, var_max_err, point_max_err  ) - hnuc(jstar, var_max_err, point_max_err+1)
                  if ( dx12*dx23 < 0.0_dbl .and. H(2 + 24*(Jstar - 1), point_max_err) < 6*CLN ) then
                     call mix_gridpoint(Jstar, point_max_err)
                     jo = 0
                     dhnuc(Jstar, :, 1:kh_nuc) = 0.0d0
                     call solver ( kr_nucsyn, ie, kt5, jo )
                     if (jo == 0) exit do_itry
                  else
                     exit do_itry
                  end if
               else
                  exit do_itry
               end if
            end do do_itry
         end if
         
         if (jo /= 0) then
            jo = 17
            exit do_jstar
         end if
      end do do_jstar  ! do Jstar = 1, ktw
      
      ! If no convergence, restart from 2 steps back, DT decreased substantially
      if (jo /= 0) then
         call backup ( dty, jo )
         if (use_quadratic_predictions) then
            call restore_previous_model()
            call predict_dh(dty, dh)
         end if
         if ( jo.eq.2 ) exit do_kp    ! Abort if timestep below limit
         goto 4
      end if
     
      ! If mass < 0, probably NaN, so exit with code 99
      ! (works for Lahey F95, but probably compiler-dependent!)
      if(h(4,1).lt.0.d0) jo = 99
     
      ! If model didn't converge, give up
      if ( jo >= 1 ) exit do_kp

      prz = rz
      call printb ( jo, jcm, rz, 1, it )
      if ( ktw.eq.2 ) then
        !average change wanted by *1
         cdd_1 = cdd  !zq(30)
         call printb ( jo, jcm, rz, 2, it )
        !average change wanted by *2
         cdd_2 = cdd  !zq(30)
        !set average change to minimum of the two components
         cdd = min(cdd_1, cdd_2)  !zq(30)
      end if

      if ( jo.eq.14 ) goto 5

      ! Check for manual termination
      read (11, *) jstop
      if ( jstop.ne.0 ) jo = 15
      close (11)
      write (11, *) 0
      close (11)

      if ( jo >= 2 ) exit do_kp
      call update ( dty )
      call update2
      if (use_quadratic_predictions) call store_current_model()
     
      ! 2nd condition saves every model, doesn't work with mod:
      if ( mod(jmod, ksv) == 1 .or. ksv == 1 ) call output ( kp, 15, jo, jf )

4     continue
      call nextdt ( dty, jo, it )
      if (use_quadratic_predictions) call predict_dh(dty, dh)
      if ( jo == 3 ) exit do_kp
      if ( jo == 5 ) exit do_kp       ! Reached age limit
      if ( jnn >= 4 ) iter = kr2
      jnn = jnn + 1
      
      ! Don't break the run just as RLOF starts!
      kr = abs(rz)/(abs(rz - prz) + 1.0d-6)
      
      if ( kp == 1 ) exit do_kp
      if (.not.( (jnn<=kp .or. kp<0) .or. jhold<3 .or. jter>2 .or. kr<=10 )) exit do_kp
   end do do_kp     ! Evolutionary loop of KP time steps:  do - ~150 lines up
3  continue

   ! Output the last converged model, unless it's at the He flash
   call output ( kp, jop, jo, jf )
   ! Save last model to file.mod as well (suggested by Selma)
   if ( .not. (mod(jmod, ksv).eq.1 .or. ksv.eq.1) )  call output ( kp, 15, jo, jf )
   if ( jo.ne.8.and.jo.ne.13 .and. isb.eq.2) write (jb, *) 'Switch to *', 3 - jb
   if ( jo.eq.8 ) write (jb, *) 'He flash'
   if ( jo.eq.13 ) write (jb,*) 'Start ZAHB'
   call flush ( jb )
end subroutine star12



!> Called from star12
!! - jip:  unit for input structure model
!! - it:   unit of 'init.dat'
!<
subroutine beginn ( jip, dty, kp, kr1, kr2, ksv, kt5, it, jo, jf )
   use real_kind
   use mesh
   use mesh_enc
   use extra_elements
   use constants
   use settings
   use model_initialiser
   use init_dat
   use control
   use test_variables
   use current_model_properties
   use binary_history
   use accretion_abundances
   use radiative_acceleration
   use nucleosynthesis, only: nucleosynthesis_enabled, Hnuc, DHnuc, Hnucpr
   use filenames
   use file_exists_module
  
   implicit none
   integer, intent(in) :: jip
   integer :: kp,kr1,kr2,ksv,kt5,it,jo,jf
   integer :: ii,kb,kh2,jch,jin,ik,ij,ir,i
   real(double) :: dty,bms,ecc,p1,tm,oa
   real(double) :: nuc_dummy(50),bm_bck,bper_bck

   logical :: status,post_he_flash
   integer :: ioerror, keq
   character(len=500) :: jip_name,fname
   
   if(debug.ge.2) write(6,'(/,A20,I12,ES20.5,8I12)')'beginn: ',jip, dty, kp, kr1, kr2, ksv, kt5, it, jo, jf


   
   h(:, 1:kh) = 0.0d0
   dh(:, 1:kh) = 0.0d0
   hspn = 0.0d0
   rlf = 0.0d0
   zet = 0.0d0
   xit = 0.0d0
   age = 0.0d0
   mc = 0.0d0
   om = 0.0d0
   
   bm_bck = bm
   bper_bck = bper
   bm = 0.0d0
   bper = 0.0d0

   sm = 0.0d0
   enc = 0.0d0

   tc = 0.0d0
   tfr = 0.0d0
   t0 = 0.0d0
   m0 = 0.0d0
   mta = 0.0d0
   om0 = 0.0d0
   omta = 0.0d0
   a0 = 0.0d0
   ata = 0.0d0
   e0 = 0.0d0
   eta = 0.0d0
   cdd = 0.0d0

   bp = 0.0d0
   horb = 0.0d0
   ro = 0.0d0
   ra2 = 0.0d0
   rs = 0.0d0
   secc = 0.0d0
   tn = 0.0d0
   wmh = 0.0d0
   wmhe = 0.0d0
   mh = 0.0d0
   mhe = 0.0d0
   mco = 0.0d0

   vmg = 0.0d0
   be = 0.0d0
   lh = 0.0d0
   lhe = 0.0d0
   lc = 0.0d0
   lnu = 0.0d0
   lth = 0.0d0
   mcb = 0.0d0
   msb = 0.0d0
   rcb = 0.0d0
   tct = 0.0d0

   sect(1:9999) = 1.0d20
   kb = 1
   ! Read miscellaneous data, usually unchanged during one evol run
   ! Try to read in init.dat; abort in case this fails
   status = read_init_dat(it, kh2, kr1, kr2, ksv, kt5, jch)
   if (status .eqv. .false.) then
      inquire(unit=it, name=jip_name)
      write(0, *) 'error reading init.dat data from "', trim(jip_name),'"'
      stop
   end if

   !     Just constructed a ZAHB model; DON'T homogenise the model,
   !     whatever the input .dat file says.
   !     Without this flag, the post-he flash model will have its
   !     composition reset to ZAMS composition if we started from a
   !     homogeneous model - which is not what we want.
   !     We should also force the accretion mode to exponential.
   post_he_flash = .false.
   if (jo == 13) then
      jch = min(3, jch)
      jo = -1
      cmi_mode = 1
      post_he_flash = .true.
   end if

   if (mutate .and. kh>0) then
     ! Don't rezone the model
      kh2 = kh
   end if

   ! Read opacity data and construct splines
   ! KOP (read from init.dat) sets which type of opacity tables
   if (kop<=1) then
     ! Iglesias & Rogers (1992), as implemented by Pols & al. 1995
      call load_opacity(20)
   else
     ! Iglesias & Rogers (1996), as implemented by Eldridge & Tout 2003
      call load_opacity_co(41)
   end if

   ! Allocate data for nucleosynthesis calculations
   call allocate_nucleosynthesis_data(kh2)

   ! Radiative levitation
   if (CRLEV > 0.0d0) then
      call set_op_directory(op_data_path)
      if (.not. init_radiative_acceleration()) then
         print *, '*** Error: cannot initialise code for radiative accelerations ***'
         print *, 'Make sure that the variable "op_data_path" in init.dat is set and points to'
         print *, 'the directory containg the OP data files.'
         print *, '(Current value: '//trim(op_data_path)//')'
         stop
      end if
   end if

   ! Abort if the requested number of meshpoints is larger than the
   ! size of the array
   if (kh2 > NM) then
      write (0, *) 'cannot rezone to ', kh2, 'meshpoints. maximum size is ', NM, 'meshpoints.'
      jo = -2
      return
   end if

   ! Autodetect if we should solve for Mg24 or not by checking if the
   ! corresponding equation is in the list of equations
   use_mg24_eqn = .false.
   do ii = 51, 100
      if (id(ii) == emg24 .or. id(ii) == esumx) then
         use_mg24_eqn = .true.
         exit
      end if
      if (id(ii) == 0) exit   ! Break loop if end of list found
   end do

   ! Detect whether rotation is treated as solid body rotation or
   ! whether we consider differential rotation. The switch is on whether
   ! or not the rotational period (var. 13) is listed as an eigenvalue.
   rigid_rotation = .true.
   keq = id(1)+id(2)    ! Total number of first+second order equations
   do ii = 11, 10 + keq
      if (id(ii) == 13) then  ! Rotational pertiod not an EV
         rigid_rotation = .false.
         exit
      end if
   end do


   ! Read data for initial model (often last model of previous run)
   ! e.g. SM = stellar mass, solar units; DTY = next timestep, years
   do ii=1, ktw
      read(jip, *, iostat=ioerror) sm, dty, age, bper, bms, ecc, p1, enc, kh, kp, jmod, jb, jin, jf
      if ( ioerror /= 0 ) then
         inquire(unit=jip, name=jip_name)
         write(0, *) 'error reading "', trim(jip_name),'"'
         stop
      end if
      if ( jmod == 0 .and. kp>0 ) kp = kp + 1
      if ( jip.eq.13 .or. jip.eq.14 ) dty = ct3*dty
      if(post_he_flash .and. ktw.eq.1) then                ! Directly after He flash, restore orbital period and binary mass
         bper = bper_bck
         bms = bm_bck/cmsn
      end if
      
      if ( ii.eq.2 ) jb = 1
      if ( jmod.lt.10 .and. ii.eq.1 ) &                    ! Write "init.dat" to file.out?
           write (jb, '(8i5, /, 7i5, /, 6i5, /, 1p, 8d8.1, 0p, /, 2(10i4, /,'//  &
           '6(20i3, /)), 3(15i3, /), i3, /, 2(20i3, /), 10f5.2, 1p, 3d8.1, /, 0p, 7f6.3, /, 1p, 5(9d9.2, /), 0p)')  &
           kh2, kr1, kr2, jch, kth, kx, ky, kz, kcl, kion, kam, kop, kcc, knuc,  &
           kcn, kt1, kt2, kt3, kt4, kt5, ksv, eps, del, dh0, cdc(1:5), id, ie, ksx,   &
           kn, kjn, ct1, ct2, ct3, ct,  &
           cc, cn, co, cne, cmg, csi, cfe,  &
           calp, cu, cos, cps, crd, cxb, cgr,   &
           cea, cet, cmt, cms, cmi, cmr, cmj,   &
           cml, chl, ctf, clt, cpa, cbr, csu, csd, cdf, cgw, cso, cmb,   &
           cq1, cth, cq
      write (jb, '(1x, 1p, 2d14.6, d17.9, 5d14.6, 0p, 6i5)') sm, dty, age, bper, bms, ecc, p1, enc, kh, kp, jmod, jb, jin, jf
      
      ! Read the initial model
      do ik = 1, kh
         read (jip, *, iostat=ioerror) (h(ij,ik), ij = 1, jin)
         if ( ioerror /= 0 )then
            inquire(unit=jip, name=jip_name)
            write(0, *) 'error reading "', trim(jip_name),'"'
            stop
         end if
      end do
      
      if (mutate) p1 = 1.0e10
      
      ! Read DH(:) if the JF flag says it's available (bit 2 is set)
      if (iand(jf, 4)==4) then
         store_changes  = .true.
         if ( jip.eq.13 .or. jip.eq.14 ) dty = dty / ct3    ! Undo scaling for timestep
         do ik = 1, kh
            read (jip, *, iostat=ioerror) (dh(ij,ik), ij = 1, jin)
            if ( ioerror /= 0 )then
               inquire(unit=jip, name=jip_name)
               write(0, *) 'error reading "', trim(jip_name),'"'
               stop
            end if
         end do
      end if
      
      tm = cmsn*sm
      
      
      ! Decrease timestep for mutation runs
      if ( mutate .and. jmod>0 ) dty = min(dty, 1.0d3)
      
      ! Make sure the timestep is small when starting ZAHB construction
      if (it == 24) then
         age = 0.0d0
         dty = min(dty, 1.0d3)
      end if
      
      ! Convert some things to `cgs' units: 10**11 cm, 10**33 gm, 10**33 erg/s
      if ( ii == 1 ) then
         dt = csy*dty
         cmi = cmi/csy
         cms = cms/csy
         cmt = cmt*1.0d-11
         if ( cea.le.0.0d0 ) cea = 1.0d-10
         
         ! Initialise orbital angular momentum
         bm = cmsn*bms
         om = bm - tm
         oa = cg1*tm*om*(cg2*bper/bm)**c3rd*dsqrt(1.0d0 - ecc*ecc)
         
         if ( jb /= 1 ) then
            ! For *2 of a binary, read in the mass-loss history of *1 from io12
            !The output file at unit 3 = io12 may not be open yet!
            !> \todo Want to do this in a more central place, e.g. filenames.f90?
            inquire(unit=3,opened=status)
            if(.not.status) then
               write (fname, '(A,i1)') 'fort.',3
               if (.not. file_exists(fname)) then
                  open(unit = 3, file=trim(basename)//'.io12')
               end if
            end if
            
            do ik = 1, 9999
               read (3, 102, end=8) ms(ik), ww, sdt(ik), ww, sect(ik), se(ik), wx, scm(ik), ww,  sang(ik), ww, ir
               if ( ik.gt.1 .and. sect(ik).gt.1.0d19 ) goto 8   ! FIXME - replace with exit? - may fail if more code is added below
            end do
102         format (3(1p, d16.9, 5d10.3, 0p, f8.3, 7f8.5, 3f8.4, /),  &
                 (1p, d16.9, 5d10.3, 0p, f8.3, 7f8.3, 3f8.4, /),  &
                 1p, d16.9, 5d10.3, 0p, 11f8.3, i6)
            
         end if
      end if
      
8     continue
      age = age - dty
      call nextdt ( dty, jo, it )
      
      ! Load the nucleosynthesis data structures:
      call set_initl_nucleosynt_abundances(ii)
      
      if (iand(jf, 8)==8) then      ! Read nucleosynthesis data
         do ik = 1, kh
            if (nucleosynthesis_enabled) then
               read (jip, *, iostat=ioerror) (Hnuc(ii,ij,ik), ij = 1, 50)
            else
               ! Nucleosynthesis datastructures not allocated, but we still need to read the data (in order to skip them)
               read (jip, *, iostat=ioerror) (nuc_dummy(ij), ij = 1, 50)
            end if
            if ( ioerror /= 0 )then
               inquire(unit=jip, name=jip_name)
               write(0, *) 'error reading nucleosynthesis data from "', trim(jip_name),'"'
               stop
            end if
         end do
      end if
      
      ! Sanity check
      if (jch < 3 .and. kh2 /= kh) then
         write(0, *) 'Warning: input mesh not equal to requested mesh, but no remesh specified.'
         write(0, *) 'Will use input mesh of ', kh, ' mesh points.'
         write(0, *) 'Change KH2 or JCH in init.dat to remove this warning.'
      end if
      
      ! optionally rezone the model, e.g. for different no. of meshpoints.
      ! also initialise some variables that were not in input model.
      if (use_smooth_remesher) then
         !   CALL REMESH ( KH, JCH, BM, H(4,1), H(13,1), ECC, OA, II, JF )
         call remesher ( kh2, jch, bm, tm, p1, ecc, oa, ii, jf )
      else
         call remesh ( kh2, jch, bm, tm, p1, ecc, oa, ii, jf )
      end if
      
      ! ZAHB construction, we need to take special care here since the ZAHB
      ! starting model will not have the correct envelope abundance. This is
      ! corrected by accreting material of the right composition, but we
      ! cannot do this immediately because we need to H burning shell to stay
      ! active and dumping pure He on top of it will make it less numerically
      ! stable.
      ! Composition accretion is switched on in printb when the stellar mass
      ! reaches the desired core mass.
      if (it == 24) ccac = 0.0
      if ( ii == 1 ) then
         ! Backup the model in case we need to go back a timestep
         ! Backup the extra variables
         if (nxvar>0) hpr(41:nvar, 1:kh) = h(41:nvar, 1:kh)
         hpr(1:16, 1:kh) = h(1:16, 1:kh)
      end if
   end do  ! do ii=1,ktw
   


   if ( ktw == 2 ) then
      h(25:40, 1:kh) = h(1:16, 1:kh)
      hpr(25:40, 1:kh) = h(1:16, 1:kh)
      h(1:16, 1:kh) = hpr(1:16, 1:kh)
      if (nxvar>0) then
         h(41+nxvstar:40+2*nxvstar, 1:kh) = h(41:40+nxvstar, 1:kh)
         hpr(41+nxvstar:40+2*nxvstar, 1:kh) = h(41:40+nxvstar, 1:kh)
         h(41:nvar, 1:kh) = hpr(41:nvar, 1:kh)
      end if
   end if

   if (nucleosynthesis_enabled) then
      DHnuc = 0.0d0
      Hnucpr = Hnuc
   end if

   ! store some numbers for possible restart with BACKUP
   jhold = 2
   prev(2:81) = (/hspn, rlf, zet, xit, age, bm, mc, om, bper, sm, enc,  &   ! zq(1:17)
       tc, tfr, t0, m0, mta, om0, omta, a0, ata, e0, eta, cdd,          &   ! zq(18:30)
       bp, horb, ro, ra2, rs, secc, tn, wmh, wmhe, mh, mhe, mco,        &   ! zq(31:43)
       vmg, be, lh, lhe, lc, lnu, lth, mcb, msb, rcb, tct/)            


   pprev(1:81) = prev(1:81)
   jm1 = jmod
   ! Determine whether I and phi are computed or not, for OUTPUT
   jf = 0
   do i = 11, 50
      if ( id(i).eq.12 .or. id(i).eq.14 ) jf = jf + 1
   end do
   if ( jf.eq.1 ) jf = 0
   ! Store derivative information
   if (store_changes) jf = jf + 4
   ! Store nucleosynthesis data
   if (nucleosynthesis_enabled) jf = jf + 8
   return
end subroutine beginn

