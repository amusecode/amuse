!---  some constant parameters to scale array sizes
      INCLUDE 'paras.inc'

!---  general MPI vars also needed in serial runs
      INTEGER   rootRank,myRank,n_proc
      PARAMETER(rootRank=0)

      COMMON /mpivars/ myRank,n_proc

!---  stdio defs
      INTEGER   STDERR,STDIN,STDOUT
      PARAMETER(STDERR=0,STDIN=5,STDOUT=6)

!---  input parameters
!           eps                  softening parameter (can be 0)
!           t_end                end time of simulation
!           dt_disk              output interval for snapshots
!           dt_contr             output interval for energy etc.
!           dt_bh                output interval for BH data
!           dt_timing            output interval for timing info
!           dt_restart           output interval for restart file
!           dt_max               maximum time step
!           eta_s                intital time step parameter
!           eta                  time step parameter
!           irestart             =1 snapshot file is a restart file
!           icmcorr              =1 correct snapshot for CofM if it 
!                                is the first snapshot (diskstep=0)
!           inp_fname            name of initial snapshot file
      INTEGER          irestart,icmcorr
      DOUBLE PRECISION eps, bh_eps, t_end,dt_disk,dt_contr,dt_bh,        &
     &     dt_timing, dt_restart,dt_max,eta_s,eta
      CHARACTER*79     inp_fname

      COMMON /inparsi/ irestart,icmcorr
      COMMON /inpars/  eps,t_end,dt_disk,dt_contr,dt_bh,dt_timing,       &
     &     dt_restart,dt_max,eta_s,eta,bh_eps
      COMMON /inpfnam/ inp_fname

!---  some useful numbers
      DOUBLE PRECISION over6,over2
      PARAMETER       (over6=1.d0/6.d0,over2=0.5d0)


!---  body data (see comments at the end for some explanations)
      INTEGER          ind(NMAX),ind_loc(NMAX_LOC),N,nbh,n_loc,n_act,    &
     &     ind_act(NMAX),ind_act_loc(NMAX_LOC),n_act_loc,                &
     &     n_act_v(NPEMAX),n3_act_v(NPEMAX),myid(NMAX),myind(NMAX),      &
     &     disp_act_v(NPEMAX),disp3_act_v(NPEMAX),inrg
      INTEGER          diskstep,Timesteps
      DOUBLE PRECISION m(NMAX),x(NDIM,NMAX),v(NDIM,NMAX),a(NDIM,NMAX),   &
     &                 adot(NDIM,NMAX),t(NMAX),dt(NMAX),pot(NMAX),MASS,  &
     &                 h(NMAX)
      DOUBLE PRECISION m_loc(NMAX_LOC),x_loc(NDIM,NMAX_LOC),             &
     &                 v_loc(NDIM,NMAX_LOC),a_loc(NDIM,NMAX_LOC),        &
     &                 adot_loc(NDIM,NMAX_LOC),t_loc(NMAX_LOC),          &
     &                 dt_loc(NMAX_LOC),pot_loc(NMAX_LOC),               &
     &                 h_loc(NMAX_LOC)
      DOUBLE PRECISION m_act(NMAX),x_act(NDIM,NMAX_LOC),                 &
     &                 v_act(NDIM,NMAX_LOC),t_act(NMAX_LOC),             &
     &                 dt_act(NMAX_LOC),pot_act(NMAX_LOC),               &
     &                 a_act(NDIM,NMAX_LOC),adot_act(NDIM,NMAX_LOC),     &
     &                 x_act_new(NDIM,NMAX),v_act_new(NDIM,NMAX),        &
     &                 body_pred_loc(2*NDIM,NMAX_LOC),                   &
     &                 body_pred(2*NDIM,NMAX),h_act(NMAX),               &
     &                 body_act_tmp(7,NMAX),body_act_new(7,NMAX)
      DOUBLE PRECISION pot_act_tmp_loc(NMAX),a_act_tmp_loc(NDIM,NMAX),   &
     &                 adot_act_tmp_loc(NDIM,NMAX)
      DOUBLE PRECISION time_cur,eps2,min_t
      DOUBLE PRECISION E_tot,E_pot,E_kin,xcm(NDIM),vcm(NDIM),mom(NDIM),  &
     &                 E_totini

!---  index parameters to identify active particle properties
      INTEGER          iaxnew,iaynew,iaznew,ijxnew,ijynew,ijznew,        &
     &                 ipnew
      PARAMETER       (iaxnew=1,iaynew=2,iaznew=3,ijxnew=4,ijynew=5,     &
     &                 ijznew=6,ipnew=7)

!---  here some arrays are set to common in parallel runs
#ifndef NOMPI
      COMMON /bodati1/ ind,ind_loc
      COMMON /bodat1/  m,x,v,a,adot,t,dt,pot,m_loc,x_loc,v_loc,a_loc,    &
     &                 adot_loc,t_loc,dt_loc,pot_loc,h,h_loc,            &
     &                 body_act_new,body_act_tmp,                        &
     &                 body_pred_loc,body_pred
#endif
#ifdef NOMPI
!---  define global and local vars to be the same which also means 
!---  only the global arrays can be in a COMMON-block
      COMMON /bodati1/ ind
      COMMON /bodat1/  m,x,v,a,adot,t,dt,pot,body_act_new,body_pred,h
      EQUIVALENCE (m_loc,m),(x_loc,x),(v_loc,v),(a_loc,a),               &
     &            (adot_loc,adot),(pot_loc,pot),(t_loc,t),               &
     &            (dt_loc,dt),(ind_loc,ind),                             &
     &            (h_loc,h),                                             &
     &            (body_act_tmp,body_act_new),                           &
     &            (body_pred_loc,body_pred)
#endif

      COMMON /bodati2/ N,nbh,n_loc,n_act_loc,ind_act,                    &
     &                 ind_act_loc,n_act,n_act_v,n3_act_v,disp_act_v,    &
     &                 disp3_act_v,myid,myind,inrg       
      COMMON /bodat2/  m_act,x_act,v_act,                                &
     &                 t_act,dt_act,pot_act,a_act,adot_act,              &
     &                 x_act_new,v_act_new,h_act,                        &
     &                 a_act_tmp_loc,adot_act_tmp_loc,pot_act_tmp_loc
      COMMON /energies/ E_tot,E_pot,E_kin,xcm,vcm,mom,MASS,E_totini


      COMMON /systemi/ diskstep,Timesteps
      COMMON /system/  time_cur,eps2,min_t

 



!---  output
      LOGICAL           lbod,ltot
      COMMON /outpvars/ lbod,ltot
      INTEGER    ubo,utot,utime,ustop,urestf,ubh
      PARAMETER (ubo=12,utot=13,utime=14,ustop=15,urestf=16,ubh=17)
      CHARACTER(LEN=79) totfile,timefile,stopfile,restfile,              &
     &                  bhfile

      PARAMETER (totfile='energy.dat',timefile='timing.dat',             &
     &           bhfile='bh.dat',                                        &
     &           stopfile='.stop',restfile='phiGRAPE.boa')


!---  time steps 
      DOUBLE PRECISION dt_min,t_timing,t_disk,t_contr,t_bh,              &
     &     t_restart,n_act_sum,n_act_sum2,nas_old,                       &
     &     n_act_step

      COMMON /times/   dt_min,t_timing,t_disk,t_contr,t_bh,              &
     &                 t_restart,n_act_sum,n_act_sum2,nas_old,           &
     &                 n_act_step

!---  GRAPE vars
      INTEGER         clusterid,npipe

      COMMON /g6varsi/ clusterid,npipe


!---  timing
      DOUBLE PRECISION cpu0,cpu_total,cpu_step,cpu_act

      COMMON /timing/  cpu0,cpu_total,cpu_step,cpu_act


!---  collisions in MUSE
      INTEGER          icollision,idprimary,idsecondary
      COMMON /colli/   icollision,idprimary,idsecondary

!--------------------------------------------------------------------
!---  body data arrays explained
!---  (NOTE: NMAX is the maximum total particle number;
!---         NMAX_LOC is the maximum per node, usually limited by 
!---         the GRAPE card memory)
!--------------------------------------------------------------------

!---  a) global body data arrays
!---       * these arrays are maintained on root only
!---       * used for data I/O and computation of energies etc.
!---       * syncronized only when needed for I/O or computation
!---
!---  arrays are:    m(1..NMAX)                mass 
!---                 x(3,1..NMAX)              position 
!---                 v(3,1..NMAX)              velocity
!---                 a(3,1..NMAX)              force
!---                 adot(3,1..NMAX)           force derivative
!---                 pot(1..NMAX)              potential
!---                 t(1..NMAX)                time
!---                 dt(1..NMAX)               time step
!---                 ind(1..NMAX)              global index starting at 1


!---  b) local body data arrays
!---       * these are the particles local on a given note
!---       * initially set by root via MPI_Scatter
!---       * updated locally every time step in update_loc_p
!---       * accessed by root via MPI_Gather to syncronize global arrays
!---
!---  arrays are:    m_loc(1..NMAX_LOC)        mass 
!---                 x_loc(3,1..NMAX_LOC)      position 
!---                 v_loc(3,1..NMAX_LOC)      velocity
!---                 a_loc(3,1..NMAX_LOC)      force
!---                 adot_loc(3,1..NMAX_LOC)   force derivative
!---                 pot_loc(1..NMAX_LOC)      potential
!---                 t_loc(1..NMAX_LOC)        time
!---                 dt_loc(1..NMAX_LOC)       time step
!---                 ind_loc(1..NMAX_LOC)      local index starting at 1


!---  c) local active body data arrays
!---       * these are local particles that require a force update, ie the active particles
!---       * selected after finding the global minimum time step
!---       * data is copied from up-to-date local arrays
!---       * used to get predicted positions and velocities
!---
!---  arrays are:    m_act(1..NMAX_LOC)        mass 
!---                 x_act(3,1..NMAX_LOC)      position 
!---                 v_act(3,1..NMAX_LOC)      velocity
!---                 a_act(3,1..NMAX_LOC)      force
!---                 adot_act(3,1..NMAX_LOC)   force derivative
!---                 pot_act(1..NMAX_LOC)      potential
!---                 t_act(1..NMAX_LOC)        time
!---                 dt_act(1..NMAX_LOC)       time step
!---                 ind_act_loc(1..NMAX_LOC)  local index 


!---  d) predicted positions and velocities
!---       * used to store predictions as computed from local active body data arrays 
!---       * in allgather_ap the predictions for local particles are distributed
!---         via MPI_allgatherv to have all predictions on all nodes
!---       * similar a global index array is distributed to all nodes
!--- 
!---  arrays are:    body_pred_loc(1..3,NMAX_LOC)  local predicted positions
!---                 body_pred_loc(4..6,NMAX_LOC)  local predicted velocities
!---                 body_pred(1..3,NMAX)          global predicted positions
!---                 body_pred(4..6,NMAX)          global predicted velocities
!---                 ind_act(1..NMAX)              global  index 


!---  e) updated forces etc.
!---       * used to store the force etc as computed on the GRAPE
!---       * particial local forces are summed up in sumforces
!---       * particial forces are also stored locally to provide a good guess 
!---         for next force update 
!--- 
!---  arrays are:    body_act_tmp(1..3,NMAX)       particial forces
!---                 body_act_tmp(4..6,NMAX)       particial force derivatives
!---                 body_act_tmp(7,NMAX)          particial potentials
!---                 a_act_tmp_loc(3,NMAX_LOC)     saving particial forces
!---                 adot_act_tmp_loc(3,NMAX_LOC)  saving particial force derivatives
!---                 pot_act_tmp_loc(NMAX_LOC)     saving particial potentials
!---                 body_act_new(1..3,NMAX)       total forces
!---                 body_act_new(4..6,NMAX)       total force derivatives
!---                 body_act_new(7,NMAX)          total potentials


!---  f) memory requirements
!---       * DP = DOUBLE PRECISION (8byte); INT = INTEGER (4byte)
!---         ie. 2 INT = 1 DP
!---
!---        global body data:   (16 DP + 1 INT) per particle
!---        local  body data:   (16 DP + 1 INT) per local particle
!---        active body data:   (16 DP + 1 INT) per local particle
!---        predictor:          (6 DP) per local particle + (6 DP + 1 INT) per particle
!---        forces:             (21 DP) per particle
!---                            --------------------------------------------------------
!---        total               (38 DP + 2 INT) per local paricle
!---                           +(43 DP + 2 INT) per global particle
!---                           = 39 DP per local particle + 44 DP per global particle
!---                            ======================================================== 
!---
!---        example:            RIT cluster 32 nodes with 128k particles (2^17) per node;
!---                            ie. a total of 4M particles (2^22)
!---                           = 39*8byte*2^17 + 44*8byte*2^22
!---                           = 39*2^20 byte  + 1408*2^20 byte
!---                           = 1447 Mbyte (+ some overhead!)
!---                             ==========

!--- (NOTE: the maximum array sizes are given, only part of the actual array may be used)
