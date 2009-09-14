      PROGRAM phiGRAPE
C=========================================================================
C===
C===           Parallel Hermite Integrator using a GRAPE hardware 
C===                         for N-body simulations
C===           --------------------------------------------------
C===                      written by S.Harfst (Aug 2006)
C===                          (harfst@astro.rit.edu)
C===
C===  Details:   the program uses a Hermite integration scheme to 
C===  --------   advance a N-body system in time. Forces are calculated
C===             using GRAPE hardware. Particles have individual time 
C===             steps. For more details see astro-ph/0608125. Also
C===             check for READMEs on how to run the code. Variables
C===             are explained in include files.
C===
C===
C===  Acknowlegdements: P.Berczik has contributed to the GRAPE 
C===  ---------------   implementation as well as the parallelisation.
C===      
C==========================================================================
      IMPLICIT NONE
#include "general.inc"
#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'
      INTEGER ierr
#endif

C---  local vars
      INTEGER istop
      DOUBLE PRECISION gflops_total,gflops_step
      LOGICAL          lsyncmode
C---  functions
      DOUBLE PRECISION cputime

C---  set some vars
      inrg       = 0
      istop      = 0
      Timesteps  = 0
      n_act_sum  = 0.d0
      n_act_sum2 = 0.d0
      cpu_total  = 0.d0
      lsyncmode  = .FALSE.

C--------------------------------------------------------------------
C---  initialise MPI
C--------------------------------------------------------------------
      CALL startmpi
      CALL initmpi
C          -------

C--------------------------------------------------------------------
C---  initialise the GRAPE
C--------------------------------------------------------------------
      CALL initgrape
C          ---------

C--------------------------------------------------------------------
C---  read parameter file & body data
C--------------------------------------------------------------------
      IF (myRank .EQ. rootRank) THEN

         CALL readpars
C             --------
      ENDIF

C---------------------------------------------------------------
C---  distribute input data
C---------------------------------------------------------------
      CALL sendinpdata
C          -----------

C---------------------------------------------------------------
C---  get body data (and more if restarting)
C---------------------------------------------------------------
      IF (irestart.NE.1) THEN 
         
         CALL readbodydata
C             ------------
      ELSE
         CALL readrestartfile
C             ---------------
      ENDIF


      IF (myRank .EQ. rootRank) THEN

C---  do center of mass correction if wanted
         IF ((icmcorr.EQ.1) .AND. (diskstep.EQ.0)) CALL cmcorr(N,m,x,v)
C                                                       ------

         OPEN(ustop,FILE=stopfile)   ! creating stopfile
         WRITE(ustop,*) 0
         CLOSE(ustop)

      ENDIF


C--------------------------------------------------------------------
C---  initialise some parameters
C--------------------------------------------------------------------
      CALL initpars
C          --------

C---------------------------------------------------------------
C---  distribute body data
C---------------------------------------------------------------
      CALL sendbodydata
C          ------------


C---  some output with parameters
#ifndef SILENT         
      IF (myRank.EQ.rootRank) THEN

         WRITE(STDOUT,*) '============================================',
     $                   '============='
         IF (irestart.NE.1) THEN
            WRITE(STDOUT,*) '===          Beginning N-body run with ph',
     $                   'iGRAPE       ==='
         ELSE
            WRITE(STDOUT,*) '===         Restarting N-body run with ph',
     $                   'iGRAPE       ==='
         ENDIF
         WRITE(STDOUT,600) N,eps
         WRITE(STDOUT,602) time_cur,t_end
         WRITE(STDOUT,604) dt_disk, dt_contr 
         WRITE(STDOUT,605) dt_bh, dt_timing
         WRITE(STDOUT,606) eta_s, eta
         WRITE(STDOUT,607) irestart,icmcorr
         WRITE(STDOUT,608) n_proc
         WRITE(STDOUT,*) '============================================',
     $                   '============='
 600     FORMAT(1P,' ===  N       = ',I12,
     $        '  eps      = ',E12.4,'  ===')
 602     FORMAT(1P,' ===  t_start = ',E12.4,
     $        '  t_end    = ',E12.4,'  ===')
 604     FORMAT(1P,' ===  dt_disk = ',E12.4,
     $        '  dt_contr = ',E12.4,'  ===')
 605     FORMAT(1P,' ===  dt_bh   = ',E12.4,
     $        '  dt_timing= ',E12.4,'  ===')
 606     FORMAT(1P,' ===  eta_s   = ',E12.4,
     $        '  eta      = ',E12.4,'  ===')
 607     FORMAT(1P,' ===  irestart= ',I12,  '  icmcorr  = ',I12,'  ===')
 608     FORMAT(1P,' ===  n_proc  = ',I12,26x,' ===')

         CALL flush(STDOUT)
      ENDIF
#endif

C--------------------------------------------------------------------
C---  opening files
C--------------------------------------------------------------------
      IF (myRank.EQ.rootRank) THEN
         IF (irestart.NE.1) THEN
            OPEN(UNIT=utot,FILE=totfile)                 ! energies etc.
            IF ((nbh.GT.0).AND.(dt_bh.LE.t_end))
     $           OPEN(UNIT=ubh,FILE=bhfile)              ! BH data
         ELSE
            OPEN(UNIT=utot,FILE=totfile,ACCESS='append') ! energies etc.
            IF ((nbh.GT.0).AND.(dt_bh.LE.t_end))
     $           OPEN(UNIT=ubh,FILE=bhfile,ACCESS='append') ! BH data
         ENDIF
      ENDIF

C--------------------------------------------------------------------
C---  Beginning to initialise the N-body system
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C---  send bodydata to GRAPE
C--------------------------------------------------------------------
      CALL sendbodies2grape
C          ----------------

      IF (irestart.EQ.1) GOTO 95   ! skip the initial force calculation 
                                   ! when restarting

C--------------------------------------------------------------------
C---  select active particles from local particles
C--------------------------------------------------------------------
      min_t = time_cur + dt_min    ! initially all particles have 
      CALL selectactive            ! same time step dt_min
C          ------------

C--------------------------------------------------------------------
C---  predictor (for local active particles)
C--------------------------------------------------------------------
      CALL predictor(1)
C          ---------         

C--------------------------------------------------------------------
C---  gather active particles on all PE's
C--------------------------------------------------------------------
      CALL allgather_ap
C          ------------

C--------------------------------------------------------------------
C---  get local forces on all PE's
C---  (as this is the first call set ifirst to 1 in this call)
C--------------------------------------------------------------------
      CALL gravity(time_cur,1)
C          -------

C--------------------------------------------------------------------
C---  sum up the partial forces
C--------------------------------------------------------------------
      CALL sumforces(1)    ! call with 1 bc all forces etc.  need initialisation
C          ---------

C--------------------------------------------------------------------
C---  get individual time steps (not when restarting)
C--------------------------------------------------------------------
      IF (irestart.NE.1) CALL timestep(0)
C                             --------

C--------------------------------------------------------------------
C---  update GRAPE
C--------------------------------------------------------------------
      CALL update_grape(1)
C          ------------

C--------------------------------------------------------------------
C---  get energy, center of mass and momentum
C--------------------------------------------------------------------
      IF (myRank.EQ.0) THEN
         CALL energy
C             ------
#ifndef SILENT
         WRITE(STDOUT,500) time_cur,E_tot,E_pot,E_kin
 500     FORMAT('time =',F10.6,1P,' E_tot =', E14.6,
     $          ' E_pot =', E14.6,' E_kin =', E14.6)
         CALL flush(STDOUT)
#endif
      ENDIF

C--------------------------------------------------------------------
C---  initial output
C--------------------------------------------------------------------
      IF ((myRank.EQ.0).AND.(diskstep.EQ.0)) THEN
         lbod = .TRUE.          ! write snapshot
         ltot = .TRUE.          ! and energy data
         CALL output(time_cur)
C             ------

         CALL writebh(time_cur,0) ! write BH data, no prediction needed
C             -------

      ENDIF

 95   CONTINUE                   ! skip until here when restarting

C---  start timing her
      cpu_total = 0.d0
      IF (myRank.EQ.rootRank) THEN
         cpu0 = cputime() - cpu_total
         IF (irestart.NE.1) THEN
            OPEN(utime,FILE=timefile)
         ELSE
            OPEN(utime,FILE=timefile,ACCESS='append')
         ENDIF
      ENDIF


C--------------------------------------------------------------------
C---  main loop begins here
C--------------------------------------------------------------------
 100  IF (time_cur.GE.t_end) GOTO 900

C---------------------------------------------------------------
C---  get the global minimum timestep 
C---------------------------------------------------------------
         CALL get_min_t(n_loc,t_loc,dt_loc, min_t, myRank)
C             ---------

C--------------------------------------------------------------------
C---  select active particles from local particles
C--------------------------------------------------------------------
         CALL selectactive
C             ------------

C--------------------------------------------------------------------
C---  predictor (for local active particles)
C--------------------------------------------------------------------
         CALL predictor(0)
C             ---------         

C--------------------------------------------------------------------
C---  gather active particles on all PE's
C--------------------------------------------------------------------
         CALL allgather_ap
C             ------------

C--------------------------------------------------------------------
C---  get local forces on all PE's
C---  (as this is not the first call set ifirst to 0 in this call)
C--------------------------------------------------------------------
         CALL gravity(min_t,0)
C             -------

C--------------------------------------------------------------------
C---  sum up the partial forces
C--------------------------------------------------------------------
         CALL sumforces(0)      ! call with 0; just sum forces
C             ---------

C--------------------------------------------------------------------
C---  make the corrector step
C---  this also updates the time steps of active particles
C--------------------------------------------------------------------
         CALL corrector
C             ---------

C--------------------------------------------------------------------
C---  update local particles from local active particles
C--------------------------------------------------------------------
         CALL update_loc_p
C             ------------

C--------------------------------------------------------------------
C---  also update local particles on GRAPE from local active particles
C--------------------------------------------------------------------
         CALL update_grape(0)
C             ------------



C--------------------------------------------------------------------
C---  synchronize system if flag is set
C--------------------------------------------------------------------
         IF (lsyncmode) THEN
            
            CALL selectsync(min_t)
C                ----------
            CALL sync
C                ----
         ENDIF


C---  write BH data if needed
         IF ((myRank.EQ.rootRank).AND.(min_t.GE.t_bh)) THEN 
            CALL writebh(min_t,1)   ! use prediction
C                -------
         ENDIF


C---  if output needs to be done collect local particles on root
C---  also compute energy (change needed when output seperate)
         IF ((min_t.GE.t_disk) .OR. (min_t.GE.t_contr)) THEN

#ifndef NOMPI
            CALL MPI_Gather(m_loc,   n_loc, MPI_DOUBLE, m,   n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(t_loc,   n_loc, MPI_DOUBLE, t,   n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(dt_loc,  n_loc, MPI_DOUBLE, dt,  n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(x_loc, 3*n_loc, MPI_DOUBLE, x, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(v_loc, 3*n_loc, MPI_DOUBLE, v, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(pot_loc, n_loc, MPI_DOUBLE, pot, n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(a_loc, 3*n_loc, MPI_DOUBLE, a, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(adot_loc, 3*n_loc, MPI_DOUBLE,adot, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
#endif


C---  get energy and make output at root
            IF (myRank.EQ.rootRank) THEN
               CALL energy
C                   ------
#ifndef SILENT
               WRITE(STDOUT,500) min_t,E_tot,E_pot,E_kin
               CALL flush(STDOUT)
#endif
               lbod = (min_t.GE.t_disk)
               ltot = (min_t.GE.t_contr)
               IF (lbod) diskstep = diskstep + 1
               CALL output(min_t)
C                   ------
            ENDIF

C---  check for external stop
            OPEN(ustop,FILE=stopfile)
            READ(ustop,*) istop
            CLOSE(ustop)

         ENDIF

C---  update time etc
         time_cur   = min_t
         Timesteps  = Timesteps  + 1
         n_act_sum  = n_act_sum  + DBLE(n_act)
         n_act_sum2 = n_act_sum2 + DBLE(n_act + 48 - MOD(n_act,48))

         IF (time_cur.GE.t_disk)   t_disk  = t_disk + dt_disk
         IF (time_cur.GE.t_contr)  t_contr = t_contr + dt_contr
         IF (time_cur.GE.t_bh)     t_bh    = t_bh + dt_bh


C---  output timing information
         IF (myRank.EQ.rootRank) THEN
            IF (time_cur.GE.t_timing) THEN
               n_act_step   = n_act_sum - nas_old
               nas_old      = n_act_sum
               cpu_act      = cputime() - cpu0
               cpu_step     = cpu_act - cpu_total
               cpu_total    = cpu_act
               gflops_total = 57.d0*N*n_act_sum/cpu_total*1d-9
               gflops_step  = 57.d0*N*n_act_step/cpu_step*1d-9

               WRITE(utime,510) Timesteps,time_cur,cpu_step,cpu_total,
     $                                        n_act_sum,n_act_sum2,
     $                                        gflops_total,gflops_step
 510           FORMAT(I12,1p,7E14.6)
               CALL flush(utime)

               t_timing = t_timing + dt_timing
            ENDIF
         ENDIF

C---  write restart file and stop run
         IF (istop.GE.1) THEN
            WRITE(STDERR,*) 'forcing stop at t =',time_cur
            IF (time_cur.GE.t_restart) t_restart = t_restart+dt_restart

            CALL writerestartfile
C                ----------------

            CALL g6_close(clusterid)
C                --------
#ifndef NOMPI
            CALL mpi_finalize(ierr)
C                ------------
#endif
            STOP
         ENDIF
         
C---  write a regular restart file but do not stop run
         IF (time_cur.GE.t_restart) THEN
            t_restart = t_restart + dt_restart

            CALL writerestartfile
C                ----------------
         ENDIF


         GOTO 100               ! go back to do next time step

 900  CONTINUE                  ! end of main loop: t_end reached

      CALL g6_close(clusterid)
C          --------

C--------------------------------------------------------------------
C---  stop MPI and end program
C--------------------------------------------------------------------
#ifndef NOMPI
      CALL mpi_finalize(ierr)
C          ------------
#endif

      STOP
      END
