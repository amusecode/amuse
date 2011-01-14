      SUBROUTINE initpars
C=========================================================================
C===  
C===  task:  initialise some parameters
C===  -----
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE
#include "general.inc"

      INTEGER i

C---  softening
      eps2   = eps*eps

C---  set minimum and maximum time step
      dt_min = 2.d0**ndtmin   ! set ndtmin in paras.inc

      dt_max = MIN(dt_disk,   dt_max)   ! dt_max is limited by output interval
      dt_max = MIN(dt_contr,  dt_max)   ! to make sure system is in sync at output

      dt_restart = MAX(dt_restart, 2.d0*dt_disk) ! make sure restart file is not written
                                                 ! more frequently than very other snapshot

C---  set individual particle times and time steps
      IF (irestart.NE.1) THEN
         DO i=1,N
            t(i)   = time_cur
            dt(i)  = dt_min
         ENDDO
      ENDIF

      n_loc = N/n_proc          ! set the local number of particles
                                ! n_proc is set to 1 in serial runs

C---  set times for output
      IF (irestart.NE.1) THEN
         t_timing  = time_cur + dt_timing
         t_disk    = time_cur + dt_disk 
         t_contr   = time_cur + dt_contr
         t_bh      = time_cur + dt_bh
         t_restart = time_cur + dt_restart

         Timesteps = 0

      ENDIF

      RETURN
      END
