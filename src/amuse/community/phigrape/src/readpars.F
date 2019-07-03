      SUBROUTINE readpars
C=========================================================================
C===
C===  task:  read parameter file
C===  -----  
C===
C===  parameter:   eps                  softening parameter (can be 0)
C===  ----------   t_end                end time of simulation
C===               dt_disk              output interval for snapshots
C===               dt_contr             output interval for energy etc.
C===               dt_bh                output interval for BH data
C===               dt_timing            output interval for timing info
C===               dt_max               maximum time step
C===               eta_s                intital time step parameter
C===               eta                  time step parameter
C===               irestart             =1 snapshot file is a restart file
C===               icmcorr              =1 correct snapshot for CofM if it 
C===                                    is the first snapshot (diskstep=0)
C===               inp_fname            name of initial snapshot file
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE
#include "general.inc"

      OPEN(10,FILE='phiGRAPE.inp')
      
      READ(10,*) eps 
      READ(10,*) t_end
      READ(10,*) dt_disk
      READ(10,*) dt_contr
      READ(10,*) dt_bh
      READ(10,*) dt_timing
      READ(10,*) dt_restart
      READ(10,*) dt_max
      READ(10,*) eta_s
      READ(10,*) eta
      READ(10,*) irestart,icmcorr
      READ(10,*) nbh
      READ(10,*) inp_fname

      CLOSE(10)

#ifdef DEBUG
         WRITE(STDERR,*) eps, t_end, dt_disk, dt_contr, dt_bh, eta_s,
     $                   eta, nbh, inp_fname
#endif


      RETURN
      END
