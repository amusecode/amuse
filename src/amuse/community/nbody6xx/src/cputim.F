      SUBROUTINE CPUTIM(TCOMP)
*
*
*       CPU time.
*       ---------
*
      INCLUDE 'common6.h'
*      REAL*8  TCOMP
      COMMON/ICPU0/ ICPU
      REAL*4  TARRAY(2),etime
*
*       Initialize timer (first call) or obtain elapsed time.
*         TCOMP = FLOAT(ITIME)/6000.0
*         TCOMP = REAL(IRTC())*6.67D-9/60.

#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif

#if MPIINIT
          TCOMP = MPI_WTIME()/60.
#else
          TCOMP = ETIME(TARRAY)/60.0
#endif
*         TCOMP = MCLOCK()/6000.
*       Elapsed CPU time in minutes on VAX, SUN or MIPS & IBM RS/6000.
*       and T3D
*
      RETURN
*
      END
