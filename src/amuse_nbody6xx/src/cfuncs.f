      SUBROUTINE CFUNCS(Z,C)
*
*
*       C-functions for two-body iteration.
*       -----------------------------------
*
*       Coded by Seppo Mikkola (11/99).
*       Reference: New Astron. 3, 309, 1998.
*       ------------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(C6=1.D0/6.D0,C132=1.D0/132.D0,C56=1.D0/56.D0,
     &   C30=1.D0/30D0,C24=1.D0/24.D0,C156=1.D0/156.D0,
     &   C90=1.D0/90.D0,C110=1.D0/110.D0,C16=1.D0/16.D0,C8=1.D0/8.D0,
     &   C72=1.D0/72.D0,C42=1.D0/42.D0,C120=1.D0/120.D0,U=1.d0)
      REAL*8 C(5)
      data large/0/
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*       Reduce the argument by 4 until < 0.07.
      h=z
      DO 1 k=0,20
          KK = K
          IF (ABS(H).lt. 0.07) go to 5
          h=0.25d0*h
    1 CONTINUE
*
      if (rank.eq.0.and.large.lt.100) then
          large=large+1
          WRITE (6,3) large, k, z, h
    3     FORMAT (' WARNING!    CFUNCS    # k Z h   ',2I4,1P,2E10.2)
      end if
*
*       Form the two highest c-functions.
    5 C(4)=(U-H*(U-H*(U-H*C90/(U+H*C132))*C56)*C30)*C24
      C(5)=(U-H*(U-H*(U-H*C110/(U+H*C156))*C72)*C42)*C120
*
*       Obtain the lower c-function by recursion.
      DO 4 I=1,KK
          C3=c6-h*C(5)
          C2=.5D0-h*C(4)
          C(5)=(C(5)+C(4)+C2*C3)*c16
          C(4)=C3*(2.d0-h*C3)*c8
          h=4.d0*h
    4 CONTINUE
*
      C(3)=c6-Z*C(5)
      C(2)=.5D0-Z*C(4)
      C(1)=1.d0-Z*C(3)
*
      RETURN
*
      END
