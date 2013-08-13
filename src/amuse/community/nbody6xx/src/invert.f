      SUBROUTINE INVERT(DT,DTAU)
*
*
*       Inversion of physical time interval.
*       ------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      LOGICAL  KSLOW,KCOLL
      REAL*8  VI(NMX3),VC(NMX3),KSCH,C(5)
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW,KCOLL
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CALLS/  TPR,TKPR,STEP,IDER,ICALL,NFN,NREG,ITER,IMCIRC
      COMMON/CLUMP/   BODYS(NMX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NMX,5),ISYS(5)
      COMMON/KSAVE/  K1,K2
      SAVE
      DATA  IWARN/0/
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*	Transform to current chain coordinates and set RINV.
      DO I=1,N-1
          L1=3*(I-1)+1
          KS1=4*(I-1)+1
          CALL KSPHYS(Q(KS1),P(KS1),XC(L1),WC(L1))
          RK = Q(KS1)**2 + Q(KS1+1)**2 + Q(KS1+2)**2 + Q(KS1+3)**2
          RINV(I) = 1.0/RK
      END DO
*
*       Obtain physical velocities (first absolute, then relative).
      L = 3*(N-2)
      DO K = 1,3
          VI(K) = -WC(K)/MC(1)
          VI(L+K+3) = WC(L+K)/MC(N)
      END DO
      DO I = 2,N-1
          L = 3*(I-1)
          DO K = 1,3
              VI(L+K) = (WC(L+K-3) - WC(L+K))/MC(I)
          END DO
      END DO  
      do j = 1,3*(n-1)
          vc(j) = vi(j+3) - vi(j)
      end do
*
*       Determine the smallest two-body distance and chain index.
      rm = 0.0
      do i = 1,n-1
          if (rinv(i).gt.rm) then
              rm = rinv(i)
              i1 = i
          end if
      end do
      rm = 1.0/rm
*
*       Form semi-major axes for the two closest distances.
      dm = 0.0
      do i = 1,n-1
          L = 3*(i - 1)
          if (i.eq.i1) then
              w2 = vc(L+1)**2 + vc(L+2)**2 + vc(L+3)**2
              mb = mc(i) + mc(i+1)
              amax = 2.0*rinv(i) - w2/mb
              a1 = 1.0/amax
              V21 = w2
*       Obtain the largest relative perturbation from nearest neighbour.
              if (i.eq.1) then
                  g1 = 2.0*mc(3)/mb*(RINV(2)/RINV(1))**3
              else
                  g1 = 2.0*mc(i-1)/mb*(RINV(i-1)/RINV(i))**3
                  if (i.lt.n-1) then
                      g1b = 2.0*mc(i+2)/mb*(RINV(i+1)/RINV(i))**3
                      g1 = max(g1b,g1)
                  end if
              end if
*       Treat the second smallest distance similarly.
          else if (rinv(i).gt.dm) then
              w2 = vc(L+1)**2 + vc(L+2)**2 + vc(L+3)**2
              mb = mc(i) + mc(i+1)
              amax = 2.0*rinv(i) - w2/mb
              a2 = 1.0/amax
              i2 = i
              dm = rinv(i)
              if (i.eq.1) then
                  g2 = 2.0*mc(3)/mb*(RINV(2)/RINV(1))**3
              else
                  g2 = 2.0*mc(i-1)/mb*(RINV(i-1)/RINV(i))**3
                  if (i.lt.n-1) then
                      g2b = 2.0*mc(i+2)/mb*(RINV(i+1)/RINV(i))**3
                      g2 = max(g2b,g2)
                  end if
              end if
          end if
      end do
      dm = 1.0/dm
*
*       Obtain reliable semi-major axis for small pericentre or large EB.
      EB = -0.5*mc(i1)*mc(i1+1)/a1
      EB1 = EB/ENERGY
      IF (rm.lt.0.2*a1.or.EB1.gt.0.9) then
*       Copy current configuration for EREL & TRANSK.
          DO 2 I = 1,N-1
              KS = 4*(I - 1)
              DO 1 J = 1,4
                  QK(KS+J) = Q(KS+J)
                  PK(KS+J) = P(KS+J)
    1         CONTINUE
    2     CONTINUE
*
*       Evaluate semi-major axis of K1 & K2 from non-singular variables.
          K1 = INAME(i1)
          K2 = INAME(i1+1)
          CALL EREL(i1,EB,SEMI)
*
          ERR = (a1 - SEMI)/SEMI
          IF (ABS(ERR).GT.1.0E-05.AND.IWARN.LT.10) THEN
              IWARN = IWARN + 1
              ZK = 1.0
              DO I = 1,N-1
                  ZK = MAX(KSCH(I),ZK)
              END DO
              if(rank.eq.0)WRITE (6,4)  NSTEP1, ZK, g1, RM, SEMI, ERR
    4         FORMAT (' WARNING!    INVERT    # K g1 RM SEMI DA/A ',
     &                                        I6,F7.2,1P,4E10.2)
          END IF
*       Replace direct value by improved determination.
          a1 = SEMI
      END IF
*
*       Restrict the interval to a fractional period or crossing time.
      IF (a1.GT.0.0) THEN
          TKS = 6.28*KSCH(i1)*a1*SQRT(a1/(mc(i1) + mc(i1+1)))
*       Adopt upper limit of TKS/2 but allow genuine small DT.
          DT = MIN(0.5*TKS,DT)
      ELSE
*       Use crossing time for hyperbolic velocity (TPR compensates small R).
          DT = MIN(rm/SQRT(V21),DT)
      END IF
*
*       Check second closest pair if bound (small period; a1 < 0 possible).
      if (a2.gt.0.0) then
          TKS = 6.28*KSCH(i2)*a2*SQRT(a2/(mc(i2) + mc(i2+1)))
          DT = MIN(0.5*TKS,DT)
      end if
*
*       Initialize time variables.
      DT0 = DT
      SUM1 = 1.0/TPR
      SUM2 = 0.0
*
*       Skip iteration for significant perturbations.
      IF (g1.GT.0.01) THEN
          DTAU = SUM1*DT
          GO TO 30
      END IF
*
*       Begin the loop (second time depends on conditions).
    5 DT = DT0
      MB = mc(i1) + mc(i1+1)
*       Form useful quantities (including EB and eccentricity if needed).
      BETA = MB/a1
      LI = 3*(i1 - 1)
      ETA = xc(LI+1)*vc(LI+1) + xc(LI+2)*vc(LI+2) + xc(LI+3)*vc(LI+3)
      ZETA = MB - BETA*RM
*     EB = -0.5*MC(I1)*MC(I1+1)/a1
*     ECC = SQRT((1.0 - RM/a1)**2 + ETA**2/(MB*a1))
*
*       Sum the constant part of the integral.
      SUM1 = SUM1 - 2.0*MC(i1)*MC(I1+1)/(KSCH(i1)*RM)
*       Reduce physical time by slow-down factor for the two-body integral.
      DT = DT/KSCH(i1)
*
*       Obtain improved initial guess by factor of 2 bisecting (Chris Tout).
      IT = 0
*       Choose a conservative starting value in case of small RM.
      Y = DT*MIN(1.0/RM,0.5/ABS(a1))
    8 Z = BETA*Y**2
      CALL CFUNCS(Z,C)
      Y0 = DT - ((ZETA*C(3)*Y + ETA*C(2))*Y + RM)*Y
      IT = IT + 1
      IF (Y0.GT.0.0) THEN
          Y = 2.0*Y
      ELSE
          Y = 0.5*Y
      END IF
      IF (IT.LT.4) GO TO 8
*
*       Perform Newton-Raphson iteration using first and second derivatives.
   10 Z = BETA*Y**2
*
*       Obtain new c-functions for argument Z = -2h*Y^2.
      CALL CFUNCS(Z,C)
*
*       Define the function Y0 and its derivatives (cf. Book eqns. 12.9).
      Y0 = DT - ((ZETA*C(3)*Y + ETA*C(2))*Y + RM)*Y
      YPR = -((ZETA*C(2)*Y + ETA*C(1))*Y + RM)
      Y2PR = -ETA*(1.0 - Z*C(2)) - ZETA*Y*C(1)
*       Adopt safety measure to avoid negative argument (Seppo Mikkola).
      D = Y0*Y2PR/YPR**2
      Y = Y - Y0/YPR/SQRT(0.01 + ABS(0.99D0 - D))
      DT1 = ((ZETA*C(3)*Y + ETA*C(2))*Y + RM)*Y
      IT = IT + 1
      IF (ABS(DT - DT1).GT.1.0D-06*ABS(DT).AND.IT.LT.25) GO TO 10
*
*       Accumulate the integral of dt/R.
      SUM2 = SUM2 + 2.0*MC(i1)*MC(i1+1)*Y
*
      IF (IT.GT.15.AND.IWARN.LT.50) THEN
          IWARN = IWARN + 1
          ECC = SQRT((1.0 - RM/a1)**2 + ETA**2/(MB*a1))
          if(rank.eq.0)WRITE (6,15) IT,KSCH(I1),ECC,RM,a1,DT-DT1,DTAU
   15     FORMAT (' WARNING!    INVERT   IT K E R A ERR DTAU ',
     &                                   I4,F7.2,F7.3,1P,5E9.1)
      END IF
*
*       Consider the second close interaction for small perturbation.
      IF (i2.ne.i1.and.g2.lt.0.01) THEN
          i1 = i2
          a1 = a2
          rm = dm
          GO TO 5
      END IF
*
*       Combine the two terms and ensure positive value.
      DTAU = SUM1*DT0 + SUM2
      IF (DTAU.LT.0.0) THEN
          DTAU = ABS(SUM1)*DT0
          IF (IWARN.LT.100) THEN
              IWARN = IWARN + 1
              if(rank.eq.0)WRITE (6,16)KSCH(i1),RM/A1,SUM1*DT0,SUM2,DTAU
   16         FORMAT (' INVERT NEGATIVE!    KSCH R/A S1*DT0 S2 DTAU ',
     &                                      F7.2,1P,4E10.2)
          END IF
      END IF
*
*     if(rank.eq.0)WRITE (6,20)  IT, a1, DT, DT-DT1, DTAU, STEP
*  20 FORMAT (' INVERT:    IT A DT ERR DTAU STEP ',I4,1P,5E10.2)
*
   30 RETURN
*
      END
