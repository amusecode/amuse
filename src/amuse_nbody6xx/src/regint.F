      SUBROUTINE REGINT(I,KLIST)
*
*
*       Regular integration.
*       --------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),DV(3),FD(3),FDR(3)
      REAL*8  FRX(3),FDX(3)
*     REAL*8  DFI(3),DFRG(3)
      INTEGER KLIST(LMAX)
*      INTEGER NLIST(NMAX)
*      LOGICAL LCHECK
*
      DO 6 K = 1,3
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
   6  CONTINUE
*
*       Set neighbour number, time-step & choice of central distance.
      NNB0 = KLIST(1)
      DTR = TIME - T0R(I)
      IRSKIP = 0
      IF (KZ(39).EQ.0) THEN
          RI2 = (XI(1) - RDENS(1))**2 + (XI(2) - RDENS(2))**2 +
     &                              (XI(3) - RDENS(3))**2
          RH2 = RSCALE**2
      ELSE
          RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
          RH2 = 9.0*RSCALE**2
      END IF
*
*       Obtain irregular & regular force and determine current neighbours.
      RS2 = RS(I)**2
*       Take volume between inner and outer radius equal to basic sphere.
*     --12/27/12 20:02-lwang-debug--------------------------------------*
***** Note:change 1.59 to 1.1------------------------------------------**
      RCRIT2 = 1.59*RS2
* 1    RCRIT2 = 1.1*RS2      
*     --12/27/12 20:02-lwang-end-debug----------------------------------*
*       Set radial velocity factor for the outer shell.
      VRFAC = -0.1*RS2/DTR
*       Start count at 2 and subtract 1 at the end to avoid ILIST(NNB+1).
      NNB = 1
*
*       Initialize scalars for forces & derivatives.
      ihost=1
#ifdef GPU
*     --12/27/12 21:45-lwang-debug--------------------------------------*
***** Note: change-----------------------------------------------------**
      if (list_g(1,I).ge.0.AND.I.LE.N) then
         DO 5 K = 1,3
            FIRR(K) = 0.0D0
            FD(K) = 0.0D0
            IF(I.GT.N) THEN
               FREG(K) = 0.0D0
               FDR(K) = 0.0D0
            ELSE
               FREG(K) = FRG(K,I)
               FDR(K) = FRGDOT(K,I)
            END IF
 5       CONTINUE
         PHII(I) = 0.D0
         PHIR(I) = PHIDBL(I)    ! PHIDBL after change to PHIREG!!!
         PHIR1(I) = 0.D0
         NNBGPU = LIST_G(1,I)
c$$$  IF(NNBGPU.LT.0)THEN
c$$$  PRINT*,' Warning Regint t,i,nnbgpu=',TIME,I,NNBGPU
c$$$  L1 = 1
c$$$  DO 7 L = 1,N-IFIRST+1
c$$$  IF(L.NE.I)THEN
c$$$  L1 = L1 + 1
c$$$  NLIST(L1) = L
c$$$  END IF
c$$$  7 CONTINUE
c$$$  NNBG = N - IFIRST
c$$$  ELSE
c$$$  NNBG = NNBGPU
c$$$  END IF
         nnbg=nnbgpu
         ihost=0
c$$$      else
c$$$         write (6,*) 'Done calculation on host'
      end if
c$$$#else
#endif
      if (ihost.eq.1) then
*     --12/27/12 21:45-lwang-end-debug----------------------------------*
         DO 205 K = 1,3
            FIRR(K) = 0.0D0
            FREG(K) = 0.0D0
            FD(K) = 0.0
            FDR(K) = 0.0
 205     CONTINUE
         PHII(I) = 0.D0
         PHIR(I) = 0.D0
         PHIR1(I) = 0.D0
*     --12/27/12 21:48-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      end if
c$$$#endif
*     --12/27/12 21:48-lwang-end----------------------------------------*
*
*       Choose appropriate force loop for single particle or c.m. body.
      IF (I.GT.N) THEN
*     See whether perturbation allows single particle approximation.
         I1 = 2*(I - N) - 1
         IF (LIST(1,I1).GT.0) THEN
*     Obtain total force on c.m. particle.
           CALL CMFREG(I,RS2,RCRIT2,VRFAC,NNB,XI,XIDOT,FIRR,FREG,FD,FDR)
            GO TO 20
         END IF
      END IF
*
      call cputim(tt1)
#ifdef GPU
*       Perform only extended neighbour loop if GPU has been used.
      if(ihost.eq.0) then
         DO 10 L = 2, NNBG + 1
*     --12/28/12 11:02-lwang-debug--------------------------------------*
***** Note:debug-------------------------------------------------------**
c$$$         IF(NNBGPU.LT.0) THEN
c$$$            J = NLIST(L)
c$$$        ELSE
c$$$            J = LIST_G(L,I)
c$$$        END IF
            j=list_g(l,i)
*     --12/28/12 11:02-lwang-end-debug----------------------------------*
*       IF (J.EQ.I) GO TO 10

            A1 = X(1,J) - XI(1)
            A2 = X(2,J) - XI(2)
            A3 = X(3,J) - XI(3)
*     Predicted coordinates avoids spurious force differences.
            DV(1) = XDOT(1,J) - XIDOT(1)
            DV(2) = XDOT(2,J) - XIDOT(2)
            DV(3) = XDOT(3,J) - XIDOT(3)
            RIJ2 = A1*A1 + A2*A2 + A3*A3
            DR2I = 1.0/RIJ2
            DR3I = BODY(J)*DR2I*SQRT(DR2I)
            DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
            DRDP = 3.0*DRDV*DR2I
*     
*       Retain particles with small steps (large derivative corrections).
*     --03/18/13 19:54-lwang-improve------------------------------------*
***** Note:suppress step(j).gt.10.0*smin-------------------------------**
*     See whether the distance exceeds the outer shell radius.
            IF (RIJ2.GT.RCRIT2)THEN
c$$$               PRINT*,
c$$$     *        ' t,i,j,regint warning rij2,rcrit2=',TIME,name(I),name(J)
c$$$     *              ,RIJ2,RCRIT2,' dr=',RIJ2-RCRIT2
c$$$               call flush(6)
              GO TO 8
            END IF
*
           IF (RIJ2.GT.RS2) THEN
               IF (DRDV.GT.VRFAC) GO TO 8
            END IF
*     --03/18/13 19:55-lwang-end-improve--------------------------------*
*     Accept member if maximum penetration factor exceeds 8 per cent.
*     
*     Increase neighbour counter and obtain current irregular force.

            NNB = NNB + 1
            ILIST(NNB) = J
            FIRR(1) = FIRR(1) + A1*DR3I
            FIRR(2) = FIRR(2) + A2*DR3I
            FIRR(3) = FIRR(3) + A3*DR3I
            FD(1)   = FD(1)   + (DV(1) - A1*DRDP)*DR3I
            FD(2)   = FD(2)   + (DV(2) - A2*DRDP)*DR3I
            FD(3)   = FD(3)   + (DV(3) - A3*DRDP)*DR3I

c     Obtain potential.

            PHII(I) = PHII(I) - DR3I*RIJ2
            GO TO 10
c
 8          CONTINUE

            FREG(1) = FREG(1) + A1*DR3I
            FREG(2) = FREG(2) + A2*DR3I
            FREG(3) = FREG(3) + A3*DR3I
            FDR(1)  = FDR(1)  + (DV(1) - A1*DRDP)*DR3I
            FDR(2)  = FDR(2)  + (DV(2) - A2*DRDP)*DR3I
            FDR(3)  = FDR(3)  + (DV(3) - A3*DRDP)*DR3I

c     Obtain potential and derivative.

            PHIR(I)  = PHIR(I)  - DR3I*RIJ2
            PHIR1(I) = PHIR1(I) + DRDV*DR3I

 10      CONTINUE               ! DO 10 L = 2, NNB

*        WRITE(6,1111) I,NNBGPU, NNB, TIME, FREG(1), FDR(1),
*    *   FRG(1,I), FRGDOT(1,I)
*1111    FORMAT(1X,' I, nnbs time fr, fdr, frgpu, frdgpu=',
*    *  3I6,1P,5D12.5)
c--------------------------------------------------------------------
c--------------------------------------------------------------------
*        IF(TIME.GT.0.01D0)STOP
*     --12/28/12 11:04-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      end if         
c$$$#else
#endif
*     --12/28/12 11:04-lwang-end----------------------------------------*
*     Perform fast force loop over single particles.
      if (ihost.eq.1) then
         DO 210 J = IFIRST,N
*RSP
            IF (J.EQ.I) GO TO 210
*RSP
            A1 = X(1,J) - XI(1)
            A2 = X(2,J) - XI(2)
            A3 = X(3,J) - XI(3)
*     Predicted coordinates avoids spurious force differences.
            DV(1) = XDOT(1,J) - XIDOT(1)
            DV(2) = XDOT(2,J) - XIDOT(2)
            DV(3) = XDOT(3,J) - XIDOT(3)
*     
            RIJ2 = A1*A1 + A2*A2 + A3*A3
            DR2I = 1.0/RIJ2
            DR3I = BODY(J)*DR2I*SQRT(DR2I)
            DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
            DRDP = 3.0*DRDV*DR2I
*
*       See whether the distance exceeds the outer shell radius.
*     --03/18/13 20:39-lwang-improve------------------------------------*
***** Note:improve-----------------------------------------------------**
            IF (RIJ2.GT.RCRIT2) GO TO 888
*     
*     Retain particles with small steps (large derivative corrections).
            IF (RIJ2.GT.RS2) THEN
*     Accept member if maximum penetration factor exceeds 8 per cent.
               IF (DRDV.GT.VRFAC) GO TO 888
            END IF
*            IF(RIJ2.GT.RS2) GO TO 888
*     --03/18/13 20:39-lwang-end-improve--------------------------------*
*     
*     Increase neighbour counter and obtain current irregular force.
            NNB = NNB + 1
            ILIST(NNB) = J
            FIRR(1) = FIRR(1) + A1*DR3I
            FIRR(2) = FIRR(2) + A2*DR3I
            FIRR(3) = FIRR(3) + A3*DR3I
            FD(1) = FD(1) + (DV(1) - A1*DRDP)*DR3I
            FD(2) = FD(2) + (DV(2) - A2*DRDP)*DR3I
            FD(3) = FD(3) + (DV(3) - A3*DRDP)*DR3I
*     Obtain potential.
            PHII(I) = PHII(I) - DR3I*RIJ2
            GO TO 210
*
*       Obtain the regular force.
 888        FREG(1) = FREG(1) + A1*DR3I
            FREG(2) = FREG(2) + A2*DR3I
            FREG(3) = FREG(3) + A3*DR3I
            FDR(1) = FDR(1) + (DV(1) - A1*DRDP)*DR3I
            FDR(2) = FDR(2) + (DV(2) - A2*DRDP)*DR3I
            FDR(3) = FDR(3) + (DV(3) - A3*DRDP)*DR3I
*     Obtain potential and derivative.
            PHIR(I) = PHIR(I) - DR3I*RIJ2
            PHIR1(I) = PHIR1(I) + DRDV*DR3I
 210     CONTINUE
*     --12/28/12 11:06-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      end if      
c$$$#endif
*     --12/28/12 11:07-lwang-end----------------------------------------*
      call cputim(tt2)
      ttfrc = ttfrc + (tt2-tt1)*60.
*
*       Add any contributions from regularized c.m. particles.
      IF (NPAIRS.GT.0) THEN
         CALL CMFREG(I,RS2,RCRIT2,VRFAC,NNB,XI,XIDOT,FIRR,FREG,FD,FDR)
      END IF
*
*       Add any contributions from regularized c.m. particles.
*       Include treatment for regularized clump.
      IF (NCH.GT.0) THEN
*       Distinguish between chain c.m. and any other particle.
          IF (NAME(I).EQ.0) THEN
              CALL CHFIRR(I,1,XI,XIDOT,FIRR,FD)
          ELSE
*       Search the chain perturber list for #I.
              NP1 = LISTC(1) + 1
              DO 15 L = 2,NP1
                  J = LISTC(L)
                  IF (J.GT.I) GO TO 20
                  IF (J.EQ.I) THEN
                      CALL FCHAIN(I,1,XI,XIDOT,FIRR,FD)
                      GO TO 20
                  END IF
   15         CONTINUE
          END IF
      END IF
*
*       See whether an external force should be added.
*
   20 CONTINUE
*       Check whether cloud forces should be included.
      IF (KZ(13).GT.0) THEN
          CALL FCLOUD(I,FREG,FDR,1)
      END IF
*
*       See whether an external force should be added.
      IF (KZ(14).GT.0) THEN
*       Save current values for deriving work done by tides (#14 = 3).
          IF (KZ(14).EQ.3) THEN
              DO 22 K = 1,3
                  FRX(K) = FREG(K)
                  FDX(K) = FDR(K)
   22         CONTINUE
          END IF
*
*       Obtain the tidal perturbation (force and first derivative).
          CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,1)
*
*       Form rate of tidal energy change during last regular step.
          IF (KZ(14).EQ.3) THEN
              WDOT = 0.0
              W2DOT = 0.0
*             W3DOT = 0.0
              DO 24 K = 1,3
                  PX = FREG(K) - FRX(K)
                  DPX = FDR(K) - FDX(K)
                  WDOT = WDOT + XIDOT(K)*PX
                  W2DOT = W2DOT + (FREG(K) + FIRR(K))*PX + XIDOT(K)*DPX
*                 W3DOT = W3DOT + 2.0*(FREG(K) + FIRR(K))*DPX +
*    &                            (FDR(K) + FD(K))*PX
   24         CONTINUE
*       Note: second-order term derived by Douglas Heggie (Aug/03).
*       Accumulate tidal energy change for general galactic potential.
*       Note: Taylor series at end of interval with negative argument.
*         ETIDE = ETIDE - BODY(I)*((ONE6*W3DOT*DTR - 0.5*W2DOT)*DTR +
*    &                                                   WDOT)*DTR
          ETIDE = ETIDE + BODY(I)*(0.5*W2DOT*DTR - WDOT)*DTR
*       Note: integral of Taylor series for V*P using final values.
          END IF
      END IF
*
      NNB = NNB - 1
*       Check for zero neighbour number.
      IF (NNB.EQ.0) THEN
*       Assume small mass at centre for distant body or no neighbours.
          R2 = XI(1)**2 + XI(2)**2 + XI(3)**2
          FIJ = 0.01*BODYM/(R2*SQRT(R2))
          RDOT = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                                 XI(3)*XIDOT(3))/R2
          DO 25 K = 1,3
              FIRR(K) = FIRR(K) - FIJ*XI(K)
              FD(K) = FD(K) - (XIDOT(K) - RDOT*XI(K))*FIJ
   25     CONTINUE
*       Modify neighbour sphere gradually but allow encounter detection.
          IF (RDOT.GT.0.0) THEN
              RS(I) = MAX(0.75*RS(I),0.1*RSCALE)
          ELSE
              RS(I) = 1.1*RS(I)
          END IF
          NBVOID = NBVOID + 1
          IRSKIP = 1
*       Skip another full N loop (necessary for GPU/OpenMP case).
          GO TO 50
*         GO TO 1
      END IF
*
*       Restrict neighbour number < NNBMAX to permit one normal addition.
      IF (NNB.LT.NNBMAX) GO TO 40
*
*       Reduce search radius by cube root of conservative volume factor.
   30 NNB2 = 0.8*NNBMAX
      A1 = FLOAT(NNB2)/FLOAT(NNB)
      IF (RS(I).GT.5.0*RSCALE) THEN
          A1 = MIN(5.0*A1,0.9D0)
          IRSKIP = 1
      END IF
      RS2 = RS2*A1**0.66667
      RCRIT2 = 1.59*RS2
      RS(I) = SQRT(RS2)
      NNB1 = 0
*
      DO 35 L = 1,NNB
          J = ILIST(L+1)
          IF (L + NNB2.GT.NNB + NNB1) GO TO 32
*       Sum of neighbours (NNB1) & those left (NNB+1-L) set to NNB2.
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
*
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
          DRDP = 3.0*DRDV*DR2I

          IF (RIJ2.GT.RCRIT2) GO TO 34

          IF (RIJ2.GT.RS2) THEN
             IF (DRDV.GT.VRFAC.AND.J.LE.N) GO TO 34
*       Retain any c.m. because of complications in force correction.
          END IF
*
   32     NNB1 = NNB1 + 1
          JLIST(NNB1+1) = J
          GO TO 35
*
*       Subtract neighbour force included above and add to regular force.
   34     FIRR(1) = FIRR(1) - A1*DR3I
          FIRR(2) = FIRR(2) - A2*DR3I
          FIRR(3) = FIRR(3) - A3*DR3I
          FD(1) = FD(1) - (DV(1) - A1*DRDP)*DR3I
          FD(2) = FD(2) - (DV(2) - A2*DRDP)*DR3I
          FD(3) = FD(3) - (DV(3) - A3*DRDP)*DR3I
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDP)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDP)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDP)*DR3I
*       Obtain potential and derivative.
          PHII(I) = PHII(I) - DR3I*RIJ2
          PHIR(I) = PHIR(I) - DR3I*RIJ2
          PHIR1(I) = PHIR1(I) + DRDV*DR3I
   35 CONTINUE
*
      DO 38 L = 2,NNB1+1
          ILIST(L) = JLIST(L)
   38 CONTINUE
      NNB = NNB1
      NBFULL = NBFULL + 1
*       See whether to reduce NNB further.
      IF (NNB.GE.NNBMAX) GO TO 30
*
*       Stabilize NNB between ZNBMIN & ZNBMAX by square root of contrast.
*       Include optional stabilization to increase neighbour number.
*       Take input parameter NNBOPT as optimal neighbour number (R.Sp.)
*       Note that it substitutes input parameter NNBMAX, which
*       is now a parameter NNBMAX=LMAX-3
   40 CONTINUE
*
      FAC = 1.D0
      IF (KZ(40).GT.0) THEN
          FAC = 1.0 + 0.1*(FLOAT(NNBOPT) - FLOAT(NNB))/FLOAT(NNB)
      END IF
*
*     A3 = ALPHA*FAC*SQRT(FLOAT(NNB)*RS(I))/RS2
      A3 = FLOAT(NNB)*FAC
      A3 = MIN(A3,0.9*ZNBMAX)
*       Reduce predicted membership slowly outside half-mass radius.
      IF (RI2.GT.RH2) THEN
          A3 = A3*RSCALE/SQRT(RI2)
      ELSE
          A3 = MAX(A3,ZNBMIN)
      END IF
      A4 = A3/FLOAT(NNB)
*
*       Include inertial factor to prevent resonance oscillations in RS.
      IF ((A3 - FLOAT(NNB0))*(A3 - FLOAT(NNB)).LT.0.0) A4 = SQRT(A4)
*
*       Modify volume ratio by radial velocity factor outside the core.
      IF (RI2.GT.RC2.AND.KZ(39).EQ.0.AND.RI2.LT.100.0*RH2) THEN
          RIDOT = (XI(1) - RDENS(1))*XIDOT(1) +
     &            (XI(2) - RDENS(2))*XIDOT(2) +
     &            (XI(3) - RDENS(3))*XIDOT(3)
          A4 = A4*(1.0 + RIDOT*DTR/RI2)
      END IF
*
*       See whether neighbour radius of c.m. body should be increased.
      IF (I.GT.N) THEN
*       Set perturber range (soft binaries & H > 0 have short duration).
          A2 = 100.0*BODY(I)/ABS(H(I-N))
*       Stabilize NNB on ZNBMAX if too few perturbers.
          IF (A2.GT.RS(I)) THEN
              A3 = MAX(1.0 - FLOAT(NNB)/ZNBMAX,0.0D0)
*       Modify volume ratio by approximate square root factor.
              A4 = 1.0 + 0.5*A3
          END IF
      END IF
*
*       Limit change of RS for small steps (RSFAC = MAX(25/TCR,1.5*VC/RSMIN).
      A5 = MIN(RSFAC*DTR,0.25D0)
*       Restrict volume ratio by inertial factor of 25 per cent either way.
      A4 = A4 - 1.0
      IF (A4.GT.A5) THEN
          A4 = A5
      ELSE
          A4 = MAX(A4,-A5)
      END IF
*
*       Modify neighbour sphere radius by volume factor.
      IF (IRSKIP.EQ.0) THEN
          A3 = ONE3*A4
          A1 = 1.0 + A3 - A3*A3
*       Second-order cube root expansion (maximum volume error < 0.3 %).
          IF (RS(I).GT.5.0*RSCALE) A1 = SQRT(A1)
*       Skip modification for small changes (avoids oscillations in RS).
          IF (ABS(A1 - 1.0D0).GT.0.003) THEN
              RS(I) = A1*RS(I)
          END IF
      END IF
*
*       Calculate the radial velocity with respect to at most 3 neighbours.
      IF (NNB.LE.3.AND.RI2.LT.100.0*RH2) THEN
          A1 = 2.0*RS(I)
*
          DO 45 L = 1,NNB
              J = ILIST(L+1)
              RIJ = SQRT((XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 +
     &                                         (XI(3) - X(3,J))**2)
              RSDOT = ((XI(1) - X(1,J))*(XIDOT(1) - XDOT(1,J)) +
     &                 (XI(2) - X(2,J))*(XIDOT(2) - XDOT(2,J)) +
     &                 (XI(3) - X(3,J))*(XIDOT(3) - XDOT(3,J)))/RIJ
*       Find smallest neighbour distance assuming constant regular step.
              A1 = MIN(A1,RIJ + RSDOT*DTR)
   45     CONTINUE
*
*       Increase neighbour sphere if all members are leaving inner region.
          RS(I) = MAX(A1,1.1*RS(I))
      END IF
*
*       Check optional procedures for adding neighbours.
      IF ((KZ(37).EQ.1.AND.LISTV(1).GT.0).OR.KZ(37).GT.1) THEN
          CALL CHECKL(I,NNB,XI,XIDOT,RS2,FIRR,FREG,FD,FDR)
      END IF
*
*       Find loss or gain of neighbours at the same time.
   50 NBLOSS = 0
      NBGAIN = 0
*
*       Check case of zero old or new membership (skip if both are zero).
      IF (NNB0.EQ.0) THEN
          NBGAIN = NNB
          GO TO 70
      END IF
*
      JMIN = 0
      L = 2
      LG = 2
*       Set termination value in ILIST(NNB+2) and save last list member.
      ILIST(NNB+2) = NTOT + 1
      ILIST(1) = KLIST(NNB0+1)
*
*       Compare old and new list members in locations L & LG.
   56 IF (KLIST(L).EQ.ILIST(LG)) GO TO 58
*
*       Now check whether inequality means gain or loss.
      IF (KLIST(L).GE.ILIST(LG)) THEN
          NBGAIN = NBGAIN + 1
          JLIST(NNB0+NBGAIN) = ILIST(LG)
*       Number of neighbour losses can at most be NNB0.
          L = L - 1
*       The same location will be searched again after increasing L below.
      ELSE
          NBLOSS = NBLOSS + 1
          J = KLIST(L)
          JLIST(NBLOSS) = J
*       Check SMIN step indicator (rare case permits fast skip below).
          IF (STEP(J).LT.SMIN) JMIN = J
          LG = LG - 1
      END IF
*
*       See whether the last location has been checked.
   58 IF (L.LE.NNB0) THEN
          L = L + 1
          LG = LG + 1
*       Last value of second search index is NNB + 2 which holds NTOT + 1.
          GO TO 56
      ELSE IF (LG.LE.NNB) THEN
          LG = LG + 1
          KLIST(L) = NTOT + 1
*       Last location of list holds termination value (saved in ILIST(1)).
          GO TO 56
      END IF
*
*       See whether any old neighbour with small step should be retained.
      IF (JMIN.EQ.0) GO TO 70
*
      K = 1
   60 IF (NNB.GT.NNBMAX.OR.I.GT.N) GO TO 70
      J = JLIST(K)
*       A single regularized component will be replaced by the c.m.
      IF (STEP(J).GT.SMIN.OR.J.LT.IFIRST.OR.J.GT.N) GO TO 68
*       Retain old neighbour inside 2*RS to avoid large correction terms.
      RIJ2 = (XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 +
     &                             (XI(3) - X(3,J))**2
      IF (RIJ2.GT.4.0*RS2) GO TO 68
*
      L = NNB + 1
   62 IF (ILIST(L).LT.J) GO TO 64
      ILIST(L+1) = ILIST(L)
      L = L - 1
      IF (L.GT.1) GO TO 62
*
*       Save index of body #J and update NNB & NBLOSS.
   64 ILIST(L+1) = J
      NNB = NNB + 1
      NBLOSS = NBLOSS - 1
*       Restore last old neighbour in case NBLOSS = 0 at end of search.
      KLIST(NNB0+1) = ILIST(1)
      NBSMIN = NBSMIN + 1
*
*       Perform correction to irregular and regular force components.
      A1 = X(1,J) - XI(1)
      A2 = X(2,J) - XI(2)
      A3 = X(3,J) - XI(3)
      DV(1) = XDOT(1,J) - XIDOT(1)
      DV(2) = XDOT(2,J) - XIDOT(2)
      DV(3) = XDOT(3,J) - XIDOT(3)
*
      RIJ2 = A1*A1 + A2*A2 + A3*A3
      DR2I = 1.0/RIJ2
      DR3I = BODY(J)*DR2I*SQRT(DR2I)
      DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
      DRDP = 3.0*DRDV*DR2I
*
      FIRR(1) = FIRR(1) + A1*DR3I
      FIRR(2) = FIRR(2) + A2*DR3I
      FIRR(3) = FIRR(3) + A3*DR3I
      FD(1) = FD(1) + (DV(1) - A1*DRDP)*DR3I
      FD(2) = FD(2) + (DV(2) - A2*DRDP)*DR3I
      FD(3) = FD(3) + (DV(3) - A3*DRDP)*DR3I
      FREG(1) = FREG(1) - A1*DR3I
      FREG(2) = FREG(2) - A2*DR3I
      FREG(3) = FREG(3) - A3*DR3I
      FDR(1) = FDR(1) - (DV(1) - A1*DRDP)*DR3I
      FDR(2) = FDR(2) - (DV(2) - A2*DRDP)*DR3I
      FDR(3) = FDR(3) - (DV(3) - A3*DRDP)*DR3I
*       Obtain potential and derivative.
      PHII(I) = PHII(I) - DR3I*RIJ2
      PHIR(I) = PHIR(I) - DR3I*RIJ2
      PHIR1(I) = PHIR1(I) + DRDV*DR3I
*
*       Remove body #J from JLIST unless it is the last or only member.
      IF (K.GT.NBLOSS) GO TO 70
      DO 66 L = K,NBLOSS
          JLIST(L) = JLIST(L+1)
   66 CONTINUE
*       Index of last body to be moved up is L = NBLOSS + 1.
      K = K - 1
*       Check the same location again since a new body has appeared.
   68 K = K + 1
*       Last member to be considered is in location K = NBLOSS.
      IF (K.LE.NBLOSS) GO TO 60
*
*       Form time-step factors and update regular force time.
   70 CONTINUE
*
      DTSQ = DTR**2
      DT6 = 6.0/(DTR*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DTR13 = ONE3*DTR
      T0R(I) = TIME
*
*       Suppress the corrector for large time-step ratios (experimental).
      IF (STEP(I).LT.5.0D0*DTMIN.AND.DTR.GT.50.0*STEP(I)) THEN
          DTR13 = 0.0
          DTSQ12 = 0.0
      END IF
*
*     LCHECK = NBLOSS+NBGAIN.EQ.0.AND.
*    *   ((STEP(I).LT.1.D-09.AND.DFIRR.GT.1.0D-9).OR.
*    *    (STEPR(I).LT.DTMIN/1.D1.AND.DFREG.GT.1.0D-9))
*     IF (LCHECK) THEN
*         DO K=1,3
*             DFI(K) = (FIRR(K)-FI(K,I))**2 
*             DFRG(K) = (FREG(K)-FR(K,I))**2
*         END DO
*         DFIRR = DSQRT(DFI(1)+DFI(2)+DFI(3))**2
*         DFREG = DSQRT(DFRG(1)+DFRG(2)+DFRG(3))**2
*     IF (rank.eq.0)
*    *WRITE (21,740)  TIME, I, NNB, DFIRR, DFREG, STEP(I), STEPR(I)
* 740 FORMAT ('-----------------------------------',/,
*    *     ' CHECK! T I NB DFI DFR DTi,r',1P,E12.4,2I6,4E12.4)
*         WRITE(21,*)NBLOSS,NBGAIN
*         DO K=1,NNB
*         J = ILIST(K)
*         IF (J.GT.N) THEN
*             JPAIR = J-N
*             J2 = 2*JPAIR-1
*         if(rank.eq.0)
*    *    WRITE (21,742) J,J2,LIST(1,J2),GAMMA(JPAIR)
* 742 FORMAT (' CHECK  J J2 L Pert         ',3I6,1P,E12.4)
*         END IF
*         END DO
*         IF (I.GT.N)THEN
*             IPAIR = I-N
*             I2 = 2*IPAIR-1
*         if(rank.eq.0)
*    *    WRITE (21,741) I,I2,LIST(1,I2),GAMMA(IPAIR)
* 741 FORMAT (' CHECK  I I2 L Pert          ',3I6,1P,E12.4)
*     CALL FLUSH(21)
*         END IF
*     END IF
*       Test for inconsistencies of nbint, regint and cmfreg Nov 2010
*       Include the corrector and save F, FI, FR & polynomial derivatives.
      DO 75 K = 1,3
*       Subtract change of neighbour force to form actual first derivative.
          DFR = FR(K,I) - (FIRR(K) - FI(K,I)) - FREG(K)
          FDR0 = FDR(K) - (FIDOT(K,I) - FD(K))
*
          FRD = FRDOT(K,I)
          SUM = FRD + FDR0
          AT3 = 2.0*DFR + DTR*SUM
          BT2 = -3.0*DFR - DTR*(SUM + FRD)
*       Use here new variables for consistency in parallel execution (R.Sp.)
          XN(K,I) = X0(K,I) + (0.6*AT3 + BT2)*DTSQ12
          XNDOT(K,I) = X0DOT(K,I) + (0.75*AT3 + BT2)*DTR13
*
          FI(K,I) = FIRR(K)
          FR(K,I) = FREG(K)
          FIDOT(K,I) = FD(K)
          FRDOT(K,I) = FDR(K)
*
          D0(K,I) = FIRR(K)
          D0R(K,I) = FREG(K)
          D1R(K,I) = FDR(K)
          D3R(K,I) = AT3*DT6
          D2R(K,I) = (3.0*AT3 + BT2)*DT2
   75 CONTINUE
*
*       Correct force polynomials due to neighbour changes (KZ(38) or I > N).
      IF (KZ(38).GT.0) THEN
          CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT,FIRR,FREG,FD,FDR,KLIST)
      ELSE IF (I.GT.N) THEN
          IF (GAMMA(I-N).GT.GMAX) THEN
          CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT,FIRR,FREG,FD,FDR,KLIST)
          END IF
      END IF
*
*       Copy new neighbour list if different from old list.
      IF (NBLOSS + NBGAIN.GT.0) THEN
          KLIST(1) = NNB
          DO 80 L = 2,NNB+1
              KLIST(L) = ILIST(L)
   80     CONTINUE
          NBFLUX = NBFLUX + NBLOSS + NBGAIN
      END IF
*
*       Obtain new regular integration step using composite expression.
*       STEPR = (ETAR*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
      FR2 = FREG(1)**2 + FREG(2)**2 + FREG(3)**2
      W1 = FDR(1)**2 + FDR(2)**2 + FDR(3)**2
      W2 = D2R(1,I)**2 + D2R(2,I)**2 + D2R(3,I)**2
      W3 = D3R(1,I)**2 + D3R(2,I)**2 + D3R(3,I)**2
*
*       Form new step by relative criterion.
      W0 = (SQRT(FR2*W2) + W1)/(SQRT(W1*W3) + W2)
      W0 = ETAR*W0
      TTMP = SQRT(W0)
      DT0 = TTMP
*
*       Winston Sweatman's suggestion
*     DVV = (XDOT(1,I)-X0DOT(1,I))**2 + (XDOT(2,I)-X0DOT(2,I))**2 +
*    &     (XDOT(3,I)-X0DOT(3,I))**2
*     FFD = FREG(1)**2 + FREG(2)**2 + FREG(3)**2
*     ETARW = ETAR
*     TTMPW = ETARW*DVV*BODY(I)/FFD
*
*     PRINT*,' Reg I=',I,' TTMP,TTMPW,RATIO=',
*    &  TTMP,TTMPW,TTMP/TTMPW
*
*     IF(TTMP.GT.TTMPW)THEN
*     IGT = IGT + 1
*     ELSE
*     ILE = ILE + 1
*     END IF
*     IF(MOD(IGT+ILE,100).EQ.0)PRINT*,' irr IGT,ILE=',IGT,ILE
*
*     TTMP = MAX(TTMPW,TTMP)
*     TTR = TSTEP(FREG,FDR,D2R(1,I),D3R(1,I),ETAR)
*
*       Adopt FAC*MIN(FREG,FIRR) (or tidal force) for convergence test.
      FAC = MIN(ETAR,0.04D0)
      IF (KZ(14).EQ.0) THEN
          FI2 = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
          DF2 = FAC**2*MIN(FR2,FI2)
*       Ignore irregular force criterion if no neighbours.
          IF (NNB.EQ.0) DF2 = FAC**2*FR2
      ELSE IF (KZ(14).EQ.1) THEN
          W0 = (TIDAL(1)*XI(1))**2
          DF2 = FAC**2*MAX(FR2,W0)
      ELSE
          DF2 = FAC**2*FR2
      END IF
*
*       Obtain regular force change using twice the predicted step.
      DTC = 2.0*TTMP
      S2 = 0.5*DTC
      S3 = ONE3*DTC
      W0 = 0.0
      DO 105 K = 1,3
          W2 = ((D3R(K,I)*S3 + D2R(K,I))*S2 + FDR(K))*DTC
          W0 = W0 + W2**2
  105 CONTINUE
*
*       See whether regular step can be increased by factor 2.
      IF (W0.LT.DF2) THEN
          TTMP = DTC
      END IF
*
*       Impose a smooth step reduction inside compact core.
      IF (NC.LT.50.AND.RI2.LT.RC2) THEN
          TTMP = TTMP*MIN(1.0D0,0.5D0*(1.0D0 + RI2*RC2IN))
      END IF
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
		IF (TTMP .GT. 2.0*STEPR(I)) THEN
			IF (DMOD(TIME,2.0*STEPR(I)) .EQ. 0.0D0) THEN
              TTMP = MIN(2.0*STEPR(I),1.D0)
          ELSE
              TTMP = STEPR(I)
          END IF
      ELSE IF (TTMP .LT. STEPR(I)) THEN
          TTMP = 0.5*STEPR(I)
      ELSE
          TTMP = STEPR(I)
      END IF
*
*       Set new regular step and reduce irregular step if STEPR < STEP.
      STEPR(I) = TTMP

*     --03/30/13 22:34-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      if (rank.eq.0.and.nstepr.gt.4314800) then
*         write (6,*) name(i), nnb,time, freg(1)
* 118     format ('name, nnb, freg, time', 2I6, 1F2.6, 1E10.2)
*      end if
*     --03/30/13 22:34-lwang-end----------------------------------------*

  110 IF (TTMP.LT.STEP(I)) THEN
          STEP(I) = 0.5D0*STEP(I)
          TIMENW(I) = T0(I) + STEP(I)
          NICONV = NICONV + 1
          GO TO 110
      END IF
*
      RETURN
*
      END
