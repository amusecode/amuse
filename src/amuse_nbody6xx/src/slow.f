      SUBROUTINE SLOW
*
*
*       Slow-down treatment of chain binary.
*       ------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      LOGICAL  KSLOW,KCOLL,KCASE
      REAL*8  KSCH,KSNEW,VI(NMX3),VC(NMX3),RC1(3),RC2(3)
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW,KCOLL
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/KSAVE/   K1,K2
      COMMON/SLOW3/  GCRIT,KZ26
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      SAVE

*
*
*       Check for switching off slow-down at start of iteration.
      if (GCRIT.eq.0.0d0) then
          ksnew = 1.0
          go to 90
      end if
*
*       Set logical variable to avoid multiple copies of QK & PK.
      KCASE = .FALSE.
*
*       Perform perturbation check if slow-down not active.
      IF (.NOT.KSLOW) THEN
*       Determine chain index and largest inverse distance.
          RM = 0.0
          DO 1 I = 1,N-1
              IF (RINV(I).GT.RM) THEN
                  RM = RINV(I)
                  I1 = I
              END IF
    1     CONTINUE
*
*       Check carefully two possible binaries (eccentricity effect).
          IF (N.GT.3) THEN
              KCASE = .TRUE.
*       Save QK & PK and copy current configuration for EREL & TRANSK.
              DO 5 I = 1,N-1
                  KS = 4*(I - 1)
                  DO 4 J = 1,4
                      QK(KS+J) = Q(KS+J)
                      PK(KS+J) = P(KS+J)
    4             CONTINUE
    5         CONTINUE
*       Evaluate first semi-major axis from non-singular variables.
              K1 = INAME(i1)
              K2 = INAME(i1+1)
              CALL EREL(i1,EB,SEMI)
*       Determine index of second smallest separation.
              RI2 = 0.0
              DO 10 I = 1,N-1
                  IF (I.NE.I1.AND.RINV(I).GT.RI2) THEN
                      I2 = I
                      RI2 = RINV(I)
                  END IF
   10         CONTINUE
*       Obtain second semi-major axis (errors only affect decision-making).
              K1 = INAME(i2)
              K2 = INAME(i2+1)
              CALL EREL(i2,EB,SEMI2)
*       Switch to #I2 if second binary is smaller or first pair is not bound.
              IF (SEMI2.GT.0.0) THEN
                  IF (SEMI2.LT.SEMI.OR.SEMI.LT.0.0) THEN
                      I1 = I2
                      SEMI = SEMI2
                      RM = RINV(I2)
                  END IF
*       Exit in case of two hyperbolic two-body motions.
              ELSE IF (SEMI.LT.0.0) THEN
                  GO TO 100
              END IF
          END IF
*
*       Sum the perturbing forces m/r^3 next to #i1 (two terms if i1 = 2).
          sum = 0.0
          do 15 i = 1,n-1
              if (iabs(i-i1).eq.1) then
                  j = i
                  if (i.gt.i1) j = i + 1
                  sum = sum + mc(j)*rinv(i)**3
              end if
   15     continue
*
*       Include one more contribution for two consecutive perturbers.
          if (i1.eq.1) then
              LJ = 3*i1
              do k = 1,3
                  RC2(k) = XC(k+LJ+3) + XC(k+LJ)
              end do
*       Add perturbation from second subsequent member (i = i1 + 2).
              RJ = SQRT(RC2(1)**2 + RC2(2)**2 + RC2(3)**2)
              j = i1 + 3
              sum = sum + mc(j)/RJ**3
          else if (i1.ge.3) then
              LJ = 3*(i1 - 2)
              do k = 1,3
                  RC2(k) = XC(k+LJ-3) + XC(k+LJ)
              end do
*       Add the previous perturbation (neglected in do 15 loop).
              RJ = SQRT(RC2(1)**2 + RC2(2)**2 + RC2(3)**2)
              j = i1 - 2
              sum = sum + mc(j)/RJ**3
          end if
*
*       Skip further search if relative perturbation exceeds limit.
          mb = mc(i1) + mc(i1+1)
          GAMMA = 2.0*sum/(mb*RM**3)
          IF (GAMMA.GT.GCRIT) THEN
              GO TO 100
          END IF
      END IF
*
*       Determine chain index for slow-down binary (otherwise from above).
      IF (KSLOW.AND.N.GT.3) THEN
          DO 20 I = 1,N-1
              IF (KSCH(I).GT.1.0D0) i1 = i
   20     CONTINUE
*       See whether a closer particle pair is present (factor of 2).
          I2 = 0
          R1 = 1.0/RINV(I1)
          DO 60 I = 1,N-1
              RI2 = 1.0/RINV(I)
              IF (I.NE.I1.AND.R1.LT.0.5*RI2) THEN
                  I2 = I
                  IF (R1.LT.0.1*SEMI.AND.RI2.LT.0.1*SEMI) THEN
                      GO TO 100
                  END IF
*       Compare closest separation with current slow-down binary.
                  IF (RI2.LT.0.5*SEMI.AND.R1.GT.0.5*SEMI) THEN
*       Evaluate semi-major axis directly (cf. small RI2 in EREL).
                      L = 3*(N-2)
                      DO 25 K = 1,3
                          VI(K) = -WC(K)/MC(1)
                          VI(L+K+3) = WC(L+K)/MC(N)
   25                 CONTINUE
                      DO 35 II = 2,N-1
                          L = 3*(II-1)
                          DO 30 K = 1,3
                              VI(L+K) = (WC(L+K-3) - WC(L+K))/MC(II)
   30                     CONTINUE
   35                 CONTINUE
                      DO 40 J = 1,3*(N-1)
                          VC(J) = VI(J+3) - VI(J)
   40                 CONTINUE
                      L = 3*(I1-1)
                      R2 = XC(L+1)**2 + XC(L+2)**2 + XC(L+3)**2
                      W2 = VC(L+1)**2 + VC(L+2)**2 + VC(L+3)**2
                      SEMI = 2.0D0/R1 - W2/(MC(I1) + MC(I1+1))
                      SEMI = 1.0D0/SEMI
                      GO TO 80
                  ELSE IF (RI2.LT.MIN(SEMI,R1).AND.R1.GT.0.5*SEMI) THEN
*       Determine second binary by regular expression (R1 not too small).
                      DO 50 II = 1,N-1
                          KS = 4*(II - 1)
                          DO 45 J = 1,4
                              QK(KS+J) = Q(KS+J)
                              PK(KS+J) = P(KS+J)
   45                     CONTINUE
   50                 CONTINUE
                      K1 = INAME(i2)
                      K2 = INAME(i2+1)
                      CALL EREL(i2,EB2,SEMI2)
                      IF (SEMI2.GT.0.0.AND.SEMI2.LT.SEMI) THEN
                          KSNEW = 1.0D0
                          GO TO 90
                      ELSE
*       Continue with the current binary (ie. small change in perturbation).
                          GO TO 80
                      END IF
                  END IF
              END IF
   60     CONTINUE
      END IF
*
*       Obtain regular semi-major axis for missing cases (including N = 3).
      IF (.NOT.KCASE) THEN
*       Save QK & PK and copy current configuration for EREL & TRANSK.
          DO 70 I = 1,N-1
              KS = 4*(I - 1)
              DO 65 J = 1,4
                  QK(KS+J) = Q(KS+J)
                  PK(KS+J) = P(KS+J)
   65         CONTINUE
   70     CONTINUE
*
*       Evaluate semi-major axis from non-singular variables.
          K1 = INAME(i1)
          K2 = INAME(i1+1)
          CALL EREL(i1,EB,SEMI)
*       Exit if no current binary (set KSLOW = .false. just in case).
          IF (SEMI.LE.0.0d0) THEN
              KSLOW = .false.
              TK2(0) = 0.0
              GO TO 100
          END IF
      END IF
*
*       Check for switching to smaller binary (exchange leads to escape).
      IF (KSLOW.AND.N.GT.3) THEN
          IF (I2.GT.0) THEN
*       Evaluate second semi-major axis (K1 & K2 can be over-written).
              K1 = INAME(i2)
              K2 = INAME(i2+1)
              CALL EREL(i2,EB2,SEMI2)
              IF (SEMI2.GT.0.0.AND.SEMI2.LT.SEMI) THEN
*       Switch off the present pair to prepare re-activation.
                  ksnew = 1.0
                  go to 90
              END IF
          END IF
      END IF
*
*       Sum the perturbations (on either side and non-dominant terms).
   80 IF (KSLOW) THEN
          sum = 0.0
          do i = 1,n-1
*
*       Include full perturbation (non-symmetrical i1 for n = 4 or n > 4).
              if ((n.eq.4.and.i1.ne.2).or.n.gt.4) then
                  do k = 1,3
                      RC1(K) = 0.0
                  end do
*       Save vector sum on either side of #i1 excluding closest neighbour.
                  if (i.lt.i1-1) then
                      do j = i+1,i1-1
                          LJ = 3*(j-1)
                          do k = 1,3
                              RC1(k) = RC1(k) + XC(k+LJ)
                          end do
                      end do
*       Check alternative case of subsequent distant members (i > i1 + 1).
                  else if (i.gt.i1+1) then
                      do j = i1+1,i-1
                          LJ = 3*(j-1)
                          do k = 1,3
                              RC1(k) = RC1(k) + XC(k+LJ)
                          end do
                      end do
                  end if
              end if
*
*       Set chain offset and mass reference index.
              L = 3*(I-1)
              j = i
              if (i.gt.i1) j = i + 1
*       Use actual separation if perturber is next to dominant binary.
              if (iabs(i-i1).eq.1) then
                  sum = sum + mc(j)*rinv(i)**3
              else if (i.ne.i1) then
*       Include chain vector to yield full distance to binary.
                  do k = 1,3
                      RC2(k) = RC1(k) + XC(k+L)
                  end do
*       Add contribution from more distant member.
                  RJ = SQRT(RC2(1)**2 + RC2(2)**2 + RC2(3)**2)
                  sum = sum + mc(j)/RJ**3
              end if
          end do
      END IF
*
*       Form relative perturbation at maximum apocentre.
      rap = 2.0*SEMI
      pert = 2.0*sum*rap**3/(mc(i1) + mc(i1+1))
*
*       Specify the slow-down factor.
      IF (pert.LT.GCRIT) THEN
          ksnew = SQRT(GCRIT/pert)
      ELSE
          ksnew = 1.0
      END IF
*
*     --------------------------------------------------------------
*       Implement discrete scheme (suppressed).
*     i = i1
*     rat = ksch(i)
*     rat = ksnew/rat
*       Check for significant changes (by 2) in slow-down factor.
*     if (rat.gt.0.5.and.rat.le.2.0) then
*         ksnew = Ksch(i)
*     else if (rat.gt.2.0) then
*         ksnew = 2*Ksch(i)
*       Allow an extra factor of 32 at initialization.
*         if (.not.KSLOW) then
*             if (rat.ge.64.0) ksnew = 32*ksnew
*         end if
*     else if (rat.le.0.5) then
*         ksnew = Ksch(i)/2
*         if (rat.le.0.125) ksnew = Ksch(i)/4
*         ksnew = max(ksnew,1.0D0)
*     end if
*     --------------------------------------------------------------
*
*       Check slow-down switch and include any change in binding energy.
   90 if (ksnew.ne.ksch(i1)) then
          if (ksnew.eq.1.0d0) then
              KSLOW = .false.
          else
              KSLOW = .true.
          end if
*       Add change in binary energy and save new slow-down index.
          eb = -0.5d0*MC(i1)*MC(i1+1)/SEMI
          DEB = eb*(1.0/ksnew - 1.0/Ksch(i1))
          Ksch(i1) = ksnew
      else
          DEB = 0.0
      end if
*
*       Accumulate total energy difference due to slow-down.
      EJUMP = EJUMP + DEB
*       Ensure zero EJUMP without slow-down to avoid small residual.
      IF (.NOT.KSLOW) EJUMP = 0.0
*
*       Form modified mass factors.
      do i = 1,n
          tk1(i) = -1.0/MC(I)
      end do
      do i = 1,n-1
          if (Ksch(i).ne.1.0d0) then
              tk1(i) = tk1(i)/Ksch(i)
              tk1(i+1) = tk1(i+1)/Ksch(i)
          end if
      end do
      DO I = 1,N-1
          TKK(I) = 0.5D0*(-tk1(i) - tk1(i+1))
          MKK(I) = MC(I)*MC(I+1)/Ksch(i)
      END DO
      do i = 1,n-1
          m12 = mc(i) + mc(i+1)
          dt12 = 0.5d0*(1.0d0 - 1.0d0/Ksch(i))/m12
          if (i.gt.1) TKK(i-1) = tkk(i-1) + dt12
          if (i.lt.n-1) TKK(i+1) = tkk(i+1) + dt12
          if (i.gt.1.and.i.lt.n-1) TK2(i) = -2.0d0*dt12
      end do
      TK2(0) = 0.0
      TK2(N) = 0.0
*
  100 RETURN
*
      END
