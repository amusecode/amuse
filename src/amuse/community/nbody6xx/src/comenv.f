***
      SUBROUTINE COMENV(M01,M1,MC1,AJ1,JSPIN1,KW1,
     &                  M02,M2,MC2,AJ2,JSPIN2,KW2,ECC,SEP,COEL)
*
* Common Envelope Evolution.
*
      INCLUDE 'common6.h'
*
      INTEGER KW1,KW2,KW
*     INTEGER IDUM
*
      REAL*8 M01,M1,MC1,AJ1,JSPIN1,R1,L1
      REAL*8 M02,M2,MC2,AJ2,JSPIN2,R2,L2,MC22
      REAL*8 TSCLS1(20),TSCLS2(20),LUMS(10),GB(10),TM1,TM2,TN
      REAL*8 EBINDI,EBINDF,EORBI,EORBF,ECIRC
      REAL*8 AA,BB,CC,EFAC
      REAL*8 ECC,SEP,SEPF,SEPL,MF,XX
      REAL*8 CONST,DELY,DERI,DELMF,MC3,FAGE1,FAGE2
      REAL*8 TB,OORB,OSPIN1,OSPIN2
      REAL*8 RC1,RC2,Q1,Q2,RL1,RL2
      REAL*8 MENV,RENV,MENVD,RZAMS
      REAL*8 LAMB1,LAMB2,K21,K22
      REAL*8 AURSUN,K3,LAMBDA,ALPHA1,MCH
      PARAMETER (AURSUN = 214.95D0,K3 = 0.21D0)
      PARAMETER (LAMBDA = 0.D0,ALPHA1 = 3.D0)
      PARAMETER (MCH = 1.44D0)
      LOGICAL COEL,ECCFLG
      REAL*8 CELAMF,RL,RZAMSF
      EXTERNAL CELAMF,RL,RZAMSF
*
* Common envelope evolution - entered only when KW1 = 2, 3, 4, 5, 6, 8 or 9.
*
* For simplicity energies are divided by -G.
*
      COEL = .FALSE.
      ECCFLG = .TRUE.
*
* Obtain the core masses and radii.
*
      KW = KW1
      CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
      IF(AJ1.GT.TN) AJ1 = TN 
      CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &            R1,L1,KW1,MC1,RC1,MENV,RENV,K21)
      OSPIN1 = JSPIN1/(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
      IF(M1.GT.MC1)THEN
         MENVD = MENV/(M1-MC1)
         RZAMS = RZAMSF(M01)
         LAMB1 = CELAMF(KW,M01,L1,R1,RZAMS,MENVD,LAMBDA)
      ELSE
         LAMB1 = 0.5D0
      ENDIF
      IF(KW.NE.KW1.and.rank.eq.0) THEN
          WRITE(38,*)' COMENV TYPE CHANGE *1'
      END IF
*
      KW = KW2
      CALL star(KW2,M02,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
      IF(AJ2.GT.TN) AJ2 = TN 
      CALL hrdiag(M02,AJ2,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS,
     &            R2,L2,KW2,MC2,RC2,MENV,RENV,K22)
      OSPIN2 = JSPIN2/(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
      IF(KW.NE.KW2.and.rank.eq.0) THEN
          WRITE(38,*)' COMENV TYPE CHANGE *2'
      END IF
*
      IF(SEP.GT.0.D0)THEN
         TB = (SEP/AURSUN)*SQRT(SEP/(AURSUN*(M1+M2)))
         OORB = TWOPI/TB
      ENDIF
*
      if(rank.eq.0) WRITE(38,66) KW1,M01,M1,MC1,AJ1,SEP
 66   FORMAT(' COMENV OLD *1:',I4,3F7.2,F9.2,F12.4)
      if(rank.eq.0) WRITE(38,67) KW2,M02,M2,MC2,AJ2
 67   FORMAT(' COMENV OLD *2:',I4,3F7.2,F9.2)
*
* Calculate the binding energy of the giant envelope (multiplied by lambda).
*
      EBINDI = M1*(M1-MC1)/(LAMB1*R1)
*
* Calculate the initial orbital energy
*
      EORBI = MC1*M2/(2.D0*SEP)
*
* Allow for an eccentric orbit.
*
*     ECIRC = EORBI/(1.D0 - ECC*ECC)
*
      IF(ECCFLG.AND.ECC.GT.0.01D0)THEN
*
* Assume 10% of the binding energy is lost on the first orbit and determine 
* whether the star just passes through. For the moment put all the 
* secondary's envelope with the primary's.
*
*        DE = MIN(EBINDI/10.D0,(ECIRC-EORBI)*ALPHA1)
*
* Adopt new procedure based on a cubic fitting (May 2003).
* (Oct 2005: YP is maximum for a head-on collision and would reduce 
*            to zero for a grazing encounter, i.e. QP=1. However, 
*            in order to be consistent with Kochanek collision detection 
*            which allows QP values up to 1.7 on entry to expel/comenv 
*            we set a non-zero YP for QP > 1. This is a temporary fix 
*            until a better collision detection scheme comes along sometime 
*            in the near future, we hope.)
*
         QP = SEP*(1.D0 - ECC)/(R1 + R2)
         YP = 0.5D0*(QP - 0.5D0)**3 - 0.375D0*(QP - 0.5D0) + 0.125D0
         YP = MAX(YP,0.05D0)
         YP = MIN(YP,0.5D0)
         DE = YP*EBINDI
* Remove all envelope mass if below 0.05 m_sun.
         IF (M1-MC1.LT.0.05) DE = EBINDI
         EBINDF = EBINDI - DE
         MF = 0.5D0*(MC1 + SQRT(MC1*MC1 + 4.D0*LAMB1*R1*EBINDF))
         EORBI = MF*M2/(2.D0*SEP)
         ECIRC = EORBI/(1.D0 - ECC*ECC)
         EORBF = EORBI + DE/ALPHA1
         IF(EORBF.LT.ECIRC)THEN
            M1 = MF
            ECC = SQRT(1.D0 - EORBF/ECIRC)
         ELSE
*           EFAC = 1.d0 - 1.d0/(ECC*ECC)
            EFAC = ECC*ECC/(1.d0 - ECC*ECC)
            AA = 1.d0/(LAMB1*R1)
*           BB = (2.d0*M1-MC1)/(LAMB1*R1) + ALPHA1*M2/(SEP*EFAC)
            BB = (2.d0*M1-MC1)/(LAMB1*R1) + ALPHA1*EFAC*M2/(2.d0*SEP)
*           CC = M1*MC1/(LAMB1*R1) + ALPHA1*M1*M2/(SEP*EFAC)
            CC = ALPHA1*EFAC*M1*M2/(2.d0*SEP)
            DM = (BB - SQRT(BB*BB - 4.d0*AA*CC))/(2.d0*AA)
            IF(DM.LE.0.0.OR.DM.GT.(M1-MC1))THEN
               if(rank.eq.0) WRITE(*,*) ' WARNING: CE DM ',DM,M1,MC1
            ENDIF
            M1 = MAX(M1-DM,MC1)
            EBINDF = M1*(M1-MC1)/(LAMB1*R1)
            DE = EBINDI - EBINDF
            EORBF = EORBI + DE/ALPHA1
            ECC = 0.001D0
         ENDIF
         SEPF = M1*M2/(2.D0*EORBF)
         IF((M1-MC1).LT.1.0D-07)THEN
            M1 = MC1
            CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &                  R1,L1,KW1,MC1,RC1,MENV,RENV,K21)
         ELSEIF(KW1.EQ.4)THEN
            AJ1 = (AJ1 - TSCLS1(2))/TSCLS1(3)
            CALL gntage(MC1,M1,KW1,ZPARS,M01,AJ1)
         ENDIF
         GOTO 30 
      ENDIF
*
* If the secondary star is also giant-like add its envelopes's energy.
*
      IF(KW2.GE.2.AND.KW2.LE.9.AND.KW2.NE.7)THEN
         IF(M2.GT.MC2)THEN
            MENVD = MENV/(M2-MC2)
            RZAMS = RZAMSF(M02)
            LAMB2 = CELAMF(KW,M02,L2,R2,RZAMS,MENVD,LAMBDA)
         ELSE
            LAMB2 = 0.5D0
         ENDIF
         EBINDI = EBINDI + M2*(M2-MC2)/(LAMB2*R2)
*
* Re-calculate the initial orbital energy
*
         EORBI = MC1*MC2/(2.D0*SEP)
      ENDIF
*
*
* Calculate the final orbital energy without coalescence.
*
      EORBF = EORBI + EBINDI/ALPHA1
*
* If the secondary is on the main sequence see if it fills its Roche lobe.
*
      IF(KW2.LE.1.OR.KW2.EQ.7)THEN
         SEPF = MC1*M2/(2.D0*EORBF)
***
* Assume the energy generated by forcing the secondary to
* co-rotate goes into the envelope (experimental).
*        TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(MC1+M2)))
*        OORB = TWOPI/TB
*        DELY = 0.5D0*M2*R2*R2*(OSPIN2*OSPIN2 - OORB*OORB)/3.91D+08
*        DELY = K22*DELY
*        EBINDI = MAX(0.D0,EBINDI - DELY)
*        EORBF = EORBI + EBINDI/ALPHA1
*        SEPF = MC1*M2/(2.D0*EORBF)
***
         Q1 = MC1/M2
         Q2 = 1.D0/Q1
         RL1 = RL(Q1)
         RL2 = RL(Q2)
         IF(RC1/RL1.GE.R2/RL2)THEN
*
* The helium core of a very massive star of type 4 may actually fill
* its Roche lobe in a wider orbit with a very low-mass secondary.
*
            IF(RC1.GT.RL1*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC1/RL1
            ENDIF
         ELSE
            IF(R2.GT.RL2*SEPF)THEN
               COEL = .TRUE.
               SEPL = R2/RL2
            ENDIF
         ENDIF
         IF(COEL)THEN
*
            KW = KTYPE(KW1,KW2) - 100
            MC3 = MC1
            IF(KW2.EQ.7.AND.KW.EQ.4) MC3 = MC3 + M2
*
* Coalescence - calculate final binding energy.
*
            EORBF = MAX(MC1*M2/(2.D0*SEPL),EORBI)
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI)
         ELSE
*
* Primary becomes a black hole, neutron star, white dwarf or helium star.
*
            MF = M1
            M1 = MC1
            CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &                  R1,L1,KW1,MC1,RC1,MENV,RENV,K21)
         ENDIF
      ELSE
*
* Degenerate or giant secondary. Check if the least massive core fills its
* Roche lobe.
*
         SEPF = MC1*MC2/(2.D0*EORBF)
***
* Assume the energy generated by forcing the secondary to
* co-rotate goes into the envelope (experimental).
*        IF(KW2.GE.10.AND.KW2.LE.14)THEN
*           TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(MC1+MC2)))
*           OORB = TWOPI/TB
*           DELY = 0.5D0*M2*R2*R2*(OSPIN2*OSPIN2 - OORB*OORB)/3.91D+08
*           DELY = K3*DELY
*           EBINDI = MAX(0.D0,EBINDI - DELY)
*           EORBF = EORBI + EBINDI/ALPHA1
*           SEPF = MC1*MC2/(2.D0*EORBF)
*        ENDIF
***
         Q1 = MC1/MC2
         Q2 = 1.D0/Q1
         RL1 = RL(Q1)
         RL2 = RL(Q2)
         IF(RC1/RL1.GE.RC2/RL2)THEN
            IF(RC1.GT.RL1*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC1/RL1
            ENDIF
         ELSE
            IF(RC2.GT.RL2*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC2/RL2
            ENDIF
         ENDIF
*
         IF(COEL)THEN
*
* If the secondary was a neutron star or black hole the outcome
* is an unstable Thorne-Zytkow object that leaves only the core.
*
            SEPF = 0.D0
            IF(KW2.GE.13)THEN
               MC1 = MC2
               M1 = MC1
               MC2 = 0.D0
               M2 = 0.D0
               KW1 = KW2
               KW2 = 15
               AJ1 = 0.D0
*
* The envelope mass is not required in this case.
*
               GOTO 30
            ENDIF
*
            KW = KTYPE(KW1,KW2) - 100
            MC3 = MC1 + MC2
*
* Calculate the final envelope binding energy.
*
            EORBF = MAX(MC1*MC2/(2.D0*SEPL),EORBI)
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI)
*
* Check if we have the merging of two degenerate cores and if so
* then see if the resulting core will survive or change form.
*
            IF(KW1.EQ.6.AND.(KW2.EQ.6.OR.KW2.GE.11))THEN
               CALL dgcore(KW1,KW2,KW,MC1,MC2,MC3,EBINDF)
            ENDIF
            IF(KW1.LE.3.AND.M01.LE.ZPARS(2))THEN
               IF((KW2.GE.2.AND.KW2.LE.3.AND.M02.LE.ZPARS(2)).OR.
     &             KW2.EQ.10)THEN
                  CALL dgcore(KW1,KW2,KW,MC1,MC2,MC3,EBINDF)
                  IF(KW.GE.10)THEN
                     KW1 = KW
                     M1 = MC3
                     MC1 = MC3
                     IF(KW.LT.15) M01 = MC3
                     AJ1 = 0.D0
                     MC2 = 0.D0
                     M2 = 0.D0
                     KW2 = 15
                     GOTO 30
                  ENDIF
               ENDIF
            ENDIF
*
         ELSE
*
* The cores do not coalesce - assign the correct masses and ages.
*
            MF = M1
            M1 = MC1
            CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &                  R1,L1,KW1,MC1,RC1,MENV,RENV,K21)
            MF = M2
            KW = KW2
            M2 = MC2
            CALL star(KW2,M02,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
            CALL hrdiag(M02,AJ2,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS,
     &                  R2,L2,KW2,MC2,RC2,MENV,RENV,K22)
         ENDIF
      ENDIF
*
      IF(COEL)THEN
         MC22 = MC2
         IF(KW.EQ.4.OR.KW.EQ.7)THEN
* If making a helium burning star calculate the fractional age 
* depending on the amount of helium that has burnt.
            IF(KW1.LE.3)THEN
               FAGE1 = 0.D0
            ELSEIF(KW1.GE.6)THEN
               FAGE1 = 1.D0
            ELSE
               FAGE1 = (AJ1 - TSCLS1(2))/(TSCLS1(13) - TSCLS1(2))
            ENDIF
            IF(KW2.LE.3.OR.KW2.EQ.10)THEN
               FAGE2 = 0.D0
            ELSEIF(KW2.EQ.7)THEN
               FAGE2 = AJ2/TM2
               MC22 = M2
            ELSEIF(KW2.GE.6)THEN
               FAGE2 = 1.D0
            ELSE
               FAGE2 = (AJ2 - TSCLS2(2))/(TSCLS2(13) - TSCLS2(2))
            ENDIF
         ENDIF
      ENDIF
*
* Now calculate the final mass following coelescence.  This requires a
* Newton-Raphson iteration.
*
      IF(COEL)THEN
*
* Calculate the orbital spin just before coalescence.
*
         TB = (SEPL/AURSUN)*SQRT(SEPL/(AURSUN*(MC1+MC2)))
         OORB = TWOPI/TB
*
         XX = 1.D0 + ZPARS(7)
         IF(EBINDF.LE.0.D0)THEN
            MF = MC3
            GOTO 20
         ELSE
            CONST = ((M1+M2)**XX)*(M1-MC1+M2-MC22)*EBINDF/EBINDI
         ENDIF
*
* Initial Guess.
*
         MF = MAX(MC1 + MC22,(M1 + M2)*(EBINDF/EBINDI)**(1.D0/XX))
   10    DELY = (MF**XX)*(MF - MC1 - MC22) - CONST
*        IF(ABS(DELY/MF**(1.D0+XX)).LE.1.0D-02) GOTO 20
         IF(ABS(DELY/MF).LE.1.0D-03) GOTO 20
         DERI = MF**ZPARS(7)*((1.D0+XX)*MF - XX*(MC1 + MC22))
         DELMF = DELY/DERI
         MF = MF - DELMF
         GOTO 10
*
* Set the masses and separation.
*
   20    IF(MC22.EQ.0.D0) MF = MAX(MF,MC1+M2)
         M2 = 0.D0
         M1 = MF
         KW2 = 15
*
* Combine the core masses.
*
         IF(KW.EQ.2)THEN
            CALL star(KW,M1,M1,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
            IF(GB(9).GE.MC1)THEN
               M01 = M1
               AJ1 = TM2 + (TSCLS2(1) - TM2)*(AJ1-TM1)/(TSCLS1(1) - TM1)
               CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            ENDIF
         ELSEIF(KW.EQ.7)THEN
            M01 = M1
            CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            AJ1 = TM1*(FAGE1*MC1 + FAGE2*MC22)/(MC1 + MC22)
         ELSEIF(KW.EQ.4.OR.MC2.GT.0.D0.OR.KW.NE.KW1)THEN
            IF(KW.EQ.4) AJ1 = (FAGE1*MC1 + FAGE2*MC22)/(MC1 + MC22)
            MC1 = MC1 + MC2
            MC2 = 0.D0
*
* Obtain a new age for the giant.
*
            CALL gntage(MC1,M1,KW,ZPARS,M01,AJ1)
            CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
         ENDIF
         CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &               R1,L1,KW,MC1,RC1,MENV,RENV,K21)
         JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
         JSPIN2 = 0.D0
         KW1 = KW
         ECC = 0.001D0
      ELSE
*
* Determine if any eccentricity remains in the orbit. 
*
*        IF(ECCFLG.AND.EORBF.LT.ECIRC)THEN
*           ECC = SQRT(1.D0 - EORBF/ECIRC)
*        ELSE
*           ECC = 0.001D0
*        ENDIF
*
* Set both cores in co-rotation with the orbit on exit of CE,
*
         TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(M1+M2)))
         OORB = TWOPI/TB
*        JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
*        JSPIN2 = OORB*(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
*
* or, leave the spins of the cores as they were on entry.
* Tides will deal with any synchronization later.
*
         JSPIN1 = OSPIN1*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
         JSPIN2 = OSPIN2*(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
      ENDIF
   30 SEP = SEPF
      AJ1 = MAX(AJ1,1.0D-10)
      AJ2 = MAX(AJ2,1.0D-10)
*     AJ1 = MAX(AJ1,1.0D-10*TM1)
*     AJ2 = MAX(AJ2,1.0D-10*TM2)
*
      if(rank.eq.0)then
      WRITE(38,68)KW1,M01,M1,MC1,AJ1
 68   FORMAT(' COMENV NEW *1:',I4,3F7.2,F9.2)
      WRITE(38,69)KW2,M02,M2,MC2,AJ2
 69   FORMAT(' COMENV NEW *2:',I4,3F7.2,F9.2)
      end if
*
      RETURN
      END
***
