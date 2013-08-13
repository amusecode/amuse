      SUBROUTINE XTRNL0
*
*
*       External force initialization.
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
*
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
*
*       Check option for cluster in circular galactic orbit.
      IF (KZ(14).NE.1) GO TO 20
*
*       Specify Oort's constants (units of km/sec/kpc).
      A = 14.4
      B = -12.0
*       Adopt local density from Gilmore & Kuijken (solar mass/pc**3).
      RHON = 0.11
*       Convert rotation constants to units of cm/sec/pc.
      A = 100.0*A
      B = 100.0*B
*
*       Define GM in cgs units and length scale in pc.
      GM = 6.67D-08*1.989D+33
      PC = 3.0857D+18
*
*       Specify the tidal term in star cluster units (solar mass & pc).
      TIDAL(1) = 4.0*A*(A - B)*(PC/GM)
*
*       Initialize the Y-component to zero.
      TIDAL(2) = 0.0
*
*       Specify the vertical force gradient.
      TIDAL(3) = -(2.0*TWOPI*RHON + 2.0*(A - B)*(A + B)*(PC/GM))
*
      FAC = 1.0E-10/(PC/GM)
      if(rank.eq.0)
     &WRITE (6,5)  ZMBAR*ZMASS, FAC*TIDAL(1), FAC*TIDAL(3), PC/GM
    5 FORMAT (/,12X,'TOTAL MASS =',F8.1,'  TIDAL(1&3) =',1P,2E10.2,
     &              '  PC/GM =',E10.2)
*
*       Adopt twice the angular velocity for Coriolis terms.
      TIDAL(4) = 2.0*(A - B)*SQRT(PC/GM)
*
*       Define time scale in seconds and velocity scale in km/sec.
      TSCALE = SQRT(PC/GM)*PC
      VSTAR = 1.0E-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSCALE = TSCALE/(3.15D+07*1.0D+06)
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      DO 10 K = 1,3
          TIDAL(K) = TIDAL(K)*RBAR**3/ZMBAR
   10 CONTINUE
      TIDAL(4) = TIDAL(4)*SQRT(RBAR**3/ZMBAR)
      TSCALE = TSCALE*SQRT(RBAR**3/(ZMASS*ZMBAR))
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
*       Consider alternatives: circular point-mass orbit or 3D galaxy model.
   20 ZMTOT = ZMASS*ZMBAR
      IF (KZ(14).EQ.2) THEN
*
*       Read galaxy mass and central distance (solar units and kpc).
          if(amusein.eq.0) READ (5,*)  GMG, RG0
#if MPIINIT
      CALL MPI_BCAST(GMG,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RG0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif

*       Set circular velocity in km/sec and angular velocity in cgs units.
          VG0 = 1.0D-05*SQRT(GMG/(1000.0*RG0))*SQRT(GM/PC)
          OMEGA = 100.0*VG0/RG0
*
*       Obtain King tidal radius in pc (eq. (9) of Fukushige & Heggie, 1995).
          RT = (ZMTOT/(3.0*GMG))**0.3333*(1000.0*RG0)
*
          IF (RTIDE.GT.0.0) THEN
*       Determine RBAR (N-body units) from RT (pc) and King model (see SCALE).
              IF(KZ(22).EQ.2) RBAR = RT/RTIDE
          ELSE
              RTIDE = RT/RBAR
          END IF
*
*       Convert from cgs to N-body units.
          OMEGA = OMEGA*SQRT(PC/GM)*SQRT(RBAR**3/ZMBAR)
*
*       Specify the galactic parameters for equations of motion.
          TIDAL(1) = 3.0*OMEGA**2
          TIDAL(2) = 0.0D0
          TIDAL(3) = -OMEGA**2
          TIDAL(4) = 2.0*OMEGA
          GMG = GMG/ZMTOT
*
*       Check re-scaling units to current RBAR (i.e. TSCALE, TSTAR & VSTAR).
          IF (KZ(22).EQ.2) THEN
              CALL UNITS
          END IF
*
          WRITE (6,35)  GMG, RG0, OMEGA, RTIDE, RBAR
*
*       Treat the general case of 3D orbit for point-mass, disk and/or halo.
      ELSE IF (KZ(14).EQ.3) THEN
*
*       Read all parameters (NB! Do not confuse with Oort's constants A, B).
          if(amusein.eq.0)then
          READ (5,*)  GMG, DISK, A, B, VCIRC, RCIRC, (RG(K),K=1,3),
     &                                               (VG(K),K=1,3)
          end if
#if MPIINIT
      CALL MPI_BCAST(GMG,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DISK,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(A,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(B,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(VCIRC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RCIRC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RG,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(VG,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif 
*
*       Specify planar motion from SEMI & ECC for no disk & halo if VZ = 0.
          IF (DISK + VCIRC.EQ.0.0.AND.VG(3).EQ.0.0D0) THEN
              RAP = RG(1)
              ECC = RG(2)
              SEMI = RAP/(1.0 + ECC)
              VG2 = GMG/(1000.0*SEMI)*(1.0 - ECC)/(1.0 + ECC)
              DO 25 K = 1,3
                  RG(K) = 0.0
                  VG(K) = 0.0
   25         CONTINUE
*       Initialize 2D orbit with given eccentricity at apocentre.
              RG(1) = RAP
              VG(2) = 1.0D-05*SQRT(VG2)*SQRT(GM/PC)
          END IF
*
*       Convert from kpc and km/sec to N-body units.
          DO 30 K = 1,3
              RG(K) = 1000.0*RG(K)/RBAR
              VG(K) = VG(K)/VSTAR
   30     CONTINUE
*
*       Define the angular velocity (z-component) and mass in N-body units.
          R02 = RG(1)**2 + RG(2)**2
          OMEGA = (RG(1)*VG(2) - RG(2)*VG(1))/R02
          TIDAL(4) = 2.0*OMEGA
          GMG = GMG/ZMTOT
*       Adopt a tidal radius of 10 units unless specified by routine SCALE.
          IF (RTIDE.EQ.0.0D0) RTIDE = 10.0
*
          if(rank.eq.0)
     &    WRITE (6,35)  GMG, SQRT(R02), OMEGA, RTIDE, RBAR
   35     FORMAT (/,12X,'POINT-MASS MODEL:    MG =',1P,E9.1,
     &                  '  RG =',E9.1,'  OMEGA =',E9.1,
     &                  '  RT =',0P,F6.2,'  RBAR =',F6.2)
*       Define disk and/or logarithmic halo parameters in N-body units.
          IF (DISK.GT.0.0D0) THEN
              DISK = DISK/ZMTOT
              A = 1000.0*A/RBAR
              B = 1000.0*B/RBAR
              WRITE (6,40)  DISK, A, B
   40         FORMAT (/,12X,'DISK MODEL:    MD =',1P,E9.1,
     &                                   '  A =',E9.1,'  B =',E9.1)
          END IF
*
*       Determine local halo velocity from total circular velocity.
          IF (VCIRC.GT.0.0D0) THEN
              VCIRC = VCIRC/VSTAR
              RCIRC = 1000.0*RCIRC/RBAR
              A2 = RCIRC**2 + (A + B)**2
              V02 = VCIRC**2 - (GMG/RCIRC + DISK*RCIRC**2/A2**1.5)
              IF (V02.LT.0.0D0) THEN
                  WRITE (6,45)  V02, 0.001*RCIRC*RBAR
   45             FORMAT (' ',' NEGATIVE HALO VELOCITY!    V02 RCIRC ',
     &                                                     1P,2E10.2)
                  STOP
              END IF
*       Specify the corresponding scale length of logarithmic halo.
              RL2 = RCIRC**2*(VCIRC**2 - V02)/V02
*       Define the asymptotic circular velocity due to halo.
              V02 = VCIRC**2
*
*       Include table of circular velocity on unit #17 (km/sec & kpc).
              RI = 1000.0/RBAR
              DR = 1000.0/RBAR
              DO 60 K = 1,30
                  RI2 = RI**2
                  A2 = RI2 + (A + B)**2
                  VCIRC2 = GMG/SQRT(RI2) + DISK*RI2/A2**1.5 +
     &                                     V02*RI2/(RL2 + RI2)
                  WRITE (17,50)  SQRT(VCIRC2)*VSTAR, RI*RBAR/1000.0
   50             FORMAT (' CIRCULAR VELOCITY:    VC R ',F7.1,F7.2)
                  RI = RI + DR
   60         CONTINUE
              CALL FLUSH(17)
*
              A2 = R02 + (A + B)**2
              VCIRC2 = GMG/SQRT(R02) + DISK*R02/A2**1.5 +
     &                                 V02*R02/(RL2 + R02)
              VCIRC = SQRT(VCIRC2)*VSTAR
              WRITE (6,62)  VCIRC, SQRT(R02)/1000.0, SQRT(RL2)/1000.0
   62         FORMAT (/,12X,'CIRCULAR VELOCITY:    VC RG RL',F7.1,2F7.2)
          ELSE
              V02 = 0.0
          END IF
*
*       Initialize F & FDOT of reference frame (point-mass galaxy is OK).
          CALL GCINIT
*
          WRITE (6,65)  (RG(K),K=1,3), (VG(K),K=1,3), SQRT(V02)
   65     FORMAT (/,12X,'SCALED ORBIT:    RG =',1P,3E10.2,
     &                                '  VG = ',3E10.2,'  V0 =',0P,F6.1)
      END IF
*
*       Include Plummer potential for 2D and 3D (set MP = 0 if not needed).
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
*       Check input for Plummer potential.
          if(amusein.eq.0)then
          READ (5,*)  MP, AP2, MPDOT, TDELAY
          end if
#if MPIINIT
      CALL MPI_BCAST(MP,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(AP2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(MPDOT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TDELAY,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          if(rank.eq.0)WRITE (6,70)  MP, AP2, MPDOT, TDELAY
   70     FORMAT (/,12X,'PLUMMER POTENTIAL:    MP =',F7.3,'  AP =',F6.2,
     &                  '  MPDOT =',F7.3,'  TDELAY =',F5.1)
          MP0 = MP
          AP2 = AP2**2
*       Rescale velocities by including the Plummer & galactic virial energy.
          IF (ZKIN.GT.0.0D0) THEN
*       Note that QVIR = Q is saved in routine SCALE and VIR < 0 with GPU.
#ifdef PARALLEL
              CALL ENERGY_MPI
#else
              CALL ENERGY
#endif
              VIR = POT - VIR
              QV = SQRT(QVIR*VIR/ZKIN)
              DO 80 I = 1,N
                  DO 78 K = 1,3
                      XDOT(K,I) = XDOT(K,I)*QV
   78             CONTINUE
   80         CONTINUE
          END IF
      ELSE
          MP = 0.0
      END IF
      RTIDE0 = RTIDE
*
*       Define tidal radius in scaled units for linearized field.
      IF (KZ(14).LE.2) THEN
          RTIDE = (ZMASS/TIDAL(1))**0.3333
          WRITE (6,90)  (TIDAL(K),K=1,4), TSCALE, RTIDE
   90     FORMAT (/,12X,'TIDAL PARAMETERS:  ',1P,4E10.2,
     &                  '  TSCALE =',E9.2,'  RTIDE =',0P,F6.2,/)
      END IF
*
      RETURN
*
      END
