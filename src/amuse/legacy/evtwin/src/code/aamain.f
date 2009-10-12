! Program EV
! An incarnation of Peter Eggleton's TWIN binary evolution code
! This version of the code can be run as a standalone program, meaning
! it doesn't need symbolic links or environment variables, although it will
! continue to use those if found.
!
! Usage:
!   ev [name] [metalicity] [evpath]
!
! Stores output files with the basename "name". If unspecified, the basename
! defaults to "mod".
! The metalicity can be specified on the command line as a pure number, or 
! through the environment variable Z. If specified on the command line, it 
! overrides the environment variable. If not specified at all, the metalicity 
! defaults to solar (0.02).
! The evpath parameter identifies the location of the code datafiles, which are
! expected in evpath/input/. Defaults to the current directory if not specified.
! Settings are read from files init.run and init.dat in the current directory.
! If these don't exist, "name".dat and "name".run will be used instead.
! An alternative starting model can be specified (by name) using the STARTFILE
! parameter in the .run file.
!
! Examples:
! Run the evolution code with default settings:
!   ev
! Run the code, saving the output file as "sun"
!   ev sun
! As above, but change the metalicity to 0.0001:
!   ev sun 0001
! 
C The working of this programme is described in a file `writeup.tex';
      USE MESH
      USE MESH_ENC
      USE SETTINGS
      USE CONSTANTS
      USE FILE_EXISTS_MODULE
      IMPLICIT REAL*8 (A-H, L-Z)
      DOUBLE PRECISION :: SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      DOUBLE PRECISION :: SM, DTY, AGE, PER, BMS, ECC, P, ENC
      DOUBLE PRECISION :: SM1, DTY1, AGE1, PER1, BMS1, ECC1, P1, ENC1
      DOUBLE PRECISION :: SM2, DTY2, AGE2, PER2, BMS2, ECC2, P2, ENC2
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /TN1/ SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      COMMON /T0 / SM, DTY, AGE, PER, BMS, ECC, P, ENC
      COMMON /T1 / SM1, DTY1, AGE1, PER1, BMS1, ECC1, P1, ENC1
      COMMON /T2 / SM2, DTY2, AGE2, PER2, BMS2, ECC2, P2, ENC2
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      REAL*8 H1(NVAR,NM), H2(NVAR,NM)
      REAL*8 DH1(NVAR,NM), DH2(NVAR,NM)
      CHARACTER*500 :: STARTFILE, ZAMSFILE
      COMMON /CINIT_RUN/ ISB, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, JMX, STARTFILE, ZAMSFILE
!MvdS Added to create different filenames for different loops
      INTEGER NML,NQL,NXL,DUMMY,NARG,IARGC,FLEN,UNITS(12),NOUT(2)
      CHARACTER FNAME*500,FNAME2*6,COMMAND*500,SUFFIX(12)*6,TMPSTR*80
      INTEGER DATEVAL(8)
      CHARACTER EVPATH*500, ZSTR*8
      INTEGER, PARAMETER :: N_INP_FILES = 13
      CHARACTER*500 :: INPUTFILENAMES(N_INP_FILES)
      INTEGER :: INPUTUNITS(N_INP_FILES)     = (/12, 24, 16, 18, 19, 20, 21, 26, 63, 22, 23, 41, 42/)
      INTEGER :: INPUT_REQUIRED(N_INP_FILES) = (/ 0,  0,  1,  1,  1,  0,  1,  1, -1,  1,  1,  0,  1/)
      SUFFIX = (/'.out1 ','.io12 ','.log  ','.out  ','.last1','.last2',
     &'.mod  ','.plt1 ','.mdl1 ','.out2 ','.plt2 ','.mdl2 '/)
      UNITS = (/1,3,8,9,13,14,15,31,33,2,32,34/)
      NOUT = (/9,12/) !N of output files for ISB=1,2
!MvdS Added to create different filenames for different loops

! Setup metalicity and evolution code directory.
! In order of precedence: command line option, environment variable, default value
      IF (IARGC() >= 2) THEN
         CALL GETARG(2,ZSTR)
      ELSE
         CALL GETENV("Z", ZSTR)
      END IF
      IF (LEN(TRIM(ZSTR)) == 0) ZSTR = "02"
      IF (IARGC() >= 3) THEN
         CALL GETARG(3,EVPATH)
      ELSE
         CALL GETENV("evpath", EVPATH)
      END IF
      IF (LEN(TRIM(EVPATH)) == 0) EVPATH = "."

! Open input files; this replaces the symbolic links to fort.nnn files
! This is actually not the best way to do this (better would be to specify
! the filenames from a configuration file), but it'll do.
! Opens the file directly by name, unless the associated fort.nnn file exists
      INPUTFILENAMES(1)=TRIM(EVPATH)//"/input/zahb"//TRIM(ZSTR)//".mod"
      INPUTFILENAMES(2)=TRIM(EVPATH)//"/input/zahb"//".dat"
      INPUTFILENAMES(3)=TRIM(EVPATH)//"/input/zams/zams"//TRIM(ZSTR)//".mod"
      INPUTFILENAMES(4)=TRIM(EVPATH)//"/input/zams/zams"//TRIM(ZSTR)//".out"
      INPUTFILENAMES(5)=TRIM(EVPATH)//"/input/zams/zams"//TRIM(ZSTR)//".mas"
      INPUTFILENAMES(6)=TRIM(EVPATH)//"/input/metals/z"//TRIM(ZSTR)//"/phys.z"//TRIM(ZSTR)
      INPUTFILENAMES(7)=TRIM(EVPATH)//"/input/lt2ubv.dat"
      INPUTFILENAMES(8)=TRIM(EVPATH)//"/input/nucdata.dat"
      INPUTFILENAMES(9)=TRIM(EVPATH)//"/input/mutate.dat"
      INPUTFILENAMES(10)='init.dat'
      INPUTFILENAMES(11)='init.run'
      INPUTFILENAMES(12)=TRIM(EVPATH)//"/input/COtables/COtables_z"//TRIM(ZSTR)
      INPUTFILENAMES(13)=TRIM(EVPATH)//"/input/physinfo.dat"

! If init.dat does not exist, try name.dat; likewise for init.run
      IF(IARGC() > 0) THEN
        CALL GETARG(1,FNAME)
        IF ( (.NOT. FILE_EXISTS(INPUTFILENAMES(10))) .AND. 
     &             (FILE_EXISTS(TRIM(FNAME)//".dat")) )
     &   INPUTFILENAMES(10)=TRIM(FNAME)//".dat"

        IF ( (.NOT. FILE_EXISTS(INPUTFILENAMES(11))) .AND. 
     &             (FILE_EXISTS(TRIM(FNAME)//".run")) )
     &   INPUTFILENAMES(11)=TRIM(FNAME)//".run"
      END IF

      DO I=1, N_INP_FILES
         WRITE (FNAME, '("fort.",I2)') INPUTUNITS(I)
         IF (FILE_EXISTS(INPUTFILENAMES(I)) .AND. .NOT. FILE_EXISTS(FNAME)) THEN
            !PRINT *, "Opening ",TRIM(INPUTFILENAMES(I))," for ", TRIM(FNAME)
            OPEN(UNIT = INPUTUNITS(I), ACTION="READ", FILE=INPUTFILENAMES(I))
         END IF
         ! Check if files exist and act appropriately
         ! If the file is required (INPUT_REQUIRED>0), abort on error
         ! If the file is optional (INPUT_REQUIRED==0), give a warning
         ! If the file is probably unneeded (INPUT_REQUIRED<0), do nothing
         IF (.NOT. (FILE_EXISTS(INPUTFILENAMES(I)) .OR. FILE_EXISTS(FNAME))) THEN
            IF (INPUT_REQUIRED(I) > 0) THEN
               WRITE (0, *) 'Required input file ', TRIM(INPUTFILENAMES(I)), ' (', TRIM(FNAME), ') not found'
               STOP
            ELSEIF (INPUT_REQUIRED(I) == 0) THEN
               WRITE (0, *) 'Warning: input file ', TRIM(INPUTFILENAMES(I)), ' (', TRIM(FNAME), ') not found'
            END IF
         END IF
      END DO
      
!MvdS Check for command line arguments and determine the output filename
      FNAME='mod'
      FLEN=3
      NARG = IARGC()
      IF(NARG.GT.0) THEN
        CALL GETARG(1,FNAME)
        FLEN = LEN(TRIM(FNAME))
      ENDIF
      
! We need to have a fort.11
      WRITE (11, *) 0
      CLOSE (11)

! We need a scratch file for output we don't want or need
      IF (.NOT. FILE_EXISTS('fort.25'))
     &   OPEN (UNIT=25, ACTION='WRITE', STATUS = 'SCRATCH')
      
      CALL SETSUP
   92 FORMAT (1X, 1P, 40D23.15, 0P)
   95 FORMAT (1X, 1P, 8D23.15, 0P, 6I6) !MvdS 5I5->5I6 to allow >10k models
   99 FORMAT (3(1P, E16.9, 5E10.3, 0P, F8.3, 7F8.5, 2F8.4, F9.5, /), 
     :         (1P, E16.9, 5E10.3, 0P, F8.3, 7F8.3, 2F8.4, F9.5, /), 
     :          1P, E16.9, 5E10.3, 0P, 11F8.3, I6, I2)
      READ (19, *) MLO, DM, MHI, KDM 
c Read data for a triple set of loops, in mass, mass ratio, and period
      CALL READ_INIT_RUN(23)
c Create short (`pruned') summary of ZAMS models from long file
      CALL PRUNER ( 18, 17, ISB )
      CLOSE (17)
      CLOSE (18)

! Open output files; this replaces the symbolic links and fort.nnn files (EG)
! Only use named files if no fort.nnn file for the same unit exists
      DO I=1, NOUT(ISB)
         IF (UNITS(I)<10) THEN
            WRITE (COMMAND, '("fort.",I1)') UNITS(I)
         ELSE
            WRITE (COMMAND, '("fort.",I2)') UNITS(I)
         END IF
         IF (.NOT. FILE_EXISTS(COMMAND)) THEN
            !PRINT *, "Opening ",FNAME(1:FLEN)//TRIM(SUFFIX(I))," for ", TRIM(COMMAND)
            OPEN(UNIT = UNITS(I), FILE=TRIM(FNAME)//SUFFIX(I))
         END IF
      END DO
      OPEN(UNIT = 29, FILE=TRIM(FNAME)//'.mas')
      IF (MUTATE) OPEN(UNIT = 65, FILE=TRIM(FNAME)//'.pmutate')

c Write out termination (JO) codes on fort.8 for later reference
      WRITE (8,*) ' -2 -- BEGINN -- requested mesh too large'
      WRITE (8,*) ' -1 -- STAR12 -- no timesteps required'
      WRITE (8,*) '  0 -- STAR12 -- finished required timesteps OK'
      WRITE (8,*) '  1 -- SOLVER -- failed; backup, reduce timestep'
c JO = 1 won't stop the code, but it may lead to JO = 2
      WRITE (8,*) '  2 -- BACKUP -- tstep reduced below limit; quit'
      WRITE (8,*) '  3 -- NEXTDT -- *2 evolving beyond last *1 model' 
      WRITE (8,*) '  4 -- PRINTB -- *1 rstar exceeds rlobe by limit'
      WRITE (8,*) '  5 -- PRINTB -- age greater than limit'
      WRITE (8,*) '  6 -- PRINTB -- C-burning exceeds limit'
      WRITE (8,*) '  7 -- PRINTB -- *2 radius exceeds rlobe by limit'
      WRITE (8,*) '  8 -- PRINTB -- close to He flash'
      WRITE (8,*) '  9 -- PRINTB -- massive (>1.2 msun) deg. C/O core'
      WRITE (8,*) ' 10 -- PRINTB -- |M1dot| exceeds limit' 
      WRITE (8,*) ' 11 -- NEXTDT -- impermissible FDT for *2' 
      WRITE (8,*) ' 14 -- PRINTB -- funny compos. distribution'
      WRITE (8,*) ' 15 -- STAR12 -- terminated by hand'
      WRITE (8,*) ' 16 -- MAIN   -- ZAHB didnt converge'
      WRITE (8,*) ' 51 -- PRINTB -- end of MS (core H below limit)'
      WRITE (8,*) ' 52 -- PRINTB -- Radius exceeds limit'
      WRITE (8,*) ' 53 -- PRINTB -- Convergence to target model reached minimum'
      WRITE (8,*) ' 12, 22, 32 -- as 2, core H, He, C non-zero, resp.' 
      WRITE (8,*) ''
      WRITE (8,*) 'CODE*1 C*2  log(Mi)   log(Qi)   log(Pi)',
     &' Modl  JMOD 27-JOP:'
      JOP = 14 
      IF ( ISB.EQ.1 ) KTW = 1
      IF ( ISB.EQ.1.OR.KTW.EQ.2 ) KP = KPT
      IF ( IP1.LE.15.OR.IP2.LE.14 ) KML = 1
      IF ( IP1.LE.14.OR.IP2.LE.15 ) KQL = 1
      IF ( IP1.LE.14.OR.IP2.LE.14 ) KXL = 1
      IF ( IP1.LE.14 ) IP2 = IP1
      
!MvdS Create a list of model properties for this run, including date and pwd
      OPEN(UNIT=40, STATUS='UNKNOWN', FORM='FORMATTED', FILE=TRIM(FNAME)//'.list')
      WRITE(40,'(A1)')' ' 
! Get date and time for list file
      CALL DATE_AND_TIME(TMPSTR, TMPSTR, TMPSTR, DATEVAL)
      WRITE(40, '(I4"-"I2"-"I2, " ", I2":"I2,":"I2, " UT"A5)') 
     &   DATEVAL(1:3), DATEVAL(5:7), TMPSTR
! Write hostname
      CALL GETENV("HOSTNAME", TMPSTR);
      WRITE(40,'(A80)') TMPSTR      
! Write current working directory
      CALL GETENV("PWD", TMPSTR);
      WRITE(40,'(A80)') TMPSTR
! Write names of files used
      DO I=1, N_INP_FILES
         IF (FILE_EXISTS(INPUTFILENAMES(I)) .AND. .NOT. FILE_EXISTS(FNAME)) THEN
            WRITE (40, *) "Using ",TRIM(INPUTFILENAMES(I))," for unit", INPUTUNITS(I)
         END IF
      END DO
      IF ( LEN(TRIM(STARTFILE))>0 ) THEN
         WRITE (40, *) "Using initial model ",TRIM(STARTFILE)," for unit", IP1
      END IF
      WRITE(40,'(A1)')' ' 
      WRITE(40,'(16X,A4,4(8X,A3))')'File','M1i','M2i',' Qi',' Pi'
      CLOSE(40)

!BEGIN LOOP      
c Cycle on *1's mass.
      NML = 0
      ML = ML1
      DO WHILE (ML < ML1 + (KML - 1)*DML + 0.1D-4)
         NML = NML + 1
c Locate the *1 model, optionally in the ZAMS (fort.16) library
         IF ( IP1==16 .AND. LEN(TRIM(ZAMSFILE))>0 ) THEN
            ! Override the ZAMS file from which the starting model is to be
            ! obtained; Necessary for stars > 200Msun
            OPEN(UNIT=40, FILE=TRIM(ZAMSFILE))
            IM1 = 1
            CALL LOAD_STAR_MODEL(40,IM1, H1, DH1, SM1,DTY1,AGE1,PER1,BMS1,ECC1,P1,ENC 1,KH1,JN1,JF1)
            CLOSE(40)
         ELSE
            IF ( IP1.EQ.16 ) IM1 = (ML - MLO)/DM + 1.501D0
            CALL LOAD_STAR_MODEL(IP1,IM1, H1, DH1, SM1,DTY1,AGE1,PER1,BMS1,ECC1,P1,ENC1,KH1,JN1,JF1)
         END IF
         IF ( JMX.EQ.0.AND.IP1.EQ.16 ) SM1 = 10.0D0**ML
         IF ( IP1.EQ.16 ) THEN
c Locate some parameters relating to *1 in the fort.18 library.
c TN = nuclear timescale, PERC = log period for RLOF on ZAMS
            DO J = 1, KDM*(IM1 - 1) + 1
               READ (17, 99) D, D, D, TN, D, PERC, (D, I = 1,79), JJ  
            END DO
            PERC = 10.0D0**PERC
            CLOSE (17)
         END IF
c Cycle on mass ratio. 
         NQL = 0
         QL = QL1
         DO WHILE (QL < QL1 + (KQL - 1)*DQL + 0.1D-4)
            NQL = NQL + 1
            SM2 = SM1/10.0D0**QL
            IF ( ISB.EQ.2 ) THEN
c Locate the *2 model, optionally in the ZAMS (fort.16) library
               IF ( IP2.EQ.16 ) IM2 = (ML - QL - MLO)/DM + 1.501D0
               CALL LOAD_STAR_MODEL(IP2,IM2, H2, DH2, SM2,DTY2,AGE2,PER2,BMS2,ECC2,P2,ENC2,KH2,JN2,JF2)
               IF ( JMX.EQ.0.AND.IP1.EQ.16 ) SM2 = SM1/10.0D0**QL
c Store total binary mass in case we're setting up the binary from previous models
               IF ( IP1 /= 16 .AND. IP2 /= 16 ) THEN
                  BMS1 = SM1 + SM2
                  BMS2 = SM1 + SM2
               END IF
            END IF
c Cycle on initial period
            NXL = 0
            XL = XL1
            DO WHILE (XL < XL1 + (KXL - 1)*DXL + 0.1D-4)
               NXL = NXL + 1

!BEGIN ACTUAL LOOP
!MvdS Make a list of models in each file	       
               WRITE(FNAME2,'(3I2.2)')NML,NQL,NXL
               OPEN(UNIT=40, STATUS='UNKNOWN', POSITION='APPEND', 
     &         FORM='FORMATTED', FILE=TRIM(FNAME)//'.list')
               IF(KML*KQL*KXL.EQ.1) THEN
                  WRITE(40,'(A20,4(2X,ES9.3))')TRIM(FNAME),
     &           10**ML,10**(ML-QL),10**QL,10**XL
               ELSE
                  WRITE(40,'(A20,4(2X,ES9.3))')TRIM(FNAME)//FNAME2,
     &           10**ML,10**(ML-QL),10**QL,10**XL
               ENDIF
               CLOSE(40)


c For ZAMS models, replace SMn, ..., JMn by correct values
               IF ( IP1 == 16 ) THEN
                  AGE1 = 0.0D0 
                  PER1 = 0.0D0 
                  BMS1 = 0.0D0 
                  ECC1 = 0.0D0 
                  P1   = H1(13,1)
                  ENC1 = 0.0D0 
                  DTY1 = 1.0D-4*TN
                  BMS1 = SM1 + SM2
                  PER1 = PERC*10.0D0**XL
                  IF ( KR.EQ.1 ) P1 = PERC*10.0D0**ROT
                  IF ( KR.EQ.2 ) P1 = DMAX1(PERC, PER1*10.0D0**ROT)
                  IF ( KR.GE.3 ) P1 = PER1 !SdM almost sync. rotation
                  ECC1 = EX
                  JM1 = 0
               END IF
               IF ( IP2 == 16 ) THEN
                  DTY2 = DTY1
                  AGE2 = AGE1
                  PER2 = PER1 
                  BMS2 = BMS1 
                  ECC2 = ECC1
                  P2   = P1
                  ENC2 = ENC1
                  JM2 = 0
               END IF
c Optionally replace SM, ..., JM by values AM, DTY, .., JMX from fort.23
               IF ( KR.EQ.3 ) PER = PER*((1.0D0 - ECC**2)/
     :                                   (1.0D0 - EX**2))**1.5D0
               IF ( KR.EQ.3 ) ECC = EX 

               IF (SM .GE.0.0D0)  SM1  = SM
               IF (DTY .GE.0.0D0) DTY1 = DTY
               IF (AGE .GE.0.0D0) AGE1 = AGE
               IF (PER .GE.0.0D0) PER1 = PER
               IF (BMS .GE.0.0D0) BMS1 = BMS
               IF (ECC .GE.0.0D0) ECC1 = ECC
               IF (P .GE.0.0D0)   P1   = P
               IF (ENC .GE.0.0D0) ENC1 = ENC

               IF (SM .GE.0.0D0)  SM2  = SM
               IF (DTY .GE.0.0D0) DTY2 = DTY
               IF (AGE .GE.0.0D0) AGE2 = AGE
               IF (PER .GE.0.0D0) PER2 = PER
               IF (BMS .GE.0.0D0) BMS2 = BMS
               IF (ECC .GE.0.0D0) ECC2 = ECC
               IF (P .GE.0.0D0)   P2   = P
               IF (ENC .GE.0.0D0) ENC2 = ENC 

               IF ( JMX.GE.0 ) JM1 = JMX
               IF ( JMX.GE.0 ) JM2 = JMX
               SM2 = BMS2 - SM1
               IF ( KTW.EQ.1 ) P2 = 1.0D6
c Store the *1 and *2 initial models in one file (fort.14)
               WRITE (JOP, 95)  SM1,DTY1,AGE1,PER1,BMS1,ECC1,P1,ENC1,
     $              KH1, KP, JM1, 1, JN1, JF1
               DO K = 1, KH1
                  WRITE (JOP, 92) (H1(J, K), J = 1, JN1)
               END DO
               IF (IAND(JF1, 4) == 4) THEN
                  DO K = 1, KH1
                     WRITE (JOP, 92) (DH1(J, K), J = 1, JN1)
                  END DO
               END IF
               IF ( ISB /= 1 ) THEN
                  KPP = 5000
                  IF ( KTW.EQ.2 ) KPP = KP
                  WRITE (JOP, 95)  SM2,DTY2,AGE2,PER2,BMS2,ECC2,P2,ENC2,
     $                 KH2, KPP, JM2, 2, JN2, JF2                
                  DO K = 1, KH2
                     WRITE (JOP, 92) (H2(J, K), J = 1, JN2)
                  END DO
                  IF (IAND(JF2, 4) == 4) THEN
                     DO K = 1, KH2
                        WRITE (JOP, 92) (DH2(J, K), J = 1, JN2)
                     END DO
                  END IF
               END IF
               CALL FLUSH (JOP)
               REWIND (JOP)
               KNT = 0
               DO
                  JO1 = -1
                  JO2 = -1
                  JIP = JOP
                  IF (MUTATE) THEN
                     ! Do mutation run of the code
                     CALL MS2BSS ( JOP, JO1 )
                     IF (JO1 /= -1) EXIT
! exit the code after constructing a starting model if that's all we want
                     IF (START_MODEL_ONLY) THEN
                        JO1 = 0
                        EXIT
                     END IF
                  END IF
c Do alternately KP steps of *1, then sufficient steps of *2 for *2 to
c catch up, then back to *1 ..., until either KPT steps are taken, or 
c one star terminates. The `termination code' JO is possibly altered 
c by various subroutines: see above.
                  CALL STAR12 ( JO1, JC1, JOP, JIP, KSV, 22 )
                  REWIND (22)
c He flash (JO=8): replace last model by ZAHB model, then continue
                  IF ( JO1.EQ.8 ) THEN
                     JO3 = -1
                     CALL FGB2HB ( JOP, 1, JO3 )
                     IF ( JO3.EQ.13 ) CYCLE
                     JO1 = 16
                  END IF
                  JM1 = JMOD
                  REWIND (1)
                  CALL PRUNER ( 1, 3, ISB )
                  CLOSE (3)
                  CLOSE (1)
                  OPEN(UNIT = 1, FILE=TRIM(FNAME)//SUFFIX(1),STATUS='OLD',POSITION='APPEND')
                  JO2 = -1
                  JMOD = 0
! Cycle for *2, but only if doing a binary in normal (ie, non-TWIN) mode
                  IF ( ISB == 2 .AND. KTW == 1 ) THEN
                     JIP = JOP
                     CALL STAR12 ( JO2, JC2, JOP, JIP, KSV, 22 )
                     CALL FLUSH (3)
                     CALL FLUSH (JOP)
                     CALL FLUSH (27 - JOP)
                     REWIND (3)
                     REWIND (JOP)
                     REWIND (27 - JOP)
                     REWIND (22)
                     KNT = KNT + KP
c Loop back if *2 has caught up with *1, etc 
                     IF ( KNT.LT.KPT .AND. JO1.EQ.0 .AND. JO2.EQ.3 ) THEN
                        JOP = 27 - JOP
                        CYCLE
                     END IF
                  END IF
                  EXIT
               END DO
c Write final answer on unit 9; initial values and termination code on unit 8
               DO I=1, ISB
                  REWIND (I)
                  CALL PRUNER ( I, 9, ISB )
                  CLOSE (1)
                  OPEN(UNIT = 1, FILE=TRIM(FNAME)//SUFFIX(1),STATUS='OLD',POSITION='APPEND')
               ENDDO
               CALL FLUSH ( 9 )
               REWIND (22)
               REWIND (JOP)
               REWIND (27 - JOP)
               IF ( JO1.EQ.2 ) JO1 = JO1 + 10*JC1
               IF ( JO2.EQ.2 ) JO2 = JO2 + 10*JC2
               WRITE (8,996) JO1, JO2, ML, QL, XL, JM1, JMOD, 27 - JOP
               CALL FLUSH ( 8 )
  996 FORMAT (2I5, 3F10.5, 7I5)
! Write description of the mass grid. Really only useful for ZAMS libraries though
               DMM = - DLOG10(1.0D0 - CMI*DTY*CSY)*KSV
               WRITE (29,97) ML, DMM, ML + DMM*KPT/KSV, KSV/KT4
 97   FORMAT (3F10.5, I5)
               CALL FLUSH (29)

!MvdS If the code was started for more than one model, rename the files for each model		       
!	       WRITE(FNAME2,'(3I2.2)')NML,NQL,NXL
!	       IF(KML*KQL*KXL.GT.1) THEN
!	       DO I=1,NOUT(ISB)
!	         CLOSE(UNITS(I))
!	         COMMAND='mv -f '//FNAME(1:FLEN)//SUFFIX(I)//' '
!     &           //FNAME(1:FLEN)//FNAME2//SUFFIX(I)
!                 DUMMY = SYSTEM(COMMAND)
!	       ENDDO
!	       ENDIF

!MvdS If the code was started for more than one model, create a directory per model and move the files
               IF(KML*KQL*KXL.GT.1) THEN
                  COMMAND='mkdir '//FNAME(1:FLEN)//FNAME2
                  DUMMY = SYSTEM(COMMAND)
                  DO I=1,NOUT(ISB)
                  CLOSE(UNITS(I))
                  COMMAND='mv -f '//FNAME(1:FLEN)//SUFFIX(I)//' '
     &             //FNAME(1:FLEN)//FNAME2
                  DUMMY = SYSTEM(COMMAND)
                  ENDDO
               ENDIF

               XL = XL + DXL
            END DO
            QL = QL + DQL
         END DO
         ML = ML + DML
      END DO

!End of a run loop
      WRITE (6,*) '  Evolution done'
      END
      
      

! Initialise constants, read in tables of physical data that are unchanged
! during the course of one run      
      SUBROUTINE SETSUP
      USE SETTINGS
      USE CONSTANTS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /UBVDAT/ TGR(34), GGR(13), TAB(34,13,5)
      
! Initialise some constants that cannot be computed at compile time on some
! compilers
      CALL INITIALISE_CONSTANTS

! Read nuclear reaction rates and neutrino loss parameters
! Used to be read from the same file as the opacity data, which is stupid
!  since this is independent of metalicity
      CALL LOAD_REACTION_AND_NEUTRINO_RATES(42)
 
! Load opacity tables - needs to be done after reading init.dat, because we need
!  to know what tables to read.
!      CALL LOAD_OPACITY(20)

c Read Bol. Corr, U-B, B-V table. 
  991 FORMAT (3(10F10.5,/), 4F10.5,/, 10F10.5,/, 3F10.5,/, 442(5f8.3,/))
      READ (21,991) TGR, GGR, (((TAB(I,J,K), K=1,5), J=1,13), I=1,34)
      CLOSE (21)

c Read nuclear reaction (QRT) and neutrino (QNT) Q values, in MeV; constants 
c for electron screening (CZA, CZB, CZC, CZD, VZ); atomic parameters (CBN, KZN),
c with masses (CAN) consistent with Q-values; ionization potentials (CHI) and 
c statistical weights (COM); molecular hydrogen parameters (CH2)
      CALL LOAD_ATOMIC_DATA(26)

      RETURN
      END
      
      
      
      SUBROUTINE READ_INIT_RUN(IR)
      USE MESH
      USE MESH_ENC
      USE SETTINGS
      USE FILE_EXISTS_MODULE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IR
      INTEGER :: S_KPT
      INTEGER :: IOERROR
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EP(3)
      DOUBLE PRECISION :: SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      DOUBLE PRECISION :: SM, DTY, AGE, PER, BMS, ECC, P, ENC
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER :: JMOD, JB, JNN, JTER, JOC, JKH, KH, KTW, KW(260)
      INTEGER :: ISB, IP1, IM1, IP2, IM2, KPT, KP, KML, KQL, KXL, KR, JMX
      DOUBLE PRECISION :: ML1, DML, QL1, DQL, XL1, DXL, ROT, EX
      CHARACTER*500 :: STARTFILE, NAME, ZAMSFILE
      COMMON H, DH, EP, KH, KTW, KW
      COMMON /TN1/ SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      COMMON /T0 / SM, DTY, AGE, PER, BMS, ECC, P, ENC
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /SAVEINIT/ S_KPT
      COMMON /CINIT_RUN/ ISB, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, JMX, STARTFILE, ZAMSFILE
! EG: Namelist for namelist I/O of init.run (fort.23)      
      NAMELIST /INIT_RUN/ ISB, KTW, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, SM, DTY, AGE, PER, BMS, ECC, 
     & P, ENC, JMX, UC, START_MODEL_ONLY, STARTFILE, ZAMSFILE,
     & START_WITH_RIGID_ROTATION
c First try to use namelist I/O; if that fails use `normal' old-fashioned IO
      ZAMSFILE = ''
      STARTFILE = ''
      START_MODEL_ONLY = .TRUE.
      READ (IR, NML=INIT_RUN, IOSTAT=IOERROR)
      REWIND (IR)
      IF (IOERROR /= 0) THEN
         READ (IR, *) ISB, KTW, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML, 
     :      QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, SM, DTY, AGE, PER, BMS,
     :      ECC, P, ENC, JMX, UC
         REWIND (IR)
      END IF
c If a starting model was named in the input file, use that instead.
c No more fort.nnn symbolic links, yay!      
      IF ( LEN(TRIM(STARTFILE))>0 ) THEN
         IF (.NOT. FILE_EXISTS(STARTFILE)) THEN
            INQUIRE(UNIT=IR, NAME=NAME)
            WRITE (0, *) 'Start file "', TRIM(STARTFILE), '" not found. Check setting in ', TRIM(NAME), '.'
            STOP
         ELSE
            IP1 = 62
            OPEN(UNIT = IP1, ACTION="READ", FILE=STARTFILE)
         END IF
      END IF
! EG: store some information we may later need so we don't have to reread the file
      S_KPT = KPT
! EG: Do some magic for ISB == -1, -2 meaning we want to do a mutation rather than
!  a regular evolution run. A mutation run is a run where we start with a
!  normal star and end up with something weird that we want to evolve next.
!  Useful for getting merger remnants, for instance, to evolve
      IF (ISB == -1 .OR. ISB == -2) THEN
         ! Make sure all parts of the code that need to know know we're mutating
         MUTATE = .TRUE.
         ! Keep the type of mutation we want to do         
         MUTATION_MODE = ISB
         ! We really do want the code to calculate things for one star
         ISB = 1
         ! Read in target model global data
         TPROF_IN = IP1
         READ (TPROF_IN, *) SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
         REWIND (TPROF_IN)
         ! To start our mutation run, we first need a suitable model from the
         !  main sequence library, so we have to setup some variables to the
         !  mass that we want.
         SM = SM0
         ML1 = LOG10(SM0)
         IP1 = 16
      END IF
      END SUBROUTINE



! Load model number IM from the file IP
      SUBROUTINE LOAD_STAR_MODEL(IP,IM, H, DH, SM,DTY,AGE,PER,BMS,ECC,P,ENC,KH,JN,JF)
      USE MESH
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP, IM
      INTEGER, INTENT(OUT) :: KH,JN,JF
      DOUBLE PRECISION, INTENT(OUT) :: SM,DTY,AGE,PER,BMS,ECC,P,ENC
      DOUBLE PRECISION, INTENT(OUT) :: H(NVAR,NM), DH(NVAR,NM)
      INTEGER :: I, J, K, JM

      DO I = 1, IM
c Skip all models upto IM
         READ (IP, *) SM, DTY, AGE, PER, BMS, ECC, P, ENC, KH, K, JM, K, JN, JF
         DO K = 1, KH
            READ (IP, *) (H(J, K), J = 1, JN)
         END DO
         IF (IAND(JF, 4) == 4) THEN
            DO K = 1, KH
               READ (IP, *) (DH(J, K), J = 1, JN)
            END DO
         END IF
      END DO
      IF (.NOT. (IAND(JF, 4) == 4)) THEN
        DH(:,:) = 0.0D0
      END IF
      REWIND (IP) 
      END SUBROUTINE
      
      
      
      SUBROUTINE PRUNER ( JN, JT, ISB )
c Reads in a long output file (eg fort.1 or 2) and prunes it down
c to 5 lines per timestep
      IMPLICIT REAL*8 (B-H, L-Z)
      DIMENSION Y(10000), Z(20000,85), IM(20000), IB(20000)
      CHARACTER*1 CHAR(150)
  100 FORMAT (150A1)
   98 FORMAT (3(1P, E16.9, 5E10.3, 0P, F8.3, 7F8.5, 2F8.4, F9.5, /), 
     :         (1P, E16.9, 5E10.3, 0P, F8.3, 7F8.3, 2F8.4, F9.5, /), 
     :          1P, E16.9, 5E10.3, 0P, 8F8.3, 24X, I6, I2)
      K=1
      DO
         READ (JN, 100, END=7) CHAR
         IF ( CHAR(8).NE.'M' ) CYCLE
         READ (JN, 100) CHAR
         READ (JN, 98, END=7) (Z(K, J), J = 1, 82), IM(K), IB(K)
         K = K + 1
      END DO
C  1-M*;    2-Porb; 3-xi ;  4-tn;   5-LH;    6-P(cr); 7-McHe;  8-CXH;   9-CXHe;
C  10-CXC; 11-CXN; 12-CXO; 13-CXNe; 14-CXMg; 15-Cpsi; 16-Crho; 17-CT; 
C  
C 18-dty;  19-Prot; 20-zet; 21-tKh; 22-LHe;  23-RCZ;  24-McCO; 25-TXH; 26-TXHe;
C 27-TXC;  28-TXN;  29-TXO; 30-TXNe; 31-TXMg; 32-Tpsi; 33-Trho; 34-TT; 
C
C 35-age;  36-ecc;  37-mdt; 38-tET; 39-LCO;  40-DRCZ; 41-McNe; 42-SXH; 43-SXHe;
C 44-SXC;  45-SXN;  46-SXO; 47-SXNe; 48-SXMg; 49-Spsi; 50-Srho; 51-ST;
C 
C 52-cM*;  53-RLF1; 54-RLF2; 55-DLT;  56-Lnu; 57-RA/R; 58-MH; 
C 59 to 66-conv. bdries;     67-logR; 68-logL
C
C 69-Horb; 70-F1; 71-DF21; 72-BE;  73-Lth;  74-Bp;  75-MHe;  76 to 81-semiconv.
C bdries;  82-k**2; 83-MV; 84- B-V; 85- U-B; then JMOD
 7    K = K - 1
      Y(1:K) = Z(1:K, 35)
c age = Y(J). Look for *small* backward jumps in age: due to convergence 
c failure and backtracking with timestep reduced.
      DO J = 1, K - 1
         IF ( DABS(Y(J)).LT.1.0D-8 ) Y(J) = -1.0D0
         IF ( DABS(Y(J + 1)).LT.1.0D-8 ) Y(J + 1) = -1.0D0
         IF ( Y(J).GT.Y(J + 1) ) IM(J) = -1
      END DO

   99 FORMAT (3(1P, E16.9, 5E10.3, 0P, F8.3, 7F8.5, 2F8.4, F9.5, /), 
     :         (1P, E16.9, 5E10.3, 0P, F8.3, 7F8.3, 2F8.4, F9.5, /), 
     :          1P, E16.9, 5E10.3, 0P, 11F8.3, I6, I2)
      DO J = 1, K
         IF ( J /= K ) THEN
            IF ( Y(J).LT.0.0D0 .AND. Y(J + 1).LT.0.0D0) CYCLE
         END IF
C add In Bol. Corrn, colours./........./........./......../........./..
         CALL LT2UBV (Z(J,68), Z(J,51), Z(J,1), Z(J,83), Z(J,84), Z(J,85))
         IF ( IM(J).GT.0 .OR. Y(J).LT.0.0D0 ) WRITE (JT, 99) 
     :    (Z(J, I), I = 1, 34), Y(J), (Z(J, I), I = 36, 85), IM(J), IB(K)
      END DO
      IF ( (JN.EQ.1 .AND. ISB.EQ.2) .OR. JT.NE.9 ) RETURN
      Z(1, 1:85) = 0.0D0
      WRITE (9,99) (Z(1, I), I = 1, 85), 0, ISB
      RETURN
      END
      
      
      SUBROUTINE LT2UBV ( LOGL, LOGT, MASS, MV, BMINV, UMINB )
c yields values of MV, B-V and U-B for given log L, log T, and mass
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /UBVDAT/ TGR(34), GGR(13), TAB(34,13,5)
      DIMENSION BC(4), UB(4), BV(4), VR(4), RI(4)
      LOGM = LOG10(MASS)
      LOGG = LOGM + 4.0D0*LOGT - LOGL - 10.6071D0
c determine values of log g to interpolate between
      ING1 = 1
      ING2 = 13
 1    IF ( ING2 - ING1.GT.1 ) THEN
         I = (ING1 + ING2)/2
         IF ( GGR(I).GT.LOGG ) THEN
            ING2 = I
         ELSE
            ING1 = I
         END IF
         GO TO 1
      END IF
      LOGG1 = GGR(ING1)
      LOGG2 = GGR(ING2)
c determine values of log T to interpolate between
      INT1 = 1
      INT2 = 34
 2    IF ( INT2 - INT1.GT.1 ) THEN
         I = (INT1 + INT2)/2
         IF ( TGR(I).GT.LOGT ) THEN
            INT2 = I
         ELSE
            INT1 = I
         END IF
         GO TO 2
      END IF
      LOGT1 = TGR(INT1)
      LOGT2 = TGR(INT2)
      DO 3 K = 1, 2
         DO 3 J = 1, 2
            K0 = (K - 1)*2 + J
            K1 = INT1 - 1 + K
            K2 = ING1 - 1 + J
            BC(K0) = TAB(K1, K2, 1)
            UB(K0) = TAB(K1, K2, 2)
            BV(K0) = TAB(K1, K2, 3)
            VR(K0) = TAB(K1, K2, 4)
    3       RI(K0) = TAB(K1, K2, 5)
      DG1 = (LOGG - LOGG1)/(LOGG2 - LOGG1)
      DG2 = 1.0D0 - DG1
      BC1 = BC(2)*DG1 + BC(1)*DG2
      UB1 = UB(2)*DG1 + UB(1)*DG2
      BV1 = BV(2)*DG1 + BV(1)*DG2
      BC2 = BC(4)*DG1 + BC(3)*DG2
      UB2 = UB(4)*DG1 + UB(3)*DG2
      BV2 = BV(4)*DG1 + BV(3)*DG2
      DT1 = (LOGT - LOGT1)/(LOGT2 - LOGT1)
      DT2 = 1.0D0 - DT1
      BCX = BC2*DT1 + BC1*DT2
      UBX = UB2*DT1 + UB1*DT2
      BVX = BV2*DT1 + BV1*DT2
      MBOL = 4.75D0 - 2.5D0*LOGL
      MV = MBOL - BCX
      BMINV = BVX
      UMINB = UBX
      RETURN
      END
