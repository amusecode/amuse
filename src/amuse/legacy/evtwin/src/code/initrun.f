      module init_run
      
      contains
      
      subroutine load_basic_init_run
      USE MESH
      USE MESH_ENC
      USE FILE_EXISTS_MODULE
      IMPLICIT NONE
      INTEGER :: S_KPT
      INTEGER :: IOERROR
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EP(3)
      DOUBLE PRECISION :: SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      DOUBLE PRECISION :: SM, DTY, AGE, PER, BMS, ECC, P, ENC
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER :: JMOD, JB, JNN, JTER, JOC, JKH, KH, KTW, KW(260)
      INTEGER :: ISB, IP1, IM1, IP2, IM2, KPT, KP, KML, KQL, KXL, KR, JMX
      DOUBLE PRECISION :: ML1, DML, QL1, DQL, XL1, DXL, ROT, EX
      CHARACTER*500 :: STARTFILE, NAME
      COMMON H, DH, EP, KH, KTW, KW
      COMMON /TN1/ SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      COMMON /T0 / SM, DTY, AGE, PER, BMS, ECC, P, ENC
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /SAVEINIT/ S_KPT
      COMMON /CINIT_RUN/ ISB, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, JMX, STARTFILE

      ISB     =           1;  KTW     =           1
      IP1     =          16;  IM1     =           1
      IP2     =          16;  IM2     =           1
      KPT     =        2000;  KP      =         200
      STARTFILE = ''

      ML1     =  0.00E+00;  DML     =  0.30E+00;  KML     =  1
      QL1     =  5.00E-02;  DQL     =  5.00E-02;  KQL     =  1
      XL1     =  6.00E+00;  DXL     =  0.30E+00;  KXL     =  1

      ROT     =  1.0;          KR      =  0
      EX      =  0.0
      SM      = -1.0
      DTY     = -1.0
      AGE     =  0.0
      PER     = -1.0
      BMS     =  1.0e3
      ECC     =  0.0
      P       = -1.0
      ENC     =  0.0
      JMX     =  1

      UC      = (/   
     &1.00E-01, 2.00E+12, 1.00E+02, 0.00E+00, 3.00E+00, 5.30E+00, 1.20E+00,
     &6.30E+00, 3.00E+02, 0.00E+00, 1.00E-06, 1.00E+06, 1.00E+03, 1.00E+03,
     &0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00 /)

      END SUBROUTINE

      end module

