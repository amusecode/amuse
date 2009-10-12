      module init_dat
      
      type init_dat_settings
      
         DOUBLE PRECISION :: EPS, DEL, DH0, WANTED_EPS

         DOUBLE PRECISION :: CDC_MS, CDC_HEC, CDC_HES, CDC_DBLSH, CDC5

         INTEGER :: KE1, KE2, KE3, KBC, KEV, KFN, KL, JH(3)
         INTEGER :: KP_VAR(40), KP_EQN(40), KP_BC(40)

         INTEGER :: KH2, KR1, KR2, JCH, KTH, KX, KY, KZ
         INTEGER :: KCL, KION, KAM, KOP, KCC, KNUC, KCN
         INTEGER :: KT1, KT2, KT3, KT4, KT5, KSV
     
         INTEGER :: KN, KJN(40)

         DOUBLE PRECISION :: CT1, CT2, CT3
         DOUBLE PRECISION :: CT(10)

         DOUBLE PRECISION :: CH, CC, CN, CO, CNE, CMG, CSI, CFE
         DOUBLE PRECISION :: CALP, CU, COS, CPS, CRD, CTH, CGRS, CCAC
         DOUBLE PRECISION :: CDSI, CSHI, CSSI, CESC, CGSF, CFMU, CFC
         DOUBLE PRECISION :: ARTMIX

         DOUBLE PRECISION :: CXB, CGR, CEA, CET
         DOUBLE PRECISION :: CMT, CMS, CML, CHL, CTF, CLT
         DOUBLE PRECISION :: CMI, CMR, CMJ, CMV, CMK, CMNL
         DOUBLE PRECISION :: CPA, CBR, CSU, CSD, CDF, CGW, CSO, CMB
         DOUBLE PRECISION :: CPHOTONTIRE
         
         INTEGER :: CONVECTION_SCHEME
         LOGICAL :: USE_FUDGE_CONTROL, ALLOW_EXTENSION, ALLOW_UNDERRELAXATION
         LOGICAL :: ALLOW_OVERRELAXATION, ALLOW_MDOTRELAXATION
         LOGICAL :: ALLOW_EGENRELAXATION, USE_PREVIOUS_MU
         LOGICAL :: ALLOW_AVMURELAXATION, USE_QUADRATIC_PREDICTIONS

         DOUBLE PRECISION :: OFF_CENTRE_WEIGHT
         
         INTEGER :: KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2(3)
         INTEGER :: KP_VAR_2(40), KP_EQN_2(40), KP_BC_2(40)

      end type

      contains
      
      subroutine push_init_dat(init_dat, KH2, KR1, KR2, KSV, KT5, JCH)
      USE MESH
      USE CONSTANTS
      USE SETTINGS
      USE FUDGE_CONTROL
      USE EXTRAPOLATE_DH
      IMPLICIT NONE
      TYPE(INIT_DAT_SETTINGS), INTENT(OUT) :: init_dat
      INTEGER, INTENT(IN) :: KH2, KR1, KR2, KSV, KT5, JCH
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      DOUBLE PRECISION :: CH_OPAC
      INTEGER :: KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH(3),KP_VAR(40), KP_EQN(40), KP_BC(40), 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2(3),KP_VAR_2(40), KP_EQN_2(40), KP_BC_2(40)
      DOUBLE PRECISION :: CDC_MS, CDC_HEC, CDC_HES, CDC_DBLSH, CDC5, 
     &      UNUSED1, UNUSED(17)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2

      init_dat%EPS = EPS
      init_dat%DEL = DEL
      init_dat%DH0 = DH0
      init_dat%WANTED_EPS = WANTED_EPS
      init_dat%CDC_MS = CDC_MS
      init_dat%CDC_HEC = CDC_HEC
      init_dat%CDC_HES = CDC_HES
      init_dat%CDC_DBLSH = CDC_DBLSH
      init_dat%CDC5 = CDC5
      init_dat%KE1 = KE1
      init_dat%KE2 = KE2
      init_dat%KE3 = KE3
      init_dat%KBC = KBC
      init_dat%KEV = KEV
      init_dat%KFN = KFN
      init_dat%KL = KL
      init_dat%JH(1:3) = JH(1:3)
      init_dat%KP_VAR(1:40) = KP_VAR(1:40)
      init_dat%KP_EQN(1:40) = KP_EQN(1:40)
      init_dat%KP_BC(1:40) = KP_BC(1:40)
      init_dat%KH2 = KH2
      init_dat%KR1 = KR1
      init_dat%KR2 = KR2
      init_dat%JCH = JCH
      init_dat%KTH = KTH
      init_dat%KX = KX
      init_dat%KY = KY
      init_dat%KZ = KZ
      init_dat%KCL = KCL
      init_dat%KION = KION
      init_dat%KAM = KAM
      init_dat%KOP = KOP
      init_dat%KCC = KCC
      init_dat%KNUC = KNUC
      init_dat%KCN = KCN
      init_dat%KT1 = KT1
      init_dat%KT2 = KT2
      init_dat%KT3 = KT3
      init_dat%KT4 = KT4
      init_dat%KT5 = KT5
      init_dat%KSV = KSV
      init_dat%KN = KN
      init_dat%KJN(1:40) = KJN(1:40)
      init_dat%CT1 = CT1
      init_dat%CT2 = CT2
      init_dat%CT3 = CT3
      init_dat%CT(1:10) = CT(1:10)
      init_dat%CH = CH
      init_dat%CC = CC
      init_dat%CN = CN
      init_dat%CO = CO
      init_dat%CNE = CNE
      init_dat%CMG = CMG
      init_dat%CSI = CSI
      init_dat%CFE = CFE
      init_dat%CALP = CALP
      init_dat%CU = CU
      init_dat%COS = COS
      init_dat%CPS = CPS
      init_dat%CRD = CRD
      init_dat%CTH = CTH
      init_dat%CGRS = CGRS
      init_dat%CCAC = CCAC
      init_dat%CDSI = CDSI
      init_dat%CSHI = CSHI
      init_dat%CSSI = CSSI
      init_dat%CESC = CESC
      init_dat%CGSF = CGSF
      init_dat%CFMU = CFMU
      init_dat%CFC = CFC
      init_dat%ARTMIX = ARTMIX
      init_dat%CXB = CXB
      init_dat%CGR = CGR
      init_dat%CEA = CEA
      init_dat%CET = CET
      init_dat%CMT = CMT
      init_dat%CMS = CMS
      init_dat%CML = CML
      init_dat%CHL = CHL
      init_dat%CTF = CTF
      init_dat%CLT = CLT
      init_dat%CMI = CMI
      init_dat%CMR = CMR
      init_dat%CMJ = CMJ
      init_dat%CMV = CMV
      init_dat%CMK = CMK
      init_dat%CMNL = CMNL
      init_dat%CPA = CPA
      init_dat%CBR = CBR
      init_dat%CSU = CSU
      init_dat%CSD = CSD
      init_dat%CDF = CDF
      init_dat%CGW = CGW
      init_dat%CSO = CSO
      init_dat%CMB = CMB
      init_dat%CPHOTONTIRE = CPHOTONTIRE
      init_dat%CONVECTION_SCHEME = CONVECTION_SCHEME
      init_dat%USE_FUDGE_CONTROL = USE_FUDGE_CONTROL
      init_dat%ALLOW_EXTENSION = ALLOW_EXTENSION
      init_dat%ALLOW_UNDERRELAXATION = ALLOW_UNDERRELAXATION
      init_dat%ALLOW_OVERRELAXATION = ALLOW_OVERRELAXATION
      init_dat%ALLOW_MDOTRELAXATION = ALLOW_MDOTRELAXATION
      init_dat%ALLOW_EGENRELAXATION = ALLOW_EGENRELAXATION
      init_dat%USE_PREVIOUS_MU = USE_PREVIOUS_MU
      init_dat%ALLOW_AVMURELAXATION = ALLOW_AVMURELAXATION
      init_dat%USE_QUADRATIC_PREDICTIONS = USE_QUADRATIC_PREDICTIONS
      init_dat%OFF_CENTRE_WEIGHT = OFF_CENTRE_WEIGHT
      init_dat%KE1_2 = KE1_2
      init_dat%KE2_2 = KE2_2
      init_dat%KE3_2 = KE3_2
      init_dat%KBC_2 = KBC_2
      init_dat%KEV_2 = KEV_2
      init_dat%KFN_2 = KFN_2
      init_dat%KL_2 = KL_2
      init_dat%JH_2(1:3) = JH_2(1:3)
      init_dat%KP_VAR_2(1:40) = KP_VAR_2(1:40)
      init_dat%KP_EQN_2(1:40) = KP_EQN_2(1:40)
      init_dat%KP_BC_2(1:40) = KP_BC_2(1:40)
      end subroutine

      subroutine pop_init_dat(init_dat, KH2, KR1, KR2, KSV, KT5, JCH)
      USE MESH
      USE CONSTANTS
      USE SETTINGS
      USE FUDGE_CONTROL
      USE EXTRAPOLATE_DH
      IMPLICIT NONE
      TYPE(INIT_DAT_SETTINGS), INTENT(IN) :: init_dat
      INTEGER, INTENT(OUT) :: KH2, KR1, KR2, KSV, KT5, JCH
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      DOUBLE PRECISION :: CH_OPAC
      INTEGER :: KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH(3),KP_VAR(40), KP_EQN(40), KP_BC(40), 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2(3),KP_VAR_2(40), KP_EQN_2(40), KP_BC_2(40)
      DOUBLE PRECISION :: CDC_MS, CDC_HEC, CDC_HES, CDC_DBLSH, CDC5, 
     &      UNUSED1, UNUSED(17)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2

      EPS = init_dat%EPS
      DEL = init_dat%DEL
      DH0 = init_dat%DH0
      WANTED_EPS = init_dat%WANTED_EPS
      CDC_MS = init_dat%CDC_MS
      CDC_HEC = init_dat%CDC_HEC
      CDC_HES = init_dat%CDC_HES
      CDC_DBLSH = init_dat%CDC_DBLSH
      CDC5 = init_dat%CDC5
      KE1 = init_dat%KE1
      KE2 = init_dat%KE2
      KE3 = init_dat%KE3
      KBC = init_dat%KBC
      KEV = init_dat%KEV
      KFN = init_dat%KFN
      KL = init_dat%KL
      JH(1:3) = init_dat%JH(1:3)
      KP_VAR(1:40) = init_dat%KP_VAR(1:40)
      KP_EQN(1:40) = init_dat%KP_EQN(1:40)
      KP_BC(1:40) = init_dat%KP_BC(1:40)
      KH2 = init_dat%KH2
      KR1 = init_dat%KR1
      KR2 = init_dat%KR2
      JCH = init_dat%JCH
      KTH = init_dat%KTH
      KX = init_dat%KX
      KY = init_dat%KY
      KZ = init_dat%KZ
      KCL = init_dat%KCL
      KION = init_dat%KION
      KAM = init_dat%KAM
      KOP = init_dat%KOP
      KCC = init_dat%KCC
      KNUC = init_dat%KNUC
      KCN = init_dat%KCN
      KT1 = init_dat%KT1
      KT2 = init_dat%KT2
      KT3 = init_dat%KT3
      KT4 = init_dat%KT4
      KT5 = init_dat%KT5
      KSV = init_dat%KSV
      KN = init_dat%KN
      KJN(1:40) = init_dat%KJN(1:40)
      CT1 = init_dat%CT1
      CT2 = init_dat%CT2
      CT3 = init_dat%CT3
      CT(1:10) = init_dat%CT(1:10)
      CH = init_dat%CH
      CC = init_dat%CC
      CN = init_dat%CN
      CO = init_dat%CO
      CNE = init_dat%CNE
      CMG = init_dat%CMG
      CSI = init_dat%CSI
      CFE = init_dat%CFE
      CALP = init_dat%CALP
      CU = init_dat%CU
      COS = init_dat%COS
      CPS = init_dat%CPS
      CRD = init_dat%CRD
      CTH = init_dat%CTH
      CGRS = init_dat%CGRS
      CCAC = init_dat%CCAC
      CDSI = init_dat%CDSI
      CSHI = init_dat%CSHI
      CSSI = init_dat%CSSI
      CESC = init_dat%CESC
      CGSF = init_dat%CGSF
      CFMU = init_dat%CFMU
      CFC = init_dat%CFC
      ARTMIX = init_dat%ARTMIX
      CXB = init_dat%CXB
      CGR = init_dat%CGR
      CEA = init_dat%CEA
      CET = init_dat%CET
      CMT = init_dat%CMT
      CMS = init_dat%CMS
      CML = init_dat%CML
      CHL = init_dat%CHL
      CTF = init_dat%CTF
      CLT = init_dat%CLT
      CMI = init_dat%CMI
      CMR = init_dat%CMR
      CMJ = init_dat%CMJ
      CMV = init_dat%CMV
      CMK = init_dat%CMK
      CMNL = init_dat%CMNL
      CPA = init_dat%CPA
      CBR = init_dat%CBR
      CSU = init_dat%CSU
      CSD = init_dat%CSD
      CDF = init_dat%CDF
      CGW = init_dat%CGW
      CSO = init_dat%CSO
      CMB = init_dat%CMB
      CPHOTONTIRE = init_dat%CPHOTONTIRE
      CONVECTION_SCHEME = init_dat%CONVECTION_SCHEME
      USE_FUDGE_CONTROL = init_dat%USE_FUDGE_CONTROL
      ALLOW_EXTENSION = init_dat%ALLOW_EXTENSION
      ALLOW_UNDERRELAXATION = init_dat%ALLOW_UNDERRELAXATION
      ALLOW_OVERRELAXATION = init_dat%ALLOW_OVERRELAXATION
      ALLOW_MDOTRELAXATION = init_dat%ALLOW_MDOTRELAXATION
      ALLOW_EGENRELAXATION = init_dat%ALLOW_EGENRELAXATION
      USE_PREVIOUS_MU = init_dat%USE_PREVIOUS_MU
      ALLOW_AVMURELAXATION = init_dat%ALLOW_AVMURELAXATION
      USE_QUADRATIC_PREDICTIONS = init_dat%USE_QUADRATIC_PREDICTIONS
      OFF_CENTRE_WEIGHT = init_dat%OFF_CENTRE_WEIGHT
      KE1_2 = init_dat%KE1_2
      KE2_2 = init_dat%KE2_2
      KE3_2 = init_dat%KE3_2
      KBC_2 = init_dat%KBC_2
      KEV_2 = init_dat%KEV_2
      KFN_2 = init_dat%KFN_2
      KL_2 = init_dat%KL_2
      JH_2(1:3) = init_dat%JH_2(1:3)
      KP_VAR_2(1:40) = init_dat%KP_VAR_2(1:40)
      KP_EQN_2(1:40) = init_dat%KP_EQN_2(1:40)
      KP_BC_2(1:40) = init_dat%KP_BC_2(1:40)
      end subroutine

      subroutine load_basic_init_dat(KH2, KR1, KR2, KSV, KT5, JCH)
      USE MESH
      USE CONSTANTS
      USE SETTINGS
      USE FUDGE_CONTROL
      USE EXTRAPOLATE_DH
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: KH2, KR1, KR2, KSV, KT5, JCH
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      DOUBLE PRECISION :: CH_OPAC
      INTEGER :: KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH(3),KP_VAR(40), KP_EQN(40), KP_BC(40), 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2(3),KP_VAR_2(40), KP_EQN_2(40), KP_BC_2(40)
      DOUBLE PRECISION :: CDC_MS, CDC_HEC, CDC_HES, CDC_DBLSH, CDC5, 
     &      UNUSED1, UNUSED(17)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2

      ! Load some reasonable default values
      KH2   =  KH
      JCH   =    2
      KL    =    1
      JH(1:3) =   (/ 0,0,0 /)

      KT1   =  100
      KT2   =    2    
      KT3   =    0    
      KT4   =    1 
      KT5   =    0
      KSV   = 1000

      KTH   = 1 
      KX    = 1 
      KY    = 1 
      KZ    = 1

      KCL   = 1
      KION  = 5
      KAM   = 1
      !KOP   = 4
      KCC   = 0
      KNUC  = 0
      KCN   = 0

      KR1   = 20
      KR2   = 5

      EPS   = 1.00E-006
      DEL   = 1.00E-002
      DH0   = 1.00E-007

      USE_QUADRATIC_PREDICTIONS = .FALSE.

      WANTED_EPS             = 1.00E-008
      USE_FUDGE_CONTROL      = .TRUE.
      ALLOW_EXTENSION        = .FALSE.
      ALLOW_UNDERRELAXATION  = .FALSE.
      ALLOW_OVERRELAXATION   = .FALSE.
      ALLOW_EGENRELAXATION   = .FALSE.
      ALLOW_MDOTRELAXATION   = .FALSE.
      ALLOW_AVMURELAXATION   = .FALSE.
      USE_PREVIOUS_MU        = .TRUE.
      OFF_CENTRE_WEIGHT      = 1.0D16

      CDC_MS    = 0.01
      CDC_HEC   = 0.25
      CDC_HES   = 1.00
      CDC_DBLSH = 4.00
      CDC5      = 1.00

      CT1   =  0.8
      CT2   =  1.2
      CT3   =  0.3

      KE1   =  6
      KE2   =  5
      KE3   =  0
      KBC   =  3
      KEV   =  1
      KFN = NFUNC    ! Just copy all functions; common error to forget. 

      KP_VAR(1:12)  =  (/ 7, 8, 4, 5, 3, 9,10,11,16, 1, 2, 6 /)
      KP_EQN(1:11)  =  (/ 1, 2, 3, 4, 5,13, 6, 7, 8, 9,10    /)
      KP_BC(1:12)   =  (/ 1, 2, 3, 4, 5,13, 6, 7, 8, 6, 7, 8 /)

      CC    =  0.176E+000
      CN    =  5.200E-002
      CO    =  0.502E+000
      CNE   =  9.200E-002
      CMG   =  3.400E-002
      CSI   =  7.200E-002
      CFE   =  7.200E-002

      CT    =  (/ 0.0E+00, 0.0E+00, 5.0E-02, 5.0E-02, 0.15E+00,
     &            2.0E-02, 0.45E+0, 1.0E-04, 1.0E+15, 2.0E+04 /)

      CONVECTION_SCHEME = 1
      CALP  =  2.000E+000
      CU    =  0.100E+000
      COS   =  0.120E+000
      CPS   =  0.120E+000
      CRD   =  1.000E-002
      CXB   =  0.150E+000
      CGR   =  1.000E-003
      CTH   =  1.000E+000
      CGRS  =  0.000E+000
      CCAC  =  0.000E+000

      ARTMIX = 0.0

      CDSI  =  0.000E+000
      CSHI  =  0.000E+000
      CSSI  =  0.000E+000
      CESC  =  0.000E+000
      CGSF  =  0.000E+000

      CEA   =  1.0E+02
      CET   =  1.0E-08

      CMI_MODE = 1
      ZSCALING_MDOT = 0.8

      CMT   =  0.0E0
      CMS   =  0.0E0
      CMI   =  0.0E0
      CMR   =  0.0E0
      CMJ   =  0.0E0
      CMV   =  0.0E0
      CMK   =  0.0E0
      CMNL  =  0.0E0
      CML   =  0.0E0
      CHL   =  0.0E0
      CTF   =  0.0E0
      CLT   =  0.0E0

      CPA   =  0.0E0
      CBR   =  0.0E0
      CSU   =  0.0E0
      CSD   =  0.0E0
      CDF   =  1.0E-02
      CGW   =  1.0E0
      CSO   =  1.0E0
      CMB   =  0.0E0

      KN      = 12
      KJN(1:12) =  (/ 1,  2,  3,  5,  6,  7, 25, 26, 27, 29, 30, 31 /)

      ! Export variables to common block
      CDC(1) = CDC_MS 
      CDC(2) = CDC_HEC
      CDC(3) = CDC_HES
      CDC(4) = CDC_DBLSH
      CDC(5) = CDC5

      end subroutine


      end module

