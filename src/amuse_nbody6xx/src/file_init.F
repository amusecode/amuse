      SUBROUTINE FILE_INIT(ISTART)
*
*
*       Opening of Files with proper names.
*       ------------------------------------
*       ISTART=0: Open all other files after reading input data.
*       ISTART=1: Open only unit 1 for mydump.
*       ISTART=2: Open only unit 2 for mydump.
*
      INCLUDE 'common6.h'
*
*
      CHARACTER*10 FILE(100)
      CHARACTER*4 FRAG(100)
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
      FRAG(1) = 'fort'
      FRAG(2) = 'fort'
      FRAG(3) = 'conf'
      FRAG(4) = 'bdat'
      FRAG(7) = 'lagr'
      FRAG(8) = 'bdat'
      FRAG(9) = 'bdat'
*
      JFMIN = 1
      JFMAX = 9
      IF (ISTART.EQ.1) THEN
          JFMIN = 1
          JFMAX = 1
      END IF
      IF (ISTART.EQ.2) THEN
          JFMIN = 2
          JFMAX = 2
      END IF
*       Initialize the file names for read and write.
      DO 100 JF = JFMIN,JFMAX
          WRITE (FILE(JF),112) FRAG(JF),JF
 100  CONTINUE
 112  FORMAT(A4,'.',I1,T4)
      FILE(10) = 'dat.10    '
      FILE(11) = 'esc.11    '
      FILE(12) = 'hia.12    '
      FILE(13) = 'hid.13    '
      FILE(15) = 'per.15    '
      FILE(82) = 'bev.82    '
      FILE(83) = 'sev.83    '
#ifdef ENSEMBLE
      DO 101 JF = JFMIN,JFMAX
      if (isize.lt.10) then
          WRITE (FILE(JF),113) FRAG(JF),rank,JF
 113  FORMAT(A4,'000',I1,'.',I1)
      else if (isize.ge.10.and.isize.lt.100) then
          WRITE (FILE(JF),115) FRAG(JF),rank,JF
 115  FORMAT(A4,'00',I2,'.',I1)
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(JF),117) FRAG(JF),rank,JF
 117  FORMAT(A4,'0',I3,'.',I1)
      else if (isize.ge.1000) then
          WRITE (FILE(JF),119) FRAG(JF),rank,JF
 119  FORMAT(A4,I4,'.',I1)
      end if
 101  CONTINUE
*
      if (isize.lt.10) then
          WRITE (FILE(10),213) rank
          WRITE (FILE(11),313) rank
          WRITE (FILE(12),413) rank
          WRITE (FILE(13),513) rank
          WRITE (FILE(15),613) rank
          WRITE (FILE(82),713) rank
          WRITE (FILE(83),813) rank
 213  FORMAT('dat000',I1,'.10')
 313  FORMAT('esc000',I1,'.11')
 413  FORMAT('hia000',I1,'.12')
 513  FORMAT('hid000',I1,'.13')
 613  FORMAT('per000',I1,'.15')
 713  FORMAT('bev000',I1,'.13')
 813  FORMAT('sev000',I1,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(10),215) rank
          WRITE (FILE(11),315) rank
          WRITE (FILE(12),415) rank
          WRITE (FILE(13),515) rank
          WRITE (FILE(15),615) rank
          WRITE (FILE(82),715) rank
          WRITE (FILE(83),815) rank
 215  FORMAT('dat00',I2,'.10')
 315  FORMAT('esc00',I2,'.11')
 415  FORMAT('hia00',I2,'.12')
 515  FORMAT('hid00',I2,'.13')
 615  FORMAT('per00',I2,'.15')
 715  FORMAT('bev00',I2,'.13')
 815  FORMAT('sev00',I2,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(10),217) rank
          WRITE (FILE(11),317) rank
          WRITE (FILE(12),417) rank
          WRITE (FILE(13),517) rank
          WRITE (FILE(15),617) rank
          WRITE (FILE(82),717) rank
          WRITE (FILE(83),817) rank
 217  FORMAT('dat0',I3,'.10')
 317  FORMAT('esc0',I3,'.11')
 417  FORMAT('hia0',I3,'.12')
 517  FORMAT('hid0',I3,'.13')
 617  FORMAT('per0',I3,'.15')
 717  FORMAT('bev0',I3,'.13')
 817  FORMAT('sev0',I3,'.15')
      else if (isize.ge.1000) then
          WRITE (FILE(10),219) rank
          WRITE (FILE(11),319) rank
          WRITE (FILE(12),419) rank
          WRITE (FILE(13),519) rank
          WRITE (FILE(15),619) rank
          WRITE (FILE(82),719) rank
          WRITE (FILE(83),819) rank
 219  FORMAT('dat',I4,'.10')
 319  FORMAT('esc',I4,'.11')
 419  FORMAT('hia',I4,'.12')
 519  FORMAT('hid',I4,'.13')
 619  FORMAT('per',I4,'.15')
 719  FORMAT('bev',I4,'.13')
 819  FORMAT('sev',I4,'.15')
      end if
#endif
#ifdef PARALLEL
      if(rank.eq.0) then
#endif
      IF (KZ(1).GT.0.AND.ISTART.EQ.1)
     &OPEN (UNIT=1,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(1))
*
      IF (KZ(2).GT.0.AND.ISTART.EQ.2)
     &OPEN (UNIT=2,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(2))
*
      IF (ISTART.GT.0) RETURN
*
      IF (KZ(3).GT.0)
     &OPEN (UNIT=3,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(3),
     &             ACCESS='APPEND')
      IF (KZ(4).GT.0)
     &OPEN (UNIT=4,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(4),
     &             ACCESS='APPEND')
      IF (KZ(7).GE.3)
     &OPEN (UNIT=7,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(7),
     &             ACCESS='APPEND')
      IF (KZ(8).GT.0.OR.NBIN0.GT.0)
     &OPEN (UNIT=8,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(8),
     &             ACCESS='APPEND')
      IF (KZ(8).GE.2.OR.NBIN0.GT.0)
     &OPEN (UNIT=9,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(9),
     &             ACCESS='APPEND')
      IF (KZ(22).GT.0)
     &OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(10))
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)
     &OPEN (UNIT=11,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(11),
     &             ACCESS='APPEND')
      IF (KZ(11).EQ.1.OR.KZ(11).EQ.3)
     &OPEN (UNIT=12,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(12),
     &             ACCESS='APPEND')
      IF (KZ(8).GT.3)
     &OPEN (UNIT=13,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(13),
     &             ACCESS='APPEND')
      IF (KZ(44).EQ.1)
     &OPEN (UNIT=15,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(15),
     &             ACCESS='APPEND')
      IF (KZ(12).GT.0)
     &OPEN (UNIT=82,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(82),
     &             ACCESS='APPEND')
      IF (KZ(12).GT.0)
     &OPEN (UNIT=83,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(83),
     &             ACCESS='APPEND')
#ifdef PARALLEL
      end if
#endif
*
      RETURN
      END
