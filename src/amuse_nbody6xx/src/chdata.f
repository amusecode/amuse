      SUBROUTINE CHDATA(XJ,VJ,BODYJ)
*
*
*       Basic data for chain members.
*       -----------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      REAL*4  XJ(3,6),VJ(3,6),BODYJ(6)
*
*
*       Save global address of all chain members in common array JLIST.
      DO 2 L = 1,NCH
          DO 1 J = IFIRST,N
*       Identify each member sequentially from NAMEC (note NAME(ICH) = 0).
              IF (NAME(J).EQ.NAMEC(L).OR.NAME(J).EQ.0) THEN
                  JLIST(L) = J
              END IF
    1     CONTINUE
    2 CONTINUE
      
*       Obtain global coordinates & velocities from current chain & c.m.
      LK = 0
      DO 10 I = 1,NCH
          BODYJ(I) = BODYC(I)
          DO 5 K = 1,3
              LK = LK + 1
              XJ(K,I) = XCH(LK) + X(K,ICH)
              VJ(K,I) = VCH(LK) + XDOT(K,ICH)
    5     CONTINUE
   10 CONTINUE
*
      RETURN
*
      END
