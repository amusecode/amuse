      SUBROUTINE CLOUD(J)
*
*
*       Initialization of interstellar cloud.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),CLDOT,VCL,SIGMA,RB2,PCL2,
     &                TCL,STEPCL,NCL,NEWCL
      REAL*8  RAN2,A(13)
*
*
*       Increase new cloud counter.
      NEWCL = NEWCL + 1
*       Generate a random position on the boundary sphere.
    1 A(1) = 0.0
      DO 2 K = 1,3
*       Note that IDUM1 is saved in COMMON6 for restarts.
          A(K+1) = 2.0*RAN2(IDUM1) - 1.0
          A(1) = A(1) + A(K+1)**2
    2 CONTINUE
      IF (A(1).GT.1.0)  GO TO 1
      DO 3 K = 1,3
          XCL(K,J) = A(K+1)*SQRT(RB2/A(1))
    3 CONTINUE
      RANPHI = TWOPI*RAN2(IDUM1)
      RANDI = SQRT(RAN2(IDUM1))
      RANDI = DASIN(RANDI)
*       Azimuth & inclination for velocity from Henon (A & A, 19, 488).
      A(1) = SIN(RANPHI)
      A(2) = COS(RANPHI)
      A(3) = SIN(RANDI)
      A(4) = COS(RANDI)
      IF (KZ(13).GE.2) THEN
          DO 10 K = 1,3
              A(5) = RAN2(IDUM1)
              A(6) = TWOPI*RAN2(IDUM1)
              A(K+10) = VCL + SIGMA*SQRT(-2.0D0*LOG(A(5)))*COS(A(6))
   10     CONTINUE
          A(10) = SQRT(A(11)**2 + A(12)**2 + A(13)**2)
      ELSE
          A(10) = VCL
      END IF
      A(11) = A(10)*A(2)*A(3)
      A(12) = A(10)*A(1)*A(3)
      A(13) = -A(10)*A(4)
      PHI = DATAN2 (-XCL(1,J),(-XCL(2,J)))
      A(1) = XCL(3,J)/SQRT(RB2)
      THETA = DACOS(A(1))
*       Eulerian angles of the transformation.
      A(1) = SIN(PHI)
      A(2) = COS(PHI)
      A(3) = SIN(THETA)
      A(4) = COS(THETA)
*       Generate isotropic velocities ignoring c.m. motion of cluster.
      XDOTCL(1,J) = -A(4)*A(1)*A(11) + A(2)*A(12) - A(3)*A(1)*A(13)
      XDOTCL(2,J) = -A(4)*A(2)*A(11) - A(1)*A(12) - A(3)*A(2)*A(13)
      XDOTCL(3,J) = -A(3)*A(11) + A(4)*A(13)
*       Convert to rotating coordinates.
      XDOTCL(1,J) = XDOTCL(1,J) + 0.5*TIDAL(4)*XCL(2,J)
      XDOTCL(2,J) = XDOTCL(2,J) - 0.5*TIDAL(4)*XCL(1,J)
      BODYCL(J) = 0.0
*       The cloud mass is increased smoothly from zero at boundary.
      DO 20 K = 1,3
          XCL(K,J) = XCL(K,J) - RDENS(K)
   20 CONTINUE
*       Coordinates transformed from density centre to inertial frame.
*
      RETURN
*
      END
