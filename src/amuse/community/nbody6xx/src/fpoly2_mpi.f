      SUBROUTINE FPOLY2_MPI(I1,I2,KCASE)
*     
*
*       Second & third force derivative.
*       --------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BARR/ ibarcount
      REAL*8  A(12),F2DOT(3),F3DOT(3)
      integer inum(maxpe),ista(maxpe)
*
*       Standard case, new c.m. or KS termination (KCASE = 0, 1, 2).
      JLAST = NTOT
*       Reduce loop size for new c.m. polynomial.
      IF (KCASE.EQ.1) JLAST = NTOT - 1
*
*       Include an initial skip for primordial binaries (large derivatives).
*     IF (NBIN0.GT.0.AND.TIME.EQ.0.0D0) GO TO 80
*
*       Loop over all bodies, pair #ICOMP & JCOMP, c.m. or one single body.
      call cputim(tt998)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ibarcount=ibarcount+1
      ttbar = ttbar + (tt999-tt998)*60
*
      nl = I2-I1+1
*
      inl = nl/isize
      idiff = nl - isize*inl
      irun = 0
*
      do 1003 ix = 1,isize
      inum(ix)=inl
      if(ix.le.idiff)inum(ix) = inum(ix) + 1
      ista(ix) = irun+1
 1003 irun = irun + inum(ix)
*
      istart = ista(rank+1)
      iend = ista(rank+1) + inum(rank+1) - 1
*
      do 70 i = istart,iend
*
*       Initialize the higher differences for body #I.
      DO 10 K = 1,3
          D2(K,I) = 0.0D0
          D3(K,I) = 0.0D0
          D2R(K,I) = 0.0D0
          D3R(K,I) = 0.0D0
   10 CONTINUE
*
      NNB = LIST(1,I)
*       Neglect F2DOT & F3DOT outside 5*RS unless high accuracy is needed.
      RCRIT2 = 25.0*RS(I)**2*(1.0 + 1.0/FLOAT(NNB+1))
*       Specify index of first neighbour to be identified.
      L = 2
      NAMEJ = LIST(L,I)
*
*       Sum over c.m. instead of regularized neighbours since F not known.
      DO 60 J = IFIRST,JLAST
*       Note IFIRST = 2*NPAIRS + 1, JCOMP + 1, ICOMP for KCASE = 0, 1, 2.
          DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
   15     CONTINUE
          RIJ2 = A(1)*A(1) + A(2)*A(2) + A(3)*A(3)
*
*       Ensure that all neighbours are considered.
          IF (RIJ2.GT.RCRIT2.AND.J.NE.NAMEJ) GO TO 60
*       Distant bodies do not contribute significantly to F2DOT & F3DOT.
*
          IF (KCASE.GT.0) THEN
              IF (J.GT.JCOMP.OR.KCASE.EQ.1) GO TO 30
          END IF
          IF (J.EQ.I) GO TO 60
*
*       See whether F & FDOT extrapolation is required.
          IF (KCASE.EQ.0.AND.TIME.GT.0.0D0) THEN
              IF (T0(J).LT.TIME) GO TO 30
*       Note that routine FCLOSE sets T0(J) = TIME for dominant bodies.
          END IF
*
*       Copy F & FDOT (all J at TIME = 0, otherwise dominant bodies only).
          DO 20 K = 1,3
              A(K+6) = F(K,J)
              A(K+9) = FDOT(K,J)
   20     CONTINUE
          GO TO 40
*
*       Obtain current force and first derivative to second order.
   30     DT = TIME - T0(J)
          DT1 = 0.5*DT
          DTR = TIME - T0R(J)
          DT1R = 0.5*DTR
*
          DO 35 K = 1,3
              A(K+6) = (D2R(K,J)*DT1R + D1R(K,J))*DTR + FR(K,J) +
     &                              (D2(K,J)*DT1 + D1(K,J))*DT + FI(K,J)
              A(K+9) = D2R(K,J)*DTR + D1R(K,J) + D2(K,J)*DT + D1(K,J)
   35     CONTINUE
*
   40     DO 45 K = 1,3
              A(K+3) = XDOT(K,J) - XDOT(K,I)
              A(K+6) = A(K+6) - F(K,I)
              A(K+9) = A(K+9) - FDOT(K,I)
   45     CONTINUE
*
          A13 = 1.0/RIJ2
          A14 = BODY(J)*A13*SQRT(A13)
          A15 = (A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A13
          A16 = A15*A15
          A17 = 3.0*A15
          A18 = 6.0*A15
          A19 = 9.0*A15
          A20 = (A(4)*A(4) + A(5)*A(5) + A(6)*A(6) + A(1)*A(7) +
     &                                  A(2)*A(8) + A(3)*A(9))*A13 + A16
          A21 = 9.0*A20
          A20 = 3.0*A20
          A22 = (9.0*(A(4)*A(7) + A(5)*A(8) + A(6)*A(9)) +
     &                 3.0*(A(1)*A(10) + A(2)*A(11) + A(3)*A(12)))*A13 +
     &                                               A17*(A20 - 4.0*A16)
*
          DO 50 K = 1,3
              F1DOTK = A(K+3) - A17*A(K)
              F2DOT(K) = (A(K+6) - A18*F1DOTK - A20*A(K))*A14
              F3DOT(K) = (A(K+9) - A21*F1DOTK - A22*A(K))*A14 -
     &                                                      A19*F2DOT(K)
   50     CONTINUE
*
*       See whether body #J is a neighbour of body #I.
          IF (J.NE.NAMEJ) THEN
              DO 52 K = 1,3
                  D2R(K,I) = D2R(K,I) + F2DOT(K)
                  D3R(K,I) = D3R(K,I) + F3DOT(K)
   52         CONTINUE
          ELSE
              DO 55 K = 1,3
                  D2(K,I) = D2(K,I) + F2DOT(K)
                  D3(K,I) = D3(K,I) + F3DOT(K)
   55         CONTINUE
*
*       Advance the neighbour list until last member is identified.
              IF (L.LE.NNB) THEN
                  L = L + 1
                  NAMEJ = LIST(L,I)
              END IF
          END IF
   60 CONTINUE
*
   70 CONTINUE
*        Distribute variables into private vectors again T3D (R.Sp.)
      isend = rank + 1
      if(isend.eq.isize)isend = 0
      irecv = rank - 1
      if(irecv.eq.-1)irecv = isize - 1
*
      do 1001 ir = 0,isize-2
*
      irank = rank - ir
      if(irank.lt.0)irank=irank+isize
*
      istart=ista(irank+1)
      icnt = inum(irank+1)
*
      if(irank.eq.0)irank=isize
      istrec = ista(irank)
      icnt2 = inum(irank)
*
*     print*,' FPOLY2: rank,irank,isend,irecv,istrec=',
*    *    rank,irank,isend,irecv,istrec
*
*     print*,' FPOLY2: bef rank,irank=',rank,irank
*     print*,' FPOLY2: bef rank ',rank,' d2(',istart,')=',d2(1,istart)
*     print*,' FPOLY2: bef rank ',rank,' d2(',istrec,')=',d2(1,istrec)
*     print*,' FPOLY2: rank=',rank,' sending ',icnt,' items to ',isend
*     print*,' FPOLY2: rank=',rank,' recving ',icnt2,
*    *  ' items from ',irecv
*
*      call mpi_barrier(MPI_COMM_WORLD,ierr)
      CALL MPI_SENDRECV(D2(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  D2(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(D2R(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  D2R(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(D3(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  D3(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(D3R(1,istart),3*icnt,MPI_REAL8,isend,rank,
     *                  D3R(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_WORLD,status,ierr)

*     print*,' FPOLY2: aft rank ',rank,' d2(',istart,')=',d2(1,istart)
*     print*,' FPOLY2: aft rank ',rank,' d2(',istrec,')=',d2(1,istrec)
*     print*,' FPOLY2: aft rank,irank=',rank,irank
*
      
      call cputim(tt998)
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       call cputim(tt999)
       ibarcount=ibarcount+1
       ttbar = ttbar + (tt999-tt998)*60

 1001 continue
*     
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNLD(I1,I2,2)
      END IF
*
*       Set new time-steps and initialize prediction variables.
      CALL STEPS(I1,I2)
*
      RETURN
*
      END
