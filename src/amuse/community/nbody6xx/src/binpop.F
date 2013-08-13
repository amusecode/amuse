      SUBROUTINE BINPOP
C
C P. Kroupa   1.04.94
C
C NOTE: This routine works for <= 100 per cent binaries as the initial
C set up. It allows choice of an
C initial period distribution and allows eigenevolution of orbital parameters.
C It also MERGES the closest binaries to one star.
C
C NOTE: A more realistic code for pre-main sequence eigenevolution 
C is contained in binpop_mardling.f. One day it will have to be 
C included in binpop_pk.f as an aditional pre-ms evoltuion sub-routine.
C
*
*
*       Initial binary population.
*       --------------------------
*
      INCLUDE 'common6.h'
      integer  pkmerge,pk_sing
      REAL*8  XORB(2),VORB(2),XREL(3),VREL(3),PX(3),QX(3),BS(NMAX)
      REAL*8  RAN2
      real*8  pkmergedist,EBpk
* Additional vars for MERGED components:
      real*8  bodypk(nmax),xpk(3,nmax),xdotpk(3,nmax)
      real*8  bodymergepk(nmax),xmergepk(3,nmax),xdotmergepk(3,nmax)
* Additional vars for initially single stars:
      real*8  s_body(nmax),s_x(3,nmax),s_xdot(3,nmax)
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
* recover random number sequence from COMMON variable IDUM1.
      KDUM = IDUM1
*
      pkmerge = 0
*
*
      if(rank.eq.0 .AND. (amusein.eq.0))
     & READ (5,*)  SEMI0, ECC0, RATIO,RANGE, NSKIP, IDORM
*
#if MPIINIT
      CALL MPI_BCAST(SEMI0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ECC0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RATIO,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RANGE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NSKIP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(IDORM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      NBIN = NBIN0
      NBIN1 = NBIN + 1
      if(rank.eq.0)
     &WRITE (6,1)  NBIN, SEMI0, ECC0, RATIO, RANGE, NSKIP, IDORM
    1 FORMAT (/,12X,'BINARIES:   NBIN =',I4,'  A =',F9.6,'  E =',F6.2,
     &              '  RATIO =',F4.1,'  RANGE =',F6.1,'  NSKIP =',I3,
     &              '  IDORM =',I2,/)
*
*       Check type of binary mass distribution (NSKIP, IMF2 or split c.m.).
      IF (NSKIP.EQ.0.OR.KZ(20).GE.2) GO TO 10
      IF (RATIO.EQ.1.0) GO TO 20
*
*       Select binaries from the most massive bodies (frequency NSKIP).
      ILAST = (1 + NSKIP)*NBIN
      JSKIP = 0
      JS = 0
      JB = 1
*
*       Transfer binary masses to first NBIN locations.
      DO 6 I = 2,ILAST
          JSKIP = JSKIP + 1
*       Copy binary mass of body #I to new global location.
          IF (JSKIP.GT.NSKIP) THEN
              JSKIP = 0
              JB = JB + 1
              BODY(JB) = BODY(I)
          ELSE
*       Save next NSKIP masses of single bodies.
              JS = JS + 1
              BS(JS) = BODY(I)
          END IF
    6 CONTINUE
*
*       Restore the single bodies in subsequent locations.
      JS = 0
      DO 8 I = NBIN1,ILAST
          JS = JS + 1
          BODY(I) = BS(JS)
    8 CONTINUE
*
*       Move main variables of all single bodies.
   10 DO 15 I = N,NBIN1,-1
          J = I + NBIN
          BODY(J) = BODY(I)
          DO 12 K = 1,3
              X(K,J) = X(K,I)
              XDOT(K,J) = XDOT(K,I)
   12     CONTINUE
   15 CONTINUE
*
*       Create space for each binary component next to primary.
   20 DO 30 I = NBIN,2,-1
          J = 2*I - 1
          BODY(J) = BODY(I)
          DO 25 K = 1,3
              X(K,J) = X(K,I)
              XDOT(K,J) = XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*=========================================================
*=======================THE LOOP==========================
*=========================================================
*       Introduce binary components from relative motion.
      DO 60 I = 1,NBIN
*
*       Randomize perihelion, node & inclination.
          PI = TWOPI*RAN2(KDUM)
          OMEGA = TWOPI*RAN2(KDUM)
          ZI = 0.5*TWOPI*RAN2(KDUM)
*
*       Set transformation elements (Brouwer & Clemence p. 35).
          PX(1) = COS(PI)*COS(OMEGA) - SIN(PI)*SIN(OMEGA)*COS(ZI)
          QX(1) =-SIN(PI)*COS(OMEGA) - COS(PI)*SIN(OMEGA)*COS(ZI)
          PX(2) = COS(PI)*SIN(OMEGA) + SIN(PI)*COS(OMEGA)*COS(ZI)
          QX(2) =-SIN(PI)*SIN(OMEGA) + COS(PI)*COS(OMEGA)*COS(ZI)
          PX(3) = SIN(PI)*SIN(ZI)
          QX(3) = COS(PI)*SIN(ZI)
*
*       Specify component masses (copy BODY0 from IMF2 or use RATIO).
          I1 = 2*I - 1
          I2 = 2*I
          IF (KZ(20).GE.2) THEN
              BODY(I1) = BODY0(I1)
              BODY(I2) = BODY0(I2)
          ELSE IF (RATIO.EQ.1.0) THEN
              BODY(I2) = BODY(I1)
          ELSE
              BODY(I1) = RATIO*BODY(I1)
              BODY(I2) = BODY(I1)*(1.0 - RATIO)/RATIO
          END IF
*
***************************
* Choose
*
* ECC0<0:     thermal distribution f(e)=2e
* 0<=ECC0<=1  fixed eccentricity
* ECC0=20     a uniform distr. between 0 and 1
* ECC0=30     an exponentially decreasing distribution akin to the data in
*             Fig.6a of Duquennoy and Mayor 1991: f(e)=0.1765/(e*e), 0.15<=e<=1
* ECC0=40:
*             The general distribution is (9.02.94):
*                              f(e) = a e^b   e0<=e<=1
*             with   a = (1+b) / (1-e0^(1+b))
*             Its generating function is
*                              e(X) = [((1+b)/a)X + e0^(1+b)]^(1/(1+b))
*             Choose some a and b below.
*             Get thermal distribution if b=1 with e0=0
*
          IF (ECC0.LT.0.0) THEN
              ECC2 = RAN2(KDUM)
              ECC = SQRT(ECC2)
          else if (ecc0.EQ.20) then
              ecc = RAN2(KDUM)
*
          else if (ecc0.EQ.30) then
              ecc = 1./0.15 - RAN2(KDUM)/0.1765
              ecc = 1./ecc
*
          else if (ecc0.EQ.40) then
              PK_e0 = 0.D0
              PK_b = 1.D0
              PK_a = (1.D0+PK_b) / (1.D0 - PK_e0**(1.D0+PK_b))

              ecc = RAN2(KDUM)
              ecc = ((PK_b+1.D0)/PK_a)*ecc + PK_e0**(1.D0+PK_b)
              ecc = ecc**(1.D0/(1.D0+PK_b))
*
*
          ELSE if (ecc0.GE.0.0.AND.ecc0.LE.1.) then
              ECC = ECC0
          END IF
*
*****************************
*
*       Select semi-major axis from uniform distribution in log(A),
*       or a uniform distribution in the range semi0 (=min. ecc.)
*       to semi0+range (=max. ecc.), or SEMI0.
*
* Use distribution in semi-major axis if
* KZ(42)=0: RANGE>0: uniform distribution in log(semi) between SEMI0 and
*                    SEMI0/RANGE
*           RANGE<0: uniform distribution in semi between SEMI0 and -1*RANGE.
*
* Use distribution in logP(days) if
* KZ(42)=1:  linearly increasing distribution function f=0.03438*logP
* KZ(42)=2:  f=3.5logP/[100+(logP)**2]
* KZ(42)=1,2 are 1st and 2nd iterations!
* KZ(42)=3:  f=2.3(logP-1)/[45+(logP-1)**2] This is a "3rd" iteration when
*            pre-ms evolution is taken into account with KZ(41).NE.0.
* KZ(42)=4:  f=2.5(logP-1)/[45+(logP-1)**2] This is a "34th" iteration when
*            pre-ms evolution is taken into account with KZ(41).NE.0, and
*            RBAR<1.5 (in runs /run34n RBAR=1.5 with KZ(42)=3 lead to
*            wrong qone distr.
* KZ(42)=5:  Duquennoy & Mayor 1991, Gaussian distr. 
*            with mean logP=4.8, SDEV in logP=2.3. Use Num.Recipes
*            routine "gasdev.f" to obtain random deviates given "idum1". 
*========= KZ(42) = 0 ========
          if (KZ(42).EQ.0) then
             IF (RANGE.GT.0.0) THEN
                 EXP1 = RAN2(KDUM)*LOG10(RANGE)
                 SEMI = SEMI0/10.0**EXP1
             else if (range.LT.0.0) then
                 exp1 = RAN2(KDUM)*(-1*range-semi0)
                 semi = semi0 + exp1
             ELSE
                 SEMI = SEMI0
             END IF
*========= KZ(42) = 1 ========
          else if (KZ(42).EQ.1) then
* exp = logP(days)
             exp1 = SQRT(RAN2(KDUM)*2.0/0.03438)
* Period in yrs:
             exp1 = 10**exp1 /365.25
* transform to semi-major axis in model units to continue prgrm from here
* as before (i.e. when semi-major axis is chosen): Use Kepler law:
* M1+M2 = a**3/P**2, where mass is in Msun, and a is in AU.
             semi = (BODY(I1)+BODY(I2))*ZMBAR*exp1*exp1
             semi = semi**(1./3.)
* semi in pc and then in modell units:
             semi = semi/206259.591
             semi = semi/RBAR
*========= KZ(42) = 2 ========
          else if (KZ(42).EQ.2) then
             exp1 = EXP(2.0*RAN2(KDUM)/3.50) - 1.
             exp1 = SQRT(exp1*100.)
* exp = logP(days)
* Period in yrs:
             exp1 = 10**exp1 /365.25
* transform to semi-major axis in model units to continue prgrm from here
* as before (i.e. when semi-major axis is chosen): Use Kepler law:
* M1+M2 = a**3/P**2, where mass is in Msun, and a is in AU.
             semi = (BODY(I1)+BODY(I2))*ZMBAR*exp1*exp1
             semi = semi**(1./3.)
* semi in pc and then in modell units:
             semi = semi/206259.591
             semi = semi/RBAR
*========= KZ(42) = 3 ========
          else if (KZ(42).EQ.3) then
             exp1 = EXP(2.D0*RAN2(KDUM)/2.3D0) - 1.D0
             exp1 = SQRT(exp1*45.D0) + 1.D0
* exp = logP(days)
* Period in yrs:
             exp1 = 10**exp1 /365.25D0
* transform to semi-major axis in model units to continue prgrm from here
* as before (i.e. when semi-major axis is chosen): Use Kepler law:
* M1+M2 = a**3/P**2, where mass is in Msun, and a is in AU.
             semi = (BODY(I1)+BODY(I2))*ZMBAR*exp1*exp1
             semi = semi**(1.D0/3.D0)
* semi in pc and then in modell units:
             semi = semi/206259.591D0
             semi = semi/RBAR
*========= KZ(42) = 4 ========
          else if (KZ(42).EQ.4) then
             exp1 = EXP(2.D0*RAN2(KDUM)/2.5D0) - 1.D0
             exp1 = SQRT(exp1*45.D0) + 1.D0
* exp = logP(days)
* Period in yrs:
             exp1 = 10**exp1 /365.25D0
* transform to semi-major axis in model units to continue prgrm from here
* as before (i.e. when semi-major axis is chosen): Use Kepler law:
* M1+M2 = a**3/P**2, where mass is in Msun, and a is in AU.
             semi = (BODY(I1)+BODY(I2))*ZMBAR*exp1*exp1
             semi = semi**(1.D0/3.D0)
* semi in pc and then in modell units:
             semi = semi/206259.591D0
             semi = semi/RBAR
*========= KZ(42) = 5 ========
          else if (KZ(42).EQ.5) then
             exp1 = -10.
C Only accept periods logP>1, i.e. longer than 10 days at birth
C AND periods shorter than 10^(10) days!
             do while (exp1.LT.-1.6522)
                exp1 = GASDEV_pk(KDUM)
                if (exp1.GT.2.2696) exp1=-10.
             end do
* exp1 is now a Gaussian deviate with mean zero and unit variance
* Transform to mean=4.8 and variance = 2.3*2.3, i.e to logP, P in days
            exp1 = 2.3*exp1 + 4.8
* exp = logP(days)
* Period in yrs:
             exp1 = 10**exp1 /365.25D0
* transform to semi-major axis in model units to continue prgrm from here
* as before (i.e. when semi-major axis is chosen): Use Kepler law:
* M1+M2 = a**3/P**2, where mass is in Msun, and a is in AU.
             semi = (BODY(I1)+BODY(I2))*ZMBAR*exp1*exp1
             semi = semi**(1.D0/3.D0)
* semi in pc and then in modell units:
             semi = semi/206259.591D0
             semi = semi/RBAR
*========= KZ(42) = ? ========
          else
             if(rank.eq.0)write(6,*)' No KZ(42) defined!!!!! '
             STOP
          end if
*
*
* At this stage we have the mases, eccentricity and period of each binary
* at "birth", i.e. prior to circularisation and "feeding". Now evolve these
* to very, very roughly take into account complete circularisation,
* partial circularisation and "feeding". Do this if option KZ(41)=1:
* (i.e. mass-exchange at proto-stellar time):
          if (KZ(41).EQ.1) then
             call proto_star_evol(
     &          ZMBAR,RBAR,BODY(I1),BODY(I2),ECC,SEMI)
          end if
*
*       Specify relative motion at apocentre.
          XORB(1) = SEMI*(1.0 + ECC)
          XORB(2) = 0.0
          VORB(1) = 0.0
          ZMBIN = BODY(I1) + BODY(I2)
          VORB(2) = SQRT(ZMBIN*(1.0D0 - ECC)/(SEMI*(1.0D0 + ECC)))
          EBpk =  -0.5*BODY(I1)*BODY(I2)/SEMI
C Moved below to account for accepted binaries only.
C          EBIN0 = EBIN0 - 0.5*BODY(I1)*BODY(I2)/SEMI
*
*       Transform to relative variables.
          DO 40 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
   40     CONTINUE
*
*       Set global variables for each component.
          DO 50 K = 1,3
              X(K,I1) = X(K,I1) + BODY(I2)*XREL(K)/ZMBIN
              X(K,I2) = X(K,I1) - XREL(K)
              XDOT(K,I1) = XDOT(K,I1) + BODY(I2)*VREL(K)/ZMBIN
              XDOT(K,I2) = XDOT(K,I1) - VREL(K)
   50     CONTINUE
*
*
* Merge binaries that have semi.LE.10 solar radii (all in AU):
* Solar radius in AU is 4.6523D-3.
          semi = semi*RBAR*206259.591D0
          pkmergedist = 10.D0 * 4.6523D-3
c          pkmergedist = 50000.D0 * 4.6523D-3
          J1 = I1 -2*pkmerge
          J2 = I2 -2*pkmerge
          if (semi.LE.pkmergedist) then
             pkmerge = pkmerge + 1
             ZMBIN = BODY(I1) + BODY(I2)
C Added 20.8.96:
             if(rank.eq.0)then
             write(6,*)
             write(6,*)' MERGED stars in binpop_pk.f:'
             write(6,*)' mass1, mass2, total mass [Msun]'
             write(6,'(3F8.3)')body(I1)*ZMBAR,body(I2)*ZMBAR,
     +                         zmbin*ZMBAR 
             end if
C End added bit. 
             do K = 1,3
                 Xmergepk(K,pkmerge) = (BODY(I1)*X(K,I1) +
     &                         BODY(I2)*X(K,I2))/ZMBIN
                 XDOTmergepk(K,pkmerge) = (BODY(I1)*XDOT(K,I1) +
     &                         BODY(I2)*XDOT(K,I2))/ZMBIN
             end do
             BODYmergepk(pkmerge) = ZMBIN
          else
C Sum binding energy for accepted binaries.
             EBIN0 = EBIN0 + EBpk
             bodypk(J1) = body(I1)
             bodypk(J2) = body(I2)
             do k=1,3
                xpk(k,J1) = x(k,I1)
                xpk(k,J2) = x(k,I2)
                xdotpk(k,J1) = xdot(k,I1)
                xdotpk(k,J2) = xdot(k,I2)
             end do
          end if
*
*=========================================================
*=================end of THE LOOP=========================
*=========================================================
   60 CONTINUE
*
C Kroupa: 5.11.96
C Take into account possible initial single stars. Before merging above 
C these occupy positions: 2*NBIN+1....N+NBIN (before merging the stars 
C in binaries occupy positions 1,2,....2*NBIN). 
      pk_sing = 0
      do i=2*NBIN+1,N+NBIN
         pk_sing = pk_sing+1
         s_body(pk_sing)=body(i)
         do k=1,3
            s_x(k,pk_sing) = x(k,i)
            s_xdot(k,pk_sing) = xdot(k,i)
         end do
      end do
*
*       Update the total particle number after primary splitting or IMF2
*       and introduce additional counters for MERGING:
      IF (RATIO.LT.1.0.OR.KZ(20).GE.2) THEN
          N = N + NBIN
          NZERO = N
          NTOT = N
          BODYM = ZMASS/FLOAT(N)
          IF (rank.eq.0.and.NSKIP.GT.0) THEN
              WRITE (6,62)  (BODY(J),J=1,10)
              WRITE (6,64)  (BODY(J),J=2*NBIN+1,2*NBIN+10)
   62         FORMAT (/,12X,'BINARY MASSES (1-10):  ',10F9.5)
   64         FORMAT (/,12X,'SINGLE MASSES (1-10):  ',10F9.5,/)
          END IF
      END IF
*
      Npk = N - pkmerge
      NZEROpk = Npk
      NTOTpk = Npk
      NBINpk = NBIN - pkmerge

c+++
CCCCCCCCCCCC for testing:
        if(rank.eq.0)then
        write(6,*)
        write(6,*)' In BINPOP_PK :'
        write(6,*)'****************************************'
        write(6,*)'N,NBIN,NBIN0,NBINpk,pkmerge,Npk,pk_sing:'
        write(6,'(7(I6))')N,NBIN,NBIN0,NBINpk,pkmerge,
     +                    Npk,pk_sing
c        write(6,*)
c        write(6,*)
c        write(6,*)' original stars'
c        do i=1,N
c           write(6,'(I4,7(F8.3))') i,body(i)*ZMBAR,
c     +     (X(k,i),xdot(k,i),k=1,3)
c        end do
c        write(6,*)
c        write(6,*)
c        write(6,*)' remaining binaries'
c        do i=1,2*NBINpk
c           write(6,'(I4,7(F8.3))') i,bodypk(i)*ZMBAR,
c     +     (Xpk(k,i),xdotpk(k,i),k=1,3)
c        end do
c        write(6,*)
c        write(6,*)' merged ones:'
c        do i=1,pkmerge
c           write(6,'(I4,7(F8.3))') i,bodymergepk(i)*ZMBAR,
c     +     (Xmergepk(k,i),xdotmergepk(k,i),k=1,3)
c        end do
c        write(6,*)'****************************************'
c        call flush(6)
         end if
CCCCCCCCCCCCCC
*
* Now overwrite the old vars with the MERGED ones:
      N = Npk
      NTOT = ntotpk
      NZERO = nzeropk
      NBIN0 = NBIN0 - pkmerge
      NBIN = NBINpk
C The stars in remaining binaries:
      do i=1,2*NBIN
         body(i) = bodypk(i)
         do k=1,3
            x(k,i) = xpk(k,i)
            xdot(k,i) = xdotpk(k,i)
         end do
      end do
C Add merged (i.e single) stars to end of the above list of stars in binaries:
      j = 2*NBIN
      do i=1,pkmerge
         j = j+1
         body(j) = bodymergepk(i)
         do k=1,3
            x(k,j) = xmergepk(k,i)
            xdot(k,j) = xdotmergepk(k,i)
         end do
      end do
C Add the initially single stars again at end of the lot:
      do i=1,pk_sing
         j=j+1
         body(j) = s_body(i)
         do k=1,3
            x(k,j)=s_x(k,i)
            xdot(k,j)=s_xdot(k,i)
         end do
      end do
*
C Added 11.08.98 with Sverre Aarseth: (body0 needed for 
C stellar evolution!)
      if (j.ne.ntot) then 
         if(rank.eq.0)
     &   write(6,*)j,ntot,' ****** j & ntot must be equal ******'
         stop
      end if
*
      zmass = 0.D0
      do i=1,ntot
         body0(i) = body(i)
         zmass = zmass + body0(i)
      end do
      BODYM = ZMASS/FLOAT(N)
      if(rank.eq.0)then
      write(6,*)
      write(6,*)' BODYM, ZMASS and BODY0(i) updated.'
      write(6,'(2(a,F10.4))')' BODYM= ',BODYM,' ZMASS= ',ZMASS
      write(6,*)
*
c+++ for testing:
c      write(6,*)
c      write(6,*)
c      write(6,*)' After merging & re-arranging:'
c      write(6,*)' N,NBIN0,NTOT: ',N,NBIN0,NTOT
c      do i=1,N
c         write(6,'(I3,7(F8.3))')i,body(i)*ZMBAR,
c     +   (x(k,i),xdot(k,i),k=1,3)
c      end do
c      call flush(6)
       end if
*
*
*       Include procedure for introducing dormant binaries.
      IF (IDORM.GT.0) THEN
          DO 66 I = 1,NBIN
              I1 = 2*I - 1
              I2 = I1 + 1
              ZMBIN = BODY(I1) + BODY(I2)
              DO 65 K = 1,3
                  X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZMBIN
                  XDOT(K,I) = (BODY(I1)*XDOT(K,I1) +
     &                         BODY(I2)*XDOT(K,I2))/ZMBIN
   65         CONTINUE
              BODY(I) = ZMBIN
   66     CONTINUE
*
*       Move the original single particles up to form compact array.
          I1 = 2*NBIN + 1
          I2 = NBIN
          DO 68 I = I1,N
              I2 = I2 + 1
              BODY(I2) = BODY(I)
              DO 67 K = 1,3
                  X(K,I2) = X(K,I)
                  XDOT(K,I2) = XDOT(K,I)
   67         CONTINUE
   68     CONTINUE
*
*       Reset particle membership and turn off binary output option (if = 1).
          N = N - NBIN
          NZERO = N
          NTOT = N
          NBIN0 = 0
          EBIN0 = 0.0
          IF (KZ(8).EQ.1) KZ(8) = 0
      END IF
*
*       Set coordinates & velocities in c.m. rest frame.
      DO 70 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   70 CONTINUE
*
      DO 80 I = 1,N
          DO 75 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   75     CONTINUE
   80 CONTINUE
*
      DO 90 I = 1,N
          DO 85 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
   85     CONTINUE
   90 CONTINUE
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
      RETURN
*
      END
C=======================================================
      subroutine proto_star_evol(ZMBAR,RBAR,mass1,mass2,ECC,SEMI)
*
      implicit none
      real*8          mass1,mass2,ecc,semi,period
      real*8          ecc_initial,period_initial
      real*8          qnew,qold,mtot,Ro,mtot_initial
      real*8          R_periastron,alpha,beta
      real*8          ZMBAR,RBAR,au,Rsun
* astr. unit, solar radius, all in AU (1pc=206259.591AU)
      parameter(au=206259.591D0,Rsun=4.6523D-3)
*
*
*==============================
* Choose these to define the model: (alpha==lambda, beta==chi)
*
c      alpha = 10.D0
c      beta = 1.D0
* best values:
      alpha = 28.D0
      beta = 0.75D0
*==============================
*
*
*
* in Msun:
      mtot = (mass1+mass2)*ZMBAR
      mtot_initial = mtot
* in AU:
      semi = semi*RBAR*au
* in years:
      period = semi*semi*semi/mtot
      period = DSQRT(period)
      ecc_initial = ecc
      period_initial = period
*
* 1) Circularisation and evolution of orbit as a function of
* periastron distance  (Ro in AU):
* Note that the algorithm used here leads to circularised orbits for
* logP<=1 approximately!! (if beta=1.5,alpha=35 approximately)
      Ro = alpha *Rsun
      R_periastron =  semi*(1.D0-ecc)
      alpha = -1.D0*(Ro/R_periastron)**beta
      if (ecc.GT.0.D0) then
         ecc = DEXP(alpha + DLOG(ecc))
      else
         ecc = ecc_initial
      end if
*
* 2) Change mass-ratio towards unity as a function of initial periastron
* distance:
*
      qold = mass1/mass2
      if (qold.GT.1.D0) qold = 1.D0/qold
      alpha = -1.D0*alpha
      if (alpha.GT.1.D0) then
         qnew = 1.D0
      else
         qnew = qold + (1.D0-qold) * alpha
      end if
*
* new masses in model units (remembering q=m1/m2<1) if mass is conserved
* NOT IMPLEMENTED!
c      mtot = mtot/ZMBAR
c      mass1 = mtot/(qnew+1.D0)
c      mass2 = mtot - mass1
*
* Keep the mass of primary fixed and adjust mass of secondary. Note that this
* algorithm leads to a gain in mass of the binary, and thus of the whole
* cluster!
*
C Added 20.06.96 write statements:
*       if(rank.eq.0)then
*       write(6,*)
*       write(6,*)' FEEDING in binpop_pk.f'
*       write(6,'(a,2F8.3)')' old masses [Msun]:', 
*    +                        mass1*ZMBAR,mass2*ZMBAR 
        mass1 = DMAX1(mass1,mass2)
        mass2 = qnew*mass1
*       write(6,'(a,2F8.3)')' new masses [Msun]:', 
*    +                        mass1*ZMBAR,mass2*ZMBAR 
*       end if
C End added bit.
*
* In Msun:
        mtot = (mass1+mass2)*ZMBAR
*
C This below is wrong as in ecc formula above constant Rperi was assumed!!!
c* Duquennoy et al. 1992 in "Binaries as tracers of stellar evolution":
c      period = period_initial * DEXP((57.D0/14.D0) *
c     & (ecc*ecc - ecc_initial*ecc_initial))
C This below is correct:
       period = period_initial*((1.D0-ecc_initial)/(1.D0-ecc))**1.5D0
       period = period * DSQRT(mtot_initial/mtot)
*
*
* New semi-major axis and back to model units:
      semi = mtot * period*period
      semi = semi**(1.D0/3.D0)
      semi = semi/(RBAR*au)
*
      return
      end
C=======================================================
      FUNCTION gasdev_pk(idum)
      INTEGER idum
      REAL*8 gasdev_pk
CU    USES ran1
C P.Kroupa: 8.11.96: changed to RAN2
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev_pk=v2*fac
        iset=1
      else
        gasdev_pk=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software >).

