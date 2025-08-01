C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
      SUBROUTINE CALC_REG_ON_GPU(IREG,istart,iend)

      INCLUDE 'common6.h'

      COMMON/BARR/ ibarcount
      PARAMETER (maxthr=1024)
      INTEGER*4 IREG(NMAX), istart, iend
      REAL*8 h2_i(maxthr)
      REAL*8 x_i(3,maxthr), v_i(3,maxthr)
      REAL*8 p_i(maxthr), a_i(3,maxthr), jerk_i(3,maxthr)
      INTEGER LISTGP(LMAX,maxthr)
      INTEGER NBFLAG(NMAX), NOFL2
*     --03/19/13 20:27-lwang-improve------------------------------------*
***** Note:add velocity criterion--------------------------------------**
      REAL*8 DTR(maxthr)
*     --03/19/13 20:27-lwang-end-improve--------------------------------*

      NOFL2= 0
      NOFL = 0
C---  calculate the pot, a & jerk on the GPU

*     WRITE(6,*) 'LOOPREG rank t istart,iend ', RANK, time,istart,iend
*     WRITE(6,*) ' IREG=',(IREG(K),k=istart,iend)
*     CALL FLUSH(6)
*      WRITE(6,*) 'RANK ISTART IEND BLOCKS',rank, istart, iend, 
*     &     iend-istart
      DO L=istart,iend,maxthr

 550     ni = maxthr
         IF( (L+ni).GT.iend ) ni = iend - L + 1

C---  fill the GPU "working" arrays

         DO ii=1,ni
            idi       = IREG(L+ii-1)
CCC   Obtain irregular & regular force and determine current neighbours.
CCC   Take volume between inner and outer radius equal to basic sphere.
            RS2      = RS(idi)**2
*     --02/02/13 18:26-lwang-debug--------------------------------------*
***** Note:change neighbor radius--------------------------------------**
***** Note:Add DTR block step------------------------------------------**
            RCRIT2   = 1.59*RS2
            h2_i(ii) = RS2
            DTR(ii)=STEPR(idi)
*     --02/02/13 18:26-lwang-end-debug----------------------------------*
*     WRITE(6,*) RANK, 'GPU UTIL h2', ii-1, idi, sqrt(h2_i(ii)), 
*     &NAME(idi), LIST(1,idi)
            DO k=1,3
               x_i(k,ii) = X(k,idi)
               v_i(k,ii) = XDOT(k,idi)
            ENDDO
         ENDDO

*     WRITE(6,*) ' call gpu L ni lmax nbmax = ',L,ni, lmax,nnbmax
*     WRITE(6,*) ' h2,x,v=',h2_i(1),x_i(1,1),v_i(1,1)
*     CALL FLUSH(6)

C---  GPU first call with first values of p_i, a_i & jerk_i
         call cputim(tt53)
         CALL gpunb_regf(ni,h2_i,DTR,x_i,v_i,a_i,jerk_i,p_i,lmax,
     &        nnbmax,LISTGP)
         call cputim(tt54)
         ttgpu = ttgpu + (tt54-tt53)*60.0

         gpu_calls = gpu_calls + 1
C---  Save new neighbour list on intermediate vector for regint
         DO 54 ii = 1,NI
            idi = IREG(L+ii-1)
            NNB = LISTGP(1,ii)
*     Warning if neighbour list overflows.
            IF (NNB.LT.0) THEN
               WRITE (6,56)  ii, NI, istart, iend, 
     *              idi, NAME(idi), LIST(1,idi), NNB, RS(idi), time
 56            FORMAT (' OVERFLOW! r ii,NI,istart,iend=',4I6,
     *              ' i NAME NB0 NBP RS t ',4I6,F8.2,1P,D12.5)
               CALL FLUSH(6)
*     Define all particles as neighbours for overflow for full force loop in regint.
*     LIST_G(1,idi) = NNB
               NOFL = NOFL + 1
               RS(idi)=0.9*RS(idi)
               NOFL2 = NOFL2 + 1
            END IF
            IF (NOFL2.GT.0) THEN
               NOFL2 = 0
               GO TO 550
            END IF
            L1 = 1
            DO 53 LL = 2,NNB+1
*     Note GPU address starts from 0 (hence add IFIRST to neighbour list).
               IF (LISTGP(LL,ii)+IFIRST.NE.idi) THEN
                  L1 = L1 + 1
                  LIST_G(L1,idi) = LISTGP(LL,ii) + IFIRST
                  NBFLAG(LIST_G(L1,idi)) = 1
               END IF
 53         CONTINUE
            LIST_G(1,idi) = L1 - 1
            PHIDBL(idi) = p_i(ii)

            DO k=1,3
               FRG(k,idi)    = a_i(k,ii)
               FRGDOT(k,idi) = jerk_i(k,ii)
            ENDDO
 54      CONTINUE
      ENDDO  

#ifdef PARALLEL
      call cputim(tt998)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ttbar = ttbar + (tt999-tt998)*60
      ibarcount=ibarcount+1
#endif
      
      RETURN
C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
      END
