      SUBROUTINE cmcorr(N,m,x,v)
      IMPLICIT NONE

C---  input vars
      INTEGER N,NDIM
      PARAMETER (NDIM=3)
      DOUBLE PRECISION m(N),x(NDIM,N),v(NDIM,N)

C---  local vars
      INTEGER          i,k
      DOUBLE PRECISION mcm,xcm(NDIM),vcm(NDIM)


C---  init local vars
      mcm = 0.d0
      DO k=1,NDIM
         xcm(k) = 0.d0
         vcm(k) = 0.d0
      ENDDO


C---  sum mass weighted positions and velocities
      DO i=1,N
         mcm = mcm + m(i)
         DO k=1,NDIM
            xcm(k) = xcm(k) + m(i)*x(k,i)
            vcm(k) = vcm(k) + m(i)*v(k,i)
         ENDDO
      ENDDO


C---  get cm data
      DO k=1,NDIM
         xcm(k) = xcm(k)/mcm
         vcm(k) = vcm(k)/mcm
      ENDDO


C---  set cm to origin 
      DO i=1,N
         mcm = mcm + m(i)
         DO k=1,NDIM
            x(k,i) = x(k,i) - xcm(k)
            v(k,i) = v(k,i) - vcm(k)
         ENDDO
      ENDDO

#ifdef DEBUG
      WRITE(0,*) 'corrected body data for center of mass'
      WRITE(0,*) 'xcm = ', xcm 
      WRITE(0,*) 'xcm = ', xcm 
      CALL flush(0)
C          -----
#endif

      RETURN
      END


