
C***********************************************************************
C
C
                        SUBROUTINE ellan
C
C
C***********************************************************************
C
C
C     Subroutine to analyze the ellipticity of the system
C
C
C=======================================================================

        INCLUDE 'common6.h'
C   Declaration of local variables.
C   -------------------------------
        DOUBLE PRECISION etoti(1:nmax),mtot,poscm(3),
     &     ti(3,3),tiwork(3,3),dwork(3),ework(3),lam(3)
        INTEGER i,j,k,ief,nbound,nstart,nnext,np,indexev(3)
        INTEGER index(1:nmax)

        INTEGER nef,nef1
        PARAMETER(nef=9,nef1=nef+1)
        DOUBLE PRECISION xf(nef1),ba(nef1),ca(nef1),taue(nef1),
     &            evec(3,3,nef1)

        PARAMETER (tiny2=1.D-30)
        DATA (xf(i),i=1,nef1) /0.01D0,0.05D0,0.1D0,0.2D0,0.5D0,
     &                        0.8D0,0.9D0,0.95D0,1.0D0,99.D0/

C=======================================================================
        
C       calculate specific energy of particles
C       --------------------------------------

        DO 100 i=ifirst,ntot
        etoti(i-ifirst+1) = 0.5D0 * (xdot(1,i)**2 + xdot(2,i)**2 +
     &                         xdot(3,i)**2) + phidbl(i)
100     CONTINUE

C      calculate number of bound particles
C      -----------------------------------

        nbound = 0
        DO 150 i=1,ntot-ifirst+1
           IF(etoti(i).LT.0.D0) nbound = nbound + 1
150     CONTINUE

C       sort for particle energy
C       ------------------------

        CALL indexx(ntot-ifirst+1,etoti,index)

C       initialize tensor of inertia
C       ----------------------------

        DO 210 i=1,3
           DO 200 k=1,3
              ti(i,k) = 0.D0
200        CONTINUE
210     CONTINUE

C       LOOP over fraction of most bound particles and all particles
C       ------------------------------------------------------------

        nstart   = 1
        mtot     = 0.D0
        poscm(1) = 0.D0
        poscm(2) = 0.D0
        poscm(3) = 0.D0
        DO 500 ief=1,nef1

           IF(ief.LE.nef) THEN
C                                  only fraction of bound particles
C                                  --------------------------------
              nnext = NINT(xf(ief) * nbound)
           ELSE
C                                   all particles
C                                   -------------
              nnext = ntot
           ENDIF

C-----------------------------------------------------------------
C--      at least two particles are required for ellipticity...
C-----------------------------------------------------------------
           IF(nnext.LT.2) THEN
              ba(ief)  = 999.
              ca(ief)  = 999.
              taue(ief) = 999.
              DO 320 k=1,3
                 DO 310 j=1,3
                    evec(k,j,ief) = 0.
310              CONTINUE
320           CONTINUE

           ELSE

C       calculate tensor of inertia
C       ---------------------------
*             print*,' nstart,nnext=',nstart,nnext
              DO 400 i=nstart,nnext
                 ipo = index(i) + ifirst - 1
                 ti(1,1) = ti(1,1) + body(ipo) *
     &                  (x(2,ipo)*x(2,ipo) + x(3,ipo)*x(3,ipo))
                 ti(2,2) = ti(2,2) + body(ipo) *
     &                  (x(1,ipo)*x(1,ipo) + x(3,ipo)*x(3,ipo))
                 ti(3,3) = ti(3,3) + body(ipo) *
     &                  (x(1,ipo)*x(1,ipo) + x(2,ipo)*x(2,ipo))
                 ti(1,2) = ti(1,2) - body(ipo) * x(1,ipo)*x(2,ipo)
                 ti(1,3) = ti(1,3) - body(ipo) * x(1,ipo)*x(3,ipo)
                 ti(2,3) = ti(2,3) - body(ipo) * x(2,ipo)*x(3,ipo)
400           CONTINUE
      
C       correct for center of mass
C       --------------------------

C   A) calculate center of mass data
C   --------------------------------

C--       remove previous correction for center of mass
              ti(1,1) = ti(1,1) + mtot * (poscm(2)**2+poscm(3)**2)
              ti(2,2) = ti(2,2) + mtot * (poscm(1)**2+poscm(3)**2)
              ti(3,3) = ti(3,3) + mtot * (poscm(1)**2+poscm(2)**2)
              ti(1,2) = ti(1,2) - mtot * poscm(1) * poscm(2)
              ti(1,3) = ti(1,3) - mtot * poscm(1) * poscm(3)
              ti(2,3) = ti(2,3) - mtot * poscm(2) * poscm(3)
              poscm(1) = poscm(1) * mtot
              poscm(2) = poscm(2) * mtot
              poscm(3) = poscm(3) * mtot
*             xav = 0.d0
*             yav = 0.d0
*             zav = 0.d0

              DO 405 i=nstart,nnext
                 ipo = index(i) + ifirst - 1
                 poscm(1) = poscm(1) + body(ipo) * x(1,ipo)
                 poscm(2) = poscm(2) + body(ipo) * x(2,ipo)
                 poscm(3) = poscm(3) + body(ipo) * x(3,ipo)
                 mtot     = mtot + body(ipo)
*                xav = xav + abs(x(1,ipo))
*                yav = yav + abs(x(2,ipo))
*                zav = zav + abs(x(3,ipo))

 405           CONTINUE
              poscm(1) = poscm(1) / mtot
              poscm(2) = poscm(2) / mtot
              poscm(3) = poscm(3) / mtot
*             print*,' av=',xav,yav,zav

 
C   B) apply correction
C   -------------------
              ti(1,1) = ti(1,1) - mtot * (poscm(2)**2+poscm(3)**2)
              ti(2,2) = ti(2,2) - mtot * (poscm(1)**2+poscm(3)**2)
              ti(3,3) = ti(3,3) - mtot * (poscm(1)**2+poscm(2)**2)
              ti(1,2) = ti(1,2) + mtot * poscm(1) * poscm(2)
              ti(1,3) = ti(1,3) + mtot * poscm(1) * poscm(3)
              ti(2,3) = ti(2,3) + mtot * poscm(2) * poscm(3)

C       set off-axis values by symmetry
C       -------------------------------

              ti(2,1) = ti(1,2)
              ti(3,1) = ti(1,3)
              ti(3,2) = ti(2,3)
*
*             print*,' mtot,poscm=',mtot,(poscm(k),k=1,3)
C=======================================================
C       determine eigenvalues and axis of inertia
C=======================================================

C------------------------------------------------------
C--          copy tensor of inertia
C------------------------------------------------------
              DO 420 i=1,3
                 DO 410 k=1,3
                    tiwork(i,k) = ti(i,k) 
410              CONTINUE
420           CONTINUE
              np = 3

C------------------------------------------------------
C--          calculate eigenvalues and eigenvectors
C------------------------------------------------------
              CALL tred2(tiwork,np,np,dwork,ework)
              CALL tqli(dwork,ework,np,np,tiwork)

C--               sort for increasing eigenvalues
              CALL indexx(np,dwork,indexev)
C--               find eigenvectors
              DO 450 i=1,np
                 lam(i) = dwork(indexev(i))
                 DO 440 k=1,np
                    evec(k,i,ief) = tiwork(k,indexev(i))
440              CONTINUE
450           CONTINUE

              xhelp    = lam(3) + lam(2) - lam(1)
              xhelp1   = lam(2) - lam(3) + lam(1)
c              IF(xhelp1.LT.0.D0) THEN
c                 PRINT*,' ellan: xhelp1 < 0',xhelp1,tnow
c                 xhelp1 = 0.D0
c             ENDIF
              ba(ief)  = SQRT(MAX(tiny2,lam(3)-lam(2)+lam(1)) / xhelp)
              ca(ief)  = SQRT(MAX(tiny2,xhelp1) / xhelp) 
              taue(ief) = (ba(ief)-ca(ief)) /MAX(tiny2,(1.D0 - ca(ief)))

              nstart = nnext + 1

           ENDIF

500     CONTINUE

C==================================================================
C==         OUTPUT of data
C===================================================================

      	DO 600 ief=1,nef1
          if(rank.eq.0)
     *     WRITE(60+ief,910) time,ba(ief),ca(ief),taue(ief),
     *     mtot,(poscm(k),k=1,3)
600     CONTINUE

*     900     FORMAT(1x,1p,e12.5,1x,0p,i9,1x,i9,1x,i2)
910     FORMAT(1x,13(f9.5,1x))

        R E T U R N
        END

