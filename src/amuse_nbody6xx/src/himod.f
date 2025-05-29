      SUBROUTINE HIMOD(evec,hvec,m3,mt,EOvec,HOvec,ICALL,dt0,dt,e,r,v)
*     
*     
*     Modification of hierarchical binary.
*     ------------------------------------
*     
*     Author:  Rosemary Mardling (3/98).

      IMPLICIT REAL*8 (A-H,M,O-Z)
      common/rksave/  coeff,HOhat(3),e0,a,hh,mb
      REAL*8  evec(3),hvec(3),EOvec(3),HOvec(3),
     &     element(6),r(3),v(3),u(6),udot(6)
*     
*     
*     Given inner Runge-Lenz vector {\bf e}=evec and {\bf h}=hvec.
*     Calculate eccentricity e, magnitude of hvec and semi-major axis a.

      mb=mt-m3
      e0=sqrt(dot(evec,evec))
      hh=sqrt(dot(hvec,hvec))
      a=hh**2/mb/(1-e0**2)

*     Given OUTER {\bf E}=Eout and {\bf H}=Hout.
*     Calculate outer eccentricity and semi-major axis.

      Eout=sqrt(dot(EOvec,EOvec))
      HHout=sqrt(dot(HOvec,HOvec))
      Aout=HHout**2/mt/(1-Eout**2)
      do i=1,3
         HOhat(i)=HOvec(i)/HHout
      enddo

*     Calculate the semi-minor axis Bout and the semi-latus rectum SLRout.

      Bout=Aout*sqrt(1-Eout**2)
      SLRout=Aout*(1-Eout**2)
      coeff=m3/(2*Aout*Bout*SLRout)
      DTSUM = 0.0

      IF (ICALL.LT.0) THEN
         do i=1,3
            u(i)=evec(i)
            u(i+3)=hvec(i)
         enddo
         CALL DERIV(u,udot)
         TAU = e0/sqrt(udot(1)**2+udot(2)**2+udot(3)**2)
         WRITE (6,5)  e0,TAU,a,Eout,Aout,(UDOT(K),K=1,3)
    5    FORMAT (' DERIV   e0 TAU a E1 A1 edot  ',
     &        F8.4,1P,2E10.2,0P,F8.4,1P,4E10.2)
         CALL FLUSH(3)
         ICALL = ICALL + 1
      END IF

*     Initialize integration variables.

 10   do i=1,3
	 u(i)=evec(i)
	 u(i+3)=hvec(i)
      enddo

*     Perform integration.

      call rkint(dt,u)

      if(e0.ge.1.0) goto 20
      eh=0.0
      do i=1,3
	 evec(i)=u(i)
	 hvec(i)=u(i+3)
         eh=eh+evec(i)*hvec(i)
      enddo

      e0=sqrt(dot(evec,evec))
      hh=sqrt(dot(hvec,hvec))

*     Apply symmetric re-orthogonalization (small eh; Seppo Mikkola 3/98).

      eps=-eh/(e0**2+hh**2) 
      eps=0.0
      do k=1,3
         evec(k)=evec(k)+eps*hvec(k)
         hvec(k)=hvec(k)+eps*evec(k)
      end do

      IF (DTSUM + DT.GT.DT0) DT = DT0 - DTSUM + 1.0D-10
      DTSUM = DTSUM + DT
      IF (DTSUM.LT.DT0) GO TO 10

*     Calculate new orbital parameters.
 20   e2=dot(evec,evec)
      h2=dot(hvec,hvec)
      a=h2/mb/(1-e2)
      e=sqrt(e2)
      hh=sqrt(h2)

*     Calculate updated vectors r and v. 
*     Note that the vectors evec and hvec give no information about
*     the orbital phase. Thus there is a kind of degeneracy. Since
*     we are updating to a later time using an approximation, we can't
*     expect to know it anyway. Hence the phase of peri is set to zero.

      call transform4(evec,hvec,mb,element)
      call transform2(element,mb,r,v)

      RETURN
      END

      
      subroutine cross(u,v,w)

*     Vectorial cross product.
*     ------------------------
      real*8 u(3),v(3),w(3)

      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)

      end


      real*8 function dot(u,v)

*     Scalar dot product.
*     ------------------
      real*8 u(3),v(3)

      dot=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)

      end


      subroutine transform2(element,mb,r,v)
*     
*     Calculates vectors r and v given 6 orbital elements.
*     ----------------------------------------------------
*     
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 r(3),v(3),element(6),mb
      real*8 tempr(3),tempv(3)
      real*8 inc,B(3,3)


      a=element(1)
      e=element(2)
      inc=element(3)
      w=element(4)
      Om=element(5)
      phi=element(6)

      cosp=cos(phi)
      sinp=sin(phi)
      cosi=cos(inc)
      sini=sin(inc)
      cosw=cos(w)
      sinw=sin(w)
      cosOm=cos(Om)
      sinOm=sin(Om)

      B(1,1) = cosw*cosOm - sinw*cosi*sinOm
      B(1,2) = cosw*sinOm + sinw*cosi*cosOm
      B(1,3) = sinw*sini
      B(2,1) = -sinw*cosOm - cosw*cosi*sinOm
      B(2,2) = -sinw*sinOm + cosw*cosi*cosOm
      B(2,3) = cosw*sini
      B(3,1) = sini*sinOm
      B(3,2) = -sini*cosOm
      B(3,3) = cosi

      h=sqrt(mb*a*(1-e**2))
      rr=a*(1-e**2)/(1+e*cosp)
      rd=e*h*sinp/(a*(1-e**2))
      phid=h/rr**2

      r(1)=rr*cosp
      r(2)=rr*sinp
      r(3)=0.0

      v(1)=rd*cosp-rr*phid*sinp
      v(2)=rd*sinp+rr*phid*cosp
      v(3)=0.0

      do i=1,3
         sum1=0.0
         sum2=0.0
         do j=1,3
            sum1=sum1 + B(j,i)*r(j)
            sum2=sum2 + B(j,i)*v(j)
         enddo
         tempr(i)=sum1
         tempv(i)=sum2
      enddo
      do i=1,3
         r(i)=tempr(i)
         v(i)=tempv(i)
      enddo

      end


      subroutine transform4(e,h,mb,element)
*     
*     Calculates 5 orbital elements given vectors evec and h and mass mb.
*     -------------------------------------------------------------------
*     
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 element(6)
      real*8 h(3),e(3),n(3),hn(3),he(3)
      real*8 eh(3)
      real*8 ii(3),jj(3),kk(3),nn,mb
      data ii/1.d0,0.d0,0.d0/
      data jj/0.d0,1.d0,0.d0/
      data kk/0.d0,0.d0,1.d0/


      hh=sqrt(dot(h,h))

      do i=1,3
         h(i)=h(i)/hh
      enddo

      call cross(kk,h,n)

      nn=sqrt(dot(n,n))
      ecc=sqrt(dot(e,e))
      a=hh**2/mb/(1-ecc**2)

      do i=1,3
         n(i)=n(i)/nn
         eh(i)=e(i)/ecc
      enddo

      call cross(h,n,hn)
      call cross(h,eh,he)

      cosOm=dot(ii,n)
      sinOm=dot(jj,n)
      cosi=dot(kk,h)
      sini=nn
      cosw=dot(n,eh)
      sinw=dot(eh,hn)

      element(1)=a
      element(2)=ecc
      element(3)=atan2(sini,cosi)
      element(4)=atan2(sinw,cosw)
      element(5)=atan2(sinOm,cosOm)
      element(6)=0.0
*     Pericentre phase (can be specified).

*     ZI = 360.0*ELEMENT(3)/(8.0*ATAN(1.0D0))
*     WRITE (6,5)  ECC,ELEMENT(3),ELEMENT(4),ELEMENT(5),ZI
*     5 FORMAT (' TRANSFORM    e i w Om IN ',F8.4,F9.3,3F9.3)
      end
