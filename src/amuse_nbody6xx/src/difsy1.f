      subroutine difsy1(n,eps,h,x,y)
c          Bulirsch-Stoer integrator.
c          --------------------------
c      For \Gamma=(H-E)/L and Y(ntime)=time
*
      PARAMETER  (NMX=80,NMX7=7*NMX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  Y(N),YA(NMX),YL(NMX),YM(NMX),DY(NMX),DZ(NMX),DT(NMX,7),
     &        D(7),X,XN,H,G,B,B1,U,V,C,TA,W,DABS
      COMMON/INTFAC/  LX,LE,LP,LV,LT,J10,NHALF2
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
*
      LOGICAL  KONV,BO,KL,GR,FYBAD
      DATA DT /NMX7*0.0D0/
*     DATA  EP/0.4D-1,0.16D-2,0.64D-4,0.256D-5/
      SAVE
*
*     nhalf2=(n/2)*2
      jti=0
      fy=1.
      reduc=0.5
       do i=1,n
       ya(i)=y(i)
       end do
      CALL DERQP(Y(1),Y(LX),Y(LE),Y(LP),Y(LV),Y(LT),
     &           DZ(1),DZ(LX),DZ(LE),DZ(LP),DZ(LV),DZ(LT))
      IF (JC.GT.0)H=DSC
   10 xn=x+h
      bo=.false.
      m=1
      jr=2
      js=3
      do  j=1,10 ! Jmax (10 -> efficient,  4 -> short steps)
       if(bo)then
       d(2)=16d0/9.d0
       d(4)=64.d0/9.d0
       d(6)=256.d0/9.d0
       else
       d(2)=2.25d0
       d(4)=9.d0
       d(6)=3.6d1
       end if
       if(j.gt.7)then
       l=7
       d(7)=6.4d1
       else
       l=j
       d(l)=m*m
       end if
      konv=l.gt.3
      m=m+m
      g=h/(m)
      b=g+g
      m=m-1
       do i=1,n
       yl(i)=ya(i)
       ym(i)=ya(i)+g*dz(i)
       end do
      do k=1,m
      CALL DERQP(YM(1),YM(LX),YM(LE),YM(LP),YM(LV),YM(LT),
     &           DY(1),DY(LX),DY(LE),DY(LP),DY(LV),DY(LT))
       do i=1,n
       u=yl(i)+b*dy(i)
       yl(i)=ym(i)
       ym(i)=u
       end do
      end do
      CALL DERQP(YM(1),YM(LX),YM(LE),YM(LP),YM(LV),YM(LT),
     &           DY(1),DY(LX),DY(LE),DY(LP),DY(LV),DY(LT))
      kl=l.lt.2
      gr=l.gt.5
      fs=0.
      do i=1,n
      v=dt(i,1)
      c=(ym(i)+yl(i)+g*dy(i))*0.5d0
      dt(i,1)=c
      ta=c
      if(.not.kl)then
      do k=2,l
      b1=d(k)*v
      b=b1-c
      w=c-v
      u=v
       if(b.ne.0.d0)then
       b=w/b
       u=c*b
       c=b1*b
       end if
      v=dt(i,k)
      dt(i,k)=u
      ta=u+ta
      end do
      is=i+n/2
      if(is.gt.nhalf2)is=i-(n/2)
      dyis=dabs(dy(is))
      if(i.eq.n)dyis=1/(abs(ya(i))+abs(ta)) ! for time
        if(konv)then
        test=dabs( (y(i)-ta)*dyis )
        if(test.gt.eps) konv=.false.
        end if
      if(.not.gr)then
      fv=dabs(w)*dyis
      if(fs.lt.fv) fs=fv
      end if
      end if
      y(i)=ta
      end do
       if(fs.ne.0.d0)then
       fa=fy
       k=l-1
       fy=(ep(k)/fs)**(1./(l+k))
       fa7=0.7*fa
       if(l.eq.2)fa7=0.0
       fybad=.not.((fa7.gt.fy).or.(fy.gt.0.7))
        if(fy bad)then
        h=h*fy
        jti=jti+1
         if(jti.gt.10)then !Seppo's suggestion Oct/06 (5 before).
         h=0.0
          do i=1,n
          y(i)=ya(i)
          end do
         return
         end if
        goto 10
        end if
       end if
      if(konv)then
      h=xn-x
      x=xn
*     if(fy.gt.10.0)fy=10. ! factor 10 may be too large
*     if(fy.gt.4.0)fy=4. ! factor 10 may be too large
      if(fy.gt.2.0)fy=2. ! factor 10 may be too large
      h=h*fy
      return
      end if
      d(3)=4.d0
      d(5)=1.6d1
      bo=.not.bo
      m=jr
      jr=js
      js=m+m
      end do
      h=reduc*h
      reduc=0.01+reduc*reduc
      goto 10
      end
