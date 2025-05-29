        function transform1new(m1,m2,r,v,element)
*
*       Calculates 6 orbital elements given vectors r and v and mass ratio q.
*       --------------------------------------------------------------------
*       updated for arbitrary m1 3-1-05
*       28/red/7
*       2/12/97     Santa Cruz
*
        implicit real*8 (a-h,m,o-z)
        real*8 r(3),v(3),element(6)
        real*8 h(3),e(3),n(3),hn(3),rh(3),he(3)
        real*8 ii(3),jj(3),kk(3),nn,dot
        data ii/1.d0,0.d0,0.d0/
        data jj/0.d0,1.d0,0.d0/
        data kk/0.d0,0.d0,1.d0/
        INTEGER transform1new

c       write(*,*) "Transform:", m1,m2,r,v
        call cross(r,v,h)
        call cross(v,h,e)

        rr=sqrt(dot(r,r))
        hh=sqrt(dot(h,h))

        do i=1,3
           e(i)=e(i)/(m1+m2)-r(i)/rr
           h(i)=h(i)/hh
           rh(i)=r(i)/rr
        enddo

        call cross(kk,h,n)

        nn=sqrt(dot(n,n))
        ecc=sqrt(dot(e,e))
        a=hh**2/(m1+m2)/(1.d0-ecc**2)


        do i=1,3
           if(nn.eq.0)then
              n(i)=ii(i)
           else
              n(i)=n(i)/nn
           endif
           e(i)=e(i)/ecc
        enddo

        call cross(h,n,hn)
        call cross(h,e,he)

        cosOm=dot(ii,n)
        sinOm=dot(jj,n)
        cosi=dot(kk,h)
        sini=nn
        cosw=dot(n,e)
        sinw=dot(e,hn)
        cosphi=dot(rh,e)
        sinphi=dot(rh,he)

        element(1)=a
        element(2)=ecc
        element(3)=atan2(sini,cosi)
        element(4)=atan2(sinw,cosw)
        element(5)=atan2(sinOm,cosOm)
        element(6)=atan2(sinphi,cosphi)

c       write(*,*) "Transform2:", element
        transform1new = 0
        end


        subroutine transform1(q,r,v,element)
*
*       Calculates 6 orbital elements given vectors r and v and mass ratio q.
*       --------------------------------------------------------------------
*       28/red/7
*       2/12/97     Santa Cruz
*
        implicit real*8 (a-h,m,o-z)
        real*8 r(3),v(3),element(6)
        real*8 h(3),e(3),n(3),hn(3),rh(3),he(3)
        real*8 ii(3),jj(3),kk(3),nn,dot
        data ii/1.d0,0.d0,0.d0/
        data jj/0.d0,1.d0,0.d0/
        data kk/0.d0,0.d0,1.d0/

        call cross(r,v,h)
        call cross(v,h,e)

        rr=sqrt(dot(r,r))
        hh=sqrt(dot(h,h))

        do i=1,3
           e(i)=e(i)/(1.d0+q)-r(i)/rr
           h(i)=h(i)/hh
           rh(i)=r(i)/rr
        enddo

        call cross(kk,h,n)

        nn=sqrt(dot(n,n))
        ecc=sqrt(dot(e,e))
        a=hh**2/(1.d0+q)/(1.d0-ecc**2)

        do i=1,3
           if(nn.eq.0)then
              n(i)=ii(i)
           else
              n(i)=n(i)/nn
           endif
           e(i)=e(i)/ecc
        enddo

        call cross(h,n,hn)
        call cross(h,e,he)

        cosOm=dot(ii,n)
        sinOm=dot(jj,n)
        cosi=dot(kk,h)
        sini=nn
        cosw=dot(n,e)
        sinw=dot(e,hn)
        cosphi=dot(rh,e)
        sinphi=dot(rh,he)

        element(1)=a
        element(2)=ecc
        element(3)=atan2(sini,cosi)
        element(4)=atan2(sinw,cosw)
        element(5)=atan2(sinOm,cosOm)
        element(6)=atan2(sinphi,cosphi)

        end


        function transform2(mtot,element,r,v)
*
*       Calculates vectors r and v given 6 orbital elements.
*
        implicit real*8 (a-h,m,o-z)

        real*8 r(3),v(3),element(6)
        real*8 tempr(3),tempv(3)
        real*8 inc,B(3,3),mtot
        INTEGER transform2

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

        h=sqrt(mtot*a*(1-e**2))
        rr=a*(1-e**2)/(1+e*cosp)
        rd=e*h*sinp/(a*(1-e**2))
        phid=h/rr**2

        r(1)=rr*cosp
        r(2)=rr*sinp
        r(3)=0

        v(1)=rd*cosp-rr*phid*sinp
        v(2)=rd*sinp+rr*phid*cosp
        v(3)=0

        do i=1,3
           sum1=0
           sum2=0
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

        transform2 = 0
        end

        subroutine cross(u,v,w)
        real*8 u(3),v(3),w(3)
        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
        end

        real*8 function dot(u,v)
        real*8 u(3),v(3)
        dot=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
        end

        subroutine transform3(q,eQD,r,v,element)
*
*       Calculates 6 orbital elements given vectors r and v and mass ratio q.
*       --------------------------------------------------------------------
*       28/red/7
*       2/12/97
*
        implicit real*8 (a-h,m,o-z)

        real*8 r(3),v(3),element(6)
        real*8 h(3),e(3),n(3),hn(3),rh(3),he(3)
        real*8 ii(3),jj(3),kk(3),nn,dot
        data ii/1.d0,0.d0,0.d0/
        data jj/0.d0,1.d0,0.d0/
        data kk/0.d0,0.d0,1.d0/

        call cross(r,v,h)
        call cross(v,h,e)

        rr=sqrt(dot(r,r))
        hh=sqrt(dot(h,h))

        do i=1,3
           e(i)=e(i)/((1+q)*(1+eQD))-r(i)/rr
           h(i)=h(i)/hh
           rh(i)=r(i)/rr
        enddo

        call cross(kk,h,n)

        nn=sqrt(dot(n,n))
        ecc=sqrt(dot(e,e))
        a=hh**2/(1+q)/(1-ecc**2)

        do i=1,3
           n(i)=n(i)/nn
           e(i)=e(i)/ecc
        enddo

        call cross(h,n,hn)
        call cross(h,e,he)

        cosOm=dot(ii,n)
        sinOm=dot(jj,n)
        cosi=dot(kk,h)
        sini=nn
        cosw=dot(n,e)
        sinw=dot(e,hn)
        cosphi=dot(rh,e)
        sinphi=dot(rh,he)

        element(1)=a
        element(2)=ecc
        element(3)=atan2(sini,cosi)
        element(4)=atan2(sinw,cosw)
        element(5)=atan2(sinOm,cosOm)
        element(6)=atan2(sinphi,cosphi)

        end


        subroutine transform4(e,h,q,element)
*
*       Calculates 5 orbital elements given vectors evec and h and mass ratio q.
*       -----------------------------------------------------------------------
*
        implicit real*8 (a-h,m,o-z)

        real*8 element(5)
        real*8 h(3),e(3),n(3),hn(3),he(3)
        real*8 eh(3),hhat(3),dot
        real*8 ii(3),jj(3),kk(3),nn,nnn
        data ii/1.d0,0.d0,0.d0/
        data jj/0.d0,1.d0,0.d0/
        data kk/0.d0,0.d0,1.d0/

        hh=sqrt(dot(h,h))

        call cross(kk,h,n)

        nn=sqrt(dot(n,n))
        if(nn.eq.0)then
           do i=1,3
              n(i)=ii(i)
           enddo
           nnn=1
        else
           nnn=nn
        endif
        ecc=sqrt(dot(e,e))
        a=hh**2/(1+q)/(1-ecc**2)

        do i=1,3
           n(i)=n(i)/nnn
           hhat(i)=h(i)/hh
           eh(i)=e(i)/ecc
        enddo

        call cross(hhat,n,hn)
        call cross(hhat,eh,he)

        cosOm=dot(ii,n)
        sinOm=dot(jj,n)
        cosi=dot(kk,hhat)
        sini=nn
        cosw=dot(n,eh)
        sinw=dot(eh,hn)

        element(1)=a
        element(2)=ecc
c       element(3)=datan2(sini,cosi)            1/3/00
        element(3)=dacos(cosi)
        element(4)=datan2(sinw,cosw)
        element(5)=datan2(sinOm,cosOm)

        end



        subroutine transform5(q,r,v,e,h)
*
*       Calculates vectors e and h given r,v and q.
*       -------------------------------------------
*
        implicit real*8 (a-h,m,o-z)

        real*8 r(3),v(3)
        real*8 h(3),e(3),dot

        call cross(r,v,h)
        call cross(v,h,e)

        rr=sqrt(dot(r,r))

        do i=1,3
           e(i)=e(i)/(1+q)-r(i)/rr
        enddo

        end



        subroutine transform6(e,h,kk,q,element)
*
*       Calculates 5 orbital elements given vectors evec and h and mass ratio q.
*       Uses Hout for kk
*       Precession Om not corrected
*       -----------------------------------------------------------------------
*
        implicit real*8 (a-h,m,o-z)

        real*8 element(5)
        real*8 h(3),e(3),n(3),hn(3),he(3)
        real*8 eh(3),dot
        real*8 ii(3),jj(3),kk(3),nn,nnn
        data ii/1.d0,0.d0,0.d0/
        data jj/0.d0,1.d0,0.d0/

        hh=sqrt(dot(h,h))

        call cross(kk,h,n)

        nn=sqrt(dot(n,n))
        if(nn.eq.0)then
           nnn=1
        else
           nnn=nn
        endif
        ecc=sqrt(dot(e,e))
        a=hh**2/(1+q)/(1-ecc**2)

        do i=1,3
           n(i)=n(i)/nnn
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
c       element(3)=atan2(sini,cosi)
        element(3)=acos(cosi)
        element(4)=atan2(sinw,cosw)
        element(5)=atan2(sinOm,cosOm)

        end



        subroutine transform7(q,r,v,ehat,qhat,hhat)
*
*       Calculates basis vectors ehat, qhat and hhat given r,v and q.
*       ------------------------------------------------------------
*
        implicit real*8 (a-h,m,o-z)

        real*8 r(3),v(3)
        real*8 h(3),e(3),dot
        real*8 ehat(3),qhat(3),hhat(3)

        call cross(r,v,h)
        call cross(v,h,e)

        rr=sqrt(dot(r,r))

        ee=0
        hh=0

        do i=1,3
           e(i)=e(i)/(1+q)-r(i)/rr
           ee=ee+e(i)**2
           hh=hh+h(i)**2
        enddo

        do i=1,3
           ehat(i)=e(i)/sqrt(ee)
           hhat(i)=h(i)/sqrt(hh)
        enddo

        call cross(hhat,ehat,qhat)
        
        end



