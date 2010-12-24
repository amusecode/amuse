*
*
      subroutine intcol(nup)
*
*      calculate collisions between stars 
*      ----------------------------------
*
*
      include 'common.h'
*
      real*8 dt,den,sm1,sm2,r1,r2,pmax,tim,ptold,ptnew,sm12,dmass,
     &       ehpot,ehkin,cor,vrol1,vrol2,vtol1,vtol2,tlog,smlog,
     &       smsss,sturn,z,zpars(20),co,zkino,zkinn,ehhhk,ehhhc,
     &       sm1s,sm2s,smtold,dmasst,ek1,ek2,ek12,eku,ekm,ekto,ektn,
     &       vvr,vvt,tform,tb3f,tb3b3,tcoll,rmax,ub,xlum,zzz,ekicks1,
     &       ekicks2,sm1so,sm2so,ptoldd,zkinoo,etotn,ek12n,ekton,
     &       ehmlevo,ehmlevn,ekickto,ekicktn,slum,zkold,ehbin3o
*
      real*4 ran2
*
      integer l,nup,lmin,lmax,i,j,im1,im2,n,is1,is2,ik1,ik2,in,icoll,
     &        iki,ii,ip,iblue,iexi,iexi1,noicoll
*
      common /timset/ tform(20),tb3f(20),tb3b3(20),tcoll(20)
*
      open(43,file='starm.dat',access='append')
*           
      n = nt
      icoll = 0
      noicoll = 0
      ehhhk = 0.d0
      ehhhc = 0.d0
      zkino = 0.d0
      zkinn = 0.d0
      dmasst = 0.d0
      ekicks1 = 0.d0
      ekicks2 = 0.d0
      smtold = smt*zmbar
      co = float(nt)/log(gamma*nt)
*
      do 10 l = 1,nup
*
         if(l.eq.1) then
           lmin = 1
         else
           lmin = nzst(l-1) + 1
         endif
*
         lmax = nzst(l)
*
         if(ltwo(l).ne.0) then
           dt = tau*2.0d0**(-ltwo(l))
         else
           dt = tau
         endif
*
         tcoll(l) = tcoll(l) + dt*tscale0*log(gamma*nt00)/
     &              log(gamma*ntnew)*float(ntnew)/float(nt00)
         tim = tcoll(l)
         print*,'intcol- l,nup,tcoll(l) = ',l,nup,tcoll(l)
         tlog = log10(1.d6*tim)
         if(tlog.lt.7.4) then
           smlog = 0.54d0*tlog**2 - 8.47*tlog + 33.95d0
         else
           smlog = 3.26d0 - tlog/3.d0
         endif
         smsss = 10.d0**smlog
         sturn = smsss
         call getZ0(1,z)
         call zcnsts(z,zpars)
         call mturn(sturn,tim,zpars)
         print*,'col- smsss,sturn,z0,tim=',smsss,sturn,z,tim
*
*       find single stars to colide
*
*       find first star to colide
*
         print*,'lmin,lmax = ',lmin,lmax
         i = lmin - 1
 20      continue
*
         i = i + 1
 
         if(i.ge.lmax) go to 10
*
         im1 = iname(i)
*
         if(abs(ikind(im1)).eq.2) go to 20
         if(abs(body(im1)).lt.tolm) go to 20
         if(nkick(im1).ge.1) go to 20
         if(r(i).gt.1.d8) go to 20
         is1 = names(im1)
         call getLum(is1,slum)
         if(isNaN(slum)) then
           print*,'i,im1,is1,lum = ',i,im1,is1,slum
           go to 20
         endif
*         
         ip = i
 25      call densi(ip,den)
         if(abs(den).lt.1.d-10) then
           print*,'--- den < 0   i,den = ',ip,den
           call flush(6)
           if(ip.gt.5) then 
             ip = ip -1
             go to 25
           else
             ip = ip + 1
             go to 25
           endif         
         endif
*
*     find second star to colide
*
         j = i
 30      continue
*
         j = j + 1
*
         if(j.gt.lmax) go to 10
*
         im2 = iname(j)
*
         if(abs(ikind(im2)).eq.2) go to 30
         if(abs(body(im2)).lt.tolm) go to 30
         if(nkick(im2).ge.1) go to 30
         if(r(j).gt.1.d8) go to 30
         is2 = names(im2)
         call getLum(is2,slum)
         if(isNaN(slum)) then
           print*,'j,im2,is2,lum = ',j,im2,is2,slum
           go to 30
         endif
*
*    compute encounter probability using rmax = R1 + R2
*    save radial and tangential velocities before collision
* 
         vrol1 = vr(im1)
         vrol2 = vr(im2)
         vtol1 = vt(im1)
         vtol2 = vt(im2)
*

         sm1 = body(im1)
         sm2 = body(im2)
         sm12 = sm1 + sm2
*                           
         ek1 = 0.5d0*sm1*(vr(im1)**2 + vt(im1)**2)
         ek2 = 0.5d0*sm2*(vr(im2)**2 + vt(im2)**2)
         ekto = ek1 + ek2
*                                 
         ii = -i         
         call relvel(ii,j,ub)
*
         is1 = names(im1)
         is2 = names(im2)
         call get_ss_type(is1,ik1)
         call get_ss_type(is2,ik2)
         call get_radius(is1,r1)
         call get_radius(is2,r2)
         call get_mass(is1,sm1so)         
         call get_mass(is2,sm2so)
*
*        R1 ad R2 are in solar units
*
*         Changed DCH to account for gravitational focusing   
*
         rmax = (r1 + r2)*rsuntopc*rtidkg/rbar
         pmax = sqrt(rmax**2 + 2*sm12*rmax/ub**2)
*
*          number density n=den/2/pi
*
         pcoll = 0.5d0*co*pmax**2*den*ub*dt
*
*Added DCH 1/8/6 to deal with the case where one star is type 15 (no remnant)
*
         if (ik1.eq.15.or.ik2.eq.15) pcoll = 0.d0
*
c      print*,'i,j,im1,im2,is1,is2,ik1,ik2,sm1,sm2,ub,pm,de,dt,pcoll=',
c     & i,j,im1,im2,is1,is2,ik1,ik2,sm1,sm2,ub,pmax,den,dt,pcoll
c         call flush(6)
*
         if(ran2(irun).le.pcoll) then
*
*   evolve two single stars to the time time - time of collision
*
           vvr = vr(im1)
           vvt = vt(im1)
           vr(im1) = vrol1
           vt(im1) = vtol1
*
           if(icoll.eq.0.and.noicoll.eq.0) then
             call energy(2)   
             ptoldd = pot      
             ptold = pot
             zkino = zkin
             zkinoo = zkin     
             ehmlevo = ehmlev
             ekickto = ekickt
             ehbin3o = ehbin3
           endif 
*
           call mloss_single(i,tim,ekicks1,iexi)
*
           call mloss_single(j,tim,ekicks2,iexi1)
*
           if(iexi.eq.0.or.iexi1.eq.0) then
             i = j
             noicoll = noicoll + 1
             print*,'intcoll-no  im1,im2,iexi,iexi1=',im1,im2,iexi,iexi1
             go to 20
           endif
*             
c           tim = tim + 2.d-14
*           
           call get_mass(is1,sm1s)         
           call get_mass(is2,sm2s)
*
*        compute new realative velocity and center of mass velocity
*        after evolution step
*
           ubo = ub
           ii = -i
           vrol1 = vr(im1)
           vtol1 = vt(im1)         
           call relvel(ii,j,ub)
           print*,'ii ubo,ub,vrol1,vtol1,vr1,vt1=',ii,ubo,ub,vrol1,
     &           vtol1,vr(im1),vt(im1)
           vvr = vr(im1)
           vvt = vt(im1)
           vr(im1) = vrol1
           vt(im1) = vtol1
*     
*         compute total potential before mergre to deal with potential
*         changes due to substituting two single stars by a merger
*
           print*,'intcoll-1  im1,im2,is1,is2,sm1,sm1so,sm1s,sm2,',
     &            'sm2so,sm2s,ekicks1,ekicks1,ekickt,tim = ',
     &            im1,im2,is1,is2,sm1*zmbar,sm1so,sm1s,sm2*zmbar,
     &            sm2so,sm2s,ekicks1,ekicks1,ekickt,tim
           sm1 = sm1s/zmbar
           sm2 = sm2s/zmbar
           sm12 = sm1 + sm2
           ekton = 0.5d0*sm1*(vr(im1)**2 + vt(im1)**2) +
     &             0.5d0*sm2*(vr(im2)**2 + vt(im2)**2)
*         
           vr(im1) = vvr
           vt(im1) = vvt
           ek12 = 0.5d0*sm12*(vr(im1)**2 + vt(im1)**2)
*
*      set the time and compute the collision output
*
           ncoll = ncoll + 1
           icoll = icoll + 1 
*
           in = 2
           print*,'in,is1,is2,sm1,sm2,sm12,tim = ',
     &             in,is1,is2,sm1*zmbar,sm2*zmbar,sm12*zmbar,tim
*           print*,'--- before collStar ---'
           call flush(6)
           call collStars(in,is1,is2,tim)
*           print*,'--- after collStar ---'
           call flush(6)
           call get_ss_updatetime(is1,uptime(im1))
*           print*,'--- update time ---'
           print*,'in,is1,uptim = ',in,is1,uptime(im1)
           call flush(6)
*
*       merger keeps the identity of the first star. The second star becomes
*       a dumy star and it is put outside of the system r(j) = 1.0d8 + im2
*
           r(j) = 1.0d8 + float(im2)
           xmax(im2) = 1.01d0*r(j)
           names(im2) = 0
           body(im2) = 0.d0
           vr(im2) = 0.d0
           vt(im2) = 0.d0
*
*       compute mass loss, heating and coooling due coliisions and 
*       changes of potential
*
           call get_mass(is1,sm1s)
           call get_ss_type(is1,iki)
           print*,'intcoll-2  im1,is1,iki,sm1s,dmass =',im1,is1,iki,
     &            sm1s,(sm12*zmbar - sm1s)
           iblue = 0
           if(ik1.le.1.and.ik2.le.1) then
             if(iki.le.1) then
               if(sm1s.gt.sturn.and.ibstra(im1).eq.0) then  
                 print*,'coll-- smsss,sturn,tim=',smsss,sturn,tim
                 ibstra(im1) = 2
                 ibsc = ibsc + 1
               endif
c               ibsc = ibsc + 1
               iblue = 1
               print*,'b. str. coll. is1,is2,ik1,ik2,iki,ibsc = ',
     &               is1,is2,ik1,ik2,iki,ibsc
               if(iki.gt.1) print*, 'collision  iki > 1  iki = ',iki
             endif
           endif
           sm1s = sm1s/zmbar
           body(im1) = sm1s
*
*  for colliding star ikind is set to 4           
*
           ikind(im1) = 4
*           
           dmass = sm12 - sm1s
           if(dmass.lt.0.d0)  then
             print*,' ---dmass < 0 ---', dmass*zmbar
             if(dmass.lt.-1.d-7) then
               dmass = 0.d0 
               sm1s = sm12*zmbar
               call setMass(is1,sm1s)
               sm1s = sm12
cFor following line see mig 24/7/6
               body(im1) = sm1s
               print*,'sm1s = sm12 ',sm12*zmbar
             endif  
           endif  
           if(abs(dmass).lt.1.d-10) then 
             dmass = 0.d0
           endif 
           dmasst = dmasst + dmass*zmbar
           slosco = slosco + dmass
           ehpot = -0.5d0*dmass*u(i)
c           ehpot = 0.d0
           ehkin = 0.5d0*dmass*(vr(im1)**2 + vt(im1)**2)
           print*,'tim,vrol1,vrol2,vr1,vtol1,vtol2,vt1,ub = ',
     &            tim,vrol1,vrol2,vr(im1),vtol1,vtol2,vt(im1),ub
           cor = -0.5d0*dmass**2/r(i)
c           cor = 0.d0
c           ehcoll = ehcoll + ehpot - ehkin + cor
           ehcoll = ehcoll - ehkin
           ehhhk = ehhhk - ehkin
           ehhhc = ehhhc - 0.5d0*sm1*sm2*ub*ub/sm12
           eku = 0.5d0*sm1*sm2*ub*ub/sm12
           ekm = ehkin
           ektn = ek12 + eku + ekm
           ek12n = 0.5d0*body(im1)*(vr(im1)**2 + vt(im1)**2)
           print*,'ek1,ek2,ekto,ekton,ek12,ek12n,eku,ekm,ektn,Dekt = ',
     &     ek1,ek2,ekto,ekton,ek12,ek12n,eku,ekm,ektn,ektn-ekton
           eccoll = eccoll + 0.5d0*sm1*sm2*ub*ub/sm12
*
           print*,'i,j,im1,im2,is1,is2,iki,sm1,sm2,sm1s,dmass,',
     & 'ehpot,ehkin,eckin,cor = ',i,j,im1,im2,is1,is2,iki,sm1*zmbar,
     & sm2*zmbar,sm1s*zmbar,dmass*zmbar,ehpot,ehkin,eku,cor
           call flush(6)
*
           call get_radius(is1,r1)
           call getLum(is1,xlum)           
*
           print*,'111',is1,im1,tim,sm1s*zmbar,sm1*zmbar,r1,xlum,
     &                  r(i),ik1,iki,ibstra(im1)
           call flush(6)
           write(43,40) is1,im1,tim,sm1s*zmbar,sm1*zmbar,r1,xlum,r(i),
     &                  ik1,iki,ibstra(im1),ikind(im1)
           print*,'222',is2,im2,tim,0.d0,sm2*zmbar,r2,0.d0,r(j),ik2,0,0
           call flush(6)
           write(43,40) is2,im2,tim,0.d0,sm2*zmbar,r2,0.d0,r(j),ik2,
     &                  0,0,0
*
 40        format(1x,'#',2i7,1p6e16.8,4i4)
           call flush(6)
*
           i = j
*
           go to 20
*           
         else
*
*       in the case of merger is not happened return to radial and 
*       tangential velocities before computing collision probability
*
           vr(im1) = vrol1
           vr(im2) = vrol2
           vt(im1) = vtol1
           vt(im2) = vtol2
*
           i = j
*
           go to 20           
*           
         endif
*
 10   continue
*
       close(43)
*       
*      compute new potential and sort stars after removing from the
*      calculation star which is now part of a merger
*
      if(icoll.eq.0.and.noicoll.eq.0) return
*
      call sort3(n,r,iname,inameo)
      nt = nt - icoll
      nzst(nsuzon) = nt
      smt = smt*zmbar - dmasst
      print*,'smtold,smt,dmasst,ntold,nt,icoll,noicoll = ',
     &  smtold,smt,dmasst,nt+icoll,nt,icoll,noicoll 
      smt = smt/zmbar
      call coepot
      call energy(2)
      ptnew = pot
      zkinn = zkin
      ehmlevn = ehmlev
      ekicktn = ekickt
      ehbin3n = ehbin3
      zzz = zkinn - zkino
      enepot = enepot + ptnew - ptold
c      enekin = enekin + zkinn - zkino
c      enekin = enekin - zzz
      etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
     &             ehmlev - ehcoll + eccoll - ekickt - enekin
      write(6,6543) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &               ekickt,enekin,escsta,ehbin3,ehb3b3
 6543 format(1x,'ontcol-f etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)   
*
      write(6,5432) zkino,zkinn,ehmlevo,ehmlevn,ehbin3o,ehbin3n,      
     &              ekickto,ekicktn,zkinn-zkino,ehmlevn-ehmlevo,
     &              ehbin3n-ehbin3o,ekicktn-ekickto,(zkinn-zkino)-
     &              (ehmlevn-ehmlevo+ehbin3n-ehbin3o+ekicktn-ekickto)
 5432 format(1x,'coll  zkold,zknew,ehmlevo,ehmlevn,ekickto,ekicktn',
     &',ehbin3o,ehbin3n,Dzk,Dhm,Dhb3,Dhki,DDcoll =',1p13e14.6)
      print*,'ico,nco,ptold,ptnew,enept=',icoll,ncoll,ptold,ptnew,enepot
      write(6,8721) zkino,zkinn,zzz,ehhhk,-ehhhc,ehhhk-ehhhc,
     &           ehmlevn-ehmlevo,zzz-(ehhhk-ehhhc)-(ehmlevn-ehmlevo)  
 8721 format(1x,'zkino,zkinn,Dzkin,ehhhk,ehhhc,Dehhh,Dehm,DDcol = ',
     &1p8e14.6)
*
      return
      end
*
*
*
   
