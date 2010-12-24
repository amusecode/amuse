*
      subroutine intb3b3b(nup)
*
*       calculate interaction between two three-body binaries according to
*       ------------------------------------------------------------------
*       Sdodolkiewicz's prescription (Acta Astronomica, 1986, 36, 19)
*       -------------------------------------------------------------
*       ib3f = 3 - Pmax used to calculate encounter probability
*       -------------------------------------------------------
*       ib3f = 4 - Pmax used to calculate encounter probability and
*       -----------------------------------------------------------
*                  numerical integration of four-body interaction
*                  ----------------------------------------------
*
*
      include 'common.h'
*
      real*8 rlmax,rbmin,rbmax,ek,pt,ap,sm1234,sig,ub,deb,deb1,
     &       deb2,deb3,deb4,dt,vek,vtek,vrek,rlmin,pbmin,pbmax,
     &       co,ptd,pb3b3d,denb(nmax),sm1,sm2,xxmt,
     &       xxmiu,vc2,vvc,pmax,denbav,timevo,tlog,rad1,rad2,smlog,
     &       smsss,sturn,z,zpars(20),deb3s,pppp,rpmax,sq,ssme,sm3,
     &       sm4,smsx,sm1s,sm2s,sm3s,xm1,xm2,x23,x13,x23x123,x3x13,
     &       x123x12,x13x123,xcof1,xcof,xform,pexch1,pexch2,yy,
     &       pexch,aold,ebold,eold,age1,ages,sm12o,timea,dmass,
     &       cor,ehpot,ehkin,dmasst,anew,sm12oo,sm3oo,den,sm3a,
     &       sm4a,a00,a01,a02,a03,a10,a11,a12,a20,a21,a30,tform,tb3f,
     &       tb3b3,tcoll,c1,anb,cfr,enew,sm12,sm123,sm12n,sm3o,sy1,
     &       sy2,xrun,pctonbu,ptold,ptnew,zkold,zknew,etotn,apnew
     &       apold,binold,rlo,rln,ss1234,ehb3b3o,ehb3b3n,ehmlevo,
     &       ehmlevn,tphys,sm1234o,sm1234n,ssm12o,ssm12n,ekickto,
     &       ekicktn,ehkin1,ssevo,sb1,sb2,sm1d,sm2d,sb11,sb12,sb21,
     &       sb22,sm1i,sm2i,sm3i,sm4i,sm1f,sm2f,sm3f,sm4f,rad1i,rad2i,
     &       rad3i,rad4i,xlum1i,xlum2i,xlum3i,xlum4i,ecc1i,ecc2i,ap1i,
     &       ap2i,dmm
*
      real*4 ran2
*
      integer l,nup,lmin,lmax,n,im1,im2,k,kl,km,kn,inum,ndist,idb,
     &        inumt,iexist,id3,id4,nb(nmax),ib1,ib2,im,jh,jhh,ik1,
     &        ik2,id1,id2,kmm,iime,iexch1,ids,id1a,id2a,iwww,nexch,
     &        id1x,imx,id3a,id4a,ixxx,nnn,iexistb1,iexistb2,ib3b3,
     &        noib3b3,ndistm,ntx,ik0n,iob,iobx,nob,nsevo,ik1o,ik2o,
     &        idb1,idb2,id11,id12,id21,id22,ik01,ik02,ik01f,ik02f,
     &        ik11,ik12,ik21,ik22,ik11f,ik12f,ik21f,ik22f,ibbb
*
      common /timset/ tform(20),tb3f(20),tb3b3(20),tcoll(20)
*
      data a00,a01,a02,a03,a10,a11,a12,a20,a21,a30 /3.70d0,7.49d0,
     &     -15.49d0,3.07d0,-1.89d0,-2.93d0,13.15d0,-2.92d0,-5.23d0,
     &     3.12d0/
                
*
*
      if((nbin3-nescb3-ndist3-ndist4-ndiste-nmerge).le.1) return
*
      open(44,file='binarym.dat',access='append')
      open(43,file='starm.dat',access='append')
*
      deb3s = 0.d0
      co = float(nt)/log(gamma*nt)
      c1 = co
      dmasst = 0.d0
      pctonbu = rtidkg/rbar
      n = nt
      ndist = 0
      ndistm = 0
      nexch = 0
      ssme = 0.d0
      ssevo = 0.d0
      ehkin1 = 0.d0
      nsevo = 0
      iime = 0
      ib3b3 = 0
      iob = 0
      noib3b3 = 0
      iexistb1 = 1
      iexistb2 = 1
*
*       probability of binary-binary interaction according to Pmax
*
*       calculate binary number density up to nzon(nup)
*
      call denbin(nup,nb,denb)
*
c      nnn = nzst(nup)
c      do 99 jh = 1,nnn
c      im = iname(jh)
c      if(ikind(im).ne.2) go to 99
c      jhh = nb(jh)
c      if (iprint.eq.0)
c     &     print*,'jh,jhh,nb(jh),denb(jhh) = ', jh,jhh,nb(jh),denb(jhh)
c 99   continue
*
      do 100 l=1,nup
*
*
*       calculate binary number density up to nzon(nup)
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
*     calculate actual time
*
         tb3b3(l) = tb3b3(l) + dt*tscale0*log(gamma*nt00)
     &              /log(gamma*ntnew)*float(ntnew)/float(nt00)
         timevo = tb3b3(l)
         print*,'intb3b3- l,nup,tb3b3(l) = ',l,nup,tb3b3(l)
         tlog = log10(1.d6*timevo)
*
         if(tlog.lt.7.4) then
           smlog = 0.54d0*tlog**2 - 8.47*tlog + 33.95d0
         else
           smlog = 3.26d0 - tlog/3.d0
         endif
         smsss = 10.d0**smlog
         sturn = smsss
         call getZ0(1,z)
         call zcnsts(z,zpars)
         call mturn(sturn,timevo,zpars)
         print*,'b3b3  smsss,sturn,z0,time=',smsss,sturn,z,timevo
         call flush(6)
*
*      find binaries to interact
*
         km = lmin - 1
*
 110     continue
*
         km = km + 1
         if(km.gt.lmax) go to 100
*
         im1 = iname(km)
         if(ikind(im1).ne.2) go to 110
         k = nwhich(im1)
*
*   do not check binaries which escaped due to relaxation or interactions or
*   were dissolved. Binaries which take part in interactions in this
*   time steap have naegative ikind. They can still take part with
*   interaction with other binaries. So ikind() ---> abs(ikins())
*
         if(abs(body(im1)).lt.tolm) go to 110
         if(abs(bin(k,1)).lt.tolm.or.
     &     abs(bin(k,2)).lt.tolm) go to 110
         if(bin(k,4).lt.0.0d0) go to 110
         if(bin(k,6).lt.0.0d0) go to 110
         if(bin(k,7).lt.0.0d0) go to 110
         if(r(km).gt.1.d8) go to 110
*
         bin(k,8) = r(km)
*
         kn = km
*
 120     continue
         kn = kn + 1
         if(kn.gt.lmax) go to 100
*
         im2 = iname(kn)
         if(ikind(im2).ne.2) go to 120            
         kl = nwhich(im2)
*
*   do not check binaries which escaped due to relaxation or interactions or
*   were dissolved. Binaries which take part in interactions in this
*   time steap have naegative ikind. They can still take part with
*   interaction with other binaries. So ikind() ---> abs(ikins())
*
         if(abs(body(im2)).lt.tolm) go to 120
         if(abs(bin(kl,1)).lt.tolm.or.
     &     abs(bin(kl,2)).lt.tolm) go to 120
         if(bin(kl,4).lt.0.0d0) go to 120
         if(bin(kl,6).lt.0.0d0) go to 120
         if(bin(kl,7).lt.0.0d0) go to 120
         if(r(kn).gt.1.d8) go to 120
*
         bin(kl,8) = r(kn)
*
*        Compute encounter probability
*
         call relvel(km,kn,ub)
         call densi(km,den)
*
*     calculate PMAX using formula given by Hut and Bahcall ApJ 268, 1983         
*     and Bacon and Sigurdsson astro-ph96/03036
*
         xxmt = body(im1) + body(im2)
         xxmiu = (body(im1)*body(im2))/xxmt
         vc2 = 2.D0/xxmiu*(bin(k,5) + bin(kl,5))
         vvc = ub/dsqrt(vc2)
*
c         if(bin(k,3).gt.bin(kl,3)) then
c           pmax = (5.D0/vvc + 0.6D0*(1.D0 + bin(k,4)))*bin(k,3)       
c           pppp = 15.d0*bin(k,3)
c         else
c           pmax = (5.D0/vvc + 0.6D0*(1.D0 + bin(kl,4)))*bin(kl,3)
c           pppp = 15.d0*bin(kl,3)
c         endif
*
*        limit Pmax - not allow to distant interactions
*
c         if(pmax.gt.pppp) pmax = pppp
*
*         limit Pmax according to maximum pericentre distance equal to
*         R_p = 10.0*a
*
         if(bin(k,3).gt.bin(kl,3)) then
c           rpmax = 3.d0*bin(k,3)
           rpmax = bin(k,3)*(1.d0 + bin(k,4))
           sq = 1.d0 + 2.d0*xxmt/ub/ub/rpmax
           pmax = rpmax*sqrt(sq)
         else
c           rpmax = 3.d0*bin(kl,3)
           rpmax = bin(kl,3)*(1.d0 + bin(kl,4))
           sq = 1.d0 + 2.d0*xxmt/ub/ub/rpmax
           pmax = rpmax*sqrt(sq)
         endif
*         
*           calculate average binary density
*
         ib1 = nb(km)
         ib2 = nb(kn)
         if(ib1.eq.0.or.ib2.eq.0) then
           print*,'im1,im2,km,kn,ib1,ib2=',im1,im2,km,kn,ib1,ib2
           call flush(6)
         endif
         denbav = 0.5d0*(denb(ib1) + denb(ib2))
         pb3b3d = pi*co*pmax**2*denbav*ub*dt
*
         kmm = km
         km = kn
*         
         if(pb3b3d.lt.ran2(irun)) go to 110
*
         print*,'nup,lmax,km,kn,k,kl,im1,im2=',nup,lmax,km,kn,
     &           k,kl,im1,im2
              
         write(6,4532) timevo,denb(ib1),denb(ib2),pmax,ub,dt,
     &   vc2,vvc,bin(k,3),bin(kl,3),bin(k,4),bin(k,4),bin(k,5),
     &   bin(kl,5),co,body(im1)*zmbar,bin(k,1)*zmbar,bin(k,2)*zmbar,
     &   body(im2)*zmbar,bin(kl,1)*zmbar,bin(kl,2)*zmbar,gamma,pb3b3d,
     &   nt,nt0,im1,im2,ib1,ib2,kmm,kn,k,kl
 4532    format(1x,'timevo,denb(ib1),denb(ib2),pmax,ub,dt,vc2,',
     &'vvc,bin(k,3),bin(kl,3),bin(k,4),bin(k,4),bin(k,5),bin(kl,5),co,',
     &'body(im1),bin(k,1),bin(k,2),body(im2),bin(kl,1),bin(kl,2),gamma',
     &'pb3b3d,nt,nt0,im1,im2,ib1,ib2,km,kn,k,kl = ', 1p23e12.4,10i8)    
         sm1234 = xxmt
*
*   evolve binary and single star to the time timevo - time of binary
*   single star interaction
*
        if(ib3b3.eq.0.and.noib3b3.eq.0) then
          call energy(2)   
          ptold = pot      
          zkold = zkin     
          ehmlevo = ehmlev
          ehb3b3o = ehb3b3
          ekickto = ekickt
        endif 
*
        ib3b3 = ib3b3 + 1
        ntx = nt
        iobx = iobt
*
        ib1 = nameb(im1)
        ib2 = nameb(im2)
        call get_loid_mass(ib1,sm1)
        call get_hiid_mass(ib1,sm2)
        call get_loid_mass(ib2,sm3)           
        call get_hiid_mass(ib2,sm4)
        write(6,7654) im1,im2,ib1,ib2,
     &                (body(im1)-bin(k,1)-bin(k,2))*zmbar,
     &                (body(im2)-bin(kl,1)-bin(kl,2))*zmbar,
     &                (body(im1)*zmbar - sm1 - sm2),
     &                (body(im2)*zmbar - sm3 - sm4)
 7654  format(1x,'im1,im2,ib1,ib2,Dsm1,Dsm2,DDsm1,DDsm2 =',
     &        4i9,1p4e14.6)
*                                                                                       *
        call mloss_single(kmm,timevo,ekickb1,iexistb1)
* 
        call mloss_single(kn,timevo,ekickb2,iexistb2)
        call flush(6)
* 
c        tphys = timevo + 2.d-14
        tphys = timevo
        sm1234o = sm1234
        sm1234 = body(im1) + body(im2)
        sm1234n = bin(k,1) + bin(k,2) + bin(kl,1) + bin(kl,2)
*
        ib1 = nameb(im1)
        ib2 = nameb(im2)
        ndistm = ndistm + (nt - ntx) - (iobt - iobx)
*
        if(iexistb1.eq.0.or.iexistb2.eq.0) then
          noib3b3 = noib3b3 + 1
          print*,'intb3b3-no  im1,im2,km,kn,ib1,ib2,iexistb1,iexistb2=',
     &            im1,im2,kmm,kn,ib1,ib2,iexistb1,iexistb2
          go to 110
        endif
*
      print*,'intb3b3-single im1,im2,ib1,ib2,ib3b3.noib3b3,ixb1,ixb2=',
     &          im1,im2,ib1,ib2,ib3b3,noib3b3,iexistb1,iexistb2 
*         
           idb = nameb(im1)
           call get_loid(idb,id1)
           call get_hiid(idb,id2)
           call get_loid_mass(idb,sm1)
           call get_hiid_mass(idb,sm2)
           call get_ss_type(id1,ik1)
           call get_ss_type(id2,ik2)
           call get_bs_type(idb,ik01)
           call getEcc(idb,ecc1i)
           call getSemi(idb,ap1i)
           call getLum(id1,xlum1i)
           call getLum(id2,xlum2i)
           call getRadius(id1,rad1i)
           call getRadius(id2,rad2i)
           ik11 = ik1
           ik12 = ik2
           id11 = id1
           id12 = id2
           idb1 = idb
           sm1i = sm1
           sm2i = sm2
           sb1 = (sm1 + sm2)/zmbar
           print*,'intb3b3-a   km,im1,idb,id1,id2,ik1,ik2,sm1,sm2=',
     &            kmm,im1,idb,id1,id2,ik1,ik2,sm1,sm2
           call flush(6)
*
           idb = nameb(im2)
           call get_loid(idb,id1)
           call get_hiid(idb,id2)
           call get_loid_mass(idb,sm1)  
           call get_hiid_mass(idb,sm2)
           call get_ss_type(id1,ik1)
           call get_ss_type(id2,ik2)
           call get_bs_type(idb,ik02)
           call getEcc(idb,ecc2i)
           call getSemi(idb,ap2i)
           call getLum(id1,xlum3i)
           call getLum(id2,xlum4i)
           call getRadius(id1,rad3i)
           call getRadius(id2,rad4i)
           ik21 = ik1
           ik22 = ik2
           id21 = id1
           id22 = id2
           idb2 = idb
           sm3i = sm1
           sm4i = sm2
           sb2 = (sm1 + sm2)/zmbar
           ss1234 = sb1 + sb2
           print*,'intb3b3-b   kn,im2,idb,id1,id2,ik1,ik2,sm1,sm2=',
     &            kn,im2,idb,id1,id2,ik1,ik2,sm1,sm2
           write(6,6543),im1,im2,sm1234o*zmbar,sm1234*zmbar,
     &                   sm1234n*zmbar,sm1234*zmbar - ss1234*zmbar
 6543      format(1x,'im1,im2,sm1234o,sm1234,sm1234n,dmass =',2i9,
     &            1p4e14.6)
           call flush(6)
* 
*       pick up outcome of binary interaction
*
         if(ran2(irun).le.0.12d0) then
*
*    FLYBY Binaries survive
*
*       compute energy generation
*
           if(bin(k,3).gt.bin(kl,3)) then
             deb = 0.4d0*bin(k,5)
             deb3 = deb
             deb3s = deb3
             deb4 = 0.0d0
             binold = bin(k,5)
             apold = bin(k,3)*rbar/rtidkg/rsuntopc
             bin(k,5) = 1.4d0*bin(k,5)
             bin(k,3) = 0.5d0*bin(k,1)*bin(k,2)/bin(k,5)
             print*,'im1,k,binold,binnew,deb,apold,apnew,sam1,sm2=',im1,
     &               k,binold,bin(k,5),apold,
     &               bin(k,3)*rbar/rtidkg/rsuntopc,deb,
     &               bin(k,1)*zmbar,bin(k,2)*zmbar
*
*     update binary parameters for binary evolution procedures
*
             idb = nameb(im1)
*
             call get_loid(idb,id1)
             call get_hiid(idb,id2)
             call getRadius(id1,rad1)
             call getRadius(id2,rad2)
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             ik1o = ik1
             ik2o = ik2
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             ssm12o = sm1 + sm2
             print*,'before evo- im1,k,idb,id1,id2,ik1,ik2,sm1,sm2 =',
     &             im1,k,idb,id1,id2,ik1,ik2,sm1,sm2
*        semi-major axis in R_sun unit
             ap = bin(k,3)*rbar/rtidkg/rsuntopc
             call setSemi(idb,ap)
             call getRL(id1,rlo)
*
*       synchronize the dynamical binary parameters (after interaction)
*       with the evolutionary status of the binary. Evolution with time
*       step equal to zero
*       
             call evBinary(idb,tphys)             
             call get_loid(idb,id1)   
             call get_hiid(idb,id2)   
             call getRadius(id1,rad1) 
             call getRadius(id2,rad2)
             call getEcc(idb,ecc)
             call getLum(id1,xlum1)
             call getLum(id2,xlum2)
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             call get_bs_type(idb,ik0)
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             sm1f = sm1
             sm2f = sm2
             sm3f = sm3i
             sm4f = sm4i
             rad1i = rad1
             rad2i = rad2
             ecc1i = ecc
             ap1i = ap
             xlum1i = xlum1
             xlum2i = xlum2
             ik11f = ik1
             ik12f = ik2
             ik01f = ik0
             idb1 = idb
             id11 = id1
             id12 = id2
             print*,'after evo- im1,k,idb,id1,id2,ik1,ik2,sm1,sm2,=',
     &              im1,k,idb,id1,id2,ik1,ik2,sm1,sm2
             bin(k,1) = sm1/zmbar
             bin(k,2) = sm2/zmbar
             bin(k,7) = 1.0d8*float(id1) + float(id2)
             body(im1) = bin(k,1) + bin(k,2)
             ssm12n = sm1 + sm2
             call getSemi(idb,apnew)
             call binary_exists(idb,iexist)
             if(iexist.eq.0) then
               rln = -1.0d0  
               print*,'intb3b3ba-merger' 
             else
               call getRL(id1,rln)
             endif
             call get_bs_type(idb,ik0n)
             dmm = ssm12o - ssm12n
             write(6,9876) im1,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,
     &                     dmm
 9876  format(1x,'1f im1,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,dmm='
     &       ,2i9,2i4,1p6e14.6)
*
*       check if semi major axis is greater than the sum of star radii
*       if not then the merger is due because of 4-body interaction.
*       Here only flag is set to distinguish between mergers due to binary
*       evolution and due to 4-body interaction.
*
             if(ik1o.le.1.and.ik2o.le.1) then
               if(rad1+rad2.gt.ap.or.iexist.eq.0) then
                 if(sm1+sm2.gt.sturn.and.(ik1.le.1.or.ik2.le.1)) then
                   print*,'33-1 sm12,sturn,t=',sm1+sm2,sturn,timevo
                   call flush(6)
                   ibstra(im1) = 4
                   ibs4 = ibs4 + 1
      print*,'b3b3-1  k,im1,idb,sm1,sm2,smsss,a,rad1,rad2,bin,ap,ibs4=',
     &     k,im1,idb,sm1,sm2,smsss,bin(k,3),rad1,rad2,bin(k,5),ap,ibs4
                   call flush(6)
                 endif
               endif
             endif
*
*        update time for the next evolution step
             call get_bs_updatetime(idb,uptime(im1))
*        check for a merger event
             nob = 0
             if(iexist.eq.0) then
               ap1i = 0.d0
               ecc1i = 0.d0
               ikind(im1) = 3
               names(im1) = idb
               nameb(im1) = 0
               call get_loid_mass(idb,sm1)
               call get_hiid_mass(idb,sm2)
               if(ik2.eq.15) then
                 bin(k,2) = 0.d0
                 bin(k,1) = sm1/zmbar
                 names(im1) = id1
                 sm2 = 0.d0
                 sm2f = 0.d0
               endif
               if(ik1.eq.15) then
                 bin(k,2) = sm2/zmbar
                 bin(k,1) = 0.d0
                 names(im1) = id2
                 sm1 = 0.d0
                 sm1f = 0.d0
               endif
               bin(k,3) = 0.d0
               bin(k,4) = 0.d0
               body(im1) = bin(k,1) + bin(k,2)
               sb1 = body(im1)
               sm1234 = sb1 + sb2
               nmerge = nmerge + 1
               iime = iime + 1
               cfr = body(im1)*zmbar
               ssme = ssme + cfr
          print*,'yyy5 timevo,iime,cfr,ssme= ',timevo,iime,cfr,ssme
          print*,'merger binary interaction intb3b3-1 sm1,sm2,nmerge=',
     &           sm1,sm2,nmerge
*
               if(sm1.le.tolm.and.sm2.le.tolm) then
*
*      mark obliterated binary (massless) to removal from the system
*
                 sm1f = 0.d0
                 sm2f = 0.d0
                 body(im1) = 0.d0
                 sb1 = 0.d0
                 sm1234 = sb1 + sb2
                 bin(k,1) = 0.d0 
                 bin(k,2) = 0.d0 
                 iob = iob + 1   
                 iobt = iobt + 1
                 nob = 1 
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
                 r(kmm) = 1.0d+15 + im1
                 xmax(im1) = 1.01d0*r(kmm)
                 print*,'intb3b3-obla iob,im1,sm1,sm2=',iob,im1,sm1,sm2
               endif
               call flush(6)
             endif
*
             if(dmm.le.1.d-10) then
               dmm = 0.d0
             else
               ssevo = dmm/zmbar
               slosev = slosev + dmm
               nmloev = nmloev + 2   
               nsevo = nsevo + 1     
               ehkin = 0.5d0*ssevo*(vr(im1)**2 + vt(im1)**2)
               ehkin1 = ehkin1 + ehkin
               ehmlev = ehmlev - ehkin
               smt = smt - ssevo
               write(6,9321) im1,nob,nsevo,dmm,ssevo,ehkin
 9321          format(1x,'intb3b3a-evo im1,nob,nsevo,dmm,ssevo,ehkin=',
     &                  3i9,1p3e14.6)
             endif
*
           else
*
             deb = 0.4d0*bin(kl,5)
             deb3 = 0.0d0
             deb4 = deb
             binold = bin(kl,5)
             apold = bin(kl,3)*rbar/rtidkg/rsuntopc
             bin(kl,5) = 1.4d0*bin(kl,5)
             bin(kl,3) = 0.5d0*bin(kl,1)*bin(kl,2)/bin(kl,5)
             print*,'im2,kl,binold,binnew,deb,apold,apnew,sm1,sm2=',im2,
     &              kl,binold,bin(kl,5),apold,
     &               bin(kl,3)*rbar/rtidkg/rsuntopc,deb,
     &               bin(kl,1)*zmbar,bin(kl,2)*zmbar
*
*     update binary parameters for binary evolution procedures
*
             idb = nameb(im2)
*
             call get_loid(idb,id1)
             call get_hiid(idb,id2)
             call getRadius(id1,rad1)
             call getRadius(id2,rad2)
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             ik1o = ik1
             ik2o = ik2
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             ssm12o = sm1 + sm2
             print*,'before evo- im2,kl,idb,id1,id2,ik1,ik2,sm1,sm2 =',
     &               im2,kl,idb,id1,id2,ik1,ik12,sm1,sm2
*        semi-major axis in R_sun unit
             ap = bin(kl,3)*rbar/rtidkg/rsuntopc
             call setSemi(idb,ap)
             call getRL(id1,rlo)
*
*       synchronize the dynamical binary parameters (after interaction)
*       with the evolutionary status of the binary. Evolution with time
*       step equal to zero
*       
             call evBinary(idb,tphys)             
             call get_loid(idb,id1)   
             call get_hiid(idb,id2)   
             call getRadius(id1,rad1) 
             call getRadius(id2,rad2)
             call getEcc(idb,ecc)
             call getLum(id1,xlum1)
             call getLum(id2,xlum2)
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             call get_bs_type(idb,ik0)
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             sm3f = sm1
             sm4f = sm2
             sm1f = sm1i
             sm2f = sm2i
             rad3i = rad1
             rad4i = rad2
             ecc2i = ecc
             ap2i = ap
             xlum3i = xlum1
             xlum4i = xlum2
             ik21f = ik1
             ik22f = ik2
             ik02f = ik0
             idb2 = idb
             id21 = id1
             id22 = id2
             print*,'after evo- im12,kl,idb,id1,id2,ik1,ik2,sm1,sm2,=',
     &              im2,kl,idb,id1,id2,ik1,ik2,sm1,sm2
             bin(kl,1) = sm1/zmbar
             bin(kl,2) = sm2/zmbar
             bin(kl,7) = 1.0d8*float(id1) + float(id2)             
             body(im2) = bin(kl,1) + bin(kl,2)
             ssm12n = sm1 + sm2
             call getSemi(idb,apnew)
             call binary_exists(idb,iexist)
             if(iexist.eq.0) then
               rln = -1.0d0  
               print*,'intb3b3bb-merger' 
             else
               call getRL(id1,rln)
             endif
             call get_bs_type(idb,ik0n)
             dmm = ssm12o - ssm12n
             write(6,9875) im2,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,
     &                     dmm
 9875  format(1x,'2f im2,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,dmm='
     &       ,2i9,2i4,1p6e14.6)
*
*       check if semi major axis is greater than the sum of star radii
*       if not then the merger is due because of 4-body interaction.
*       Here only flag is set to distinguish between mergers due to binary
*       evolution and due to 4-body interaction.
*
             if(ik1o.le.1.and.ik2o.le.1) then
               if(rad1+rad2.gt.ap.or.iexist.eq.0) then
                 if(sm1+sm2.gt.sturn.and.(ik1.le.1.or.ik2.le.1)) then
                   print*,'33-2 sm12,sturn,t=',sm1+sm2,sturn,timevo
                   call flush(6)
                   ibstra(im2) = 4
                   ibs4 = ibs4 + 1
      print*,'b3b3-2  k,im1,idb,sm1,sm2,smsss,a,rad1,rad2,bin,ap,ibs4=',
     &     k,im1,idb,sm1,sm2,smsss,bin(k,3),rad1,rad2,bin(k,5),ap,ibs4
                   call flush(6)
                 endif
               endif
             endif
*
*        update time for the next evolution step
*
             call get_bs_updatetime(idb,uptime(im1))
*
*        check for a merger event
*
             nob = 0
             if(iexist.eq.0) then
               ap2i = 0.d0
               ecc2i = 0.d0
               ikind(im2) = 3
               names(im2) = idb
               nameb(im2) = 0
               call get_loid_mass(idb,sm1)
               call get_hiid_mass(idb,sm2)
               if(ik2.eq.15) then
                 bin(kl,2) = 0.d0
                 bin(kl,1) = sm1/zmbar
                 names(im2) = id1
                 sm2 = 0.d0
                 sm4f = 0.d0
               endif
               if(ik1.eq.15) then
                 bin(kl,2) = sm2/zmbar
                 bin(kl,1) = 0.d0
                 names(im2) = id2
                 sm1 = 0.d0
                 sm3f = 0.d0
               endif
               bin(kl,3) = 0.d0
               bin(kl,4) = 0.d0
               body(im2) = bin(kl,1) + bin(kl,2)
               sb2 = body(im2)
               sm1234 = sb1 + sb2
               nmerge = nmerge + 1
               iime = iime + 1
               cfr = body(im2)*zmbar
               ssme = ssme + cfr
          print*,'yyy6 timevo,iime,cfr,ssme= ',timevo,iime,cfr,ssme
          print*,'merger binary interaction intb3b3-2 sm1,sm2,nmerge=',
     &            sm1,sm2,nmerge
*
               if(sm1.le.tolm.and.sm2.le.tolm) then
*
*      mark obliterated binary (massless) to removal from the system
*
                 sm3f = 0.d0
                 sm4f = 0.d0
                 body(im2) = 0.d0
                 sb2 = 0.d0
                 sm1234 = sb1 + sb2
                 bin(kl,1) = 0.d0 
                 bin(kl,2) = 0.d0 
                 iob = iob + 1   
                 iobt = iobt + 1
                 nob = 1 
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
                 r(kn) = 1.0d+15 + im2
                 xmax(im2) = 1.01d0*r(kn)
                 print*,'intb3b3-oblb iob,im2,sm1,sm2=',iob,im1,sm1,sm2
               endif
                 call flush(6)
             endif
*
             if(dmm.le.1.d-10) then
               dmm = 0.d0
             else
               ssevo = dmm/zmbar
               slosev = slosev + dmm
               nmloev = nmloev + 2   
               nsevo = nsevo + 1     
               ehkin = 0.5d0*ssevo*(vr(im2 )**2 + vt(im2)**2)
               ehkin1 = ehkin1 + ehkin
               ehmlev = ehmlev - ehkin
               smt = smt - ssevo
               write(6,9322) im2,nob,nsevo,dmm,ssevo,ehkin
 9322          format(1x,'intb3b3b-evo im2,nob,nsevo,dmm,ssevo,ehkin=',
     &                  3i9,1p3e14.6)
               call flush(6)
             endif
           endif
*
           deb1 = body(im2)*deb/sm1234
           deb2 = body(im1)*deb/sm1234
*
*   set negative ikind for escape procedure. Negative ikind marks that
*   star/binary was taken place in strong interaction
*
           ikind(im1) = -ikind(im1)
           ikind(im2) = -ikind(im2)
*
*       calculate kick velocities
*
           write(6,5467) im1,im2,body(im1),body(im2),sm1234,deb,
     &                   deb1,deb2
 5467      format(1x,'b3b3<0.12-1 im1,im2,body(im1),body(im2),',
     &'sm1234,deb,deb1,deb2 =',2i9,1p6e14.6)
           call kick(kmm,deb1)
           call kick(kn,deb2)
           write(6,*) 'b3b3<0.12 km,kn,im1,im2,deb,deb1,deb2=',
     &                 kmm,kn,im1,im2,deb,deb1,deb2
           call flush(6)
*
           nb3b3 = nb3b3 + 1
           ehb3b3 = ehb3b3 + deb
*
*      store information about binary binary interaction - first binary
*
           iinte3(k) = iinte3(k) + 1
           inum = iinte3(k)
*
           if(inum.gt.50) then
             iinte3(k) = 1
             inum = 1
           endif
*
           inumt = 50*(k - 1) + inum
           binin(inumt,1) = r(kmm)/pctonbu
           binin(inumt,2) = deb3s
           binin(inumt,3) = body(im1)*zmbar
           binin(inumt,4) = body(im2)*zmbar
           binin(inumt,5) = timevo
           binin(inumt,6) = bin(k,5)
           if(iexist.eq.0) then
             bin(k,5) = 0.d0
             binin(inumt,7) = 5
           else
             binin(inumt,7) = 2
           endif
*
           ibbb = binin(inumt,7)
           write(44,420) idb1,id11,id21,im1,timevo,sm1i,sm1f,sm2i,
     &                 sm2f,rad1i,rad2i,xlum1i,xlum2i,ap1i,ecc1i,r(kmm),
     &                 ik01,ik01f,ik11,ik11f,ik12,ik12f,iexist,
     &                 ibstra(im1),ibbb
*
 420            format(1x,'#',4i7,1p12e16.8,9i4)
*
*      store information about binary binary interaction - second binary
*
           iinte3(kl) = iinte3(kl) + 1
           inum = iinte3(kl)
*
           if(inum.gt.50) then
             iinte3(kl) = 1
             inum = 1
           endif
*
           inumt = 50*(kl - 1) + inum
           binin(inumt,1) = r(kn)/pctonbu
           binin(inumt,2) = deb4
           binin(inumt,3) = body(im1)*zmbar
           binin(inumt,4) = body(im2)*zmbar
           binin(inumt,5) = timevo
           binin(inumt,6) = bin(kl,5)
           if(iexist.eq.0) then
             bin(kl,5) = 0.d0
             binin(inumt,7) = 5
           else
             binin(inumt,7) = 2
           endif
*
           ibbb = binin(inumt,7)
           write(44,420) idb2,id21,id22,im2,timevo,sm3i,sm3f,sm4i,
     &                 sm4f,rad3i,rad4i,xlum3i,xlum4i,ap2i,ecc2i,r(kn),
     &                 ik02,ik02f,ik21,ik21f,ik22,ik22f,iexist,
     &                 ibstra(im2),ibbb
                          
*
         else
*
*       BINARY DISRUPTION
*       
*      compute energy generation
*
           deb = 0.516d0*(bin(k,5) + bin(kl,5))
*
           nb3b3 = nb3b3 + 1
           ehb3b3 = ehb3b3 + deb
           ndist = ndist + 1
*
*      compute potential energy to deal with potential changes due to
*      binary disruption
*
*      Modyfication is needed in the case of primordial binaries. Binary 
*      dynamicaly formed consists of two single stars which were already in 
*      the system. In the case of primordial binaries, binary dissolution 
*      create a new single star. The primary star is put into the binary cm.
*
           nt0 = nt0 + 1
           nt = nt + 1
*
*    FIRST BINARY IS DISRUPTED
*    
           if(bin(k,3).gt.bin(kl,3)) then
             deb1 = body(im2)*deb/sm1234
             deb2 = body(im1)*deb/sm1234
cChanged DCH 2/8/6
             if (body(im1).gt.0.d0) then
                deb3 = bin(k,2)*deb1/body(im1)
                deb4 = bin(k,1)*deb1/body(im1)
             else
                deb3 = 0.d0
                deb4 = 0.d0
             endif
*
             idb = nameb(im1)
             call get_loid(idb,id3)
             call get_hiid(idb,id4)
             call get_loid_mass(idb,sy1)
             call get_hiid_mass(idb,sy2)
             call getRadius(id3,rad1i)
             call getRadius(id4,rad2i)
             call getLum(id3,xlum1i)
             call getLum(id4,xlum2i)
             call get_ss_type(id3,ik11f)
             call get_ss_type(id4,ik21f)
             call get_bs_type(idb,ik01f)
             sm1f = sy1
             sm2f = sy2
             sm3f = sm3i
             sm4f = sm4i
             ap1i = -1000.d0
             ecc1i = -1000.d0
             idb1 = idb
             id11 = id3
             id12 = id4
             sm3 = sy1
             sm4 = sy2
*                                             
             id3a = bin(k,7)/1.0d8
             write(6,2445) im1,k,bin(k,7),id3,id3a
 2445        format(1x,'first binary disrupted: im1,k,bin,id3,id3a =',
     &              2i7,1pe22.14,2i7)
             id4a = bin(k,7) - 1.0d8*id3a
             write(6,2446) id4,id4a
 2446        format(1x,'id4,id4a =',2i7)
             call flush(6)
*
             iname(nt0) = nt0
             ikind(im1) = -1
             ikind(nt0) = -1
             r(nt0) = r(kmm)*1.0000000001d0
             xmax(nt0) = 1.01d0*r(nt0)
             vr(nt0) = -vr(im1)
             vt(nt0) = -vt(im1)
*
*       set the vector ibstra for binary components (now single stars)
*
             ibstra(im1) = 0
             ibstra(nt0) = 0
*
*               remove binary from the list nameb - zero means not binary
*               add star nt0 with id id3 to the names list 
*
             nameb(im1) = 0
             names(im1) = id3
             names(nt0) = id4
             body(im1) = sy1/zmbar
             body(nt0) = sy2/zmbar
             sb11 = body(im1)
             sb12 = body(nt0)
             sb1 = sb11 + sb12
             print*,'k,nt0,im1,names1-2,nameb,id3,id4,id3a,id4a=',
     &             k,nt0,im1,names(im1),names(nt0),nameb(im1),id3,id4,
     &             id3a,id4a
             print*,'sy1,sy2,sm1b,sm2b = ',sy1,sy2,bin(k,1)*zmbar,
     &                   bin(k,2)*zmbar
             call flush(6)
*
             write(43,440) id3,im1,timevo,sy1,sy1,rad1i,xlum1i,r(kmm),
     &                     ik11f,ik11f,ibstra(im1),ikind(im1)
             write(43,440) id4,nt0,timevo,sy2,sy2,rad2i,xlum2i,r(nt0),
     &                     ik21f,ik21f,ibstra(nt0),ikind(nt0)
*
 440             format(1x,'@',2i7,1p6e16.8,4i4)
                       
*
*        compute the new binding energy of survived binary
*
             binold = bin(kl,5)
             apold = bin(kl,3)*rbar/rtidkg/rsuntopc
             bin(kl,5) = bin(kl,5) + deb + bin(k,5)
             bin(kl,3) = 0.5d0*bin(kl,1)*bin(kl,2)/bin(kl,5)
             print*,'im2,kl,binold,binnew,deb,apold,apnew,m1,m2=',im2,
     &               kl,binold,bin(kl,5),deb,apold,
     &               bin(kl,3)*rbar/rtidkg/rsuntopc,bin(kl,1)*zmbar,
     &               bin(kl,2)*zmbar
*
*     update binary parameters for binary evolution procedures
*
             idb = nameb(im2)
             call get_loid(idb,id1)
             call get_hiid(idb,id2)
             call getRadius(id1,rad1)
             call getRadius(id2,rad2)
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             ik1o = ik1
             ik2o = ik2
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             ssm12o = sm1 + sm2
             ebold = bin(kl,5)
             print*,'before diss- im2,kl,idb,id1,id2,ik1,ik2,sm1,sm2 =',
     &               im2,kl,idb,id1,id2,ik1,ik2,sm1,sm2
*        semi-major axis in R_sun unit
             ap = bin(kl,3)*rbar/rtidkg/rsuntopc
             call setSemi(idb,ap)
             call getRL(id1,rlo)
*
*       synchronize the dynamical binary parameters (after interaction)
*       with the evolutionary status of the binary. Evolution with time
*       step equal to zero
*       
             call evBinary(idb,tphys)             
             call get_loid(idb,id1)   
             call get_hiid(idb,id2)   
             call getRadius(id1,rad1) 
             call getRadius(id2,rad2) 
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             print*,'after diss- im2,kl,idb,id1,id2,ik1,ik2,sm1,sm2 =',
     &               im2,kl,idb,id1,id2,ik1,ik2,sm1,sm2
             bin(kl,1) = sm1/zmbar
             bin(kl,2) = sm2/zmbar
             bin(kl,7) = 1.0d8*float(id1) + float(id2)
             ssm12n = sm1 + sm2
             call getSemi(idb,apnew)  
             call binary_exists(idb,iexist)
             if(iexist.eq.0) then
               rln = -1.0d0     
               print*,'intb3b3bc-merger'
             else
               call getRL(id1,rln)
             endif
             call get_bs_type(idb,ik0n)
             dmm = ssm12o - ssm12n
             write(6,9974) im2,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,
     &                     dmm
 9974  format(1x,'1d im2,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,dmm='
     &       ,2i9,2i4,1p6e14.6)
             call flush(6)
*
*       check for exchange interaction
*
*           special treatment for exchange case 
*           (Heggie et al. 1996, ApJ, 467, 359)
*
             if(iexist.eq.0) go to 53
*
             iexch1 = 0
             ixxx = 0
             if(iexch.eq.2.and.iexist.eq.1) then
               if(sm3.ge.sm4) then
                 ids = id3
                 sm3 = sm3/zmbar
                 sm1 = sm1/zmbar
                 sm2 = sm2/zmbar
                 sm12 = sm1 + sm2
                 sm123 = sm12 + sm3
                 ixxx = 1
               else
                 ids = id4
                 sm3a = sm3
                 sm4a = sm4
                 sm3 = sm4/zmbar
                 sm1 = sm1/zmbar
                 sm2 = sm2/zmbar
                 sm12 = sm1 + sm2
                 sm123 = sm12 + sm3
                 ixxx = 2
               endif              
*
               yy = sm123/bin(kl,3)
               if(ub*ub.gt.yy) go to 53
*                  
*      compute probability for exchange
* 
               id1a = bin(kl,7)/1.0d8
               id2a = bin(kl,7) - 1.0d8*id1a
               idb = nameb(im2)
               call get_loid(idb,id1)
               call get_hiid(idb,id2)
               call get_loid_mass(idb,sm1s)
               call get_hiid_mass(idb,sm2s)
               ids = names(im1)
               if(ixxx.eq.2) ids = names(nt0)
               call get_mass(ids,sm3s)
*
      print*,'3b3b- id1a,id1,id2a,id2,ids,sm1,sm1s,sm2,sm2s,am3,sm3s,'
     &,'ub2,yy,ixxx=',id1a,id1,id2a,id2,ids,sm1*zmbar,sm1s,sm2*zmbar,
     &sm2s,sm3*zmbar,sm3s,ub*ub,yy,ixxx                    
*
*     make sure that id of stars and their masses are the same in the
*     evolution procedures (get_loid(idb,id1), get_loid_mass(idb,sm1s)) 
*     and in the bin(k,1), bin(k,2) and bin(k,7)
* 
               if(id1.ne.id1a) then
                 if(sm1s.lt.0.99999d0*sm1*zmbar.or.
     &             sm1s.gt.1.00001*sm1*zmbar) then
                   smsx = sm1
                   sm1 = sm2
                   sm2 = smsx
                   id1x = id1a
                   id1a = id2a
                   id2a = id1x
                   bin(kl,1) = sm1
                   bin(kl,2) = sm2
                   bin(kl,7) = 1.0d8*float(id1a) + float(id2a)
*
                   print*,'3b3b- smsx,sm1,sm2,id1x,id1a,id2a=',
     &                   smsx*zmbar,sm1*zmbar,sm2*zmbar,id1x,id1a,id2a
*
                 endif
               endif
*                      
*       first star is exchanged
* 
               xm1 = sm1/sm12
               xm2 = sm3/sm123
               x23 = sm2 + sm3
               x13 = sm1 + sm3
               x23x123 = x23/sm123
               x3x13 = sm3/x13
               x123x12 = sm123/sm12
               x13x123 = x13/sm123
               xcof1=x23x123**one6*x3x13**3.5*x123x12**one3*x13x123
               xcof = 0.5d0*pi*sm123*bin(kl,3)/ub/ub
               xform = a00 + a01*xm1 + a10*xm2 + a02*xm1**2 +
     &                 a11*xm1*xm2 + a20*xm2**2 + a03*xm1**3 +
     &                 a12*xm1**2*xm2 + a21*xm1*xm2**2 + a30*xm2**3         
*                  
               xform = exp(xform)
               pexch1 = 0.5d0*co*xcof*xcof1*xform*den*ub*dt/pi
               print*,'3b3b- sm1,sm2,sm12,sm3,sm123,ub,a,sm123=',
     &         sm1*zmbar,sm2*zmbar,sm12*zmbar,sm3*zmbar,sm123*zmbar,
     &         ub,bin(kl,3),sm123
               print*,'3b3b- xcof,xcof1,xform,c1,den,ub,dt,iexch1=',
     &                      xcof,xcof1,xform,c1,den,ub,dt,iexch1
*
*       second star is exchanged
*
               xm1 = sm2/sm12
               xm2 = sm3/sm123
               x23 = sm1 + sm3
               x13 = sm2 + sm3
               x23x123 = x23/sm123
               x3x13 = sm3/x13
               x123x12 = sm123/sm12
               x13x123 = x13/sm123
               xcof1=x23x123**one6*x3x13**3.5*x123x12**one3*x13x123
               xcof = 0.5d0*pi*sm123*bin(kl,3)/ub/ub
               xform = a00 + a01*xm1 + a10*xm2 + a02*xm1**2 +
     &                 a11*xm1*xm2 + a20*xm2**2 + a03*xm1**3 +
     &                 a12*xm1**2*xm2 + a21*xm1*xm2**2 + a30*xm2**3
*                    
               xform = exp(xform)
               pexch2 = 0.5d0*co*xcof*xcof1*xform*den*ub*dt/pi
               pexch = pexch1 + pexch2
*
               print*,'3b3b- pexch1,pexch2,pexch=',pexch1,pexch2,pexch
               print*,'3b3b- xcof,xcof1,xform,c1,den,ub,dt = ',xcof,
     &                      xcof1,xform,c1,den,ub,dt
*                                 
               xrun = ran2(irun)
*
               if(xrun.gt.pexch) go to 53
*                    
               iexch1 = 1
*
               nexchang2 = nexchang2 + 1
               nexch = nexch + 1
               print*,'nexchang2,nexch,iexch1',nexchang2,nexch,iexch1
*
               ebold = bin(kl,5)
               aold = bin(kl,3)
               eold = bin(kl,4)
               sm12oo = sm1 + sm2
               sm3oo = sm3
               inexch2(im2) = inexch2(im2) + 1                 
*
               write(6,17) im2,im1,idb,id1,id2,ids,kl,inexch2(im2)
 17            format(1x,'im2,imx,idb,id1,id2,ids,kl,inexch =',9i8)
*
*      check the ages of binary and single star. If they are different
*      synhronize them,
*
               call getAge(id1,age1)
               call getAge(ids,ages)
               if(age1.ne.ages) then
                  print*,'b3b3 wrong ages-1  im1,im2,age1,ages =',
     &                 im1,im2,age1,ages               
               endif
               timea = timevo
*
               sm12 = sm1 + sm2
               sm123 = sm12 + sm3
*                     
               imx = im1
               sm4 = sm4/zmbar
               if(ixxx.eq.2) then
                 imx = nt0
                 sm4 = sm3a/zmbar
               endif 
               xrun = pexch*ran2(irun)
               if(xrun.le.pexch1) then
*
*       the first star is exchanged
*
                 bin(kl,1) = sm3
                 bin(kl,2) = sm2
                 body(im2) = sm3 + sm2
                 body(imx) = sm1
                 sm12n = sm3 + sm2
                 sb2 = sm12n
                 sb11 = body(imx)
                 sb12 = sm4
                 sb1 = sb11 + sb12
                 eold = bin(kl,4)
                 bin(kl,3) = 0.5d0*bin(kl,1)*bin(kl,2)/ebold
                 anew = bin(kl,3)
                 bin(kl,5) = ebold
                 enew = 1.d0 - (sm1/sm3)*(sm123/sm12)**one3
                 if(enew.lt.0.d0) enew = 0.d0
                 bin(kl,4) = enew
                 vr(im2) = sqrt(vr(im2)*vr(im2)*sm12/sm12n)
                 vt(im2) = sqrt(vt(im2)*vt(im2)*sm12/sm12n)
                 vr(imx) = sqrt(vr(imx)*vr(imx)*sm3/sm1)
                 vt(imx) = sqrt(vt(imx)*vt(imx)*sm3/sm1)
                 bin(kl,7) = 1.0d8*float(ids) + float(id2)
                 names(imx) = id1
                 if(names(imx).eq.nameb(im2))then
                   nameb(im2) = ids
                   idb = ids
                 endif
*
                 anb = bin(kl,3)*rbar/rtidkg/rsuntopc
*
       print*,'3b3b- exchange-1:kl,kn,km,im1,im2,idb,id1,id2,ids,names,'
     &,'nameb,sm3,sm2,sm1,sm12o,sm12n,ebold,aold,anew,eold,enew,age1,'
     &,'ages,timea,,anb,anb-NB = ',kl,kn,km,imx,im2,idb,id1,id2,ids,
     &names(imx),nameb(im2),sm3*zmbar,sm2*zmbar,sm1*zmbar,
     &(sm1+sm2)*zmbar,sm12n*zmbar,ebold,aold,anew,eold,enew,age1,
     &ages,timea,anb,anb*rsuntopc*rtidkg/rbar                      
*
*      create a new binary which formed after an exchange interaction
*
                 call create_binary(idb,ids,id2,anb,enew,timevo)
                 call get_loid_mass(idb,sm1d)
                 call get_hiid_mass(idb,sm2d)
                 call binary_exists(idb,iexist)     
                 dmm = sm12n*zmbar - (sm1d + sm2d)
                 write(6,7865) im2,idb,iexist,sm12n*zmbar,
     &                            sm1d + sm2d,dmm
 7865        format(1x,'3b3b-creat1- im2,idb,iexist,sm12o,sm12n,dmm,',
     &        2i9,i4,1p5e14.6)
                 if(iexist.eq.0) print*,'intb3b3-1 iexist=0' 
*                                                
                 go to 53
               endif
*
               if(xrun.gt.pexch1) then
*
*       the second star is exchanged
*
                 bin(kl,1) = sm1
                 bin(kl,2) = sm3
                 body(im2) = sm1 + sm3
                 body(imx) = sm2
                 sm12n = sm1 + sm3
                 sb2 = sm12n
                 sb11 = body(imx)
                 sb12 = sm4
                 sb1 = sb11 + sb12
                 eold = bin(kl,4) 
                 bin(kl,3) = 0.5d0*bin(kl,1)*bin(kl,2)/ebold
                 anew = bin(kl,3)
                 bin(kl,5) = ebold
                 enew = 1.d0 - (sm2/sm3)*(sm123/sm12)**one3
                 if(enew.lt.0.d0) enew = 0.d0
                 bin(kl,4) = enew
                 vr(im2) = sqrt(vr(im2)*vr(im2)*sm12/sm12n)
                 vt(im2) = sqrt(vt(im2)*vt(im2)*sm12/sm12n)
                 vr(imx) = sqrt(vr(imx)*vr(imx)*sm3/sm2)
                 vt(imx) = sqrt(vt(imx)*vt(imx)*sm3/sm2)
                 bin(kl,7) = 1.0d8*float(id1) + float(ids)
                 names(imx) = id2
                 if(names(imx).eq.nameb(im2)) then
                   nameb(im2) = ids
                   idb = ids
                 endif  
*
                 anb = bin(kl,3)*rbar/rtidkg/rsuntopc
*
       print*,'3b3b- exchange-2:kl,kn,km,imx,im2,idb,id1,id2,ids,names,'
     &,'nameb,sm3,sm2,sm1,sm12o,sm12n,ebold,aold,anew,eold,enew,age1,'
     &,'ages,timea,,anb,anb-NB = ',kl,kn,km,im1,im2,idb,id1,id2,ids,
     &names(imx),nameb(im2),sm3*zmbar,sm2*zmbar,sm1*zmbar,
     &(sm1+sm2)*zmbar,sm12n*zmbar,ebold,aold,anew,eold,enew,age1,
     &ages,timea,anb,anb*rsuntopc*rtidkg/rbar
*
*      create a new binary which formed after an exchange interaction
*
                 call create_binary(idb,id1,ids,anb,enew,timea)
                 call get_loid_mass(idb,sm1d)
                 call get_hiid_mass(idb,sm2d)
                 call binary_exists(idb,iexist)
                 dmm = sm12n*zmbar - (sm1d + sm2d)
                 write(6,7864) im2,idb,iexist,sm12n*zmbar,
     &                            sm1d + sm2d,dmm
 7864        format(1x,'3b3b-creat2- im2,idb,iexist,sm12o,sm12n,dmm,',
     &        2i9,i4,1p5e14.6)
                 if(iexist.eq.0) print*,'intb3b3-2 iexist=0'
*
                 go to 53
*
               endif
*
             endif
*
 53          continue
*
*       check if semi major axis is greater than the sum of star radii
*       if not then the merger is due because of 4-body interaction.
*       Here only flag is set to distinguish between mergers due to binary
*       evolution and due to 4-body interaction.
*
             sm1 = bin(kl,1)*zmbar
             sm2 = bin(kl,2)*zmbar
             call get_loid(idb,id1)
             call get_hiid(idb,id2)
             call getRadius(id1,rad1)
             call getRadius(id2,rad2)  
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             call get_loid_mass(idb,sm1d)
             call get_hiid_mass(idb,sm2d)
             call getEcc(idb,ecc2i)
             call getSemi(idb,ap2i)
             call getLum(id1,xlum3i)
             call getLum(id2,xlum4i)
             call get_bs_type(idb,ik02f)
             sm3f = sm1
             sm4f = sm2
             rad3i = rad1
             rad4i = rad2
             ik21f = ik1
             ik22f = ik2
             idb2 = idb
             id21 = id1
             id22 = id2
             if(abs(sm1+sm2-sm1d-sm2d).gt.1.d-10) then
               write(6,5668) im2,idb,id1,id2,ik1,ik2,sm1,sm2,sm1d,sm2d,
     &                        sm1 + sm2 - (sm1d + sm2d)
 5668          format(1x,'wrong-diss-mass-2 im2,idb,id1,id2,ik1,ik2,',
     &'sm1,sm2,sm1d,sm2d,dmmdis=',7i9,1p5e14.6)
             endif                                                                                           
*
             call binary_exists(idb,iexist)
             if(ik1o.le.1.and.ik2o.le.1) then
               if(rad1+rad2.gt.ap.or.iexist.eq.0) then
                 if(sm1+sm2.gt.sturn.and.(ik1.le.1.or.ik2.le.1)) then
                   print*,'33-3 sm12,sturn,t=',sm1+sm2,sturn,timevo
                   call flush(6)
                   ibstra(im2) = 4
                   ibs4 = ibs4 + 1
      print*,'b3b3-3  k,im2,idb,sm1,sm2,smsss,a,rad1,rad2,bin,ap,ibs4=',
     &     k,im2,idb,sm1,sm2,smsss,bin(kl,3),rad1,rad2,bin(kl,5),ap,ibs4
                   call flush(6)
                 endif
               endif
             endif
*
*        update time for the next evolution step
*
             call get_bs_updatetime(idb,uptime(im2))
*
*        check for a merger event
*
             nob = 0
             if(iexist.eq.0) then
               ecc2i = -1000.d0
               ap2i = -1000.d0
               ikind(im2) = 3
               names(im2) = idb
               nameb(im2) = 0
               call get_loid_mass(idb,sm1)
               call get_hiid_mass(idb,sm2)
               if(ik2.eq.15) then
                 bin(kl,2) = 0.d0
                 bin(kl,1) = sm1/zmbar
                 names(im2) = id1
                 sm2 = 0.d0
                 sm4f = 0.d0
               endif
               if(ik1.eq.15) then
                 bin(kl,2) = sm2/zmbar
                 bin(kl,1) = 0.d0
                 names(im2) = id2
                 sm1 = 0.d0
                 sm3f = 0.d0
               endif
               bin(kl,3) = 0.d0
               bin(kl,4) = 0.d0
               body(im2) = bin(kl,1) + bin(kl,2)
               sb2 = body(im2)
               nmerge = nmerge + 1
               sm1234 = sb1 + sb2
               iime = iime + 1
               cfr = body(im2)*zmbar
               ssme = ssme + cfr
            print*,'yyy7 timevo,iime,cfr,ssme= ',timevo,iime,cfr,ssme
            print*,'merger after interaction intb3b3-3 sm1,sm2,nmerge',
     &              sm1,sm2,nmerge
               call flush(6)
*
               if(sm1.le.tolm.and.sm2.le.tolm) then
*
*      mark obliterated binary (massless) to removal from the system
*
                 sm3f = 0.d0
                 sm4f = 0.d0
                 body(im2) = 0.d0
                 sb2 = 0.d0
                 sm1234 = sb1 + sb2
                 bin(kl,1) = 0.d0
                 bin(kl,2) = 0.d0
                 iob = iob + 1  
                 iobt = iobt + 1
                 nob = 1 
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
                 r(kn) = 1.0d+15 + im2
                 xmax(im2) = 1.01d0*r(kn)
                 print*,'intb3b3-oblc iob,im2,sm1,sm2=',iob,im2,sm1,sm2
               endif
               call flush(6)
             endif
*
             if(dmm.le.1.d-10) then
               dmm = 0.d0
             else
               ssevo = dmm/zmbar
               slosev = slosev + dmm 
               nmloev = nmloev + 2   
               nsevo = nsevo + 1     
               ehkin = 0.5d0*ssevo*(vr(im2)**2 + vt(im2)**2)
               ehkin1 = ehkin1 + ehkin
               ehmlev = ehmlev - ehkin
               smt = smt - ssevo
               write(6,9399) im2,nob,nsevo,dmm,ssevo,ehkin
 9399          format(1x,'intb3b3c-evo im2,nob,nsevo,dmm,ssevo,ehkin=',
     &                  3i9,1p3e14.6)
             endif
*
             ikind(im2) = -ikind(im2)
*
             if(sb2.gt.0.d0) then
               deb1 = sb2*deb/sm1234
               deb2 = sb1*deb/sm1234
               deb3 = sb12*deb1/sb1
               deb4 = sb11*deb1/sb1
               call kick(kn,deb2)
               write(6,*) 'b3b3>0.12 kl,kn,im2,deb2,ebin=',
     &                   kl,kn,im2,deb2,bin(kl,5)
               call flush(6)
             else
               deb2 = ebold
               deb3 = sb12*deb/sm1234
               deb4 = sb11*deb/sm1234
             endif
*
*     compute kick velocities
*
              write(6,5468) im1,im2,sb1,sb2,sb11,sb12,sm1234,deb,
     &                      deb1,deb2,deb3,deb4
 5468      format(1x,'b3b3>0.12-1 im1,im2,sb1,sb2,sb11,sb12,sm1234,',  
     &'deb,deb1,deb2,deb3,deb4 =',2i9,1p10e14.6)
             print*,'sb2,deb,deb1-4 b=',sb2,deb,deb1,deb2,deb3,deb4
             call kick(kmm,deb3)
             write(6,*) 'b3b3>0.12 k,km,im1,id3,deb3,ebin=',
     &                   k,kmm,im1,names(im1),deb3,bin(k,5)
             write(6,*) 'nt,nn,nt0,im1,id3,id4 = ',nt,n+ndist+ndistm,
     &                   nt0,im1,id3,id4
             call flush(6)
             call kick(nt0,deb4)
             write(6,*) 'b3b3>0.12 nt0,id4,deb4=',nt0,names(nt0),deb4
             call flush(6)
*
*      store information about binary binary interaction - second binary
*
             iinte3(kl) = iinte3(kl) + 1
             inum = iinte3(kl)
*
             if(inum.gt.50) then
               iinte3(kl) = 1
               inum = 1
             endif
*
             inumt = 50*(kl - 1) + inum
             binin(inumt,1) = r(kn)/pctonbu
             binin(inumt,2) = deb2
             binin(inumt,3) = sb1*zmbar
             binin(inumt,4) = sb2*zmbar
             binin(inumt,5) = timevo
             binin(inumt,6) = bin(kl,5)
             iwww = 0
             if(iexist.eq.0) then
               binin(inumt,7) = 5
               bin(kl,5) = 0.d0
               iwww = 6
               if(iexch1.eq.1) binin(inumt,7) = 9 + iwww
             else  
               binin(inumt,7) = 2
               if(iexch1.eq.1) binin(inumt,7) = 9
             endif
*
             ibbb = binin(inumt,7)
             write(44,420) idb2,id21,id22,im2,timevo,sm3i,sm3f,sm4i,
     &                 sm4f,rad3i,rad4i,xlum3i,xlum4i,ap2i,ecc2i,r(kn),
     &                 ik02,ik02f,ik21,ik21f,ik22,ik22f,iexist,
     &                 ibstra(im2),ibbb
*
*      store information about binary binary interaction - first binary
*
             iinte3(k) = iinte3(k) + 1
             inum = iinte3(k)
*
             if(inum.gt.50) then
               iinte3(k) = 1
               inum = 1
             endif
*
             inumt = 50*(k - 1) + inum
             binin(inumt,1) = r(kmm)/pctonbu
             binin(inumt,2) = -bin(k,5)
             binin(inumt,3) = sb1*zmbar
             binin(inumt,4) = sb2*zmbar
             binin(inumt,5) = timevo
             binin(inumt,6) = bin(k,5)
             binin(inumt,7) = 7
*
             ibbb = binin(inumt,7)
             write(44,420) idb1,id11,id21,im1,timevo,sm1i,sm1f,sm2i,
     &                 sm2f,rad1i,rad2i,xlum1i,xlum2i,ap1i,ecc1i,r(kmm),
     &                 ik01,ik01f,ik11,ik11f,ik12,ik12f,iexist,
     &                 ibstra(im1),ibbb
*
*        update the list of binaries in the system to distinguish
*        between binaries removed due to escape (relaxation - bin(k,6),
*        interaction bin(k,7)) and dissolved bin(k,4)
*
             bin(k,4) = -bin(k,4)
*
           else
*
*     SECOND BINARY IS DISRUPTED
*     
             deb1 = body(im2)*deb/sm1234
             deb2 = body(im1)*deb/sm1234
cChanged DCH 2/8/6
             if(body(im2).ne.0.d0) then
                deb3 = bin(kl,2)*deb2/body(im2)
                deb4 = bin(kl,1)*deb2/body(im2)
             else
                deb3 = 0.d0
                deb4 = 0.d0
             endif
*
             idb = nameb(im2)
             call get_loid(idb,id3)
             call get_hiid(idb,id4)
             call get_loid_mass(idb,sy1)
             call get_hiid_mass(idb,sy2)
             call getRadius(id3,rad3i)
             call getRadius(id4,rad4i)
             call getLum(id3,xlum3i)
             call getLum(id4,xlum4i)
             call get_ss_type(id3,ik21f)
             call get_ss_type(id4,ik22f)
             call get_bs_type(idb,ik02f)
             sm3f = sy1
             sm4f = sy2
             sm1f = sm1i
             sm2f = sm2i
             ap2i = -1000.d0
             ecc2i = -1000.d0
             idb2 = idb
             id21 = id3
             id22 = id4             
             sm3 = sy1
             sm4 = sy2
*                                       
             id3a = bin(kl,7)/1.0d8
             write(6,2447) im2,kl,bin(kl,7),id3,id3a
 2447        format(1x,'second binary disrupted: im2,kl,bin,id3,id3a =',
     &             2i7,1pe22.14,2i7)
             id4a = bin(kl,7) - 1.0d8*id3a
             write(6,2448) id4,id4a
 2448        format(1x,'id4,id4a = ',2i7)
             call flush(6)
*
             iname(nt0) = nt0
             ikind(im2) = -1
             ikind(nt0) = -1
             r(nt0) = r(kn)*1.0000000001d0
             xmax(nt0) = 1.01d0*r(nt0)
             vr(nt0) = vr(im2)
             vt(nt0) = vt(im2)
*
*       set the vector ibstra for binary components (now single stars)
*
             ibstra(im2) = 0
             ibstra(nt0) = 0
*
*               remove binary from the list nameb - zero means not binary
*               add star nt0 with id id3 to the names list
*
             nameb(im2) = 0
             names(im2) = id3
             names(nt0) = id4
             body(im2) = sy1/zmbar
             body(nt0) = sy2/zmbar
             sb11 = body(im2)
             sb12 = body(nt0)
             sb2 = sb11 + sb12
             print*,'kl,nt0,im2,names1-2,nameb,id3,id4,id3a,id4a=',
     &            kl,nt0,im2,names(im2),names(nt0),nameb(im2),id3,id4,
     &            id3a,id4a
             print*,'sy1,sy2,sm1b,sm2b = ',sy1,sy2,bin(kl,1)*zmbar,
     &                   bin(kl,2)*zmbar
             call flush(6)
*
             write(43,440) id3,im2,timevo,sy1,sy1,rad3i,xlum3i,r(kn),
     &                     ik21f,ik21f,ibstra(im2),ikind(im2)
             write(43,440) id4,nt0,timevo,sy2,sy2,rad4i,xlum4i,r(nt0),
     &                     ik22f,ik22f,ibstra(nt0),ikind(nt0)
*
*        compute the new binding energy of survived binary
*
             binold = bin(k,5)
             apold = bin(k,3)*rbar/rtidkg/rsuntopc
             bin(k,5) = bin(k,5) + deb + bin(kl,5)
             bin(k,3) = 0.5d0*bin(k,1)*bin(k,2)/bin(k,5)
             print*,'im1,k,binold,binnew,deb,apold,apnew,m1,m2=',im1,
     &               k,binold,bin(k,5),deb,apold,
     &               bin(k,3)*rbar/rtidkg/rsuntopc,bin(k,1)*zmbar,
     &               bin(k,2)*zmbar
*
*     update binary parameters for binary evolution procedures
*
             idb = nameb(im1)
             call get_loid(idb,id1)
             call get_hiid(idb,id2)
             call getRadius(id1,rad1)
             call getRadius(id2,rad2)
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             ik1o = ik1
             ik2o = ik2
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             ssm12o = sm1 + sm2
             ebold = bin(k,5)
             print*,'before diss- im1,k,idb,id1,id2,ik1,ik2,sm1,sm2 =',
     &               im1,k,idb,id1,id2,ik1,ik2,sm1,sm2
*        semi-major axis in R_sun unit
             ap = bin(k,3)*rbar/rtidkg/rsuntopc
             call setSemi(idb,ap)
             call getRL(id1,rlo)
*
*       synchronize the dynamical binary parameters (after interaction)
*       with the evolutionary status of the binary. Evolution with time
*       step equal to zero
*       
             call evBinary(idb,tphys)            
             call get_loid(idb,id1)  
             call get_hiid(idb,id2)  
             call getRadius(id1,rad1)
             call getRadius(id2,rad2)
             call get_ss_type(id1,ik1)
             call get_ss_type(id2,ik2)
             call get_loid_mass(idb,sm1)
             call get_hiid_mass(idb,sm2)
             print*,'after diss- im1,k,idb,id1,id2,ik1,ik2,sm1,sm2 =',
     &               im1,k,idb,id1,id2,ik1,ik2,sm1,sm2
             bin(k,1) = sm1/zmbar
             bin(k,2) = sm2/zmbar
             bin(k,7) = 1.0d8*float(id1) + float(id2)
             ssm12n = sm1 + sm2
             call getSemi(idb,apnew)
             call getRL(id1,rln)
             call binary_exists(idb,iexist)
             if(iexist.eq.0) then
               rln = -1.0d0      
               print*,'intb3b3bd-merger'
             else
               call getRL(id1,rln)
             endif
             call get_bs_type(idb,ik0n)
             dmm = ssm12o - ssm12n
             write(6,9877) im1,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,
     &                     dmm
 9877  format(1x,'2d im1,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,dmm='
     &       ,2i9,2i4,1p6e14.6)
*
*       check for exchange interaction
*
*           special treatment for exchange case 
*           (Heggie et al. 1996, ApJ, 467, 359)
*
             if(iexist.eq.0) go to 54
*             
             iexch1 = 0
             ixxx = 0
             if(iexch.eq.2) then
               if(sm3.ge.sm4) then
                 ids = id3
                 sm3 = sm3/zmbar
                 sm1 = sm1/zmbar
                 sm2 = sm2/zmbar
                 sm12 = sm1 + sm2
                 sm123 = sm12 + sm3
                 ixxx = 1
               else
                 ids = id4
                 sm4a = sm4
                 sm3a = sm3
                 sm3 = sm4/zmbar
                 sm1 = sm1/zmbar
                 sm2 = sm2/zmbar
                 sm12 = sm1 + sm2
                 sm123 = sm12 + sm3
                 ixxx = 2
               endif              
*
               yy = sm123/bin(k,3)
               if(ub*ub.gt.yy) go to 54
*                 
*      compute probability for exchange
* 
               id1a = bin(k,7)/1.0d8
               id2a = bin(k,7) - 1.0d8*id1a
               idb = nameb(im1)
               call get_loid(idb,id1)
               call get_hiid(idb,id2)
               call get_loid_mass(idb,sm1s)
               call get_hiid_mass(idb,sm2s)
               ids = names(im2)
               if(ixxx.eq.2) ids = names(nt0)
               call get_mass(ids,sm3s)
*
      print*,'3b3b-- id1a,id1,id2a,id2,ids,sm1,sm1s,sm2,sm2s,am3,sm3s,'
     &,'ub2,yy,ixxx=',id1a,id1,id2a,id2,ids,sm1*zmbar,sm1s,sm2*zmbar,
     &sm2s,sm3u*zmbar,sm3s,ub*ub,yy,ixxx                    
*
*     make sure that id of stars and their masses are the same in the
*     evolution procedures (get_loid(idb,id1), get_loid_mass(idb,sm1s)) 
*     and in the bin(k,1), bin(k,2) and bin(k,7)
* 
               if(id1.ne.id1a) then
                 if(sm1s.lt.0.99999d0*sm1*zmbar.or.
     &             sm1s.gt.1.00001*sm1*zmbar) then
                   smsx = sm1
                   sm1 = sm2
                   sm2 = smsx
                   id1x = id1a
                   id1a = id2a
                   id2a = id1x
                   bin(k,1) = sm1
                   bin(k,2) = sm2
                   bin(k,7) = 1.0d8*float(id1a) + float(id2a)
*
                   print*,'3b3b-- smsx,sm1,sm2,id1x,id1a,id2a=',
     &                 smsx*zmbar,sm1*zmbar,sm2*zmbar,id1x,id1a,id2a
*
                 endif
               endif
*                      
*       first star is exchanged
* 
               xm1 = sm1/sm12
               xm2 = sm3/sm123
               x23 = sm2 + sm3
               x13 = sm1 + sm3
               x23x123 = x23/sm123
               x3x13 = sm3/x13
               x123x12 = sm123/sm12
               x13x123 = x13/sm123
               xcof1=x23x123**one6*x3x13**3.5*x123x12**one3*x13x123
               xcof = 0.5d0*pi*sm123*bin(k,3)/ub/ub
               xform = a00 + a01*xm1 + a10*xm2 + a02*xm1**2 +
     &                 a11*xm1*xm2 + a20*xm2**2 + a03*xm1**3 +
     &                 a12*xm1**2*xm2 + a21*xm1*xm2**2 + a30*xm2**3         
*                  
               xform = exp(xform)
               pexch1 = 0.5d0*co*xcof*xcof1*xform*den*ub*dt/pi
               print*,'3b3b-- sm1,sm2,sm12,sm3,sm123,ub,a,sm123=',
     &         sm1*zmbar,sm2*zmbar,sm12*zmbar,sm3*zmbar,sm123*zmbar,
     &         ub,bin(k,3),sm123
               print*,'3b3b-- xcof,xcof1,xform,c1,den,ub,dt,iexch1=',
     &                      xcof,xcof1,xform,c1,den,ub,dt,iexch1
*
*       second star is exchanged
*
               xm1 = sm2/sm12
               xm2 = sm3/sm123
               x23 = sm1 + sm3
               x13 = sm2 + sm3
               x23x123 = x23/sm123
               x3x13 = sm3/x13
               x123x12 = sm123/sm12
               x13x123 = x13/sm123
               xcof1=x23x123**one6*x3x13**3.5*x123x12**one3*x13x123
               xcof = 0.5d0*pi*sm123*bin(k,3)/ub/ub
               xform = a00 + a01*xm1 + a10*xm2 + a02*xm1**2 +
     &                 a11*xm1*xm2 + a20*xm2**2 + a03*xm1**3 +
     &                 a12*xm1**2*xm2 + a21*xm1*xm2**2 + a30*xm2**3
*                  
               xform = exp(xform)
               pexch2 = 0.5d0*co*xcof*xcof1*xform*den*ub*dt/pi
               pexch = pexch1 + pexch2
*
               print*,'3b3b-- pexch1,pexch2,pexch=',pexch1,pexch2,pexch
               print*,'3b3b-- xcof,xcof1,xform,c1,den,ub,dt = ',xcof,
     &                      xcof1,xform,c1,den,ub,dt
*                                 
               xrun = ran2(irun)
*
               if(xrun.gt.pexch) go to 54
*                    
               iexch1 = 1
*
               nexchang2 = nexchang2 + 1
               nexch = nexch + 1
               print*,'nexchang2,nexch,iexch1',nexchang2,nexch,iexch1
*
               ebold = bin(k,5)
               aold = bin(k,3)
               eold = bin(k,4)
               sm12oo = sm1 + sm2
               sm3oo = sm3
               inexch2(im1) = inexch2(im1) + 1                 
*
               write(6,17) im2,im1,idb,id1,id2,ids,kl,inexch(im1)
*
*      check the ages of binary and single star. If they are different
*      synhronize them,
*
               call getAge(id1,age1)
               call getAge(ids,ages)
               if(age1.ne.ages) then
                 print*,'b3b3 wrong ages-2  im1,im2,age1,ages =',
     &                 im1,im2,age1,ages               
               endif
               timea = timevo
*
                sm12 = sm1 + sm2
                sm123 = sm12 + sm3
*         
                imx = im2
                sm4 = sm4/zmbar
                if(ixxx.eq.2) then
                  imx = nt0
                  sm4 = sm3a/zmbar
                endif             
                xrun = pexch*ran2(irun)
                if(xrun.le.pexch1) then
*
*       the first star is exchanged
*
                  bin(k,1) = sm3
                  bin(k,2) = sm2
                  body(im1) = sm3 + sm2
                  body(imx) = sm1
                  sm12n = sm3 + sm2
                  sb1 = sm12n
                  sb11 = body(imx)
                  sb12 = sm4
                  sb2 = sb11 + sb12
                  eold = bin(k,4)
                  bin(k,3) = 0.5d0*bin(k,1)*bin(k,2)/ebold
                  anew = bin(k,3)
                  bin(k,5) = ebold
                  enew = 1.d0 - (sm1/sm3)*(sm123/sm12)**one3
                  if(enew.lt.0.d0) enew = 0.d0
                  bin(k,4) = enew
                  vr(im1) = sqrt(vr(im1)*vr(im1)*sm12/sm12n)
                  vt(im1) = sqrt(vt(im1)*vt(im1)*sm12/sm12n)
                  vr(imx) = sqrt(vr(imx)*vr(imx)*sm3/sm1)
                  vt(imx) = sqrt(vt(imx)*vt(imx)*sm3/sm1)
                  bin(k,7) = 1.0d8*float(ids) + float(id2)
                  names(imx) = id1
                  if(names(imx).eq.nameb(im1))then
                    nameb(im1) = ids
                    idb = ids
                  endif
*
                  anb = bin(k,3)*rbar/rtidkg/rsuntopc
*
       print*,'3b3b--exchange-1:kl,kn,km,im1,im2,idb,id1,id2,ids,names,'
     &,'nameb,sm3,sm2,sm1,sm12o,sm12n,ebold,aold,anew,eold,enew,age1,'
     &,'ages,timea,,anb,anb-NB = ',kl,kn,km,im1,imx,idb,id1,id2,ids,
     &names(imx),nameb(im1),sm3*zmbar,sm2*zmbar,sm1*zmbar,
     &(sm1+sm2)*zmbar,sm12n*zmbar,ebold,aold,anew,eold,enew,age1,
     &ages,timea,anb,anb*rsuntopc*rtidkg/rbar                      
*
*      create a new binary which formed after an exchange interaction
*
                  call create_binary(idb,ids,id2,anb,enew,timea)
                  call get_loid_mass(idb,sm1d)
                  call get_hiid_mass(idb,sm2d)
                  call binary_exists(idb,iexist)
                  dmm = sm12n*zmbar - (sm1d + sm2d)
                  write(6,7863) im1,idb,iexist,sm12n*zmbar,
     &                            sm1d + sm2d,dmm
 7863        format(1x,'3b3b-creat3- im1,idb,iexist,sm12o,sm12n,dmm,',
     &        2i9,i4,1p5e14.6)
                  if(iexist.eq.0) print*,'intb3b3-1 iexist=0'
*                                                
                  go to 54
                endif
*
                if(xrun.gt.pexch1) then
*
*       the second star is exchanged
*
                bin(k,1) = sm1
                bin(k,2) = sm3
                body(im1) = sm1 + sm3
                body(imx) = sm2
                sm12n = sm1 + sm3
                sb1 = sm12n
                sb11 = body(imx)
                sb12 = sm4
                sb2 = sb11 + sb12
                eold = bin(k,4) 
                bin(k,3) = 0.5d0*bin(k,1)*bin(k,2)/ebold
                anew = bin(k,3)
                bin(k,5) = ebold
                enew = 1.d0 - (sm2/sm3)*(sm123/sm12)**one3
                if(enew.lt.0.d0) enew = 0.d0
                bin(k,4) = enew
                vr(im1) = sqrt(vr(im1)*vr(im1)*sm12/sm12n)
                vt(im1) = sqrt(vt(im1)*vt(im1)*sm12/sm12n)
                vr(imx) = sqrt(vr(imx)*vr(imx)*sm3/sm2)
                vt(imx) = sqrt(vt(imx)*vt(imx)*sm3/sm2)
                bin(k,7) = 1.0d8*float(id1) + float(ids)
                names(imx) = id2
                if(names(imx).eq.nameb(im1)) then
                  nameb(im1) = ids
                  idb = ids
                endif  
*
                anb = bin(k,3)*rbar/rtidkg/rsuntopc
*
      print*,'3b3b-- exchange-2:kl,kn,km,im1,im2,idb,id1,id2,ids,names,'
     &,'nameb,sm3,sm2,sm1,sm12o,sm12n,ebold,aold,anew,eold,enew,age1,'
     &,'ages,timea,,anb,anb-NB = ',kl,kn,km,im1,imx,idb,id1,id2,ids,
     &names(imx),nameb(im1),sm3*zmbar,sm2*zmbar,sm1*zmbar,
     &(sm1+sm2)*zmbar,sm12n*zmbar,ebold,aold,anew,eold,enew,age1,
     &ages,timea,anb,anb*rsuntopc*rtidkg/rbar
*
*      create a new binary which formed after an exchange interaction
*
                call create_binary(idb,id1,ids,anb,enew,timea)
                call get_loid_mass(idb,sm1d)
                call get_hiid_mass(idb,sm2d)
                call binary_exists(idb,iexist)
                dmm = sm12n*zmbar - (sm1d + sm2d)
                write(6,7862) im1,idb,iexist,sm12n*zmbar,
     &                            sm1d + sm2d,dmm
 7862        format(1x,'3b3b-creat4- im1,idb,iexist,sm12o,sm12n,dmm,',
     &        2i9,i4,1p5e14.6)
                if(iexist.eq.0) print*,'intb3b3-2 iexist=0'
*
                go to 54
*
              endif
*
            endif
*
 54         continue
*
*       check if semi major axis is greater than the sum of star radii
*       if not then the merger is due because of 4-body interaction.
*       Here only flag is set to distinguish between mergers due to binary
*       evolution and due to 4-body interaction.
*
            sm1 = bin(k,1)*zmbar
            sm2 = bin(k,2)*zmbar
            call get_loid(idb,id1)
            call get_hiid(idb,id2)
            call getRadius(id1,rad1)
            call getRadius(id2,rad2)
            call get_ss_type(id1,ik1)
            call get_ss_type(id2,ik2)
            call get_loid_mass(idb,sm1d)
            call get_hiid_mass(idb,sm2d)
            call getEcc(idb,ecc1i)
            call getSemi(idb,ap1i)
            call getLum(id1,xlum1i)
            call getLum(id2,xlum2i)
            call get_bs_type(idb,ik01f)
            sm1f = sm1
            sm2f = sm2
            rad1i = rad1
            rad2i = rad2
            ik11f = ik1
            ik12f = ik2
            idb1 = idb
            id11 = id1
            id12 = id2
            if(abs(sm1+sm2-sm1d-sm2d).gt.1.d-10) then
              write(6,6688) im1,idb,id1,id2,ik1,ik2,sm1,sm2,sm1d,sm2d,
     &                       sm1 + sm2 - (sm1d + sm2d)
 6688         format(1x,'wrong-diss-mass-2 im1,idb,id1,id2,ik1,ik2,',
     &'sm1,sm2,sm1d,sm2d,dmmdis=',7i9,1p5e14.6)
            endif  
*
            call binary_exists(idb,iexist)
            if(ik1o.le.1.and.ik2o.le.1) then
              if(rad1+rad2.gt.ap.or.iexist.eq.0) then
                if(sm1+sm2.gt.sturn.and.(ik1.le.1.or.ik2.le.1)) then
                  print*,'33-4 sm12,sturn,t=',sm1+sm2,sturn,timevo
                  call flush(6)
                  ibstra(im1) = 4
                  ibs4 = ibs4 + 1
      print*,'b3b3-4  k,im1,idb,sm1,sm2,smsss,a,rad1,rad2,bin,ap,ibs4=',
     &     k,im1,idb,sm1,sm2,smsss,bin(k,3),rad1,rad2,bin(k,5),ap,ibs4
                  call flush(6)
                endif
              endif
            endif
*
*        update time for the next evolution step
*
            call get_bs_updatetime(idb,uptime(im1))
*
*        check for a merger event
*
            nob = 0
            if(iexist.eq.0) then
              ecc1i = -1000.d0
              ap1i = -1000.d0
              ikind(im1) = 3
              names(im1) = idb
              nameb(im1) = 0
              call get_loid_mass(idb,sm1)
              call get_hiid_mass(idb,sm2)
              if(ik2.eq.15) then
                bin(k,2) = 0.d0
                bin(k,1) = sm1/zmbar
                names(im1) = id1
                sm2 = 0.d0
                sm2f = 0.d0
              endif
              if(ik1.eq.15) then 
                bin(k,2) = sm2/zmbar
                bin(k,1) = 0.d0
                names(im1) = id2
                sm1 = 0.d0
                sm1f = 0.d0
              endif
              bin(k,3) = 0.d0
              bin(k,4) = 0.d0
              body(im1) = bin(k,1) + bin(k,2)
              sb1 = body(im1)
              nmerge = nmerge + 1
              sm1234 = sb1 + sb2
              iime = iime + 1
              cfr = body(im1)*zmbar
              ssme = ssme + cfr
      print*,'yyy8 timevo,iime,cfr,ssme= ',timevo,iime,cfr,ssme
      print*,'merger binary interaction intb3b3-4 k,sm1,sm2,nmerge,',
     &'im1,idb,iexist,ikind = ',k,sm1,sm2,nmerge,im1,idb,iexist,
     &ikind(im1)
              call flush(6)
*
              if(sm1.le.tolm.and.sm2.le.tolm) then
*
*      mark obliterated binary (massless) to removal from the system
*
                sm1f = 0.d0
                sm2f = 0.d0
                body(im1) = 0.d0
                sb1 = 0.d0
                sm1234 = sb1 + sb2
                bin(k,1) = 0.d0  
                bin(k,2) = 0.d0  
                iob = iob + 1     
                iobt = iobt + 1   
                nob = 1 
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
                r(kmm) = 1.0d+15 + im1
                xmax(im1) = 1.01d0*r(kmm)
                print*,'intb3b3-obld iob,im1,sm1,sm2=',iob,im1,sm1,sm2
              endif
              call flush(6)
            endif
*
            if(dmm.le.1.d-10) then
              dmm = 0.d0
            else
              ssevo = dmm/zmbar
              slosev = slosev + dmm
              nmloev = nmloev + 2  
              nsevo = nsevo + 1    
              ehkin = 0.5d0*ssevo*(vr(im1)**2 + vt(im1)**2)
              ehkin1 = ehkin1 + ehkin
              ehmlev = ehmlev - ehkin
              smt = smt - ssevo
              write(6,9328) im1,nob,nsevo,dmm,ssevo,ehkin
 9328         format(1x,'intb3b3d-evo im1,nob,nsevo,dmm,ssevo,ehkin=',
     &                 3i9,1p3e14.6)
            endif
*
            ikind(im1) = -ikind(im1)
*
            if(sb1.gt.0.d0) then
              deb1 = sb2*deb/sm1234
              deb2 = sb1*deb/sm1234
              deb3 = sb12*deb2/sb2 
              deb4 = sb11*deb2/sb2 
              call kick(kmm,deb1)   
              write(6,*) 'b3b3>0.12 k,kmm,im1,deb1,ebin=',
     &                   k,kmm,im1,deb1,bin(k,5)
              call flush(6)
            else
              deb1 = ebold
              deb3 = sb12*deb/sm1234
              deb4 = sb11*deb/sm1234
            endif
*
*         compute kick velocity
*
            print*,'sb1,deb,deb1-4 b=',sb1,deb,deb1,deb2,deb3,deb4
            call kick(kn,deb3)
            write(6,*) 'b3b3>0.12 ** kl,kn,im2,id3,deb3,ebin=',
     &                   kl,kn,im2,names(im2),deb3,bin(kl,5)
            write(6,*) 'nt,nn,nt0,im2,id3,id4 = ',nt,n+ndist+ndistm,
     &                  nt0,im2,id3,id4
            call flush(6)
            call kick(nt0,deb4)
            write(6,*) 'b3b3>0.12 ** nt0,id4,deb4=',nt0,names(nt0),deb4
            call flush(6)
*
*      store information about binary binary interaction - first binary
*
            iinte3(k) = iinte3(k) + 1
            inum = iinte3(k)
*
            if(inum.gt.50) then
              iinte3(k) = 1
              inum = 1
            endif
*
            inumt = 50*(k - 1) + inum
            binin(inumt,1) = r(kmm)/pctonbu
            binin(inumt,2) = deb1
            binin(inumt,3) = sb1*zmbar
            binin(inumt,4) = sb2*zmbar
            binin(inumt,5) = timevo
            binin(inumt,6) = bin(k,5)
            iwww = 0
            if(iexist.eq.0) then
              binin(inumt,7) = 5
              bin(k,5) = 0.d0
              iwww = 5
              if(iexch1.eq.1) binin(inumt,7) = 9 + iwww
            else
              binin(inumt,7) = 2
              if(iexch1.eq.1) binin(inumt,7) = 9
            endif
            if(iexch1.eq.1) binin(inumt,7) = 9 + iwww
*
            ibbb = binin(inumt,7)
            write(44,420) idb1,id11,id21,im1,timevo,sm1i,sm1f,sm2i,
     &                 sm2f,rad1i,rad2i,xlum1i,xlum2i,ap1i,ecc1i,r(kmm),
     &                 ik01,ik01f,ik11,ik11f,ik12,ik12f,iexist,
     &                 ibstra(im1),ibbb
*
*      store information about binary binary interaction - second binary
*
            iinte3(kl) = iinte3(kl) + 1
            inum = iinte3(kl)
*
            if(inum.gt.50) then
              iinte3(kl) = 1
              inum = 1
            endif
*
            inumt = 50*(kl - 1) + inum
            binin(inumt,1) = r(kn)/pctonbu
            binin(inumt,2) = -bin(kl,5)
            binin(inumt,3) = sb1*zmbar
            binin(inumt,4) = sb2*zmbar
            binin(inumt,5) = timevo
            binin(inumt,6) = bin(kl,5)
            binin(inumt,7) = 7
*
            ibbb = binin(inumt,7)
            write(44,420) idb2,id21,id22,im2,timevo,sm3i,sm3f,sm4i,
     &                 sm4f,rad3i,rad4i,xlum3i,xlum4i,ap2i,ecc2i,r(kn),
     &                 ik02,ik02f,ik21,ik21f,ik22,ik22f,iexist,
     &                 ibstra(im2),ibbb
*        
*        update the list of binaries in the system to distinguish
*        between binaries removed due to escape (relaxation bin(kl,6),
*        interaction bin(kl,7)) and dissolved bin(kl,4)
*
            bin(kl,4) = -bin(kl,4)
*
          endif
*
         endif
*
         go to 110
*
 100  continue
*
      close(43)
      close(44)
*      
      if(ib3b3.eq.0.and.noib3b3.eq.0) return   
*
*     deal with change of number of single stars due to binary disruption
*     compute new potential and sort stars after destroing binaries in
*     binary--binary interactions. One component of destroied binary
*     is substituded for centre of binary mass, for another component
*     name of one of stars which previeously created binary is taken
*
      ndist4 = ndist4 + ndist
      n = n + ndist + ndistm
      nt = nt - iob
c      call sort2(nt0,r,iname)
      call sort3(nt0,r,iname,inameo)
      nzst(nsuzon) = nt
*
      call coepot
      call energy(2)
      ptnew = pot
      zknew = zkin
      ehmlevn = ehmlev
      ehb3b3n = ehb3b3
      ekicktn = ekickt
      write(6,5432) zkold,zknew,ehmlevo,ehmlevn,ehbin3o,ehbin3n,      
     &              ekickto,ekicktn,zknew-zkold,ehmlevn-ehmlevo,
     &              ehbin3n-ehbin3o,ekicktn-ekickto,(zknew-zkold)-
     &              (ehmlevn-ehmlevo+ehbin3n-ehbin3o+ekicktn-ekickto)
 5432 format(1x,'binb3b3  zkold,zknew,ehmlevo,ehmlevn,ekickto,ekicktn',
     &',ehbin3o,ehbin3n,Dzk,Dhm,Dhb3,Dhki,DDb3b3 =',1p13e14.6)
*                                                                             
      enepot = enepot + ptnew - ptold
      etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
     &         ehmlev - ehcoll + eccoll - ekickt - enekin
       write(6,6541) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &               ekickt,enekin,escsta,ehbin3,ehb3b3
 6541 format(1x,'bin3b3-f etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12) 
*
      print*,'potold,potnew,zkold,zknew,n,nt,ndist,ndistm,nexch,',
     &'nexchang2,ib3b3,noib3b3,iob,iobx,iobt = ',ptold,ptnew,zkold,
     &zknew,n,nt,ndist,ndistm,nexch,nexchang2,ib3b3,noib3b3,iob,iobx,
     &iobt
*                                                                             
      call flush(6)
*
      return
*
      end
*
*
*
