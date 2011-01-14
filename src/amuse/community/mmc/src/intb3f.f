*
*
      subroutine intb3f(nup)
*
*      calculate interaction between t binaries and field 
*      --------------------------------------------------
*      stars according to Spitzer's or Heggie's formulae.
*      --------------------------------------------------
*             ib3f = 1   -  Spitzer
*             ---------------------
*             ib3f = 2   -  Heggie
*             --------------------
*             ib3f = 3   - Pmax for interaction probability
*             ---------------------------------------------
*             ib3f = 4   - numerical integration
*             ----------------------------------
*
*     intb3f is the second subroutine in the interactions sequence
*     (intb3b3 - intb3bf - formb3) so ikind can be negative
*
*
      include 'common.h'
*
      real*8 dt,co,velave,smave,den,pb3f,zz,deb,sm1,sm2,sm3,sm12,debb,
     &       deb3,zc,zd,ze,zq,ebdeb,xx1,xx2,ap,ub,sm123,xxmiu,vc2,vvc,
     &       pmax,rad1,rad2,timevo,tlog,smlog,smsss,sturn,z,zpars(20),
     &       c1,pppp,a00,a01,a02,a03,a10,a11,a12,a20,a21,a30,yy,xm1,
     &       xm2,x23x123,x3x13,x123x12,x13x123,x23,x13,xcof,xform,
     &       pexch1,pexch2,pexch,aold,ebold,eold,sm3s,smsx,age1,
     &       ages,sm12o,ptold,potnew,timea,dmass,cor,ehpot,ehkin,
     &       xcof1,anew,sm12oo,sm3oo,tform,tb3f,tb3b3,tcoll,
     &       anb,enew,rpmax,sm12n,sm1s,sm2s,sm3o,sq,xrun,pctonbu,sma,
     &       ecc,ekickb,ekicks,vrol1,vtol1,potold,zkino,zkinn,sss,
     &       tphys,ehmlevo,ehmlevn,ehbin3o,ehbin3n,sm123o,sm123n,
     &       sm123nn,ss123,rlo,rln,dm3f,ehkin1,ekickto,ekicktn,r1,r2,
     &       xlum1,xlum2,sm1oo,sm2oo,slum
*
      real*4 ran2
*
      integer l,nup,lmin,lmax,i,im1,im2,k,idb,inum,inumt,id1,id2,
     &        ii2,n,iexist,ii,isingl,im,ij,ik1,ik2,iexch1,ids,ik0
     &        id1a,id2a,iwww,nexch,id1x,iexi,ifint,nfint,ntx,ndistm,
     &        ik0n,iob,nsevo,nob,iobx,ik1o,ik2o,ik1oo,ik2oo,ibbb
*
      common /timset/ tform(20),tb3f(20),tb3b3(20),tcoll(20)     
*
      data a00,a01,a02,a03,a10,a11,a12,a20,a21,a30 /3.70d0,7.49d0,
     &     -15.49d0,3.07d0,-1.89d0,-2.93d0,13.15d0,-2.92d0,-5.23d0,
     &     3.12d0/
*
      if((nbin3-nescb3-ndist3-ndist4-ndiste-nmerge).eq.0) return
*
      open(44,file=trim(datadir)//'/binarym.dat',access='append')
      open(43,file=trim(datadir)//'/starm.dat',access='append')
*
      n = nt
      c1 = float(nt)/log(gamma*nt)
      co = sqrt(3.0d0)*c1
      nexch = 0
      ndistm = 0
      ifint = 0
      nfint = 0
      nsevo = 0
      iexi = 0
      iob = 0
      ekiskb = 0.d0
      ekicks = 0.d0
      ehkin1 = 0.d0
      pctonbu = rtidkg/rbar
*
      if(ib3f.eq.1) then
*
*       co = 5*pi*As/(8*4*pi/3)*2/3
*
        co = co*21.0d0*10.0d0/32.0d0
      endif
*
      if(ib3f.eq.2) then
*
*       co = 2*Ah*/(7*4*pi/3)*2/3
*
        co = co*30.0d0*2.0d0/(14.0d0*pi)
      endif
*
      if(alphal.eq.-1.0d0) then
*
*          equal mass case
*
        co = co
      else
*
*          unequal mass case
*
        co = co*sqrt(2.0d0/3.0d0)
      endif       
*
      do 10 l = 1,nup
*
         if(l.eq.1) then
           lmin = 1
         else
           lmin = nzst(l-1)+1
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
         isingl = 0
*
*     calculate actual time
*
         tb3f(l) = tb3f(l)+dt*tscale0*log(gamma*nt00)
     &             /log(gamma*ntnew)*float(ntnew)/float(nt00)
         timevo = tb3f(l)
         print*,'intb3f- l,nup,tb3f(l) = ',l,nup,tb3f(l)
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
         print*,'b3f- smsss,sturn,z0,tim=',smsss,sturn,z,timevo
*                                                                      
         do 20 i = lmin,lmax
*
            ssevo = 0.d0
            if(i.gt.nt) go to 20
*
            im1 = iname(i)
*
            if(abs(ikind(im1)).eq.2) then
              k = nwhich(im1)
*
*   do not check binaries which escaped due to relaxation or interactions or
*   were dissolved. Binaries which take part in interactions in this
*   time steap have naegative ikind. They can still take part with
*   interaction with other binaries. So ikind() ---> abs(ikins())
*
              if(abs(body(im1)).lt.tolm) go to 20
              if(abs(bin(k,1)).lt.tolm.or.
     &           abs(bin(k,2)).lt.tolm)  go to 20
              if(r(i).gt.1.d8) go to 20
              if(bin(k,4).lt.0.0d0) go to 20
              if(bin(k,6).lt.0.0d0) go to 20
              if(bin(k,7).lt.0.0d0) go to 20
*
              call veloci(i,velave,smave)
              call densi(i,den)
*
*     find single star to interact with binary if isingl = 0 (no single
*     stars "below" binary
*
              do 21 ij = i+1,lmax
                 im = iname(ij)
                 if(abs(ikind(im)).eq.2) go to 21
                 if(abs(body(im)).lt.tolm) go to 21
                 if(r(ij).gt.1.d8) go to 21
                 if(nkick(im).ge.1) go to 21
                 ids = names(im)
                 call getLum(ids,slum)
                 if(isNaN(slum)) then
                   print*,'21- i,im,ids,lum = ',ij,im,ids,slum
                   go to 21
                 endif
                 isingl = ij
                 go to 22
 21           continue
 22           continue
*
              if(isingl.eq.0) then
                do 23 ij = i-1,lmin,-1
                   im = iname(ij)
                   if(abs(ikind(im)).eq.2) go to 23
                   if(abs(body(im)).lt.tolm) go to 23
                   if(nkick(im).ge.1) go to 23
                   if(r(ij).gt.1.d8) go to 23
                   ids = names(im)
                   call getLum(ids,slum)
                   if(isNaN(slum)) then 
                     print*,'23- i,im,ids,lum = ',ij,im,ids,slum
                     go to 23
                   endif
                   isingl = ij
                   go to 24
 23             continue
              endif
 24           continue
*             
*
              if(isingl.eq.0) go to 20
              im2 = iname(isingl)
              ii = isingl
              isingl = 0
*
*          decide if there is an interaction between binary and field star
*
              sm12 = body(im1)
              sm1 = bin(k,1)
              sm2 = bin(k,2)
*                                                       
              if(ib3f.le.2) then
                if(alphal.eq.-1.0d0) then
*
*                  equal mass case
*
                  pb3f = co*smave**3*den*dt/sqrt(velave)/bin(k,5)
                else
*
*                  unequal mass case
*
                  xx1 = (sm1*sm2)**2
                  xx2 = sm12/smave
                  xx2 = 1.0d0 + xx2
                  xx2 = xx2/smave/sm12
                  xx2 = sqrt(xx2)
                  xx1 = xx1*xx2
                  pb3f = co*xx1*den*dt/sqrt(velave)/bin(k,5)
                endif
*
              else
*
*    compute encounter probability using Pmax 
*
                call relvel(i,ii,ub)
*
                sm3 = body(im2)
                sm123 = body(im1) + sm3
*
                rpmax = 3.d0*bin(k,3)
                sq = 1.d0 + 2.d0*sm123/ub/ub/rpmax
                pmax = rpmax*sqrt(sq)
                                          
*
*          number density n=de/2/pi
*
                pb3f = 0.5d0*c1*pmax**2*den*ub*dt
*
              endif                      
                      
              bin(k,8) = r(i)
              xrun = ran2(irun)
*
cChanged DCH 2/8/6
              if(xrun.le.pb3f.and.body(im1)*body(im2).ne.0.d0)then
*
*   evolve binary and single star to the time timevo - time of binary
*   single star interaction
*
                if(ifint.eq.0.and.nfint.eq.0) then
                  call energy(2)   
                  potold = pot
                  zkino = zkin
                  ehmlevo = ehmlev
                  ehbin3o = ehbin3
                  ekickto = ekickt
                endif 
*
                ifint = ifint + 1
                ntx = nt
                iobx = iobt
                idb = nameb(im1)
                ids = names(im2)
                call get_loid_mass(idb,sm1s)
                call get_hiid_mass(idb,sm2s)
                call get_mass(ids,sm3s)
                sm123o = sm1s + sm2s + sm3s
                write(6,7654) n,nt,nt0,im1,im2,idb,ids,k,i,ii,
     &                (body(im1)-bin(k,1)-bin(k,2))*zmbar,
     &                (body(im1)*zmbar - sm1s - sm2s),
     &                (body(im2)*zmbar - sm3s)
 7654           format(1x,'n,nt,nt0,im1,im2,idb,ids,k,i,ii,Dsmbd,',
     &'Dsmb,Dsms =',10i9,1p3e14.6)
                call flush(6)                                        
*                                                                                       *
                call mloss_single(i,timevo,ekickb,iexist)
                call mloss_single(ii,timevo,ekicks,iexi)
*                
c                timevo = timevo + 2.d-14
                tphys = timevo
                idb = nameb(im1)
                ids = names(im2)
                sm123n = (body(im1) + body(im2))*zmbar
                sm123nn = (bin(k,1) + bin(k,2) + body(im2))*zmbar
                sm3 = body(im2)
                sm1 = bin(k,1)
                sm2 = bin(k,2)
                sm12 = sm1 + sm2
                sm123 = body(im1) + sm3
                ndistm = ndistm + (nt - ntx) - (iobt - iobx) 
c                ekickt = ekickt + ekickb + ekicks
       print*,'intb3f-mloss_single n,nt,nt0,im1,im2,idb,ids,iexist,',
     &'iexi,ifint=',n,nt,nt0,im1,im2,idb,ids,iexist,iexi,ifint 
                if(iexist.eq.0.or.iexi.eq.0) then 
                  nfint = nfint + 1
                  print*,'intb3f-no  im1,im2,iexist,iexi,nfint =',
     &                    im1,im2,iexist,iexi,nfint
                
                  go to 20
                endif
*
                call get_loid(idb,id1)
                call get_hiid(idb,id2)
                call get_loid_mass(idb,sm1s)
                call get_hiid_mass(idb,sm2s)
                call get_ss_type(id1,ik1o)
                call get_ss_type(id2,ik2o)
                call get_mass(ids,sm3s)
                call get_bs_type(idb,ik0)
                sm1oo = sm1s
                sm2oo = sm2s
                ik1oo = ik1o
                ik2oo = ik2o
                ss123 = sm1s + sm2s + sm3s
                if(sm3s.lt.tolm*zmbar) go to 20
                call getSemi(idb,ap)
                call getEcc(idb,ecc)
                write(6,6543),im1,im2,sm123o,sm123*zmbar,
     &                   sm123n,sm123nn,sm123nn - ss123
 6543           format(1x,'im1,im2,sm123o,sm123,sm123n,sm123nn,dm3bf =',
     &                 2i9,1p5e14.6)
                call flush(6)
*
*           compute energy generated by binary and kick velocities of
*           the binary and field star
*
*           special treatment for exchange case 
*           (Heggie et al. 1996, ApJ, 467, 359)
*
                iexch1 = 0
                if(iexch.ge.1) then
                  yy = sm123/bin(k,3)
                  if(ub*ub.gt.yy) go to 51
*                  
*      compute probability for exchange
* 
                    id1a = bin(k,7)/1.0d8
                    id2a = bin(k,7) - 1.0d8*id1a
*
      print*,'id1a,id1,id2a,id2,ids,sm1,sm1s,sm2,sm2s,am3,sm3s,ub2,yy=',
     & id1a,id1,id2a,id2,ids,sm1*zmbar,sm1s,sm2*zmbar,sm2s,sm3*zmbar,
     & sm3s,ub*ub,yy                    
*
*     make sure that id of stars and their masses are the same in the
*     evolution procedures (get_loid(idb,id1), get_loid_mass(idb,sm1s)) 
*     and in the bin(k,1), bin(k,2) and bin(k,7)
* 
                    if(id1.ne.id1a) then
                      if(sm1s.lt.0.99999d0*sm1*zmbar.or.
     &                  sm1s.gt.1.00001*sm1*zmbar) then
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
                        print*,'smsx,sm1,sm2,id1x,id1a,id2a=',
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
                    xcof = 0.5d0*pi*sm123*bin(k,3)/ub/ub
                    xform = a00 + a01*xm1 + a10*xm2 + a02*xm1**2 +
     &                      a11*xm1*xm2 + a20*xm2**2 + a03*xm1**3 +
     &                      a12*xm1**2*xm2 + a21*xm1*xm2**2 + a30*xm2**3         
*                  
                    xform = exp(xform)
                    pexch1 = 0.5d0*c1*xcof*xcof1*xform*den*ub*dt/pi
                    print*,'sm1,sm2,sm12,sm3,sm123,ub,a,sm123=',
     &         sm1*zmbar,sm2*zmbar,sm12*zmbar,sm3*zmbar,sm123*zmbar,
     &         ub,bin(k,3),sm123
                    print*,'xcof,xcof1,xform,c1,den,ub,dt,iexch1= ',
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
     &                      a11*xm1*xm2 + a20*xm2**2 + a03*xm1**3 +
     &                      a12*xm1**2*xm2 + a21*xm1*xm2**2 + a30*xm2**3
*                    
                    xform = exp(xform)
                    pexch2 = 0.5d0*c1*xcof*xcof1*xform*den*ub*dt/pi
                    pexch = pexch1 + pexch2
*
                    print*,'pexch1,pexch2,pexch=',pexch1,pexch2,pexch
                    print*,'xcof,xcof1,xform,c1,den,ub,dt = ',xcof,
     &                      xcof1,xform,c1,den,ub,dt
*                                 
                    xrun = ran2(irun)
*
                    if(xrun.gt.pexch) go to 51
*                    
                    iexch1 = 1
*
                    nexchang = nexchang + 1
                    nexch = nexch + 1
                    print*,'nexchang,nexch,iexch1',nexchang,nexch,iexch1
*
                    ebold = bin(k,5)
                    aold = bin(k,3)
                    sm12oo = sm1 + sm2
                    sm3oo = sm3
                    inexch(im1) = inexch(im1) + 1                 
*
                    write(6,17) im1,im2,idb,id1,id2,ids,k,inexch(im1)
 17                 format(1x,'im1,im2,idb,id1,id2,ids,k,inexch =',9i8)
*
*      check the ages of binary and single star. They should be the same
*
                    call getAge(id1,age1)
                    call getAge(ids,ages)
                    if(age1.ne.ages) then
                      print*,'intb3f wrong ages  im1,im2,age1,ages =',
     &                     im1,im2,age1,ages               
                    endif
                    timea = timevo
*
                    xrun = pexch*ran2(irun)
                    if(xrun.le.pexch1) then
*
*       the first star is exchanged
*
                      bin(k,1) = sm3
                      bin(k,2) = sm2
                      body(im1) = sm3 + sm2
                      body(im2) = sm1
                      sm12n = sm3 + sm2
                      aold = ap*rsuntopc*pctonbu
                      ebold = 0.5d0*sm1*sm2/aold
                      eold = ecc
                      bin(k,3) = 0.5d0*bin(k,1)*bin(k,2)/ebold
                      anew = bin(k,3)
                      bin(k,5) = ebold
                      enew = 1.d0 - (sm1/sm3)*(sm123/sm12)**one3
                      if(enew.lt.0.d0) enew = 0.d0
                      bin(k,4) = enew
                      deb = 0.d0
                      sss = sqrt(sm12/sm12n)
                      vr(im1) = vr(im1)*sss
                      vt(im1) = vt(im1)*sss
                      sss = sqrt(sm3/sm1)
                      vr(im2) = vr(im2)*sss
                      vt(im2) = vt(im2)*sss
                      bin(k,7) = 1.0d8*float(ids) + float(id2)
                      names(im2) = id1
                      if(names(im2).eq.nameb(im1))then
                        nameb(im1) = ids
                        idb = ids
                      endif
*
                      anb = bin(k,3)*rbar/rtidkg/rsuntopc
*
*      create a new binary which formed after an exchange interaction
*
                      call create_binary(idb,ids,id2,anb,enew,timea)
                      call get_loid_mass(idb,sm1s)
                      call get_hiid_mass(idb,sm2s)              
                      call binary_exists(idb,iexist)
                      dm3f = (sm3 + sm2)*zmbar - (sm1s + sm2s)
                      write(6,7875) im1,idb,iexist,sm1*zmbar,sm2*zmbar,   
     &                sm1s,sm2s,dm3f
 7875  format(1x,'3bf-creat1- im1,idb,iexist,sm1,sm2,sm1s,sm2s,dm3f=',        
     &        2i9,i4,1p5e14.6) 
                      if(iexist.eq.0) print*,'intb3f-1 iexist=0'
*
         print*,'exchange-1:k,i,ii,im1,im2,idb,id1,id2,ids,names,nameb'
     &,',sm3,sm2,sm1,sm12o,sm12n,ebold,aold,anew,eold,enew,age1,ages'
     &,',timea,,anb,anb-NB,iexist = ',k,i,ii,im1,im2,idb,id1,id2,ids,
     &names(im2),nameb(im1),sm3*zmbar,sm2*zmbar,sm1*zmbar,
     &(sm1+sm2)*zmbar,sm12n*zmbar,ebold,aold,anew,eold,enew,age1,
     &ages,timea,anb,anb*rsuntopc*rtidkg/rbar,iexist                      
*
                      go to 51
                    endif
*
                    if(xrun.gt.pexch1) then
*
*       the second star is exchanged
*
                      bin(k,1) = sm1
                      bin(k,2) = sm3
                      body(im1) = sm1 + sm3
                      body(im2) = sm2
                      sm12n = sm1 + sm3
                      aold = ap*rsuntopc*pctonbu
                      ebold = 0.5d0*sm1*sm2/aold
                      eold = ecc 
                      bin(k,3) = 0.5d0*bin(k,1)*bin(k,2)/ebold
                      anew = bin(k,3)
                      bin(k,5) = ebold
                      enew = 1.d0 - (sm2/sm3)*(sm123/sm12)**one3
                      if(enew.lt.0.d0) enew = 0.d0
                      bin(k,4) = enew
                      deb = 0.d0
                      sss = sqrt(sm12/sm12n)
                      vr(im1) = vr(im1)*sss
                      vt(im1) = vt(im1)*sss
                      sss = sqrt(sm3/sm2)
                      vr(im2) = vr(im2)*sss
                      vt(im2) = vt(im2)*sss
                      bin(k,7) = 1.0d8*float(id1) + float(ids)
                      names(im2) = id2
                      if(names(im2).eq.nameb(im1)) then
                        nameb(im1) = ids
                        idb = ids
                      endif  
*
                      anb = bin(k,3)*rbar/rtidkg/rsuntopc
*
*      create a new binary which formed after an exchange interaction
*
                      call create_binary(idb,id1,ids,anb,enew,timea)
                      call get_loid_mass(idb,sm1s)
                      call get_hiid_mass(idb,sm2s)
                      call binary_exists(idb,iexist)
                      dm3f = (sm1 + sm3)*zmbar - (sm1s + sm2s)
                      write(6,7876) im1,idb,iexist,sm1*zmbar,sm2*zmbar,
     &                sm1s,sm2s,dm3f                 
 7876  format(1x,'3bf-creat2- im1,idb,iexist,sm1,sm2,sm1s,sm2s,dm3f=',
     &         2i9,i4,1p5e14.6)
                      if(iexist.eq.0) print*,'intb3f-2 iexist=0'
*
         print*,'exchange-2:k,i,ii,im1,im2,idb,id1,id2,ids,names,nameb'
     &,',sm1,sm3,sm2,sm12o,sm12n,ebold,aold,anew,eold,enew,age1,ages'  
     &,',timea,,anb,anb-NB,iexist = ',k,i,ii,im1,im2,idb,id1,id2,ids,
     &names(im2),nameb(im1),sm1*zmbar,sm3*zmbar,sm2*zmbar,
     &(sm1+sm2)*zmbar,sm12n*zmbar,ebold,aold,anew,eold,enew,age1,
     &ages,timea,anb,anb*rsuntopc*rtidkg/rbar,iexist
*
                      go to 51
*
                    endif                                                                                                                                                                                                                   
*
                endif
*
 51             continue                    
*
                if(iexch1.eq.1) go to 54
                                                      
                if(ib3f.le.3) then
*                
 50               continue
                  xrun = ran2(irun)
                  zz = 0.5d0*pi*xrun     
                  zq = ran2(irun)
                  zc = cos(zz)
* 
                  if(ib3f.eq.1) then
                    zc = zc**6
                    if(zq.gt.zc) go to 50
                    deb = bin(k,5)*tan(zz)**2
                  else
                    zd = sin(zz)
                    ze = zc**2.5/(zc + zd)**4.5
                    if(zq.gt.ze) go to 50
                    deb = bin(k,5)*tan(zz)
                  endif
*
                  print*,'zz,zq,zc,zd,ze,xrun,irun =',zz,zq,zc,zd,ze,
     &                    xrun,irun 
                  ebdeb = bin(k,5)
                  bin(k,5) = bin(k,5) + deb
                  bin(k,3) = 0.5d0*sm1*sm2/bin(k,5)
*
                  call getSemi(idb,sss)
                  ap = bin(k,3)/pctonbu/rsuntopc
                  call setSemi(idb,ap)
                  call getRL(id1,rlo)
                  call get_ss_type(id1,ik1o)
                  call get_ss_type(id2,ik2o)
*
*       synchronize the dynamical binary parameters (after interaction)
*       with the evolutionary status of the binary. Evolution with time
*       step equal to zero
*
                  call evBinary(idb,tphys)
                  call get_loid_mass(idb,sm1s)
                  call get_hiid_mass(idb,sm2s)
                  call get_mass(ids,sm3s)
                  call getSemi(idb,zd)
                  call binary_exists(idb,iexist)
                  if(iexist.eq.0) then
                    rln = -1.0d0
                    print*,'intb3f-merger'
                  else
                    call getRL(id1,rln)
                  endif
                  call get_bs_type(idb,ik0n)
                  dm3f = (sm1 + sm2)*zmbar - (sm1s + sm2s)
                  write(6,9876) im1,idb,iexist,ik0n,sss,ap,zd,rlo,rln,
     &                          dm3f
 9876  format(1x,'3bf- im1,idb,iexist,ik0n,apold,ap,apnew,rlo,rln,dm3f='
     &        ,2i9,2i4,1p6e14.6)
                  write(6,8765) im1,im2,idb,ids,tphys,timevo,sm1s,sm2s,
     &                          sm3s,deb,ebdeb,bin(k,5),k
 8765  format(1x,'im1,im2,idb,ids,tphys,timevo,sm1s,sm2s,sm3s,deb,',
     &'ebdeb,ebnew,k=',4i9,1p8e14.6,i9)
                  ap = zd 
                else
*
*
*           numerical integration
*
*
                  ap = bin(k,3)/pctonbu/rsuntopc
                  call setSemi(idb,ap)
                  call get_ss_type(id1,ik1o)
                  call get_ss_type(id2,ik2o)
*
*       synchronize the dynamical binary parameters (after interaction)
*       with the evolutionary status of the binary. Evolution with time
*       step equal to zero
*
                  call evBinary(idb,tphys)
*
                endif
*
 54             continue                                
*
*     update binary parameters for binary evolution procedures
*
                call get_loid(idb,id1)
                call get_hiid(idb,id2)
                call getRadius(id1,rad1)
                call getRadius(id2,rad2)
                call getLum(id1,xlum1)
                call getLum(id2,xlum2)
                call get_ss_type(id1,ik1)
                call get_ss_type(id2,ik2)
                call get_bs_type(idb,ik0n)
                call getEcc(idb,ecc)
                call getSemi(idb,ap)
                call get_loid_mass(idb,sm1s)
                call get_hiid_mass(idb,sm2s)                
                print*,'im1,idb,id1,id2,ik1,ik2,ik1o,ik2o,sm1,sm2=',
     &          im1,idb,id1,id2,ik1,ik2,ik1o,ik2o,sm1s,sm2s
*
*       final set of bin(k,...) arrey - just to be sure that masses and
*       ids are the same for bin(k,...) and evolution
*
                bin(k,1) = sm1s/zmbar
                bin(k,2) = sm2s/zmbar
                body(im1) = (sm1s + sm2s)/zmbar
                bin(k,7) = 1.0d8*float(id1) + float(id2)                
*
*       check if semi major axis is greater than the sum of star radii
*       if not then the merger is due because of 3-body interaction.
*       Here only flag is set to distinguish between mergers due to binary
*       evolution and due to 3-body interaction.
*
                call binary_exists(idb,iexist)
                if(ik1o.le.1.and.ik2o.le.1) then
                  if(rad1+rad2.gt.ap.or.iexist.eq.0) then
                  if(sm1s+sm2s.gt.sturn.and.(ik1.le.1.or.ik2.le.1)) then
                      print*,'b3f  smsss,sturn,tim=',smsss,sturn,timevo
                      ibstra(im1) = 3
                      ibs3 = ibs3 + 1
                    endif  
                  endif
                endif
*
                print*,'k,im1,ik1,ik2,a,bin,ap,rad1,rad2,ibstra,ibs3=',
     &    k,im1,ik1,ik2,bin(k,3),bin(k,5),ap,rad1,rad2,ibstra(im1),ibs3
*
*        update time for the next evolution step
*
                call get_bs_updatetime(idb,uptime(im1))
*
*        check for a merger event
*                call binary_exists(idb,iexist)
*
                nob = 0
                if(iexist.eq.0) then
                  ikind(im1) = -3
                  names(im1) = idb
                  nameb(im1) = 0
                  call get_loid_mass(idb,sm1)
                  call get_hiid_mass(idb,sm2)
                  ap = 0.d0
                  ecc = 0.d0
                  if(ik2.eq.15) then
                    bin(k,2) = 0.d0
                    bin(k,1) = sm1/zmbar
                    names(im1) = id1
                    sm2 = 0.d0
                    sm2s = 0.d0
                    rad2 = 0.d0
                    xlum2 = 0.d0
                  endif
                  if(ik1.eq.15) then
                    bin(k,2) = sm2/zmbar
                    bin(k,1) = 0.d0
                    names(im1) = id2
                    sm1 = 0.d0
                    sm1s = 0.d0
                    rad1 = 0.d0
                    xlum1 = 0.d0
                  endif
                  bin(k,3) = 0.d0
                  bin(k,4) = 0.d0
                  body(im1) = bin(k,1) + bin(k,2)
                  nmerge = nmerge + 1
                  print*,'merger after binary interaction intb3f-1'
                  if(sm1.le.tolm.and.sm2.le.tolm) then
*
*      mark obliterated binary (massless) to removal from the system
*
                    body(im1) = 0.d0
                    bin(k,1) = 0.d0
                    bin(k,2) = 0.d0
                    iob = iob + 1
                    iobt = iobt + 1
                    nob = 1
                    deb = 0.d0
                    sm1s = 0.d0
                    sm2s = 0.d0
                    rad1 = 0.d0
                    rad2 = 0.d0
                    xlum1 = 0.d0
                    xlum2 = 0.d0
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
                    r(i) = 1.0d+15 + im1
                    xmax(im1) = 1.01d0*r(i)
                    print*,'intb3f-obl iob,im1,sm1,sm2=',iob,im1,sm1,sm2
                  endif
                endif
*
                if(dm3f.le.1.d-10) then
                  dmf3 = 0.d0
                else
                  ssevo = ssevo + dm3f/zmbar
                  slosev = slosev + dm3f
                  nmloev = nmloev + 2
                  nsevo = nsevo + 1
                  ehkin = 0.5d0*ssevo*(vr(im1)**2 + vt(im1)**2)
                  ehkin1 = ehkin1 + ehkin
                  ehmlev = ehmlev - ehkin
                  smt = smt - ssevo
                  write(6,9321) im1,nob,nsevo,dm3f,ssevo,ehkin
 9321             format(1x,
     &                 'intb3f-evo im1,nob,nsevo,dm3f,ssevo,ehkin=',
     &                  3I9,1p3e14.6)
                endif
*                  
                if(deb.eq.0.d0) go to 52 
*                
                sm3 = body(im1) + body(im2)
                debb = body(im2)*deb/sm3
                deb3 = body(im1)*deb/sm3
                ikind(im1) = -abs(ikind(im1))
                ikind(im2) = -abs(ikind(im2))
*
                call kick(i,debb)
                write(6,*) 'kick-b i,im1,debb,vr,vt,ebin,ikind = ',
     &          i,im1,debb,vr(im1),vt(im1),bin(k,5),ikind(im1)
                ii2 = ii
                call kick(ii2,deb3)
                write(6,*) 'kick-s i,im2,deb3,vr,vt,ikind,body = ',
     &          ii2,im2,deb3,vr(im2),vt(im2),ikind(im2),body(im2)
*
 52             continue
* 
                nb3fin = nb3fin + 1
                ehbin3 = ehbin3 + deb
*
*      store information about binary field star interaction
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
                binin(inumt,1) = r(i)/pctonbu
                binin(inumt,2) = deb
                binin(inumt,3) = body(im1)*zmbar
                binin(inumt,4) = body(im2)*zmbar
                binin(inumt,5) = timevo
                binin(inumt,6) = bin(k,5)
                iwww = 0
                if(iexist.eq.0) then
                  bin(k,5) = 0.d0
                  binin(inumt,7) = 4
                  iwww = 4
                else
                  binin(inumt,7) = 1
                endif
                if(iexch1.eq.1) binin(inumt,7) = 9 + iwww
*
                ibbb = int(binin(inumt,7) + 0.1d0)
                write(44,420) idb,id1,id2,im1,timevo,sm1oo,sm1s,sm2oo,
     &                 sm2s,rad1,rad2,xlum1,xlum2,ap,ecc,r(i),ik0,
     &                 ik0n,ik1oo,ik1,ik2oo,ik2,iexist,ibstra(im1),
     &                 ibbb
*
420             format(1x,'#',4i7,1p12e16.8,9i4)
*
              endif
            endif
*
 20      continue
*
 10   continue
*
      close(43)
      close(44)
*      
      if(ifint.eq.0.and.nfint.eq.0) return
*
      smt = smt
      n = n + ndistm
      nt = nt - iob
      call sort3(nt0,r,iname,inameo)
      nzst(nsuzon) = nt
      call coepot
      call energy(2)
      potnew = pot
      zkinn = zkin
      ehmlevn = ehmlev
      ehbin3n = ehbin3
      ekicktn = ekickt
      write(6,5432) zkino,zkinn,ehmlevo,ehmlevn,ehbin3o,ehbin3n,
     &              ekickto,ekicktn,zkinn-zkino,ehmlevn-ehmlevo,
     &              ehbin3n-ehbin3o,ekicktn-ekickto,(zkinn-zkino)-
     &              (ehmlevn-ehmlevo+ehbin3n-ehbin3o+ekicktn-ekickto)
 5432 format(1x,'intb3f- zkold,zknew,ehmlevo,ehmlevn,ekickto,ekicktn',
     &'ehbin3o,ehbin3n,Dzk,Dhm,Dhb3,Dhki,DDb3f =',1p13e14.6)
      enepot = enepot + potnew - potold
*
      print*,'intb3f nexch,nb3fin,ifint,nfint,iob,ndistm,n,nt,nt0,',
     &'potold,potnew,enepot,ssevo,ehkin1,iobx,iobt=',nexch,nb3fin,
     &ifint,nfint,iob,ndistm,n,nt,nt0,potold,potnew,enepot,
     &ssevo*zmbar,ehkin1,iobx,iobt
      etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
     &             ehmlev - ehcoll + eccoll - ekickt - enekin
      write(6,7321) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &              ekickt,enekin,escsta,ehbin3,ehb3b3
 7321 format(1x,'3bf- etotn,zkin,pot,ehmlev,enepot,ehcoll,eccol,',
     &'ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)
*   
      return
*
      end
*
*
*
