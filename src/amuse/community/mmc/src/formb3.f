*
*
      subroutine formb3(nup)
*
*      compute formation of binaries due to three-body interactions
*      ------------------------------------------------------------
*      Goodman & Hut 1994, .....
*      -------------------------
*
*     abs(ikind(im)) - formb3 is the last routine in the interaction sequence
*     (intb3b3 - intb3f - formb3) therefore same objects could take place 
*     already in interactions in intb3b3 or intb3f and have negative ikind
*
*
      include 'common.h'
*
      real*8 dt,co,vol,velave,smave,den,vb,f2,wrr,ptold,ptnew,
     &       xmm,xmm1,xmm2,age1,age2,sm1s,ehpot,ehkin,cor,timea,anb,
     &       body12,tform,tb3f,tb3b3,tcoll,tim,sm2s,ekicks1,ekicks2,
     &       ekicks3,sm1so,sm2so,zkold,zknew,ehmlevo,ehmlevn,ehbin3o,
     &       ehbin3n,ekickto,ekicktn,r1,r2,xlum1,xlum2,slum
*
      real*4 ran2
*
      integer nup,l,im1,im2,im3,ibin,lmin,lmax,n,inum,inumt,idb,
     &        id1,id2,is1,is2,is3,ik1,ik2,ik3,im,k1,k11,k2,k3,iexi,
     &        iexi1,iexi2,noibin,ik0
*
      common /timset/ tform(20),tb3f(20),tb3b3(20),tcoll(20)
*
      open(44,file=trim(datadir)//'/binarym.dat',access='append')
      open(43,file=trim(datadir)//'/starm.dat',access='append')

      n = nt
      sm2s = 0.d0
      f2 = 0.9d0
      ekisks1 = 0.d0
      ekicks2 = 0.d0
      ekicks3 = 0.d0
      pctonbu = rtidkg/rbar
      iexi = 0
      noibin = 0
*
*     Stodolkiewicz's density definition
*     cos = 4/3*pi*3^4.5/(4/3*pi)^3*f2*float(n)/log(gamma*n)*(2/3)^3
*
*     Henon's density definition
*     coh = 4/3*pi*3^4.5/(2*pi)^3*float(n)/log(gamma*n) = cos
*
      if(alphal.eq.-1.0d0) then
*
*     coefficient for single-mass case
*
        co = 2.36916d0*f2*float(nt)/log(gamma*nt)
      else
*
*     coefficient for multi-mass case
*
*     co = 3^6*f2*(m1*m2)^4*m3^2.5*nt0/(2^1.5*pi^2*((m1+m2)*(m1+m2+m3))^0.5*
*          <m>^4.5*(v^2)^4.5*log(gamma*nt0))
*
        co = 26.114565d0*f2*float(nt)/log(gamma*nt)
      endif
*
      ibin = 0
*
      do 10 l = 1,nup
*
         if(l.eq.1) then
           lmin = 1
         else
           lmin  = nzst(l-1) + 1
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
         tform(l) = tform(l) + dt*tscale0*log(gamma*nt00)
     &              /log(gamma*ntnew)*float(ntnew)/float(nt00)
         print*,'formb3- l,nup,tform(l) = ',l,nup,tform(l)
         call flush(6)
         tim = tform(l)
*
*      acomulate probability 'pbin3' for formation of 3-body binaries
*      find stars to form a binary
*
         k1 = lmin - 1
*
 20      continue
*
         k1 = k1 + 1
c         if(k1.gt.int(0.8*n)) go to 10
         if(k1.gt.lmax) go to 10 
*
         im1 = iname(k1)
         if(abs(ikind(im1)).eq.2.or.names(im1).eq.0) go to 20
         if(abs(body(im1)).lt.tolm) go to 20
         if(nkick(im1).ge.1) go to 20
         if(r(k1).gt.1.d8) go to 20
         is1 = names(im1)
         call getLum(is1,slum)
         if(isNaN(slum)) then
           print*,'k1,im1,is1,lum = ',k1,im1,is1,slum
           go to 20
         endif
*
         k2 = k1
*
 30      continue
*
         k2 = k2 + 1         
c         if(k2.gt.int(0.8*n)) go to 10
         if(k2.gt.lmax) go to 10
*
         im2 = iname(k2)
         if(abs(ikind(im2)).eq.2.or.names(im2).eq.0) go to 30
         if(abs(body(im2)).lt.tolm) go to 30
         if(nkick(im2).ge.1) go to 30
         if(r(k2).gt.1.d8) go to 30
         is1 = names(im2)
         call getLum(is1,slum)
         if(isNaN(slum)) then
           print*,'k2,im2,is2,lum = ',k2,im2,is1,slum
           go to 30
         endif
*                                    
         k3 = k2
*
 40      continue
*
         k3 = k3 + 1
c         if(k3.gt.int(0.8*n)) go to 10
         if(k3.gt.lmax) go to 10
*
         im3 = iname(k3)
         if(abs(ikind(im3)).eq.2.or.names(im3).eq.0) go to 40
         if(abs(body(im3)).lt.tolm) go to 40
         if(nkick(im3).ge.1) go to 40
         if(r(k3).gt.1.d8) go to 40
         is1 = names(im3)
         call getLum(is1,slum)
         if(isNaN(slum)) then
           print*,'k3,im3,is3,lum = ',k3,im3,is1,slum
           go to 40
         endif
*
*       take care for multi-mass systems
*
         if(alphal.eq.-1.0d0) then
           if(k1.eq.1) then
             vol = co*r(k1+2)**3
           else
             vol = co*(r(k1+2)**3 - r(k1-1)**3)
           endif
           call veloci(k1,velave,smave)
           call densi(k1,den)
*
           pbin3 = vol*smave**5*den**3*dt/velave**4.5
*
         else
           call veloci(k1,velave,smave)
           call densi(k1,den)
cAdded DCH 28/8/6 to deal with a central triple of black holes
           if(k3.ge.3.and.k3.le.5) then
             is1 = names(im1)
             call get_ss_type(is1,ik1)
             if(ik1.ge.10.and.ik1.lt.15) then
               is2 = names(im2)
               call get_ss_type(is2,ik2)
               if(ik2.ge.10.and.ik2.lt.15) then
                 is3 = names(im3)
                 call get_ss_type(is3,ik3)
                 if(ik3.ge.10.and.ik3.lt.15) then
                   velave = body(im1)*(vr(im1)**2 + vt(im1)**2) +
     &                      body(im2)*(vr(im2)**2 + vt(im2)**2) +
     &                      body(im3)*(vr(im3)**2 + vt(im3)**2)      
                   smave = body(im1) + body(im2) + body(im3)
                   velave = velave/smave
                   smave = smave/3.d0
                   den = 3.d0/(r(k3)**3 - r(k1)**3)
                   print*,'Three/or more inmost stars are WD-NS-BH',
     &                       body(im2),body(im2),body(im3),den,smave,
     &                       velave,co,dt,k3
                 endif
               endif
             endif
           endif
           if(body(im1)*body(im2)*body(im3).ne.0.d0) then
             xmm = (body(im1)*body(im2))**4*body(im3)**2.5
             xmm1 = body(im1) + body(im2)
             xmm2 = xmm1*(xmm1 + body(im3))
             xmm = xmm/(sqrt(xmm2)*(smave*velave)**4.5)
             pbin3 = co*xmm*dt*den**2
           else
             pbin3 = 0.d0
           endif
*
         endif
*
*      new check for binary formations starts from k3 + 1
*
         k11 = k1
         k1 = k3
         wrr = ran2(irun)
         if(wrr.gt.pbin3) go to 20
*
*         compute total potential before binary formation to deal with
*         potential changes due to substituting two single stars by binary
*
         if(ibin.eq.0.and.noibin.eq.0) then
           call energy(2)
           ptold = pot
           zkold = zkin
           ehmlevo = ehmlev
           ehbin3o = ehbin3
           ekickto = ekickt
         endif
*
*   evolve single stars to the time tim - time of binary formation
*
         call mloss_single(k11,tim,ekicks1,iexi)
*
         call mloss_single(k2,tim,ekicks2,iexi1)
*
         call mloss_single(k3,tim,ekicks3,iexi2)
*
*        check for obliterated stars
*        
c         if(r(k11).gt.1.d8.or.r(k2).gt.1.d8.or.r(k3).gt.1.d8) then
c           print*,'formb3-obl  k1,k2,k3,im1,im2,im3,sm1,sm2,sm3 =',
c     &     k1,k2,k3,im1,im2,im3,body(im1)*zmbar,body(im2)*zmbar,
c     &     body(im3)*zmbra
c           go to 20
c         endif  
*            
         if(iexi.eq.0.or.iexi1.eq.0.or.iexi2.eq.0) then
           noibin = noibin + 1
           print*,'formb3-no  im1,im2,im3,iexi,iexi1 iexi2 = ',
     &             im1,im2,im3,iexi,iexi1,iexi2
           go to 20
         endif
*                                                             
c         tim = tim + 2.d-14
         pbin3 = 0.0d0
         nbin3 = nbin3 + 1
         ik0 = -1
*
         print*,'k1,velave,smave,den =',k11,velave,smave,den
         call veloc(k11,velave,smave)
         print*,'k1,velave,smave,den =',k11,velave,smave,den
*
*       define parameters of binary and compute kick velocities
*       for binary and field star
*
         vb = smave*velave
*
         ibin = ibin + 1
         ikind(im1) = 2
*
         id1 = names(im1)
         id2 = names(im2)
*
*      Update names and nameb and prepare the new binary for binary 
*      evolution
*      nwhich contains the name of binary in the binary list.
*      if ibx=im1-nss0 is negative then binary index can be get from the
*      nwhich(im1)
*
*     set ibstra for newly formed binary. If ibstra(im1) or/and ibstra(im2)
*     is greater than 0 then ibstra(im1) is set to ibstra(im1) or/and
*     ibstra(im2)
*
         if(ibstra(im1).ge.ibstra(im2)) then
           ibstra(im1) = ibstra(im1)
           ibstra(im2) = 0
         else
           ibstra(im1) = ibstra(im2)
           ibstra(im2) = 0
         endif
*
         nameb(im1) = id1
         nameb(im2) = 0
         nwhich(im1) = nbin3
         nbinar(nbin3) = im1
         call get_mass(id1,sm1s)
         call get_mass(id2,sm2s)
         idb = nameb(im1)
         body12 = sm1s + sm2s
         call getAge(id1,age1)
         call getAge(id2,age2)
*
         print*,'k1,k2,k3,im1,im2,im3,id1,id2,idb,nameb,names1,',
     &   'names2,sm1b,sm2b,sm1s,sm2s body12,age1,age2 =',k11,k2,k3,im1,
     &   im2,im3,id1,id2,idb,nameb(im1),names(im1),names(im2),
     &   body(im1)*zmbar,body(im2)*zmbar,sm1s,sm2s,body12,age1,age2
*
         names(im1) = 0
         names(im2) = 0
*
*       change binary parameters in the case when one or two stars
*       change their masses because of stellar evolution (age1=age2)
*       
         body(im1) = sm1s/zmbar
         body(im2) = sm2s/zmbar
         body12 = body(im1) + body(im2)
*                      
         bin(nbin3,3) = 0.5d0*body(im1)*body(im2)/vb
         bin(nbin3,4) = sqrt(ran2(irun))
         bin(nbin3,5) = vb
         bin(nbin3,6) = time
         bin(nbin3,7) = 1.0d8*float(id1) + float(id2)
*             
         anb = bin(nbin3,3)*rbar/rtidkg/rsuntopc
*
      print*,'formb3-for im1,im2,id1,id2,idb,nb3,sm1,sm2,anb,ecc,tim=',
     &       im1,im2,id1,id2,idb,nbin3,sm1s,sm2s,anb,bin(nbin3,4),tim
*
*      create a new binary
*
         call create_binary(idb,id1,id2,anb,bin(nbin3,4),tim)
*
*       check if binary exist after formation
*       
         call binary_exists(idb,iexist)
         if(iexist.eq.0) then
           print*,'formb3-noexist im1,idb,iexist',im1,idb,iexist
         endif
*           
         call get_loid(idb,id1)
         call get_hiid(idb,id2)
         call get_loid_mass(idb,sm1so)
         call get_hiid_mass(idb,sm2so)
         call get_radius(id1,r1)
         call getLum(id1,xlum1) 
         call get_radius(id2,r2)
         call getLum(id2,xlum2)
         call get_bs_type(idb,ik0n)
         call get_ss_type(id1,ik1n)
         call get_ss_type(id2,ik2n)
         ik0n = -1
         ik1 = -1
         ik2 = -1
*
         write(44,420) idb,id1,id2,im1,tim,sm1so,sm1s,sm2so,sm2s,
     &                 r1,r2,xlum1,xlum2,anb,bin(nbin3,4),r(k11),ik0,
     &                 ik0n,ik1,ik1n,ik2,ik2n,iexist,ibstra(im1),0
*
 420               format(1x,'#',4i7,1p12e16.8,9i4)
*                         
         bin(nbin3,1) = sm1so/zmbar
         bin(nbin3,2) = sm2so/zmbar
         bin(nbin3,7) = 1.0d8*float(id1) + float(id2)      
c
      print*,'im1,im2,id1,id2,idb,ikind,sm1so,sm1s,sm2so,sm2,body(im1)'
     &,',body(im2),bin(nb,1),bin(nb,2),body12 =',im1,im2,id1,id2,idb,
     &ikind(im1),sm1so,sm1s,sm2so,sm2s,body(im1)*zmbar,body(im2)*zmbar,
     &sm1s,sm2s,body12*zmbar
*
*      store information about binary field star interaction
*
         iinte3(nbin3) = iinte3(nbin3) + 1
         inum = iinte3(nbin3)
*
         if(inum.gt.50) then
           iinte3(nbin3) = 1
           inum = 1
         endif
*
         r(k11) = (body(im1)*r(k11) + body(im2)*r(k2))/body12
         bin(nbin3,8) = r(k11)
*
         inumt = 50*(nbin3 - 1) + inum
         binin(inumt,1) = r(k11)/pctonbu
         binin(inumt,2) = vb
         binin(inumt,3) = body12*zmbar
         binin(inumt,4) = body(im3)*zmbar
         binin(inumt,5) = tim
         binin(inumt,6) = bin(nbin3,5)
         binin(inumt,7) = 0
*
*      determine energy generated during binary formation and split 
*      this energy for the binary and the star
*
         print*,'vr(im1),vt(im1),vr(im2),vt(im2) = ',vr(im1),
     &          vt(im1),vr(im2),vt(im2)
         vb = 0.5d0*(body(im1)*(vr(im1)**2 + vt(im1)**2) + 
     &        body(im2)*(vr(im2)**2 + vt(im2)**2) + smave*velave)
         body(im1) = body(im1) + body(im2)
c            vr(im1) = sqrt(smave*velave/(3.0d0*body12))
c            vt(im1) = sqrt(2.0d0*smave*velave/(3.0d0*body12))
         vr(im1) = sqrt(smave*velave/(3.0d0*body(im1)))
         vt(im1) = sqrt(2.0d0*smave*velave/(3.0d0*body(im1)))
c            nameb(im1) = nbin3 
         r(k2) = 1.0d8 + float(im2)
         xmax(im2) = 1.01d0*r(k2)
*
         print*,'vb,smave,velave,body(im1),body(im2),body(im3),',
     &   'vr(im1),vt(im1),vr(im2),vt(im2) = ',vb,smave,velave,
     &   body(im1)*zmbar,body(im2)*zmbar,body(im3)*zmbar,vr(im1),
     &   vt(im1),vr(im2),vt(im2)
            vb = vb*body(im3)/(body(im1) + body(im3))
         call kick(k11,vb)
         write(6,*) 'kick-b-k1 k1,im1,debb,m12=',k11,im1,vb,
     &               body(im1)*zmbar
         vb = vb*body(im1)/body(im3)
         call kick(k3,vb)
         write(6,*) 'kick-s-k3  k3,im3,deb3,m3=',k3,im3,vb,
     &              body(im3)*zmbar
*
         ehbin3 = ehbin3 + bin(nbin3,5)
*
         go to 20
*
 10   continue
*
*      compute new potential and sort stars after removing from the
*      calculation star which is now one of the binary components (centre 
*      of mass of newly formed binary is substitute for position of 
*      another star which form binary)
*
      close(43)
      close(44)
*
      if(ibin.eq.0.and.noibin.eq.0) return
*
      call sort3(n,r,iname,inameo)
      nt = nt - ibin
      smt = smt
      nzst(nsuzon) = nt
      call coepot
      call energy(2)
      ptnew = pot
      zknew = zkin
      ehmlevn = ehmlev
      ehbin3n = ehbin3
      ekicktn = ekickt
      enepot = enepot + ptnew - ptold
      write(6,5432) zkold,zknew,ehmlevo,ehmlevn,ehbin3o,ehbin3n,   
     &              ekickto,ekicktn,zknew-zkold,ehmlevn-ehmlevo,
     &              ehbin3n-ehbin3o,ekicktn-ekickto,(zknew-zkold)-
     &              (ehmlevn-ehmlevo+ehbin3n-ehbin3o+ekicktn-ekickto)
 5432 format(1x,'formb3  zkold,zknew,ehmlevo,ehmlevn,ekickto,ekicktn',
     &',ehbin3o,ehbin3n,Dzk,Dhm,Dhb3,Dhki,DDform =',1p13e14.6) 
*
      etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -  
     &         ehmlev - ehcoll + eccoll - ekickt - enekin
      write(6,6543) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &               ekickt,enekin,escsta,ehbin3,ehb3b3
 6543 format(1x,'formb3-f etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)
*                                            
      enekin = enekin + zknew - zkold
*
      print*,'formb3-fin  ibin,noibin,zko,zkn,poto,potn, = ', 
     &  ibin,noibin,zkold,zknew,ptold,ptnew,enepot,enekin
      enekin = enekin - zknew + zkold
*
      return
      end
*
*

      logical function isnan(x)
      double precision x
      if (.not.x.le.0.d0.and..not.x.gt.0.d0) then
         isnan = .true.
      else
         isnan = .false.
      endif
      return
      end
