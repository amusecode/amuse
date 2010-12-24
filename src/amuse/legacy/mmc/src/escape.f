      subroutine escape(lmax,tim)
*
*
*       determine escapers and save their parameters, compute energy
*       -------------------------------------------------------------
*       changes due to escapers
*       -------------------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 e1,a1,sm1,desc,escr(100000),zpars(20),sscal,e1r,aursun,
     &       ecc,ecc0,epoch1,epoch2,rad1,rad2,radc1,radc2,rade1,rade2,
     &       rl1,rl2,s0mass1,s0mass2,scmass1,scmass2,semass1,semass2,
     &       semi0,slum1o,slum2o,sm1so,sm2so,smass1,smass2,spin1,spin2,
     &       ss1,ss2,t1,t2,tb0,tmstim1,tmstim2,tphys1,tphys2,tphysf1,
     &       tphysf2,vs(3),vss,yeardy,z,pctonbu,timevo,tim
*
      integer i,k,kk,im1,l,izesc,lmax,id,iks,id1,id2,ik0,ik1,ik2,kh,
     &        idb
*
*
*       check for escapers  -   energy > 0/-smt/rtid   or/and r > rtid
*
*
      print*,'lmax,dt = ',lmax,tim
      open(37,file='escape_binary.dat',access='append')
*
      aursun = 214.95d0
      pctonbu = rtidkg/rbar
      sscal = rbar/rtidkg/rsuntopc
*
      timevo = tim
*
      do 10 k = 1,nt
*
*    gard for cases in which the total number of objects in the system has 
*    changed because of obliterated objects formed in interactions
*
      if(lmax.gt.nt) then 
        print*,'lmax > nt   tim=', lmax,nt,tim
        go to 10
      endif
*
*     determine total energy and angular momentum of stars inside 'lmax' shell
*
      im1 = iname(k)
*
      if(xmax(im1).lt.r(k).and.xmax(im1).gt.0.d0
     &  .and.time.gt.0.d0) then
        write(6,1000) k,im1,xmax(im1),r(k),tim
 1000   format(1x,'xmaxxx k,im1,xmax(im1),r(k),time',2i8,1p3e16.8)
      endif
*
*     gard for cases in which second object taken part in interactions
*     intb3b3 or intb3f was choosen outside lmax (this can happen because
*     second binary has to have orbit which intersect with the orbin of
*     the first binary. Its position in the system can be elswere. For
*     intb3f the star for interaction can be choosen outside lmax if binary
*     is at lmax or there are a lot of binaries close to lmax). For 
*     interacting objects ikind is negative
*
c
c      Now every call of escape all stars are checked if they are due to
c      escape. Then stars which can be removed because of smaller Rtid
c      (evolution mass loss) will be properly treated
c
c      if(k.gt.lmax.and.ikind(im1).ge.0) go to 10
c
*        
      izesc = 0
      sm1 = body(im1)
      e1 = u(k) + 0.5d0*(vr(im1)**2 + vt(im1)**2)
      e1r = u(k) + 0.5d0*vr(im1)**2
      a1 = r(k)*vt(im1)
*
cChanged for M67:
c      if(imodel.eq.3.or.imodel.eq.4) then
      if(imodel.eq.3.or.imodel.eq.4.or.imodel.eq.5) then
        if(xmax(im1).gt.rtid.or.r(k).gt.rtid) izesc = 1
c        if(e1.gt.-xtid*smt/rtid.or.body(im1).eq.0.d0) izesc = 1
c        if(e1.gt.-1.25d0*smt/rtid.or.body(im1).eq.0.d0) izesc = 1
c        if(e1.gt.-1.5d0*smt/rtid.or.body(im1).eq.0.d0) izesc = 1
c        if(e1.gt.-smt/rtid.or.body(im1).eq.0.d0) izesc = 1
c        if(e1.gt.0.0d0.or.body(im1).eq.0.d0) izesc = 1
        if(xtid.eq.0) then
          if(e1.gt.0.0d0.or.abs(body(im1)).lt.tolm) izesc = 1
        else
c          xtid1 = 1.5d0 - 2.5d0*(log(gamma*nt)/nt)**0.25
c          if(e1r.gt.-xtid*smt/rtid.or.e1.gt.0.d0.or.
c     &       body(im1).eq.0.d0) izesc = 1
          if(e1.gt.-xtid*smt/rtid.or.abs(body(im1)).lt.tolm) izesc = 1
c          if(e1.gt.(-xtid1*smt/rtid + 0.5*(a1/rtid)**2).or.
c     &      e1.gt.0.d0.or. body(im1).eq.0.d0) izesc = 1               
        endif       
      else
        if(e1.gt.0.0d0.or.abs(body(im1)).lt.tolm) izesc = 1
      endif
*
      if(izesc.eq.1) then
*
        write(6,6712) k,im1,ikind(im1),lmax,nt,body(im1)*zmbar,xtid,
     &                rtid,smt,e1,r(k),xmax(im1),timevo 
 6712   format(1x,'k,im1,ikind,lmax,nt,body,xt,rtd,mt,e1,r,xm,time= '
     &        ,5i8,1p8e12.4)
*
        if(abs(body(im1)).lt.tolm) then
          nescm = nescm + 1
          print*,'k,im1,nescm,m,r = ',k,im1,nescm,body(im1),r(k)
        endif
*
        if(xescp(im1).lt.1.0d-10) then
          xescp(im1) = r(k)
          xesct(im1) = timevo
        endif
*
        nescst = nescst + 1
        escsta = escsta + e1*sm1 + 0.5d0*sm1*sm1/r(k)
        sloses = sloses + sm1
*
        if(abs(ikind(im1)).eq.1.or.abs(ikind(im1)).ge.3) then
*
*      check how many stars after SN natal kick are removed from the system
*
          if(nkick(im1).eq.1) ntsn1 = ntsn1 + 1
          if(nkick(im1).eq.2) ntsn2 = ntsn2 + 1
*                
          id = names(im1)
          call get_ss_type(id,iks)
          if(iks.ge.13) then
            mnsbh = mnsbh + 1
            print*,'im1,id,iks,ikind,mnsbh=',im1,id,iks,ikind(im1),mnsbh
          endif
        endif
*
        if(abs(ikind(im1)).eq.2) then
*
*      check how many binaries after SN natal kick are removed from the system
*
          if(nkick(im1).eq.1) ntsnb = ntsnb + 1
*
*      compute binary parameters needed to restart its evolution
*
          idb = nameb(im1)
          call get_loid_mass(idb,sm1so)
          call get_hiid_mass(idb,sm2so)
          call get_loid(idb,id1)
          call get_hiid(idb,id2)
          call get_bs_type(idb,ik0)
          call getLum(id1,slum1o)
          call getLum(id2,slum2o)
          call get_ss_type(id1,ik1)
          call get_ss_type(id2,ik2)
          call getEcc(idb,ecc)
          call get_mass(id1,ss1)
          call get_mass(id2,ss2)
          call getAge(id1,t1)
          call getAge(id2,t2)
          tphys1 = t1
          tphysf1 = t1
          tphys2 = t2
          tphysf2 = t2
          ik0 = ik0
          ik1 = ik1
          ik2 = ik2
          call getZ0(id1,z)
          call zcnsts(z,zpars)
          call getVel(id1,vs)
          vss = 0.d0
          call getEpoch(id1,epoch1)
          call getEpoch(id2,epoch2)
          call getMStime(id1,tmstim1)
          call getMStime(id2,tmstim2)
          call getMass0(id1,s0mass1)
          call getMass0(id2,s0mass2)
          smass1 = ss1
          smass2 = ss2
          call getMc(id1,scmass1)
          call getMc(id2,scmass2)
          call getMenv(id1,semass1)
          call getMenv(id2,semass2)
          call getRadius(id1,rad1)
          call getRadius(id2,rad2)
          call getRc(id1,radc1)
          call getRc(id2,radc2)
          call getRenv(id1,rade1)
          call getRenv(id2,rade2)
          slum1o = slum1o
          slum2o = slum2o
          call getSpin(id1,spin1)
          call getSpin(id2,spin2)
          call getRL(id1,rl1)
          call getRL(id2,rl2)
          ecc0 = ecc
          call getSemi(idb,semi0)
          aursun = 214.95d0
          yeardy = 365.24d0
          tb0 = (semi0/aursun)*SQRT(semi0/(aursun*
     &          (smass1 + smass2)))*yeardy
*
*      end of computation of binary parameters
*
          kk = nwhich(im1)
*          
          write(37,1010) im1,kk,idb,id1,id2,ikind(im1),timevo,z,
     &                   (zpars(kh),kh=1,20),vss,ecc0,semi0,tb0,ik0,
     &                   ik1,tphys1,tphysf1,epoch1,tmstim1,smass1,
     &                   s0mass1,scmass1,semass1,rad1,radc1,rade1,
     &                   slum1o,spin1,rl1,ik2,tphys2,tphysf2,epoch2,
     &                   tmstim2,smass2,s0mass2,scmass2,semass2,rad2,
     &                   radc2,rade2,slum2o,spin2,rl2,bin(kk,3)*sscal,
     &                   bin(kk,4)
 1010     format(1x,6i8,1p26e14.6,2i4,1p14e14.4,i4,1p16e14.4)
*
        endif
*
        if(ikind(im1).eq.2) then       
*
          if(bin(kk,6).eq.0.d0) then
*
*  negative time marks that binary escape due to relaxation process
*
            bin(kk,6) = -timevo
          else
            bin(kk,6) = -bin(kk,6)
          endif
*
          write(6,*) ' relaxation escape of binary ',kk,im1
          call flush(6)
*
          nescb3 = nescb3 + 1
          erelb3 = erelb3 + e1*sm1
          erb3in = erb3in + bin(kk,5)
        endif
*
        if(ikind(im1).lt.0) then
*
          write(6,*) 'interaction escape',im1,ikind(im1)
          call flush(6)
*
          if(ikind(im1).eq.-1.or.ikind(im1).le.-3) then
            nesb3s = nesb3s + 1
            escb3s = escb3s + e1*sm1
            slob3s = slob3s + sm1
          endif
*
          if(ikind(im1).eq.-2) then
            nescb3 = nescb3 + 1
            escbi3 = escbi3 + e1*sm1
            escbb3 = escbb3 + bin(kk,5)
            slob3b = slob3b + sm1
*
*     negative number on the global list marks that object escape
*     due to interaction processes
*
            bin(kk,7) = -bin(kk,7)
          endif
        endif
*
        iesc = iesc + 1
        iesc1 = iesc1 + 1
        iename(iesc) = im1
        ienam1(iesc1) = im1
*
*     check if a blue straggler is escaping
*
        if(ibstra(im1).gt.0.and.ibstra(im1).le.4) ibstra(im1) = 0
*
        if(abs(ikind(im1)).eq.2.or.abs(ikind(im1)).eq.3) then
           ienam1(iesc1) = nameb(im1)
        endif
*
        iekind(iesc1) = ikind(im1)
        escmas(iesc1) = sm1
        escene(iesc1) = e1*sm1
        escdis(iesc1) = xescp(im1)
        escang(iesc1) = a1*sm1
        esctim(iesc1) = xesct(im1)
        escr(iesc) = r(k)
        r(k) = 1.0d8 + float(im1)
      else
        if(ikind(im1).lt.0) ikind(im1) = -ikind(im1)
        if(xescp(im1).gt.1.0d-10) then
          xescp(im1) = 0.0d0
          xesct(im1) = 0.0d0
        endif
        if(nkick(im1).ge.1) nkick(im1) = 0
      endif
*
 10   continue
*
      close(37)
*
*     reduce the total mass and the total number of stars
*
      desc = 0.0d0
      do 80 i=1,iesc
         k = iename(i)
         smt = smt - body(k)
         nt = nt - 1
         if(iesc.eq.1) go to 80
         do 90 l = i+1,iesc
            kk = iename(l)
            desc = desc + body(k)*body(kk)/escr(l)
 90      continue
 80   continue
      print *,'escape: time,iesc,desc ',time,iesc,desc
*
      escsta = escsta + desc
      nzst(nsuzon) = nt
cChanged for M67:
      if(imodel.ge.3) rtid = rtidkg*smt**0.33333333
*
      return
*
      end
*
*
*
