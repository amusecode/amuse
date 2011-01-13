*
      subroutine mloss_single(i,tphys,ekick,iexist)
*
*
*       Mass loss from evolving a single star or a binary.
*       ------------------------------
*     In this version, stellar and binary evolution are carried out via the
*     interface of McScatter. Modified 26/4/6 so that mass is lost only if
*     the update time is passed.
*
*     mloss_single is called inside interactions subroutines, so ikind
*     can be negative. 
*     
*     ekick - kinetic energy added to the kinetic energy of an object 
*     because of natal kick due to spernove explosion. In interaction
*     subroutines ekick is sumed up to ekickt. 
*     Remember that at the and of each interactn/relaxation loop escape is
*     called, so all kicked object are removed from the system. Kinetic
*     energy of escaping stars is properly taken into account.
*
*     After each call of mloss_single in the interaction subroutine
*     one has to check if binary exist. If not the interaction has to be
*     aborded.       
*
      include 'common.h'
*
*
      real*8  timevo,sm1,sm2,sm1s,sm2s,sm1o,sm2o,a,anb,ssevo,cor,ebb,
     &        abb,pctonbu,ehpot,ehkin,r1,r2,xlum1,xlum2,tlog,smlog,
     &        smsss,sturn,z,zpars(20),ecc,sm1b,sm2b,sm1so,sm2so,slum1o,
     &        slum2o,xxm1,xxm2,xxl1,xxl2,potold,tphys,ssme,tphys1,
     &        tphys2,tphysf1,tphysf2,epoch1,epoch2,tmstim1,tmstim2,
     &        smass1,smass2,s0mass1,s0mass2,scmass1,scmass2,semass1,
     &        semass2,rad1,rad2,radc1,radc2,spin1,spin2,rl1,rl2,tb0,
     &        semi0,ecc0,rade1,rade2,vs(3),xxn,vscal,ekick,ekick1,
     &        ekick2,vst2,vst,vrold,vtold,vss,aursun,yeardy,cfr,potnew,
     &        ss1,ss2,t1,t2,upt,ehmlev1,semi,ekickbs,ekickbd,smbb,
     &        vrold2,vtold2
*
      integer n,i,k,im1,ibx,ik1,ik2,kk,ik0,ipoi,id3,id4,ievo,idb,ids,
     &        id1,id2,ikio,ikin,iexist,iblue,ik0n,inum,inumt,iob,ik1n,
     &        ik2n,ndist,iime,ikick,ixkick,ikickbs,ixkickbs,iexistbs,
     &        ikickbd,ixkickbd,iexistbd
*
      real*8 vs2(3)
      common /companion/ vs2 
*
*     scale velocity from km/s to N-body units
*
      r2 = 0.d0
      r1 = 0.d0
      xlum1 = 0.d0
      xlum2 = 0.d0
      smbb = 0.d0
      xxn = zmbar*rtidkg/rbar
      vscal = 0.06558d0*sqrt(xxn)
      vscal = 1.d0/vscal            
      pctonbu = rtidkg/rbar
      ekick = 0.d0
      ekickbs = 0.d0
      ekickbd = 0.d0
      ikick = 0
      ikickb = 0
      iexist = 1
      iexistbs = 1
      iexistbd = 1
      cor = 0.d0
      n = nt
      ssme = 0.d0
      iime = 0
      iob = 0
      ndist = 0
*
*       check the stellar evolution time for stars in a superzone
*
*       upgrade time for each zone
*
      timevo = tphys
      tlog = log10(1.d6*timevo)
*
      print*,'mloss_single i,im1,timevo = ',i,iname(i),timevo
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
      print*,'smsss,sturn,z0=',smsss,sturn,z
*
      ixkick = 0
      ixkickbs = 0
      ixkickbd = 0
*
      vs(1) = 0.d0
      vs(2) = 0.d0
      vs(3) = 0.d0
*
      ievo = 0
      ssevo = 0.d0
      im1 = iname(i)
*
      if(abs(ikind(im1)).eq.2) then
*
*          take care for ojects which are binaries
*
        oldtime(im1) = timevo
        ibx = nwhich(im1)
        idb = nameb(im1)
        ebb = bin(ibx,5)
        abb = bin(ibx,3)
*                                  
        if(bin(ibx,7).lt.0.0d0.or.bin(ibx,6).lt.0.0d0) then
          write (6,*) 'bin(7) or bin(6) < 0 in mloss'
          stop
        endif
*
        sm1b = bin(ibx,1)
        sm2b = bin(ibx,2)
        call get_loid_mass(idb,sm1so)
        call get_hiid_mass(idb,sm2so)
        sm1 = sm1so/zmbar
        sm2 = sm2so/zmbar
        call get_loid(idb,id1)
        call get_hiid(idb,id2)
        xxm1 = abs((sm1 - sm1b)/sm1) 
        xxm2 = abs((sm2 - sm2b)/sm2)
        if(xxm1.gt.1.0d-5.or.xxm2.gt.1.d-5) 
     & print*,'wrong-b-sing im1,idb,id1,id2,sm1,',
     &'sm1b,sm2,sm2b,sm1+sm2,sm1b+sm2b= ',im1,idb,id1,id2,sm1*zmbar,
     &sm1b*zmbar,sm2*zmbar,sm2b*zmbar,(sm1+sm2)*zmbar,(sm1b+sm2b)*zmbar
        call get_bs_type(idb,ik0)
        call getLum(id1,slum1o)
        call getLum(id2,slum2o)
        call get_ss_type(id1,ik1)
        call get_ss_type(id2,ik2)
        call getEcc(idb,ecc)
        call get_mass(id1,ss1)
        call get_mass(id2,ss2)
*
*      compute binary parameters needed to restart its evolution
*
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
        call setVel(id1,vs)
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
     &        (smass1 + smass2)))*yeardy
*
*      end of computation of binary parameters
*
        call ev_binary(idb,timevo)
*
        call get_bs_type(idb,ik0n)
        call binary_exists(idb,iexist)
*
        call get_loid(idb,id1)
        call get_hiid(idb,id2)
        call get_loid_mass(idb,sm1s)
        call get_ss_type(id1,ik1n)
        call get_hiid_mass(idb,sm2s)
        call get_ss_type(id2,ik2n)
        call get_mass(id1,ss1)
        call get_mass(id2,ss2)
        call getEcc(idb,ecc)
        call getSemi(idb,semi)
        call getVel(idb,vs)
        write(6,7654) im1,idb,id1,id2,iexist,ik0n,ik0,ik1,ik1n,ik2,
     &                ik2n,sm1so,sm1s,sm2so,sm2s,ecc0,ecc,semi0,semi
 7654   format(1x,'im1,idb,id1,id2,iexist,ik0n,ik0,ik1,ik1n,ik2,',
     &'ik2n,sm1so,sm1s,sm2so,sm2s,ecc0,ecc,semi0,semi=',4i9,7i4,
     &1p8e14.6)
*                                                   
        ievo = 1
*
*      first take care about binaries which survive
*
        if(iexist.ne.0) then
*
*     check for supernove formation and NS/BH kick
*
          smbb = (sm1s + sm2s)/zmbar
*          
          if((ik1n.ge.13.and.ik1n.lt.15.and.ik1.lt.13).or.
     &       (ik2n.ge.13.and.ik2n.lt.15.and.ik2.lt.13)) then
            if((ik1n.eq.13.or.ik2n.eq.13).and.iflagns.le.1) then
              vs(1) = 0.d0
              vs(2) = 0.d0
              vs(3) = 0.d0
              go to 140
            endif
            if((ik1n.eq.14.or.ik2n.eq.14).and.iflagbh.le.1) then
               vs(1) = 0.d0
               vs(2) = 0.d0
               vs(3) = 0.d0
               go to 140
            endif          
            nkick(im1) = 1
            iexistbs = 0
            ikickbs = ikickbs + 1
            ikicktbs = ikicktbs + 1
            vrold = vr(im1)
            vtold = vt(im1)
            vst2 = vs(2)**2 + vs(3)**2
            vst = sqrt(vst2)
            ixkickbs = 1
            ekick1 = 0.5d0*smbb*vscal**2*(vs(1)**2 + vst2)
            ekick2 = smbb*vscal*(vrold*vs(1) + vtold*vst)
            ekickbs = ekickbs + ekick1 + ekick2
            ekicktbs = ekicktbs + ekick1 + ekick2
            vr(im1) = vr(im1) + vscal*vs(1)
            vt(im1) = vt(im1)  + vscal*vst
            write(6,8712) timevo,im1,ik1n,ik2n,ikickbs,ikicktbs,ekick1,
     &                    ekick2,ekickbs,ekicktbs,vscal,vs(1),vs(2),
     &                  vs(3),vr(im1),vt(im1),iflagns,iflagbh
 8712       format(1x,'timevo-sin-bs,im1,ik1n,ik2n,ikickbs,ikicktbs,',
     &       'ekick1,ekick2,ekickbs,ekicktbs,vscale,vs(1),vs(2),vs(3),',
     &       'vr(im1),vt(im1),iflagns,iflagbh = ',1pe12.4,5i8,1p10e12.4,
     &        2i4)     
          endif
*
 140      continue
*
          sm1o = sm1
*
          sm1 = sm1s/zmbar
*
          nmloev = nmloev + 1
          ssevo = sm1o - sm1
          slosev = slosev + ssevo
*
          sm2o = sm2
*
          sm2 = sm2s/zmbar
*
          nmloev = nmloev + 1
          slosev = slosev + (sm2o - sm2)
          ssevo = ssevo + (sm2o - sm2)
*
*       modify binary parameters and center of mass parameters
*
          call get_bs_updatetime(idb,uptime(im1))
*                 
          if(ssevo.gt.1.0d-10) then
            print*,'evo-bin: ,im1,idb,id1,id2,ibx,nt0,nameb,',
     &'ik0,ik0n,ik1,ik1n,ik2n,sm1so,sm1s,sm2so,sm2s,sevo,uptim= ',
     & im1,idb,id1,id2,ibx,nt0,nameb(im1),ik0,ik0n,ik1,ik1n,ik2,
     & ik2n,sm1so,sm1s,sm2so,sm2s,ssevo*zmbar,uptime(im1)
          endif
*
          bin(ibx,7) = 1.0d8*float(id1) + float(id2)
          body(im1) = sm1 + sm2
          bin(ibx,1) = sm1
          bin(ibx,2) = sm2
          call get_sma(idb,a)
          call getEcc(idb,ecc)
          bin(ibx,4) = ecc
*
*            Here a is in solar radii
*
          anb = a*rsuntopc*pctonbu
          bin(ibx,3) = anb
          bin(ibx,5) = sm1*sm2/2.d0/bin(ibx,3)
          call get_radius(id1,r1)
          call getLum(id1,xlum1)
          call get_radius(id2,r2)
          call getLum(id2,xlum2)
*
*      mark obliterated binary (massless) to removal from the system
*
          if(sm1.lt.tolm.and.sm2.lt.tolm) then
            body(im1) = 0.d0
            bin(ibx,1) = 0.d0
            bin(ibx,2) = 0.d0
            iob = iob + 1
            iobt = iobt + 1
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
            r(i) = 1.0d+15 + im1
            xmax(im1) = 1.01d0*r(i)
            print*,'iob-ms,im1,sm1,sm2=',iob,im1,sm1,sm2
          endif
*
*     print only when binary component masses and luminosities change
*     because of binary evolution by more than 1%
*
          xxm1 = abs((sm1s - sm1so)/sm1so)
          xxm2 = abs((sm2s - sm2so)/sm2so)
          xxl1 = abs((xlum1 - slum1o)/slum1o)
          xxl2 = abs((xlum2 - slum2o)/slum2o)
c          if(xxm1.gt.0.01d0.or.xxm2.gt.0.01d0.or.xxl1.gt.0.01d0
c     &      .or.xxl2.gt.0.01d0) then
*
            write(44,420) idb,id1,id2,im1,timevo,sm1*zmbar,
     &                    sm1o*zmbar,sm2*zmbar,sm2o*zmbar,r1,r2,
     &                    xlum1,xlum2,a,ecc,r(i),ik0,ik0n,ik1,ik1n,
     &                    ik2,ik2n,iexist,ibstra(im1),-1
*
 420        format(1x,'#',4i7,1p12e16.8,9i4)
*
          go to 120
*
        endif
*
*      first take care about binary mergers
*
        if((iexist.eq.0.and.ik0n.eq.4).or.
     &     (iexist.eq.0.and.ik0n.eq.5)) then
*
      write(6,4390) im1,idb,id1,id2,ik0,ik1,ik2,sm1*zmbar,sm2*zmbar,
     & semi0,ecc0,rad1,rad2,rl1,rl2,t1,t2,epoch1,epoch2,ik0n,ik1n,
     & ik2n,sm1s,sm2s
 4390 format(1x,'evo-merged-init im1,idb,id1,id2,ik0,ik1,ik2,sm1,sm2,',
     &'a,e,r1,r2,rl1,rt2,t1,t2,ep1,ep2 = ',7i8,1p12e14.6,2x,
     &'ik0n,ik1n,ik2n,sm1s,sm2s = ',3i6,1p2e14.6)      
*
          ikind(im1) = 3
          ibx = nwhich(im1)
*
*      print binary parameters which merge and form a black hole
*      initial types of binary components: bh - wd, bh - ns, bh - bh
*
          if(ik1.eq.14.and.(ik2.ge.10.and.ik2.le.14)) then
            open(27,file=trim(datadir)//'/bhmergers.dat',
     &  access='append')
            write(27,1000) im1,ibx,idb,id1,id2,abs(ikind(im1)),
     &  timevo,z,(zpars(kk),kk=1,20),vss,ecc0,semi0,tb0,ik0,ik1,tphys1,
     &  tphysf1,epoch1,tmstim1,smass1,s0mass1,scmass1,semass1,rad1,
     &  radc1,rade1,slum1o,spin1,rl1,ik2,tphys2,tphysf2,epoch2,tmstim2,
     &  smass2,s0mass2,scmass2,semass2,rad2,radc2,rade2,slum2o,spin2,rl2
 1000 format(1x,6i8,1p26e14.6,2i4,1p14e14.4,i4,1p14e14.4)
            close(27)
          endif 
*
          if(ik2.eq.14.and.(ik1.ge.10.and.ik1.le.13)) then
            open(27,file=trim(datadir)//'/bhmergers.dat',
     &  access='append')
            write(27,1000) im1,ibx,idb,id1,id2,abs(ikind(im1)),
     &  timevo,z,(zpars(kk),kk=1,20),vss,ecc0,semi0,tb0,ik0,ik1,tphys1,
     &  tphysf1,epoch1,tmstim1,smass1,s0mass1,scmass1,semass1,rad1,
     &  radc1,rade1,slum1o,spin1,rl1,ik2,tphys2,tphysf2,epoch2,tmstim2,
     &  smass2,s0mass2,scmass2,semass2,rad2,radc2,rade2,slum2o,spin2,rl2
            close(27)
          endif                 
*
*      print binary parameters which merge and form a black hole
*      initial types of binary components: bh - all types
*
          if(ik1.eq.14.or.ik2.eq.14) then
            open(28,file=trim(datadir)//'/bhmergers-all.dat',
     &  access='append')
            write(28,1000) im1,ibx,idb,id1,id2,abs(ikind(im1)),
     &  timevo,z,(zpars(kk),kk=1,20),vss,ecc0,semi0,tb0,ik0,ik1,tphys1,
     &  tphysf1,epoch1,tmstim1,smass1,s0mass1,scmass1,semass1,rad1,
     &  radc1,rade1,slum1o,spin1,rl1,ik2,tphys2,tphysf2,epoch2,tmstim2,
     &  smass2,s0mass2,scmass2,semass2,rad2,radc2,rade2,slum2o,spin2,rl2
            close(28)
          endif 
*
*      print binary parameters which merge and form a black hole
*      initial types of binary components: all types - all types
*
          if(ik1n.eq.14.or.ik2n.eq.14) then
            open(29,file=trim(datadir)//'/bhmergers-all-all.dat',
     &  access='append')
            write(29,1010) im1,ibx,idb,id1,id2,abs(ikind(im1)),   
     &  timevo,z,(zpars(kk),kk=1,20),vss,ecc0,semi0,tb0,ik0,ik1,tphys1,
     &  tphysf1,epoch1,tmstim1,smass1,s0mass1,scmass1,semass1,rad1,
     &  radc1,rade1,slum1o,spin1,rl1,ik2,tphys2,tphysf2,epoch2,tmstim2,
     &  smass2,s0mass2,scmass2,semass2,rad2,radc2,rade2,slum2o,spin2,
     &  rl2,ik1n,ik2n
 1010 format(1x,6i8,1p26e14.6,2i4,1p14e14.4,i4,1p14e14.4,2i4)
            close(29)
          endif 
*
cAltered DCH 5/9/6 so as not to assume that the first component is the survivor
cFirst, kind may change because of call to ev_binary
*
          if(ik1n.eq.15) then
            names(im1) = id2
            ids = id2
            ipoi = 2
            sm1s = 0.d0
            write(6,4321) ik1n,ik2n,sm1s,sm2s,im1,id2
 4321       format(1x,'ipo1=2, ik1n,ik2n,sm1s,sm2s,im1,id2 = ',2i6,
     &             1p2e14.6,2i9)
          else
            names(im1) = id1
            ids = id1
            ipoi = 1
            sm2s = 0.d0
            write(6,4322) ik1n,ik2n,sm1s,sm2s,im1,id2
 4322       format(1x,'ipo1=1, ik1n,ik2n,sm1s,sm2s,im1,id2 = ',2i6,
     &             1p2e14.6,2i9)
          endif
*
          nmerge = nmerge + 1
          if(ipoi.eq.2) then
            call get_radius(id2,r1)
            call getLum(id2,xlum1)
            r2 = 0.d0
            xlum2 = 0.d0
          endif
          if(ipoi.eq.1) then
            call get_radius(id1,r1)
            call getLum(id1,xlum1)
            r2 = 0.d0
            xlum2 = 0.d0
          endif 
*
          a = 0.d0
          ecc = 0.d0
*
*      store information about binary merger
*
          iinte3(ibx) = iinte3(ibx) + 1
          inum = iinte3(ibx)
          if(inum.gt.50) then
            iinte3(ibx) = 1
            inum = 1
          endif
*
          inumt = 50*(ibx - 1) + inum
          binin(inumt,1) = r(i)/pctonbu
          binin(inumt,2) = bin(ibx,5)
          binin(inumt,3) = sm1s
          binin(inumt,4) = sm2s
          binin(inumt,5) = timevo
          binin(inumt,6) = 1.d6
          binin(inumt,7) = 3
*
          sm1o = sm1
          sm1 = sm1s/zmbar
*
          nmloev = nmloev + 1
          ssevo = sm1o - sm1
          slosev = slosev + ssevo
*
          sm2o = sm2
          sm2 = sm2s/zmbar
*
          nmloev = nmloev + 1
          slosev = slosev + (sm2o - sm2)
          ssevo = ssevo + (sm2o - sm2)
          cfr = ((sm1o - sm1) + (sm2o - sm2))*zmbar
          ssme = ssme + cfr
          body(im1) = sm1 + sm2
*      
          print*,'ddd timevo,nmerge,iime,ssme=',timevo,nmerge,iime,cfr
*
*       check for blue stragglers for binary mergers
*
          iblue = 0
          if(ik1.le.1.and.ik2.le.1.) then
            if(ik1n.le.1) then
              iblue = 1
              if(sm1s.gt.sturn.and.ibstra(im1).eq.0) then
                ibstra(im1) = 1
                ibsm = ibsm + 1
              endif
      print*,'b. str. merg1.id1,id2,ik1,ik2,ikn1,ikn2,ibsm,sm1s,sturn=',
     &               id1,id2,ik1,ik2,ik1n,ik2n,ibsm,sm1s,sturn
            endif
            if(ik2n.le.1.and.iblue.eq.0) then
              if(sm2s.gt.sturn.and.ibstra(im1).eq.0) then
                ibstra(im1) = 1
                ibsm = ibsm + 1
              endif
      print*,'b. str. merg2.id1,id2,ik1,ik2,ikn1,ikn2,ibsm,sm2s,sturn=',
     &                id1,id2,ik1,ik2,ik1n,ik2n,ibsm,sm2s,sturn
            endif
          endif
*
*       modify binary parameters and center of mass parameters
*
          call get_bs_updatetime(idb,uptime(im1))
          call get_ss_updatetime(ids,upt)
*                                  
          if(ssevo.gt.1.0d-10) then
      print*,'evo-merged: im1,idb,id1,id2,ibx,nt0,names,nameb,',
     &'ik0,ik0n,ik1,ik1n,ik2,ik2n,sm1so,sm1s,sm2so,sm2s,ss1,ss2,ssevo,',
     &'uptimeb,uptimes,nmerge = ',im1,idb,id1,id2,ibx,nt0,
     &names(im1),nameb(im1),ik0,ik0n,ik1,ik1n,ik2,ik2n,sm1so,sm1s,
     &sm2so,sm2s,ss1,ss2,ssevo*zmbar,uptime(im1),upt,nmerge
          endif
*
           nameb(im1) = 0
           bin(ibx,1) = body(im1)
           bin(ibx,2) = 0.d0
           bin(ibx,3) = 0.d0
           bin(ibx,4) = 0.d0
           bin(ibx,5) = 0.d0
           bin(ibx,6) = 0.d0
           bin(ibx,7) = 0.d0
*
*      mark obliterated binary (massless) to removal from the system
*
          if(sm1.lt.tolm.and.sm2.lt.tolm) then
            body(im1) = 0.d0
            bin(ibx,1) = 0.d0
            bin(ibx,2) = 0.d0
            iob = iob + 1
            iobt = iobt + 1
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
            r(i) = 1.0d+15 + im1
            xmax(im1) =1.01d0*r(i)          
            print*,'merge-sm  iob,im1,sm1,sm2=',iob,im1,sm1,sm2
          endif
*
          write(44,420) idb,id1,id2,im1,timevo,ss1,                     
     &                   sm1o*zmbar,ss2,sm2o*zmbar,r1,r2,
     &                   xlum1,xlum2,a,ecc,r(i),ik0,ik0n,ik1,ik1n,
     &                   ik2,ik2n,iexist,ibstra(im1),-1
*
          go to 120
*
        endif
*
*   deal with disrupted binaries due to supernovea explosion
*
        if(iexist.eq.0.and.ik0n.eq.6) then
*
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
*
          id3 = bin(ibx,7)/1.0d8
          id4 = bin(ibx,7) - 1.0d8*id3
          iname(nt0) = nt0
          ikind(im1) = 1
          ikind(nt0) = 1
          r(nt0) = r(i)*1.0000000001d0
          xmax(nt0) = 1.01d0*r(nt0)
          vr(nt0) = -vr(im1)
          vt(nt0) = -vt(im1)
          print*,'idb,id1,id2,id3,id4=',idb,id1,id2,id3,id4
*
*     check for supernove formation and NS/BH kick
*
          sm1o = sm1
          sm2o = sm2
          sm1 = sm1s/zmbar
          sm2 = sm2s/zmbar
*          
          if((ik1n.ge.13.and.ik1.lt.13).or.
     &       (ik2n.ge.13.and.ik2.lt.13)) then
            if((ik1n.eq.13.or.ik2n.eq.13).and.iflagns.le.1) then
              vs(1) = 0.d0
              vs(2) = 0.d0
              vs(3) = 0.d0
              go to 145
            endif
            if((ik1n.eq.14.or.ik2n.eq.14).and.iflagbh.le.1) then
               vs(1) = 0.d0
               vs(2) = 0.d0
               vs(3) = 0.d0
               go to 145
            endif          
            nkick(im1) = 1
            nkick(nt0) = 2
            iexistbd = 0
            ikickbd = ikickbd + 1
            ikicktbd = ikicktbd + 1
            vrold = vr(im1)
            vtold = vt(im1)
            vrold2 = vr(nt0)
            vtold2 = vt(nt0)
            vst2 = vs(2)**2 + vs(3)**2
            vst = sqrt(vst2)
            ixkickbd = 1
            ekick1 = 0.5d0*sm1*vscal**2*(vs(1)**2 + vst2)
            ekick2 = sm1*vscal*(vrold*vs(1) + vtold*vst)
            ekickbd = ekickbd + ekick1 + ekick2
            ekicktbd = ekicktbd + ekick1 + ekick2
            vr(im1) = vr(im1) + vscal*vs(1)
            vt(im1) = vt(im1)  + vscal*vst
            vst2 = vs2(2)**2 + vs2(3)**2
            vst = sqrt(vst2)
            ekick12 = 0.5d0*sm2*vscal**2*(vs2(1)**2 + vst2)
            ekick22 = sm2*vscal*(vrold2*vs2(1) + vtold2*vst) 
            ekickbd = ekickbd + ekick12 + ekick22
            ekicktbd = ekicktbd + ekick12 + ekick22
            ekicktb2 = ekicktb2 + ekick12 + ekick22            
            vr(nt0) = vr(nt0) + vscal*vs2(1)
            vt(nt0) = vt(nt0)  + vscal*vst
            write(6,8717) timevo,im1,ik1n,ik2n,ikickbd,ikicktbd,ekick1,
     &                    ekick2,ekick12,ekick22,ekickbd,ekicktbd,
     &                    vscal,vs(1),vs(2),vs(3),vr(im1),vt(im1),
     &                    vs2(1),vs2(2),vs2(3),vr(nt0),vt(nt0),nt0,
     &                    iflagns,iflagbh,ekicktb2
 8717       format(1x,'timevo-sin-bd,im1,ik1n,ik2n,ikickbd,ikicktbd,',
     &       'ekick1,ekick2,ekick12,ekick22,ekickbd,ekicktbd,vscale,',
     &       'vs(1),vs(2),vs(3),vr(im1),vt(im1),vs2(1),vs2(2),vs2(3),',
     &       'vr(nt0),vt(nt0),nt0,iflagns,iflagbh,ekicktb2 =',1pe12.4,
     &       5i8,1p17e12.4,3i8,1pe12.4)     
          endif
*
 145      continue
*
*       set the vector ibstra for binary components (now single stars)
*
          ibstra(im1) = 0
          ibstra(nt0) = 0
*
*          remove binary from the list nameb - zero means not binary
*          add star nt0 with id id3 to the names list
*
          nameb(im1) = 0
          names(im1) = id1
          names(nt0) = id2
*
          nmloev = nmloev + 1
          ssevo = sm1o - sm1
          slosev = slosev + ssevo
*
          nmloev = nmloev + 1
          slosev = slosev + (sm2o - sm2)
          ssevo = ssevo + (sm2o - sm2)
*
          call get_ss_updatetime(id1,uptime(im1))
          call get_ss_updatetime(id2,uptime(nt0))
*                                  
          if(ssevo.gt.1.0d-10) then
      print*,'evo-disrup: im1,idb,id1,id2,ibx,nt0,names1-2,nameb,',
     &'ik0,ik0n,ik1,ik1n,ik2n,sm1so,sm1s,sm2so,sm2s,ss1,ss2,ssevo,',
     &'uptimes1,uptimes2,ndist,id3,id4= ',im1,idb,id1,id2,ibx,nt0,
     &names(im1),names(nt0),nameb(im1),ik0,ik0n,ik1,ik1n,ik2,ik2n,sm1so,
     &sm1s,sm2so,sm2s,ss1,ss2,ssevo*zmbar,uptime(im1),uptime(nt0),ndist,
     &id3,id4
          endif
*                                                                                                              
          body(im1) = sm1
          body(nt0) = sm2
*
          nameb(im1) = 0
*
      print*,'im1,nt0,idb,id1,id2,names(im1),names(nt0),nameb(im1),',
     & 'sm1b,sm2b,id3,id4 =',im1,nt0,idb,id1,id2,names(im1),names(nt0),
     & nameb(im1),sm1*zmbar,sm2*zmbar,id3,id4
*
*       modify binary parameters and center of mass parameters
*
          bin(ibx,1) = body(im1)
          bin(ibx,2) = body(nt0)
          a = -1000.d0
          ecc = -1000.d0
          bin(ibx,4) = -bin(ibx,4)
          call get_radius(id1,r1)
          call getLum(id1,xlum1)
          call get_radius(id2,r2)
          call getLum(id2,xlum2)
*
*      store information about binary dissolution
*
          iinte3(ibx) = iinte3(ibx) + 1
          inum = iinte3(ibx)
          if(inum.gt.50) then
            iinte3(ibx) = 1
            inum = 1
          endif
*
          inumt = 50*(ibx - 1) + inum
          binin(inumt,1) = r(i)/pctonbu
          binin(inumt,2) = bin(ibx,5)
          binin(inumt,3) = sm1s
          binin(inumt,4) = sm2s/zmbar
          binin(inumt,5) = timevo
          binin(inumt,6) = 1.d6
          binin(inumt,7) = 6
*
          bin(ibx,3) = 1.0d-12
          bin(ibx,5) = 0.d0
*
*      mark obliterated binary (massless) to removal from the system
*
          if(sm1.lt.tolm.and.sm2.lt.tolm) then
            body(im1) = 0.d0
            body(nt0) = 0.d0
            bin(ibx,1) = 0.d0
            bin(ibx,2) = 0.d0
            iob = iob + 2
            iobt = iobt + 2
            print*,'dissolu-sm iob,im1,sm1,sm2=',iob,im1,sm1,sm2
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
            r(i) = 1.0d+15 + im1
            xmax(im1) = 1.01d0*r(i)
            r(nt0) = 1.0d+15 + nt0
            xmax(nt0) = 1.01d0*r(nt0)
          endif
*
          write(44,420) idb,id1,id2,im1,timevo,sm1*zmbar,
     &                  sm1o*zmbar,sm2*zmbar,sm2o*zmbar,r1,r2,
     &                  xlum1,xlum2,a,ecc,r(i),ik0,ik0n,ik1,ik1n,
     &                  ik2,ik2n,iexist,ibstra(im1),-1
*
          go to 120
*
      endif
*
      else if(abs(ikind(im1)).eq.1.or.abs(ikind(im1)).ge.3) then
*
        oldtime(im1) = timevo
*
        idb = nameb(im1)
        ids = names(im1)
*
        ievo = 1
*
        call setVel(ids,vs)
        call get_mass(ids,sm1so)
        sm1b = sm1so/zmbar
        sm1 = body(im1)
        sm1o = sm1
        xxm1 = abs((sm1-sm1b)/sm1)
        if(xxm1.gt.1.0d-5) then
          print*,'wrong-ss im1,ids,idb,sm1o,sm1b,ikind=',
     &            im1,ids,idb,sm1o*zmbar,sm1b*zmbar,abs(ikind(im1))
        endif
        call getLum(ids,slum1o)
*
        call get_ss_type(ids,ikio)
*
        call ev_star(ids,timevo)
*
        call get_mass(ids,sm1s)
        call get_radius(ids,r1)
        call getLum(ids,xlum1)
        call getVel(ids,vs)
*
        sm1 = sm1s/zmbar
*
        call get_ss_type(ids,ikin)
*
*     check for supernove formation and NS/BH kick
*
        if(ikin.ge.13.and.ikin.lt.15.and.ikio.lt.13) then
          if(ikin.eq.13.and.iflagns.eq.0) then
            vs(1) = 0.d0
            vs(2) = 0.d0
            vs(3) = 0.d0
            go to 150
          endif
          if(ikin.eq.14.and.iflagbh.eq.0) then
             vs(1) = 0.d0
             vs(2) = 0.d0
             vs(3) = 0.d0
             go to 150
          endif          
          nkick(im1) = 1
          iexist = 0
          ikick = ikick + 1
          ikickt = ikickt + 1
          vrold = vr(im1)
          vtold = vt(im1)
          vst2 = vs(2)**2 + vs(3)**2
          vst = sqrt(vst2)
          ixkick = 1
          ekick1 = 0.5d0*sm1*vscal**2*(vs(1)**2 + vst2)
          ekick2 = sm1*vscal*(vrold*vs(1) + vtold*vst)
          ekick = ekick + ekick1 + ekick2
          ekickt = ekickt + ekick1 + ekick2
          vr(im1) = vr(im1) + vscal*vs(1)
          vt(im1) = vt(im1)  + vscal*vst
          write(6,8719) timevo,im1,ikin,ikick,ikickt,ekick1,
     &                  ekick2,ekick,ekickt,vscal,vs(1),vs(2),
     &                  vs(3),vr(im1),vt(im1),iflagns,iflagbh
 8719     format(1x,'timevo-sing-s,im1,ikin,ikick,ikickt,ekick1,',
     &        'ekick2,ekick,ekickt,vscale,vs(1),vs(2),vs(3),vr(im1)',
     &         ',vt(im1),iflagns,iflagbh = ',1pe12.4,4i8,1p10e12.4,2i4)     
        endif
*
 150    continue
*
*       check if a blue straggler leaves the main sequence
*
        if(ibstra(im1).gt.0.and.ibstra(im1).le.4.and.
     &    ikin.gt.1) ibstra(im1) = 0
*
        nmloev = nmloev + 1
        ssevo = sm1o - sm1
        slosev = slosev + ssevo
*
        call get_ss_updatetime(ids,uptime(im1))
*
        if(ssevo.gt.1.d-10) then
      print*,'single: im1,ids,idb,ikio,ikin,sm1so,sm1s,ssevo,ut =',
     & im1,ids,idb,ikio,ikin,sm1so,sm1s,ssevo*zmbar,uptime(im1)
        endif
*        
        body(im1) = sm1
*
*      mark massless star to removal from the system
*
        if(sm1.lt.tolm) then
          body(im1) = 0.d0
          iob = iob + 1
          iobt = iobt + 1
          print*,'iob-sm,im1,sm1=',iob,im1,sm1
*
*   for obliterated objects set r(i)=1.d15+im1 to mark them for remmoval
*   from the system
*
          r(i) = 1.0d+15 + im1
          xmax(im1) = 1.01d0*r(i)
        endif
*
*     print only when star mass and luminosities change
*     because of binary evolution by more than 1%
*
        xxm1 = (sm1s - sm1so)/sm1so
        xxl1 = (xlum1 - slum1o)/slum1o
c        if(xxm1.gt.0.01d0.or.xxl1.gt.0.01d0) then
c        print*,'test',ids,im1,timevo,sm1*zmbar,sm1o*zmbar,r1,
c     &                xlum1,r(i),ikio,ikin,ibstra(im1)
        call flush(6)
        write(43,440) ids,im1,timevo,sm1*zmbar,sm1o*zmbar,r1,
     &                xlum1,r(i),ikio,ikin,ibstra(im1),ikind(im1)
*
 440    format(1x,'@',2i7,1p6e16.8,4i4)
*
        go to 120
*
*      Calculate masslose and heating due to stelar evolution
*      For neutron stars assigne proper velocity - to be done
*      Update the potential
*
      else
*
        write (6,*) 'ikind != 1, 2 or 3 in mloss_single'
        stop
*
      endif
*
 120  continue
*
      if(imodel.ge.3) then
        cor = -0.5d0*ssevo**2/r(i)
*         ehpot = -ssevo*(-smt/rtid - u(i))
        ehpot = -ssevo*u(i)
        if(ixkick.eq.1.or.ixkickbs.eq.1.or.ixkickbd.eq.1) then
          ehkin = 0.5d0*ssevo*(vrold**2 + vtold**2)
        else  
          ehkin = 0.5d0*ssevo*(vr(im1)**2 + vt(im1)**2)
        endif
        ehmlev = ehmlev + ehpot - ehkin + cor
      else
        cor = -0.5d0*ssevo**2/r(i)
        ehpot = -ssevo*u(i)
        if(ixkick.eq.1.or.ixkickbs.eq.1.or.ixkickbd.eq.1) then
          ehkin = 0.5d0*ssevo*(vrold**2 + vtold**2)
        else
          ehkin = 0.5d0*ssevo*(vr(im1)**2 + vt(im1)**2)
        endif
        ehmlev = ehmlev + ehpot - ehkin + cor
      endif
*
*      update total mass, tidal radius and potential
*
      smt = smt - ssevo
c      if(imodel.ge.3) rtid = rtidkg*smt**0.3333333
c      if(ndist.eq.0) call coepot
*
*     Only changes in the kinetic energy are used in the energy balance
*     The changes of the potrntial energy are taken ibto account by 
*     calling before and after interaction the coepot subroutine
*
      ehmlev1 = ehmlev
      ehmlev = ehmlev - ehpot - cor
*
      if(iexistbs.eq.0) iexist = 0
      if(iexistbd.eq.0) iexist = 0          
*
      print*,'im1,i,ixkick,ixkickbs,ixkickbd,iime,ssme,ssevo,cor,',
     & 'ehpot,ehkin,ehmlev1,ehmlev =',im1,i,ixkick,ixkickbs,ixkickbd,
     & iime,ssme,ssevo*zmbar,cor,ehpot,ehkin,ehmlev1,ehmlev
*
      if(ndist.eq.0.and.iob.eq.0) return
*
*     deal with change of number of single stars due to binary disruption
*     compute new potential and sort stars after destroing binaries in
*     binary--binary interactions. One component of destroied binary
*     is substituded for centre of binary mass, for another component
*     name of one of stars which previeously created binary is taken
*
      ndiste = ndiste + ndist
      n = n + ndist - iob
      nt = nt + ndist - iob
c      call sort2(nt0,r,iname)
      nzst(nsuzon) = nt
*
      return
*
      end
*
*
*
