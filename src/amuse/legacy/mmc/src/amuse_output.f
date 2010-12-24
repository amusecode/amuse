* this code does the same as the original output subroutine, except
* it does not stop, instead returns an error (or stopping) code

      function amuse_output()
*
*
*       generate information about the system for each model
*       ----------------------------------------------------------
*
*     
      include 'common.h'
*
      real*8 velo(nlagra),de,etotn,ro5,eintb3,xktc,xxn,smvc,tphys,
     &       avms10,avmsh,avmst,avwd10,avwdh,avwdt,avns10,avnsh,avnst,
     &       ams10,amsh,amst,awd10,awdh,awdt,ans10,ansh,anst,rms10,
     &       rmsh,rmst,rwd10,rwdh,rwdt,rns10,rnsh,rnst,smms,smwd,
     &       smnsbh,avgs10,avgsh,avgst,ags10,agsh,agst,rgs10,rgsh,rgst,
     &       bm10,bmh,bmt,ba10,bah,bat,br10,brh,brt,xbar10,xbarh,xbar,
     &       xbet10,xbeth,xbet,xrs10,xrsh,xrst,amas10,amash,amast,smb,
     &       sms,ehbint,vscal,slum,slum1,sm1,sm2,ap,teff,epoch1,epoch2,
     &       temp,smsl,smslh,smsl10,smsl10h,rmsl,rmslh,rmsl10h,txxx,
     &       ekick,zkink,spin1,spin2,ecc,rad1,rad2,slum2,smgs,pctonbu,
     &       csb,rtid0,rtid1,xrtid,smt00,ekickbs,ekickbd,rcob,ttime
*
      real*4 runtim,cpp,xxx,timerun
*
      integer i,ii,ibirel,ibiesc,nms,nwd,nnsbh,n,ns,nb,inum,inumt,
     &        ids,ik1,ik2,ikb,idb,idstar,k,ngs,kbsm,kbsc,kbs3,kbs4,l,
     &        im1,idestr,imerge,lmin,lmax,nup,nmsl,nmsl10,nmsll,idd,
     &        iii,nbh,nbb,irtid,int,int0,lwd,lns,lbh,lhg,lgb,lch,lfag,
     &        lsag,lhms,lhhg,lhgb,l2wd,l2ns,l2bh,lwdms,lwdns,lwdbh,
     &        lwdot,lnsms,lnsbh,lnsot,lbhms,lbhot,lbh3,lbh2,idb1,idb2,
     &        kebs,kebb,kdbs,kdbb,n1000
*
      common /colla/ avms10,avmsh,avmst,avwd10,avwdh,avwdt,avns10,
     &               avnsh,avnst,ams10,amsh,amst,awd10,awdh,awdt,
     &               ans10,ansh,anst,rms10,rmsh,rmst,rwd10,rwdh,rwdt,
     &               rns10,rnsh,rnst,smms,smwd,smnsbh,avgs10,avgsh,
     &               avgst,ags10,agsh,agst,rgs10,rgsh,rgst,bm10,bmh,
     &               bmt,ba10,bah,bat,br10,brh,brt,smb,xbar10,xbarh,
     &               xbar,xbet10,xbeth,xbet,xrs10,xrsh,xrst,sms,
     &               smgs,amas10,amash,smsl,smslh,smsl10,smsl10h,
     &               rmsl,rmslh,rmsl10h,nms,ngs,nwd,nnsbh,nb,ns,
     &               nmsl,nmsl10,nmsll
cAdded DCH 4/1/7      
      double precision mb,mv,mi,mb1,mv1,mb2,mv2,mi1,mi2

      integer amuse_output

*
*
*       open output files
*
*                 f77
*
*
      cpp = cpu
      write(6,*)   ' in  output'
      call flush(6)
      pctonbu = rtidkg/rbar
*
      open(10,file='lagrangi.dat',access='append')
      open(11,file='velocity.dat',access='append')
      open(12,file='anisotro.dat',access='append')
      open(13,file='system.dat',access='append')
      open(14,file='core.dat',access='append')
      open(20,file='collabo.dat',access='append')
*
*                 f90
*
*      open(10,file='lagrangi.dat',position='append') 
*      open(11,file='velocity.dat',position='append')
*      open(12,file='anisotro.dat',position='append')
*      open(13,file='system.dat',position='append')
*      open(14,file='core.dat',position='append')
*      open(15,file='escape.dat',position='append')
*      open(20,file='collabo.dat',position='append')
*
*
*     determine the physical time
*
      ekick = 0.d0
      ekickbs = 0.d0
      ekickbd = 0.d0
      tphys = timet
      if(tcrit.lt.1000.d0) then
        ttime = time
      else
        ttime = tphys
      endif
      txxx = tscale0*log(gamma*nt00)/log(gamma*nt)*
     &       float(nt)/float(nt00)
      n = nt
*
*     change tau and bmin/bmax for tphys > tcrevo
*
      if(tphys.le.tcrevo) then
        bmin = bmin0
        bmax = 2.d0*bmin0
        tau = tau0
      else
        bmin = ybmin*bmin0
        bmax = 2.d0*bmin
        tau = ytau*tau0
      endif  
*
*     call mloss for all stars between lmin=1 and lmax=nt and force
*     stellar evolution up to the time tphys
*     nup = -100 for this case
*
      if((tphys.gt.tte).or.(ttime.gt.tcrit)) then
        nup = -100
        if(tphys.le.tcrevo) then
          tte = tte + dtte0
        else
          tte = tte + dtte
        endif
*cello  print*,'tphys,timet,time,tcrit,tcrevo,tte,dtte0,dtte ',
*     &          tphys,timet,time,tcrit,tcrevo,tte,dtte0,dtte
      else
        nup = 1
      endif
      lmin = 1
      lmax = nt
*      
      irtid = 0
      int0 = nt
      n = nt
      rtid0 = rtid
      if(tphys.gt.0.d0) then
        call mloss(nup,lmin,lmax,tphys,ekick,ekickbs,ekickbd)
        rtid1 = rtid
        int = nt
        etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
     &         ehmlev - ehcoll + eccoll - ekickt - enekin -
     &         ekicktbs - ekicktbd - ekicktwd    
        write(6,6543) tphys,etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &                ekickt,ekicktbs,ekicktbd,ekicktwd,enekin,escsta,
     &                ehbin3,ehb3b3        
 6543   format(1x,'output-m tphys,etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     & 'eccol,ekickt,ekicktbs,ekicktbd,ekicktwd,enekin,escsta,ehbin3,',
     & 'ehb3b3 =',1p16e20.12)
*
*   iteration of the tidal redius after stellar evolution
*
        if(nitesc.eq.0) go to 8717
 8712   continue       
        iesc = 0
        smt00 = smt
        call escape(-1,tphys)
        nescrt = nescrt + iesc
        sescrt = sescrt + (smt00 - smt)
        xrtid = (rtid0 - rtid)/rtid0
        irtid = irtid + 1
        write(6,9154) irtid,int0,int,nt,rtid0,rtid1,rtid,xrtid,iesc,
     &                nescrt,sescrt
 9154   format(1x,'irt,int0,int,nt,rt0,rt1,rt,xrt,iesc,nescrt,sescrt=',
     &         i4,3i9,1p4e14.6,2i9,1pe14.6)
        rtid0 = rtid
        call sort2(nt0,r,iname)
        n = nt
        int = nt
        if(xrtid.gt.0.003d0) go to 8712
      endif
*
      if(irtid.gt.0) then
        nzst(nsuzon) = nt
        call coepot
      endif
*
 8717 continue
*       
      call energy(2)
c      print*,'after escape in mloss  int,nt=',int,nt
      n = nt
*
      etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
     &        ehmlev - ehcoll + eccoll - ekickt - enekin -
     &        ekicktbs - ekicktbd - ekicktwd     
      write(6,6544) tphys,etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &              ekickt,ekicktbs,ekicktbd,ekicktwd,enekin,escsta,
     &              ehbin3,ehb3b3                 
 6544 format(1x,'output-e tphys,etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     & 'eccol,ekickt,ekicktbs,ekicktbd,ekicktwd,enekin,escsta,ehbin3,',
     & 'ehb3b3 =',1p16e20.12)                
*celloprint*,'tphys,ltot: ',tphys,ltot
*
*       calculate conservation of the total energy
*
c      call energy(2)
*
      etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 - 
     &        ehmlev - ehcoll + eccoll - ekickt - enekin -
     &        ekicktbs - ekicktbd - ekicktwd
c      de = (etotn - etot)/zkin
      de = etotn - etot
      error = error + de
      etot = etotn
*
cChanged for M67
c
      if(imodel.gt.3.and.rtid.le.0.d0) then
        write(6,*) ' rtid = 0.0',rtidkg,smt
        call flush(6)
        call mydump(1)
      endif
*
*       determine lagrangian radii and average radial and tangential 
*       velocities between lagrangian radii
*
c
      write(6,*)  ' call lagr '
      call flush(6)
c
      call lagrad
c
      write(6,*)  ' out lagr'
      call flush(6)
*
      write(10,100) iseed,time,tphys,(rlag(i),i=1,nlagra)
*
      do 20 i=1,nlagra
   20    velo(i) = sqrt(v2rl(i) + v2tl(i))
*
      write(11,100) iseed,time,tphys,(velo(i),i=1,nlagra)
*
      write(12,100) iseed,time,tphys,(ani(i),i=1,nlagra)
*
  100 format(1x,i5,1p102e12.4,:)
*
      close(10)
      close(11)
      close(12)
*
*     total binary binding energy - energy 'generated' + 'primordial'
*
      ehbint = ehbin3 + ehbi3p
*
*     compute curent number of blue stragglers
*
      kbsm = 0
      kbsc = 0
      kbs3 = 0
      kbs4 = 0
      kebs = 0
      kebb = 0
      kdbs = 0
      kdbb = 0
      csb = 0.0d0
      nbb = 0
      lwd = 0
      lns = 0
      lbh = 0
      lbh2 = 0
      lbh3 = 0
      lwdms = 0
      lwdns = 0
      lwdbh = 0
      lwdot = 0
      lnsms = 0
      lnsbh = 0
      lnsot = 0
      lbhms = 0
      lbhot = 0
      l2wd = 0
      l2ns = 0
      l2bh = 0
      lhg = 0
      lgb = 0
      lch = 0
      lfag = 0
      lsag = 0
      lhms = 0
      lhhg = 0
      lhgb = 0
      do 60 l = 1,n
         im1 = iname(l)
         if(ikind(im1).eq.2) then 
           idb = nameb(im1)
           call get_loid(idb,idb1)
           call get_hiid(idb,idb2)
           if(ikind(im1).eq.2) nbb = nbb + 1
           if(ibstra(idb1).eq.1) kbsm = kbsm + 1
           if(ibstra(idb1).eq.2) kbsc = kbsc + 1
           if(ibstra(idb1).eq.3) kbs3 = kbs3 + 1
           if(ibstra(idb1).eq.4) kbs4 = kbs4 + 1
           if(ibstra(idb1).eq.5) kebs = kebs + 1
           if(ibstra(idb1).eq.6) kebb = kebb + 1
           if(ibstra(idb1).eq.7) kdbs = kdbs + 1
           if(ibstra(idb1).eq.8) kdbb = kdbb + 1
           if(ibstra(idb2).eq.1) kbsm = kbsm + 1
           if(ibstra(idb2).eq.2) kbsc = kbsc + 1
           if(ibstra(idb2).eq.3) kbs3 = kbs3 + 1
           if(ibstra(idb2).eq.4) kbs4 = kbs4 + 1
           if(ibstra(idb2).eq.5) kebs = kebs + 1
           if(ibstra(idb2).eq.6) kebb = kebb + 1
           if(ibstra(idb2).eq.7) kdbs = kdbs + 1
           if(ibstra(idb2).eq.8) kdbb = kdbb + 1
         else
           ids = names(im1)
           if(ibstra(ids).eq.1) kbsm = kbsm + 1
           if(ibstra(ids).eq.2) kbsc = kbsc + 1
           if(ibstra(ids).eq.3) kbs3 = kbs3 + 1
           if(ibstra(ids).eq.4) kbs4 = kbs4 + 1
           if(ibstra(ids).eq.5) kebs = kebs + 1
           if(ibstra(ids).eq.6) kebb = kebb + 1
           if(ibstra(ids).eq.7) kdbs = kdbs + 1
           if(ibstra(ids).eq.8) kdbb = kdbb + 1
         endif
         slum = 0.d0
         if(ikind(im1).eq.1.or.ikind(im1).ge.3) then
           ids = names(im1)
           call getlum(ids,slum)
           call get_ss_type(ids,ik1)
           if(ik1.eq.2) lhg = lhg + 1
           if(ik1.eq.3) lgb = lgb + 1
           if(ik1.eq.4) lch = lch + 1
           if(ik1.eq.5) lfag = lfag + 1
           if(ik1.eq.6) lsag = lsag + 1
           if(ik1.eq.7) lhms = lhms + 1
           if(ik1.eq.8) lhhg = lhhg + 1
           if(ik1.eq.9) lhgb = lhgb + 1
           if(ik1.gt.9.and.ik1.lt.13) lwd = lwd + 1
           if(ik1.eq.13) lns = lns + 1
           if(ik1.eq.14) lbh = lbh + 1
           if(ik1.eq.14.and.ikind(im1).eq.3) lbh3 = lbh3 + 1
           if(ik1.eq.14.and.ikind(im1).eq.4) lbh2 = lbh2 + 1
         elseif (ikind(im1).eq.2) then
           idb = nameb(im1)
           call get_loid(idb,ids)
           call get_ss_type(ids,ik1)
           call getlum(ids,slum1)
           call get_hiid(idb,ids)
           call get_ss_type(ids,ik2)
           call getlum(ids,slum2)
           slum = slum1 + slum2
           if(ik1.eq.14.and.ik2.eq.14) l2bh = l2bh + 1
           if(ik1.eq.13.and.ik2.eq.13) l2ns = l2ns + 1
           if((ik1.gt.9.and.ik1.lt.13).and.(ik2.gt.9.and.ik2.lt.13)) 
     &        l2wd = l2wd + 1
           if((ik1.eq.14.and.ik2.le.1).or.(ik2.eq.14.and.ik1.le.1)) 
     &        lbhms = lbhms + 1
           if((ik1.eq.14.and.(ik2.ge.2.and.ik2.le.9)).or.
     &        (ik2.eq.14.and.(ik1.ge.2.and.ik1.le.9))) 
     &        lbhot = lbhot + 1
           if((ik1.eq.13.and.ik2.le.1).or.(ik2.eq.13.and.ik1.le.1)) 
     &        lnsms = lnsms + 1          
           if((ik1.eq.13.and.ik2.eq.14).or.(ik2.eq.13.and.ik1.eq.14)) 
     &        lnsbh = lnsbh + 1
           if((ik1.eq.13.and.(ik2.ge.2.and.ik2.le.9)).or.
     &        (ik2.eq.13.and.(ik1.ge.2.and.ik1.le.9)))
     &        lnsot = lnsot + 1
           if(((ik1.ge.10.and.ik1.le.12).and.ik2.le.1).or.
     &        ((ik2.ge.10.and.ik2.le.12).and.ik1.le.1)) 
     &        lwdms = lwdms + 1
           if(((ik1.ge.10.and.ik1.le.12).and.ik2.eq.13).or.
     &        ((ik2.ge.10.and.ik2.le.12).and.ik1.eq.13))
     &        lwdns = lwdns + 1
           if(((ik1.ge.10.and.ik1.le.12).and.ik2.eq.14).or.
     &        ((ik2.ge.10.and.ik2.le.12).and.ik1.eq.14))          
     &        lwdbh = lwdbh + 1
           if(((ik1.ge.10.and.ik1.le.12).and.(ik2.ge.2.and.ik2.le.9))
     &     .or.((ik2.ge.10.and.ik2.le.12).and.(ik1.ge.2.and.ik1.le.9)))          
     &        lwdot = lwdot + 1    
         endif
           if(slum.gt.0.0.and.r(l).gt.0.d0)
     &        csb = csb + slum/(r(l)*rbar/rtidkg)**2
 60   continue
      csb = csb/(2.d0*pi)
      slum = 0.d0
*                                
*celloprint*,' ekick,ekickbs,ekickbd = ',ekick,ekickbs,ekickbd
      if(ekick.gt.0.d0.or.ekickbs.gt.0.d0.or.ekickbd.gt.0.d0) then
        zkink = zkin - ekick - ekickbs - ekickbd
      else
        zkink = zkin
      endif
*
      xxx = runtim(cpp)/60.0
      if(timeold.lt.xxx) then
        timeold = xxx
      else
        itime = itime + 1
        timeold = xxx
      endif
      timerun = itime*1440.0 + xxx                                      
*
*       determine and print central parameters
*
      write(6,*)  ' in core'
      call flush(6)          
c
      call core
c
      write(6,*)  ' out core'
      call flush(6)          
*
*     compute the "observational" core radius (rcob = distance where the
*     surface brightnes is half of the central surface brightnes)
*
      call profiles(tphys,vscal,0,rcob)
* 
      write(13,110) iseed,time,tphys,smt,etot,zkink,pot,escsta,ehbint,
     &              ehb3b3,ehmlev,ehcoll,eccoll,error,enepot,sloses,
     &              slosev,slosco,rc,rcob,rlag(11),rtid,nt,nescst,
     &              nescm,nmloev,ncoll,nbin3,nescb3,ndist3,ndist4,
     &              ndiste,nmerg3,nmerg4,nmerge,nbin3-nescb3-ndist3-
     &              ndist4-ndiste-nmerg3-nmerg4-nmerge,kbsm,kbsc,kbs3,
     &              kbs4,kebs,kebb,kdbs,kdbb,kbsm+kbsc+kbs3+kbs4+kebs+
     &              kebb+kdbs+kdbb,ibsm+ibsc+ibs3+ibs4,ibsm,ibsc,ibs3,
     &              ibs4,nt0,ivnewg,ivrr,enrad,timerun,ehbin3,txxx,
     &              ekickt,ekicktbs,ekicktbd,ekicktb2,ikickt,ikicktbs,
     &              ikicktbd,mnsbh,nexchang,nexchang2,nbb,de,csb,nescrt,
     &              sescrt,lwd,lns,lbh,lhg,lgb,lch,lfag,lsag,lhms,lhhg,
     &              lhgb,l2wd,l2ns,l2bh,lwdms,lwdns,lwdbh,lwdot,lnsms,
     &              lnsbh,lnsot,lbhms,lbhot,lbh2,lbh3,ntsn1,ntsn2,ntsnb,
     &              ikicktwd,ekicktwd,ntwd1
*
  110 format(1x,i5,1p7e12.4,1pe15.7,1p13e12.4,29i10,2i12,1p8e12.4,7i7,
     &       1p2e12.4,i9,1pe12.4,29i9,1pe12.4,i9)
*
      close(13)
*
*     prepare output for collaboration project
*
      write(6,*)   '  in collabr '
*
cAltered DCH 3/8/6:
      call colabr(nbh)
*
      write(6,*)   '  out collabr  '
*      
      amast = smt/float(nt)
      write(20,190) iseed,time,tphys,smt,nt,rlag(7),rlag(11),rtid,
     &              amas10,amash,amast,ani(7),ani(11),ani(16),sms,ns,
     &              xrs10,xrsh,xrst,xbar10,xbarh,xbar,xbet10,xbeth,
     &              xbet,smms,nms,rms10,rmsh,rmst,avms10,avmsh,avmst,
     &              ams10,amsh,amst,smgs,ngs,rgs10,rgsh,rgst,avgs10,
     &              avgsh,avgst,ags10,agsh,agst,smwd,nwd,rwd10,rwdh,
     &              rwdt,avwd10,avwdh,avwdt,awd10,awdh,awdt,smnsbh,
     &              nnsbh,rns10,rnsh,rnst,avns10,avnsh,avnst,ans10,
     &              ansh,anst,smb,nb,br10,brh,brt,bm10,bmh,bmt,ba10,
     &              bah,bat,nbh,smsl,smslh,rmsl,rmslh,nmsl,smsl10,
     &              smsl10h,rmsl10h,nmsl10,nmsll
*
 190  format(1x,i5,1p3e12.4,i7,1p10e12.4,i7,1p10e12.4,i7,1p10e12.4,
     &       i7,1p10e12.4,i7,1p10e12.4,i7,1p10e12.4,i7,1p9e12.4,i5,
     &       1p4e12.4,i7,1p3e12.4,2i7)
*
      close(20)
*
*      print core parameters
*
      xktc = smc*vc*vc/(3.0d0*nc)
      smvc = 0.d0
      do 111 i = 1,5
      ii = iname(i)
 111  smvc = smvc + body(ii)
*
      ro5 = 1.5d0*smvc/twopi/r(5)**3
***************
*Added DCH 8/2/10
      nddb = 0
      nbc = 0
      sumr2rho2 = 0.d0
      sumrho2 = 0.d0
      do i = 1,n
         call densi(i,den)
         sumr2rho2 = sumr2rho2 + r(i)**2*den**2
         sumrho2 = sumrho2 + den**2
      enddo
      rdens = sqrt(sumr2rho2/sumrho2)
      do i = 1,n
         k = iname(i)
         if (ikind(k).eq.2.and.r(i).lt.rdens) then
            nbc = nbc + 1
            idb = nameb(k)
            call get_loid(idb,idstar)           
            call get_ss_type(idstar,ik1)
            call get_hiid(idb,idstar)
            call get_ss_type(idstar,ik2)
            if (ik1.ge.10.and.ik2.ge.10) nddb = 
     &           nddb + 1
         endif
      enddo
      rdens = rdens*rbar/rtidkg
*******************************************
*
      write(14,120) iseed,time,tphys,smc,rc,vc,roc,u(1),
c     &              -3.d0*u(1)/vc/vc,r(5),ro5,nc,rcob,r(11),rtid
     &              -3.d0*u(1)/vc/vc,r(5),ro5,nc,rcob,r(11),rtid,rdens,
     &     nddb,nbc,n
*
c  120 format(1x,i5,1p10e12.4,i9,1p3e12.4)
  120 format(1x,i5,1p10e12.4,i9,1p4e12.4,3i8)
*
      close(14)
*   
*     print escapers if there are any
*
      if(iesc1.gt.0) then
*      
        open(15,file='escape.dat',access='append')       
*
        do 30 i=1,iesc1
   30      write(15,130) iseed,time,tphys,ienam1(i),iekind(i),
     &                   escmas(i)*zmbar,escdis(i)/pctonbu,
     &                   escene(i),escang(i),esctim(i),xktc,
     &                   rc/pctonbu
  130      format(1x,i5,1p2e12.4,2i8,1p7e12.4)
        iesc1 = 0
        close(15)
      endif
*
*     print information about  binaries if any
*
      if(nbin3.gt.0) then
*
        open(17,file='bin3glo.dat',access='append')
c        open(18,file='bin3inf.dat',access='append')
c        open(19,file='bin3int.dat')
c        open(19,file='bin3int.dat',access='append')
*        open(17,file='bin3glo.dat',position='append')
*        open(18,file='bin3inf.dat',position='append')
*        open(19,file='bin3int.dat',position='append')
*
        eintb3 = 0.0d0
        ibiesc = 0
        ibirel = 0
        idestr = 0
        imerge = 0
*
        do 150 i=1,nbin3
*
           if(bin(i,7).lt.0.0d0) ibiesc = ibiesc + 1
           if(bin(i,6).lt.0.0d0) ibirel = ibirel + 1
           if(bin(i,4).lt.0.0d0) idestr = idestr + 1
           if(bin(i,3).eq.0.0d0) imerge = imerge + 1
*
*        find the name of binary i
*
           im1 = nbinar(i)
*
           if(bin(i,7).ge.0.d0.and.bin(i,6).ge.0.d0.and.
     &       bin(i,4).ge.0.d0.and.bin(i,3).gt.0.d0.and.
     &       bin(i,1).ge.0.d0.and.bin(i,2).ge.0.d0) then
             eintb3 = eintb3 + bin(i,5)
*
c             if (iprint.eq.0) then
c                write(18,170) iseed,time,tphys,im1,i,
c     &               (bin(i,iii),iii=1,6),bin(i,7),bin(i,8),xktc,rc
c             endif
           endif
*
c           inum = iinte3(i)
c           if(inum.eq.0) go to 150
c           do 160 ii = 1,inum
c              inumt = 50*(i-1) + ii
c              write(19,180) iseed,time,tphys,im1,i,ii,
c     &                      (binin(inumt,iii),iii=1,7)
c 160       continue
 150    continue
*
        write(17,140) iseed,time,tphys,nbin3,nbin3-nescb3-ndist3-
     &                ndist4-ndiste-nmerg3-nmerg4-nmerge,nb3b3,ndist3,
     &                ndist4,ndiste,nmerg3,nmerg4,nmerge,nb3fin,nesb3s, 
     &                nescb3,ibiesc,ibirel,idestr,imerge,ehbin3,eintb3,
     &                escb3s,escbi3,escbb3,erelb3,erb3in,ehb3b3,slob3s,
     &                slob3b,nexchang,nexchang2
  140   format(1x,i5,1p2e12.4,16i7,1p10e12.4,2i7)
        close(17)
*
c 170    format(1x,i5,1p2e12.4,2i8,1p6e12.4,1pe24.16,1pe12.4,1p2e12.4)
c 180    format(1x,i5,1p2e12.4,3i8,1p7e12.4)
c        close(18)
c        close(19)
*
      endif
c
*
*     snapshot, profile
*
      write(6,*)  ' snapshot '
      call flush(6)
*
c      if((tphys.eq.0.d0).or.(tphys.ge.ttp).or.(ttime.gt.tcrit)) then
c        n1000 = n
*       save restart data
c        call mydump(1)
c      else
c        n1000 = 1000 
c      endif                
      n1000 = 1000
      if(tphys.gt.tcrit) n1000 = n
*
      open(21,file='snapshot.dat',access='append')
*      open(21,file='snapshot.dat',position='append')
*
      ttp = ttp + dttp
*celloprint*,'tphys,ttp,dttp = ',tphys,ttp,dttp
      xxn = zmbar*rtidkg/rbar
      vscal = 0.06558d0*sqrt(xxn)

cIn the following (DCH) the temperature is in units of (4*pi*sigma)**-0.25
      write(21,210) n,tphys
 210  format(1x,i9,1pe12.4,'  *****')
*
      ltot = 0.d0
      do 200 i = 1,n1000
         k = iname(i)
         if(r(i).gt.rtid) then
*cello     print*,'r>rt tim,r,rt=',tphys,r(i),rtid
           go to 200
         endif
         if(abs(body(k)).lt.tolm) go to 200
         if(ikind(k).lt.0) then
*cello     print*,'output: i,k,ikind = ',i,k,ikind(k)
         endif
         if(ikind(k).eq.1.or.ikind(k).ge.3) then
           ids = names(k)
           idd = ids
cAdded DCH 9/9/6
           call getLum(ids,slum1)
           ltot = ltot + slum1
           call getSpin(ids,spin1)
           spin2 = 0.d0
           slum = slum1
           if(isNaN(slum)) then
*cello       print*,'i,im1,ids,r,lum = ',i,k,ids,r(i),slum1
             go to 200
           endif
           call get_ss_type(ids,ik1)
           call getRadius(ids,rad1)
           call getEpoch(ids,epoch1)
           if (rad1.gt.0.d0) then
              teff = slum1**0.25d0/sqrt(rad1)
           else
              teff = 0.d0
           endif
           ik2 = 0
           ikb = 0
           sm2 = 0.d0
           ap = 0.d0
           ecc = 0.d0
           epoch2 = 0.d0
           rad2 = 0.d0
           slum2 = 0.d0
           call get_mass(ids,temp)
           sm1 = temp
           call mycolour(ik1,sm1,slum,rad1,mv,mb,mi)
           mv1 = mv
           mv2 = 0.d0
           mb1 = mb
           mb2 = 0.d0
           mi1 = mi
           mi2 = 0.d0
         else
           idb = nameb(k)
           idd = idb
cAdded DCH 9/9/6
           call get_bs_type(idb,ikb)
           call get_sma(idb,ap)           
           call get_loid(idb,idstar)           
           call getEpoch(idstar,epoch1)
           call getSpin(idstar,spin1)
           call getEcc(idb,ecc)
           call get_loid_mass(idb,sm1)
           temp = sm1
           call get_ss_type(idstar,ik1)
           call getLum(idstar,slum1)
           ltot = ltot + slum1
           call getRadius(idstar,rad1)
           call get_hiid(idb,idstar)
           call getEpoch(idstar,epoch2)
           call getSpin(idstar,spin2)
           call get_hiid_mass(idb,sm2)
           call get_ss_type(idstar,ik2)
           call getLum(idstar,slum2)
           ltot = ltot + slum2
           call getRadius(idstar,rad2)
           if(rad1.gt.0.d0.and.rad2.gt.0.d0) then
             if(ik1.le.9.and.ik2.le.9) then
               teff = (slum1**1.25d0/sqrt(rad1) + 
     &         slum2**1.25d0/sqrt(rad2))/(slum1 + slum2)
               go to 159
             endif
             if(ik1.gt.9.and.ik2.le.9) then
               teff = slum2**0.25d0/sqrt(rad2)
               go to 159
             endif
             if(ik2.gt.9.and.ik1.le.9) then
               teff = slum1**0.25d0/sqrt(rad1)
               go to 159
             endif
             if(ik1.gt.9.and.ik2.gt.9) then
               teff = (slum1**1.25d0/sqrt(rad1) +
     &         slum2**1.25d0/sqrt(rad2))/(slum1 + slum2)
               go to 159
             endif
           else
             teff = 0.d0
           endif
 159       slum = slum1 + slum2
           call mycolour(ik1,sm1,slum1,rad1,mv1,mb1,mi1)
           call mycolour(ik2,sm2,slum2,rad2,mv2,mb2,mi2)           
           mv = -2.5d0*log10(10.d0**(-0.4d0*mv1)+10.d0**(-0.4d0*mv2))
           mb = -2.5d0*log10(10.d0**(-0.4d0*mb1)+10.d0**(-0.4d0*mb2))
           mi = -2.5d0*log10(10.d0**(-0.4d0*mi1)+10.d0**(-0.4d0*mi2))
         endif   
*       
         if(ikind(k).eq.1.or.ikind(k).ge.3) ikb = 0
*
         write(21,220) i,k,idd,body(k)*zmbar,temp,sm1,sm2,slum,slum1,
     &                 rad1,rad2,r(i)*rbar/rtidkg,
     &                 vr(k)*vscal,vt(k)*vscal,ap,ik1,ik2,ikb,ikind(k),
     &                 nwhich(k),teff,epoch1,epoch2,ibstra(k),mv,mb-mv,
     &                 mi,mv1,mb1-mv1,mi1,mv2,mb2-mv2,mi2,spin1,spin2,
     &                 ecc,inexch(k),inexch2(k)
 220     format(1x,3i8,1p12e12.4,4i5,i8,1p3e12.4,i5,1p12e12.4,2i6)
         if(ikind(k).ge.3.and.abs(body(k)*zmbar - temp).gt.0.01) then
*cello     print*,'output k,body(k)*zmbar,temp,time: ',
*cello&            k,body(k)*zmbar,temp,time
         endif
 200  continue
*
      close(21)
*
      write(6,*)  ' in profiles '
      call flush(6)
*            
      call profiles(tphys,vscal,1,rcob)
*
      write(6,*)  ' out profiles '
      call flush(6)
*
      write(6,*)  ' out output',timerun
      call flush(6)
*            
      if(abs(de).gt.qe) then
        write(6,*)'energy conservation is not satisfied! ',time,etot,de
*       stop 'energy error'
        amuse_output = -1
        return
      endif
*
*       check if cpu is greater than tcomp
*
c      if(tcomp.lt.timerun) then
c        write(6,*) 'cpu time is greater than tcomp   ',time
c        stop 'cpu > tcomp'
c      endif
*
*       check if time is greater than tcrit
*
c      if(ttime.gt.tcrit) then
      if(tphys.gt.tcrit) then
        write(6,*) 'time/tphys is greater than tcrit   ',time,tphys
*        stop 'time/tphys > tcrit'
        amuse_output = -2
        return
      endif
*
*     stop when the number of objects is too small
*
      if(nt*gamma.le.3.d0) then
        write(6,*) '     N*gamma < 3.0    ',nt,time
*        stop 'N*gamam < 3.0 '
        amuse_output = -3
        return
        
      endif
*
      amuse_output = 0
      return
*
      end function
