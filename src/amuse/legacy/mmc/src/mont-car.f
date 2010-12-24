*
*
*
*         Program Monte-Carlo  MONTCARL
*         -----------------------------
*
*
*
      include 'common.h'
*
      integer iphase
*
*
*
*
*
*             read initial parameters
*
*
      print*,'calling input'
      call flush(6)
      call input
*     
      if(istart.eq.1) then
*
*
*             compute initial model and save it
*
*
         call start
*
         iphase = 1
*
      else
*
*
*             read restart model
*
*
         call mydump(2)
*            
*            read parameters which will be changed
*
         call input
*
         go to 40
*
      endif
*
*
 10   continue
*
      go to (20, 30) iphase
*
*
*             devision on zones and super-zones and output
*
*
 20   call output
 40   call zone
      iphase = 2
      go to 10
*
*
*            main flow rutine: calculation of cluster evolution
*
*
 30   call relaxt
      iphase = 1
      go to 10
*
      end
*
*
*
*
 
      subroutine coepot
*
*
*         compute smooth potential for all particles
*         ------------------------------------------
*
*
      include 'common.h'
*
      real*8 ux,uy,smx,rx
*
      integer i,n,im
*
*
      n = nt
*
      smx = smt
      uy = 0.0d0
      rx = r(n)
      ux = -smx/rx
      u(n) = ux
*
      do 10 i=n-1,1,-1
*    
         im = iname(i+1)
         smx = smx - body(im)
         ux = -smx/r(i)
         uy = uy - body(im)/rx
         rx = r(i)
         u(i) = ux + uy
*
 10   continue
*
*     
      return
*
      end
*
*
*
*    
      subroutine colabr(nbh)
* 
*
*       lagrangian radii & mean mass & radial and tangential velocities
*       --------------------------------------------------------------- 
*       & anisotropy in the lagrangian shells 
*       -------------------------------------
*
*
      include 'common.h'
*
      real*8 zm,zm1,zm2,vvr,vvt,zmi,vr2,vt2,xbar,xbarh,xbar10,xbet,
     &       xbeth,xbet10,frac,axx,sxx,rms,rgs,rwd,rnsbh,rb,rs,
     &       r13,r23,rr,rxx,avms10,avmsh,avmst,avwd10,avwdh,avwdt,
     &       avns10,avnsh,avnst,ams10,amsh,amst,awd10,awdh,awdt,ans10,
     &       ansh,anst,rms10,rmsh,rmst,rwd10,rwdh,rwdt,rns10,rnsh,
     &       rnst,smms,smwd,smnsbh,smgs,smb,sms,avgs10,avgsh,avgst,
     &       ags10,agsh,agst,rgs10,rgsh,rgst,bm10,bmh,bmt,ba10,bah,
     &       bat,br10,brh,brt,xrs10,xrsh,xrst,amas10,amash,r10pc,
     &       rmsl,smsl,smsl10,rmslh,rmsl10h,zmilh,zmil10h,sm1,sm2,
     &       zm3,smsl10h,smslh

*
      integer i,j,n,im,imk,nxx,nms,nwd,nnsbh,ngs,nb,ns,imm,imo,imn,
     &        ids,nmsl,nmsl10,nmsll,idb,id1,id2,ik1,ik2,ibir,nbh
*
      dimension axx(3),sxx(3),frac(3),nxx(3),rxx(3)
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
*
*       determine the lagrangian radii and mean radial and tangential
*       velocities in the lagrangian shells
*
      n = nt
*
      frac(1) = 0.1d0
      frac(2) = 0.5d0
      frac(3) = 1.0d0
*
*     calculate the taotal masses for main sequence, giants, white dwarf, 
*     neutron + black holes stars
*
      r10pc = 10.d0
      smnsbh = 0.d0
      smwd = 0.d0
      smms = 0.d0
      smgs = 0.d0
      sms = 0.d0
      smb = 0.d0
      smsl = 0.d0
      smsl10 = 0.d0
      rs = 0.d0
      rms = 0.d0
      rgs = 0.d0
      rwd = 0.d0
      rnsbh = 0.d0
      rb = 0.d0
      rmsl = 0.d0
      nnsbh = 0
      nwd = 0
      nms = 0
      ngs = 0
      ns = 0
      nb = 0
      nmsl = 0
      nmsll = 0
      nmsl10 = 0
cAdded DCH 3/8/6
      nbh = 0
*
      do 10 i = 1,n
         im = iname(i)
         if(r(i).gt.rtid) go to 10
cchanged following ftncheck
c         if(ikind(im1).lt.0) then
c            print*,'collabr: i,im1,ikind = ',i,im1,ikind(im1)
         if(ikind(im).lt.0) then
            print*,'collabr: i,im,ikind = ',i,im,ikind(im)
            ikind(im) = -ikind(im)
         endif
         if(ikind(im).eq.1.or.ikind(im).ge.3) then
           ids = names(im)
           sms = sms + body(im)
           ns = ns + 1
           rs = r(i)
           call get_ss_type(ids,imk)         
           if(imk.le.1) then
	     smms = smms + body(im)
             nms = nms + 1
	     rms = r(i)
           endif
           if(imk.le.9) then
	     if(body(im)*zmbar.ge.0.5d0) then
	       smsl = smsl + body(im)
	       nmsl = nmsl + 1
	       rmsl = r(i)
	       if(r(i)*rbar/rtidkg.le.r10pc) then
	         smsl10 = smsl10 + body(im)
	         nmsl10 = nmsl10 + 1
	       endif
	     else
	        nmsll = nmsll + 1
	     endif
	   endif
	   if(imk.ge.2.and.imk.le.9) then
	     smgs = smgs + body(im)
	     ngs = ngs + 1
	     rgs = r(i)
	   endif
           if(imk.ge.10.and.imk.le.12) then
             smwd = smwd + body(im)
             nwd = nwd + 1
             rwd = r(i)
           endif
	   if(imk.ge.13.and.imk.le.14) then
	     smnsbh = smnsbh + body(im)
	     nnsbh = nnsbh + 1
	     rnsbh = r(i)
             if (imk.eq.14) nbh = nbh + 1
	   endif
	 else
           smb = smb + body(im)
           nb = nb + 1
           rb = r(i)
           idb = nameb(im)
           call get_loid(idb,id1)
           call get_hiid(idb,id2)
           call get_ss_type(id1,ik1)
           call get_ss_type(id2,ik2)
           call get_loid_mass(idb,sm1)
           call get_hiid_mass(idb,sm2)
c           print*,'i,im,idb,id1,id2,ik1,ik2,sm1,sm2=',i,im,idb,id1,
c     &             id2,ik1,ik2,sm1,sm2
           call flush(6)
           if(ik1.le.9) then
             if(sm1.ge.0.5d0) then
               smsl = smsl + sm1/zmbar
               nmsl = nmsl + 1
               rmsl = r(i)
               if(r(i)*rbar/rtidkg.le.r10pc) then
                 smsl10 = smsl10 + sm1/zmbar
                 nmsl10 = nmsl10 + 1
               endif
             else
               nmsll = nmsll + 1
             endif
           endif
           if(ik2.le.9) then
             if(sm2.ge.0.5d0) then
               smsl = smsl + sm2/zmbar
               nmsl = nmsl + 1
               rmsl = r(i)
               if(r(i)*rbar/rtidkg.le.r10pc) then
                 smsl10 = smsl10 + sm2/zmbar
                 nmsl10 = nmsl10 + 1
               endif
             else
               nmsll = nmsll + 1
             endif
           endif
         endif
* 	
 10   continue
* 
*
*     main sequence stars
*
*     Mass of main sequence stars > 0.5 Mo
*
      zmilh = 0.5d0*smsl
      zmil10h = 0.5d0*smsl10      
      smslh = zmilh
      smsl10h = zmil10h
      imo = 0
      imn = 1
      im = 0
      zm = 0.d0
 12   im = im + 1
      i = iname(im)
      if(ikind(i).ne.2) then
        if(r(im).gt.rtid) go to 12
        if(body(i)*zmbar.lt.0.5d0) go to 12
        ids = names(i)
        call get_ss_type(ids,imk)
        if(imk.gt.9) go to 12
        imo = imn
        imn = im
        zm2 = body(i)
        zm = zm + zm2
*
      else
*
        ibir = 0
        if(r(im).gt.rtid) go to 12
        idb = nameb(i)
        call get_loid(idb,id1)
        call get_hiid(idb,id2)
        call get_ss_type(id1,ik1)
        call get_ss_type(id2,ik2)
        call get_loid_mass(idb,sm1)
        call get_hiid_mass(idb,sm2)
        if(ik1.le.9) then
          zm3 = 0.d0 
          if(sm1.lt.0.5d0) go to 13
          imo = imn
          imn = im
          zm2 = sm1/zmbar
          zm = zm + zm2
          ibir = 1
        endif
 13     continue       
        if(ik2.le.9) then
          if(sm2.lt.0.5d0) go to 14
          if(ibir.eq.1) zm3 = zm2
          imo = imn
          imn = im
          zm2 = sm2/zmbar
          zm = zm + zm2
        endif
 14     continue 
      endif
*      
      if(zm.lt.zmilh) go to 12
      if(ibir.eq.1) zm2 = zm2 + zm3
      zm1 = zm - zm2
      if(imo.gt.0) then
        r13 = r(imo)**3
      else
        r13 = 0.d0
      endif
      if(imn.gt.imo) then  
        r23 = r(imn)**3
      else
        r13 = 0.d0
      endif  
      rr = (zmilh - zm1)/(zm - zm1)
      rr = r13 + rr*(r23 - r13)
      rmslh = rr**one3
*
*     Mass of main sequence stars > 0.5Mo inside 10pc
*
      imo = 0
      imn = 1
      im = 0
      zm = 0.d0
 15   im = im + 1
      i = iname(im)
      if(r(im)*rbar/rtidkg.gt.r10pc) go to 15
      if(ikind(i).ne.2) then
        if(r(im).gt.rtid) go to 15
        if(body(i)*zmbar.lt.0.5d0) go to 15
        ids = names(i)
        call get_ss_type(ids,imk)
        if(imk.gt.9) go to 15
        imo = imn
        imn = im
        zm2 = body(i)
        zm = zm + zm2
*
      else
*
        ibir = 0
        if(r(im).gt.rtid) go to 15
        idb = nameb(i)
        call get_loid(idb,id1)
        call get_hiid(idb,id2)
        call get_ss_type(id1,ik1)
        call get_ss_type(id2,ik2)
        call get_loid_mass(idb,sm1)
        call get_hiid_mass(idb,sm2)
        if(ik1.le.9) then 
          zm3 = 0.d0
          if(sm1.lt.0.5d0) go to 16
          imo = imn
          imn = im
          zm2 = sm1/zmbar
          zm = zm + zm2
          ibir = 1
        endif
 16     continue       
        if(ik2.le.9) then
          if(sm2.lt.0.5d0) go to 17
          if(ibir.eq.1) zm3 = zm2
          imo = imn
          imn = im
          zm2 = sm2/zmbar
          zm = zm + zm2
        endif
 17     continue 
      endif
*      
      if(zm.lt.zmil10h) go to 15
      if(ibir.eq.1) zm2 = zm2 + zm3
      zm1 = zm - zm2
      if(imo.gt.0) then
        r13 = r(imo)**3
      else
        r13 = 0.d0
      endif
      if(imn.gt.imo) then
        r23 = r(imn)**3  
      else
        r13 = 0.d0
      endif  
      rr = (zmil10h - zm1)/(zm - zm1)
      rr = r13 + rr*(r23 - r13)
      rmsl10h = rr**one3
*
*  All main sequence stars
*
      imo = 0
      imn = 1
      im = 0
      imm = 0
      zm = 0.0d0
      vvr = 0.0d0
      vvt = 0.0d0
*         
      do 20 j=1,3
         sxx(j) = 0.d0
         axx(j) = 0.d0
         rxx(j) = 0.d0
         nxx(j) = 0
         zmi = frac(j)*smms
   30    im = im + 1
         i = iname(im)
	 if(ikind(i).eq.2) go to 30
	 if(r(im).gt.rtid) go to 30
         ids = names(i)
         call get_ss_type(ids,imk)
         if(imk.gt.1) go to 30
*	 
	 imm = imm + 1
         imo = imn
	 imn = im
         zm2 = body(i)
         zm = zm + zm2
         vr2 = vr(i) * vr(i)
         vt2 = vt(i) * vt(i)
         vvr = vvr + vr2
         vvt = vvt + vt2
*
         if(imm.eq.nms) then
           vr2 = 0.0d0
           vt2 = 0.0d0
           zm2 = 0.d0
	   imm = imm + 1
           go to 40
         endif
*
         if (zm.lt.zmi) go to 30
   40    zm1 = zm - zm2
         if(j.le.2) then
           if(imo.gt.0) then
             r13 = r(imo)**3
           else
             r13 = 0.d0
           endif
           if(imn.gt.imo) then
             r23 = r(imn)**3  
           else
             r13 = 0.d0
           endif  
           rr = (zmi - zm1)/(zm - zm1)
           rr = r13 + rr*(r23 - r13)
           rxx(j) = rr**one3
         else
           rxx(j) = r(im)
         endif
*
         sxx(j) = zm1
         vvr = vvr - vr2
         vvt = vvt - vt2
         axx(j) = 1.0d0 - 0.5d0*vvt/vvr
         im = im - 1
         nxx(j) = imm - 1
         imm = imm - 1
         zm = zm1
*
 20   continue
*
*       calculate average lagrangian radii, masses and anisotropies 
*       for 10%, 50% and 100% lagrangian radii for main sequence stars
*
      avms10 = sxx(1)/float(nxx(1))
      avmsh = sxx(2)/float(nxx(2))
      avmst = smms/float(nxx(3))
      ams10 = axx(1)
      amsh = axx(2)
      amst = axx(3)
      rms10 = rxx(1)
      rmsh = rxx(2)
      rmst = rms
*
*
*       giant  stars
*
      if(ngs.gt.20) then
*
c
c      write(6,*) ' giant stars'
c 
        imn = 1
        imo = 0
        im = 0
	imm = 0
        zm = 0.0d0
        vvr = 0.0d0
        vvt = 0.0d0
*
        do 45 j=1,3
           sxx(j) = 0.d0
           axx(j) = 0.d0
           rxx(j) = 0.d0
           nxx(j) = 0
           zmi = frac(j)*smgs
 46        im = im + 1
cAdded DCH 1/8/6
           if (im.gt.nt) then
              print*,'im > nt in collabr.f(5)',im,nt,zm,zmi,j,smgs
              goto 47
           endif
           i = iname(im)
           if(ikind(i).eq.2) go to 46
           if(r(im).gt.rtid) go to 46
           ids = names(i)
           call get_ss_type(ids,imk)
           if(imk.lt.2.or.imk.gt.9) go to 46
*
	   imm = imm + 1
           imo = imn
           imn = imm
           zm2 = body(i)
           zm = zm + zm2
           vr2 = vr(i) * vr(i)
           vt2 = vt(i) * vt(i)
           vvr = vvr + vr2
           vvt = vvt + vt2
*
           if(imm.eq.ngs) then
             vr2 = 0.0d0
             vt2 = 0.0d0
             zm2 = 0.d0
             imm = imm  + 1
             go to 47
           endif
*
           if (zm.lt.zmi) go to 46
 47        zm1 = zm - zm2
           sxx(j) = zm1
           if(j.le.2) then
             if(imo.gt.0) then
               r13 = r(imo)**3
             else
               r13 = 0.d0
             endif
             if(imn.gt.imo) then
               r23 = r(imn)**3  
             else
               r13 = 0.d0
             endif
             rr = (zmi - zm1)/(zm - zm1)
             rr = r13 + rr*(r23 - r13)
             rxx(j) = rr**one3
           else
             rxx(j) = r(im)
           endif
           vvr = vvr - vr2
           vvt = vvt - vt2
           axx(j) = 1.0d0 - 0.5d0*vvt/vvr
           im = im - 1
           nxx(j) = imm - 1
           zm = zm1
           imm = imm - 1
*  
 45     continue
*
*       calculate average lagrangian radii, masses and anisotropies 
*       for 10%, 50% and 100% lagrangian radii for giant stars
*
        avgs10 = sxx(1)/float(nxx(1))
        avgsh = sxx(2)/float(nxx(2))
        avgst = smgs/float(nxx(3))
        ags10 = axx(1)
        agsh = axx(2)
        agst = axx(3)
        rgs10 = rxx(1)
        rgsh = rxx(2)
        rgst = rgs
*
      else
*      
        avgs10 = 0.d0
        avgsh = 0.d0
        avgst = 0.d0
        ags10 = 0.d0
        agsh = 0.d0
        agst = 0.d0
        rgs10 = 0.d0
        rgsh = 0.d0
        rgst = rgs
*        
      endif

*
*      white dwarfs
*
      if(nwd.gt.20) then
*
c
c      write(6,*) ' wwwwwwwwwwwwwddddddddddddd'
c 
        imn = 1
        imo = 0
        im = 0
	imm = 0
        zm = 0.0d0
        vvr = 0.0d0
        vvt = 0.0d0
*
        do 50 j=1,3
           sxx(j) = 0.d0
           axx(j) = 0.d0
           rxx(j) = 0.d0
           nxx(j) = 0
           zmi = frac(j)*smwd
 60        im = im + 1
cAdded DCH 1/8/6
           if (im.gt.nt) then
              print*,'im > nt in collabr.f(4)',im,nt,zm,zmi,j,smwd
              goto 70
           endif
           i = iname(im)
           if(ikind(i).eq.2) go to 60
           if(r(im).gt.rtid) go to 60
           ids = names(i)
           call get_ss_type(ids,imk)
           if(imk.lt.10.or.imk.gt.12) go to 60
*
	   imm = imm + 1
           imo = imn
           imn = imm
           zm2 = body(i)
           zm = zm + zm2
           vr2 = vr(i) * vr(i)
           vt2 = vt(i) * vt(i)
           vvr = vvr + vr2
           vvt = vvt + vt2
*
           if(imm.eq.nwd) then
             vr2 = 0.0d0
             vt2 = 0.0d0
             zm2 = 0.d0
             imm = imm  + 1
             go to 70
           endif
*
           if (zm.lt.zmi) go to 60
 70        zm1 = zm - zm2
           sxx(j) = zm1
           if(j.le.2) then
             if(imo.gt.0) then
               r13 = r(imo)**3
             else
               r13 = 0.d0
             endif
             if(imn.gt.imo) then
               r23 = r(imn)**3  
             else
               r13 = 0.d0
             endif  
             rr = (zmi - zm1)/(zm - zm1)
             rr = r13 + rr*(r23 - r13)
             rxx(j) = rr**one3
c             print*,'imm,imo,imn,im,j,r13,r23,rr,zm,zm1,zm2,zmi = ',
c     &              imm,imo,imn,im,j,r13,r23,rr,zm,zm1,zm2,zmi
           else
             rxx(j) = r(im)
           endif
           vvr = vvr - vr2
           vvt = vvt - vt2
           axx(j) = 1.0d0 - 0.5d0*vvt/vvr
           im = im - 1
           nxx(j) = imm - 1
           zm = zm1
           imm = imm - 1
*  
 50     continue
*
*       calculate average lagrangian radii, masses and anisotropies 
*       for 10%, 50% and 100% lagrangian radii for white dwarfs
*
        avwd10 = sxx(1)/float(nxx(1))
        avwdh = sxx(2)/float(nxx(2))
        avwdt = smwd/float(nxx(3))
        awd10 = axx(1)
        awdh = axx(2)
        awdt = axx(3)
        rwd10 = rxx(1)
        rwdh = rxx(2)
        rwdt = rwd
*
      else
*      
        avwd10 = 0.d0
        avwdh = 0.d0
        avwdt = 0.d0
        awd10 = 0.d0
        awdh = 0.d0
        awdt = 0.d0
        rwd10 = 0.d0
        rwdh = 0.d0
        rwdt = rwd
*        
      endif
*
*       neutron stars + black holes
*
      if(nnsbh.gt.20) then
*
c
c      write(6,*) ' nnnnnnnnnnnnnnnnnnsssssssssssssssss'
c
        imn = 1
        imo = 0
        im = 0
	imm = 0
        zm = 0.0d0
        vvr = 0.0d0
        vvt = 0.0d0
*
        do 80 j=1,3
           sxx(j) = 0.d0
           axx(j) = 0.d0
           rxx(j) = 0.d0
           nxx(j) = 0
           zmi = frac(j)*smnsbh
 90        im = im + 1
cAdded DCH 1/8/6
           if (im.gt.nt) then
              print*,'im > nt in collabr.f(3)',im,nt,zm,zmi,j,smnsbh
              goto 100
           endif
           i = iname(im)
           if(ikind(i).eq.2) go to 90
           if(r(im).gt.rtid) go to 90
           ids = names(i)
           call get_ss_type(ids,imk)
           if(imk.lt.13.or.imk.gt.14) go to 90
*
	   imm = imm + 1 
	   imo = imn
           imn = imm
           zm2 = body(i)
           zm = zm + zm2
           vr2 = vr(i) * vr(i)
           vt2 = vt(i) * vt(i)
           vvr = vvr + vr2
           vvt = vvt + vt2
*
           if(imm.eq.nnsbh) then
             vr2 = 0.0d0
             vt2 = 0.0d0
             zm2 = 0.d0
             imo = imm + 1
             go to 100
           endif
*
           if (zm.lt.zmi) go to 90
 100       zm1 = zm - zm2
           sxx(j) = zm1
           if(j.le.2) then
             if(imo.gt.0) then
               r13 = r(imo)**3
             else
               r13 = 0.d0
             endif
             if(imn.gt.imo) then
               r23 = r(imn)**3  
             else
               r13 = 0.d0
             endif  
             rr = (zmi - zm1)/(zm - zm1)
             rr = r13 + rr*(r23 - r13)
             rxx(j) = rr**one3
           else
             rxx(j) = r(im)
           endif
	   rxx(j) = r(imo)
           vvr = vvr - vr2
           vvt = vvt - vt2
           axx(j) = 1.0d0 - 0.5d0*vvt/vvr
           im = im - 1
           nxx(j) = imm - 1
           zm = zm1
           imm = imm - 1
*  
 80     continue
*
*       calculate average lagrangian radii, masses and anisotropies 
*       for 10%, 50% and 100% lagrangian radii for neutron stars and
*       black holes
*
        avns10 = sxx(1)/float(nxx(1))
        avnsh = sxx(2)/float(nxx(2))
        avnst = smnsbh/float(nxx(3))
        ans10 = axx(1)
        ansh = axx(2)
        anst = axx(3)
        rns10 = rxx(1)
        rnsh = rxx(2)
        rnst = rnsbh
*
      else
*      
        avns10 = 0.d0
        avnsh = 0.d0
        avnst = 0.d0
        ans10 = 0.d0
        ansh = 0.d0
        anst = 0.d0
        rns10 = 0.d0
        rnsh = 0.d0
        rnst = rnsbh
*        
      endif
*
*     all stars
*
c
c      write(6,*) 'alllllll starssssssss'
c
      im = 0
      imo = 0
      imn = 1
      imm = 0
      zm = 0.0d0
      vvr = 0.0d0
      vvt = 0.0d0
*
      do 110 j=1,3
         sxx(j) = 0.d0
         axx(j) = 0.d0
         rxx(j) = 0.d0
         nxx(j) = 0
         zmi = frac(j)*sms
 120     im = im + 1
cAdded DCH 1/8/6
         if (im.gt.nt) then
            print*,'im > nt in collabr.f',im,nt,zm,zmi,j,sms
            goto 130
         endif
         i = iname(im)
         if(ikind(i).eq.2) go to 120
         if(r(im).gt.rtid) go to 120
         imm = imm + 1
         imo = imn
         imn = im
         zm2 = body(i)
         zm = zm + zm2
         vr2 = vr(i) * vr(i)
         vt2 = vt(i) * vt(i)
         vvr = vvr + vr2
         vvt = vvt + vt2
*
         if(imm.eq.ns) then
	   vr2 = 0.0d0
           vt2 = 0.0d0
           zm2 = 0.d0
           imm = imm + 1
           go to 130
         endif
*
         if (zm.lt.zmi) go to 120
 130     zm1 = zm - zm2
         if(j.le.2) then
           if(imo.gt.0) then
             r13 = r(imo)**3
            else
              r13 = 0.d0
            endif
            if(imn.gt.imo) then
              r23 = r(imn)**3  
            else
              r13 = 0.d0
           endif  
           rr = (zmi - zm1)/(zm - zm1)
           rr = r13 + rr*(r23 - r13)
           rxx(j) = rr**one3
         else
           rxx(j) = r(im)
         endif
         sxx(j) = zm1
         vvr = vvr - vr2
         vvt = vvt - vt2
         axx(j) = 1.0d0 - 0.5d0*vvt/vvr
         im = im - 1
         nxx(j) = imm - 1
         zm = zm1
         imm = imm - 1
*
 110  continue
*
*       calculate average masses and anisotropies for 10%, 50% and 100%
*       lagrangian radii
*
      xbar10 = sxx(1)/float(nxx(1)) 
      xbarh = sxx(2)/float(nxx(2))
      xbar = sms/float(nxx(3))
      xbet10 = axx(1)
      xbeth = axx(2)
      xbet = axx(3)
      xrs10 = rxx(1)
      xrsh = rxx(2)
      xrst = rs
*
*     only binaries
*
      if(nb.gt.20) then
*
c
c      write(6,*) 'only binaries'
c
        im = 0
        imn = 1
        imm = 0
        imo = 0 
        zm = 0.0d0
        vvr = 0.0d0
        vvt = 0.0d0
*
        do 140 j=1,3
           sxx(j) = 0.d0
           axx(j) = 0.d0
           rxx(j) = 0.d0
           nxx(j) = 0
           zmi = frac(j)*smb
 150       im = im + 1
cAdded by DCH 1/8/6
           if (im.gt.nt) then
              print*,'im > nt in collabr.f(2)',im,nt,zm,zmi,j,smb
              goto 160
           endif
           i = iname(im)
cChanged by DCH 1/8/6
c           if(ikind(i).eq.2) go to 150
           if(ikind(i).ne.2) go to 150
           if(r(im).gt.rtid) go to 150
           imm = imm + 1
           imo = imn
           imn = imm
           zm2 = body(i)
           zm = zm + zm2
           vr2 = vr(i) * vr(i)
           vt2 = vt(i) * vt(i)
           vvr = vvr + vr2
           vvt = vvt + vt2
*
           if(imm.eq.nb) then
             vr2 = 0.0d0
             vt2 = 0.0d0
             zm2 = 0.d0
             imm = imm + 1
             go to 160
           endif
*
           if (zm.lt.zmi) go to 150
 160       zm1 = zm - zm2
           sxx(j) = zm1
           if(j.le.2) then
             if(imo.gt.0) then
               r13 = r(imo)**3
             else
               r13 = 0.d0
             endif
             if(imn.gt.imo) then
               r23 = r(imn)**3  
             else
               r13 = 0.d0
             endif  
             rr = (zmi - zm1)/(zm - zm1)
             rr = r13 + rr*(r23 - r13)
             rxx(j) = rr**one3
           else
              rxx(j) = r(im)
           endif
           vvr = vvr - vr2
           vvt = vvt - vt2
           axx(j) = 1.0d0 - 0.5d0*vvt/vvr
           im = im - 1
           nxx(j) = imm - 1
           zm = zm1
           imm = imm - 1
*
 140    continue
*
*       calculate average masses and anisotropies for 10%, 50% and 100%
*       lagrangian radii
*
        bm10 = sxx(1)/float(nxx(1)) 
        bmh = sxx(2)/float(nxx(2))
        bmt = smb/float(nxx(3))
        ba10 = axx(1)
        bah = axx(2)
        bat = axx(3)
        br10 = rxx(1)
        brh = rxx(2)
        brt = rb
      else
        bm10 = 0.d0
        bmh = 0.d0
        bmt = 0.d0
        ba10 = 0.d0
        bah = 0.d0
        bat = 0.d0
        br10 = 0.d0
        brh = 0.d0
        brt = rb
      endif 
*
*     whole system
*
      zm = 0.d0
      amas10 = 0.d0
      amash = 0.d0
*      
      do 200 i = 1,nt
         im = iname(i)
         zm = zm + body(im)
         if(zm.le.0.1d0*smt) amas10 = 0.1d0*smt/float(i)
         if(zm.le.0.5d0*smt) amash = 0.5d0*smt/float(i)
 200  continue
*
*
      return
*
      end
*
*
*
*
       subroutine core
*
*
*     calculate core parameters: rc,vc,nc,smc,roc
*     -------------------------------------------
*
*
      include 'common.h'
*
      real*8 sx,v2x
*
      integer i,n,nx,im,ncut
*
*
      n = nt
*     
*       determination number of stars to calculate the central parameters
*
      ncut = 0.01*n
      if(ncor.le.ncut) then
        nx = ncor
      else
        if(ncut.gt.nmin) then
          nx = ncut
        else
          nx = nmin
        endif
      endif    
*
*
*       determination of the central parameters
*
      sx = 0.0d0
      v2x = 0.0d0
*
      do 10 i=1,nx
         im = iname(i)
         sx = sx + body(im)
   10    v2x = v2x + body(im)*(vr(im)*vr(im) + vt(im)*vt(im))
*
      v2x = v2x/sx
      roc = 3.0d0*sx/(4.0d0*pi*r(nx)**3)
      rc = 3.0d0*v2x/(4.0d0*pi*roc)
      vc = sqrt(v2x)
      rc = sqrt(rc)
*
      i = 0
      sx = 0.0d0
*
   20 continue
      i = i + 1
      im = iname(i)
      sx = sx + body(im)
      if(rc.gt.r(i)) go to 20
      smc = sx - body(im)
      nc = i-1
*
*
      return
*
      end
*
*
*
*
       subroutine data
*
*
*       initial conditions.
*       -------------------
*
      include 'common.h'
*
*
      real*8  b(8),alpha1,const,fmnn,fm1,fmn,fmi,zmbar1,cm(3),
     &        cmdot(3),ri,g1,g2,g3,g4,range,lambda,chi,rrm,bbmax,
     &        ssmax,delta,eta,ebirth,ein,pbirth,pminl,rperi,sm1,
     &        sm2,rho,rhop,qbirth,qin,sm2in,sum1,sum2,sum2p,sm12,
     &        sm12in,smtp,rperi1,auro,cmin,cmax,rpauro,ftot,g5,
     &        g6,g7,g8,pin,q,qmin,smts,xx,zm,zm1,zm2,fbrake,
     &        standalone_rzamsf
* 
      real*4 ran2
*
      integer i,k,n,ntt,nss,nbb,ibb,ii,ixb,ibx,im
*
      data g1,g2,g3,g4 /0.19d0,1.55d0,0.05d0,0.6d0/
      DATA g5,g6,g7,g8 /0.75,0.04,0.25,1.04/
      data lambda,chi,delta,eta,pminl/28.d0,0.75d0,45.d0,2.5d0,1.d0/
*
*     scalling factor from AU to Ro, AU = 1.496E+13 cm, Ro = 6.955E+10 cm       
*     auro = AU/Ro = 2.151E+02                                                  
*
      auro = 2.151e+02
*      
      if(fracb.eq.0.d0) then
        ntt = nt
        nss = nt
        nbb = 0
        ixb = 0
        range = 0.d0
      else
        ntt = (1.0d0 + fracb)*nt + 0.5d0
        nss = (1.0d0 -fracb)*nt + 0.5d0
        nbb = fracb*nt + 0.5d0
        if(amin.gt.0.d0) then
          range = abs(amax)/amin
        else
          if(amin.eq.0.d0) then
*  
*     M67 case amin=2*(R1+R2). R1=R2 for 0.1Mo ==> amin=4*R1 for m=0.1Mo
*       
            amin = 4.d0*standalone_rzamsf(0.1d0)
            range = abs(amax)/amin
          endif  
        endif
      endif
*
      nss0 = nss
      n = nt
*
      if(iprint.eq.0)
     & write (6,20) alphal,alphah,fbrake,body1,bodyn,fracb,amin,amax,
     &              ntt,nt,nss,nbb
   20 format (1x,'mass function:','  alphal,alphah,fbrake =',1p3d10.4,
     &        '  body1 =',1pd10.4,'  bodyn =',1pd10.4,'  primordial' 
     &        ,' binaries: ','  fracb =',1pd10.4,'  amin =',1pd10.4,
     &        '  amax =',1pd10.4,'  ntt =',i8,'  nt =',i8,'  nss =',
     &        i8,'  nbb =',i8)
*
*       include the case of equal masses (defined by alphal = -1)
*
      if (alphal.eq.-1.0d0) then
*
          smt = 0.d0
          ibb = 0
          ii = 0
*
          do 10 i = 1,ntt
             body(i) = 1.0d0
             smt = smt + body(i)
             if(i.le.nss) then 
               ii = ii + 1
               ikind(i) = 1
               names(i) = i - 1
               nameb(i) = 0
               iname(i) = i
             else
               ibb = ibb + 1
               if(ibb.eq.1) then
                 go to 10
               else
                 ii = ii + 1
                 ikind(ii) = 2
                 nameb(ii) = i - 2
                 names(ii) = 0
                 iname(ii) = ii 
                 ixb = ii - nss0
                 nwhich(ii) = ixb
                 nbinar(ixb) = ii
                 bin(ixb,1) = body(i-1)
                 bin(ixb,2) = body(i)
                 fmi = ran2(irun)
                 fmn = fmi*log10(range)
                 bin(ixb,3) = abs(amax)/10.d0**fmn
                 fmi = ran2(irun)
                 bin(ixb,4) = sqrt(fmi)
                 bin(ixb,6) = 0.d0
                 bin(ixb,7) = 1.0d8*float(i-2) + float(i-1)
                 bin(ixb,8) = 0.d0
                 body(ii) = body(i-1) + body(i)
                 ibb = 0
               endif
             endif
*
 10       continue
*
          nbin3 = nbb
          if(ixb.ne.nbb) then
            write(6,*) 'wrong number of binaries:  nbb,ixb = ',nbb,ixb
          endif
*
          smt = float(ntt)
          go to 60
*
      else
*
*       generate IMF
*
        if(alphal.eq.alphah) then
*
* generate a power-law mass function with exponent alphal=alphah
* IMF is not a broken power-law
*
          alpha1 = alphal - 1.0d0
          const=alpha1/(bodyn**(-alpha1)-body1**(-alpha1))
          fmnn=const*bodyn**(-alphal)
*
          sum1 = 0.d0
          sum2 = 0.d0
          smts = 0.d0
          ssmax = 0.d0
          bbmax = 0.d0
          smt = 0.d0
          ibb = 0
          ii = 0
*
          do 30 i = 1,ntt
   40        fm1=ran2(irun)
             body(i)=bodyn+(body1-bodyn)*fm1
             fmi=const*body(i)**(-alphal)
             fmn=ran2(irun)*fmnn
             if(fmn.gt.fmi) go to 40
             smt = smt + body(i)
             if(i.le.nss) then
               sum1 = sum1 + body(i)
               smts = smts + body(i)
               if(body(i).gt.ssmax) ssmax = body(i) 
               ii = ii + 1
               ikind(i) = 1
               names(i) = i - 1
               nameb(i) = 0
               iname(i) = i
             else
               ibb = ibb + 1
               sum2 = sum2 + body(i)
               if(ibb.eq.1) then
                 go to 30
               else
                 ii = ii + 1
                 ikind(ii) = 2
                 nameb(ii) = i - 2
                 names(ii) = 0
                 iname(ii) = ii 
                 ixb = ii - nss0
                 nwhich(ii) = ixb
                 nbinar(ixb) = ii
                 if(body(i-1).ge.body(i)) then
                   bin(ixb,1) = body(i-1)
                   bin(ixb,2) = body(i)
                   bin(ixb,7) = 1.0d8*float(i-2) + float(i-1)
                 else
                   bin(ixb,1) = body(i)
                   bin(ixb,2) = body(i-1)
                   bin(ixb,7) = 1.0d8*float(i-1) + float(i-2)
                 endif
                 fmi = ran2(irun)
                 fmn = fmi*log10(range)
                 bin(ixb,3) = abs(amax)/10.d0**fmn
                 fmi = ran2(irun)
                 bin(ixb,4) = sqrt(fmi)
                 bin(ixb,5) = body(i-1)*body(i)/2.d0/bin(ixb,3)
                 bin(ixb,6) = 0.d0
                 bin(ixb,8) = 0.d0
                 body(ii) = body(i-1) + body(i)
                 if(body(ii).gt.bbmax) bbmax = body(ii)
                 ibb = 0
               endif
             endif
*
 30       continue
*
          print *,'Power-law mass function - no break'
          print *,'Finished singles and binaries; total an average mass'
          print *,'total Ms+Mb, <ms+mb> = ', smt,smt/float(ntt)
          print *,'singles Ms, <ms> = ',sum1,sum1/float(nss)  
          print *,'binaries Mb, <mb>  = ',sum2, sum2/float(nbb)
          print*,'max mass star =',ssmax,'  maxx mass bin =',bbmax
*
          nbin3 = nbb
          if(ixb.ne.nbb) then
            write(6,*) 'wrong number of binaries:  nbb,ixb = ',nbb,ixb
          endif
          zmbar = smt
*
        else
*
*       generate Kroupa mass function (Kroupa et al MN 262. 545)
*
          sum1 = 0.d0
          sum2 = 0.d0
          sum2p = 0.d0
          smt = 0.d0
          smtp = 0.d0
          smts = 0.d0
          ssmax = 0.d0
          bbmax = 0.d0
          cmin = 1.d30
          cmax = 0.d0
          ibb = 0
          ii = 0
*
* IMF is a broken power-law: lower index is alphal, higher index is alphah,
* breake is at mass breakm.
*
          do 35 i = 1,ntt
 37          continue
             fm1 = ran2(irun)
             if(brakem.lt.bodyn) then
               zm = 0.08d0+(g1*fm1**g2 + g3*fm1**g4)/(1.d0 - fm1)**0.58
             else
               fbreak = -(brakem/(alphal - 1.d0))*(1.d0 - (bodyn/brakem)
     &                  **(-alphal + 1.d0))
               ftot = fbreak - (brakem/(alphah - 1.d0))*((body1/brakem)
     &                **(-alphah + 1.d0) - 1.d0)
               xx = fm1*ftot
               if(xx.gt.fbreak) then
                 zm = brakem*(1.d0 - ((alphah - 1.d0)/brakem)*
     &                (xx - fbreak))**(-1.d0/(alphah - 1.d0))
               else
                 zm = brakem*((bodyn/brakem)**(-alphal + 1.d0) + 
     &                (1.d0 - alphal)*xx/brakem)**
     &                (-1.d0/(alphal - 1.d0))
               endif
             endif
*
* masses of the binary components are picked up from the same IMF as
* single stars
* 
             if((zm.ge.bodyn).and.(zm.le.body1)) then
               body(i) = zm
               smt = smt + zm
               smtp = smtp + zm
               if(i.le.nss) then 
                 sum1 = sum1 + zm
                 smts = smts + zm
                 if(body(i).gt.ssmax) ssmax = body(i)
                 ii = ii + 1
                 ikind(i) = 1
                 names(i) = i - 1
                 nameb(i) = 0
                 iname(i) = i
               else
                 ibb = ibb + 1
                 sum2 = sum2 + zm
                 sum2p = sum2p + zm
                 if(ibb.eq.1) then
                   go to 35
                 else
                   ii = ii + 1
                   ikind(ii) = 2
                   nameb(ii) = i - 2
                   names(ii) = 0
                   iname(ii) = ii 
                   ixb = ii - nss0
                   nwhich(ii) = ixb
                   nbinar(ixb) = ii
                   if(body(i-1).ge.body(i)) then
                     bin(ixb,1) = body(i-1)
                     bin(ixb,2) = body(i)
                     bin(ixb,7) = 1.0d8*float(i-2) + float(i-1)
                   else
                     bin(ixb,1) = body(i)
                     bin(ixb,2) = body(i-1)
                     bin(ixb,7) = 1.0d8*float(i-1) + float(i-2)
                   endif
*
                   if(ikroupa.eq.0) then
*
*   eigenevolution and a feeding algorithm Kroupa 1995 (MNRAS, 277, 1507)                 
*
                     ibb = 0
*
*                 always sm1 > sm2
*                   
                     sm1 = bin(ixb,1)
                     sm2 = bin(ixb,2)
                     sm12 = sm1 + sm2
                     qbirth = sm2/sm1
 355                 fmi = ran2(irun)
                     pbirth = sqrt(delta*(exp(2.d0*fmi/eta) - 1.d0))
                     pbirth = 10.d0**(pbirth + pminl)                                   
                     fmi = ran2(irun)
                     ebirth = sqrt(fmi)
                     rperi1 = (sm12*(pbirth/365.25d0)**2.d0)**one3
                     rperi = rperi1*(1.d0 - ebirth)
                     rho = (lambda/auro/rperi)**chi
c                     print*,'rho = ',rho
                     ein = ebirth*exp(-rho)
                     if(rho.le.1.d0) then
                       rhop = rho
                     else
                       rhop = 1.d0
                     endif
                     qin = qbirth + (1.d0 - qbirth)*rhop
                     sm2in = sm1*qin
                     sm12in = sm1 + sm2in
                     bin(ixb,2) = sm2in
                     sum2p = sum2p - sm2 + sm2in
                     smtp = smtp - sm2 + sm2in
                     pin = pbirth*sqrt(sm12/sm12in)
                     pin = pin*((1.d0 - ebirth)/(1.d0 - ein))**1.5d0
                     rperi1 = (sm12in*(pin/365.25d0)**2.d0)**one3
c                     if(rperi1.gt.50.d0) then
c                       smtp = smtp - sm2in + sm2
c                       sum2p = sum2p - sm2in + sm2
c                       go to 355
c                     endif
                     rperi = rperi1*(1.d0 - ein)
                     rpauro = rperi*auro
                     if(rpauro.lt.5.d0*standalone_rzamsf(sm1)) then
                       sum2p = sum2p - sm2in + sm2
                       smtp = smtp - sm2in + sm2
                       go to 355
                     endif
c                     write(6,987) i,sm1,sm2,sm2in,qbirth,qin,ebirth,
c     &                            ein,pbirth,pin,rperi1,rperi,rho,rhop
c 987                 format(1x,i9,1p13e12.4)
                     bin(ixb,4) = ein
                     bin(ixb,3) = rperi1*auro 
                     bin(ixb,6) = 0.d0
                     bin(ixb,8) = 0.d0
                     body(ii) = sm12in
                     if(rperi1*auro.le.cmin) cmin = rperi1*auro
                     if(rperi1*auro.ge.cmax) cmax = rperi1*auro
                     if(body(ii).gt.bbmax) bbmax = body(ii)
*
                   else                               
*
* binary parameters as for Hurley's M67 model
*
 988                 continue
                     fmi = ran2(irun)
                     fmn = fmi*log10(range)
                     bin(ixb,3) = abs(amax)/10.d0**fmn
                     fmi = ran2(irun)
                     bin(ixb,4) = sqrt(fmi)
                     bin(ixb,6) = 0.d0
                     bin(ixb,8) = 0.d0
                     body(ii) = body(i-1) + body(i)
                     if(bin(ixb,3).le.cmin) cmin = bin(ixb,3)
                     if(bin(ixb,3).ge.cmax) cmax = bin(ixb,3)
                     if(body(ii).gt.bbmax) bbmax = body(ii)
                     ibb = 0
                     sm1 = bin(ixb,1)
                     sm2 = bin(ixb,2)
                     sm2in = sm2
                     qin = sm2/sm1
                     qbirth = qin
                     ein = bin(ixb,4)
                     rperi1 = bin(ixb,3)
                     rperi = rperi1*(1.d0 - ein)
                     if (rperi.lt.(standalone_rzamsf(sm1) +
     &                   standalone_rzamsf(sm2))) go to 988 
                   
                     rhop = rho
                     pin = 365.25d0*(rperi1/auro)**1.5/sqrt(sm1 + sm2)
                     pbirth = pin
                     if(rho.gt.1.d0) rhop = 1.d0
c                     write(6,987) i,sm1,sm2,sm2in,qbirth,qin,ebirth,
c     &                            ein,pbirth,pin,rperi1,rperi,rho,rhop
                   endif
                 endif
               endif
             else
               go to 37
             endif
*
 35       continue
*
          print *,'Broken Power-law mass function'
          print *,'Finished singles and binaries; total an average mass'
          print *,'total Ms+Mb, <ms+mb> = ', smt,smt/float(ntt)
          print *,'singles Ms, <ms> = ',sum1,sum1/float(nss)  
          print *,'binaries Mb, <mb>  = ',sum2, sum2/float(nbb)
          print*,'max mass star =',ssmax,'  maxx mass bin =',bbmax
          if(ikroupa.eq.0) then
            print *,'Corections for eigenevolution and feeding algorit.'
            print *,'new total Ms+Mb, <ms+mb> = ', smtp,smtp/float(ntt)
            print *,'the same singles Ms, <ms> = ',sum1,sum1/float(nss)
            print *,'new binaries Mb, <mb>  = ',sum2p, sum2p/float(nbb)
            print *,'min and max semi major axis = ',cmin, cmax
            print*,'max mass star =',ssmax,'  maxx mass bin =',bbmax
            smt = smtp
          endif
*
* For M67, throw away these binaries and pick up new - binary parameters as
* for Hurley's M67 model
*
          if(imodel.eq.4.and.ikroupa.eq.1) then
            sum1 = 0.d0
            bbmax = 0.d0
            smt = smts
            cmin = 1.d30
            cmax = 0.d0
            do ii = nss+1,nss+nbb
 36            continue
               fm1 = ran2(irun)
               zm = 0.33d0*((1.D0/(fm1**g5 + g6*fm1**g7)) - (fm1**2/g8))
               if (zm.gt.2*body1.or.zm.lt.2*bodyn) go to 36
               qmin = max(bodyn/(zm-bodyn),(zm-body1)/body1)
               q = qmin + (1 - qmin)*ran2(irun)
               zm1 = q*zm/(1+q)
               zm2 = zm - zm1
               if (zm1.gt.body1.or.zm2.gt.body1.
     &             or.zm1.lt.bodyn.or.zm2.lt.bodyn) goto 36
               sum1 = sum1 + zm
               ikind(ii) = 2
               i = nss + (ii-nss)*2 
               nameb(ii) = i - 2
               names(ii) = 0
               iname(ii) = ii
               ixb = ii - nss0
               nwhich(ii) = ixb
               nbinar(ixb) = ii
               body(i-1) = zm1
               body(i) = zm2
               smt = smt + zm
               if(body(i-1).ge.body(i)) then
                 bin(ixb,1) = body(i-1)
                 bin(ixb,2) = body(i)
                 bin(ixb,7) = 1.0d8*float(i-2) + float(i-1)
               else
                 bin(ixb,1) = body(i)
                 bin(ixb,2) = body(i-1)
                 bin(ixb,7) = 1.0d8*float(i-1) + float(i-2)
               endif              
cSelect a; at present this is in solar radii
 375           continue
               fmi = ran2(irun)
               fmn = fmi*log10(range)
cSelect e
               fmi = ran2(irun)
               bin(ixb,4) = sqrt(fmi)
               ebirth = bin(ixb,4)
               bin(ixb,3) = abs(amax)/10.d0**fmn
               if (bin(ixb,3).lt.2.d0*(standalone_rzamsf(body(i-1)) +
     &              standalone_rzamsf(body(i)))) go to 375
               rperi = bin(ixb,3)*(1.d0 - bin(ixb,4))
               if (rperi.lt.(standalone_rzamsf(body(i-1)) +
     &              standalone_rzamsf(body(i)))) go to 375
               rho = (lambda/rperi)**chi
               bin(ixb,4) = exp(-(lambda/rperi)**chi)*bin(ixb,4)
               if(bin(ixb,4).lt.1.d-3) bin(ixb,4) = 1.d-3*ran2(irun)
               bin(ixb,6) = 0.d0
               bin(ixb,8) = 0.d0
               body(ii) = body(i-1) + body(i)
               if(bin(ixb,3).le.cmin) cmin = bin(ixb,3)
               if(bin(ixb,3).ge.cmax) cmax = bin(ixb,3)
               if(body(ii).gt.bbmax) bbmax = body(ii)
               sm1 = bin(ixb,1)
               sm2 = bin(ixb,2)
               sm2in = sm2
               qin = sm2/sm1
               qbirth = qin
               ein = bin(ixb,4)
               rperi1 = bin(ixb,3)
               rperi = rperi1*(1.d0 - ein)
               rhop = rho
               pin = 365.25d0*(rperi1/auro)**1.5/sqrt(sm1 + sm2)
               pbirth = pin
               if(rho.gt.1.d0) rhop = 1.d0
c               write(6,987) sm1,sm2,sm2in,qbirth,qin,ebirth,ein,
c     &                          pbirth,pin,rperi1,rperi,rho,rhop
            enddo
            print *,'Finished generating binaries'
            print *,  'binaries Mb, <mb>  = ',sum1, sum1/float(nbb)
            print*,'max mass star =',ssmax,'  maxx mass bin =',bbmax
            print *,'min and max semi major axis = ',cmin, cmax
          endif
*
          nbin3 = nbb
          if(ixb.ne.nbb) then
            write(6,*) 'wrong number of binaries:  nbb,ixb = ',nbb,ixb
          endif
          zmbar = smt
          write(6,38) nt,ixb,smt,smt/nt
 38       format(1x,' number of objects = ',i8,' number of binaries = '
     &          ,i8,' total mass in sollar units = ',1pd13.6,
     &          ' average mass in sollar units = ',1pd13.6)
*
        endif
      endif
*
*       first scale the masses to <m> = 1
*
      zmbar1 = smt/float(nt)
      do 50 i = 1,n
         im = iname(i)
         if(ikind(im).eq.2) then
           ibx = nwhich(im)
            bin(ibx,1) = bin(ibx,1)/zmbar1
            bin(ibx,2) = bin(ibx,2)/zmbar1
         endif
         body(i) = body(i)/zmbar1
   50 continue
*
      smt = float(nt)
*
   60 if (imodel.ne.1) go to 130
*
*       set up a uniform spherical system
*
      do 70 i = 1,n
   80     b(1) = 0.0
          do 90 k = 1,3
              b(k+1) = 2.0d0*ran2(irun) - 1.0d0
              b(1) = b(1) + b(k+1)**2
   90     continue
          if(b(1).gt.1.0d0) go to 80
  100     b(5) = 0.0d0
          do 110 k = 1,3
              b(k+5) = 2.0d0*ran2(irun) - 1.0d0
              b(5) = b(5) + b(k+5)**2
  110     continue
          if(b(5).gt.1.0d0) go to 100
          do 120 k = 1,3
*
*         density proportional to 1/r**2
*
*             x(i,k) = b(1)*b(k+1)
*
*         constant density
*
              x(i,k) = b(k+1)
*
*         isotropic velocities (magnitude randomized; no radial dependence)
*
              xdot(i,k) = b(k+5)
*
  120     continue
*
      
   70 continue
*
      go to 180
*
cModified for M67
  130 if(imodel.ne.2.and.imodel.ne.4) go to 170
*
*       generate initial conditions from plummer model (a & a 37, 183)
*       Take care for the M67 case  (imodel = 4)
*
         rrm = 20.d0
         b(1) = 1.5d0*twopi/16.0d0
         b(2) = sqrt(smt/b(1))
         if (imodel.eq.4) then
            rtidkg = rplum*b(1)/sqrt(2.d0**(2.d0/3.d0) - 1.d0)
c            rtidkg = rplum*b(1)
            rtid = rtidkg
            rrm = rtid
         endif
*
      do 140 i = 1,n
  150     b(1) = ran2(irun)
          if(b(1).lt.1.0d-10) go to 150
          ri = (b(1)**(-0.6666667) - 1.0d0)**(-0.5)
*
*       reject distant particles
*
          if (ri.gt.rrm) go to 150
          b(2) = ran2(irun)
          b(3) = ran2(irun)
          x(i,3) = (1.0d0 - 2.0d0*b(2))*ri
          x(i,1) = sqrt(ri**2 - x(i,3)**2)*cos(twopi*b(3))
          x(i,2) = sqrt(ri**2 - x(i,3)**2)*sin(twopi*b(3))
  160     b(4) = ran2(irun)
          b(5) = ran2(irun)
          b(6) = b(4)**2*(1.0d0 - b(4)**2)**3.5
          if (0.1d0*b(5).gt.b(6)) go to 160
          b(8) = b(4)*sqrt(2.0d0)/(1.0d0 + ri**2)**0.25
          b(6) = ran2(irun)
          b(7) = ran2(irun)
          xdot(i,3) = (1.0d0 - 2.0d0*b(6))*b(8)
          xdot(i,1) = sqrt(b(8)**2 - xdot(i,3)**2)*cos(twopi*b(7))
          xdot(i,2) = sqrt(b(8)**2 - xdot(i,3)**2)*sin(twopi*b(7))
*
  140 continue
*
*       scale coordinates & velocities to analytical expectation values
*
      b(1) = 1.5d0*twopi/16.0d0
      b(2) = sqrt(smt/b(1))
*
      go to 180
*
 170  continue
*
*       generate initial conditions from king model 
*
      b(1)=1.0d0
      b(2)=1.0d0
*
      call king
*
 180  continue
*
*       accumulate the centre of mass terms
*
      do 190 k = 1,3
          cm(k) = 0.0d0
          cmdot(k) = 0.0d0
  190 continue
*
      do 200 i=1,n
      do 210 k=1,3
           cm(k) = cm(k) + body(i)*x(i,k)
           cmdot(k) = cmdot(k) + body(i)*xdot(i,k)
  210 continue
  200 continue
*
      do 220 i = 1,n
      do 230 k = 1,3
           x(i,k) = x(i,k) - cm(k)/smt
           xdot(i,k) = xdot(i,k) - cmdot(k)/smt
           x(i,k) = b(1)*x(i,k)
           xdot(i,k) = b(2)*xdot(i,k)
  230 continue
  220 continue
*
*     sort all particles according to increase radius
*
      do 240 i=1,n
         r(i)= sqrt(x(i,1)**2 + x(i,2)**2 + x(i,3)**2)
 240  continue
*
      call sort2(n,r,iname)
*
      return
*
      end
*
*
*
*
      real*8 FUNCTION standalone_rzamsf(m)
      implicit none
      real*8 m,mx,a(200)
*
* A function to evaluate Rzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      a(8) = 1.715359d0
      a(9) = 6.597788d0
      a(10) = 10.08855d0
      a(11) = 1.012495d0
      a(12) = 0.07490166d0
      a(13) = 0.01077422d0
      a(14) = 3.082234d0
      a(15) = 17.84778d0
      a(16) = 0.00022582
      mx = SQRT(m)
      standalone_rzamsf = ((a(8)*m**2 + a(9)*m**6)*mx + a(10)*m**11 +
     &          (a(11) + a(12)*mx)*m**19)/
     &         (a(13) + a(14)*m**2 + 
     &          (a(15)*m**8 + m**18 + a(16)*m**19)*mx)
*
      return
      end
      subroutine densi(k,den)
*
*
*       determine density in the vicinity of the 'k' star
*       --------------------------------------------------
*
*
      include 'common.h'
*
      real*8 den,co
*
      integer k,nspan,n,i1,i2
*
*
      n = nt
c      nspan = nminzo/2
c      nspan = 7
      nspan = 4
 10   co = 3.0d0*float(nspan)

      if(k.gt.nspan) then
        i1 = k - nspan
        i2 = k + nspan
*
        if(i2.gt.n) then
          i2 = n
 20       i1 = n - 2*nspan -1 
          if(r(i1).gt.1.d8) then
            nspan = nspan + 1
            co = 3.0d0*float(nspan)
            print*,'dens  nspan,r1,r2,i1,i2=',nspan,r(i1),r(i2),i1,i2
            go to 20
          else  
            den = co/(r(i2)**3 - r(i1)**3)
            return
          endif
        endif  
*           
        if(r(i1).gt.1.d8.or.r(i2).gt.1.d8) then
          nspan = nspan + 1
          print*,'density  nspan,r1,r2,i1,i2=',nspan,r(i1),r(i2),i1,i2
          go to 10
        else
          den = co/(r(i2)**3 - r(i1)**3)
          return
        endif  
*        
      else
*
 30     i2 = 2*nspan + 1
*
        if(r(i2).gt.1.d8) then
          nspan = nspan + 1
          print*,'density  nspan,r2,i2=',nspan,r(i2),i2
          co = 3.0d0*float(nspan)
          go to 30
        else
          den = co/r(i2)**3
          return
        endif
      endif  
*
*
      return
*
      end
*
*
*
*

*
*
      subroutine denbin(nup,nb,denb)
*
*     calculate binary number density in the case of primordial binares
*     --------------------------------------------------------------------
*     density is calculated for binaries which position is inside r(lmax).
*     --------------------------------------------------------------------
*     lmax = nzon(nup)
*     ----------------
*
*
      include 'common.h'
*
cChanged following ftncheck
c      real*8 rb(nmax),denb(nmax),co
      real*8 rb(nmax),denb(nmax),co
*
      integer nup,lmax,i,im1,ib,nspan,nspan2,nmin1,nmax1,nb(nmax)
*
      nspan = 2
      nspan2 = 2*nspan
      co = 0.75d0*(nspan2 + 1)/pi
      lmax = nzst(nup)
      ib = 0
*
      do 10 i = 1,lmax
*
         nb(i) = 0
         im1 = iname(i)
         if(ikind(im1).eq.2) then
           ib = ib + 1
           rb(ib) = r(i)
           nb(i) = ib
           denb(ib) = 0.d0
         endif
*
 10   continue
*
      do 20 i = 1,ib
*
         if(i.le.nspan) then
           nmin1 = 1
           nmax1 = nspan2 + 1
           denb(i) = co/(rb(nmax1)**3 - rb(nmin1)**3)
*
         else
*
           if(i+nspan.gt.ib) then
c             nmin1 = 2*i - ib - nspan
c             denb(i) = co/(rb(ib)**3 - rb(nmin1)**3)
             nmin1 = ib - 2*nspan
             if(nmin1.le.0) then
               print*,' wrong-nmin0 i,ib,nmin,nspan =',i,ib,nmin,nspan
               call flush(6)
               nmin1 = 1
             endif
             nmax1 = ib
             denb(i) = co/(rb(nmax1)**3 - rb(nmin1)**3)
*
           else
*
             nmin1 = i - nspan
             if(nmin1.le.0) then 
               print*,' wrong-nmin1 i,ib,nmin,nspan =',i,ib,nmin,nspan
               call flush(6)
               nmin1 = 1
             endif
             nmax1 = i + nspan
             denb(i) = co/(rb(nmax1)**3 - rb(nmin1)**3)
*
           endif
*
         endif
*
 20      continue
*
*       maybe is good to smooth the binary density
* 
         return
*
         end
*
*      subroutine energy(ienerg)
*
*
*       total energy.
*       -------------
*
      include 'common.h'
*
      real*8 zzz
      integer i,n,ienerg,im
*
*
*     ----------------------------------------------------------
*
*         ienerg = 1    energies calculated from x and xdot
*         ienerg = 2    energies calculated from r, vr and vt
*
*     ----------------------------------------------------------
*
      n=nt
      zzz = 0.d0
*
*       calculate the potential and kinetic energy
*
*   
      if(ienerg.eq.1) then
*
        call coepot
*
        zkin = 0.0d0
        pot = 0.0d0
*
        do 20 i=1,n
           im = iname(i)
*
*      sum the potential energy
* 
           pot = pot - 0.5d0*body(im)*u(i)
*
*      sum the kinetic energy (include c.m. bodies but not components)
*
           zkin = zkin + body(im)*(xdot(im,1)**2 + xdot(im,2)**2 + 
     &            xdot(im,3)**2)
*
   20   continue
*
        zkin = 0.5d0*zkin
*
        return
*
      endif
*
*
      if(ienerg.eq.2) then
*
        if(time.eq.0.0d0) call coepot
*
        zkin = 0.0d0
        pot = 0.0d0
*
        do 30 i=1,n
           im = iname(i)
           zzz = zzz + body(im)
*
*      sum the potential energy
* 
           pot = pot - 0.5d0*body(im)*u(i)
*
*      sum the kinetic energy (include c.m. bodies but not components)
*
           zkin = zkin + body(im)*(vr(im)**2 + vt(im)**2)
*
   30   continue
*
      endif
*
      zkin = 0.5*zkin
*
        print*,'energy smt,zkin,pot,n,nt =',zzz,zkin,pot,n,nt
*
      return
*
      end
*
*
*       total energy = zkin - pot + etide + ebin + esub + emerge + ecoll.
*
*
*
*
*
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
      if(imodel.eq.3.or.imodel.eq.4) then
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
      open(44,file='binarym.dat',access='append')
      open(43,file='starm.dat',access='append')

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

      subroutine input
*
*
*       parameter input.
*       ----------------
*
      include 'common.h'
*
      integer i,ixx1,ixx2,ixx3
*
      real*8 xxx4
*
*
*              open input file and read initial parameters
*
*
      open(1,file='mont.run')
*
      read(1,*) ixx1,ixx2,istart,ncor,nmin,ixx3,nzonc,nminzo,ntwo,
     &          imodel,iprint,ib3f,iexch,tcrit,tcomp,qe,alphal,
     &          alphah,brakem,body1,bodyn,fracb,amin,amax,qvir,rbar,
     &          xxx4,w0,bmin0,bmax,tau0,gamma,xtid,rplum,dttp,dtte,
     &          dtte0,tcrevo,xtau,ytau,ybmin,zini,ikroupa,iflagns,
     &          iflagbh,nitesc
*
      read(1,*) (flagr(i),i=1,nlagra)
*
      close(1)
      write(6,*)ixx1,ixx2,istart,ncor,nmin,ixx3,nzonc,nminzo,ntwo,
     &          imodel,iprint,ib3f,iexch,tcrit,tcomp,qe,alphal,
     &          alphah,brakem,body1,bodyn,fracb,amin,amax,qvir,rbar,
     &          xxx4,w0,bmin0,bmax,tau0,gamma,xtid,rplum,dttp,dtte,
     &          dtte0,tcrevo,xtau,ytau,ybmin,zini,ikroupa,iflagns,
     &          iflagbh,nitesc
      call flush(6)
*
      bmin = bmin0
      bmax = 2.d0*bmin0 
      if(istart.eq.1) then
        irun = ixx1
        nt = ixx2
        nt0 = nt
        iseed = ixx1
        nz0 = ixx3
        zmbar = xxx4
        nt00 = nt0
      endif
*
*
*     ------------------------------------------------------------------------
*
*                  INPUT PARAMETERS
*                  ----------------
*
*     irun    -  initial sequence of random numbers
*     nt      -  total number of objects (stars and binaries) at T=0
*                ns - number of single stars, nb - number of binaries
*                (nt = ns + nb), nss - number of stars (nss = nt + nb)
*     istart  -  1 - initial model,    .ne.1 - restart
*     ncor    -  number of stars to calculate the central parameters
*     nmin    -  minimum number of stars to calculate the central
*                parameters
*     nz0     -  number of stars in each zone at T=0
*     nzonc   -  minimum number of zones in the core
*     nminzo  -  minimum number of stars in a zone
*     ntwo    -  maximum index of 2
*     imodel  -  initial model: 1- uniform & isotropic, 2- Plummer,
*                3- King, 4 - M67
*     iprint  -  0- full diagnostic information, 1- diagnostic info.
*                suppressed
*     ib3f    -  1 - Spitzer's, 2 - Heggie's formula for three-body binary
*                interaction with field stars, 3 - use Pmax for interaction
*                probability  4 - three- and four-body numerical integration
*     iexch   -  0 - no exchange in any interactions, 1 - exchange only in 
*                binary field star interacions, 2 - exchange in all
*                interactions (binary - field and binary - binary)
*     tcrit   -  termination time in units of the crossing time
*     tcomp   -  maximum computing time in hours
*     qe      -  energy tolerance
*     alphal  -  power-law index for initial mass function for masses
*                smaller than breake mass: -1 - equal mass case
*     alphah  -  power-law index for initial mass function for masses
*                greater than breake mass. If alphal=alphah the IMF does
*                not have a break
*     brakem  -  the mass in which the IMF is broken. If brakem is smaller
*                than the minimum mass (bodyn) than the break mass is as for 
*                the Kroupa mass function (brakem = 0.5 Mo) 
*     body1   -  maximum particle mass before scaling (solar mass)
*     bodyn   -  minimum particle mass before scaling (solar mass)
*     fracb   -  primordial binary fraction by number. nb = fracb*nt,
*                ns = (1 - fracb)*nt, nss = (1 + fracb)*nt
*                fracb > 0 - primordial binaries
*                fracb = 0 - only dynamical binaries
*     amin    -  minimum semi-major axis of binaries (in sollar units)
*                = 0 then amin = 2*(R1+R2), > 0 then amin = amin
*     amax    -  maximum semi-major axis of binaries (in sollar units)
*     qvir    -  virial ratio  (qvir = 0.5 for equilibrium)
*     rbar    -  tidal radius in pc, halfmass radius in pc for isolated
*                cluster. No scaling - rbar = 1
*     zmbar   -  total mass of the cluster in sollar mass, 
*                no scaling zmbar = 1
*     w0      -  king model parameter
*     bmin    -  minimum value of sin(beta^2/2)
*     bmax    -  maximum value of sin(beta^2/2)
*     tau0    -  time step for a complite cluster model
*     gamma   -  parameter in the Coulomb logarithm (standard value = 0.11)
*     xtid    -  coeficient in the front of cluster tidal energy:
*                                                            -xtid*smt/rtid
*     rplum   -  for M67 rtid = rplum*rsplum (rsplum - scale radius for
*                plummer model)
*     dttp    -  time step (Myr) for profile output
*     dtte    -  time step (Myr) for mloss call for all objects
*                greater then tcrevo. For tphys less then tcrevo time step
*                is eqiual to dtte0
*     dtte0   -  time step (Myr) for mloss call for all objects for tphys
*                less then tcrevo. For tphys greater then tcrevo time step
*                is eqiual to dtte
*     tcrevo  -  critical time for which time step for mloss call changes from
*                dtte0 to dtte
*     xtau    -  call mloss for a particlular object when
*                (uptime(im1) - olduptime(im1))/tau/tscale < xtau
*     ytau    -  multiplication of tau0 (tau = ytau*tau0) after time
*                greater than tcrevo
*     ybmin   -  multiplication of bmin0 (bmin = ybmin*bmin0) after time
*                greater than tcrevo
*     zini    -  initial metalicity (solar z = 0.02, globular clusters
*                M4 - z = 0.002,  NGC6397 - z = 0.0002)
*     ikroupa -  0 - the initial binary parameters are picked up
*                according Kroupa's eigenevolution and feeding algorithm
*                (Kroupa 1995, MNRAS 277, 1507)
*                1 - the initial binary parameters are picked as for M67
*                model (Hurley et al. 2005)
*     iflagns -  0 - no SN natal kiks for NS formation, 1 - SN natal kicks
*                only for single NS formation, 2 - SN natal kick for single
*                NS formation and NS formation in binaries
*     iflagbh -  0 - no SN natal kiks for BH formation, 1 - SN natal kicks
*                only for single BH formation, 2 - SN natal kick for single
*                BH formation and BH formation in binaries
*     nitesc  -  0 - no iteration of the tidal radius and induced mass loss
*                due to stellar evolution, 1 - iteration of the tidal radius
*                and induced mass loss due to stellar evolution
*
*     ------------------------------------------------------------------------
*
*

      return
*
      end
*
*
*
*
*
      subroutine inb3b3(nup)
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
      integer nup
*
      if((nbin3-nescb3-ndist3-ndist4-ndiste-nmerge).le.1) return
*
*     if number of binaries presented in the system is less then 100
*     compute probability for binary-binary interaction using
*     Stodolkiewicz procedure (Acta Astronomica 1986), otherwise
*     compute probability for binary-binary interaction using Pmax
*     (Hut and Bahcal ApJ 268, 1983 and Bacon and Sigurdsson 
*     astro-ph96/03036)
*
      if((nbin3-nescb3-ndist3-ndist4-ndiste-nmerge).le.100) then
*
        call intb3b3a(nup)
*
      else
*
        call intb3b3b(nup)
*
      endif
*
      return
*
      end
*
*
*
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
      open(44,file='binarym.dat',access='append')
      open(43,file='starm.dat',access='append')
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
 9321             format(1x'intb3f-evo im1,nob,nsevo,dm3f,ssevo,ehkin=',
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
*
*
      subroutine kick(k,vb)
*
*      compute components of the 'kick' velocity (gained during
*      --------------------------------------------------------
*      interaction between objects) of object and new radial
*      ------------------------------------------------------- 
*      and tangencial velocities
*      -------------------------

      include 'common.h'
*
      real*8 q1,q2,cteta,steta,cosffi,sinffi,vb,b1,c1,dv,zz,xvr,xvt
*
      real*4 ran2
*
      integer im1,k
*
*      used formulae from Stodolkiewicz Acta Astronomica 1986, 36, p19
*
      q1 =ran2(irun) - 0.5d0
      cteta = -2.0d0*q1
      steta = sqrt(1.0d0 - cteta**2)
      im1 = iname(k)
      xvr = vr(im1)
      xvt = vt(im1)
*
      q2 = 2.0d0*pi*ran2(irun)
      sinffi = sin(q2)
      cosffi = cos(q2)
      b1 = 2.0d0*(vr(im1)*cteta + vt(im1)*steta*cosffi)
      c1 = -2.0d0*vb/body(im1)
      dv = 0.5d0*(-b1 + sqrt(b1*b1 - 4.0d0*c1))
      vr(im1) = vr(im1) + dv*cteta
      b1 = vt(im1) + dv*steta*cosffi
      c1 = dv*steta*sinffi
      vt(im1) = sqrt(b1*b1 + c1*c1)
*
      zz = vt(im1)**2 + vr(im1)**2 - xvr**2 - xvt**2
      zz = 0.5d0*body(im1)*zz - vb
*
c      if(vb.gt.1.0d-4) then
c        write(6,*) 'k,im1,r,xvr,xvt,vr,vt,zz =',k,im1,r(k),xvr,xvt,
c     &              vr(im1),vt(im1),zz
c      endif
*
      return
      end
*
*
 
       
      subroutine king
*
*
*       king model as a initial model of a system
*       -----------------------------------------
*
*
      include 'common.h'
*
      integer n,npo,nsp,nxy,i,iqw,k,ispr,nok,nbad,ixter,iro,ilo
*
      real*4 ran2
*
      real*8 xa,ya,y2a,row0,z1,z2,z3,wlgr,wlw,ystart,s1,s2,g,
     &       dw,ss1,ss2,aa,b,ssd,b2,yy,yp1,xx0,z5,x2,xend,
     &       dx,y1old,y2old,x1,x3,x4,zz,row,z6,z66,z7,z8,z4,dr,r1,
     &       smass,beta,rr0,rrr0,rr1,sma,bet,yk,ri,wr,fmax,func1,
     &       func3,func4
*
*
      common /spl/ xa(100000),ya(100000),y2a(100000),row0,nsp
      common /fu3/ z1(100000),z3(100000),wlgr(100000),npo
      common /fu4/ z2(100000),wlw(100000)
*
      dimension ystart(2),s1(100000),s2(100000),g(8)
*
      external func3,func1,func4
*
*
      n = nt
      dw=0.005d0
*
      nxy=sqrt(w0+1.0d0)/dw
      ss1=0.0d0
      ss2=0.0d0
      xa(1)=0.0d0
      ya(1)=0.0d0
*
      do 10 i=1,nxy
      aa=(i-1)*dw
      b=i*dw
*
      call qromb(func1,aa,b,ssd)
      ss1=ss1+ssd
      b2=b*b
      xa(i+1)=b2
      ya(i+1)=ss1*exp(b2)
*
 10   continue
*
      nsp=nxy+1
      yy=xa(nsp)
      yp1=ya(nsp)+0.5d0*yy**1.5
      call spline(xa,ya,nsp,1.0d30,yp1,y2a)
*
      call splint(xa,ya,y2a,nsp,w0,row0)
*
      xx0=0.01
      ystart(1)=w0-1.5d0*xx0**2+27.0d0*xx0**4/40.0d0/w0
      ystart(2)=-3.0d0*xx0+2.7d0*xx0**3/w0
      z5=ystart(2)
*
      x2=log10(xx0)
      xend=5.0d0
      dx=0.02d0
      z1(1)=xx0
      z2(1)=ystart(1)
      iqw=0
      yy=0.0d0
*
      do 20 k=2,1000
*
      y1old=ystart(1)
      y2old=ystart(2)
 30   x1=x2
      x2=x1+dx
      x3=10.0d0**x1
      x4=10.0d0**x2
      ispr=0
      call odeint(ystart,2,x3,x4,1.d-8,1.d-5,0.0d0,nok,nbad)
*
      if(ystart(1).lt.0.0d0) then
        if(abs(y1old-yy).lt.1.d-11.and.ispr.eq.1) go to 40
        if(y1old.lt.1.d-8) go to 40
        dx=dx/2.0d0
        x2=x1
        ystart(1)=y1old
        ystart(2)=y2old
        yy=y1old
        go to 30
      endif
*
      z1(k)=x4
      z2(k)=ystart(1)
      ispr=1
*
      if(x2.gt.xend) go to 40
*
      iqw=iqw+1
*
 20   continue
*
 40   continue
*
      npo=iqw+1
      do 50 k=1,npo
      zz=z2(k)
      call splint(xa,ya,y2a,nsp,zz,row)
      if(k.eq.1) z6=row
      if(k.eq.npo) z66=row
      z3(k)=row/row0
*      write(1,1000) log10(z1(k)),z2(k),log10(z3(k))
 50   continue
*
* 1000 format(1x,1p3e12.4)
*
      z7=z5*(z6+0.5d0*z2(1)**1.5)/row0
      z8=z4*(z66+0.5d0*z2(npo)**1.5)/row0
      call spline(z1,z3,npo,z7,z8,wlgr)
      call spline(z1,z2,npo,z5,z4,wlw)
      rtid=log10(z1(npo))
      dr=0.01d0
      r1=-2.0d0
      iro=1
      smass=0.0d0
      beta=0.0d0
      ixter=0
      s1(1)=0.0d0
      s2(1)=0.0d0
*
      do 60 k=1,1000
      rr0=r1
      r1=rr0+dr
      if(r1.gt.rtid) then
        r1=rtid
        ixter=1
      endif
      rrr0=10.0d0**rr0
      rr1=10.0d0**r1
      if(k.eq.1) then
        rrr0=1.0d-5
        rr1=10.0d0**rr0
        r1=rr0
      endif
*
      call qromb(func3,rrr0,rr1,sma)
      call qromb(func4,rrr0,rr1,bet)
      smass=smass+sma
      beta=beta+bet
      iro=iro+1
      s1(iro)=rr1
      s2(iro)=smass
*      write(2,1000) float(iro),rr1,smass
      if(ixter.eq.1) go to 70
 60   continue
*
 70   continue
*
      yk=beta*rr1/18.0d0/smass**2
      rtidkg=2.0d0*yk+1.0d0
*
      rtid = rtidkg
*
      do 80 k=1,iro
      s2(k)=s2(k)/s2(iro)
 80   continue
*
*
      do 90 i = 1,n
 100      g(1) = ran2(irun)
          if (g(1).lt.1.0d-10) go to 100
          call hunt(s2,iro,g(1),ilo)
          if((ilo.eq.0).or.(ilo.eq.iro)) go to 100
          ri=(g(1)-s2(ilo))*(s1(ilo+1)-s1(ilo))
          ri=ri/(s2(ilo+1)-s2(ilo))+s1(ilo)
          g(2) = ran2(irun)
          g(3) = ran2(irun)
          x(i,3) = (1.0d0 - 2.0d0*g(2))*ri
          x(i,1) = sqrt(ri**2 - x(i,3)**2)*cos(twopi*g(3))
          x(i,2) = sqrt(ri**2 - x(i,3)**2)*sin(twopi*g(3))
          call splint(z1,z2,wlw,npo,ri,wr)
          fmax=exp(wr-1d0)/wr
          if(wr.lt.2.0d0) then
		fmax=0.55d0*(exp(0.5d0*wr)-1.0d0)
	  endif
 110      g(4) = ran2(irun)
          g(5) = ran2(irun)*fmax
          g(6) = g(4)**2*(exp(-wr*(g(4)**2-1.0d0))-1.0d0)
          if (g(5).gt.g(6)) go to 110
          g(8)=g(4)*sqrt(2.0d0*wr)
          g(6) = ran2(irun)
          g(7) = ran2(irun)
          xdot(i,3) = (1.0d0 - 2.0d0*g(6))*g(8)
          xdot(i,1) = sqrt(g(8)**2 - xdot(i,3)**2)*cos(twopi*g(7))
          xdot(i,2) = sqrt(g(8)**2 - xdot(i,3)**2)*sin(twopi*g(7))
*
*      write(3,1010) (x(i,k),k=1,3),(xdot(i,k),k=1,3)
* 1010 format(1x,1p6e12.4)
   90 continue

      return
      end
*
*
*
*
      double precision function func3(x)
*
      integer npo
*
      real*8 z1,z3,wlgr,x,yy
*
      common /fu3/ z1(100000),z3(100000),wlgr(100000),npo
*
*
      call splint(z1,z3,wlgr,npo,x,yy)
      func3=yy*x*x
*
      return
      end
*
*
*
*
      double precision function func4(x)   
*
      integer npo
*
      real*8 z1,z2,z3,wlgr,wlw,x,yy,ww
*
      common /fu3/ z1(100000),z3(100000),wlgr(100000),npo
      common /fu4/ z2(100000),wlw(100000) 
*
*
      call splint(z1,z3,wlgr,npo,x,yy)
      call splint(z1,z2,wlw,npo,x,ww)
      func4=ww*yy*x*x
*
      return
      end
*
*
*
*
      double precision function func1(x)
*
      real*8 x
*
*
      func1=exp(-x*x)*x**4
*
      return
      end
*
*
*
      subroutine derivs(x,y,dy)
*
      integer nsp
*
      real*8 xa,ya,y2a,row0,y,dy,x,ro
*
      common /spl/ xa(100000),ya(100000),y2a(100000),row0,nsp
      dimension y(2),dy(2)
*
*
      dy(1)=y(2)
      call splint(xa,ya,y2a,nsp,y(1),ro)
      dy(2)=-2.0d0/x*y(2)-9.0d0*ro/row0
*
      return
      end
*
*
*
*
      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad)
*
      integer nvar,nok,nbad,kmax,kount,i,maxstp,nmax1,nstp
*
      real*8 ystart,x1,x2,eps,h1,hmin,two,zero,tiny,dxsav,xp,yp,
     &       yscal,y,dydx,x,h,xsav,hdid,hnext
*
      common /path/ kmax,kount,dxsav,xp(200),yp(10,200)
* 
      parameter (maxstp=10000,nmax1=10,two=2.d0,zero=0.d0,tiny=1.d-30)
*
      dimension ystart(nvar),yscal(nmax1),y(nmax1),dydx(nmax1)
*
*
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      xsav=x-dxsav*two
      do 16 nstp=1,maxstp
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.zero) h=x2-x
        call rkqc(y,dydx,nvar,x,h,eps,yscal,hdid,hnext)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.zero)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause 'stepsize smaller than minimum.'
           h=hnext
16    continue
      pause 'too many steps.'
*
      return
      end
*
*
*
      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext)
*
      integer nmax1,n,i
*
      real*8 y,dydx,x,htry,eps,yscal,hdid,hnext,ytemp,ysav,
     &       dysav,pgrow,pshrnk,xsav,h,hh,errmax,fcor,one,safety,
     &       errcon

      parameter (nmax1=10,fcor=.0666666667d0,
     &          one=1.0d0,safety=0.9d0,errcon=6.d-4)
*
      dimension y(n),dydx(n),yscal(n),ytemp(nmax1),ysav(nmax1),
     &          dysav(nmax1)
*
*
      pgrow=-0.20d0
      pshrnk=-0.25d0
      xsav=x
      do 11 i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
11    continue
      h=htry
1     hh=0.5d0*h
      call rk4(ysav,dysav,n,xsav,hh,ytemp)
      x=xsav+hh
      call derivs(x,ytemp,dydx)
      call rk4(ytemp,dydx,n,x,hh,y)
      x=xsav+h
      if(x.eq.xsav)pause 'stepsize not significant in rkqc.'
      call rk4(ysav,dysav,n,xsav,h,ytemp)
      errmax=0.0d0
      do 12 i=1,n
        ytemp(i)=y(i)-ytemp(i)
        errmax=max(errmax,abs(ytemp(i)/yscal(i)))
12    continue
      errmax=errmax/eps
      if(errmax.gt.one) then
        h=safety*h*(errmax**pshrnk)
        goto 1
      else
        hdid=h
        if(errmax.gt.errcon)then
          hnext=safety*h*(errmax**pgrow)
        else
          hnext=4.0d0*h
        endif
      endif
      do 13 i=1,n
        y(i)=y(i)+ytemp(i)*fcor
13    continue
*
      return
      end
*
*
*
*
      subroutine rk4(y,dydx,n,x,h,yout)
*
      integer n,nmax1,i
      parameter (nmax1=10)
*
*
      real*8 y,dydx,x,h,yout,yt,dyt,dym,hh,h6,xh
*
      dimension y(n),dydx(n),yout(n),yt(nmax1),dyt(nmax1),dym(nmax1)
*
*
      hh=h*0.5d0
      h6=h/6.0d0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0*dym(i))
14    continue
*
      return
      end
*
*
*
*
      subroutine splint(xa,ya,y2a,n,x,y)
*
      integer n,klo,khi,k
*
      real*8 xa,ya,y2a,x,y,h,a,b
*
      dimension xa(n),ya(n),y2a(n)
*
*
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.0d0) pause 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
*
      return
      end
*
*
*
*
      subroutine qromb(func,a,b,ss)
*
      integer jmax,jmaxp,k,km,j,l
*
      real*8 func,a,b,ss,s,h,dss,eps
*
      parameter(eps=1.0d-6,jmax=20,jmaxp=jmax+1,k=5,km=4)
*
      dimension s(jmaxp),h(jmaxp) 
*
      external func
*
*
      h(1)=1.0d0
      do 11 j=1,jmax
        call trapzd(func,a,b,s(j),j)
        if (j.ge.k) then
	  l=j-km
          call polint(h(l),s(l),k,0.0d0,ss,dss)
          if (abs(dss).lt.eps*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
*
      pause 'too many steps.'
*
      end
*
*
*
      subroutine polint(xa,ya,n,x,y,dy)
*
      integer n,nmax1,ns,i,m
      parameter (nmax1=10) 
*
*
      real*8 xa,ya,x,y,dy,dif,dift,c,d,ho,hp,w,den
*
      dimension xa(n),ya(n),c(nmax1),d(nmax1)
*
*
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n 
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.0d0)pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
*
      return
      end
*
*
*
*
      subroutine trapzd(func,a,b,s,n)
*
      save
*      
      integer n,it,j
*
      real*8 func,a,b,s,tnm,del,x,sum
*
*
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
        it=1
      else
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.0d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
        it=2*it
      endif
*
      return
      end
*
*
*
      subroutine spline(x,y,n,yp1,ypn,y2)
*
      integer nmax2,n,i,k
      parameter (nmax2=1000)
*
*
      real*8 x,y,yp1,ypn,y2,u,sig,p,qn,un
*
      dimension x(n),y(n),y2(n),u(nmax2)
*
*
      if (yp1.gt.0.99d30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0d0
        y2(i)=(sig-1.0d0)/p
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt.0.99d30) then
        qn=0.9d0
        un=0.0d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
*
      return
      end
*
*
*
      subroutine hunt(xx,n,x,jlo)
*
      integer n,jlo,jhi,inc,jm
*
      real*8 xx,x
*
      dimension xx(n)
*
      logical ascnd
*
*
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
*
      end
*
*
*
      subroutine lagrad
* 
*
*       lagrangian radii & mean radial and tangential velocities in the 
*       ---------------------------------------------------------------
*       lagrangian shells
*       -----------------
*
*
      include 'common.h'
*
      real*8 zm,zm1,zm2,vvr,vvt,zmi,vr2,vt2,r21,r22,rla
*
      integer i,j,n,im
*
*
*       determine the lagrangian radii and mean radial and tangential
*       velocities in the lagrangian shells
*
      n = nt
*
      im = 0
      zm = 0.0d0
      zm2 = 0.0d0
*
      do 10 j=1,nlagra
         vvr = 0.0d0
         vvt = 0.0d0
         zmi = flagr(j)*smt
   20    im = im + 1
         i = iname(im)
         zm = zm + body(i)
         vr2 = vr(i) * vr(i)*body(i)
         vt2 = vt(i) * vt(i)*body(i)
         vvr = vvr + vr2
         vvt = vvt + vt2
*
         if(im.eq.n) then
           vr2 = 0.0d0
           vt2 = 0.0d0
           go to 30
         endif
*
         if (zm.le.zmi) go to 20
   30    zm1 = zm - body(i)
         vvr = vvr - vr2
         vvt = vvt - vt2
         if(im.eq.1) then
           r21 = 0.d0
         else
           r21 = r(im - 1)**3
         endif
         v2rl(j) = vvr/(zm1 - zm2)
         v2tl(j) = vvt/(zm1 - zm2)
         ani(j) = 2.0d0 - v2tl(j)/v2rl(j)
         r22 = r(im)**3
         rla = (zmi - zm) * (r22 - r21)/(zm - zm1) + r22
         rlag(j) = rla**one3
         im = im - 1
         zm = zm1
         zm2  =  zm1
*
 10   continue
*
*
      return
*
      end
*
*
*
*
*
      subroutine mloss(nup,lmin,lmax,tphys,ekick,ekickbs,ekickbd)
*
*
*       Mass loss from evolving stars.
*       ------------------------------
*     In this version, stellar and binary evolution are carried out via the
*     interface of McScatter. Modified 26/4/6 so that mass is lost only if
*     the update time is passed.
*
*     mloss is the first subroutine called in relaxall (escape is the last
*     subroutine) so ikind cannot be negative
*
*
      include 'common.h'
*
*
      real*8  timevo,sm1,sm2,sm1s,sm2s,sm1o,sm2o,a,anb,ssevo,ssevot,
     &        cor,desc,escr(nmax),escb(nmax),ebb,abb,pctonbu,ehpot,
     &        ehkin,r1,r2,xlum1,xlum2,tlog,smlog,smsss,sturn,z,
     &        zpars(20),ecc,sm1b,sm2b,sm1so,sm2so,slum1o,slum2o,xxm1,
     &        xxm2,xxl1,xxl2,potold,tphys,dtt,xdt,ssme,tphys1,tphys2,
     &        tphysf1,tphysf2,epoch1,epoch2,tmstim1,tmstim2,smass1,
     &        smass2,s0mass1,s0mass2,scmass1,scmass2,semass1,semass2,
     &        rad1,rad2,radc1,radc2,spin1,spin2,rl1,rl2,tb0,semi0,ecc0,
     &        rade1,rade2,vs(3),xxn,vscal,ekick,ekick1,ekick2,vst2,vst,
     &        vrold,vtold,vss,aursun,yeardy,cfr,dttot,potnew,ss1,ss2,
     &        t1,t2,upt,ddd,ekickbs,ekickbd,smbb,vrold2,vtold2
*
      integer nup,n,i,k,lmin,lmax,im1,ibx,ik1,ik2,kk,ik0,ipoi,id3,id4,
     &        ievo,kevo,idb,ids,id1,id2,ikio,ikin,iexist,iblue,
     &        ik0n,inum,inumt,iob,ik1n,ik2n,ndist,iime,ikick,ixkick,
     &        ikickbs,ixkickbs,ikickbd,ixkickbd
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
      kevo = 0
      desc = 0.d0
      cor = 0.d0
      n = nt
      ssevot = 0.d0
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
      dttot = tau*tscale0*log(gamma*nt00)/log(gamma*nt)*
     &        float(nt)/float(nt00)
      tlog = log10(1.d6*timevo)
*
      print*,'mloss  lmin,lmax,nup,timevo = ',lmin,lmax,nup,timevo
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
      call energy(2)
      potold = pot
*
      open(44,file='binarym.dat',access='append')
      open(43,file='starm.dat',access='append')          
*
      do 10 i = lmin,lmax
*
         ixkick = 0
         ixkickbs = 0
         ixkickbd = 0
*
         vs(1) = 0.d0
         vs(2) = 0.d0
         vs(3) = 0.d0
*
         if(i.gt.nt) go to 10
*
           ievo = 0
           ssevo = 0.d0
           im1 = iname(i)
*
*      guard against massless particles
*
           if((abs(body(im1)).lt.tolm).or.(r(i).gt.1.d8)) go to 10
*
           if(ikind(im1).lt.0) then
             ikind(im1) = -ikind(im1)
           endif
*
           dtt = uptime(im1) - oldtime(im1)
           xdt = dtt/dttot
*
           if(ikind(im1).eq.2) then
*
*               take care for ojects which are binaries
*
             if(uptime(im1).lt.timevo.or.xdt.le.xtau.or.
     &          nup.eq.-100) then
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
               ddd = (sm1b + sm2b)*zmbar - (sm1so + sm2so)
               write(6,9123) im1,ibx,sm1b*zmbar,sm1so,sm2b*zmbar,sm2so,
     &                       ddd
 9123          format(1x,'im1,ibx,sm1b,sm1so,sm2b,sm2so,dmassdm-1 =',
     &                2i9,1p5e16.8)
               if(xxm1.gt.1.0d-5.or.xxm2.gt.1.d-5)
     &         print*,'wrong-b-sing im1,idb,id1,id2,sm1,',
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
     &              (smass1 + smass2)))*yeardy
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
     &              (ik2n.ge.13.and.ik2n.lt.15.and.ik2.lt.13)) then
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
 8712       format(1x,'timevo-all-bs,im1,ik1n,ik2n,ikickbs,ikicktbs,',
     &       'ekick1,ekick2,ekickbs,ekicktbs,vscale,vs(1),vs(2),vs(3),',
     &       'vr(im1),vt(im1),iflagns,iflagbh = ',1pe12.4,5i8,1p10e12.4,
     &        2i4)     
                 endif
*
 140             continue
*
                 sm1o = sm1
*
                 sm1 = sm1s/zmbar
*
                 nmloev = nmloev + 1
                 slosev = slosev + (sm1o - sm1)
                 ssevo = sm1o - sm1
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
                   kevo = kevo + 1
                   escr(kevo) = r(i)
                   escb(kevo) = ssevo
      print*,'evo-bin-m: kevo,im1,idb,id1,id2,ibx,nt0,nameb,',
     &'ik0,ik0n,ik1,ik1n,ik2n,sm1so,sm1s,sm2so,sm2s,sevo,uptim= ',
     &kevo,im1,idb,id1,id2,ibx,nt0,nameb(im1),ik0,ik0n,ik1,ik1n,ik2,
     &ik2n,sm1so,sm1s,sm2so,sm2s,ssevo*zmbar,uptime(im1)
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
                   print*,'iob,im1,sm1,sm2=',iob,im1,sm1,sm2
                 endif
*
*     print only when binary component masses and luminosities change
*     because of binary evolution by more than 1%
*
                 ddd = (body(im1) + ssevo)*zmbar - (sm1so + sm2so)
                 write(6,9124) im1,ibx,body(im1)*zmbar,ssevo*zmbar,sm1so
     &                         ,sm2so,ddd
 9124            format(1x,'im1,ibx,smb,ssevo,sm1so,sm2so,dmassdm-2 =',
     &                2i9,1p5e16.8)
*                          
                 xxm1 = abs((sm1s - sm1so)/sm1so)
                 xxm2 = abs((sm2s - sm2so)/sm2so)
                 xxl1 = abs((xlum1 - slum1o)/slum1o)
                 xxl2 = abs((xlum2 - slum2o)/slum2o)
c                 if(xxm1.gt.0.01d0.or.xxm2.gt.0.01d0.or.xxl1.gt.0.01d0
c     &              .or.xxl2.gt.0.01d0) then
*
                   write(44,420) idb,id1,id2,im1,timevo,sm1*zmbar,
     &                       sm1o*zmbar,sm2*zmbar,sm2o*zmbar,r1,r2,
     &                       xlum1,xlum2,a,ecc,r(i),ik0,ik0n,ik1,ik1n,
     &                       ik2,ik2n,iexist,ibstra(im1),-1
*
 420               format(1x,'#',4i7,1p12e16.8,9i4)
*
c                 endif
*
                 go to 120
*
               endif
*
*      first take care about binary mergers
*
               if((iexist.eq.0.and.ik0n.eq.4).or.
     &            (iexist.eq.0.and.ik0n.eq.5)) then
*
      write(6,4390) im1,idb,id1,id2,ik0,ik1,ik2,sm1*zmbar,sm2*zmbar,
     & semi0,ecc0,rad1,rad2,rl1,rl2,t1,t2,epoch1,epoch2,ik0n,ik1n,
     & ik2n,sm1s,sm2s
 4390 format(1x,'mevo-merged-init im1,idb,id1,id2,ik0,ik1,ik2,sm1,sm2,',
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
                   open(27,file='bhmergers.dat',access='append')
                   write(27,1000) im1,ibx,idb,id1,id2,ikind(im1),
     &  timevo,z,(zpars(kk),kk=1,20),vss,ecc0,semi0,tb0,ik0,ik1,tphys1,
     &  tphysf1,epoch1,tmstim1,smass1,s0mass1,scmass1,semass1,rad1,
     &  radc1,rade1,slum1o,spin1,rl1,ik2,tphys2,tphysf2,epoch2,tmstim2,
     &  smass2,s0mass2,scmass2,semass2,rad2,radc2,rade2,slum2o,spin2,rl2
 1000 format(1x,6i8,1p26e14.6,2i4,1p14e14.4,i4,1p14e14.4)
                   close(27)
                  endif 
*
                 if(ik2.eq.14.and.(ik1.ge.10.and.ik1.le.13)) then
                   open(27,file='bhmergers.dat',access='append')
                   write(27,1000) im1,ibx,idb,id1,id2,ikind(im1),
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
                   open(28,file='bhmergers-all.dat',access='append')
                   write(28,1000) im1,ibx,idb,id1,id2,ikind(im1),
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
                   open(29,file='bhmergers-all-all.dat',access='append')
                   write(29,1010) im1,ibx,idb,id1,id2,ikind(im1),   
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
 4321              format(1x,'mipo1=2, ik1n,ik2n,sm1s,sm2s,im1,id2 = ',
     &2i6,1p2e14.6,2i9)
                 else
                   names(im1) = id1
                   ids = id1
                   ipoi = 1
                   sm2s = 0.d0
                   write(6,4322) ik1n,ik2n,sm1s,sm2s,im1,id2
 4322              format(1x,'mipo1=1, ik1n,ik2n,sm1s,sm2s,im1,id2 = ',
     &2i6,1p2e14.6,2i9)
                 endif
*
                 nmerge = nmerge + 1
                 iime = iime + 1
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
                 slosev = slosev + (sm1o - sm1)
                 ssevo = sm1o - sm1
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
      print*,'ddd-m timevo,nmerge,iime,ssme=',timevo,nmerge,iime,cfr
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
      print*,'b.str.merg1m.id1,id2,ik1,ik2,ikn1,ikn2,ibsm,sm1s,sturn=',
     &               id1,id2,ik1,ik2,ik1n,ik2n,ibsm,sm1s,sturn
                   endif
                   if(ik2n.le.1.and.iblue.eq.0) then
                     if(sm2s.gt.sturn.and.ibstra(im1).eq.0) then
                       ibstra(im1) = 1
                       ibsm = ibsm + 1
                     endif
      print*,'b.str.merg2m.id1,id2,ik1,ik2,ikn1,ikn2,ibsm,sm2s,sturn=',
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
                   kevo = kevo + 1
                   escr(kevo) = r(i)
                   escb(kevo) = ssevo
      print*,'evo-mergedm: kevo,im1,idb,id1,id2,ibx,nt0,names,nameb,',
     &'ik0,ik0n,ik1,ik1n,ik2,ik2n,sm1so,sm1s,sm2so,sm2s,ss1,ss2,ssevo,',
     &'uptimeb,uptimes,nmerge = ',kevo,im1,idb,id1,id2,ibx,nt0,
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
*      from the system
*
                   r(i) = 1.0d+15 + im1
                   xmax(im1) = 1.01d0*r(i)
                   print*,'merge-m  iob,im1,sm1,sm2=',iob,im1,sm1,sm2
                 endif
*
                 ddd = (body(im1) + ssevo)*zmbar - (sm1so + sm2so)
                 write(6,9125) im1,ibx,body(im1)*zmbar,ssevo*zmbar,sm1so
     &                         ,sm2so,ddd
 9125            format(1x,'im1,ibx,smb,ssevo,sm1so,sm2so,dmassdm-3 =',
     &                2i9,1p5e16.8)
*
                 write(44,420) idb,id1,id2,im1,timevo,ss1,
     &                         sm1o*zmbar,ss2,sm2o*zmbar,r1,r2,
     &                         xlum1,xlum2,a,ecc,r(i),ik0,ik0n,ik1,ik1n,
     &                         ik2,ik2n,iexist,ibstra(im1),-1
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
                 print*,'idb-m,id1,id2,id3,id4=',idb,id1,id2,id3,id4
*
*     check for supernove formation and NS/BH kick
*
                 sm1o = sm1
                 sm2o = sm2
                 sm1 = sm1s/zmbar
                 sm2 = sm2s/zmbar
*          
                 if((ik1n.ge.13.and.ik1n.lt.15.and.ik1.lt.13).or.
     &             (ik2n.ge.13.and.ik2n.lt.15.and.ik2.lt.13)) then
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
            write(6,8719) timevo,im1,ik1n,ik2n,ikickbd,ikicktbd,ekick1,          ,
     &                    ekick2,wkick12,ekick22,ekickbd,ekicktbd,
     &                    vscal,vs(1),vs(2),vs(3),vr(im1),vt(im1),
     &                    vs2(1),vs2(2),vs2(3),vr(nt0),vt(nt0),nt0,
     &                    iflagns,iflagbh,ekicktb2
 8719       format(1x,'timevo-all-bd,im1,ik1n,ik2n,ikickbd,ikicktbd,',
     &       'ekick1,ekick2,ekick12,ekick22,ekickbd,ekicktbd,vscale,',
     &       'vs(1),vs(2),vs(3),vr(im1),vt(im1),vs2(1),vs2(2),vs2(3),',
     &       'vr(nt0),vt(nt0),nt0,iflagns,iflagbh,ekicktb2 =',1pe12.4,
     &       5i8,1p17e12.4,3i8,1pe12.4)     
                 endif
*
 145             continue
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
                 slosev = slosev + (sm1o - sm1)
                 ssevo = sm1o - sm1
*
                 nmloev = nmloev + 1
                 slosev = slosev + (sm2o - sm2)
                 ssevo = ssevo + (sm2o - sm2)
*
                 call get_ss_updatetime(id1,uptime(im1))
                 call get_ss_updatetime(id2,uptime(nt0))
*                                  
                 if(ssevo.gt.1.0d-10) then
                   kevo = kevo + 1
                   escr(kevo) = r(i)
                   escb(kevo) = ssevo
      print*,'evo-disrum: kevo,im1,idb,id1,id2,ibx,nt0,names1-2,nameb,',
     &'ik0,ik0n,ik1,ik1n,ik2n,sm1so,sm1s,sm2so,sm2s,ss1,ss2,ssevo,',
     &'uptimes1,uptimes2,ndist,id3,id4= ',kevo,im1,idb,id1,id2,ibx,nt0,
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
      print*,'m :im1,nt0,idb,id1,id2,names(im1),names(nt0),nameb(im1),',
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
                 binin(inumt,4) = sm2s
                 binin(inumt,5) = timevo
                 binin(inumt,6) = 1.d6
                 binin(inumt,7) = 6

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
                   print*,'dissol m: iob,im1,sm1,sm2=',iob,im1,sm1,sm2
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
                 ddd = (body(im1) + body(nt0) + ssevo)*zmbar - 
     &                 (sm1so + sm2so)
                 write(6,9126) im1,ibx,(body(im1)+body(nt0))*zmbar,
     &                         ssevo*zmbar,sm1so,sm2so,ddd
 9126            format(1x,'im1,ibx,smb,ssevo,sm1so,sm2so,dmassdm-4 =',
     &                2i9,1p5e16.8)
*
                 write(44,420) idb,id1,id2,im1,timevo,sm1*zmbar,
     &                         sm1o*zmbar,sm2*zmbar,sm2o*zmbar,r1,r2,
     &                         xlum1,xlum2,a,ecc,r(i),ik0,ik0n,ik1,ik1n,
     &                         ik2,ik2n,iexist,ibstra(im1),-1
*
                 go to 120
*
               endif
*
             endif
*
           else if(ikind(im1).eq.1.or.ikind(im1).ge.3) then
*
             if(uptime(im1).lt.timevo.or.xdt.le.xtau.or.
     &         nup.eq.-100) then
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
                 print*,'wrong-ssm im1,ids,idb,sm1o,sm1b,ikind=',
     &                   im1,ids,idb,sm1o*zmbar,sm1b*zmbar,ikind(im1)
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
               if((ikin.ge.13.and.ikin.lt.15).and.(ikio.lt.13)) then
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
                 write(6,8717) timevo,im1,ikin,ikick,ikickt,ekick1,
     &                         ekick2,ekick,ekickt,vscal,vs(1),vs(2),
     &                         vs(3),vr(im1),vt(im1),iflagns,iflagbh
 8717            format(1x,'timevo-mloss,im1,ikin,ikick,ikickt,ekick1',
     &           'ekick2,ekick,ekickt,vscale,vs(1),vs(2),vs(3),vr(im1)',
     &           'vt(im1),iflagns,iflagbh = ',1pe12.4,4i8,1p10e12.4,2i4)     
               endif
*
 150           continue               
*
*       check if a blue straggler leaves the main sequence
*
               if(ibstra(im1).gt.0.and.ibstra(im1).le.4.and.
     &            ikin.gt.1) ibstra(im1) = 0
*
               nmloev = nmloev + 1
               slosev = slosev + (sm1o - sm1)
               ssevo = sm1o - sm1
*
               call get_ss_updatetime(ids,uptime(im1))
*
               if(ssevo.gt.1.0d-10) then
                 kevo = kevo + 1
                 escr(kevo) = r(i)
                 escb(kevo) = ssevo
      print*,'single: kevo,im1,ids,idb,ikio,ikin,sm1so,sm1s,ssevo,ut =',
     &kevo,im1,ids,idb,ikio,ikin,sm1so,sm1s,ssevo*zmbar,uptime(im1)
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
                  print*,'iob-m,im1,sm1=',iob,im1,sm1
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
*
               ddd = (body(im1) + ssevo)*zmbar - sm1o*zmbar
               write(6,9127) im1,ids,body(im1)*zmbar,ssevo*zmbar,
     &                       sm1o*zmbar,sm1so,ddd
 9127            format(1x,'im1,ibx,sms,ssevo,sm1so,dmassdm-5 =',
     &                2i9,1p5e16.8)
*
               xxm1 = (sm1s - sm1so)/sm1so
               xxl1 = (xlum1 - slum1o)/slum1o
c               if(xxm1.gt.0.01d0.or.xxl1.gt.0.01d0) then
c                 print*,'000', ids,im1,timevo,sm1*zmbar,sm1o*zmbar,r1,
c     &                       xlum1,r(i),ikio,ikin,ibstra(im1)
                 call flush(6)
                 write(43,440) ids,im1,timevo,sm1*zmbar,sm1o*zmbar,r1,
     &                     xlum1,r(i),ikio,ikin,ibstra(im1),ikind(im1)
*
 440             format(1x,'@',2i7,1p6e16.8,4i4)
*
               go to 120
*
             endif
*
*      Calculate masslose and heating due to stelar evolution
*      For neutron stars assigne proper velocity - to be done
*      Update the potential
*
           else
*
              write (6,*) 'ikind != 1, 2,3  or 4 in mloss1'
              stop
*
           endif
*
 120       continue
*
           if(ievo.eq.1) then
*
              ssevot = ssevot + ssevo
              if(imodel.ge.3) then
c                 cor = -0.5d0*ssevo**2/r(i)
c                 ehpot = -ssevo*(-smt/rtid - u(i))
c                 ehpot = -ssevo*u(i)
                 if(ixkick.eq.1.or.ixkickbs.eq.1.or.ixkickbd.eq.1) then
                   ehkin = 0.5d0*ssevo*(vrold**2 + vtold**2)
                 else  
                   ehkin = 0.5d0*ssevo*(vr(im1)**2 + vt(im1)**2)
                 endif
                 ehmlev = ehmlev - ehkin
c                 ehmlev = ehmlev + ehpot - ehkin + cor
              else
c                 cor = -0.5d0*ssevo**2/r(i)
c                 ehpot = -ssevo*u(i)
                 if(ixkick.eq.1.or.ixkickbs.eq.1.or.ixkickbd.eq.1) then
                   ehkin = 0.5d0*ssevo*(vrold**2 + vtold**2)
                 else
                   ehkin = 0.5d0*ssevo*(vr(im1)**2 + vt(im1)**2)
                 endif
                 ehmlev = ehmlev - ehkin
c                 ehmlev = ehmlev + ehpot - ehkin + cor
              endif
           endif
*
 10   continue
*
      close(43)
      close(44) 
*
*      update total mass, tidal radius and potential
*
c      do 130 k = 1,kevo
c         if(kevo.eq.1) go to 130
c         do 140 kk = k+1,kevo
c            desc = desc - escb(k)*escb(kk)/escr(kk)
c 140     continue
c 130  continue
*
c      ehmlev = ehmlev + desc
      smt = smt - ssevot
      if(imodel.ge.3) rtid = rtidkg*smt**0.3333333
      if(ndist.eq.0.and.iob.eq.0) then
        call coepot
        call energy(2)
        potnew = pot
        enepot = enepot + potnew - potold
        print*,'potold,potnew,n,nt,ndist,nt0,nup=',potold,potnew,n,nt,
     &  ndist,nt0,nup
      endif                       
*
      print*,'nup,kevo,desc,ehmlev =',nup,kevo,desc,ehmlev
      print*,'fff timevo,iime,ssme =',timevo,iime,ssme
      print*,'ixkick,ixkickbs,ixkickbd=',ikick,ikickbs,ikickbd
      print*,'ixkick,ixkickbs,ixkickbd=',ekick,ekickbs,ekickbd
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
      call sort2(nt0,r,iname)
      nzst(nsuzon) = nt
*
      call coepot
      call energy(2)
      potnew = pot
      enepot = enepot + potnew - potold
      write(6,8321) potold,potnew,n,nt,ndist,iob,nt0,nup,iobt
 8321 format(1x,'potold,potnew,n,nt,ndist,iob,nt0,nup,iobt=',1p2e14.6,
     &7i9)
*
*
      return
*
      end
*
*
*
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
            open(27,file='bhmergers.dat',access='append')
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
            open(27,file='bhmergers.dat',access='append')
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
            open(28,file='bhmergers-all.dat',access='append')
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
            open(29,file='bhmergers-all-all.dat',access='append')
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
      subroutine mydump(i)
*
*
*       common save or read.
*       --------------------
*
      include 'params.h'
c      include 'zdata.h'
*
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i14,i15,i16,i17,
     &        i18,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,k1,
     &        i,nco,i32,i33,i34
*
      parameter (i1=32,i2=nlagra+10,i3=8*nmax,
     &           i4=8*nbmax3,i5=350*nbmax3,i6=57,i7=26,
     &           i8=18,i9=4,i10=120,i11=7*nmax+nbmax3+nzonma+2*nsupzo+3,
     &           i12=48+4*nmax,i14=2*nmax,i15=24*nmax+24,i16=3*nmax+3,
     &           i17=nmax+1,i18=3*nmax+3,i20=3,i21=2,i22=1,i23=1,
     &           i24=6,i25=3,i26=2,i27=2,i28=1,i29=34,i30=225,i31=5,
     &           i32=200,i33=200,i34=2*nmax,k1=5*nmax)
*
      real*8  y1,y2,y3,y4,y5,y6,y7,y8,y9,y14,y15,rtidkg,smto,
     &        ys1,ys2,ys3,ys4,ys5,ys6,ys7,ys8,ys9,ys10,ys11,z1,zini
*
      real*4 timeold
*      
      integer iy10,iy11,iy12,iy13,nto,iys1,iys2,iys3,iys4,iys5,
     &        iys6,iys7,iflagns,iflagbh,itime
*
*     monte carlo commons
*
      common /param/  y1(i1)
      common /coefi/  y2(i2)
      common /body/   y3(i3)
      common /binar/  y4(i4)
      common /bininf/ y5(i5)
      common /system/ y6(i6)
      common /proba/  y8(i8)
      common /corep/  y9(i9),nco
      common /iparam/ iy10(i10)
      common /ibody/  iy11(i11)
      common /isyste/ iy12(i12)
      common /kingm/  rtidkg
      common /uptime/ y14(i14)
      common /oldpot/ smto,nto
      common /integral/ z1(k1)
      common /randx/  iy13(i29)
      common /zset/ zini
      common /fflags/ iflagns,iflagbh
      common /runtime/ timeold
      common /iruntime/ itime
*
*     stellar evolution commons
*
      common /value1/ ys1(i20)
      common /value2/ ys2(i21)
      common /value4/ ys3(i22),iys1(i23)
      common /value5/ ys4(i24)
      common /points/ ys5(i25)
      common /tstepc/ ys6(i26)
      common /params/ ys7(i27)
      common /stellar/ ys8(i15),iys6(i17)
      common /binary/ ys9(i16),iys7(i18)
      common /mscff/ ys10(i32)
      common /gbcff/ ys11(i33)
*
      common /value3/ iys2(i28)
      common /rand3/ iys3(i29)
      common /types/ iys4(i30,i30)
      common /flags/ iys5(i31)
*
*       open restart file   -  restart.fil 
*
c       open(1,file='restart.fil',status='unknown',form='formatted')
       open(1,file='restart.fil',status='unknown',form='unformatted')
*
*
*       read all common variables saved for restart
*
      if(i.eq.2) then
*
c         read (1,*)   (y1(k),k=1,i1)
c         read (1,*)   (y2(k),k=1,i2)
c         read (1,*)   (y3(k),k=1,i3)
c         read (1,*)   (y4(k),k=1,i4)
c         read (1,*)   (y5(k),k=1,i5)
c         read (1,*)   (y6(k),k=1,i6)
c         read (1,*)   (y7(k),k=1,i7)
c         read (1,*)   (y8(k),k=1,i8)
c         read (1,*)   (y9(k),k=1,i9)
c         read (1,*)   (iy10(k),k=1,i10),nco
c         read (1,*)   (iy11(k),k=1,i11)
c         read (1,*)   (iy12(k),k=1,i12)
c         read (1,*)   (iy13(k),k=1,34)
*
         read (1)   y1
         read (1)   y2
         read (1)   y3
         read (1)   y4
         read (1)   y5
         read (1)   y6
c         read (1)   y7
         read (1)   y8
         read (1)   y9
         read (1)   y14
c         read (1)   y15
         read (1)   ys1
         read (1)   ys2
         read (1)   ys3
         read (1)   ys4
         read (1)   ys5
         read (1)   ys6
         read (1)   ys7
         read (1)   ys8
         read (1)   ys9
         read (1)   ys10
         read (1)   ys11
         read (1)   rtidkg
         read (1)   smto,nto
         read (1)   z1
         read (1)   zini
         read (1)   timeold
         read (1)   iy10,nco
         read (1)   iy11
         read (1)   iy12
         read (1)   iy13
         read (1)   iys1
         read (1)   iys2
         read (1)   iys3
         read (1)   iys4
         read (1)   iys5
         read (1)   iys6
         read (1)   iys7
         read (1)   iflagns,iflagbh
         read (1)   itime
*
      endif
*
*       save all common variables needed for restart
      if(i.eq.1) then
*
c         write (1,*)   (y1(k),k=1,i1)
c         write (1,*)   (y2(k),k=1,i2)
c         write (1,*)   (y3(k),k=1,i3)
c         write (1,*)   (y4(k),k=1,i4)
c         write (1,*)   (y5(k),k=1,i5)
c         write (1,*)   (y6(k),k=1,i6)
c         write (1,*)   (y7(k),k=1,i7)
c         write (1,*)   (y8(k),k=1,i8)
c         write (1,*)   (y9(k),k=1,i9)
c         write (1,*)   (iy10(k),k=1,i10),nco
c         write (1,*)   (iy11(k),k=1,i11)
c         write (1,*)   (iy12(k),k=1,i12)
c         write (1,*)   (iy13(k),k=1,34)
*
         write (1)   y1
         write (1)   y2
         write (1)   y3
         write (1)   y4
         write (1)   y5
         write (1)   y6
c         write (1)   y7
         write (1)   y8
         write (1)   y9
         write (1)   y14
c         write (1)   y15
         write (1)   ys1
         write (1)   ys2
         write (1)   ys3
         write (1)   ys4
         write (1)   ys5
         write (1)   ys6
         write (1)   ys7
         write (1)   ys8
         write (1)   ys9
         write (1)   ys10
         write (1)   ys11
         write (1)   rtidkg
         write (1)   smto,nto
         write (1)   z1
         write (1)   zini
         write (1)   timeold
         write (1)   iy10,nco
         write (1)   iy11 
         write (1)   iy12
         write (1)   iy13
         write (1)   iys1
         write (1)   iys2
         write (1)   iys3
         write (1)   iys4
         write (1)   iys5
         write (1)   iys6
         write (1)   iys7
         write (1)   iflagns,iflagbh
         write (1)   itime
*
      endif
*
*
      close(1)
*
*
      return
*
      end
*
*
*
**
*
c      subroutine newpos(k,nup,rnew1,rnew2,iruntemp)
      subroutine newpos(k,nup)
*
*
*       determine new positions of two interacting objects after encounter
*       ------------------------------------------------------------------
*       M. HENON 1971, Astrophysics and Space Science Vol. 14, 151
*       ----------------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 e1,e2,a1,a2,q1,dz,b1,c1,gmin1,gmax1,gmin2,gmax2,
     &       rmin1,rmin2,rmax1,rmax2,rmaxz,s,u1,u2,v1,drds,q0,
     &       rtot,rles,xcut,ycut,xrr,rnn1,rnn2,ytot,qmax,rnew1,
     &       rnew2
c     &       rtot,rles,xcut,ycut,xrr,vrp,vtp,rnn1,rnn2
*
      real*4 ran2
*
      integer k,im1,im2,n,i,ipoin1,ipoin2,irm1,irm2,ipp,nup,m,imi1,
     &        imi2,ipo,kmin1,kmax1,kmin2,kmax2,imk1,imk2,maxk,
     &        mmin,mmax
*
c      common /timep/ vrp(nmax),vtp(nmax)
*
*
      n = nt
      i = 0
      im1 = iname(k)
      im2 = iname(k+1)
      gmin1 = 0.0d0
      gmin2 = 0.0d0
      gmax1 = 0.0d0
      gmax2 = 0.0d0
      u1 = 0.d0
      u2 = 0.d0
      v1 = 0.d0
      q0 = 0.d0
      ipoin1 = 0
      ipoin2 = 0
      irm1 = 0
      irm2 = 0
      imi1 = 0
      imi2 = 0
      imk1 = 0
      imk2 = 0
c      xmin(im1) = 0.0d0
c      xmin(im2) = 0.0d0
c      xmax(im1) = 0.0d0
c      xmax(im2) = 0.0d0
      xmaxz(im1) = 0.0d0
      xmaxz(im2) = 0.0d0
      xgmin(im1) = 0.0d0
      xgmin(im2) = 0.0d0
      xgmax(im1) = 0.0d0
      xgmax(im2) = 0.0d0
      rnew1 = r(k)
      rnew2 = r(k+1)
      vrp(im1) = vr(im1)
      vtp(im1) = vt(im1)
      vrp(im2) = vr(im2)
      vtp(im2) = vt(im2)
*
*       find maximal distance for new position of objects in super-zone
*
      ipp = nzst(nup)
      maxk = ipp + 1
      
*
      if(ipp.eq.n) then
        rmaxz = 1.d12
        maxk = n
      else
        rmaxz = 0.01d0*r(ipp) + 0.99d0*r(ipp+1)
      endif
*     
*       determine total energy and angular momentum for interacting objects
*
      e1 = u(k) + 0.5d0*(vr(im1)**2 + vt(im1)**2)
      e2 = u(k+1) + 0.5d0*(vr(im2)**2 + vt(im2)**2)
      a1 = r(k)*vt(im1)
      a2 = r(k+1)*vt(im2)
*
*       check for escapers  -   energy > 0
*       do not find a new possition for massles particles
*
      if(e1.gt.0.0d0.or.body(im1).eq.0.d0) then

*
c          xmin(im1) = rnew1
c          xmax(im1) = rnew1
          xmaxz(im1) = rmaxz
          ipoin1 = 1
*
       endif
*
*       do not find a new possition for massles particles
*
      if(e2.gt.0.0d0.or.body(im2).eq.0.d0) then
*
c          xmin(im2) = rnew2
c          xmax(im2) = rnew2
          xmaxz(im2) = rmaxz
          ipoin2 = 1
*
      endif
*
*
*       find rmin and rmax for interacting objects
*
*       first star
*
      if(ipoin1.eq.1) go to 30
*
cThis call attempts to find m for a simple case (the majority).  Otherwise it returns
czero, and then the standard procedure is followed.
      m = mmin(k,e1,a1,n)
      if (m.gt.0) goto 13
      i = -1
   10 i = i + 1
      m = k - i
*
      if(m.eq.0) then
        b1 = 0.0d0
        c1 = u(1)
        go to 15
      endif
*
      q1 = 2.0d0*(e1 - u(m)) - a1**2/r(m)**2
      if(q1.ge.0.0d0) go to 10
*
      if(m.eq.k) then
        imk1 = 1
        rmin1 = 0.9999d0*r(k)
        write(6,*) '< r(k)   k,im1,r(k),u(k),e1,vr,vt,body = ',
     &  k,im1,r(k),u(k),e1,vr(im1),vt(im1),body(im1)
	call flush(6)
        kmin1 = 1
        go to 17
      endif
*
      if(m.eq.n) then
        b1 = u(m)*r(m)
        c1 = 0.0d0
        go to 15
      endif
*
 13   continue
      dz = 1.0d0/r(m+1) - 1.0d0/r(m)
      b1 = (u(m+1) - u(m))/dz
      c1 = (u(m)/r(m+1) - u(m+1)/r(m))/dz
   15 rmin1 = (-b1 + sqrt(b1*b1 - 2.0d0*a1*a1*(c1 - e1)))/a1/a1
      rmin1 = 1.0d0/rmin1
*
      if(rmin1.gt.r(k)) then 
*
        xrr = 1.d0 - rmin1/r(k)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then
          write(6,*) '    rmin1  >  r(k)  iseed = ',iseed
          write(6,*) 'm,k,rmin1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1 =',
     &                m,k,rmin1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1
          call flush(6)
*
          rmin1 = 0.9999d0*r(k)
          imi1 = 1
        else
          rmin1 = r(k)
          imi1 = 1
        endif
*
      endif
*
      gmin1 = 2.0d0*a1*a1/rmin1**3 + 2.0d0*b1/rmin1**2
c      print*,gmin1,rmin1
*
      if(gmin1.lt.0.0d0) then
*
        write(6,*) '    gmin1  <   0.0   iseed = ',iseed
        write(6,*) 'm,k,rm1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1,b1,gm1=',
     &            m,k,rmin1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1,b1,gmin1

        call flush(6)
*
        gmin1 = 1.d10
      endif
*
      if(m.le.1) then
        kmin1 = 1
      else
        kmin1 = m - 1
      endif
*
      if((m.eq.n).or.(imi1.eq.1)) kmin1 = 1
*
 17   continue
*
      m = mmax(k,e1,a1,n)
      if (m.gt.0) goto 23

      i = -1
   20 i = i + 1
      m = k + i
*
      if(m.gt.n) then
        b1 = u(m-1)*r(m-1)
        c1 = 0.0d0
        go to 25
      endif
*
      q1 = 2.0d0*(e1 - u(m)) - a1**2/r(m)**2
c      print*,i,m,q1,e1,u(m),a1,r(m)
      if(q1.ge.0.0d0) go to 20
*
      if(m.eq.k) then
        if(imk1.eq.1) then
          rnew1 = r(k)
          imk1 = 2
          write(6,*) '-- O --  k,im1,r(k),u(k),e1= ',k,im1,r(k),u(k),e1
	  call flush(6)
          rmin1 = rnew1
          rmax1 = rnew1
          go to 30
        endif
*
        imk1 = 1
        rmax1 = 1.0001d0*r(k)
        write(6,*) '> r(k)    k,im1,r(k),u(k),e1 = ',k,im1,r(k),u(k),e1
	call flush(6)
        kmax1 = maxk
        go to 30
      endif
*
      if(m.eq.1) then
*
        if(imi1.eq.1) then
*
        write(6,*) ' m = 1   rnew1 = r(k)  iseed = ',iseed
        write(6,*) 'm,k,r(k),r(m),u(m),e1,a1 =',      
     &              m,k,r(k),r(m),u(m),e1,a1      
        call flush(6)
*
          imi1 = 2
          rnew1 = r(k)
          go to 30
        endif
*
        write(6,*) ' m = 1   rmax1 = 1.0001*r(k)   iseed = ',iseed
        write(6,*) 'm,k,r(k),r(m),u(m),e1,a1 =',      
     &              m,k,r(k),r(m),u(m),e1,a1      
        call flush(6)
*
        rmax1 = 1.0001d0*r(k)
        gmax1 = -1.d10
        go to 27
      endif
*
 23   continue
      dz = 1.0d0/r(m) - 1.0d0/r(m-1)
      b1 = (u(m) - u(m-1))/dz
      c1 = (u(m-1)/r(m) - u(m)/r(m-1))/dz
   25 continue
c      print*,c1,u(m-1),r(m),u(m),r(m-1),dz,m,k
      rmax1 = 0.5d0*(-b1 + sqrt(b1*b1 - 2.0d0*a1*a1*(c1 - e1)))
     &        /(c1 - e1)
*
      if(rmax1.lt.r(k)) then 
*
        if(imi1.eq.1) then
          write(6,*) 'circular orbit     iseed,k = ',iseed,k
	  call flush(6)
          imi1 = 2
          rnew1 = r(k)
          rmin1 = rnew1
          rmax1 = rnew1
          go to 30
        endif
*        
        xrr = 1.d0 - rmax1/r(k)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then
          write(6,*) ' rmax1 < r(k)  iseed =',iseed
          write(6,*) 'm,k,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1 =',
     &                m,k,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1          
          call flush(6)
*          
          rmax1 = 1.0001d0*r(k)
          imi1 = 3
        else
          rmax1 = r(k)
          imi1 = 3
        endif
*
      endif
*
      gmax1 = 2.0d0*a1*a1/rmax1**3 + 2.0d0*b1/rmax1**2
*
      if(gmax1.gt.0.0d0) then
*
        write(6,*) ' gmax1  >  0.0   iseed,gmax1 = ',iseed,gmax1
        write(6,*) 'm,k,rmax1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1 =',
     &              m,k,rmax1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1          
        call flush(6)
*          
        gmax1 = -1.d10
      endif
*
  27  if(rmax1.lt.rmin1) then
*
        write(6,*) 'rmax1 < rmin1  iseed = ',iseed
        write(6,*) 'r(k),a1,e1,rmax1,rmin1,m = ',
     &             r(k),a1,e1,rmax1,rmin1,m
        call flush(6)
*
      endif
*
      gmax1 = sqrt(-3.0d0*(rmax1 - rmin1)/gmax1)
      gmin1 = sqrt(3.0d0*(rmax1 - rmin1)/gmin1)
*
      if(rmax1.gt.rmaxz) irm1 = 1
*     
      if(m.gt.n) then
        kmax1 = n
      else
        kmax1 = m + 1
      endif
*
      if((m.eq.1).or.(imi1.eq.3)) kmax1 = maxk
*
C       print*,ipoin1,k,e1,a1,maxk,rmaxz,gmax1,gmin1,istar,rmin1,
C      &     rmax1,imi1,imk1,kmax1,kmin1,irm1
c      stop
*       second star
*
   30 continue
*
      xgmin(im1) = gmin1
      xgmax(im1) = gmax1
*
      if(ipoin2.eq.1) go to 60
*
      m = mmin(k+1,e2,a2,n)
      if (m.gt.0) goto 43
      i = -1
   40 i = i + 1
      m = k + 1 - i
*
      if(m.eq.0) then
        b1 = 0.0d0
        c1 = u(1)
        go to 45
      endif
*
      q1 = 2.0d0*(e2 - u(m)) - a2**2/r(m)**2
      if(q1.ge.0.0d0) go to 40
*
      if(m.eq.k+1) then
        imk2 = 1
        rmin2 = 0.9999d0*r(k+1)
        write(6,*) '< r(k+1) k+1,im2,r,u,e2,vr,vt=',
     &  k+1,im2,r(k+1),u(k+1),e2,vr(im2),vt(im2)
	call flush(6)
        kmin2 = 1
        go to 47
      endif
*
      if(m.eq.n) then
        b1 = u(m)*r(m)
        c1 = 0.0d0
        go to 45
      endif
*
 43   continue
      dz = 1.0d0/r(m+1) - 1.0d0/r(m)
      b1 = (u(m+1) - u(m))/dz
      c1 = (u(m)/r(m+1) - u(m+1)/r(m))/dz
   45 rmin2 = (-b1 + sqrt(b1*b1 - 2.0d0*a2*a2*(c1 - e2)))/a2/a2
      rmin2 = 1.0d0/rmin2
*
      if(rmin2.gt.r(k+1)) then
*
        xrr = 1.d0 - rmin2/r(k+1)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then
          write(6,*) '    rmin2  >  r(k+1)  iseed = ',iseed
          write(6,*) 'm,k1,rm2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2 =',
     &                m,k+1,rmin2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2
          call flush(6)
*
          rmin2 = 0.9999d0*r(k+1)
          imi2 = 1
        else
          rmin2 = r(k+1)
          imi2 = 1
        endif
*
      endif
*
      gmin2 = 2.0d0*a2*a2/rmin2**3 + 2.0d0*b1/rmin2**2
*
      if(gmin2.lt.0.0d0) then
*         
        write(6,*) '    gmin2 < 0.0   iseed,gmin2 = ',iseed,gmin2
        write(6,*) 'm,k1,rmin2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2 =', 
     &              m,k+1,rmin2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2
        call flush(6)
*         
        gmin2 = 1.d10   
      endif
*
      if(m.le.1) then
        kmin2 = 1
      else
        kmin2 = m - 1
      endif
*
      if((m.eq.n).or.(imi2.eq.1)) kmin2 = 1
*
 47   continue
*
      m = mmax(k+1,e2,a2,n)
      if (m.gt.0) go to 53

      i = -1
*
   50 i = i + 1
      m = k + 1 + i
*
      if(m.gt.n) then
        b1 = u(m-1)*r(m-1)
        c1 = 0.0d0
        go to 55
      endif
*
      q1 = 2.0d0*(e2 - u(m)) - a2**2/r(m)**2
      if(q1.ge.0.0d0) go to 50
*
      if(m.eq.k+1) then
        if(imk2.eq.1) then
          rnew2 = r(k+1)
          imk2 = 2
          write(6,*) '-- O -- k+1,im2,r,u,e2=',k+1,im2,r(k+1),u(k+1),e2
	  call flush(6)
          rmin2 = rnew2
          rmax2 = rnew2
          go to 60
        endif

        imk2 = 1
        rmax2 = 1.0001d0*r(k+1)
        write(6,*) '> r(k+1) k+1,im2,r,u,e2 = ',k+1,im2,r(k+1),u(k+1),e2
	call flush(6)
        kmax2 = maxk
        go to 60
      endif
*
 53   continue
      dz = 1.0d0/r(m) - 1.0d0/r(m-1)
      b1 = (u(m) - u(m-1))/dz
      c1 = (u(m-1)/r(m) - u(m)/r(m-1))/dz
   55 rmax2 = 0.5d0*(-b1 + sqrt(b1*b1 - 2.0d0*a2*a2*(c1 - e2)))
     &        /(c1 - e2)
*
      if(rmax2.lt.r(k+1)) then
*
        if(imi2.eq.1) then
*         
          write(6,*) ' circular orbit   iseed,k+1 =',iseed,k+1
	  call flush(6)
          imi2 = 2
          rnew2 = r(k+1)
          rmin2 = r(k+1)
          rmax2 = r(k+1)
          go to 60
        endif
*      
        xrr = 1.0d0 - rmax2/r(k+1)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then     
          write(6,*) 'rmax2 < r(k+1), rmax2=1.0001*r(k+1),iseed= ',iseed
          write(6,*) 'm,k1,rmax2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2=',
     &              m,k+1,rmax2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2            
          call flush(6)
*                   
          rmax2 = 1.0001d0*r(k+1)
          imi2 = 3
        else
          rmax2 = r(k+1)
          imi2 = 3
        endif
*
      endif
*
      gmax2 = 2.0d0*a2*a2/rmax2**3 + 2.0d0*b1/rmax2**2
*
      if(gmax2.gt.0.0d0) then
*                   
        write(6,*) '    gmax2  >  0.0   iseed,gmax2 = ',iseed,gmax2
        write(6,*) 'm,k1,rmax2,r(k+1),r(m),r(m-1),u(m),u(m-1),e2,a2 =',
     &              m,k+1,rmax2,r(k+1),r(m),r(m-1),u(m),u(m-1),e2,a2            
        call flush(6)
*                   
        gmax2 = -1.d10
      endif
*
      if(rmax2.lt.rmin2) then
*
        write(6,*) 'rmax2 < rmin2  iseed = ',iseed
        write(6,*) 'r(k+1),a2,e2,rmax2,rmin2,m = ',
     &             r(k+1),a2,e2,rmax2,rmin2,m
        call flush(6)
      endif
*
      gmax2 = sqrt(-3.0d0*(rmax2 - rmin2)/gmax2)
      gmin2 = sqrt(3.0d0*(rmax2 - rmin2)/gmin2)
*
      if(m.gt.n) then
        kmax2 = n
      else
        kmax2 = m + 1
      endif
*
      if((m.eq.2).or.(imi2.eq.3)) kmax2 = maxk
*
      if(rmax2.gt.rmaxz) irm2 = 1
*
   60 continue
*
      xgmin(im2) = gmin2
      xgmax(im2) = gmax2
*
C       print*,ipoin2,k,e2,a2,maxk,rmaxz,gmax2,gmin2,istar,rmin2,
C      &     rmax2,imi2,imk2,kmax2,kmin2,irm2
c      stop

*       determination of the new position
*
*       first star
*
      if(ipoin1.eq.1) go to 80
      xmin(im1) = rmin1
      xmax(im1) = rmax1
      xmaxz(im1) = rmaxz
      if(imi1.eq.2) go to 80
      if(imk1.eq.2) go to 80
*
      if(kmax1-kmin1.lt.2) then
        kmax1 = n
        kmin1 = 1
      endif  
*
      ipp = 0
      qmax = 1.9d0
      xcut = 0.01d0
      ycut = 1.d0 - xcut
      rles = rmin1
      rtot = rmax1
      if(irm1.eq.1) then
        rtot = rmaxz
        ycut = 1.0d0
      endif
*
   70 s = ran2(irun)
*
      if(s.lt.xcut) s = xcut
      if(s.gt.ycut) s = ycut
*
      s = -1.0d0 + 2.0d0*s
*
      rnew1 = 0.5d0*(rmin1 + rmax1) + 
     &                    0.25d0*(rmax1 - rmin1)*(3.0d0*s - s**3)
*
      ipp = ipp + 1
*
      if(ipp.lt.10000) then
*
        if((rnew1.gt.rtot).or.(rnew1.lt.rmglob)) go to 70
        if((k.eq.1).and.(rnew1.lt.r(1))) go to 70
        if(ikind(im1).lt.0) then
          print*,'newpos-1: i,im1,ikind = ',i,im1,ikind(im1)
        endif
        if(ikind(im1).eq.2.and.rnew1.lt.r(1)) go to 70
c        if(rnew1.lt.r(1)) go to 70
*
        q0 = qmax*ran2(irun)*dmax1(gmin1,gmax1)
c        print*,'q0,s ',q0,s
        ipo = 0
        rnn1 = rnew1
        call potent(ipo,kmin1,kmax1,rnn1,u1)
        v1 = 2.0d0*(e1 - u1) - a1*a1/rnew1/rnew1
*
      else
*
        write(6,*) 'ipp>10000 isee,time,rmz = ',iseed,time,rmaxz
        write(6,*) 'k,im1,rx1,ri1,v1,q0,e1,u1,a1,rnew1 = ',
     &              k,im1,rmax1,rmin1,v1,q0,e1,u1,a1,rnew1
        call flush(6)
        rnew1 = r(k)
        vrp(im1) = vr(im1)
        vtp(im1) = vt(im1)
        go to 80
*       
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1234) time,ipp,k,im1,xcut,ycut,s,rmin1,rnew1,
c     &               rmax1,rmaxz,e1,a1,u1,v1
c 1234 format(1x,'-1- ',1pe12.4,3i6,1p13e12.4)
c      close(33)
c
      if(v1.le.0.0d0) go to 70
      v1 = sqrt(v1)
      drds = 0.75d0*(rmax1 - rmin1)*(1.0d0 - s*s)
      q1 = drds/v1
*
      if(q1.gt.qmax*dmax1(gmin1,gmax1)) then
c        write(6,*) 'q1 > q  iseed,q1,gmin1,gmax1=',iseed,q1,gmin1,gmax1
        qmax = 2.0d0*qmax
        go to 70
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1235) time,ipp,k,im1,qmax,gmin1,gmax1,v1,drds,q0,q1
c 1235 format(1x,'-1a- ',1pe12.4,3i6,1p7e12.4)
c      close(33)
c      
      if(q0.gt.q1) go to 70
*
*       check for a nearly radial orbit if the new position is not
*       pick up too close to the apocentre
*
      q1 = a1/rnew1
      vrp(im1) = v1
      vtp(im1) = q1
*
*
   80 continue
*
*       second star
*
      if(ipoin2.eq.1) go to 100
      xmin(im2) = rmin2
      xmax(im2) = rmax2
      xmaxz(im2) = rmaxz
      if(imi2.eq.2) go to 100
      if(imk2.eq.2) go to 100
*
      if(kmax2-kmin2.lt.2) then
        kmax2 = n
        kmin2 = 1
      endif
*
      ipp = 0
      qmax = 1.9d0
      xcut = 0.01d0
      ycut = 1.0d0 - xcut
      rles = rmin2
      rtot = rmax2
      if(irm2.eq.1) then
        rtot = rmaxz
        ytot = 1.0d0
      endif
*
   90 s = ran2(irun)
*
      if(s.lt.xcut) s = xcut
      if(s.gt.ycut) s = ycut
*
      s = -1.0d0 + 2.0d0*s
*
      rnew2 = 0.5d0*(rmin2 + rmax2) + 
     &                       0.25d0*(rmax2 - rmin2)*(3.0d0*s - s**3)
*
      ipp = ipp + 1
*
      if(ipp.lt.10000) then
*
        if((rnew2.gt.rtot).or.(rnew2.lt.rmglob)) go to 90
        if((k+1.eq.2).and.(rnew2.lt.r(1))) go to 90
        if(ikind(im2).lt.0) then
          print*,'newpos-2: i,im2,ikind = ',i,im2,ikind(im2)
        endif
        if(ikind(im2).eq.2.and.rnew2.lt.r(1)) go to 90
c        if(rnew2.lt.r(1)) go to 90
*
        q0 = qmax*ran2(irun)*dmax1(gmin2,gmax2)
        ipo = 0
        rnn2 = rnew2
        call potent(ipo,kmin2,kmax2,rnn2,u2)
        v1 = 2.0d0*(e2 - u2) - a2*a2/rnew2/rnew2
*
      else
*
        write(6,*) 'ipp>10000 isee,t,rmz = ',iseed,time,rmaxz
        write(6,*) 'k+1,im2,rx2,ri2,v1,q0,e2,u2,a2,rnew2 = ',
     &              k+1,im2,rmax2,rmin2,v1,q0,e2,u2,a2,rnew2
        call flush(6)
        rnew2 = r(k+1)
        vrp(im2) = vr(im2)
        vtp(im2) = vt(im2)
        go to 100
*       
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1236) time,ipp,k+1,im2,xcut,ycut,s,rmin2,rnew2,
c     &               rmax2,rmaxz,e2,a2,u2,v1
c 1236 format(1x,'-2- ',1pe12.4,3i6,1p13e12.4)
c      close(33)
c
      if(v1.le.0.0d0) go to 90
      v1 = sqrt(v1)
      drds = 0.75d0*(rmax2 - rmin2)*(1.0d0 - s*s)
      q1 = drds/v1
*
      if(q1.gt.qmax*dmax1(gmin2,gmax2)) then
c        write(6,*) 'q1 > q  iseed,q1,gmin2,gmax2 =',iseed,q1,gmin2,gmax2
        qmax = 2.0d0*qmax
        go to 90
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1237) time,ipp,k+1,im2,qmax,gmin2,gmax2,v1,drds,q0,q1
c 1237 format(1x,'-2a- ',1pe12.4,3i6,1p7e12.4)
c      close(33)
c
      if(q0.gt.q1) go to 90
*
*       check for a nearly radial orbit if the new position is not 
*       pick up too close to the apocentre
*
      q1 = a2/rnew2
      vrp(im2) = v1
      vtp(im2) = q1
*
*
  100 continue
*
*
      rn1 = rnew1
      rn2 = rnew2
c      print*,rnew1,rnew2
c      stop
      return
*
      end
*
*
*
*
c      subroutine newton(x,r,r2,r1)
*
*       newton iteration for finding minimum and maximum values 
*       -------------------------------------------------------
*       of random number s
*       ------------------
*
*
c      real*8 x,r,r2,r1,ri,f,fp,xs,xs1,dx
*
c      integer iter
*
*
c      iter = 0
c      xs = 0.5d0
*
c      ri = (r - r1)/(r2 - r1)
*
c 10   continue
*
c      f = (3.0d0 - 2.0d0*xs)*xs*xs - ri
c      fp = 6.0d0*xs*(1.0d0 - xs)
*
c      xs1 = xs - f/fp
*
c      dx = abs(xs1 - xs)
c      dx = dx/xs
c      xs = xs1
*
c      iter = iter + 1
c      if(iter.gt.1000) stop  'iter > 1000'
*
c      if(dx.gt.1.0d-6) go to 10
*
c      x = xs
*
c      return
c      end
*
*
*
      function mmin(k,e,a,n)
cThis function returns the m-value for rmin in simple case, otherwise zero
      include 'common.h'
      double precision e,a,qk,q1,qstar,q2
      integer k,n,k1,k2,kstar,nit,mmin
      q1 = 2.0d0*(e - u(1)) - a**2/r(1)**2
      qk = 2.0d0*(e - u(k)) - a**2/r(k)**2
      if(q1.lt.0.d0.and.qk.gt.0.d0.and.k.gt.1.and.k.le.n) then
         k1 = 1
         k2 = k
         q2 = qk
         nit = 1
 10      continue
         if (k2.eq.k1 + 1) then
            mmin = k1
            return
         else
            nit = nit + 1
            if (nit.gt.1000) then
               print*,'nit>1000',k,e,a,n,k1,k2,q1,q2
               mmin = 0
               return
            endif
            kstar = (k1+k2)/2
            qstar = 2.0d0*(e - u(kstar)) - a**2/r(kstar)**2
            if (qstar.lt.0.d0) then
               k1 = kstar
               q1 = qstar
               go to 10
            else if (qstar.gt.0.d0) then
               k2 = kstar
               q2 = qstar
               go to 10
            else
               mmin = 0
               print*,'qstar = 0',k,e,a,n,k1,k2,q1,q2
               return
            endif
         endif
      else
         print*,'non-standard',k,e,a,n,q1,qk
         mmin = 0
         return
      endif
      end

      function mmax(k,e,a,n)
cThis function returns the m-value for rmax in simple case, otherwise zero
cNote that it is the index of the star above the maximu
      include 'common.h'
      double precision e,a,qk,qn,qstar,q1,q2
      integer k,n,k1,k2,nit,kstar,mmax
      q1 = 0.d0
      q2 = 0.d0
      qk = 2.0d0*(e - u(k)) - a**2/r(k)**2
      qn = 2.0d0*(e - u(n)) - a**2/r(n)**2
      if(qn.lt.0.d0.and.qk.gt.0.d0.and.k.ge.1.and.k.lt.n) then
         k1 = k
         k2 = n
         q1 = qk
         q2 = qn
         nit = 1
 10      continue
         if (k2.eq.k1 + 1) then
            mmax = k2
            return
         else
            nit = nit + 1
            if (nit.gt.1000) then
               print*,'nit>1000',k,e,a,n,k1,k2,q1,q2
               mmax = 0
               return
            endif
            kstar = (k1+k2)/2
            qstar = 2.0d0*(e - u(kstar)) - a**2/r(kstar)**2
            if (qstar.gt.0.d0) then
               k1 = kstar
               q1 = qstar
               go to 10
            else if (qstar.lt.0.d0) then
               k2 = kstar
               q2 = qstar
               go to 10
            else
               mmax = 0
               print*,'qstar = 0',k,e,a,n,k1,k2,q1,q2
               return
            endif
         endif
      else
         print*,'non-standard',k,e,a,n,q1,qk
         mmax = 0
         return
      endif
      end

      subroutine newvel(k,beta,den,dt)
*
*
*       determine new velocities of two interacting objects after encounter
*       -------------------------------------------------------------------
*       M. HENON 1971, Astrophysics and Space Science Vol. 14, 151
*       ----------------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 vr1,vr2,vtx1,vtx2,vty1,vty2,beta,wx,wy,
     &       wz,w,wx1,wy1,wz1,wx2,wy2,wz2,wp,psi,sinpsi,cospsi,
     &       sbeta,cbeta,wnx,wny,wnz,vnx1,vnx2,vny1,vny2,sm1,sm2,
     &       sm12,sm112,sm212,den,dt
*
      real*4 ran2
*
      integer k,im1,im2
*
      common /chesc/ sm12,sm112,sm212,w,wx,wy,wz,wx1,wy1,wz1,wx2,
     &               wy2,wz2,vtx1,vtx2,vty1,vty2,vr1,vr2,
     &               sinpsi,cospsi
*
*
      im1 = iname(k)
      im2 = iname(k+1)
*
      sm1 = body(im1)
      sm2 = body(im2)
      sm12 = sm1 + sm2
      sm112 = sm1/sm12
      sm212 = sm2/sm12
      vr1 = vr(im1)
      vr2 = vr(im2)
*
      vtx1 = vt(im1)
      vtx2 = vt(im2)*cosfi
      vty1 = 0.0d0
      vty2 = vt(im2)*sinfi
*
      wx = vtx2 - vtx1
      wy = vty2 - vty1
      wz = vr2 - vr1
      wp = wx*wx + wy*wy
      w = wp + wz*wz
      w = sqrt(w)
      wp = sqrt(wp)
*
      wx1 = wy*w/wp
      wy1 = -wx*w/wp
      wz1 = 0.0d0
      wx2 = -wx*wz/wp
      wy2 = -wy*wz/wp
      wz2 = wp
*
      psi = twopi*ran2(irun)
      sinpsi = sin(psi)
      cospsi = cos(psi)
      sbeta = sin(beta)
      cbeta = cos(beta)
*
      wnx = wx*cbeta + wx1*sbeta*cospsi + wx2*sbeta*sinpsi
      wny = wy*cbeta + wy1*sbeta*cospsi + wy2*sbeta*sinpsi
      wnz = wz*cbeta + wz1*sbeta*cospsi + wz2*sbeta*sinpsi
*
      vnx1 = vtx1 - sm212*(wnx - wx)
      vnx2 = vtx2 + sm112*(wnx - wx)
      vny1 = vty1 - sm212*(wny - wy)
      vny2 = vty2 + sm112*(wny - wy)
      vrn1 = vr1 - sm212*(wnz - wz)
      vrn2 = vr2 + sm112*(wnz - wz)
      vtn1 = sqrt(vnx1*vnx1 + vny1*vny1)
      vtn2 = sqrt(vnx2*vnx2 + vny2*vny2)
c
c      write(6,*) 'k,im1,im2,vrn1,vtn1,vrn2,vtn2 = ',k,im1,im2,vrn1,
c     &            vtn1,vrn2,vtn2
   
*
*     check for the true escapers. Method suggested by D.C. Heggie - 1996
*     do not check for tidaly limited clouster
*
cChanged for M67:
c      if(imodel.ne.3) call checke(k,dt,den)
      if(imodel.ne.3.and.imodel.ne.4) call checke(k,dt,den)
*
      return
*
      end
*
*
*
*
      
      subroutine checke(k,dt,den)
*
*
*       determine if escaping stars can realy escape from the system.
*       -----------------------------------------------------------------
*       Method proposed by D.C. Heggie - 1996.  
*       -----------------------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 vr1,vr2,vtx1,vtx2,vty1,vty2,beta,wx,wy,
     &       wz,w,wx1,wy1,wz1,wx2,wy2,wz2,sinpsi,cospsi,
     &       sbeta,cbeta,wnx,wny,wnz,vnx1,vnx2,vny1,vny2,
     &       sm12,sm112,sm212,dt,den,lambda,cden,yran,pp,po,pppo,
     &       e1,e2,vrn1n,vrn2n,vtn1n,vtn2n,sin2b2
*
      real*4 ran2
*
      integer k
*
      common /chesc/ sm12,sm112,sm212,w,wx,wy,wz,wx1,wy1,wz1,wx2,
     &               wy2,wz2,vtx1,vtx2,vty1,vty2,vr1,vr2,
     &               sinpsi,cospsi
*
*
*     check if stars escape after relaxation
*
      e1 = u(k) + 0.5d0*(vrn1**2 + vtn1**2)
      e2 = u(k+1) + 0.5d0*(vrn2**2 + vtn2**2)
*
      if((e1.gt.0.0d0).or.(e2.gt.0.0d0)) then
*
*     determine deflection angle beta
*
      cden = 0.5d0/pi
      lambda = cden*pi*den*w*dt*nt/log(gamma*nt)
      yran = ran2(irun) - 1.d0
      yran = abs(yran)
      pp = sqrt(-log(yran)/lambda)
      po = sm12/w/w
      pppo = (pp/po)**2
      sin2b2 = 1.0d0/(1.0d0 + pppo)
      beta = 2.0d0*atan(sqrt(sin2b2/(1.0d0 - sin2b2)))
*
      sbeta = sin(beta)
      cbeta = cos(beta)
*
      wnx = wx*cbeta + wx1*sbeta*cospsi + wx2*sbeta*sinpsi
      wny = wy*cbeta + wy1*sbeta*cospsi + wy2*sbeta*sinpsi
      wnz = wz*cbeta + wz1*sbeta*cospsi + wz2*sbeta*sinpsi
*
      vnx1 = vtx1 - sm212*(wnx - wx)
      vnx2 = vtx2 + sm112*(wnx - wx)
      vny1 = vty1 - sm212*(wny - wy)
      vny2 = vty2 + sm112*(wny - wy)
      vrn1n = vr1 - sm212*(wnz - wz)
      vrn2n = vr2 + sm112*(wnz - wz)
      vtn1n = sqrt(vnx1*vnx1 + vny1*vny1)
      vtn2n = sqrt(vnx2*vnx2 + vny2*vny2)
*
*     check if any star has positive energy - escapers
*
c      if((e1.gt.0.0d0).or.(e2.gt.0.0d0)) then
        vrn1 = vrn1n   
        vtn1 = vtn1n
        vrn2 = vrn2n
        vtn2 = vtn2n
      endif
*
c      e1 = u(k) + 0.5d0*(vrn1n**2 + vtn1n**2)
c      e2 = u(k+1) + 0.5d0*(vrn2n**2 + vtn2n**2)
*
c      if((e1.gt.0.0d0).or.(e2.gt.0.0d0)) then
c        vrn1 = vrn1n
c        vtn1 = vtn1n
c        vrn2 = vrn2n
c        vtn2 = vtn2n
c      endif
*
      return
*
      end
*
*
*
*
      
      subroutine output
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
     &       csb,rtid0,rtid1,xrtid,smt00,ekickbs,ekickbd
*
      real*4 runtim,cpp,xxx,timerun
*
      integer i,ii,ibirel,ibiesc,nms,nwd,nnsbh,n,ns,nb,inum,inumt,
     &        ids,ik1,ik2,ikb,idb,idstar,k,ngs,kbsm,kbsc,kbs3,kbs4,l,
     &        im1,idestr,imerge,lmin,lmax,nup,nmsl,nmsl10,nmsll,idd,
     &        iii,nbh,nbb,irtid,int,int0,lwd,lns,lbh,lhg,lgb,lch,lfag,
     &        lsag,lhms,lhhg,lhgb,l2wd,l2ns,l2bh,lwdms,lwdns,lwdbh,
     &        lwdot,lnsms,lnsbh,lnsot,lbhms,lbhot,lbh3,lbh2
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
      if((tphys.gt.tte).or.(time.gt.tcrit)) then
        nup = -100
        if(tphys.le.tcrevo) then
          tte = tte + dtte0
        else
          tte = tte + dtte
        endif
        print*,'tphys,timet,time,tcrit,tcrevo,tte,dtte0,dtte ',
     &          tphys,timet,time,tcrit,tcrevo,tte,dtte0,dtte
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
     &         ekicktbs - ekicktbd     
        write(6,6543) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &                ekickt,ekicktbs,ekicktbd,enekin,escsta,ehbin3,
     &                ehb3b3        
 6543    format(1x,'output-m etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     & 'eccol,ekickt,ekicktbs,ekicktbd,enekin,escsta,ehbin3,ehb3b3 =',
     & 1p14e20.12)
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
     &        ekicktbs - ekicktbs     
      write(6,6544) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,        
     &              ekickt,ekicktbs,ekicktbd,enekin,escsta,ehbin3,
     &              ehb3b3                 
 6544    format(1x,'output-e etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     & 'eccol,ekickt,ekicktbs,ekicktbd,,enekin,escsta,ehbin3,ehb3b3 =',
     & 1p14e20.12)                
*
*       calculate conservation of the total energy
*
c      call energy(2)
*
      etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 - 
     &        ehmlev - ehcoll + eccoll - ekickt - enekin -
     &        ekicktbs - ekicktbd
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
         if(ikind(im1).eq.2) nbb = nbb + 1
         if(ibstra(im1).eq.1) kbsm = kbsm + 1
         if(ibstra(im1).eq.2) kbsc = kbsc + 1
         if(ibstra(im1).eq.3) kbs3 = kbs3 + 1
         if(ibstra(im1).eq.4) kbs4 = kbs4 + 1
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
      print*,' ekick,ekickbs,ekickbd = ',ekick,ekickbs,ekickbd
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
      write(13,110) iseed,time,tphys,smt,etot,zkink,pot,escsta,ehbint,
     &              ehb3b3,ehmlev,ehcoll,eccoll,error,enepot,sloses,
     &              slosev,slosco,rlag(11),rtid,nt,nescst,nescm,nmloev,
     &              ncoll,nbin3-nescb3-ndist3-ndist4-ndiste-nmerge,kbsm,
     &              kbsc,kbs3,kbs4,kbsm+kbsc+kbs3+kbs4,ibsm+ibsc+ibs3+
     &              ibs4,ibsm,ibsc,ibs3,ibs4,nt0,ivnewg,ivrr,enrad,
     &              timerun,ehbin3,txxx,ekickt,ekicktbs,ekicktbd,
     &              ekicktb2,ikickt,ikicktbs,ikicktbd,mnsbh,nexchang,
     &              nexchang2,nbb,de,csb,nescrt,sescrt,lwd,lns,lbh,lhg,
     &              lgb,lch,lfag,lsag,lhms,lhhg,lhgb,l2wd,l2ns,l2bh,
     &              lwdms,lwdns,lwdbh,lwdot,lnsms,lnsbh,lnsot,lbhms,
     &              lbhot,lbh2,lbh3,ntsn1,ntsn2,ntsnb
*
  110 format(1x,i5,1p7e12.4,1pe15.7,1p11e12.4,17i10,2i12,1p8e12.4,7i7,
     &       1p2e12.4,i9,1pe12.4,28i9)
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
      xktc = smc*vc*vc/(3.0d0*nc)
      smvc = 0.d0
      do 111 i = 1,5
      ii = iname(i)
 111  smvc = smvc + body(ii)
*
      ro5 = 1.5d0*smvc/twopi/r(5)**3
*
      write(14,120) iseed,time,tphys,smc,rc,vc,roc,u(1),
     &              -3.d0*u(1)/vc/vc,r(5),ro5,nc
*
  120 format(1x,i5,1p10e12.4,i6)
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
        open(18,file='bin3inf.dat',access='append')
        open(19,file='bin3int.dat')
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
             if (iprint.eq.0) then
                write(18,170) iseed,time,tphys,im1,i,
     &               (bin(i,iii),iii=1,6),bin(i,7),bin(i,8),xktc,rc
             endif
           endif
*
           inum = iinte3(i)
           if(inum.eq.0) go to 150
           do 160 ii = 1,inum
              inumt = 50*(i-1) + ii
              write(19,180) iseed,time,tphys,im1,i,ii,
     &                      (binin(inumt,iii),iii=1,7)
 160       continue
 150    continue
*
        write(17,140) iseed,time,tphys,nbin3,nbin3-nescb3-ndist3-
     &          ndist4-ndiste-nmerge,nb3b3,ndist3,ndist4,nmerge,ndiste,
     &                nb3fin,nesb3s,nescb3,ibiesc,ibirel,idestr,imerge,
     &                ehbin3,eintb3,escb3s,escbi3,escbb3,erelb3,erb3in,
     &                ehb3b3,slob3s,slob3b,nexchang,nexchang2
  140   format(1x,i5,1p2e12.4,14i7,1p10e12.4,2i7)
        close(17)
*
 170    format(1x,i5,1p2e12.4,2i8,1p6e12.4,1pe24.16,1pe12.4,1p2e12.4)
 180    format(1x,i5,1p2e12.4,3i8,1p7e12.4)
        close(18)
        close(19)
*
      endif
c
*
*     snapshot, profile
*
      write(6,*)  ' snapshot '
      call flush(6)
*
      if((tphys.eq.0.d0).or.(tphys.ge.ttp).or.time.gt.tcrit) then
*
      open(21,file='snapshot.dat',access='append')
*      open(21,file='snapshot.dat',position='append')
*
      ttp = ttp + dttp
      print*,'tphys,ttp,dttp = ',tphys,ttp,dttp
      xxn = zmbar*rtidkg/rbar
      vscal = 0.06558d0*sqrt(xxn)

cIn the following (DCH) the temperature is in units of (4*pi*sigma)**-0.25
      write(21,210) n,tphys
 210  format(1x,i9,1pe12.4,'  *****')
*
      do 200 i = 1,n
         k = iname(i)
         if(r(i).gt.rtid) then
           print*,'r>rt tim,r,rt=',tphys,r(i),rtid
           go to 200
         endif
         if(abs(body(k)).lt.tolm) go to 200
         if(ikind(k).lt.0) then
           print*,'output: i,k,ikind = ',i,k,ikind(k)
         endif
         if(ikind(k).eq.1.or.ikind(k).ge.3) then
           ids = names(k)
           idd = ids
cAdded DCH 9/9/6
           call getLum(ids,slum1)
           call getSpin(ids,spin1)
           spin2 = 0.d0
           slum = slum1
           if(isNaN(slum)) then
             print*,'i,im1,ids,r,lum = ',i,k,ids,r(i),slum1
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
           call getRadius(idstar,rad1)
           call get_hiid(idb,idstar)
           call getEpoch(idstar,epoch2)
           call getSpin(idstar,spin2)
           call get_hiid_mass(idb,sm2)
           call get_ss_type(idstar,ik2)
           call getLum(idstar,slum2)
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
           mv = -2.5d0*log10(10.d0**(-0.4d0*mv1) + 10.d0**(-0.4d0*mv2))
           mb = -2.5d0*log10(10.d0**(-0.4d0*mb1) + 10.d0**(-0.4d0*mb2))
           mi = -2.5d0*log10(10.d0**(-0.4d0*mi1) + 10.d0**(-0.4d0*mi2))
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
           print*,'output k,body(k)*zmbar,temp,time: ',
     &            k,body(k)*zmbar,temp,time
         endif
 200  continue
*
      close(21)
*
      write(6,*)  ' in profiles '
      call flush(6)
*            
      call profiles(tphys,vscal)
*
      write(6,*)  ' out profiles '
      call flush(6)
*
*       save restart data
*
      call mydump(1)
*
      endif
*
      write(6,*)  ' out output',timerun
      call flush(6)
*            
      if(abs(de).gt.qe) then
        write(6,*)'energy conservation is not satisfied! ',time,etot,de
        stop 'energy error'
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
      if(time.gt.tcrit) then
        write(6,*) 'time is greater than tcrit   ',time
        stop 'time > tcrit'
      endif
*
*     stop when the number of objects is too small
*
      if(nt*gamma.le.3.d0) then
        write(6,*) '     N*gamma < 3.0    ',nt,time
        stop 'N*gamam < 3.0 '
      endif
*
      return
*
      end
*
*
*
*
      subroutine profiles(tphys,vscal)
*
*     Code for generating surface brightness and other profiles
*      --------------------------------------------------------
*     DCH 24-25/4/6
*     -------------
*
*
      include 'common.h'
*
      integer i,j,k,ngrid,ids,idb,idstar,n
      parameter (ngrid = 31)
      real*8  rlog(ngrid),lum,mass,sb(ngrid),sumrhovr2(ngrid),
     &     sumrho(ngrid),rgrid(ngrid),lumstar,massstar,radius,
     &     v_radial,v_trans,sinth,costh,tphys,vscal,rtt
c      data rlog,sumrhovr2,sumrho /ngrid*0.d0,ngrid*0.d0,ngrid*0.d0/
c      data sb /ngrid*0.d0/
*
*    First create the grid in radius [(0.005  0)*rtid in code units)
*
      do i = 1,ngrid
         rlog(i) = (-2.3 + 2.3*float((i-1))/float((ngrid - 1)))*rtid
         rgrid(i) = 10.d0**rlog(i)*rbar/rtidkg
         sb(i) = 0.d0
         sumrhovr2(i) = 0.d0
         sumrho(i) = 0.d0
      enddo
*
*     Loop through every object, and calculate its contribution to each 
*     radius
*
cAdded DCH 1/8/6
      n = nt
      rtt = rtid*rbar/rtidkg
      print *,'profiles: n = ',n
      do 10 i = 1,n
         k = iname(i)
         lum = 0.d0
         if(ikind(k).lt.0) then
           print*,'profiles: i,k,ikind = ',i,k,ikind(k)
         endif
         if (ikind(k).eq.1.or.ikind(k).ge.3) then
            ids = names(k)
            call getLum(ids,lum)
            call getMass(ids,mass)
         else if(ikind(k).eq.2) then
            idb = nameb(k)
            call getIprim(idb,idstar)
            call getLum(idstar,lum)
            call getMass(idstar,mass)
            call getIsec(idb,idstar)
            call getLum(idstar,lumstar)
            call getMass(idstar,massstar)
            lum = lum + lumstar
            mass = mass + massstar
         else 
            print*,'profiles: i,k,ikind = ',i,k,ikind(k)
         endif
         radius = r(i)*rbar/rtidkg
         v_radial = vr(k)*vscal
         v_trans = vt(k)*vscal
         do  20 j = 1,ngrid
            if(rgrid(j).gt.rtt) go to 20
            if (rgrid(j).lt.radius) then
               sinth = rgrid(j)/radius
               costh = sqrt(1.d0 - sinth**2)
               if(isNaN(lum)) then
            print*,'i,k,j,lum,rad,sin,cos=',i,k,j,lum,radius,sinth,costh
                 go to 20
               endif
               sb(j) = sb(j) + lum/(2*pi*radius**2*costh)
               sumrhovr2(j) = sumrhovr2(j) + 
     &         mass*(0.5*v_trans**2*sinth**2 + v_radial**2*costh**2)/
     &         (2*pi*radius**2*costh)
               sumrho(j) = sumrho(j) + mass/(2*pi*radius**2*costh)
            endif
 20      continue
 10   continue
*
      open(22,file='profile.dat',access='append')
*      open(22,file='profile.dat',position='append')
*      
      write (22,30) tphys
      do j = 1,ngrid
         write (22,40) rgrid(j),sb(j),sumrhovr2(j)/sumrho(j)
      enddo
*   
 30   format(1x,1pe12.4)
 40   format(1x,1p3e12.4)
*
      close(22)
*      
      return
*
      end
*
*
*
      subroutine mycolour(kw,sm,slum,rad,mv,mb,mi)           
      implicit double precision (a-h,o-z)
      double precision mv,mb,logz,logte,logr,logl,massi,mi,bminv,
     &                 bolc,rad,slum,sm,uminb,vmini,zini,z
      integer ifirst,kw
      common /zset/ zini
      data ifirst/0/
      save ifirst
       z = zini
c      z = 0.002d0
c      z = 0.02d0
      logz = log10(z/0.02)
*
* Read in the bolometric correction tables. 
*
      if (ifirst.eq.0) then
         CALL ubvtab
         ifirst = 1
      endif
      logr = log10(rad)
      logl = log10(slum)
      logte = 3.762d0 + 0.25d0*logl - 0.5d0*logr
      massi = sm
      if(kw.ge.10.and.kw.le.12)then
         CALL wd2ubv(logl,logte,massi,bolc,mv,
     &               uminb,bminv,vmini)
      else
         CALL lt2ubv(logl,logte,massi,logz,bolc,mv,
     &               uminb,bminv,vmini)
      endif
      mb = mv + bminv
      mi = mv - vmini
*      print*,'mycolour: ',mb,mv
      return
      end
*
*
      subroutine potent(ipot,k1,k2,ri,ui)
*
*
*         compute potential 'ui' for radius 'ri'
*         --------------------------------------
*
*
      include 'common.h'
*
      real*8 ri,ui,z,z1,z2,z3,z4,z5
*
      integer i,n,ipot,k1,k2
*
*
      z = 1.d0/ri
*
      if(ipot.eq.0) then
        n = nt
*
        call locate(r,n,k1,k2,ri,i)
*
        if(i.eq.n) then
          ui = -smt*z
          return
        endif
*
        if(i.eq.0) then
          ui = u(1)
          return
        endif
*
        z1 = 1.0d0/r(i)
        z2 = 1.0d0/r(i+1)
        z3 = z1 - z2
        z4 = (u(i+1) - u(i))/z3
        z5 = z1 - z
        ui = u(i) + z5*z4
        return
*
      else
*
        n = nto
*
        call locate(ro,n,k1,k2,ri,i)
*
        if(i.eq.n) then
          ui = -smto*z
          return
        endif
*
        if(i.eq.0) then
          ui = uo(1)
          return
        endif
*
        z1 = 1.0d0/ro(i)
        z2 = 1.0d0/ro(i+1)
        z3 = z1 - z2
        z4 = (uo(i+1) - uo(i))/z3
        z5 = z1 - z
        ui = uo(i) + z5*z4
        return
*
      endif
*        
*
      end
*
*
*
*    
      subroutine locate(xx,n,k1,k2,z,j)
*
      include 'params.h'
*
      real*8 xx(nmax),z
*
      integer n,k1,k2,j,jl,jm,ju
*      
      jl = k1 - 1
      ju = k2 + 1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(z.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl
      return
      end
*
*
*
*

      function ran2(idum)
*
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
*
      real ran2,am,eps,rnmx
*
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *    ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     *    ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
*
      integer idum2,j,k,iv,iy
*
      common /randx/ idum2,iy,iv(ntab)
*
      if(idum.le.0)then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
           k=idum/iq1
           idum=ia1*(idum-k*iq1)-k*ir1
           if(idum.lt.0)idum=idum+im1
           if(j.le.ntab)iv(j)=idum
 11     continue
        iy=iv(1)
      end if

      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
*
      if(idum.lt.0)idum=idum+im1
*
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
*
      if(idum2.lt.0)idum2=idum2+im2
*
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
*
      if(iy.lt.1)iy=iy+imm1
*
      ran2=min(am*iy,rnmx)
*
      return
      end
*
*
*

      subroutine relax(nup,lmin,lmax,dt)
*
*
*
*       determine relaxation process for stars between lmin and lamx
*       ------------------------------------------------------------
*       with time-step dt
*       -----------------
*
*
      include 'common.h'
*
      real*8 sm1,sm2,w,den,a12,sin2b2,dt1,beta,dt,dt0,rx1,rx2,ek,ek0,
     &       ep,ep1,ep2,voo,dt2
*
      real*4 ran2
*
      integer i,im1,im2,lmin,lmax,nup,nc2,ibound,m,k,ibx
*
*
      ibound = 0
      ep = 0.0d0
      ek = 0.0d0
      ep1 = 0.0d0
      ep2 = 0.0d0
      dt0 = dt
      nc2 = nc/2
      if(nc2.gt.nminzo) nc2 = nminzo
*
*      determine deflection angle and new velocities
*
      do 10 i=lmin,lmax,2
*
         im1 = iname(i)
*
         vro(im1) = vr(im1)
         vto(im1) = vt(im1)
         ro(i) = r(i)
         rn(i) = r(i)
         if(i+1.gt.lmax) go to 10
*
         im2 = iname(i+1)
         sm1 = body(im1)
         sm2 = body(im2)
cAdded DCH 2/8/6
         if (sm1*sm2.eq.0.d0) go to 10
*
*       do not relax pair of stars which apocentre distances are 
*       smaller then actual position of the next star
*
c        if(i.lt.nc2) then
        if(i.lt.0) then
          if(xmax(im1).eq.0.0d0) go to 15
          rx1 = xmax(im1)
          rx2 = xmax(im2)
          if(rx1.gt.xmaxz(im1)) rx1 = xmaxz(im1)
          if(rx2.gt.xmaxz(im2)) rx2 = xmaxz(im2)
          if((rx1.lt.r(i+2)).and.(rx2.lt.r(i+2))) go to 30
        endif
*
 15     continue
*
*       determine relative velocity of stars i and i+1
*
         if(ran2(irun).lt.0.5) vr(im1) = -vr(im1)
         if(ran2(irun).lt.0.5) vr(im2) = -vr(im2)
*
         call relvel(i,0,w)
*
*       determine density in the vicinity of stars i and i+1
*
         call densi(i,den) 
*
*       determine the deflection angle
*
         a12 = den*(sm1 + sm2)**2/w**3*float(nt)
         sin2b2 = a12*dt0
* 
*       if beta greater than PI form a sequence of interaction
*
         dt1 = dt0
   20    continue
         if(sin2b2.gt.1) then
*
           dt2 = 0.5d0/a12
c           dt2 = 1.0d0/a12
           dt1 = dt1 - dt2
*
*       find new velocities after interaction with deflection angle = PI/2
*
           beta = 0.5d0*pi
c           beta = pi
           call newvel(i,beta,den,dt2)
*
           vr(im1) = vrn1
           vr(im2) = vrn2
           vt(im1) = vtn1
           vt(im2) = vtn2
*
           call relvel(i,0,w)
*
*       determine new deflection angle
*
           a12 = den*(sm1 + sm2)**2/w**3*float(nt)
           sin2b2 = a12*dt1
c           dt1 = dt2
           go to 20
         endif
*
         beta = 2.0d0*atan(sqrt(sin2b2/(1.0d0-sin2b2)))
*
*       determine velocities after interaction
*
         call newvel(i,beta,den,dt1)
*
         vr(im1) = vrn1
         vr(im2) = vrn2
         vt(im1) = vtn1
         vt(im2) = vtn2
*
*       compute the self-energy of up to nminzo innermost stars
*       if it is negative devide total kinetic energy of these stars
*       on two equal parts and isotropise the velocities
*
 30     continue
        ibound = 1
        if(ibound.eq.0) then
          if(i.le.nminzo) then
            ek0 = 0.5d0*(sm1*(vr(im1)**2 + vt(im1)**2) +
     &                   sm2*(vr(im2)**2 + vt(im2)**2))
            ek = ek + ek0
            ep1 = ep1 + 0.5d0*(sm1**2/r(i) + sm2**2/r(i+1))
            do 40 k=2,i+1
            do 50 m=k,i+1
               ep2 = ep2 + body(iname(k-1))*body(iname(m))/r(m)
 50         continue
 40         continue
            ep = ep1 + ep2
            ep2 = 0.0d0
            if(ep.lt.ek) then
              ibound = 1
              if((i+1).eq.2) go to 10
              ek0 = (ek - ek0)/float(i-1)
              do 60 k = 1,i-1
                 im1 = iname(k)
cChanged by DCH 2/8/6
                 if (body(im1).ne.0.d0) then
                    voo = 2.0d0*ek0/3.0d0/body(im1)
                    vt(im1) = sqrt(2.0d0*voo)
                    vr(im1) = sqrt(voo)
                 endif
 60           continue
*
              write(6,*) 'bound num,ep,ek,lmi,lma =',i-1,ep,ek,lmin,lmax
*
              go to 10
            endif
          endif
        endif
*
 10   continue  
*
*
      do 70 i=lmin,lmax,2
*
*       define velocities of stars after interaction as old velocities
*
         if(i+1.gt.lmax) go to 70
*
         im1 = iname(i)
         im2 = iname(i+1)         
*
         vro(im1) = vr(im1) 
         vro(im2) = vr(im2)
         vto(im1) = vt(im1)
         vto(im2) = vt(im2)
*
*       determine new positions of two interacting stars
*
c
c         rrroo1 = r(i)
c         rrroo2 = r(i+1)

         call newpos(i,nup)
*
*       define positions of stars before interaction as old position
*
         ro(i) = r(i)
         ro(i+1) = r(i+1)
         rn(i) = rn1
         rn(i+1) = rn2
*
*        set new position for binaries (binary array)
*
         if(ikind(im1).ge.2) then
c           ibx = im1 - nss0
c           if(ibx.lt.0) then
           ibx = nwhich(im1)
c           endif
           bin(ibx,8) = rn1
           bin(ibx,6) = time + dt
         endif
*
         if(ikind(im2).ge.2) then
c           ibx = im2 - nss0
c           if(ibx.lt.0) then
           ibx = nwhich(im2)
c           endif
           bin(ibx,8) = rn2
           bin(ibx,6) = time + dt
         endif
*                                                                                    *
c
c         if(rn1.gt.rtid.or.rn2.gt.rtid) then
c      print*,'rrr: im1,im2,i,i+1,ro1,ro2,rn1,rn2,2(rmax,rmin,rmaxz)=',
c     &im1,im2,i,i+1,rrroo1,rrroo2,rn1,rn2,xmin(im1),xmax(im1),
c     &xmaxz(im1),xmin(im2),xmax(im2),xmaxz(im2)
c         endif
*
 70   continue
*
*
      return
*
      end
*
*
*
*
       subroutine relaxt
*
*
*       compute relaxation process for all sub-models which consist of 
*       ---------------------------------------------------------------
*       one model of the whole cluster
*       ------------------------------
*
*
      include 'common.h'
*
      real*8 dt,eko,xx,vrrt,tform,tb3f,tb3b3,tcoll,tescp(20),tim
*
      integer ld,nmodel,k,m,l,j,nup,lmin,lmax,n,mm,i,ll,nm
*
      common /timset/ tform(20),tb3f(20),tb3b3(20),tcoll(20)
*
*       determine which super-zone is due to relaxation and calculate
*       in sequence all partial models which consist of the whole
*       cluster model
*      
*
      n = nt
      ld = ltwo(1) - ltwo(nsuzon)
      nmodel = 2**ld
*
*     set the time for later use in the interaction procedures. The
*     vector txx(j) are set separatly for each iteraction procedure
*
       do j = 1,20
          tform(j) = timet
          tb3f(j) = timet
          tb3b3(j) = timet
          tcoll(j) = timet
          tescp(j) = timet
       enddo
*
      do 10 i = 1, nmodel
*
      write(6,*) ' --- relaxt ---  i,nmodel = ',i,nmodel
      call flush(6)
*
         k = i
         m = mod(k,2)
*
         if(m.gt.0) then
           nup = 1
*
         else
*
           j = 0
   20      k = k/2
           j = j + 1
           m = mod(k,2)
           if(m.eq.0) go to 20
           nup = j + 1
         endif
*
         do 30 l=1,nup
*
            if(l.eq.1) then
              lmin = 1
            else
              lmin = nzst(l-1) + 1
            endif
*
            lmax = nzst(l)
*
          write(6,*) 'l,nup,lmin,lmax = ',l,nup,lmin,lmax
          call flush(6)
*
*        lmin and lmax are the first and the last stars in the 
*        super zone which is evoluated with time-step dt        
*
            if(ltwo(l).ne.0) then
              dt = tau*2.0d0**(-ltwo(l))
            else
              dt = tau
            endif
*
            tescp(l) = tescp(l) + dt*tscale0*log(gamma*nt00)
     &              /log(gamma*ntnew)*float(ntnew)/float(nt00)
            print*,'escape- nup,tescp(l) = ',l,nup,tescp(l)
            call flush(6)
            tim = tescp(l) 
*
*        calculate relaxation process for stars between lmin and lmax
*
            write(6,*) ' ------ in  relax ------'
            call flush(6)
*
            call relax(nup,lmin,lmax,dt)
*
            write(6,*) ' ------ out  relax ------'
            call flush(6)
*
   30    continue
*
         do 40 l=1,n
         inameo(l) = l
   40    uo(l)=u(l)
*
*        sort paritcles according to increase distsnce      
*
         do 50 l=1,lmax
   50    r(l) = rn(l)
*   
         do 60 l=lmax+1,n
         mm = iname(l)
         vro(mm) = vr(mm)
         vto(mm) = vt(mm)
   60    ro(l) = r(l)
*
         call sort3(n,r,iname,inameo)
*
         smto = smt
         nto = nt
*
*        compute new potential
*
            write(6,*) ' ------ in  coepot ------'
            call flush(6)
*
c         call energy(2)
c         etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
c     &   ehmlev - ehcoll + eccoll - ekickt - enekin
c         write(6,6542) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
c     &               ekickt,enekin,escsta,ehbin3,ehb3b3
c 6542    format(1x,'relax-1 etotn,zkin,pot,ehmlev,enepot,ehcoll,',
c     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)
*
         call coepot
c         call energy(2)
c         etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -  
c     &   ehmlev - ehcoll + eccoll - ekickt - enekin               
c         write(6,6541) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,            
c     &               ekickt,enekin,escsta,ehbin3,ehb3b3                       
c 6541    format(1x,'relax-2 etotn,zkin,pot,ehmlev,enepot,ehcoll,',            
c     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)                
*
            write(6,*) ' ------ out  coepot ------'
            call flush(6)
*
*        determine new velocities from changes of the potential
*
            write(6,*) ' ------ in  timepot ------'
            call flush(6)
*
         call timpot(lmax)
*
c         call energy(2)
c         etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -  
c     &   ehmlev - ehcoll + eccoll - ekickt - enekin               
c         write(6,6540) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,            
c     &               ekickt,enekin,escsta,ehbin3,ehb3b3                       
c 6540    format(1x,'timepot-1 etotn,zkin,pot,ehmlev,enepot,ehcoll,',
c     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)                
*
            write(6,*) ' ------ out  timepot ------'
            call flush(6)
*
*        deal with energy errors caused by the potential adjustment procedure
*
         ivrr = ivrr + lmax
         vrrt = 0.0d0
         eko = 0.0d0
         do 70 ll = 1,n
            nm = iname(ll)
            if(vrr(nm).ge.0.d0) go to 75
              ivnewg = ivnewg + 1
              enrad = enrad + vrr(nm)
 75         continue
            vrrt = vrrt + vrr(nm)
            vrr(nm) = 0.0d0
            eko = eko + 0.5d0*body(nm)*(vr(nm)**2 + vt(nm)**2)
 70      continue
         xx = (eko + vrrt)/eko
         xx = sqrt(xx)
         do 80 ll = 1,n
            nm = iname(ll)
            vr(nm) = vr(nm)*xx
            vt(nm) = vt(nm)*xx
 80      continue
*
         call energy(2)
         etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
     &   ehmlev - ehcoll + eccoll - ekickt - enekin
         write(6,6539) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &               ekickt,enekin,escsta,ehbin3,ehb3b3       
 6539    format(1x,'timepot-2 etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)
*
*        determine whether binary is due to form or generate
*        energy in interactions with field stars or binary-binary
*        interactions. 
*
            write(6,*) ' ------ in  inb3b3 ------'
            call flush(6)
*
         call inb3b3(nup)
*
            write(6,*) ' ------ out  inb3b3 ------'
            call flush(6)
*
         n = nt
*
            write(6,*) ' ------ in  intb3f ------'
            call flush(6)
*
         call intb3f(nup)
*
            write(6,*) ' ------ out  intb3f ------'
            call flush(6)
*
         n = nt
*
         write(6,*) ' ------ in  intcol ------'
         call flush(6)
*
         call intcol(nup)
*
         write(6,*) ' ------ out  intcol ------'
         call flush(6)
*
         n = nt
*
            write(6,*) ' ------ in  formb3 ------'
            call flush(6)
*
*        determine formation of binaries in three-body interactions
*
         call formb3(nup)
*
            write(6,*) ' ------ out  formb3 ------'
            call flush(6)
*
         n = nt
         iesc = 0
*
*        deal with escapers after strong interactions
*
            write(6,*) ' ------ in  escape after interactions ------'
            call flush(6)
*
         call escape(lmax,tim)
*
            write(6,*) ' ------ out  escape after interactions ------'
            call flush(6)
*
*        compute new potential after escapers removal
*
         if(iesc.ne.0) then
           call sort3(n,r,iname,inameo)
           n = nt
           nzst(nsuzon) = n
           call coepot
         endif
*
         call energy(2)
         etotn = zkin - pot + escsta - ehbin3 + enepot - ehb3b3 -
     &         ehmlev - ehcoll + eccoll - ekickt - enekin
         write(6,6543) etotn,zkin,pot,ehmlev,enepot,ehcoll,eccoll,
     &               ekickt,enekin,escsta,ehbin3,ehb3b3
 6543    format(1x,'escape-f etotn,zkin,pot,ehmlev,enepot,ehcoll,',
     &'eccol,ekickt,enekin,escsta,ehbin3,ehb3b3 =',1p12e20.12)
*
   10    continue
*
*        determine the evolution time
*
         time = time + tau
         timet = timet + tau*tscale0*log(gamma*nt00)/log(gamma*nt)*
     &           float(nt)/float(nt00)
         ntnew = nt
*
         write(6,*) '  time,tau,timet,ntnew = ',time,tau,timet,ntnew
         call flush(6)
*
*
         return
*
         end
*
*
*
*
      subroutine relvel(kn,kk,w)
*
*
*       determine relative velocity of two interacting stars
*       ----------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 w,vr1,vr2,vtx1,vtx2,vty1,vty2,wx,wy,wz,vx,vy,vz,
     &       sm1, sm2, sm12
*
      real*4 ran2
*
      integer k,kk,kn,im1,im2
*
*
      k = abs(kn)
      im1 = iname(k)
*
      if(kk.eq.0) then
        im2 = iname(k+1)
      else
        im2 = iname(kk)
      endif
*
      vr1 = vr(im1)
      vr2 = vr(im2)
      sm1 = body(im1)
      sm2 = body(im2)
      sm12 = sm1 + sm2
*
      fi = twopi*ran2(irun)
      sinfi = sin(fi)
      cosfi = cos(fi)
*
      vtx1 = vt(im1)
      vtx2 = vt(im2)*cosfi
      vty1 = 0.0d0
      vty2 = vt(im2)*sinfi
*
      wx = vtx2 - vtx1
      wy = vty2 - vty1
      wz = vr2 - vr1
      w = wx*wx + wy*wy + wz*wz
      w = sqrt(w)
*
*     compute ceter of mass velocity components for mergrer
*     k < 0 means merger 
*
      if(kn.lt.0) then
        vx = (sm1*vtx1 + sm2*vtx2)/sm12
        vy = (sm1*vty1 + sm2*vty2)/sm12
        vz = (sm1*vr1 + sm2*vr2)/sm12
        vr(im1) = vz
        vt(im1) = sqrt(vx*vx + vy*vy)
      endif
*
      return
*
      end
*
*
*
*
      
      real function runtim(x)
*
*
*       calculte elapsed time in hours
*       ------------------------------
*
*        for g77 and gfortran
*
      real*4 x,secnds
*
      runtim = secnds(x)
*
      return
*
      end
*
*
*
*
      subroutine scale0
*
*
*       scale to the N-body units.
*       ---------------------
*
      include 'common.h'
*
*
      real*8 qv,e0,sx,v2,vdr,vd,a1,xxn,trh,zm,zm1,aRsun,ecc,rscal
*
      real ran2
*
      integer i,k,n,ibx,nm,id1,id2,in,ib,im,iexist,in1,in2
*
*
      n=nt
*
*       scale masses to standard units of <m> = 1/n
*
      do 10 i = 1,n
          im = iname(i)
          if(ikind(im).eq.2) then
             body(im) = body(im)/smt
             ibx = nwhich(im)
             bin(ibx,1) = bin(ibx,1)/smt
             bin(ibx,2) = bin(ibx,2)/smt
             bin(ibx,3) = bin(ibx,3)*rtid*rsuntopc/rbar
             bin(ibx,5) = bin(ibx,1)*bin(ibx,2)/2.d0/bin(ibx,3)
             ehbi3p = ehbi3p + bin(ibx,5)
          else
             body(im) = body(im)/smt
          endif
   10 continue
*
      smt = 1.0d0
*
*       obtain the total kinetic & potential energy
*
      call energy(1)
*
*       scale non-zero velocities by virial theorem ratio
*
      if (zkin.gt.0.0d0) then
*
          qv = sqrt(qvir*pot/zkin)
          do 20 i = 1,n
             im = iname(i)
             do 30 k = 1,3
                xdot(im,k) = xdot(im,k)*qv
   30        continue
   20     continue
*
      end if
*
*       scale total energy to standard units (e = -0.25 for q < 1)
*
      e0 = -0.25d0
      etot = (qvir - 1.0d0)*pot
      sx = e0/etot
*
*
      if(iprint.eq.0)
     & write (6,65)  sx, etot, body(1), body(n), smt/float(n)
   65 format (//'scaling:   sx =',1pe11.3,'  e =',1pe11.3,
     &       '  m(1) =',1pe11.3,'  m(n) =',1pe11.3,'  <m> =',1pe11.3)
*
*       scale coordinates & velocities to the new units
*
      do 40 i = 1,n
         im = iname(i)
         do 50 k = 1,3
            x(im,k) = x(im,k)/sx
            xdot(im,k) = xdot(im,k)*sqrt(sx)
   50    continue
   40 continue
*
*
*       compute radius, radial and tangential velocities
*
        do 60 i=1,n
           im = iname(i)
           r(i)=sqrt(x(im,1)**2 + x(im,2)**2 + x(im,3)**2)
           if(ikind(im).eq.2) then
             ibx = nwhich(im)
             bin(ibx,8) = r(i)
           endif
           v2=xdot(im,1)**2 + xdot(im,2)**2 + xdot(im,3)**2
           vdr=x(im,1)*xdot(im,1)+x(im,2)*xdot(im,2)+x(im,3)*xdot(im,3)
           vd = vdr*vdr/r(i)/r(i)
           vr(im) = sqrt(vd)
           if(ran2(irun).gt.0.5) vr(im) = -vr(im)
           vd = v2 - vd
           vt(im) = sqrt(vd)
*
*       check option for writing the initial conditions on unit 6
*
   60     continue
*
*       compute the potential for all particles
*
      call energy(2)
      etot = zkin - pot
*
*       set initial crossing time in scaled units
*
      tcr = smt**2.5/(2.0d0*abs(etot))**1.5
*
*       obtain approximate half-mass radius after scaling
*
      rscale = 0.25d0*smt**2/abs(etot)
*
*       form equilibrium rms velocity (temporarily defined as vc)
*
      vc = sqrt(2.0d0*abs(etot)/smt)
*
*       print scale radius equilibrium rms velocity 
*       half-mass relaxation time & equilibrium crossing time
*
      a1 = float(n)
      trh = 4.0d0*twopi/3.0d0*(vc*rscale)**3
     &      /(15.4d0*smt**2*log(a1)/a1)
      if(iprint.eq.0) then
         write (6,80)  rscale,vc,trh, tcr, 2.0*rscale/vc,zkin,pot
   80    format (/,1x,'rscale = ',1pe10.2,'   vc = ',1pe10.2,'   trh ='
     &       ,1pe10.2,'  tcr =',1pe10.2,'  2<r>/<v> =',1pe10.2,
     &       '  zkin, pot = ',1p2e12.4)
      endif
*
*       determine core parameters
*
      call core
*
      if(iprint.eq.0) then
        write(6,90) smc,nc,rc,vc,roc,u(1)
  90    format(/,' core:  mass = ',1pe12.4,'  number of stars = ',i6,
     &        '  radius = ',1pe12.4,'  velocity = ',1pe12.4,
     &        '  density = ',1pe12.4,'  potential = ',1pe12.4)
      endif
*
*    set the time scale factor for the physical units
*
      xxn = float(nt00)/log(gamma*nt00)
*
      if(imodel.lt.3) then
*
         tscale = 20.837*rbar**1.5*xxn/sqrt(zmbar)
*
         rscal = rbar
         rtidkg = 1.d0
      else
*
         tscale = 14.90976*(rbar/rtid)**1.5*xxn/sqrt(zmbar)
*
         rscal = rbar/rtidkg
         print*,'rbar,rtidkg,rscal,zmbar,tscale = ',
     &           rbar,rtidkg,rscal,zmbar,tscale    
      endif
      tscale0 = tscale
      tolm = 1.0d-5/zmbar
*
*    single star and binary initialization (McScatter interface)
*
      in = 1
      do i = 1,n
         nm = iname(i)
         if (ikind(nm).eq.1) then
            id1 = names(nm)
c            in = 1
            zm = body(nm)*zmbar
            call init_stars(in,id1,zm)
            call flush(6)
            call get_ss_updatetime(id1,uptime(nm))
c            print*,'scale0: i,nm,id1,zm upt= ',i,nm,id1,zm,uptime(nm)
c            call flush(6)
         else if (ikind(nm).eq.2) then
            id2 = nameb(nm)
c            in = 1
            ib = nwhich(nm)
            zm = bin(ib,1)*zmbar
            zm1 = bin(ib,2)*zmbar
            ecc = bin(ib,4)
            aRsun = bin(ib,3)*rbar/rtid/rsuntopc
            call init_binaries(in,id2,aRsun,ecc,zm,zm1)
            call get_loid(id2,in1)
            call get_hiid(id2,in2)
             call get_bs_updatetime(id2,uptime(nm))
            call binary_exists(id2,iexist)
            if(iexist.eq.0) then
               ikind(nm) = 3
            endif
c          print*,'i,nm,binary_id,loid,hiid,id2,zm,zm1,ecc,upt,iexist',
c     &           '= ',i,nm,ib,in1,in2,id2,zm,zm1,ecc,uptime(nm),iexist
c            call flush(6)
         else
            write (6,*) 'ikind not 1 or 2 in data.f'
            stop
         endif
      enddo
*
      call flush(6)
*
      return
*
      end
*
*
*
*      subroutine sort2(n,ra,rb)
*
*
*       heapsort method (press p. 231).
*       -------------------------------
*
      real*8 ra,rra
*
      integer rb,rrb,l,n,ir,i,j
*
      dimension  ra(n),rb(n)
*
*
      l = n/2+1
      ir=n
   10 continue
*
      if(l.gt.1)then
	  l=l-1
	  rra=ra(l)
	  rrb=rb(l)
      else
          rra=ra(ir)
	  rrb=rb(ir)
	  ra(ir)=ra(1)
	  rb(ir)=rb(1)
          ir=ir-1
*
	  if(ir.eq.1)then
	      ra(1)=rra
	      rb(1)=rrb
	      return
          endif
*
      endif
*
      i=l
      j=l+l
*
   20 if(j.le.ir)then
*
	  if(j.lt.ir)then
	      if(ra(j).lt.ra(j+1))j=j+1
          endif
*
	  if(rra.lt.ra(j))then
	       ra(i)=ra(j)
	       rb(i)=rb(j)
	       i=j
	       j=j+j
           else
	       j=ir+1
           endif
*
           go to 20
      endif
*
      ra(i)=rra
      rb(i)=rrb
*
      go to 10
*
      end
*
*
*
*
      subroutine sort3(n,ra,rb,rc)
*
*
*       heapsort method (press p. 231).
*       -------------------------------
*
      real*8 ra,rra
*
      integer rb,rrb,rc,rrc,l,n,ir,i,j
*
      dimension  ra(n),rb(n),rc(n)
*
*
      l = n/2+1
      ir=n
   10 continue
*
      if(l.gt.1)then
	  l=l-1
	  rra=ra(l)
	  rrb=rb(l)
          rrc=rc(l)
      else
          rra=ra(ir)
	  rrb=rb(ir)
          rrc=rc(ir)
	  ra(ir)=ra(1)
	  rb(ir)=rb(1)
          rc(ir)=rc(1)
          ir=ir-1
*
	  if(ir.eq.1)then
	      ra(1)=rra
	      rb(1)=rrb
              rc(1)=rrc
	      return
          endif
*
      end if
*
      i=l
      j=l+l
*
   20 if(j.le.ir)then
*
	  if(j.lt.ir)then
	      if(ra(j).lt.ra(j+1))j=j+1
          endif
*
	  if(rra.lt.ra(j))then
	       ra(i)=ra(j)
	       rb(i)=rb(j)
               rc(i)=rc(j)
	       i=j
	       j=j+j
           else
	       j=ir+1
           endif
*
           go to 20
      endif
*
      ra(i)=rra
      rb(i)=rrb
      rc(i)=rrc
*
      go to 10
*
      end
*
*
*
*
      subroutine start
*
*
*       initialization of initial model
*       -------------------------------
*
      include 'common.h'
*
      real*4 runtim,cpp      
*
      integer n
*
*
*       inicialize cpu time
*
      cpp = runtim(0.0)
      cpu = cpp
      itime = 0
      timeold = 0.0
*
      n = nt
      ntnew = nt
*
*       initialize global scalars,arrays, counters & useful constants
*
      call zero
*
*       set initial conditions: body(i), r(i), vr(i), vt(i)
*
      call data
*
*       scale initial conditions to N-body units
*
      call scale0
*
*       define mean mass in scaled units 
*
      bodym = smt/float(n)
*
*
      return
*
      end
*
*
*
*
*
*
      subroutine timpot(lmax)
*
*
*       determine changes to the system due to changes of mass distribution
*       -------------------------------------------------------------------
*
*     timpot is called before interaction sequence (intb3b3 - intb3f - formb3)
*     so ikind cannot be negative
*
*
      include 'common.h'
*
c      real*8 uoo,uon,uno,unn,vold,vnew,dd,v2t,vrp,vtp
      real*8 uoo,uon,uno,unn,vold,vnew,dd,v2t
*  
      integer lmax,i,k1,k2,im,imo,ipo,iex,im1
*
c      common /timep/ vrp(nmax),vtp(nmax)
*
      do 10 i=1,lmax
         im = iname(i)
         imo = inameo(i)
*
*     determine old potentials for old and new positions
*
         uoo = uo(imo)
         ipo = 1
         k1 = 1
         k2 = lmax
         call potent(ipo,k1,k2,r(i),uon)
*
*     determine new potentials for old and new positions
*
         ipo = 0
         k1 = 1
         k2 = lmax
         call potent(ipo,k1,k2,ro(imo),uno)
         unn = u(i)
*
*     determine old and new velocities
*
         vold = vro(im)*vro(im) + vto(im)*vto(im)
         vnew = vold + uoo - uon - unn + uno
* 
*     determine new radial and tangential velocities
*
         vt(im) = vto(im)*ro(imo)/r(i)
         v2t = vt(im)*vt(im)
         dd = vnew - v2t
         if(body(im).lt.0.d0) then
           write(6,*) ' m < 0  i,im,imo,m = ',i,im,imo,body(im)
           stop ' mass < 0 '
         endif
*
         if(dd.gt.0.0d0) then
           vr(im) = sqrt(dd)
           vrr(im) = 0.0d0
         else   
*
           if(vnew.lt.0.0d0) then
             vr(im) = 0.0d0
             vrr(im) = 0.5d0*dd*body(im)
           else
             dd = xmax(im)/r(i) - 1.0d0
             v2t = r(i)/xmin(im) - 1.0d0
             if((dd.lt.0.03d0).or.(v2t.lt.0.03d0)) then
               vr(im) = 0.0d0
               vt(im) = sqrt(vnew)
               vrr(im) = 0.0d0
             else
               dd = vrp(im)**2 + vtp(im)**2
               dd = sqrt(dd/vnew)
               vr(im) = 0.0d0
               vt(im) = sqrt(vnew)
               vr(im) = vrp(im)/dd
               vt(im) = vtp(im)/dd
               vrr(im)= 0.0d0
             endif
           endif
         endif
*
*     for the nminzo stars the velocities are artificialy isotropised
*
c 20      continue
c  20     iex = 3*nminzo
c  20     iex = 0.1*nt
         iex = -10
         if(i.le.iex) then
         im1 = iname(i)
         if(ikind(im1).lt.0) then
           print*,'timpot: i,im1,ikind = ',i,im1,ikind(im1)
         endif
           if(ikind(im1).eq.2) go to 10
           dd = (vr(im)**2 + vt(im)**2)/3.0d0
           vr(im) = sqrt(dd)
           vt(im) = sqrt(2.0d0*dd)
         endif
*
 10   continue
*
*
      return
*
      end
*
*
*
*
*
      subroutine veloci(k,velave,smave)
*
*       determine average velocity and mass of stars in the vicinity
*       ------------------------------------------------------------
*       of 'k' star
*       -----------
*
*
      include 'common.h'
*
      real*8 velave,smave,v2
*
      integer k,i,i1,i2,im,n,nspan,inum
*
*
      velave = 0.0d0
      smave = 0.0d0
      n = nt
c      nspan = int(1.5*nminzo/2.0)
c      nspan = 10 
       nspan = 4
*
      if(k.gt.nspan) then
        i1 = k-nspan
        i2 = k+nspan
*
        if(i2.ge.n) then
          i1 = n - 2*nspan - 1
          i2 = n
        endif
*
      else
        i1 = 1
        i2 = 2*nspan + 1
      endif
*
*      compute average velocity and mass of stars from 2*nspan + 1 stars
*
      inum = 0
      do 10 i = i1,i2
         im = iname(i)
         inum = inum + 1
c         print*,'i,im,inum,vr,vt,sms=',i,im,inum,vr(im),vt(im),body(im)
         v2 = vr(im)**2 + vt(im)**2
         velave = velave + v2*body(im)
         smave = smave + body(im)
 10   continue
*
      velave = velave/smave
      smave = smave/float(inum)
*
      return
      end
*
*
*
*
      subroutine veloc(k,velave,smave)
*
*       determine average velocity and mass of stars in the vicinity
*       ------------------------------------------------------------
*       of 'k' star
*       -----------
*
*
      include 'common.h'
*
      real*8 velave,smave,v2
*
      integer k,i,i1,i2,im,n,nspan,inum
*
*
      velave = 0.0d0
      smave = 0.0d0
      n = nt
c      nspan = int(1.5*nminzo/2.0)
c      nspan = 10 
       nspan = 4
*
      if(k.gt.nspan) then
        i1 = k-nspan
        i2 = k+nspan
*
        if(i2.ge.n) then
          i1 = n - 2*nspan - 1
          i2 = n
        endif
*
      else
        i1 = 1
        i2 = 2*nspan + 1
      endif
*
*      compute average velocity and mass of stars from 2*nspan + 1 stars
*
      inum = 0
      do 10 i = i1,i2
         im = iname(i)
         inum = inum + 1
         print*,'i,im,inum,vr,vt,sms=',i,im,inum,vr(im),vt(im),body(im)
         v2 = vr(im)**2 + vt(im)**2
         velave = velave + v2*body(im)
         smave = smave + body(im)
 10   continue
*
      velave = velave/smave
      smave = smave/float(inum)
*
      return
      end
*
*
      subroutine zero
*
*
*       initialization of global scalars. 
*       ---------------------------------
*
      include 'common.h'
*
      integer i,j
*
*      integer idum2,iv,iy,i
*
c      common /randx/ idum2,iy,iv(32)
*
*      data idum2/123456789/, iv/32*0/, iy/0/
*
*
*       initialize parameters, arrays, counters & set useful constants.
*
      time    = 0.0d0
      timet   = 0.0d0
      ttp     = 0.0d0
      tte     = 0.0d0
      rtid    = 0.0d0
      ekickt  = 0.0d0
      ekicktbs  = 0.0d0
      ekicktbd  = 0.0d0
      ekicktb2  = 0.0d0
      escsta  = 0.0d0
      escbi2  = 0.0d0
      escbi3  = 0.0d0
      escbi4  = 0.0d0
      escbip  = 0.0d0
      escbb2  = 0.0d0
      escbb3  = 0.0d0
      escbb4  = 0.0d0
      escbbp  = 0.0d0
      escb2s  = 0.0d0
      escb3s  = 0.0d0
      escb4s  = 0.0d0
      escbps  = 0.0d0
      ehbin2  = 0.0d0
      ehbin3  = 0.0d0
      ehbi3p  = 0.0d0
      ehb3b3  = 0.0d0
      ehbin4  = 0.0d0
      ehbinp  = 0.0d0
      ecbin2  = 0.0d0
      ehmlev  = 0.0d0
      erelb3  = 0.0d0
      erb3in  = 0.0d0
      enepot  = 0.0d0
      enekin  = 0.0d0
      error   = 0.0d0
      enrad   = 0.0d0
      sloses  = 0.0d0
      slob3s  = 0.0d0
      slob3b  = 0.0d0
      slosev  = 0.0d0
      slosco  = 0.0d0
      sescrt  = 0.0d0
      eccoll  = 0.0d0
      ehcoll  = 0.0d0
      pbin2   = 0.0d0
      pbin3   = 0.0d0
      pbin4   = 0.0d0
      pbin2s  = 0.0d0
      pbin3s  = 0.0d0
      pbin4s  = 0.0d0
      pbinps  = 0.0d0
      pb2b2   = 0.0d0
      pb2b3   = 0.0d0
      pb2b4   = 0.0d0
      pb2bp   = 0.0d0
      pb3b3   = 0.0d0
      pb3b4   = 0.0d0
      pb3bp   = 0.0d0
      pb4b4   = 0.0d0
      pb4bp   = 0.0d0
      pbpbp   = 0.0d0
      pcoll   = 0.0d0
*
      mnsbh   = 0
      ikickt  = 0
      ikicktbs  = 0
      ikicktbd  = 0
      nexchang = 0
      ntsn1 = 0
      ntsn2 = 0
      ntsnb = 0
      nmerge  = 0
      nescst  = 0
      nescb2  = 0
      nescb3  = 0
      nescb4  = 0
      nescrt  = 0
      nescbp  = 0
      nesb2s  = 0
      nesb3s  = 0
      nesb4s  = 0
      nesbps  = 0
      nmloev  = 0
      nmloco  = 0
      iobt = 0
      ncoll   = 0
      necb2   = 0
      nehb2   = 0
      necb3   = 0
      nehb3   = 0
      necb4   = 0
      nehb4   = 0
      necbp   = 0
      nehbp   = 0
      nbin2   = 0
      nbin3   = 0
      nb3b3   = 0
      nb3fin  = 0
      nbin4   = 0
      nbinp   = 0
      ndist3  = 0
      ndist4  = 0
      iesc1   = 0
      ivnewg  = 0
      ivrr    = 0
      ibsm    = 0
      ibsc    = 0
      ibs3    = 0
      ibs4    = 0
*
      do 20 i=1,32
         iv(i) = 0
 20   continue
*
      idum2 = 123456789
      iy = 0
*
*       set fractional constants pi & two pi
*
      one2  = 1.0d0/2.0d0
      one3  = 1.0d0/3.0d0
      one4  = 1.0d0/4.0d0
      one5  = 1.0d0/5.0d0
      one6  = 1.0d0/6.0d0
      one9  = 1.0d0/9.0d0
      one12 = 1.0d0/12.0d0
      pi = 4.0d0*atan(1.0d0)
      twopi = 2.0d0*pi
      rmglob = 1.8d0/float(nt)
      rsuntopc = 6.9599e10/3.085678e18
*
      do 30 i = 1,nmax
         nkick(i) = 0
         inexch(i) = 0
         iname(i) = 0
         ikind(i) = 0
         names(i) = 0
         nameb(i) = 0
         nwhich(i) = 0
         nbinar(i) = 0
         ibstra(i) = 0
         body(i)= 0.d0
         vr(i) = 0.d0
         vt(i) = 0.d0
         r(i) = 0.d0
         u(i) = 0.d0
         xescp(i) = 0.d0
         xesct(i) = 0.d0
         vrr(i) = 0.d0
         ro(i) = 0.d0
         vro(i) = 0.d0
         vto(i) = 0.d0
         uo(i) = 0.d0
         rn(i) = 0.d0
         xmin(i) = 0.d0
         xmax(i) = 0.d0
         xmaxz(i) = 0.d0
         xgmin(i) = 0.d0
         xgmax(i) = 0.d0
         x(i,1) = 0.d0
         x(i,2) = 0.d0
         x(i,3) = 0.d0
         xdot(i,1) = 0.d0
         xdot(i,2) = 0.d0
         xdot(i,3) = 0.d0                                                                
         uptime(i) = 0.d0
         oldtime(i) = 0.d0
         vrp(i) = 0.d0
         vtp(i) = 0.d0
 30   continue
*
      do 40 i = 1,nbmax3
         iinte3(i) = 0
         do 41 j = 1,8
             bin(i,j) = 0.d0
 41      continue
 40   continue
*
      do 50 i =1,50*nbmax3
         do 51 j = 1,7
            binin(i,j) = 0.d0
 51      continue
 50   continue
*
      return
*
      end
*
*
*
*
      subroutine zone
*
*
*       compute devision of the system on spherical zones and grup
*       ----------------------------------------------------------
*       them in super-zones. Compute individual time steps for
*       ------------------------------------------------------
*       each super-zone.
*       ----------------
*
*
*
      include 'common.h'      
*
*
      real*8 smz,w2,dt,beta,xld,az,tref,tccc,xxl,bmin1,bmax1,vv,rr
*
      integer n,i,nzon,num,nco,ibound,nud,lo,ln,ld,k,im,itwo
*
      dimension az(nzonma),vv(nzonma),rr(nzonma)
*
*
      n = nt
*
*       determine number of zones and stars in each zone
*
      if(time.eq.0.0d0) then
        nzon0 = nt/nz0
        izon = 0
        nzon = nzon0
      else
        nzon = n/nz0
        if(izon.eq.1) go to 10
*
        if(nzon.lt.0.9*nzon0) then
          nz0 = nz0 - 2
          izon = 0
          nzon = n/nz0
*
          if(nz0.lt.nminzo) then
            nz0 = nminzo
            nzon = n/nz0
            izon = 1
          endif
*
        endif
*
      endif
*
   10 continue
*
*       determine the name of the last star in each zone
*
      i = 0
   20 i = i + 1
      if(i.gt.nzon-1) go to 30
      nzs(i) = nz0*i
      go to 20
*
   30 nzs(i) = n 
*
*
*       check if in the core is at least nzonc zones
*
      i = 0
   40 i =i + 1
      if(nzs(i).lt.nc) go to 40
*
      if(i.ge.nzonc) go to 90
*
      ibound = i
      nco = nzs(ibound)
      num = nco/nzonc
*
        if(num.gt.nminzo) then
*
          nud = nzonc - ibound
          k = mod(num,2)
*
            if(k.gt.0) then
              num = num - 1
            endif
*
          do 50 i=nzon,ibound+1,-1
   50        nzs(i+nud) = nzs(i)
*
          do 60 i=1,nzonc-1
   60        nzs(i) = num*i
*
          nzs(nzonc) = nco
          nzon = nzon + nud
*
        else
* 
          num = nzs(ibound)/nminzo
          nud = num - ibound
*
          do 70 i=nzon,ibound+1,-1
   70        nzs(i+nud)=nzs(i)
*
          do 80 i=1,num-1
   80        nzs(i)= nminzo*i
*
          nzs(num) = nco
          nzon = nzon + nud
*
        endif
*     
*
   90 continue
*
*       determine for each zone coefficient in formulae for sin^2(beta/2)
*       Stodolkiewicz 1982, A.A. Vol.32, 69, eq.(20)
*
      lo = 1   
*
      do 200 i=1,nzon
         ln = nzs(i)
         ld = ln + 1 - lo
         smz = 0.0d0
         w2 = 0.0d0
*
      do 210 k=lo,ln
         im=iname(k)
         smz = smz + body(im)
  210    w2 = w2 + body(im)*(vr(im)*vr(im) + vt(im)*vt(im))
*
         xld = float(ld)
         smz = smz/xld
         w2 = 2.0d0*w2/smz/xld
         vv(i) = sqrt(0.5d0*w2)
         rr(i) = 0.5d0*(r(lo) + r(ln))
*
         az(i) = 6.0d0*nt*xld*smz*smz/w2**1.5/(r(ln)**3 - r(lo)**3)
*
         if(i.eq.2) then
           if(az(1).gt.2.d0*az(2)) az(1) = 1.9d0*az(2)
         endif
*
         if(i.gt.3) then
           if((az(i-2).lt.az(i-1)).and.(az(i).lt.az(i-1))) 
     &                     az(i-1) = 0.5d0*(az(i-2) + az(i))
         endif
*
         if(i.gt.3) then
           if((az(i-2).gt.az(i-1)).and.(az(i).gt.az(i-1))) 
     &                     az(i-1) = 0.5d0*(az(i-2) + az(i))
         endif
*
         lo = ln + 1
* 
  200 continue
*
*       determine zones for which time-step is the same and collect them
*       in a one super-zone.
*       determine time step for each super-zone
*
c      tau = tau0
      xxl = log(gamma*nt)/nt
c      tccc = rc*xxl/vc
      bmin1 = bmin
      bmax1 = bmax
*
  300 k = 0
*
      do 310 i=1,nzon
         itwo = 0
*
  320    itwo = itwo + 1
         if(itwo.gt.ntwo) go to 370
         dt = az(i)*2.0d0**(-itwo)
         beta = dt*tau
*
         if((itwo.eq.1).and.(beta.lt.bmin1)) then
           itwo = 0
           go to 325
         endif
*
         if(beta.gt.bmax1) go to 320
*
  325    continue
*
c         if(itwo.gt.0) then
c           tref = 2.0d0**(-itwo)*tau
c           tccc = rr(i)*xxl/vv(i)
c           if(tref.lt.1.5d0*tccc) then
c           if((k.le.3).and.(tref.lt.1.5d0*tccc)) then
c	      bmin1 = bmax1
c	      bmax1 = 2.0d0*bmax1
c	      open(23,file='crostime.dat',access='append')
c              write(23,*) time,tref,tccc,bmax,tau,bmax1,itwo,k
c	      close(23)
c	      tau = tau0
c	      go to 300
c            endif
c         endif
*
         ln = nzs(i)
         k = k + 1
         nzst(k) = ln
         ltwo(k) = itwo
*
         if(k.eq.1) go to 310
*
         ld = ltwo(k) - ltwo(k-1)
*
         if(ld) 330,340,350
*
  330    if(-ld.eq.1) go to 310
         ltwo(k) = ltwo(k-1) - 1
         go to 310
*
  340    nzst(k-1) = nzst(k)
         k = k - 1
         go to 310
*
  350    continue
* 
         nzst(k-1) = nzst(k)
         k = k - 1
         go to 310
*
  370    tau = tau/2.0d0
         go to 300
*
  310 continue
*
      nsuzon = k
c      nsuzon = ltwo(1)
*
*       print control data
*
c      if(iprint.eq.0) then
        write(6,*) ' nzon = ',nzon,'  nsuzon = ',nsuzon,'  tau = ',tau
*
c        do 400 k=1,nzon
c           bbb = az(k)*tau*2.d0**(-nsuzon+1)
c  400      write(6,*) k,nzs(k),az(k),bbb
*
        do 410 k=1,nsuzon
  410      write(6,*) k,nzst(k),ltwo(k)
*
c      endif
*
*
      return
*
      end
*
*
*
************************************************************************
*
*      INTERFACE_BSE.F 
*
*      STELLAR AND BINARY EVOLUTION MODULE 
*      [comprised of interface_bse.f, interface_bse.h, const_bse.h, 
*       & libstr.a]
*      Version 0.0 March 20, 2003 by JH
c     Version 3.1 June 14, 2004 by DCH
*
*      Requires libstr.a - BSE library of stellar and binary evolution 
*      functions/prescriptions described in the papers: 
*
*      "Comprehensive Analytic Formulae for Stellar Evolution as a 
*       Function of Mass and Metallicity", 
*       Hurley J.R., Pols O.R. and Tout C.A. 2000, MNRAS, 315, 543. 
*       (the SSE paper). 
*
*      "Evolution of Binary Stars and the Effect of Tides on Binary Populations", 
*       Hurley J.R., Tout C.A. and Pols O.R. 2002, MNRAS, 329, 897
*       (the BSE paper). 
*
*      [note that libstr.a was compiled with either g77 or fort 
*       - the Compaq Fortran compiler - on Linux (or Mac OSX) 
*       and using the -O2 option in all cases]
*      [contact jhurley@amnh.org regarding any problems with the module]
*
*      The subroutines evStar and evBinary provide the link between the 
*      main program and the SSE/BSE library (or module). These routines 
*      evolve the star (or binary) forward by a user specified interval. 
*      The library routines return a recommended update timestep based on 
*      the evolution stage of the star, and the parameters dmmax (maximum 
*      allowed change in mass) and drmax (maximum allowed change in radius). 
*      The main program may or may not utilise this timestep. 
*
*      Input options for the SSE/BSE library are set in the subroutine initPar 
*      and these may be modified by the user (see the routine for an 
*      explanation of these options). The common blocks that convey these 
*      options to the library are declared in const_bse.h. 
*      Single stars are initialized by calling initStar and binaries by 
*      calling initBinary (which in turn initializes two single stars). 
*
*      The following quantities are available for star i: 
*
*          Age     - time to which star has been evolved (Myr)
*          MStime  - main sequence lifetime (Myr)
*          Tstep   - recommended time between updates (Myr)
*          Epoch   - effective zero-age of star (Myr)(age = Age - Epoch)
*          Mass    - current stellar mass (Msun)
*          Mass0   - initial stellar mass (Msun)
*          Radius  - current stellar radius (Rsun)
*          Lum     - current stellar luminosity (Lsun) 
*          Y0      - initial helium mass fraction (0->1)(Y = 0.24 + 2*Z)
*          Z0      - initial metallicity mass fraction (0->1)(0.02 is solar)
*          Y       - current helium mass fraction (0->1)
*          Z       - current metallicity mass fraction (0->1)
*          Mc      - core mass (Msun)
*          Rc      - core radius (Rsun)
*          Menv    - mass of the convective envelope (Msun)
*          Renv    - radius of the convective envelope (Msun)
*          Spin    - stellar rotation (1/yr)
*          Rl      - Roche-lobe radius (Rsun) 
*          Pos     - x,y,z co-ordinates (?)
*          Vel     - velocity in x,y,z (km/s)
*          Type    - index of stellar type (0->15, e.g. MS star = 1)
*
*      and these additional quantities are available for binaries: 
*
*          Btype   - index of binary type (1->6, e.g. detached = 1)
*          Iprim   - array index of the primary star (0->nmax)
*          Isec    - array index of the primary star (0->nmax)
*          Mprim   - mass of the primary (primary has greater radius/rl)
*          Msec    - mass of the secondary 
*          Semi    - semi-major axis of the orbit (Rsun) 
*          Ecc     - eccentricity of the orbit 
*          Tb      - orbital period (day)
*
*      Note that quantity X may be obtained by a call to the subroutine getX 
*      and may be set by a call to setX. The arrays that store these quantities 
*      are declared in interface_bse.h (where the user may choose to alter the 
*      size of these arrays, i.e. nmax). 
*
*      Additional subroutines included in this interface: 
*
*         ssupdatetime - provides next update time for star i 
*         bsupdatetime - calls ssupdatetime using index of primary star 
*         binaryexists - determines if binary remains bound based on Btype 
*         getLabel     - text label associated with index of stellar type 
*         getLabelb    - text label associated with index of binary type 
*         printstar    - formatted output for a single star 
*         printbinary  - formatted output for a binary 
* 
************************************************************************
*
       SUBROUTINE initPar
       implicit none
*
* Input options for the SSE/BSE library functions 
* (see BSE paper for details on most of the options). 
*
       include "const_bse.h"
       include "interface_bse.h"
       integer i,iflagns,iflagbh
       common /fflags/ iflagns,iflagbh
*
* Random number seed used to determine supernovae velocity kicks. 
* [NOTE THAT THIS SHOULD IDEALLY BE COMMUNICATED FROM MAIN PROGRAM]
       idum = -999
*
* Flag for use in common-envelope routine (not currently utilised).
       ceflag = 0
*
* Flag to activate tidal circularisation (0=off; 1=on).
       tflag = 1
*
* Flag to choose which WD IFMR to use (0=SSE; 1=HPE, 1995, MN, 272, 800).
       ifflag = 0
*
* Flag to determine NS/BH mass (0=SSE; 1=Belczynski et al. 2002, ApJ, 572, 407).
       nsflag = 1
*
* Flag to choose WD cooling track (0=simple-Mestel; 1="modified-Mestel").
       wdflag = 1
*
* Flag to allow kicks at birth for BHs (0=no; 1=yes).
       bhflag = iflagbh
*
* Maximum NS mass (set at 1.8 Msun for nsflag=0; 3.0 Msun for nsflag=1).
       mxns = 3.d0
*
* Reimers mass-loss coefficient (neta*4x10^-13; =0.5 normally).
       neta = 0.5d0
*
* Binary enhanced mass-loss parameter (inactive for single stars).
       bwind = 0.d0
*
* Common-envelope efficiency parameter (1.0).
       alpha1 = 1.d0
*
* Binding energy factor for common-envelope evolution (0.5).
       lambda = 0.5d0
*
* Dispersion in the Maxwellian for the SN kick speed (190 km/s).
       if (iflagns.eq.0) then
          sigma = 0.d0
       else
          sigma = 190.d0
       endif
       print*,'initpar: iflagns,sigma',iflagns,sigma
*

* Wind velocity parameter (proportional to V_wind**2; 0.125).
       beta = 0.125d0
*
* Wind accretion efficiency factor for momentum transfer (1.0).
       xi = 1.d0
*
* Bondi-Hoyle accretion factor (1.5).
       acc2 = 1.5d0
*
* Fraction of matter retained in a nova eruption (0.001).
       epsnov = 0.001d0
*
* Eddington limit factor for mass transfer (1.0; set large to allow super-Edd).
       eddfac = 1.d0
*
* Parameter to determine angular momentum change during RLOF mass-loss 
* (>0: lost material carries with it a fraction gamma of orbital angular momentum;
*  -1: material carries with it specific angular momentum of the primary; 
*  -2: material is lost from system as if a wind from the secondary).
       gamma = -1.d0
*
* Parameters to determine the timestep based on the evolution phase, 
* i.e. dt = pts*t_phase (pts1=MS,HG,HeMS; pts2=GB,CHeB,HeGB; pts3=AGB).
       pts1 = 0.05d0
       pts2 = 0.01d0
       pts3 = 0.02d0
*
* Maximum allowed fractional change in mass per timestep.
       dmmax = 0.05d0
*
* Maximum allowed fractional change in radius per timestep.
       drmax = 0.1d0
*
       aursun = 214.95d0
       yeardy = 365.24d0
*
* Set the collision matrix (ktype). 
*
       CALL instar
*
* Initialize the stellar index array to negative value. 
*
       do i = 0,nmax
          CALL setType(i,-1)
       enddo
*
       RETURN
       END
***
       SUBROUTINE printStar(idstar)
       integer idstar,kw,iend
       real*8 mass,rad
       character*20 label
*
       CALL getType(idstar,kw)
       CALL getLabel(kw,label,iend)
       CALL getMass(idstar,mass)
       mass = MIN(mass,999.d0)
       CALL getRadius(idstar,rad)
       rad = MIN(rad,9999.d0)
       WRITE(6,500)label(1:iend),mass,rad
  500  FORMAT(a,',  mass = ',f12.8,' radius = ',f13.8)
*
       RETURN
       END
***
       SUBROUTINE printBinary(idbin)
       implicit none
       integer idbin,idstar,ip,is
       integer k,kw,iend1,iend2,iend3
       real*8 time,epoch1,epoch2
       real*8 m1,r1,m2,r2,a,ecc,v(3),vd
       character*20 label1,label2,label3
*
       CALL getIprim(idbin,idstar)
       ip = idstar
       CALL getType(idstar,kw)
       CALL getLabel(kw,label1,iend1)
       CALL getMass(idstar,m1)
       m1 = MIN(m1,999.d0)
       CALL getRadius(idstar,r1)
       r1 = MIN(r1,9999.d0)
*
       CALL getAge(idstar,time)
       CALL getVel(idstar,v)
       call getepoch(idstar,epoch1)	
       vd = 0.d0
       do k = 1,3
          vd = vd + v(k)**2
       enddo
       if(vd.gt.0.d0) vd = SQRT(vd)
*
       CALL getIsec(idbin,idstar)
       is = idstar
       CALL getType(idstar,kw)
       CALL getLabel(kw,label2,iend2)
       CALL getMass(idstar,m2)
       m2 = MIN(m2,999.d0)
       CALL getRadius(idstar,r2)
       call getepoch(idstar,epoch2)
       r2 = MIN(r2,9999.d0)
*
       CALL getBtype(idbin,kw)
       CALL getLabelb(kw,label3,iend3)
       CALL getSemi(idbin,a)
       a = MIN(a,9.999d+07)
       CALL getEcc(idbin,ecc)
       ecc = MIN(ecc,99.d0)
       vd = MIN(vd,99999.d0)
*
       WRITE(6,600)time
       WRITE(6,601)label3(1:iend3),label1(1:iend1),label2(1:iend2)
cThis is arranged so that the component of lowest id is printed first
       if (ip.lt.is) then
          WRITE(6,602)m1,r1
          WRITE(6,603)m2,r2
       elseif (ip.gt.is) then
          WRITE(6,602)m2,r2
          WRITE(6,603)m1,r1
       else
          write (6,*) 'ip = is 1, stopping'
          stop
       endif
       WRITE(6,604)a,ecc,vd
       print*,'epochs ',epoch1,epoch2
*
  600  FORMAT(' status at time = ',f12.4)
  601  FORMAT(a,' (',a,', ',a,')')
  602  FORMAT(' M = ',f12.8,' R = ',f13.8)
  603  FORMAT(' m = ',f12.8,' r = ',f13.8)
  604  FORMAT(' a = ',f14.4,' e = ',f7.3,' v = ',f8.1)
*
       RETURN
       END
***
       SUBROUTINE initStar(in,idstar,m_init)
       implicit none
       integer in,idstar,k,kw
       real*8 m_init,x(3),v(3),zini
       common /zset/ zini
*
       if(in.eq.1)then
          CALL initPar
          in = 2
       endif
*
       CALL getType(idstar,kw)
       if(kw.ge.0)then
          WRITE(6,*)' WARNING: id already in use ',idstar
       endif
*
       CALL setType(idstar,1)
       CALL setAge(idstar,0.d0)
       CALL setEpoch(idstar,0.d0)
       CALL setMStime(idstar,1.0d+10)
       CALL setTstep(idstar,0.d0)
       CALL setMass0(idstar,m_init)
       CALL setMass(idstar,m_init)
       CALL setMc(idstar,0.d0)
       CALL setMenv(idstar,0.d0)
       CALL setRadius(idstar,m_init)
       CALL setRc(idstar,0.d0)
       CALL setRenv(idstar,0.d0)
       CALL setLum(idstar,m_init)
       CALL setSpin(idstar,0.d0)
       CALL setRl(idstar,1.0d+10)
*
* Assume solar abundance for now. 
* (Note that chemical evolution is not yet accounted for).
*
       CALL setZ0(idstar,zini)
       CALL setY0(idstar,0.28d0)
       CALL setZ(idstar,zini)
       CALL setY(idstar,0.28d0)
*
       if(idstar.lt.10) print*,'idstar,zini =',idstar,zini
*                            
c       CALL setZ0(idstar,0.002d0)
c       CALL setY0(idstar,0.28d0)
c       CALL setZ(idstar,0.002d0)
c       CALL setY(idstar,0.28d0)
*
       do k = 1,3
          x(k) = 0.d0
          v(k) = 0.d0
       enddo
       CALL setPos(idstar,x)
       CALL setVel(idstar,v)
*
* Obtain initial values for stellar parameters. 
*
       CALL evStar(idstar,0.d0)
*
       RETURN
       END
***
       SUBROUTINE initBinary(in,idbin,a,e,m1,m2)
       implicit none
       integer in,idbin,idstar
       real*8 a,e,tb0,m1,m2,mx,aursun,yeardy
       COMMON /PARAMS/ aursun,yeardy
*
       if(m2.gt.m1)then
          mx = m1
          m1 = m2
          m2 = mx
       endif
*
       idstar = idbin
       CALL initStar(in,idstar,m1)
       CALL setIprim(idbin,idstar)
       idstar = idbin + 1
       CALL initStar(in,idstar,m2)
       CALL setIsec(idbin,idstar)
*
       CALL setSemi(idbin,a)
       CALL setEcc(idbin,e)
       tb0 = (a/aursun)*SQRT(a/(aursun*(m1+m2)))*yeardy
       CALL setTb(idbin,tb0)
       CALL setBtype(idbin,1)
*
       CALL evBinary(idbin,0.d0)
*
       RETURN
       END
***
       SUBROUTINE getLabel(kw,labelx,i2)
       implicit none
       include "interface_bse.h"
       integer kw,i1,i2
       character*20 labelx
       labelx = label(kw)
       CALL strip(labelx,i1,i2)
*     i2 = i2 + 3
*     if(kw.ge.10) i2 = i2 + 1
       RETURN
       END
***
       SUBROUTINE getLabelb(kw,labelx,i2)
       implicit none
       include "interface_bse.h"
       integer kw,i1,i2
       character*20 labelx
       labelx = labelb(kw)
       CALL strip(labelx,i1,i2)
       RETURN
       END
***
       SUBROUTINE strip(word,i1,i2)
       implicit none
       integer i,i1,i2,n
       parameter(n=80)
       character*(n) word
*
       do i = 1,n
          i1 = i
          if(word(i:i).ne.' ') goto 1
       enddo
  1    do i = i1,n
          i2 = i
          if(word(i:i).eq.' ') goto 2
       enddo
  2    i2 = i2 - 1
*
       RETURN
       END
***
       SUBROUTINE ssupdatetime(idstar,time)
       implicit none
       integer idstar
       real*8 time,age,dt
       CALL getAge(idstar,age)
       CALL getTstep(idstar,dt)
       time = age + dt
       RETURN
       END
***
       SUBROUTINE bsupdatetime(idbin,time)
       implicit none
       integer idbin,idstar
       real*8 time
       CALL getIprim(idbin,idstar)
       CALL ssupdatetime(idstar,time)
       RETURN
       END
***
       SUBROUTINE binaryexists(idbin,iexist)
       implicit none
       integer idbin,iexist,kw
       CALL getBtype(idbin,kw)
       iexist = 1
       if(kw.gt.3) iexist = 0
       RETURN
       END
***
       SUBROUTINE setAge(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       age(idstar) = time
       RETURN
       END
***
       SUBROUTINE getAge(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = age(idstar)
       RETURN
       END
***
       SUBROUTINE setMStime(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       ms_lifetime(idstar) = time
       RETURN
       END
***
       SUBROUTINE getMStime(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = ms_lifetime(idstar)
       RETURN
       END
***
       SUBROUTINE setTstep(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       standard_timestep(idstar) = time
       RETURN
       END
***
       SUBROUTINE getTstep(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = standard_timestep(idstar)
       RETURN
       END
***
       SUBROUTINE setEpoch(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       epoch(idstar) = time
       RETURN
       END
***
       SUBROUTINE getEpoch(idstar,time)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 time
       time = epoch(idstar)
       RETURN
       END
***
       SUBROUTINE setMass(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       zmass(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMass(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = zmass(idstar)
       RETURN
       END
***
       SUBROUTINE setMass0(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       zmass0(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMass0(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = zmass0(idstar)
       RETURN
       END
***
       SUBROUTINE setRadius(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       radius(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRadius(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = radius(idstar)
       RETURN
       END
***
       SUBROUTINE setLum(idstar,lum)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 lum
       zlum(idstar) = lum
       RETURN
       END
***
       SUBROUTINE getLum(idstar,lum)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 lum
       lum = zlum(idstar)
       RETURN
       END
***
       SUBROUTINE setY0(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       yinit(idstar) = y
       RETURN
       END
***
       SUBROUTINE getY0(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       y = yinit(idstar)
       RETURN
       END
***
       SUBROUTINE setZ0(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       zinit(idstar) = z
       RETURN
       END
***
       SUBROUTINE getZ0(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       z = zinit(idstar)
       RETURN
       END
***
       SUBROUTINE setY(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       ycurr(idstar) = y
       RETURN
       END
***
       SUBROUTINE getY(idstar,y)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 y
       y = ycurr(idstar)
       RETURN
       END
***
       SUBROUTINE setZ(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       zcurr(idstar) = z
       RETURN
       END
***
       SUBROUTINE getZ(idstar,z)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 z
       z = zcurr(idstar)
       RETURN
       END
***
       SUBROUTINE setMc(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       zmassc(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMc(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = zmassc(idstar)
       RETURN
       END
***
       SUBROUTINE setRc(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       radc(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRc(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = radc(idstar)
       RETURN
       END
***
       SUBROUTINE setMenv(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       menv(idstar) = mass
       RETURN
       END
***
       SUBROUTINE getMenv(idstar,mass)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 mass
       mass = menv(idstar)
       RETURN
       END
***
       SUBROUTINE setRenv(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       renv(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRenv(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = renv(idstar)
       RETURN
       END
***
       SUBROUTINE setSpin(idstar,sp)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 sp
       spin(idstar) = sp
       RETURN
       END
***
       SUBROUTINE getSpin(idstar,sp)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 sp
       sp = spin(idstar)
       RETURN
       END
***
       SUBROUTINE setRl(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rlof(idstar) = rad
       RETURN
       END
***
       SUBROUTINE getRl(idstar,rad)
       implicit none
       include "interface_bse.h"
       integer idstar
       real*8 rad
       rad = rlof(idstar)
       RETURN
       END
***
       SUBROUTINE setPos(idstar,x)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 x(3)
       do k = 1,3
          xstar(k,idstar) = x(k)
       enddo
       RETURN
       END
***
       SUBROUTINE getPos(idstar,x)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 x(3)
       do k = 1,3
          x(k) = xstar(k,idstar)
       enddo
       RETURN
       END
***
       SUBROUTINE setVel(idstar,v)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 v(3)
       do k = 1,3
          vstar(k,idstar) = v(k)
       enddo
       RETURN
       END
***
       SUBROUTINE getVel(idstar,v)
       implicit none
       include "interface_bse.h"
       integer idstar,k
       real*8 v(3)
       do k = 1,3
          v(k) = vstar(k,idstar)
       enddo
       RETURN
       END
***
       SUBROUTINE setType(idstar,ktype1)
       implicit none
       include "interface_bse.h"
       integer idstar,ktype1
       kstar(idstar) = ktype1
       RETURN
       END
***
       SUBROUTINE getType(idstar,ktype1)
       implicit none
       include "interface_bse.h"
       integer idstar,ktype1
       ktype1 = kstar(idstar)
       RETURN
       END
***
       SUBROUTINE setBtype(idbin,ktype1)
       implicit none
       include "interface_bse.h"
       integer idbin,ktype1
       kstarb(idbin) = ktype1
       RETURN
       END
***
       SUBROUTINE getBtype(idbin,ktype1)
       implicit none
       include "interface_bse.h"
       integer idbin,ktype1
       ktype1 = kstarb(idbin)
       RETURN
       END
***
       SUBROUTINE setIprim(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       iprim(idbin) = idstar
       RETURN
       END
***
       SUBROUTINE getIprim(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       idstar = iprim(idbin)
       RETURN
       END
***
       SUBROUTINE setIsec(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       isec(idbin) = idstar
       RETURN
       END
***
       SUBROUTINE getIsec(idbin,idstar)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       idstar = isec(idbin)
       RETURN
       END
***
       SUBROUTINE getMprim(idbin,mass)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       real*8 mass
       idstar = iprim(idbin)
       mass = zmass(idstar)
       RETURN
       END
***
       SUBROUTINE getMsec(idbin,mass)
       implicit none
       include "interface_bse.h"
       integer idbin,idstar
       real*8 mass
       idstar = isec(idbin)
       mass = zmass(idstar)
       RETURN
       END
***
       SUBROUTINE setSemi(idbin,arsun)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 arsun
       sep(idbin) = arsun
       RETURN
       END
***
       SUBROUTINE getSemi(idbin,arsun)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 arsun
       arsun = sep(idbin)
       RETURN
       END
***
       SUBROUTINE setEcc(idbin,e)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 e
       ecc(idbin) = e
       RETURN
       END
***
       SUBROUTINE getEcc(idbin,e)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 e
       e = ecc(idbin)
       RETURN
       END
***
       SUBROUTINE setTb(idbin,tbday)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 tbday
       tb(idbin) = tbday
       RETURN
       END
***
       SUBROUTINE getTb(idbin,tbday)
       implicit none
       include "interface_bse.h"
       integer idbin
       real*8 tbday
       tbday = tb(idbin)
       RETURN
       END
***
       SUBROUTINE evStar(idstar,time)
       implicit none
       integer idstar,kw
       real*8 time,tphys,tphysf,dmmax,drmax
       COMMON /TSTEPC/ dmmax,drmax
       real*8 z,zpars(20)
       real*8 mass,mt,mc,me
       real*8 rad,rc,re
       real*8 epch,tm,lum,ospin,vs(3)
*
       CALL getAge(idstar,tphys)
       tphysf = time
       if(tphysf.ge.tphys.or.tphys.le.1.0d-10)then
*
          CALL getZ0(idstar,z)
          CALL zcnsts(z,zpars)
*
          CALL getType(idstar,kw)
          CALL getEpoch(idstar,epch)
          CALL getMStime(idstar,tm)
          CALL getMass0(idstar,mass)
          CALL getMass(idstar,mt)
          CALL getMc(idstar,mc)
          CALL getMenv(idstar,me)
          CALL getRadius(idstar,rad)
          CALL getRc(idstar,rc)
          CALL getRenv(idstar,re)
          CALL getLum(idstar,lum)
          CALL getSpin(idstar,ospin)
          CALL getVel(idstar,vs)
c          print*,rad, 'before evolvlb'
*
          CALL evolv1b(kw,mass,mt,rad,lum,mc,rc,me,re,ospin,
     &                epch,tm,tphys,tphysf,z,zpars,dmmax,drmax,vs)
*
c         if (kw.ne.13) then
c          print*,rad, 'after evolvlb'
          CALL setType(idstar,kw)
          CALL setAge(idstar,tphys)
          CALL setEpoch(idstar,epch)
          CALL setMStime(idstar,tm)
          CALL setTstep(idstar,tphysf-tphys)
          CALL setMass0(idstar,mass)
          CALL setMass(idstar,mt)
          CALL setMc(idstar,mc)
          CALL setMenv(idstar,me)
          CALL setRadius(idstar,rad)
          CALL setRc(idstar,rc)
          CALL setRenv(idstar,re)
          CALL setLum(idstar,lum)
          CALL setSpin(idstar,ospin)
          CALL setVel(idstar,vs)
c         else
          if (kw.eq.13) then
              write (6,*) 'interface ',idstar,kw,tphys,epch,tm,tphysf,
     &            mass,
     &           mt,mc,me,rad,rc,re,lum,ospin,vs
          endif
*
       endif
*
       RETURN
       END
***
       SUBROUTINE evBinary(idbin,time)
       implicit none
       integer idbin,idstar,id1,id2
       integer kw1,kw2,kw(2),k
       real*8 time,tphys,tphysf,dmmax,drmax,aursun,yeardy
       COMMON /TSTEPC/ dmmax,drmax
       COMMON /PARAMS/ aursun,yeardy
       real*8 z,zpars(20)
       real*8 mass(2),mt(2),mc(2),me(2)
       real*8 rad(2),rc(2),re(2),rol(2)
       real*8 epch(2),tm(2),lum(2),ospin(2),dmdt(2)
       real*8 semi,tb,ecc,vs(3)
       logical iprint
*
       CALL getIprim(idbin,idstar)
       CALL getAge(idstar,tphys)
       tphysf = time
       iprint = .false.
c      iprint = .true.
*
       if(tphysf.ge.tphys.or.tphys.le.1.0d-10)then
*
          CALL getZ0(idstar,z)
          CALL zcnsts(z,zpars)
          CALL getVel(idstar,vs)
*
          id1 = idstar
          do k = 1,2
             if(k.eq.2) CALL getIsec(idbin,idstar)
             CALL getType(idstar,kw(k))
             CALL getEpoch(idstar,epch(k))
             CALL getMStime(idstar,tm(k))
             CALL getMass0(idstar,mass(k))
             CALL getMass(idstar,mt(k))
             CALL getMc(idstar,mc(k))
             CALL getMenv(idstar,me(k))
             CALL getRadius(idstar,rad(k))
             CALL getRc(idstar,rc(k))
             CALL getRenv(idstar,re(k))
             CALL getLum(idstar,lum(k))
             CALL getSpin(idstar,ospin(k))
             CALL getRl(idstar,rol(k))
          enddo
          id2 = idstar
          kw1 = kw(1)
          kw2 = kw(2)
*
          CALL getSemi(idbin,semi)
          CALL getEcc(idbin,ecc)
          tb = (semi/aursun)*SQRT(semi/(aursun*(mt(1)+mt(2))))*yeardy
*
          CALL evolv2b(kw,mass,mt,rad,lum,mc,rc,me,re,ospin,rol,
     &               dmdt,epch,tm,tphys,tphysf,z,zpars,dmmax,drmax,
     &                tb,ecc,vs)
*
c          print*,'evbin period',tb
          if(rad(2)/rol(2).gt.rad(1)/rol(1))then
             CALL setIprim(idbin,id2)
             CALL setIsec(idbin,id1)
             iprint = .true.
          elseif(kw1.ne.kw(1).or.kw2.ne.kw(2))then
             iprint = .true.
          endif
*
          idstar = id1
          do k = 1,2
             if(k.eq.2) idstar = id2
             CALL setType(idstar,kw(k))
             CALL setAge(idstar,tphys)
             CALL setEpoch(idstar,epch(k))
             CALL setMStime(idstar,tm(k))
             CALL setTstep(idstar,tphysf-tphys)
             CALL setMass0(idstar,mass(k))
             CALL setMass(idstar,mt(k))
             CALL setMc(idstar,mc(k))
             CALL setMenv(idstar,me(k))
             CALL setRadius(idstar,rad(k))
             CALL setRc(idstar,rc(k))
             CALL setRenv(idstar,re(k))
             CALL setLum(idstar,lum(k))
             CALL setSpin(idstar,ospin(k))
             CALL setRl(idstar,rol(k))
             CALL setVel(idstar,vs)
          enddo
*
* Set binary type and semi-major axis. 
*
          CALL getBtype(idbin,id1)
          if(tb.le.0.d0)then
             if(kw(1).lt.15.and.kw(2).lt.15)then
                id2 = 6
             elseif(kw(1).lt.15.or.kw(2).lt.15)then
                id2 = 4
             else
                id2 = 5
             endif
             semi = -1.d0
          else
             if(rad(1)/rol(1).lt.1.d0.and.rad(2)/rol(2).lt.1.d0)then
                id2 = 1
             else
                if(rad(1)/rol(1).ge.1.d0.and.rad(2)/rol(2).ge.1.d0)then
                   id2 = 3
                else
                   id2 = 2
                endif
             endif
             semi = aursun*((mt(1)+mt(2))*(tb/yeardy)**2)**(1.d0/3.d0)
          endif
          CALL setBtype(idbin,id2)
          CALL setSemi(idbin,semi)
          CALL setEcc(idbin,ecc)
          CALL setTb(idbin,tb)
          if(id2.ne.id1) iprint = .true.
*
          if(iprint)then
             CALL printBinary(idbin)
             call print_roche_data(idbin)
          endif
*
       endif
*
       RETURN
       END
*
************************************************************************
cAdditions for v3.1:
cAlso, an extra call to print_roche_data was also added to evBinary
       subroutine init_binaries (in, id, aRsun, e, m1, m2)
       implicit none
       integer in,id
       double precision aRsun,e,m1,m2,a

C        if (in.ne.1) then
C           write (6,*) 'in != 1 in init_binaries not implemented,',
C      &        ' stopping'
C           stop
C        endif
       a = aRsun
       call initBinary(in,id,a,e,m1,m2)
       open (7,file='BSE.data')
       write (7,*) '      TIME    M1     M2  KW1 KW2    SEP    ECC R1/',
     &     'ROL1 R2/ROL2  TYPE'
       return
       end

       subroutine out_binary(id)
       implicit none
       integer id
       call printBinary(id)
       call print_roche_data(id)
       return
       end

       subroutine get_bs_updatetime(id, updt)
       implicit none
       double precision updt
       integer id
       call bsupdatetime(id,updt)
       return
       end

       subroutine morph_binary(idb, id3, a_factor, e, outcome)
       implicit none
       integer idb,id3,outcome
       double precision a_factor,e
       if (outcome.eq.0) then
          a_factor = 1
          call getEcc(idb,e)
          return
       endif
       if (outcome.eq.1) then
          call exchange_hiid(idb,id3,a_factor,e)
          return
       endif
       if (outcome.eq.2) then
          call exchange_loid(idb,id3,a_factor,e)
          return
       endif
       write (6,*) 'outcome ',outcome,
     &     ' unexpected in morph_binary, stopping'
       stop
       end

       subroutine out_star(ids)
       implicit none
       integer ids,kw,iend2
       character*20 label2
       call getType(ids,kw)
       CALL getLabel(kw,label2,iend2)
       write (6,*) 'type (',label2,'), '
       return
       end

       subroutine ev_binary(id,tMyr)
       implicit none
       integer id
       real*8 tMyr,updatetime
c 10    continue
c       call bsupdatetime(id,updatetime)
c       if (updatetime.ge.tMyr) then
          call evBinary(id,tMyr)
c          print*,'binary',id,tMyr
c       else
c          call evBinary(id,updatetime)
c          print*,'binary',id,updatetime
c          goto 10
c       endif
       return
       end

       subroutine binary_exists(id, iexist)
       implicit none
       integer id,iexist
       call binaryexists(id,iexist)
       return
       end

       subroutine get_sma(id, aRsun)
       implicit none
       integer id
       real*8 aRsun
       call getSemi(id,aRsun)
       return
       end

       subroutine get_ecc(id, e)
       implicit none
       integer id
       real*8 e
       call getEcc(id,e)
       return
       end

       subroutine get_loid_mass(id, m1)
       implicit none
       integer id,ip,is
       real*8 m1
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.lt.is) then
          call getMprim(id,m1)
       elseif (ip.gt.is) then
          call getMsec(id,m1)
       else
          write (6,*) 'ip = is 2, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine get_hiid_mass(id, m2)
       implicit none
       integer id,ip,is
       real*8 m2
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.gt.is) then
          call getMprim(id,m2)
       elseif (ip.lt.is) then
          call getMsec(id,m2)
       else
          write (6,*) 'ip = is 3, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine out_scatter(id,id3,tMyr,write_text)
       implicit none
       integer id,id3
       real*8 tMyr
       logical write_text
       return
       end

       subroutine init_stars(in, id3, m3)
       implicit none
       integer in,id3
       real*8 m3
       call initStar(in,id3,m3)
       return
       end

       subroutine ev_star(id3, tmyr)
       implicit none
       integer id3
       real*8 tmyr
       call evStar(id3,tmyr)
       return
       end

       subroutine get_mass(id3, m3)
       implicit none
       integer id3
       real*8 m3
       call getMass(id3,m3)
       return
       end

       subroutine get_radius(id3, r3)
       implicit none
       integer id3
       real*8 r3
       call getRadius(id3,r3)
       return
       end

       subroutine exchange_loid(idb,id3,a_factor,e)
       implicit none
       integer idb,id3,ip,is
       real*8 a_factor,e
       CALL getIprim(idb,ip)
       CALL getIsec(idb,is)
       if (ip.lt.is) then
          call setIprim(idb,id3)
       elseif (ip.gt.is) then
          call setIsec(idb,id3)
       else
          write (6,*) 'ip = is 4, stopping...'
          stop
       endif
       call update_binary_parameters(idb,a_factor,e)
       return
       end

       subroutine exchange_hiid(idb,id3,a_factor,e)
       implicit none
       integer idb,id3,ip,is
       real*8 a_factor,e
       CALL getIprim(idb,ip)
       CALL getIsec(idb,is)
       if (ip.lt.is) then
          call setIsec(idb,id3)
       elseif (ip.gt.is) then
          call setIprim(idb,id3)
       else
          write (6,*) 'ip = is 5, stopping...'
          stop
       endif
       call update_binary_parameters(idb,a_factor,e)
       return
       end

       subroutine update_binary_parameters(idb,a_factor,e)
       implicit none
       integer idb
       real*8 a_factor,e,arsun
       call getSemi(idb,arsun)
       arsun = arsun*a_factor
       call setSemi(idb,arsun)
       call setEcc(idb,e)
       return
       end

       subroutine print_roche_data(idbin)
       implicit none
       integer idbin,idstar
       integer k,kw,iend1,iend2,iend3
       real*8 time
       real*8 m1,r1,m2,r2,a,ecc,v(3),vd
       real*8 rl1,rl2
       integer kw1,kw2,id1,id2,id10,id20,kw10,kw20
       character*20 label1,label2,label3,type,label30
       logical first,swap
       data first/.true./
       save id10,id20,kw10,kw20,swap,label30
*
       CALL getIprim(idbin,idstar)
       id1 = idstar
       CALL getType(idstar,kw)
       kw1 = kw
       call getRl(idstar,rl1)
       CALL getLabel(kw,label1,iend1)
       CALL getMass(idstar,m1)
       m1 = MIN(m1,999.d0)
       CALL getRadius(idstar,r1)
       r1 = MIN(r1,9999.d0)
*
       CALL getAge(idstar,time)
       CALL getVel(idstar,v)
       vd = 0.d0
       do k = 1,3
          vd = vd + v(k)**2
       enddo
       if(vd.gt.0.d0) vd = SQRT(vd)
*
       CALL getIsec(idbin,idstar)
       id2 = idstar
       CALL getType(idstar,kw)
       kw2 = kw
       call getRl(idstar,rl2)
       CALL getLabel(kw,label2,iend2)
       CALL getMass(idstar,m2)
       m2 = MIN(m2,999.d0)
       CALL getRadius(idstar,r2)
       r2 = MIN(r2,9999.d0)
*
       CALL getBtype(idbin,kw)
       CALL getLabelb(kw,label3,iend3)
       CALL getSemi(idbin,a)
       a = MIN(a,9.999d+07)
       CALL getEcc(idbin,ecc)
       ecc = MIN(ecc,99.d0)
       vd = MIN(vd,99999.d0)
*
       if (first) then
          first = .false.
          type = 'INITIAL'
          swap = .false.
       else
          if (id1.eq.id10) then
             if (id2.eq.id20) then
cBoth components still in place; could be kw change or change of binary type 
c(e.g. no longer detached) or end of run
                if (kw1.ne.kw10.or.kw2.ne.kw20) then
                   type = 'KW CHNGE'
                else
                   if (label3.ne.label30) then
                      type = label3
                   else
                      type = 'MAX TIME'
                   endif
                endif
             else
cStar 2 must have been exchanged
                type = 'EXCHANGE'
             endif
          else
cStar 1 no longer in place
             if (id1.eq.id20) then
c..but star 1 is still in the binary; swop printing order
                swap = .not.swap
                if (id2.eq.id10) then
cStars 1 and 2 are interchanged
                   if (kw1.ne.kw20.or.kw2.ne.kw10) then
                      type = 'KW CHNGE'
                   else
                      if (label3.ne.label30) then
                         type = label3
                      else
cThis should not happen
                         type = 'ERROR 1'
                      endif
                   endif
                else
cStar 1 has swapped, but star 2 is no longer in the system
                   type = 'EXCHANGE'
                endif
             else
cStar 1 is no longer in the sytem...
                if (id2.eq.id20) then
c...but star 2 is still in place
                   type = 'EXCHANGE'
                else if (id2.eq.id10) then
c...but star 2 has been swapped
                   type = 'EXCHANGE'
                   swap = .not.swap
                else
cNeither star surviving
                   type = 'ERROR 2'
                endif
             endif
          endif
       endif
       if (.not.swap) then
          WRITE(7,100)time,m1,m2,kw1,kw2,a,ecc,r1/rl1,r2/rl2,type,label3
       else
          WRITE(7,100)time,m2,m1,kw2,kw1,a,ecc,r2/rl2,r1/rl1,type,label3
       endif
  100  format (f11.4,2f7.3,2i3,f10.3,f6.2,2f7.3,1x,2a10)
*
       id10 = id1
       id20 = id2
       kw10 = kw1
       kw20 = kw2
       label30 = label3
c
       RETURN
       END
c===========================================================================
cNew wrappers added Edinburgh 1/5/6
cThe following two corrected 13/11/6 following Mirek's email of 12/11/6

       subroutine get_loid(id,k)
       implicit none
       integer id,k
       integer ip,is
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.lt.is) then
          k = ip
       elseif (ip.gt.is) then
          k = is
       else
          write (6,*) 'ip = is 2, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine get_hiid(id,k)
       implicit none
       integer id,k
       integer ip,is
       CALL getIprim(id,ip)
       CALL getIsec(id,is)
       if (ip.lt.is) then
          k = is
       elseif (ip.gt.is) then
          k = ip
       else
          write (6,*) 'ip = is 3, stopping... id.ip,is =',id,ip,is
          stop
       endif
       return
       end

       subroutine get_ss_type(i,k)
       implicit none
       integer i,k
       call getType(i,k)
       return
       end

       subroutine get_bs_type(i,k)
       implicit none
       integer i,k
       call getBtype(i,k)
       return
       end

       subroutine get_ss_updatetime(i,t)
       implicit none
       integer i
       double precision t
       call ssupdatetime(i,t)
       return
       end
c===========================================================================
cNew subroutine added Edinburgh 4/5/6
cModelled on initBinary
       subroutine create_binary(idbin,id1,id2,sma,e,time)
       implicit none
c      include "interface_bse.h"
       double precision mx,m1,m2,sma,t1,t2,tb0,e
       integer idbin,idstar,id1,id2,idx,ix1,ix2
       integer kw1,kw2,kw(2),k
       real*8 time,tphys,tphysf,dmmax,drmax,aursun,yeardy
       COMMON /TSTEPC/ dmmax,drmax
       COMMON /PARAMS/ aursun,yeardy
       real*8 z,zpars(20)
       real*8 mass(2),mt(2),mc(2),me(2)
       real*8 rad(2),rc(2),re(2),rol(2)
       real*8 epch(2),tm(2),lum(2),ospin(2),dmdt(2)
       real*8 semi,tb,ecc,vs(3)
       logical iprint
c      
       call getAge(id1,t1)
       call getAge(id2,t2)
cPut the star of higher mass first
       call getMass(id1,m1)
       call getMass(id2,m2)
cCheck that both stars have been evolved to the same time
       if (t1.ne.t2) then
          print*,'create_binary: stars have different ages',id1,id2,  
     &            idbin,t1,t2,m1,m2,sma,e,time
c          stop
       endif
                                
       if(m2.gt.m1)then
          mx = m1
          m1 = m2
          m2 = mx
          idx = id1
          id1 = id2
          id2 = idx
       endif
*
       print*,'id1,id2,m1,m2,idbin =',id1,id2,m1,m2,idbin
       CALL setIprim(idbin,id1)
       CALL setIsec(idbin,id2)
       CALL setSemi(idbin,sma)
       CALL setEcc(idbin,e)
       tb0 = (sma/aursun)*SQRT(sma/(aursun*(m1+m2)))*yeardy
       CALL setTb(idbin,tb0)
cIt would be nice to distinguish primordial from three-body binaries, but still....
       CALL setBtype(idbin,1)
cNow we have to set all the parameters that would be needed in a call to evBinary
       tphys = t1
       tphysf = t1
       idstar = id1
       CALL getZ0(idstar,z)
       CALL zcnsts(z,zpars)
       CALL getVel(idstar,vs)
*
       iprint = .false.
       do k = 1,2
          if(k.eq.2) CALL getIsec(idbin,idstar)
          CALL getType(idstar,kw(k))
          CALL getEpoch(idstar,epch(k))
          CALL getMStime(idstar,tm(k))
          CALL getMass0(idstar,mass(k))
          CALL getMass(idstar,mt(k))
          CALL getMc(idstar,mc(k))
          CALL getMenv(idstar,me(k))
          CALL getRadius(idstar,rad(k))
          CALL getRc(idstar,rc(k))
          CALL getRenv(idstar,re(k))
          CALL getLum(idstar,lum(k))
          CALL getSpin(idstar,ospin(k))
          CALL getRl(idstar,rol(k))
       enddo
       kw1 = kw(1)
       kw2 = kw(2)
*
       CALL getSemi(idbin,semi)
       CALL getEcc(idbin,e)
       tb = (semi/aursun)*SQRT(semi/(aursun*(mt(1)+mt(2))))*yeardy
*
       CALL evolv2b(kw,mass,mt,rad,lum,mc,rc,me,re,ospin,rol,
     &    dmdt,epch,tm,tphys,tphysf,z,zpars,dmmax,drmax,
     &     tb,e,vs)
*
       if(rad(2)/rol(2).gt.rad(1)/rol(1))then
          CALL setIprim(idbin,id2)
          CALL setIsec(idbin,id1)
          iprint = .true.
       elseif(kw1.ne.kw(1).or.kw2.ne.kw(2))then
          iprint = .true.
       endif
*
       idstar = id1
       do k = 1,2
          if(k.eq.2) idstar = id2
          CALL setType(idstar,kw(k))
          CALL setAge(idstar,tphys)
          CALL setEpoch(idstar,epch(k))
          CALL setMStime(idstar,tm(k))
          CALL setTstep(idstar,tphysf-tphys)
          CALL setMass0(idstar,mass(k))
          CALL setMass(idstar,mt(k))
          CALL setMc(idstar,mc(k))
          CALL setMenv(idstar,me(k))
          CALL setRadius(idstar,rad(k))
          CALL setRc(idstar,rc(k))
          CALL setRenv(idstar,re(k))
          CALL setLum(idstar,lum(k))
          CALL setSpin(idstar,ospin(k))
          CALL setRl(idstar,rol(k))
          CALL setVel(idstar,vs)
       enddo
*
* Set binary type and semi-major axis. 
*
       CALL getBtype(idbin,ix1)
       if(tb.le.0.d0)then
          if(kw(1).lt.15.and.kw(2).lt.15)then
             ix2 = 6
          elseif(kw(1).lt.15.or.kw(2).lt.15)then
             ix2 = 4
          else
             ix2 = 5
          endif
          semi = -1.d0
       else
          if(rad(1)/rol(1).lt.1.d0.and.rad(2)/rol(2).lt.1.d0)then
             ix2 = 1
          else
             if(rad(1)/rol(1).ge.1.d0.and.rad(2)/rol(2).ge.1.d0)then
                ix2 = 3
             else
                ix2 = 2
             endif
          endif
          semi = aursun*((mt(1)+mt(2))*(tb/yeardy)**2)**(1.d0/3.d0)
       endif
       CALL setBtype(idbin,ix2)
       CALL setSemi(idbin,semi)
       CALL setEcc(idbin,e)
       CALL setTb(idbin,tb)
       if(id2.ne.id1) iprint = .true.
*     
       if(iprint)then
          CALL printBinary(idbin)
          call print_roche_data(idbin)
       endif
*
       return
       end
*
*
      SUBROUTINE collStars(in,id1,id2,time)
      implicit none
*
* Developed by J. Hurley (May 1, 2006) for use 
* with the SSE/BSE McScatter interface. 
*
* Determines collision product of Star 1 and Star 2 
* based on Hurley et al. 2002, MNRAS, 329, 897. 
* Assumes that product is to be placed in position 
* occupied by Star 1. 
* Any cases that potentially lead to common-envelope 
* (subgiant/giant + any other star) are treated here 
* as mergers without mass-loss. If instead it is 
* determined that the stars form a close binary prior 
* to contact, then the binary can be evolved with the 
* routine evBinary which contains the common-envelope 
* treatment.  
*
      include "const_bse.h"
*
      integer oldkw
      integer in,id1,id2,idstar
      integer kw1,kw2,kwstar
      real*8 time
      real*8 m1,m2,mstar
      real*8 m01,m02,m0star
      real*8 mc1,mc2,mcstar
      real*8 aj1,aj2,epch,ajstar,f1,f2
      real*8 tms1,tms2,tms
      real*8 xs(3),vs(3)
      real*8 z,zpars(20),tscls(20),lums(10),gb(10),tn
      logical tzo
*
cSee email from jhurley 15/x/8
c      tzo = .true.
      tzo = .false.
*
*     tn is not deffined so is set to 0
*
      tn = 0.d0      
*
      CALL getType(id1,kw1)
      CALL getType(id2,kw2)
      print*,'id1,id2,kw1,kw2,time = ',id1,id2,kw1,kw2,time
      if(kw1.lt.0) kw1 = abs(kw1)
      if(kw2.lt.0) kw2 = abs(kw2)
      kwstar = ktype(kw1,kw2)
      print*,'id1,id2,kw1,kw2,kw,time = ',id1,id2,kw1,kw2,kwstar,time
      if(kwstar.ge.113) kwstar = kwstar - 100
*
      CALL getMass(id1,m1)
      CALL getMass(id2,m2)
      mstar = m1 + m2
      CALL getMass0(id1,m01)
      CALL getMass0(id2,m02)
      m0star = mstar
*
      CALL getMc(id1,mc1)
      CALL getMc(id2,mc2)
      mcstar = mc1 + mc2
*
      CALL getEpoch(id1,epch)
      aj1 = time - epch
      CALL getEpoch(id2,epch)
      aj2 = time - epch
      ajstar = 0.d0
      CALL getMStime(id1,tms1)
      CALL getMStime(id2,tms2)
      print*,'kw1,kw2,kwstar,m1,m2,m01,m02,mc1,mc2,aj1,aj2,time= ',
     &        kw1,kw2,kwstar,m1,m2,m01,m02,mc1,mc2,aj1,aj2,time
      print*,'mstar,m0star,mcstar = ',mstar,m0star,mcstar
*
      CALL getPos(id1,xs)
      CALL getVel(id1,vs)
*
      idstar = id1
      CALL setType(idstar,-1)
      print*,' before initS '
      call flush(6)
      CALL initStar(in,idstar,mstar)
      print*,' after initS '
      call flush(6)
      CALL getZ0(idstar,z)
      print*,'z of merger',z
      CALL zcnsts(z,zpars)
      CALL getMStime(idstar,tms)
      print*,'z, tms',z,tms
*
      if(kwstar.eq.1)then
         if(kw1.eq.7) kwstar = 7
         ajstar = 0.1d0*tms*(aj1*m1/tms1+aj2*m2/tms2)/mstar
      elseif(kwstar.eq.4)then
         if(kw1.le.1)then
            mcstar = m2
            ajstar = aj2/tms2
         else
            mcstar = m1
            ajstar = aj1/tms1
         endif
         print*,'a- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
         CALL gntage(mcstar,mstar,kwstar,zpars,m0star,ajstar)
         print*,'b- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
      elseif(kwstar.eq.7)then
         if(kw1.lt.10)then
            ajstar = (tms/mstar)*(aj1*m1/tms1)
         else
            ajstar = (tms/mstar)*(aj2*m2/tms2)
         endif
      elseif(kwstar.le.9)then
         print*,'c- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
         CALL gntage(mcstar,mstar,kwstar,zpars,m0star,ajstar)
         print*,'d- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
      elseif(kwstar.le.12)then
         if(kwstar.lt.12.and.mstar.ge.1.44)then
            mstar = 0.d0
            kwstar = 15
         endif
      elseif(kwstar.eq.13.or.kwstar.eq.14)then
*       Absorb all mass into NS/BH unless taking the 
*       unstable Thorne-Zytkow object option. 
         if(tzo)then
            mstar = 0.d0
            if(kw1.ge.10) mstar = mstar + m1
            if(kw2.ge.10) mstar = mstar + m2
            m0star = mstar
            mcstar = mstar
            print*,'TZO kw1,kw2,mstar,m1,m2=',kw1,kw2,mstar,m1,m2
         endif
         if(kwstar.eq.13.and.mstar.gt.mxns) kwstar = 14
      elseif(kwstar.eq.15)then
         mstar = 0.d0
      elseif(kwstar.gt.100)then
*       Common envelope cases.  
*       In the absence of CE treatment assume that 
*       merger proceeds with no mass lost from the system. 
         kwstar = kwstar - 100
         if(kwstar.eq.4.or.kwstar.eq.7)then
            f1 = mc1
            if(kw1.le.3.or.kw1.eq.10) f1 = 0.d0
            if(kw1.eq.7) f1 = m1*aj1/tms1
            if(kw1.eq.4.or.kw1.eq.5)then
        print*,'e- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               CALL star(kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars)
        print*,'f- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               f1 = f1*(aj1 - tscls(2))/(tscls(13) - tscls(2))
            endif
            f2 = mc2
            if(kw2.le.3.or.kw2.eq.10) f2 = 0.d0
            if(kw2.eq.7) f2 = m2*aj2/tms2
            if(kw2.eq.4.or.kw2.eq.5)then
        print*,'g- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               CALL star(kw2,m02,m2,tms2,tn,tscls,lums,gb,zpars)
        print*,'h- kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars =',
     & kw1,m01,m1,tms1,tn,tscls,lums,gb,zpars
               f2 = f2*(aj2 - tscls(2))/(tscls(13) - tscls(2))
            endif
            ajstar = (f1+f2)/mcstar
         endif
         if(kwstar.eq.7)then
            ajstar = ajstar*tms
cccc         elseif(kwstar.eq.7)then
         else
         print*,'i- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
            oldkw = kwstar
            CALL gntage(mcstar,mstar,kwstar,zpars,m0star,ajstar)
         print*,'j- kw1,kwa,mcstar,mstar,m0star=',kw1,kwstar,mcstar,
     &    mstar,m0star
cThe following two lines are really meant to deal with black hole outcomes
cThe previous version gave NaN stellar radii
            print*,'collision routine kw,ajstar',kwstar,ajstar,
     &        'changed to original value',oldkw,'and 0.d0'
c            kwstar = oldkw
c            ajstar = 0.d0
         endif
         print*,'kwstar,mstar,m0star,mcstar,ajstar = ',
     &           kwstar,mstar,m0star,mcstar,ajstar
      else
*       This should not be reached.
        print*,' *************** '
        kwstar = 1
      endif
*
      call getz0(idstar,z)
c      print*,z
      call getz(idstar,z)
c      print*,z
      epch = time - ajstar
      ajstar = MAX(0.d0,time-1.0d-14)
      CALL setAge(idstar,time)       
      CALL setEpoch(idstar,epch)
      CALL setType(idstar,kwstar)
      CALL setMass0(idstar,m0star)
      CALL setMass(idstar,mstar)
      CALL setMc(idstar,mcstar)
      CALL setPos(idstar,xs)
      CALL setVel(idstar,vs)
*
      call getz0(idstar,z)
      print*,z
      print*,'before evSTAR kwstar,time,epoch,mstar,m0star,mcstar =',
     &        kwstar,time,epch,mstar,m0star,mcstar
      CALL evStar(idstar,time)
      call getz0(idstar,z)
      print*,z
      call getMass(idstar,mstar)
      call getMass0(idstar,m0star)
      call getMc(idstar,mcstar)
      call getEpoch(idstar,epch)
      call getType(idstar,kwstar)
      print*,'after evStar kwstar,time,epoch,mstar,m0star,mcstar =',
     &        kwstar,time,epch,mstar,m0star,mcstar
           
      CALL printStar(idstar)
*
      RETURN
      END
*
***
*
* Developed by J. Hurley (August 29, 2006) for use
* with the SSE/BSE McScatter interface.
*
*   TURN  - an initial guess at the turn-off mass (set large if unsure)
*   TPHYS - current age in Myr
*   ZPARS - the metallicity parameter array set by an earlier call
*            to ZCNSTS.
*
***
      SUBROUTINE MTURN(TURN,TPHYS,ZPARS)
*
*
*       Current MS turn-off mass.
*       --------------------------------------
*
      IMPLICIT NONE
      INTEGER  I,II,IMAX
      PARAMETER(IMAX=30)
      REAL*8 TURN,TPHYS,ZPARS(20)
      REAL*8 TM,TURN2,DM,FMID,TACC
      PARAMETER(TACC=0.001D0)
      REAL*8 THOOKF,TBGBF
      EXTERNAL THOOKF,TBGBF
*
      TURN2 = 100.D0
      TM = MAX(ZPARS(8),THOOKF(TURN2))*TBGBF(TURN2)
      IF(TM.GT.TPHYS)THEN
         TURN = TURN2
         GOTO 40
      ENDIF
*
      II = 0
 25   TM = MAX(ZPARS(8),THOOKF(TURN))*TBGBF(TURN)
      IF(TM.GT.TPHYS)THEN
         IF(TPHYS.LE.0.D0.OR.TURN.GT.98.D0) GOTO 40 
         TURN = 2.D0*TURN
         II = II + 1
         GOTO 25
      ENDIF
      TURN2 = TURN
      DM = TURN
      DO 30 , I = 1,IMAX
         DM = 0.5D0*DM
         TURN = TURN2 - DM
         TM = MAX(ZPARS(8),THOOKF(TURN))*TBGBF(TURN)
         FMID = TM - TPHYS
         IF(FMID.LT.0.0) TURN2 = TURN
         IF(DM.LT.TACC.OR.ABS(FMID).LT.1.0D-14) GOTO 40
         IF(I.EQ.IMAX)THEN
            GOTO 40
         ENDIF
 30   CONTINUE
 40   CONTINUE
*
      RETURN
*
      END
***
        



***
      SUBROUTINE COMENV(M01,M1,MC1,AJ1,JSPIN1,KW1,
     &                  M02,M2,MC2,AJ2,JSPIN2,KW2,
     &                  ZPARS,ECC,SEP,JORB,COEL)
*
* Common Envelope Evolution.
*
*     Author : C. A. Tout
*     Date :   18th September 1996
*
*     Redone : J. R. Hurley
*     Date :   7th July 1998
*
      IMPLICIT NONE
*
      INTEGER KW1,KW2,KW
      INTEGER KTYPE(0:14,0:14)
      COMMON /TYPES/ KTYPE
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
*
      REAL*8 M01,M1,MC1,AJ1,JSPIN1,R1,L1,K21
      REAL*8 M02,M2,MC2,AJ2,JSPIN2,R2,L2,K22,MC22
      REAL*8 TSCLS1(20),TSCLS2(20),LUMS(10),GB(10),TM1,TM2,TN,ZPARS(20)
      REAL*8 EBINDI,EBINDF,EORBI,EORBF,ECIRC,SEPF,SEPL,MF,XX
      REAL*8 CONST,DELY,DERI,DELMF,MC3,FAGE1,FAGE2
      REAL*8 ECC,SEP,JORB,TB,OORB,OSPIN1,OSPIN2,TWOPI
      REAL*8 RC1,RC2,Q1,Q2,RL1,RL2,LAMB1,LAMB2
      REAL*8 MENV,RENV,MENVD,RZAMS,VS(3)
      REAL*8 AURSUN,K3,ALPHA1,LAMBDA
      PARAMETER (AURSUN = 214.95D0,K3 = 0.21D0) 
      COMMON /VALUE2/ ALPHA1,LAMBDA
      LOGICAL COEL
      REAL*8 CELAMF,RL,RZAMSF
      EXTERNAL CELAMF,RL,RZAMSF
*
* Common envelope evolution - entered only when KW1 = 2, 3, 4, 5, 6, 8 or 9.
*
* For simplicity energies are divided by -G.
*
      TWOPI = 2.D0*ACOS(-1.D0)
      COEL = .FALSE.
*
* Obtain the core masses and radii.
*
      KW = KW1
      CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
      CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &            R1,L1,KW1,MC1,RC1,MENV,RENV,K21)
      OSPIN1 = JSPIN1/(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
      MENVD = MENV/(M1-MC1)
      RZAMS = RZAMSF(M01)
      LAMB1 = CELAMF(KW,M01,L1,R1,RZAMS,MENVD,LAMBDA)
*     IF(KW1.NE.KW)THEN
*        WRITE(66,*)' COMENV TYPE CHANGE1 ',KW,KW1
*     ENDIF
      KW = KW2
      CALL star(KW2,M02,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
      CALL hrdiag(M02,AJ2,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS,
     &            R2,L2,KW2,MC2,RC2,MENV,RENV,K22)
      OSPIN2 = JSPIN2/(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
*     IF(KW2.NE.KW)THEN
*        WRITE(66,*)' COMENV TYPE CHANGE2 ',KW,KW2
*     ENDIF
*
* Calculate the binding energy of the giant envelope (multiplied by lambda).
*
      EBINDI = M1*(M1-MC1)/(LAMB1*R1)
*
* If the secondary star is also giant-like add its envelopes's energy.
*
      IF(KW2.GE.2.AND.KW2.LE.9.AND.KW2.NE.7)THEN
         MENVD = MENV/(M2-MC2)
         RZAMS = RZAMSF(M02)
         LAMB2 = CELAMF(KW,M02,L2,R2,RZAMS,MENVD,LAMBDA)
         EBINDI = EBINDI + M2*(M2-MC2)/(LAMB2*R2)
*
* Calculate the initial orbital energy
*
         EORBI = MC1*MC2/(2.D0*SEP)
      ELSE
         EORBI = MC1*M2/(2.D0*SEP)
      ENDIF
*
* Allow for an eccentric orbit.
*
      ECIRC = EORBI/(1.D0 - ECC*ECC)
*
* Calculate the final orbital energy without coalescence.
*
      EORBF = EORBI + EBINDI/ALPHA1
*
* If the secondary is on the main sequence see if it fills its Roche lobe.
*
      IF(KW2.LE.1.OR.KW2.EQ.7)THEN
         SEPF = MC1*M2/(2.D0*EORBF)
***
* Assume the energy generated by forcing the secondary to 
* co-rotate goes into the envelope (experimental). 
*        IF(CEFLAG.GT.5)THEN
*           TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(MC1+M2)))
*           OORB = TWOPI/TB
*           DELY = 0.5D0*M2*R2*R2*(OSPIN2*OSPIN2 - OORB*OORB)/3.91D+08
*           DELY = K22*DELY
*           EBINDI = MAX(0.D0,EBINDI - DELY)
*           EORBF = EORBI + EBINDI/ALPHA1
*           SEPF = MC1*M2/(2.D0*EORBF)
*        ENDIF
***
         Q1 = MC1/M2
         Q2 = 1.D0/Q1
         RL1 = RL(Q1)
         RL2 = RL(Q2)
         IF(RC1/RL1.GE.R2/RL2)THEN
*
* The helium core of a very massive star of type 4 may actually fill
* its Roche lobe in a wider orbit with a very low-mass secondary.
*
            IF(RC1.GT.RL1*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC1/RL1
            ENDIF
         ELSE
            IF(R2.GT.RL2*SEPF)THEN
               COEL = .TRUE.
               SEPL = R2/RL2
            ENDIF
         ENDIF
         IF(COEL)THEN
*
            KW = KTYPE(KW1,KW2) - 100
            MC3 = MC1
            IF(KW2.EQ.7.AND.KW.EQ.4) MC3 = MC3 + M2
*
* Coalescence - calculate final binding energy.
*
            EORBF = MAX(MC1*M2/(2.D0*SEPL),EORBI)
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI)
         ELSE
*
* Primary becomes a black hole, neutron star, white dwarf or helium star.
*
            MF = M1
            M1 = MC1
            CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &                  R1,L1,KW1,MC1,RC1,MENV,RENV,K21)
            IF(KW1.GE.13)THEN
               print*,'kick 1a- kw1,mf,m1,m2,ecc,sepf,jorb,vs=',
     &                 kw1,mf,m1,m2,ecc,sepf,jorb,vs(1),vs(2),vs(3)
               CALL kickv(KW1,MF,M1,M2,ECC,SEPF,JORB,VS)
               print*,'kick 1b- kw1,mf,m1,m2,ecc,sepf,jorb,vs=',
     &                 kw1,mf,m1,m2,ecc,sepf,jorb,vs(1),vs(2),vs(3)               
               IF(ECC.GT.1.D0) GOTO 30
            ENDIF
         ENDIF
      ELSE
*
* Degenerate or giant secondary. Check if the least massive core fills its
* Roche lobe.
*
         SEPF = MC1*MC2/(2.D0*EORBF)
***
* Assume the energy generated by forcing the secondary to 
* co-rotate goes into the envelope (experimental). 
*        IF(KW2.GE.10.AND.KW2.LE.14.AND.CEFLAG.GT.5)THEN
*           TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(MC1+MC2)))
*           OORB = TWOPI/TB
*           DELY = 0.5D0*M2*R2*R2*(OSPIN2*OSPIN2 - OORB*OORB)/3.91D+08
*           DELY = K3*DELY
*           EBINDI = MAX(0.D0,EBINDI - DELY)
*           EORBF = EORBI + EBINDI/ALPHA1
*           SEPF = MC1*MC2/(2.D0*EORBF)
*        ENDIF
***
         Q1 = MC1/MC2
         Q2 = 1.D0/Q1
         RL1 = RL(Q1)
         RL2 = RL(Q2)
         IF(RC1/RL1.GE.RC2/RL2)THEN
            IF(RC1.GT.RL1*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC1/RL1
            ENDIF
         ELSE
            IF(RC2.GT.RL2*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC2/RL2
            ENDIF
         ENDIF
*
         IF(COEL)THEN
*
* If the secondary was a neutron star or black hole the outcome
* is an unstable Thorne-Zytkow object that leaves only the core.
*
            SEPF = 0.D0
            IF(KW2.GE.13)THEN
               MC1 = MC2
               M1 = MC1
               MC2 = 0.D0
               M2 = 0.D0
               KW1 = KW2
               KW2 = 15
               AJ1 = 0.D0
*
* The envelope mass is not required in this case.
*
               GOTO 30
            ENDIF
*
            KW = KTYPE(KW1,KW2) - 100
            MC3 = MC1 + MC2
*
* Calculate the final envelope binding energy.
*
            EORBF = MAX(MC1*MC2/(2.D0*SEPL),EORBI)
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI)
*
* Check if we have the merging of two degenerate cores and if so
* then see if the resulting core will survive or change form.
*
            IF(KW1.EQ.6.AND.(KW2.EQ.6.OR.KW2.GE.11))THEN
               CALL dgcore(KW1,KW2,KW,MC1,MC2,MC3,EBINDF)
            ENDIF
            IF(KW1.LE.3.AND.M01.LE.ZPARS(2))THEN
               IF((KW2.GE.2.AND.KW2.LE.3.AND.M02.LE.ZPARS(2)).OR.
     &             KW2.EQ.10)THEN
                  CALL dgcore(KW1,KW2,KW,MC1,MC2,MC3,EBINDF)
                  IF(KW.GE.10)THEN
                     KW1 = KW
                     M1 = MC3
                     MC1 = MC3
                     IF(KW.LT.15) M01 = MC3
                     AJ1 = 0.D0
                     MC2 = 0.D0
                     M2 = 0.D0
                     KW2 = 15
                     GOTO 30
                  ENDIF
               ENDIF
            ENDIF
*
         ELSE
*
* The cores do not coalesce - assign the correct masses and ages.
*
            MF = M1
            M1 = MC1
            CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &                  R1,L1,KW1,MC1,RC1,MENV,RENV,K21)
            IF(KW1.GE.13)THEN
               print*,'kick 2a- kw1,mf,m1,m2,ecc,sepf,jorb,vs=',
     &                 kw1,mf,m1,m2,ecc,sepf,jorb,vs(1),vs(2),vs(3)
               CALL kickv(KW1,MF,M1,M2,ECC,SEPF,JORB,VS)
               print*,'kick 2b- kw1,mf,m1,m2,ecc,sepf,jorb,vs=',
     &                 kw1,mf,m1,m2,ecc,sepf,jorb,vs(1),vs(2),vs(3)               
               IF(ECC.GT.1.D0) GOTO 30
            ENDIF
            MF = M2
            KW = KW2
            M2 = MC2
            CALL star(KW2,M02,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
            CALL hrdiag(M02,AJ2,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS,
     &                  R2,L2,KW2,MC2,RC2,MENV,RENV,K22)
            IF(KW2.GE.13.AND.KW.LT.13)THEN
               print*,'kick 3a- kw1,mf,m1,m2,ecc,sepf,jorb,vs=',
     &                 kw2,mf,m2,m1,ecc,sepf,jorb,vs(1),vs(2),vs(3)            
               CALL kickv(KW2,MF,M2,M1,ECC,SEPF,JORB,VS)
               print*,'kick 3b- kw1,mf,m1,m2,ecc,sepf,jorb,vs=',
     &                 kw2,mf,m2,m1,ecc,sepf,jorb,vs(1),vs(2),vs(3)               
               IF(ECC.GT.1.D0) GOTO 30
            ENDIF
         ENDIF
      ENDIF
*
      IF(COEL)THEN
         MC22 = MC2
         IF(KW.EQ.4.OR.KW.EQ.7)THEN
* If making a helium burning star calculate the fractional age 
* depending on the amount of helium that has burnt.
            IF(KW1.LE.3)THEN
               FAGE1 = 0.D0
            ELSEIF(KW1.GE.6)THEN
               FAGE1 = 1.D0
            ELSE
               FAGE1 = (AJ1 - TSCLS1(2))/(TSCLS1(13) - TSCLS1(2))
            ENDIF
            IF(KW2.LE.3.OR.KW2.EQ.10)THEN
               FAGE2 = 0.D0
            ELSEIF(KW2.EQ.7)THEN
               FAGE2 = AJ2/TM2
               MC22 = M2
            ELSEIF(KW2.GE.6)THEN
               FAGE2 = 1.D0
            ELSE
               FAGE2 = (AJ2 - TSCLS2(2))/(TSCLS2(13) - TSCLS2(2))
            ENDIF
         ENDIF
      ENDIF
*
* Now calculate the final mass following coelescence.  This requires a
* Newton-Raphson iteration.
*
      IF(COEL)THEN
*
* Calculate the orbital spin just before coalescence. 
*
         TB = (SEPL/AURSUN)*SQRT(SEPL/(AURSUN*(MC1+MC2)))
         OORB = TWOPI/TB
*
         XX = 1.D0 + ZPARS(7)
         IF(EBINDF.LE.0.D0)THEN
            MF = MC3
            GOTO 20
         ELSE
            CONST = ((M1+M2)**XX)*(M1-MC1+M2-MC22)*EBINDF/EBINDI
         ENDIF
*
* Initial Guess.
*
         MF = MAX(MC1 + MC22,(M1 + M2)*(EBINDF/EBINDI)**(1.D0/XX))
   10    DELY = (MF**XX)*(MF - MC1 - MC22) - CONST
*        IF(ABS(DELY/MF**(1.D0+XX)).LE.1.0D-02) GOTO 20
         IF(ABS(DELY/MF).LE.1.0D-03) GOTO 20
         DERI = MF**ZPARS(7)*((1.D0+XX)*MF - XX*(MC1 + MC22))
         DELMF = DELY/DERI
         MF = MF - DELMF
         GOTO 10
*
* Set the masses and separation.
*
   20    IF(MC22.EQ.0.D0) MF = MAX(MF,MC1+M2)
         M2 = 0.D0
         M1 = MF
         KW2 = 15
*
* Combine the core masses.
*
         IF(KW.EQ.2)THEN
            CALL star(KW,M1,M1,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
            IF(GB(9).GE.MC1)THEN
               M01 = M1
               AJ1 = TM2 + (TSCLS2(1) - TM2)*(AJ1-TM1)/(TSCLS1(1) - TM1)
               CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            ENDIF
         ELSEIF(KW.EQ.7)THEN
            M01 = M1
            CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            AJ1 = TM1*(FAGE1*MC1 + FAGE2*MC22)/(MC1 + MC22)
         ELSEIF(KW.EQ.4.OR.MC2.GT.0.D0.OR.KW.NE.KW1)THEN
            IF(KW.EQ.4) AJ1 = (FAGE1*MC1 + FAGE2*MC22)/(MC1 + MC22)
            MC1 = MC1 + MC2
            MC2 = 0.D0
*
* Obtain a new age for the giant.
*
            CALL gntage(MC1,M1,KW,ZPARS,M01,AJ1)
            CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
         ENDIF
         CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &               R1,L1,KW,MC1,RC1,MENV,RENV,K21)
         JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
         KW1 = KW
         ECC = 0.D0
      ELSE
*
* Check if any eccentricity remains in the orbit by first using 
* energy to circularise the orbit before removing angular momentum. 
* (note this should not be done in case of CE SN ... fix).  
*
         IF(EORBF.LT.ECIRC)THEN
            ECC = SQRT(1.D0 - EORBF/ECIRC)
         ELSE
            ECC = 0.D0
         ENDIF
*
* Set both cores in co-rotation with the orbit on exit of CE, 
*
         TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(M1+M2)))
         OORB = TWOPI/TB
         JORB = M1*M2/(M1+M2)*SQRT(1.D0-ECC*ECC)*SEPF*SEPF*OORB
*        JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
*        JSPIN2 = OORB*(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
*
* or, leave the spins of the cores as they were on entry.
* Tides will deal with any synchronization later.
*
         JSPIN1 = OSPIN1*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
         JSPIN2 = OSPIN2*(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
      ENDIF
   30 SEP = SEPF
      RETURN
      END
***
***
      REAL*8 FUNCTION CORERD(KW,MC,M0,MFLASH)
*     
* A function to determine the radius of the core of a giant-like star.
* NOTE: this is out of date so rc should be obtained using HRDIAG!
* It is still OK to use but bear in mind that the core radius calculated
* for non-degenerate giant cores is only a rough estimate.
*
*     Author : C. A. Tout
*     Date :   26th February 1997
*     Updated 6/1/98 by J. Hurley
*
      IMPLICIT NONE
      INTEGER KW
      REAL*8 MC,MCH,M0,MFLASH
      PARAMETER (MCH = 1.44d0)
*
* First do the black holes and neutron stars.
*
      IF(KW.EQ.14)THEN
         CORERD = 4.24d-06*MC
      ELSEIF(KW.EQ.13)THEN
         CORERD = 1.4d-05
*
* Main sequence stars.
*
      ELSEIF(KW.LE.1.OR.KW.EQ.7)THEN
         CORERD = 0.d0
*
* Core-helium-burning stars, FAGB stars and non-degenerate giant cores.
*
      ELSEIF(KW.EQ.4.OR.KW.EQ.5.OR.(KW.LE.3.AND.M0.GT.MFLASH))THEN
         CORERD = 0.2239d0*MC**0.62d0
*
* The degenerate giants and white dwarfs.
*
      ELSE
         CORERD = 0.0115d0*SQRT(MAX(1.48204d-06,(MCH/MC)**(2.d0/3.d0)
     &                                        - (MC/MCH)**(2.d0/3.d0)))
* 
* Degenerate giants have hot subdwarf cores.
*
         IF(KW.LE.9) CORERD = 5.d0*CORERD
      ENDIF
*
      RETURN
      END
***
***
      SUBROUTINE deltat(kw,age,tm,tn,tscls,dt,dtr)
      implicit none
*
      INTEGER kw
      REAL*8 age,tm,tn,tscls(20)
      REAL*8 dt,dtr
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
*
*     Base new time scale for changes in radius & mass on stellar type.
*
      if(kw.le.1)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.2)then
         dt = pts1*(tscls(1) - tm)
         dtr = tscls(1) - age
      elseif(kw.eq.3)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = MIN(tscls(2),tn) - age
      elseif(kw.eq.4)then
         dt = pts2*tscls(3)
         dtr = MIN(tn,tscls(2) + tscls(3)) - age
      elseif(kw.eq.5)then
         if(age.lt.tscls(9))then
            dt = pts3*(tscls(7) - age)
         else
            dt = pts3*(tscls(8) - age)
         endif
         dtr = MIN(tn,tscls(13)) - age
      elseif(kw.eq.6)then
         if(age.lt.tscls(12))then
            dt = pts3*(tscls(10) - age)
         else
            dt = pts3*(tscls(11) - age)
         endif
         dt = MIN(dt,0.005d0)
         dtr = tn - age
      elseif(kw.eq.7)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.8.or.kw.eq.9)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = tn - age
      else
*        dt = MAX(0.1d0,age*10.d0)
         dt = MAX(0.1d0,dt*10.d0)
         dt = MIN(dt,5.0d+02)
         dtr = dt
      endif
*
      RETURN
      END
***
***
      SUBROUTINE dgcore(kw1,kw2,kw3,m1,m2,m3,ebinde)
*
* A routine to determine the outcome of a collision or coalescence
* of two degenerate cores.
* Entered with kw1,kw2 = 2 or 3 with M <= Mflash, 6, 10, 11 or 12
*
      implicit none
*
      integer kw1,kw2,kw3
*
      real*8 m1,m2,m3,ebinde
      real*8 r1,r2,r3,mhe,mc,mne,ebindi,ebindf,deleb,de,enuc
      real*8 temp,x,y,m0,mflash
      real*8 cvhe,cvc,cvne
      parameter(cvhe=3.1d+07,cvc=8.27d+06,cvne=7.44d+06)
      real*8 ehe,ec,ene
      parameter(ehe=5.812d+17,ec=2.21d+17,ene=2.06d+17)
      real*8 the,tc,gmr,mch
      parameter(the=1.0d+08,tc=1.0d+09,gmr=1.906d+15,mch=1.44d0)
*
      real*8 corerd
      external corerd
*
* Calculate the core radii setting m0 < mflash using dummy values as we
* know it to be true if kw = 2 or 3.
      m0 = 1.d0
      mflash = 2.d0
      r1 = corerd(kw1,m1,m0,mflash)
      r2 = corerd(kw2,m2,m0,mflash)
      r3 = corerd(kw3,m3,m0,mflash)
* Calculate the initial binding energy of the two seperate cores and the
* difference between this and the final binding energy of the new core.
      ebindi = m1*m1/r1 + m2*m2/r2
      ebindf = m3*m3/r3
      deleb = ABS(ebindi - ebindf)
* If an envelope is present reduce its binding energy by the amount
* of energy liberated by the coalescence.
      ebinde = MAX(0.d0,ebinde - deleb)
      if(kw1.gt.3) goto 90
* Distribute the mass into core mass groups where mhe represents Helium 
* core mass, mc represents Carbon core mass and mne represents a Neon
* core mass or any mass that is all converted Carbon.
      mhe = 0.d0
      mc = 0.d0
      mne = 0.d0
      if(kw1.le.3.or.kw1.eq.10)then
         mhe = mhe + m1
      elseif(kw1.eq.12)then
         mne = mne + m1
      else
         mc = mc + m1
      endif
      if(kw2.le.3.or.kw2.eq.10)then
         mhe = mhe + m2
      elseif(kw2.eq.12)then
         mne = mne + m2
      else
         mc = mc + m2
      endif
* Calculate the temperature generated by the merging of the cores.
      temp = (deleb/(cvhe*mhe+cvc*mc+cvne*mne))*gmr
*
* To decide if He is converted to C we use:
*    3He4 -> C12 , T > 10^8 K , 7.274 Mev released, 
* to decide if C is converted further we use:
*    2C12 -> Ne20 + alpha , T > 10^9 K , 4.616 Mev released.
* and to decide if O is converted further we use:
*    2O16 -> P31 + p , T > 10^9 K , 7.677 Mev released.
* To obtain the heat capacity of an O/Ne WD and to gain an idea of the 
* energy released from the further processing of an O/Ne WD we use:
*    2Ne20 + gamma -> O16 + Mg24 +gamma , T > 10^9 K , 4.583 Mev released.
* For a CO core the composition is assumed to be 20% C, 80% O and for 
* an ONe core 80% O, 20% Ne.
*
* Decide if conversion of elements can take place.
*     if(temp.gt.the)then
         x = 1.d0
*     else
*        x = 0.d0
*     endif
*     if(temp.gt.tc)then
*        y = 1.d0
*     else
         y = 0.d0
*     endif
* Calculate the nuclear energy generated from any element conversion.
      enuc = (x*ehe*mhe + y*(ec*mc + ene*mne))/gmr
* Calculate the difference between the binding energy of the star
* (core + envelope) and the nuclear energy. The new star will be
* destroyed in a SN if dE < 0.
      de = (ebindf + ebinde) - enuc
* If the star survives and an envelope is present then reduce the 
* envelope binding energy by the amount of liberated nuclear energy.
* The envelope will not survive if its binding energy is reduced to <= 0.
      if(de.ge.0.d0) ebinde = MAX(0.d0,ebinde - enuc)
* Now determine the final evolution state of the merged core.
      if(de.lt.0.d0) kw3 = 15 
      if(kw3.eq.3)then
         if(x.gt.0.5d0)then
            kw3 = 6
         elseif(ebinde.le.0.d0)then
            kw3 = 10
         endif
      elseif(kw3.eq.4)then
         if(x.gt.0.5d0)then
            kw3 = 6
         elseif(ebinde.le.0.d0)then
            kw3 = 7
         endif
      endif
      if(kw3.eq.6.and.y.lt.0.5d0)then
         if(ebinde.le.0.d0) kw3 = 11
      elseif(kw3.eq.6.and.y.gt.0.5d0)then
         if(ebinde.le.0.d0) kw3 = 12
      endif
      if(kw3.eq.10.and.x.gt.0.5d0) kw3 = 11
      if(kw3.eq.11.and.y.gt.0.5d0) kw3 = 12
      if(kw3.ge.10.and.kw3.le.12.and.m3.ge.mch) kw3 = 15
*
      if(kw3.eq.15) m3 = 0.d0
 90   continue
*
      return
      end
***
***
      SUBROUTINE evolv1b(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &                   epoch,tm,tphys,tphysf,z,zpars,
     &                   dmmax,drmax,vd)
c
c-------------------------------------------------------------c
c
c     Evolves a single star with mass loss.
c     The timestep is not constant but determined by certain criteria.
c
c     Note: no chemical evolution in SSE at this time.  
c
c     Written by Jarrod Hurley 05/10/02 at AMNH, NY. 
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
*
      integer i,j,k,it,kw,kwold,nv
      parameter(nv=50000)
*
      real*8 mass,mt,z,aj,tm,tn,epoch
      real*8 tphys,tphysf,tphys2,tmold,tbgold
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,rc,menv,renv,vd(3),vs(3)
      real*8 ospin,jspin,djt,djmb,k2,k3
      parameter(k3=0.21d0)
      real*8 m0,mt2,mc1,ajhold,rm0
      real*8 dt,dtm,dtr,dtm0,dr,dms,dml,rl
      real*8 tiny
      parameter(tiny=1.0d-14)
      real*8 mlwind,rzamsf,vrotf
      external mlwind,rzamsf,vrotf
      real*8 neta,bwind,mxns
      common /value1/ neta,bwind,mxns
      real*8 dmmax,drmax
*
      dtm = 0.d0
      djt = 0.d0
      mc1 = 0.d0
      rl = 0.d0
      if(ospin.le.0.d0)then
         ospin = 1.0d-10
         jspin = 1.0d-10
      endif
      k2 = 0.15d0
*
      if(tphys.lt.tiny)then
         r = rzamsf(mt)
c         print*,r
         ospin = 45.35d0*vrotf(mt)/r
      endif
      aj = tphys - epoch
*
      CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
c      print*,'hrdiag input 1',mass,aj,mt,tm,tn,r,
c     &        lum,kw,mc,rc,menv,renv,k2

      CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &            r,lum,kw,mc,rc,menv,renv,k2)
c         print*,'hrdiag output 1',mass,aj,mt,tm,tn,r,
c     &        lum,kw,mc,rc,menv,renv,k2

      jspin = ospin*(k2*r*r*(mt-mc)+k3*rc*rc*mc)
*
      do 10 , j = 1,nv
*
         CALL deltat(kw,aj,tm,tn,tscls,dtm,dtr)
         dtm = MIN(dtm,dtr)
*
* Choose minimum of time-scale and remaining interval (> 100 yrs).
*
         dtm = MAX(dtm,1.0d-07*aj)
         dtm = MIN(dtm,tphysf-tphys)
*
         rm0 = r
         ajhold = aj
         tmold = tm
*
* Calculate mass loss.
*
         if(kw.lt.10.and.neta.gt.tiny)then
            dt = 1.0d+06*dtm
            dms = mlwind(kw,lum,r,mt,mc,rl,z)*dt
            dml = mt - mc
            if(dml.lt.dms)then
               dtm = (dml/dms)*dtm
               dms = dml
            endif
         else
            dms = 0.d0
         endif
*
* Limit to 1% mass loss.
*
         if(dms.gt.0.01d0*mt)then
            dtm = 0.01d0*mt*dtm/dms
            dms = 0.01d0*mt
         endif
*
* Calculate the rate of angular momentum loss due to magnetic braking
* and/or mass loss.
*
         if(dtm.gt.tiny)then
            djt = (2.d0/3.d0)*(dms/(1.0d+06*dtm))*r*r*ospin
            if(mt.gt.0.35d0.and.kw.lt.10)then
               djmb = 5.83d-16*menv*(r*ospin)**3/mt
               djt = djt + djmb
            endif
         endif
*
* Update mass and time and reset epoch for a MS (and possibly a HG) star.
*
         if(dms.gt.tiny)then
            mt = mt - dms
            if(kw.le.2.or.kw.eq.7)then
               m0 = mass
               mc1 = mc
               mass = mt
               tbgold = tscls(1)
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(kw.eq.2)then
                  if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                     mass = m0
                  else
                     epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                            (tbgold - tmold)
                     epoch = tphys - epoch
                  endif
               else
                  epoch = tphys - ajhold*tm/tmold
               endif
            endif
         endif
         tphys2 = tphys
         tphys = tphys + dtm
*
* Find the landmark luminosities and timescales as well as setting
* the GB parameters.
*
         aj = tphys - epoch
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
* Find the current radius, luminosity, core mass and stellar type
* given the initial mass, current mass, metallicity and age
*
         kwold = kw
         m0 = mass
c         print*,'hrdiag input',mass,aj,mt,tm,tn,r,
c     &        lum,kw,mc,rc,menv,renv,k2
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r,lum,kw,mc,rc,menv,renv,k2)
c         print*,'hrdiag output',mass,aj,mt,tm,tn,r,
c     &        lum,kw,mc,rc,menv,renv,k2
c         print*,'2',r
*
* Adjust the spin of the star and reset epoch.
*
         jspin = MAX(1.0d-10,jspin - djt*1.0d+06*dtm)
         ospin = jspin/(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         epoch = tphys - aj
*
* Force new NS or BH to have a one second period. 
* 
         if(kw.ne.kwold.and.(kw.eq.13.or.kw.eq.14))then
            ospin = 2.0d+08
            jspin = k3*rc*rc*mc*ospin
            CALL kickv(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,vs)
            do k = 1,3
               vd(k) = vd(k) + vs(k)
            enddo
         endif
*
* Check exit conditions. 
*
         if(tphys.ge.(tphysf-tiny)) goto 90
         if(kw.eq.15)then
            r = 0.d0
            goto 90
         endif
*
 10   continue
*
 90   continue
      tphys = tphysf
*
* Set next timestep. 
*
      dtm0 = 1.0d+10
      if(kw.lt.10)then
         if(neta.gt.tiny)then
            dms = mlwind(kw,lum,r,mt,mc,rl,z)
            if(dms.gt.tiny)then
               dtm0 = dmmax*mt/(1.0d+06*dms)
            endif
         endif
         if(dtm.gt.tiny)then
            dtr = drmax*r*dtm/ABS(r-rm0)
            dtm0 = MIN(dtm0,dtr)
         endif
      endif
      aj = tphys - epoch
      CALL deltat(kw,aj,tm,tn,tscls,dtm,dtr)
      dtm = MIN(dtm,dtr)
      dtm = MIN(dtm,dtm0)
      dtm = MAX(dtm,1.0d-07*aj)
      tphysf = tphys + dtm
*
      RETURN
      END
***
***
      SUBROUTINE evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &                  epoch,tm,tphys,tphysf,dtp,z,zpars)
c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c     Plots the HRD and variables as a function of time.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
*
      integer kw,it,ip,jp,j,kwold,rflag
      integer nv
      parameter(nv=50000)
*
      real*8 mass,z,aj
      real*8 epoch,tphys,tphys2,tmold,tbgold
      real*8 mt,tm,tn,tphysf,dtp,tsave
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,teff,rc,menv,renv,vs(3)
      real*8 ospin,jspin,djt,djmb,k2,k3
      parameter(k3=0.21d0)
      real*8 m0,r1,lum1,mc1,rc1,menv1,renv1,k21
      real*8 dt,dtm,dtr,dr,dtdr,dms,dml,mt2,rl
      real*8 tol,tiny
      parameter(tol=1.0d-10,tiny=1.0d-14)
      real*8 ajhold,rm0,eps,alpha2
      parameter(eps=1.0d-06,alpha2=0.09d0)
      real*8 mlwind,vrotf
      external mlwind,vrotf
      logical iplot,isave
      REAL*8 neta,bwind,mxns
      COMMON /VALUE1/ neta,bwind,mxns
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL scm(50000,10),spp(20,3)
      COMMON /SINGLE/ scm,spp
*
      dtm = 0.d0
      djt = 0.d0
      r = 0.d0
      lum = 0.d0
      mc = 0.d0
      mc1 = 0.d0
      rc = 0.d0
      rl = 0.d0
      if(ospin.le.0.d0)then
         ospin = 1.0d-10
         jspin = 1.0d-10
      endif
      k2 = 0.15d0
      rflag = 0
*
* Setup variables which control the output (if it is required).
*
      ip = 0
      jp = 0
      tsave = tphys
      isave = .true.
      iplot = .false.
      if(dtp.le.0.d0)then
         iplot = .true.
         isave = .false.
         tsave = tphysf
      elseif(dtp.gt.tphysf)then
         isave = .false.
         tsave = tphysf
      endif
* 
      do 10 , j = 1,nv
*
         if(neta.gt.tiny.and.j.gt.1)then
*
* Calculate mass loss from the previous timestep.
*
            dt = 1.0d+06*dtm
            dms = mlwind(kw,lum,r,mt,mc,rl,z)*dt
            if(kw.lt.10)then
               dml = mt - mc
               if(dml.lt.dms)then
                  dtm = (dml/dms)*dtm
                  dms = dml
               endif
            endif
         else
            dms = 0.d0
         endif
*
* Limit to 1% mass loss.
*
         if(dms.gt.0.01d0*mt)then
            dtm = 0.01d0*mt*dtm/dms
            dms = 0.01d0*mt
         endif
*
* Calculate the rate of angular momentum loss due to magnetic braking 
* and/or mass loss.
*
         if(j.gt.1)then
            djt = (2.d0/3.d0)*(dms/(1.0d+06*dtm))*r*r*ospin
            if(mt.gt.0.35d0.and.kw.lt.10)then
               djmb = 5.83d-16*menv*(r*ospin)**3/mt
               djt = djt + djmb
            endif
         endif
*
* Update mass and time and reset epoch for a MS (and possibly a HG) star.
*
         if(dms.gt.0.d0)then
            mt = mt - dms
            if(kw.le.2.or.kw.eq.7)then
               m0 = mass
               mc1 = mc
               mass = mt
               tmold = tm
               tbgold = tscls(1)
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(kw.eq.2)then
                  if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                     mass = m0
                  else
                     epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                            (tbgold - tmold)
                     epoch = tphys - epoch
                  endif
               else
                  epoch = tphys - ajhold*tm/tmold
               endif
            endif
         endif
         tphys2 = tphys
         tphys = tphys + dtm
*
* Find the landmark luminosities and timescales as well as setting
* the GB parameters.
*
         aj = tphys - epoch
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
* Find the current radius, luminosity, core mass and stellar type
* given the initial mass, current mass, metallicity and age
*
         kwold = kw
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r,lum,kw,mc,rc,menv,renv,k2)
*
* If mass loss has occurred and no type change then check that we
* have indeed limited the radius change to 10%.
*
         if(kw.eq.kwold.and.dms.gt.0.d0.and.rflag.ne.0)then
            mt2 = mt + dms
            dml = dms/dtm
            it = 0
 20         dr = r - rm0
            if(ABS(dr).gt.0.1d0*rm0)then
               it = it + 1
               if(it.eq.20.and.kw.eq.4) goto 30
               if(it.gt.30)then
                  WRITE(99,*)' DANGER1! ',it,kw,mass,dr,rm0
                  WRITE(*,*)' STOP: EVOLV1 FATAL ERROR '
                  CALL exit(0)
                  STOP 
               endif
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(rm0,r)*dtdr
               if(it.ge.20) dtm = 0.5d0*dtm
               if(dtm.lt.1.0d-07*aj) goto 30
               dms = dtm*dml
               mt = mt2 - dms
               if(kw.le.2.or.kw.eq.7)then
                  mass = mt
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
                  if(kw.eq.2)then
                     if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                        mass = m0
                     else
                        epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                               (tbgold - tmold)
                        epoch = tphys2 - epoch
                     endif
                  else
                     epoch = tphys2 - ajhold*tm/tmold
                  endif
               endif
               tphys = tphys2 + dtm
               aj = tphys - epoch
               mc = mc1
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                     r,lum,kw,mc,rc,menv,renv,k2)
               goto 20
            endif
 30         continue
         endif
*
* Initialize or adjust the spin of the star.
*
         if(j.eq.1)then
            if(tphys.lt.tiny.and.ospin.lt.0.001d0)then
               ospin = 45.35d0*vrotf(mt)/r
            endif
            jspin = ospin*(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         else
            jspin = MAX(1.0d-10,jspin - djt*1.0d+06*dtm)
            ospin = jspin/(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         endif
*
* Test for changes in evolution type.
*
         if(j.eq.1.or.kw.ne.kwold)then
*
* Force new NS or BH to have a one second period. 
* 
            if(kw.eq.13.or.kw.eq.14)then
               ospin = 2.0d+08
               jspin = k3*rc*rc*mc*ospin
               CALL kickv(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,vs)
            endif
            jp = jp + 1
            spp(jp,1) = tphys
            spp(jp,2) = float(kw)
            if(kw.eq.15)then
               spp(jp,3) = mass 
               goto 90
            else
               spp(jp,3) = mt
            endif
         endif
*
* Record values for plotting and reset epoch.
*
         epoch = tphys - aj
         if((isave.and.tphys.ge.tsave).or.iplot)then
            ip = ip + 1
            scm(ip,1) = tphys
            scm(ip,2) = float(kw)
            scm(ip,3) = mass
            scm(ip,4) = mt
            scm(ip,5) = log10(lum)
            scm(ip,6) = log10(r)
            teff = 1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))
            scm(ip,7) = log10(teff)
            scm(ip,8) = mc
            scm(ip,9) = epoch
            scm(ip,10) = ospin
            if(isave) tsave = tsave + dtp
            if(tphysf.lt.tiny)then
               ip = ip + 1
               do 35 , it = 1,10
                  scm(ip,it) = scm(ip-1,it)
 35            continue
            endif
         endif
*
         if(tphys.ge.tphysf)then
            jp = jp + 1
            spp(jp,1) = tphys
            spp(jp,2) = float(kw)
            spp(jp,3) = mt
            goto 90
         endif
*
* Record radius and current age.
*
         rm0 = r
         ajhold = aj
         if(kw.ne.kwold) kwold = kw
         CALL deltat(kw,aj,tm,tn,tscls,dtm,dtr)
*
* Check for type change.
*
         it = 0
         m0 = mass
         if((dtr-dtm).le.tol.and.kw.le.9)then
*
* Check final radius for too large a jump.
*
            aj = MAX(aj,aj*(1.d0-eps)+dtr)
            mc1 = mc 
            CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r1,lum1,kw,mc1,rc1,menv1,renv1,k21)
            dr = r1 - rm0
            if(ABS(dr).gt.0.1d0*rm0)then
               dtm = dtr - ajhold*eps
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(r1,rm0)*dtdr
               goto 40
            else
               dtm = dtr
               goto 50
            endif
         endif
*
* Limit to a 10% increase in radius assuming no further mass loss
* and thus that the pertubation functions due to small envelope mass
* will not change the radius.
*
 40      aj = ajhold + dtm
         mc1 = mc 
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r1,lum1,kw,mc1,rc1,menv1,renv1,k21)
         dr = r1 - rm0
         it = it + 1
         if(it.eq.20.and.kw.eq.4) goto 50
         if(it.gt.30)then
            WRITE(99,*)' DANGER2! ',it,kw,mass,dr,rm0
            WRITE(*,*)' STOP: EVOLV1 FATAL ERROR '
            CALL exit(0)
            STOP 
         endif
         if(ABS(dr).gt.0.1d0*rm0)then
            dtdr = dtm/ABS(dr)
            dtm = alpha2*MAX(rm0,r1)*dtdr
            if(it.ge.20) dtm = 0.5d0*dtm
            goto 40
         endif
*
 50      continue
*
* Ensure that change of type has not occurred during radius check. 
* This is rare but may occur for HG stars of ZAMS mass > 50 Msun. 
*
         if(kw.ne.kwold)then
            kw = kwold
            mass = m0
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         endif
*
* Choose minimum of time-scale and remaining interval (> 100 yrs).
*
         dtm = MAX(dtm,1.0d-07*aj)
         dtm = MIN(dtm,tsave-tphys)
*
 10   continue
*
 90   continue
*
      tphysf = tphys
      scm(ip+1,1) = -1.0
      spp(jp+1,1) = -1.0
      if(ip.ge.nv)then
         WRITE(99,*)' EVOLV1 ARRAY ERROR ',mass
         WRITE(*,*)' STOP: EVOLV1 ARRAY ERROR '
         CALL exit(0)
         STOP
      endif
*
      RETURN
      END
***
***
      SUBROUTINE evolv2b(kstar,mass0,mass,rad,lumin,massc,radc,
     &                   menv,renv,ospin,rol,dmdt,epoch,tms,
     &                   tphys,tphysf,z,zpars,dmmax,drmax,
     &                   tb,ecc,vd)
      implicit none
*
* Alternative version of evolv2. 
*
      INTEGER loop,iter,intpol,k,kk,ip,jp,j1,j2
      PARAMETER(loop=20000)
      INTEGER kstar(2),kw,kst,kw1,kw2,kmin,kmax
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
*
      REAL*8 km,km0,tphys,tphys0,dtm0,tphys00
      REAL*8 tphysf
      REAL*8 aj(2),aj0(2),epoch(2),tms(2),tbgb(2),tkh(2),dtmi(2)
      REAL*8 mass0(2),mass(2),massc(2),menv(2),mass00(2),mcxx(2)
      REAL*8 rad(2),rol(2),rol0(2),rdot(2),radc(2),renv(2),radx(2)
      REAL*8 lumin(2),k2str(2),q(2),dms(2),dmr(2),dmt(2),dmdt(2)
      REAL*8 dml,vorb2,vwind2,omv2,ivsqm,lacc,vd(3),vs(3)
      REAL*8 sep,dr,tb,dme,tdyn,taum,dm1,dm2,dmchk,qc,dt,pd,rlperi
      REAL*8 m1ce,m2ce,mch,tmsnew,dm22,mew
      PARAMETER(mch=1.44d0)
      REAL*8 yeardy,aursun
      PARAMETER(yeardy=365.24d0,aursun=214.95d0)
      REAL*8 acc1,tiny
      PARAMETER(acc1=3.920659d+08,tiny=1.0d-14)
      REAL*8 ecc,ecc1,tc,tcirc,ttid,ecc2,omecc2,sqome2,sqome3,sqome5
      REAL*8 f1,f2,f3,f4,f5,f,raa2,raa6,eqspin,rg2,tcqr
      REAL*8 k3,mr23yr,twopi
      PARAMETER(k3=0.21d0,mr23yr=0.4311d0)
      REAL*8 jspin(2),ospin(2),jorb,oorb,jspbru,ospbru
      REAL*8 delet,delet1,dspint(2),djspint(2),djtx(2)
      REAL*8 dtj,djorb,djgr,djmb,djt,djtt,rmin,rdisk
      REAL*8 neta,bwind,mxns
      COMMON /VALUE1/ neta,bwind,mxns
      REAL*8 beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
*
      REAL*8 z,tm,tn,m0,mt,rm,lum,mc,rc,me,re,k2,age,dtm,dtr
      REAL*8 tscls(20),lums(10),GB(10),zpars(20)
      REAL*8 dmmax,drmax,dtmx
      REAL*8 mt2,mcx
      LOGICAL coel,com,prec,inttry,snova,sgl,change,check
      LOGICAL supedd,novae,disk
      REAL*8 rl,mlwind,vrotf,corerd
      EXTERNAL rl,mlwind,vrotf,corerd
*
      twopi = 2.d0*ACOS(-1.d0)
*
* Initialize the parameters.
*
      kmin = 1
      kmax = 2
      sgl = .false.
      mt2 = MIN(mass(1),mass(2))
      kst = 0
*
      if(mt2.lt.tiny.or.tb.le.0.d0)then
         sgl = .true.
         if(mt2.lt.tiny)then
            mt2 = 0.d0
            if(mass(1).lt.tiny)then
               if(tphys.lt.tiny)then
                  mass0(1) = 0.01d0
                  mass(1) = mass0(1)
                  kst = 1
               else
                  kmin = 2
               endif
               ospin(1) = 1.0d-10
               jspin(1) = 1.0d-10
               dmdt(1) = 0.d0
            else
               if(tphys.lt.tiny)then
                  mass0(2) = 0.01d0
                  mass(2) = mass0(2)
                  kst = 2
               else
                  kmax = 1
               endif
               ospin(2) = 1.0d-10
               jspin(2) = 1.0d-10
               dmdt(2) = 0.d0
            endif
         endif
         ecc = -1.d0
         tb = 0.d0
         sep = 1.0d+10
         oorb = 0.d0
         jorb = 0.d0
         if(ospin(1).lt.0.d0) ospin(1) = 1.0d-10
         if(ospin(2).lt.0.d0) ospin(2) = 1.0d-10
         q(1) = 1.0d+10
         q(2) = 1.0d+10
         rol(1) = 1.0d+10
         rol(2) = 1.0d+10
      else
         tb = tb/yeardy
         sep = aursun*(tb*tb*(mass(1) + mass(2)))**(1.d0/3.d0)
         oorb = twopi/tb
         jorb = mass(1)*mass(2)/(mass(1)+mass(2))
     &          *SQRT(1.d0-ecc*ecc)*sep*sep*oorb
         if(ospin(1).lt.0.d0) ospin(1) = oorb
         if(ospin(2).lt.0.d0) ospin(2) = oorb
      endif
*
      do 500 , k = kmin,kmax
         age = tphys - epoch(k)
         CALL star(kstar(k),mass0(k),mass(k),tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(mass0(k),age,mass(k),tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kstar(k),mc,rc,me,re,k2)
         aj(k) = age
         epoch(k) = tphys - age
         rad(k) = rm
         lumin(k) = lum
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
*
         if(tphys.lt.tiny.and.ospin(k).le.0.001d0)then
            ospin(k) = 45.35d0*vrotf(mass(k))/rm
         endif
         jspin(k) = ospin(k)*(k2*rm*rm*(mass(k)-mc)+k3*rc*rc*mc)
         if(.not.sgl)then
            q(k) = mass(k)/mass(3-k)
            rol(k) = rl(q(k))*sep
         endif
         rol0(k) = rol(k)
         dmr(k) = 0.d0
         dmt(k) = 0.d0
         dmdt(k) = 0.d0
         djspint(k) = 0.d0
*
         dt = 0.01d0
         CALL deltat(kstar(k),age,tm,tn,tscls,dt,dtr)
         dtmi(k) = MIN(dt,dtr)
*
 500  continue
*
      if(mt2.lt.tiny)then
         sep = 0.d0
         if(kst.gt.0)then
            mass0(kst) = 0.d0
            mass(kst) = 0.d0
            kmin = 3 - kst
            kmax = kmin
         endif
      endif
*
* On the first entry the previous timestep is zero to prevent mass loss.
*
      dtm = 0.d0
      delet = 0.d0
      djorb = 0.d0
      dtmx = 1.0d+10
      change = .false.
      check = .false.
      if(tphysf.lt.0.d0)then
         tphysf = -1.d0*tphysf
         check = .true.
      endif
*
      if(tphys.ge.tphysf) goto 140
*
 4    iter = 0
      intpol = 0
      inttry = .false.
      prec = .false.
      snova = .false.
      coel = .false.
      com = .false.
      tphys0 = tphys
      ecc1 = ecc
      j1 = 1
      j2 = 2
      if(kstar(1).ge.10.and.kstar(1).le.14) dtmi(1) = 0.01d0
      if(kstar(2).ge.10.and.kstar(2).le.14) dtmi(2) = 0.01d0
      dmr(1) = 0.d0
      dmt(1) = 0.d0
      dmr(2) = 0.d0
      dmt(2) = 0.d0
*
 5    kw1 = kstar(1)
      kw2 = kstar(2)
*
      dt = 1.0d+06*dtm
      eqspin = 0.d0
      djtt = 0.d0
      dtmx = 1.0d+10
*
      if(intpol.eq.0.and.ABS(dtm).gt.tiny.and..not.sgl)then
         vorb2 = acc1*(mass(1)+mass(2))/sep
         ivsqm = 1.d0/SQRT(1.d0-ecc*ecc)
         do 501 , k = 1,2
*
* Calculate wind mass loss from the previous timestep.
*
            if(neta.gt.tiny)then
               rlperi = rol(k)*(1.d0-ecc)
               dmr(k) = mlwind(kstar(k),lumin(k),rad(k),mass(k),
     &                         massc(k),rlperi,z)
*
* Calculate how much of wind mass loss from companion will be
* accreted (Boffin & Jorissen, A&A 1988, 205, 155).
*
               vwind2 = 2.d0*beta*acc1*mass(k)/rad(k)
               omv2 = (1.d0 + vorb2/vwind2)**(3.d0/2.d0)
               dmt(3-k) = ivsqm*acc2*dmr(k)*((acc1*mass(3-k)/vwind2)**2)
     &                    /(2.d0*sep*sep*omv2)
               dmt(3-k) = MIN(dmt(3-k),0.8d0*dmr(k))
            else
               dmr(k) = 0.d0
               dmt(3-k) = 0.d0
            endif
 501     continue
*
* Calculate orbital angular momentum change due to wind mass loss.
*
         ecc2 = ecc*ecc
         omecc2 = 1.d0 - ecc2
         sqome2 = SQRT(omecc2)
*
         djorb = ((dmr(1)+q(1)*dmt(1))*mass(2)*mass(2) + 
     &            (dmr(2)+q(2)*dmt(2))*mass(1)*mass(1))*
     &           sep*sep*sqome2*oorb/(mass(1)+mass(2))**2
         delet = ecc*(dmt(1)*(0.5d0/mass(1) + 1.d0/(mass(1)+mass(2))) +
     &                dmt(2)*(0.5d0/mass(2) + 1.d0/(mass(1)+mass(2))))
*
* For very close systems include angular momentum loss owing to 
* gravitational radiation. 
*
         if(sep.le.10.d0)then
            djgr = 8.315d-10*mass(1)*mass(2)*(mass(1)+mass(2))/
     &             (sep*sep*sep*sep)
            f1 = (19.d0/6.d0) + (121.d0/96.d0)*ecc2
            sqome5 = sqome2**5
            delet1 = djgr*ecc*f1/sqome5
            djgr = djgr*jorb*(1.d0+0.875d0*ecc2)/sqome5
            djorb = djorb + djgr
            delet = delet + delet1
         endif
*
         do 502 , k = 1,2
*
* Calculate change in the intrinsic spin of the star.
*
            djtx(k) = (2.d0/3.d0)*xi*dmt(k)*rad(3-k)*rad(3-k)*ospin(3-k)
            djspint(k) = (2.d0/3.d0)*(dmr(k)*rad(k)*rad(k)*ospin(k)) -
     &                   djtx(k)
*
* Include magnetic braking for stars that have appreciable convective 
* envelopes. This includes MS stars with M < 1.25, HG stars near the GB 
* and giants. MB is not allowed for fully convective MS stars. 
*
            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
               djspint(k) = djspint(k) + djmb
*
* Limit to a 3% angular momentum change for the star owing to MB. 
* This is found to work best with the maximum iteration of 20000, 
* i.e. does not create an excessive number of iterations, while not 
* affecting the evolution outcome when compared with a 2% restriction.  
*
               if(djmb.gt.tiny)then
                  dtj = 0.03d0*jspin(k)/ABS(djmb)
                  dt = MIN(dt,dtj)
               endif
            endif
*
* Calculate circularization, orbital shrinkage and spin up.
*
            dspint(k) = 0.d0
            if(((kstar(k).le.9.and.rad(k).ge.0.01d0*rol(k)).or.
     &         (kstar(k).ge.10.and.k.eq.j1)).and.tflag.gt.0)then
*
               raa2 = (rad(k)/sep)**2
               raa6 = raa2**3
*
* Hut's polynomials.
*
               f5 = 1.d0+ecc2*(3.d0+ecc2*0.375d0)
               f4 = 1.d0+ecc2*(1.5d0+ecc2*0.125d0)
               f3 = 1.d0+ecc2*(3.75d0+ecc2*(1.875d0+ecc2*7.8125d-02))
               f2 = 1.d0+ecc2*(7.5d0+ecc2*(5.625d0+ecc2*0.3125d0))
               f1 = 1.d0+ecc2*(15.5d0+ecc2*(31.875d0+ecc2*(11.5625d0
     &                  +ecc2*0.390625d0)))
*
               if((kstar(k).eq.1.and.mass(k).ge.1.25d0).or.
     &            kstar(k).eq.4.or.kstar(k).eq.7)then
*
* Radiative damping (Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329).
*
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*rad(k)*rad(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k)
               elseif(kstar(k).le.9)then
*
* Convective damping (Hut, 1981, A&A, 99, 126).
*
                  tc = mr23yr*(menv(k)*renv(k)*(rad(k)-0.5d0*renv(k))/
     &                 (3.d0*lumin(k)))**(1.d0/3.d0)
                  ttid = twopi/(1.0d-10 + ABS(oorb - ospin(k)))
                  f = MIN(1.d0,(ttid/(2.d0*tc)**2))
                  tcqr = 2.d0*f*q(3-k)*raa6*menv(k)/(21.d0*tc*mass(k))
                  rg2 = (k2str(k)*(mass(k)-massc(k)))/mass(k)
               else
*
* Degenerate damping (Campbell, 1984, MNRAS, 207, 433)
*
                  f = 7.33d-09*(lumin(k)/mass(k))**(5.d0/7.d0)
                  tcqr = f*q(3-k)*q(3-k)*raa2*raa2/(1.d0+q(3-k))
                  rg2 = k3
               endif
*
* Circularization.
*
               sqome3 = sqome2**3
               delet1 = 27.d0*tcqr*(1.d0+q(3-k))*raa2*(ecc/sqome2**13)*
     &                  (f3 - (11.d0/18.d0)*sqome3*f4*ospin(k)/oorb)
               tcirc = ecc/(ABS(delet1) + 1.0d-20)
               delet = delet + delet1
*
* Spin up of star.
*
               dspint(k) = (3.d0*q(3-k)*tcqr/(rg2*omecc2**6))*
     &                     (f2*oorb - sqome3*f5*ospin(k))
*
* Calculate the equilibrium spin at which no angular momentum 
* can be transferred.
*
               eqspin = oorb*f2/(sqome3*f5)
*
* Calculate angular momentum change for the star owing to tides. 
*
               djt = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               if(kstar(k).le.6.or.ABS(djt)/jspin(k).gt.0.1d0)then
                  djtt = djtt + djt
               endif
            endif
 502     continue
*
* Limit to 2% orbital angular momentum change.
*
         djtt = djtt + djorb 
         if(ABS(djtt).gt.tiny)then
            dtj = 0.02d0*jorb/ABS(djtt)
            dt = MIN(dt,dtj)
         endif
         dtm = dt/1.0d+06
*
      elseif(ABS(dtm).gt.tiny.and.sgl)then
         do 503 , k = kmin,kmax
            if(neta.gt.tiny)then
               rlperi = 0.d0
               dmr(k) = mlwind(kstar(k),lumin(k),rad(k),mass(k),
     &                         massc(k),rlperi,z)
            else
               dmr(k) = 0.d0
            endif
            dmt(k) = 0.d0
            djspint(k) = (2.d0/3.d0)*dmr(k)*rad(k)*rad(k)*ospin(k)
 503     continue
      endif
*
      do 504 , k = kmin,kmax
*
         dms(k) = (dmr(k) - dmt(k))*dt
         if(kstar(k).lt.10)then
            dml = mass(k) - massc(k)
            if(dml.lt.dms(k))then
               dtm = (dml/dms(k))*dtm
               if(k.eq.2) dms(1) = dms(1)*dml/dms(2)
               dms(k) = dml
               dt = 1.0d+06*dtm
            endif
*
* Limit to 1% mass loss.
*
            if(dms(k).gt.0.01d0*mass(k))then
               dtm = 0.01d0*mass(k)*dtm/dms(k)
               if(k.eq.2) dms(1) = dms(1)*0.01d0*mass(2)/dms(2)
               dms(k) = 0.01d0*mass(k)
               dt = 1.0d+06*dtm
            endif
*
            if(dmr(k).gt.tiny)then
               dtmx = MIN(dtmx,dmmax*mass(k)/(1.0d+06*dmr(k)))
            endif
*
         endif
*
 504  continue
*
* Update mass and intrinsic spin (checking that the star is not spun 
* past the equilibrium) and reset epoch for a MS (and possibly a HG) star. 
*
      do 505 , k = kmin,kmax
*
         if(eqspin.gt.0.d0.and.ABS(dspint(k)).gt.tiny)then
            if(intpol.eq.0)then
               if(dspint(k).ge.0.d0)then
                  dspint(k) = MIN(dspint(k),(eqspin-ospin(k))/dt)
               else
                  dspint(k) = MAX(dspint(k),(eqspin-ospin(k))/dt)
               endif
               djt = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               djorb = djorb + djt
               djspint(k) = djspint(k) - djt
            endif
         endif
*
         jspin(k) = MAX(1.0d-10,jspin(k) - djspint(k)*dt)
*
* Ensure that the star does not spin up beyond break-up.
*
         ospbru = twopi*SQRT(mass(k)*aursun**3/rad(k)**3)
         jspbru = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &             k3*massc(k)*radc(k)*radc(k))*ospbru
         if(jspin(k).gt.jspbru.and.ABS(dtm).gt.tiny)then
            mew = 1.d0
            if(djtx(k).gt.0.d0)then
               mew = MIN(mew,(jspin(k) - jspbru)/djtx(k))
            endif
            jspin(k) = jspbru
* If excess material should not be accreted, activate next line.
*           dms(k) = dms(k) + (1.d0 - mew)*dmt(k)*dt
         endif
*
         if(ABS(dms(k)).gt.tiny)then
            mass(k) = mass(k) - dms(k)
            if(kstar(k).le.2.or.kstar(k).eq.7)then
               m0 = mass0(k)
               mass0(k) = mass(k)
               CALL star(kstar(k),mass0(k),mass(k),tm,tn,tscls,
     &                   lums,GB,zpars)
               if(kstar(k).eq.2)then
                  if(GB(9).lt.massc(k).or.m0.gt.zpars(3))then
                     mass0(k) = m0
                  else
                     epoch(k) = tm + (tscls(1) - tm)*(aj(k)-tms(k))/
     &                               (tbgb(k) - tms(k))
                     epoch(k) = tphys - epoch(k)
                  endif
               else
                  epoch(k) = tphys - aj(k)*tm/tms(k)
               endif
            endif
         endif
*
 505  continue
*
      if(.not.sgl)then
*
         ecc1 = ecc1 - delet*dt
         ecc = MAX(ecc1,0.d0)
         if(ecc.lt.1.0d-10) ecc = 0.d0
*
         if(ecc.ge.1.d0) goto 135
*
         jorb = jorb - djorb*dt
         sep = (mass(1) + mass(2))*jorb*jorb/
     &         ((mass(1)*mass(2)*twopi)**2*aursun**3*(1.d0-ecc*ecc))
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
      endif
*
* Advance the time.
*
      if(intpol.eq.0)then
         tphys0 = tphys
         dtm0 = dtm
      endif
      tphys = tphys + dtm
*
      do 6 , k = kmin,kmax
*
* Acquire stellar parameters (M, R, L, Mc & K*) at apparent evolution age.
*
         age = tphys - epoch(k)
         aj0(k) = age
         kw = kstar(k)
         m0 = mass0(k)
         mt = mass(k)
         mc = massc(k)
         if(intpol.eq.0) mcxx(k) = mc
         if(intpol.gt.0) mc = mcxx(k)
         mass00(k) = m0
*
*
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(m0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kw,mc,rc,me,re,k2)
*
* At this point there may have been a supernova.
*
         if(kw.ne.15)then
            ospin(k) = jspin(k)/(k2*(mt-mc)*rm*rm+k3*mc*rc*rc)
         endif
         if(kw.ne.kstar(k).and.kstar(k).le.12.and.
     &      (kw.eq.13.or.kw.eq.14))then
            mt2 = mass(3-k)
            if(sgl) mt2 = 0.d0
            print*,'kick 4a- kw,k,mf,m1,m2,ecc,sepf,jorb,vs=',
     &              kw,k,mass(k),mt,mt2,ecc,sep,jorb,vs(1),vs(2),vs(3)
            CALL kickv(kw,mass(k),mt,mt2,ecc,sep,jorb,vs)
            print*,'kick 4b- kw,k,mf,m1,m2,ecc,sepf,jorb,vs=',  
     &              kw,k,mass(k),mt,mt2,ecc,sep,jorb,vs(1),vs(2),vs(3) 
            print*, 'after kickv in evolv2b, ecc ',ecc
            do kk = 1,3
               vd(kk) = vd(kk) + vs(kk)
            enddo
            if(.not.sgl)then
               if(ecc.gt.1.d0)then
                  ospin(k) = 2.0d+08
                  jspin(k) = k3*rc*rc*mc*ospin(k)
                  kstar(k) = kw
                  mass(k) = mt
                  epoch(k) = tphys - age
                  goto 135
               endif
               tb = (sep/aursun)*SQRT(sep/(aursun*(mt+mass(3-k))))
               oorb = twopi/tb
               snova = .true.
            endif
         endif
         if(kw.ne.kstar(k))then
            change = .true.
            mass(k) = mt
            dtmi(k) = 0.01d0
            if((kw.gt.6.and.kstar(k).le.6).or.
     &         (kw.gt.9.and.kstar(k).le.9))then
               dmr(k) = 0.d0
               dmt(3-k) = 0.d0
            endif
            if(kw.eq.15)then
               kstar(k) = kw
               goto 135
            endif
            mass0(k) = m0
            epoch(k) = tphys - age
*
* Force new NS or BH to have a second period.
*
            if(kw.eq.13.or.kw.eq.14)then
               ospin(k) = 2.0d+08
               jspin(k) = k3*rc*rc*mc*ospin(k)
            endif
         endif
*
* Set radius derivative for later interpolation.
*
         if(ABS(dtm).gt.tiny)then
            rdot(k) = ABS(rm - rad(k))/dtm
            if(intpol.eq.0.and.kw.eq.kstar(k))then
               if(rdot(k).gt.0.d0) dtmx = MIN(dtmx,drmax*rm/rdot(k))
            endif
         else
            rdot(k) = 0.d0
         endif
*
*     Base new time scale for changes in radius & mass on stellar type.
*
         dt = dtmi(k)
         CALL deltat(kw,age,tm,tn,tscls,dt,dtr)
*
* Choose minimum of time-scale and remaining interval.
*
         dtmi(k) = MIN(dt,dtr)
*
* Save relevent solar quantities.
*
         aj(k) = age
         kstar(k) = kw
         rad(k) = rm
         lumin(k) = lum
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
         dmdt(k) = dmt(k) - dmr(k)
*
 6    continue
*
      if(.not.sgl)then
*
* Determine the mass ratios.
*
         do 506 , k = 1,2
            q(k) = mass(k)/mass(3-k)
 506     continue
*
* Determine the Roche lobe radii and adjust the radius derivative.
*
         do 507 , k = 1,2
            rol(k) = rl(q(k))*sep
            if(ABS(dtm).gt.tiny)then
               rdot(k) = rdot(k) + (rol(k) - rol0(k))/dtm
               rol0(k) = rol(k)
            endif
 507     continue
      else
         do 508 , k = kmin,kmax
            rol(k) = 10000.d0*rad(k)
 508     continue
      endif
*
      if(snova)then
         dtm = 0.d0
         goto 4
      endif
*
* If not interpolating set the next timestep.
*
      if(intpol.eq.0)then
         dtm = MAX(1.0d-07*tphys,MIN(dtmi(1),dtmi(2)))
         dtm = MIN(dtm,tphysf-tphys)
         if(iter.eq.0) dtm0 = dtm
      endif
      if(sgl) goto 98
*
* Set j1 to the donor - the primary
* and j2 to the accretor - the secondary.
*
      if(intpol.eq.0)then
         if(rad(1)/rol(1).ge.rad(2)/rol(2))then
            j1 = 1
            j2 = 2
         else
            j1 = 2
            j2 = 1
         endif
      endif
*
* Test whether Roche lobe overflow has begun.
*
      if(rad(j1).gt.rol(j1))then
*
* Interpolate back until the primary is just filling its Roche lobe.
*
         if(rad(j1).ge.1.002d0*rol(j1))then
            if(intpol.eq.0) tphys00 = tphys
            intpol = intpol + 1
            if(iter.eq.0) goto 7
            if(inttry) goto 7
            if(intpol.ge.100)then
               ecc = 2.d0
               goto 140 
            endif
            dr = rad(j1) - 1.001d0*rol(j1)
            if(ABS(rdot(j1)).lt.tiny.or.prec)then
               goto 7
            endif
            dtm = -dr/ABS(rdot(j1))
            if(ABS(tphys0-tphys).gt.tiny) dtm = MAX(dtm,tphys0-tphys)
            if(kstar(1).ne.kw1)then
               kstar(1) = kw1
               mass0(1) = mass00(1)
               epoch(1) = tphys - aj0(1)
            endif
            if(kstar(2).ne.kw2)then
               kstar(2) = kw2
               mass0(2) = mass00(2)
               epoch(2) = tphys - aj0(2)
            endif
         else
*
* Enter Roche lobe overflow
*
            if(tphys.ge.tphysf) goto 140
            goto 7
         endif
      else
*
* Check if already interpolating.
*
         if(intpol.gt.0)then
            intpol = intpol + 1
            if(intpol.ge.80)then
               inttry = .true.
            endif
            if(ABS(rdot(j1)).lt.tiny)then
               prec = .true.
               dtm = 1.0d-07*tphys
            else
               dr = rad(j1) - 1.001d0*rol(j1)
               dtm = -dr/ABS(rdot(j1))
            endif
            if((tphys+dtm).ge.tphys00)then
*
* If this occurs then most likely the star is a high mass type 4
* where the radius can change very sharply or possibly there is a 
* discontinuity in the radius as a function of time and HRDIAG
* needs to be checked!
*
               dtm = 0.5d0*(tphys00 - tphys0)
               dtm = MAX(dtm,1.0d-10)
               prec = .true.
            endif
            tphys0 = tphys
         endif
      endif
*
* Check for collision at periastron.
*
      pd = sep*(1.d0 - ecc)
      if(pd.lt.(rad(1)+rad(2)).and.intpol.eq.0) goto 130
*
* Go back for the next step or interpolation.
*
 98   continue
      if(tphys.ge.tphysf.and.intpol.eq.0) goto 140
      if(check.and.change) goto 140
*
      iter = iter + 1
*
      if(iter.ge.loop)then
         ecc = 2.d0
         goto 140
      endif
      goto 5
*
* Set the nuclear timescale in years and slow-down factor.
*
 7    km0 = dtm0*1.0d+03/tb
      if(km0.lt.tiny) km0 = 0.5d0
*
* Force co-rotation of primary and orbit to ensure that the tides do not
* lead to unstable Roche (not currently used).
*
*     if(ospin(j1).gt.1.05d0*oorb)then
*        ospin(j1) = oorb
*        jspin(j1) = (k2str(j1)*rad(j1)*rad(j1)*(mass(j1)-massc(j1))+
*    &                k3*radc(j1)*radc(j1)*massc(j1))*ospin(j1)
*     endif
*
      coel = .false.
      radx(j1) = rol(j1)
      radx(j2) = rad(j2)
*
      if(iter.gt.0.or.tphys.le.tiny)then
         change = .true.
         dtm = 0.d0
*
* Make an estimate of starting mass transfer rate for new RLOF. 
*
         dm1 = 3.0d-06*(LOG(rad(j1)/rol(j1))**3)*
     &         MIN(mass(j1),5.d0)**2
         if(kstar(j1).ge.10)then
            dm1 = dm1*1.0d+03*mass(j1)/MAX(rad(j1),1.0d-04)
         endif
         if(kstar(j1).ge.2.and.kstar(j1).le.9.and.kstar(j1).ne.7)then
            tkh(j1) = 1.0d+07*mass(j1)/(rad(j1)*lumin(j1))
            tkh(j1) = tkh(j1)*(mass(j1) - massc(j1))
            dm1 = MIN(dm1,mass(j1)/tkh(j1))
         else
            tdyn = 5.05d-05*SQRT(rad(j1)**3/mass(j1))
            dm1 = MIN(dm1,mass(j1)/tdyn)
         endif
         dmdt(j1) = -1.d0*dm1
         dmdt(j2) = 0.d0
*
         goto 85
      endif
      iter = 0
*
* Eddington limit for accretion on to the secondary in one orbit.
*
 8    dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(j2)*tb
      supedd = .false.
      novae = .false.
      disk = .false.
      dtmx = 1.0d+10
*
* Determine whether the transferred material forms an accretion 
* disk around the secondary or hits the secondary in a direct 
* stream, by using eq.(1) of Ulrich & Burger (1976, ApJ, 206, 509) 
* fitted to the calculations of Lubow & Shu (1974, ApJ, 198, 383). 
*
*     if(kstar(j2).ge.10) disk = .true.
      rmin = 0.0425d0*sep*(q(j2)*(1.d0+q(j2)))**(1.d0/4.d0)
      if(rmin.gt.rad(j2)) disk = .true.
*
* Kelvin-Helmholtz time from the modified classical expression.
*
      do 13 , k = 1,2
         tkh(k) = 1.0d+07*mass(k)/(rad(k)*lumin(k))
         if(kstar(k).le.1.or.kstar(k).eq.7.or.kstar(k).ge.10)then
            tkh(k) = tkh(k)*mass(k)
         else
            tkh(k) = tkh(k)*(mass(k) - massc(k))
         endif
         dmdt(k) = 0.d0
 13   continue
*
* Dynamical timescale for the primary.
*
      tdyn = 5.05d-05*SQRT(rad(j1)**3/mass(j1))
*
* Identify special cases.
*
      if(kstar(j1).eq.2)then
         qc = 4.d0
      elseif(kstar(j1).eq.3.or.kstar(j1).eq.5.or.kstar(j1).eq.6)then
         qc = (1.67d0-zpars(7)+2.d0*(massc(j1)/mass(j1))**5)/2.13d0
*
* Alternatively use condition of Hjellming & Webbink, 1987, ApJ, 318, 794. 
*        qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
      elseif(kstar(j1).eq.8.or.kstar(j1).eq.9)then
         qc = 0.784d0
      else
         qc = 3.d0
      endif
*
      if(kstar(j1).eq.0.and.q(j1).gt.0.695d0)then
*
* This will be dynamical mass transfer of a similar nature to
* common-envelope evolution.  The result is always a single
* star placed in *2.
*
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         if(kstar(j2).le.1)then
*
* Restrict accretion to thermal timescale of secondary.
*
            dm2 = taum/tkh(j2)*dm1
            mass(j2) = mass(j2) + dm2
*
* Rejuvenate if the star is still on the main sequence.
*
            mass0(j2) = mass(j2)
            CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,
     &                tscls,lums,GB,zpars)
* If the star has no convective core then the effective age decreases,
* otherwise it will become younger still.
            if(mass(j2).lt.0.35d0.or.mass(j2).gt.1.25d0)then
               aj(j2) = tmsnew/tms(j2)*aj(j2)*(mass(j2) - dm2)/mass(j2)
            else
               aj(j2) = tmsnew/tms(j2)*aj(j2)
            endif
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.6)then
*
* Add all the material to the giant's envelope.
*
            dm2 = dm1
            mass(j2) = mass(j2) + dm2
            if(kstar(j2).eq.2)then
               mass0(j2) = mass(j2)
               CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &                   lums,GB,zpars)
               aj(j2) = tmsnew + tscls(1)*(aj(j2)-tms(j2))/tbgb(j2)
               epoch(j2) = tphys - aj(j2)
            endif
         elseif(kstar(j2).le.12)then
*
* Form a new giant envelope.
*
            dm2 = dm1
            kst = ktype(kstar(j1),kstar(j2))
            if(kst.gt.100) kst = kst - 100
            if(kst.eq.4)then
               aj(j2) = aj(j2)/tms(j2)
               massc(j2) = mass(j2)
            endif
*
* Check for planets or low-mass WDs.
*
            if((kstar(j2).eq.10.and.mass(j2).lt.0.05d0).or.
     &         (kstar(j2).ge.11.and.mass(j2).lt.0.5d0))then
               kst = kstar(j1)
               mass(j1) = mass(j2) + dm2
               mass(j2) = 0.d0
            else
               mass(j2) = mass(j2) + dm2
               CALL gntage(massc(j2),mass(j2),kst,zpars,
     &                     mass0(j2),aj(j2))
               epoch(j2) = tphys - aj(j2)
            endif
            kstar(j2) = kst
         else
*
* The neutron star or black hole simply accretes at the Eddington rate.
*
            dm2 = MIN(dme*taum/tb,dm1)
            if(dm2.lt.dm1) supedd = .true. 
            mass(j2) = mass(j2) + dm2
         endif
         coel = .true.
         if(mass(j2).gt.0.d0)then
            mass(j1) = 0.d0
            kstar(j1) = 15
         else
            kstar(j1) = kstar(j2)
            kstar(j2) = 15
            mass(j2) = 0.d0
         endif
         goto 135
      elseif(((ABS(ABS(2*kstar(j1)-11)-3).eq.2.or.kstar(j1).eq.9).
     &        and.(q(j1).gt.qc.or.radx(j1).le.radc(j1))).or.
     &        (kstar(j1).eq.2.and.q(j1).gt.qc).or.
     &        (kstar(j1).eq.4.and.q(j1).gt.qc))then
*
* Common-envelope evolution.
*
         m1ce = mass(j1)
         m2ce = mass(j2)
         CALL comenv(mass0(j1),mass(j1),massc(j1),aj(j1),jspin(j1),
     &               kstar(j1),mass0(j2),mass(j2),massc(j2),aj(j2),
     &               jspin(j2),kstar(j2),zpars,ecc,sep,jorb,coel)
*
         epoch(j1) = tphys - aj(j1)
         if(coel)then
            com = .true.
            goto 135
         endif
         epoch(j2) = tphys - aj(j2)
         if(ecc.gt.1.d0)then
            if(kstar(1).ge.13)then
               rc = corerd(kstar(1),mass(1),mass(1),zpars(2))
               ospin(1) = jspin(1)/(k3*rc*rc*mass(1))
            endif
            if(kstar(2).ge.13)then
               rc = corerd(kstar(2),mass(2),mass(2),zpars(2))
               ospin(2) = jspin(2)/(k3*rc*rc*mass(2))
            endif
            goto 135
         endif
*
* Next step should be made without changing the time.
*
         dm1 = m1ce - mass(j1)
         dm2 = mass(j2) - m2ce
         dm22 = dm2
         dtm = 0.d0
*
* Reset orbital parameters as separation may have changed.
*
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
         change = .true.
      elseif(kstar(j1).ge.10.and.kstar(j1).le.12.and.
     &       q(j1).gt.0.628d0)then
*
* Dynamic transfer from a white dwarf.  Secondary will have KW > 9.
*
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         dm2 = MIN(dme*taum/tb,dm1)
         if(dm2.lt.dm1) supedd = .true. 
         mass(j2) = mass(j2) + dm2
*
         if(kstar(j1).eq.10.and.kstar(j2).eq.10)then
*
* Assume the energy released by ignition of the triple-alpha reaction 
* is enough to destroy the star. 
*
            kstar(j2) = 15
            mass(j2) = 0.d0
         elseif(kstar(j1).eq.10.or.kstar(j2).eq.10)then
*
* Should be helium overflowing onto a CO or ONe core in which case the 
* helium swells up to form a giant envelope so a HeGB star is formed. 
* Allowance for the rare case of CO or ONe flowing onto He is made. 
*
            kst = 9
            if(kstar(j2).eq.10) massc(j2) = dm2
            CALL gntage(massc(j2),mass(j2),kst,zpars,mass0(j2),aj(j2))
            kstar(j2) = kst
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.12)then
            mass0(j2) = mass(j2)
            if(kstar(j1).eq.12.and.kstar(j2).eq.11)then
*
* Mixture of ONe and CO will result in an ONe product. 
*
               kstar(j2) = 12
            endif
         endif
         kstar(j1) = 15
         mass(j1) = 0.d0
*
* Might be a supernova that destroys the system.
*
         if(kstar(j2).le.11.and.mass(j2).gt.mch)then
            kstar(j2) = 15
            mass(j2) = 0.d0
         endif
         coel = .true.
         goto 135
      elseif(kstar(j1).eq.13)then
*
* Gamma ray burster?
*
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
         kstar(j2) = 14
         coel = .true.
         goto 135
      elseif(kstar(j1).eq.14)then
*
* Both stars are black holes.  Let them merge quietly.
*
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
         coel = .true.
         goto 135
      else
*
* Mass transfer in one Kepler orbit.
*
         dm1 = 3.0d-06*tb*(LOG(rad(j1)/rol(j1))**3)*
     &         MIN(mass(j1),5.d0)**2
         if(kstar(j1).eq.2)then
            mew = (mass(j1) - massc(j1))/mass(j1)
            dm1 = MAX(mew,0.01d0)*dm1
         elseif(kstar(j1).ge.10)then
*           dm1 = dm1*1.0d+03/MAX(rad(j1),1.0d-04)
            dm1 = dm1*1.0d+03*mass(j1)/MAX(rad(j1),1.0d-04)
         endif
         kst = kstar(j2)
*
* Possibly mass transfer needs to be reduced if primary is rotating 
* faster than the orbit (not currently implemented). 
*
*        spnfac = MIN(3.d0,MAX(ospin(j1)/oorb,1.d0))
*        dm1 = dm1/spnfac**2
*
* Limit mass transfer to the thermal rate for remaining giant-like stars
* and to the dynamical rate for all others.
*
         if(kstar(j1).ge.2.and.kstar(j1).le.9.and.kstar(j1).ne.7)then
***
* JH_temp ... this may be good for HG RLOF??
*           if(kstar(j1).eq.2)then
*              mew = rad(j1)/rol(j1) - 1.d0
*              mew = 2.d0*mew
*              dm1 = dm1*10.d0**mew
*           endif
***
            dm1 = MIN(dm1,mass(j1)*tb/tkh(j1))
         elseif(rad(j1).gt.10.d0*rol(j1).or.(kstar(j1).le.1.and.
     &          kstar(j2).le.1.and.q(j1).gt.qc))then
*
* Allow the stars to merge with the product in *1.
*
            m1ce = mass(j1)
            m2ce = mass(j2)
            CALL mix(mass0,mass,aj,kstar,zpars)
            dm1 = m1ce - mass(j1)
            dm2 = mass(j2) - m2ce
*
* Next step should be made without changing the time.
*
            dtm = 0.d0
            epoch(1) = tphys - aj(1)
            coel = .true.
            goto 135
         else
            dm1 = MIN(dm1,mass(j1)*tb/tdyn)
         endif
*
* Calculate wind mass loss from the stars during one orbit.
*
         vorb2 = acc1*(mass(1)+mass(2))/sep
         ivsqm = 1.d0/SQRT(1.d0-ecc*ecc)
         do 14 , k = 1,2
            if(neta.gt.tiny)then
               rlperi = rol(k)*(1.d0-ecc)
               dmr(k) = mlwind(kstar(k),lumin(k),radx(k),
     &                         mass(k),massc(k),rlperi,z)
               vwind2 = 2.d0*beta*acc1*mass(k)/radx(k)
               omv2 = (1.d0 + vorb2/vwind2)**(3.d0/2.d0)
               dmt(3-k) = ivsqm*acc2*dmr(k)*((acc1*mass(3-k)/vwind2)**2)
     &                    /(2.d0*sep*sep*omv2)
               dmt(3-k) = MIN(dmt(3-k),dmr(k))
            else
               dmr(k) = 0.d0
               dmt(3-k) = 0.d0
            endif
 14      continue
*
         do 15 , k = 1,2
            dms(k) = (dmr(k)-dmt(k))*tb
 15      continue
*
* Increase time-scale to relative mass loss of 0.5% but not more than twice.
* KM is the number of orbits for the timestep.
*
         km = MIN(2.d0*km0,5.0d-03/
     &            MAX(ABS(dm1+dms(j1))/mass(j1),dms(j2)/mass(j2)))
         km0 = km
*
*       Modify time-step & mass loss terms by speed-up factor.
*
         dt = km*tb
         dtm = dt/1.0d+06
*
* Take the stellar evolution timestep into account but don't let it 
* be overly restrictive for long lived phases. 
*
         if(iter.le.1000) dtm = MIN(dtm,dtmi(1),dtmi(2)) 
         dtm = MIN(dtm,tphysf-tphys)
         dt = dtm*1.0d+06
         km = dt/tb
*
* Decide between accreted mass by secondary and/or system mass loss.
*
         taum = mass(j2)/dm1*tb
         if(kstar(j2).le.2.or.kstar(j2).eq.4)then 
*
* Limit according to the thermal timescale of the secondary.
*
            dm2 = MIN(1.d0,10.d0*taum/tkh(j2))*dm1
         elseif(kstar(j2).ge.7.and.kstar(j2).le.9)then
*
* Naked helium star secondary swells up to a core helium burning star
* or SAGB star unless the primary is also a helium star.
*
            if(kstar(j1).ge.7)then
               dm2 = MIN(1.d0,10.d0*taum/tkh(j2))*dm1
            else
               dm2 = dm1
               dmchk = dm2 - 1.05d0*dms(j2)
               if(dmchk.gt.0.d0.and.dm2/mass(j2).gt.1.0d-04)then
                  kst = MIN(6,2*kstar(j2)-10)
                  if(kst.eq.4)then
                     aj(j2) = aj(j2)/tms(j2)
                     mcx = mass(j2)
                  else
                     mcx = massc(j2)
                  endif
                  mt2 = mass(j2) + km*(dm2 - dms(j2))
                  CALL gntage(mcx,mt2,kst,zpars,mass0(j2),aj(j2))
                  epoch(j2) = tphys + dtm - aj(j2)
*
               endif
            endif            
         elseif(kstar(j1).le.6.and. 
     &           (kstar(j2).ge.10.and.kstar(j2).le.12))then
*
* White dwarf secondary.
*
            if(dm1/tb.lt.2.71d-07)then
               if(dm1/tb.lt.1.03d-07)then
*
* Accrete until a nova explosion blows away most of the accreted material.
*
                  novae = .true. 
                  dm2 = MIN(dm1,dme)
                  if(dm2.lt.dm1) supedd = .true. 
                  dm22 = epsnov*dm2
               else
*
* Steady burning at the surface
*
                  dm2 = dm1
               endif
            else
*
* Make a new giant envelope.
*
               dm2 = dm1
*
* Check for planets or low-mass WDs.
*
               if((kstar(j2).eq.10.and.mass(j2).lt.0.05d0).or.
     &            (kstar(j2).ge.11.and.mass(j2).lt.0.5d0))then
                  kst = kstar(j2)
               else
                  kst = MIN(6,3*kstar(j2)-27)
                  mt2 = mass(j2) + km*(dm2 - dms(j2))
                  CALL gntage(massc(j2),mt2,kst,zpars,mass0(j2),aj(j2))
                  epoch(j2) = tphys + dtm - aj(j2)
                  change = .true. 
               endif
*
            endif
         elseif(kstar(j2).ge.10)then
*
* Impose the Eddington limit.
*
            dm2 = MIN(dm1,dme)
            if(dm2.lt.dm1) supedd = .true. 
*
         else
*
* We have a giant whose envelope can absorb any transferred material.
*
            dm2 = dm1
         endif
         if(.not.novae) dm22 = dm2
*
         if(kst.ge.10.and.kst.le.12)then
            mt2 = mass(j2) + km*(dm22 - dms(j2))
            if(kstar(j1).le.10.and.kst.eq.10.and.mt2.ge.0.7d0)then
*
* HeWD can only accrete helium-rich material up to a mass of 0.7 when
* it is destroyed in a possible Type 1a SN.
*
               mass(j1) = mass(j1) - km*(dm1 + dms(j1))
               mass(j2) = 0.d0
               kstar(j2) = 15
               goto 135
            elseif(kstar(j1).le.10.and.kst.ge.11)then
*
* CO and ONeWDs accrete helium-rich material until the accumulated 
* material exceeds a mass of 0.15 when it ignites. For a COWD with 
* mass less than 0.95 the system will be destroyed as an ELD in a 
* possible Type 1a SN. COWDs with mass greater than 0.95 and ONeWDs 
* will survive with all the material converted to ONe (JH 30/09/99). 
*
** Now changed to an ELD for all COWDs when 0.15 accreted (JH 11/01/00).  
*
               if((mt2-mass0(j2)).ge.0.15d0)then
                  if(kst.eq.11)then
                     mass(j1) = mass(j1) - km*(dm1 + dms(j1))
                     mass(j2) = 0.d0
                     kstar(j2) = 15
                     goto 135
                  endif
                  mass0(j2) = mt2
               endif
            else
               mass0(j2) = mt2
            endif
*
* If the Chandrasekhar limit is exceeded for a white dwarf then destroy
* the white dwarf in a supernova. If the WD is ONe then a neutron star
* will survive the supernova and we let HRDIAG take care of this when
* the stars are next updated.
*
            if(kst.eq.10.or.kst.eq.11)then
               if(mt2.ge.mch)then
                  dm1 = mch - mass(j2) + km*dms(j2)
                  mass(j1) = mass(j1) - dm1 - km*dms(j1)
                  mass(j2) = 0.d0
                  kstar(j2) = 15
                  goto 135
               endif
            endif
         endif
*
*       Modify mass loss terms by speed-up factor.
*
         dm1 = km*dm1
         dm2 = km*dm2
         dm22 = km*dm22
         dme = km*dme
*
* Calculate orbital angular momentum change due to system mass loss.
*
         djorb = ((dmr(1)+q(1)*dmt(1))*mass(2)*mass(2) +
     &            (dmr(2)+q(2)*dmt(2))*mass(1)*mass(1))/
     &           (mass(1)+mass(2))**2
         djorb = djorb*dt
*
* For super-Eddington mass transfer rates, for gamma = -2.0, 
* and for novae systems, assume that material is lost from  
* the system as if a wind from the secondary. 
* If gamma = -1.0 then assume the lost material carries with it 
* the specific angular momentum of the primary and for all 
* gamma > 0.0 assume that it takes away a fraction gamma of 
* the orbital angular momentum. 
*
         if(supedd.or.novae.or.gamma.lt.-1.5d0)then
            djorb = djorb + (dm1 - dm22)*mass(j1)*mass(j1)/
     &              (mass(1)+mass(2))**2
         elseif(gamma.ge.0.d0)then
            djorb = djorb + gamma*(dm1 - dm2)
         else
            djorb = djorb + (dm1 - dm2)*mass(j2)*mass(j2)/
     &              (mass(1)+mass(2))**2
         endif
*
         ecc2 = ecc*ecc
         omecc2 = 1.d0 - ecc2
         sqome2 = SQRT(omecc2)
*
         djorb = djorb*sep*sep*sqome2*oorb
         delet = 0.d0
*
* For very close systems include angular momentum loss mechanisms.
*
         if(sep.le.10.d0)then
            djgr = 8.315d-10*mass(1)*mass(2)*(mass(1)+mass(2))/
     &             (sep*sep*sep*sep)
            f1 = (19.d0/6.d0) + (121.d0/96.d0)*ecc2
            sqome5 = sqome2**5
            delet1 = djgr*ecc*f1/sqome5
            djgr = djgr*jorb*(1.d0+0.875d0*ecc2)/sqome5
            djorb = djorb + djgr*dt
            delet = delet + delet1*dt
         endif
*
         do 602 , k = 1,2
*
            dms(k) = km*dms(k)
            if(kstar(k).lt.10) dms(k) = MIN(dms(k),mass(k) - massc(k))
*
* Calculate change in the intrinsic spin of the star.
*
            djspint(k) = (2.d0/3.d0)*(dmr(k)*radx(k)*radx(k)*ospin(k) -
     &                   xi*dmt(k)*radx(3-k)*radx(3-k)*ospin(3-k))
            djspint(k) = djspint(k)*dt
*
            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
               djspint(k) = djspint(k) + djmb*dt
            endif
*
 602     continue
*
* Adjust the spin angular momentum of each star owing to mass transfer 
* and conserve total angular momentum. 
*
         djt = dm1*radx(j1)*radx(j1)*ospin(j1)
         djspint(j1) = djspint(j1) + djt
         djorb = djorb - djt
         if(disk)then
*
* Alter spin of the degenerate secondary by assuming that material
* falls onto the star from the inner edge of a Keplerian accretion
* disk and that the system is in a steady state.
*
            djt = dm2*twopi*aursun*SQRT(aursun*mass(j2)*radx(j2)) 
            djspint(j2) = djspint(j2) - djt 
            djorb = djorb + djt
*
         else
*
* No accretion disk. 
* Calculate the angular momentum of the transferred material by 
* using the radius of the disk (see Ulrich & Burger) that would 
* have formed if allowed. 
*
            rdisk = 1.7d0*rmin
            djt = dm2*twopi*aursun*SQRT(aursun*mass(j2)*rdisk) 
            djspint(j2) = djspint(j2) - djt
            djorb = djorb + djt
*
         endif
         djtx(2) = djt
*
* Adjust the secondary spin if a nova eruption has occurred. 
*
         if(novae)then
            djt = (dm2 - dm22)*radx(j2)*radx(j2)*ospin(j2) 
            djspint(j2) = djspint(j2) + djt 
            djtx(2) = djtx(2) - djt
         endif
*
* Calculate circularization, orbital shrinkage and spin up.
*
         do 603 , k = 1,2
*
            dspint(k) = 0.d0
            if(((kstar(k).le.9.and.rad(k).ge.0.01d0*rol(k)).or.
     &         (kstar(k).ge.10.and.k.eq.j1)).and.tflag.gt.0)then
*
               raa2 = (radx(k)/sep)**2
               raa6 = raa2**3
*
               f5 = 1.d0+ecc2*(3.d0+ecc2*0.375d0)
               f4 = 1.d0+ecc2*(1.5d0+ecc2*0.125d0)
               f3 = 1.d0+ecc2*(3.75d0+ecc2*(1.875d0+ecc2*7.8125d-02))
               f2 = 1.d0+ecc2*(7.5d0+ecc2*(5.625d0+ecc2*0.3125d0))
               f1 = 1.d0+ecc2*(15.5d0+ecc2*(31.875d0+ecc2*(11.5625d0
     &                  +ecc2*0.390625d0)))
*
               if((kstar(k).eq.1.and.mass(k).ge.1.25d0).or.
     &            kstar(k).eq.4.or.kstar(k).eq.7)then
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*radx(k)*radx(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k)
               elseif(kstar(k).le.9)then
                  renv(k) = MIN(renv(k),radx(k)-radc(k))
                  renv(k) = MAX(renv(k),1.0d-10)
                  tc = mr23yr*(menv(k)*renv(k)*(radx(k)-0.5d0*renv(k))/
     &                 (3.d0*lumin(k)))**(1.d0/3.d0)
                  ttid = twopi/(1.0d-10 + ABS(oorb - ospin(k)))
                  f = MIN(1.d0,(ttid/(2.d0*tc)**2))
                  tcqr = 2.d0*f*q(3-k)*raa6*menv(k)/(21.d0*tc*mass(k))
                  rg2 = (k2str(k)*(mass(k)-massc(k)))/mass(k)
               else
                  f = 7.33d-09*(lumin(k)/mass(k))**(5.d0/7.d0)
                  tcqr = f*q(3-k)*q(3-k)*raa2*raa2/(1.d0+q(3-k))
                  rg2 = k3
               endif
               sqome3 = sqome2**3
               delet1 = 27.d0*tcqr*(1.d0+q(3-k))*raa2*(ecc/sqome2**13)*
     &                  (f3 - (11.d0/18.d0)*sqome3*f4*ospin(k)/oorb)
               tcirc = ecc/(ABS(delet1) + 1.0d-20)
               delet = delet + delet1*dt
               dspint(k) = (3.d0*q(3-k)*tcqr/(rg2*omecc2**6))*
     &                     (f2*oorb - sqome3*f5*ospin(k))
               eqspin = oorb*f2/(sqome3*f5)
               if(dt.gt.0.d0)then
                  if(dspint(k).ge.0.d0)then
                     dspint(k) = MIN(dt*dspint(k),eqspin-ospin(k))/dt
                  else
                     dspint(k) = MAX(dt*dspint(k),eqspin-ospin(k))/dt
                  endif
               endif
               djt = (k2str(k)*(mass(k)-massc(k))*radx(k)*radx(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               djorb = djorb + djt*dt
               djspint(k) = djspint(k) - djt*dt
*
            endif
*
            jspin(k) = MAX(1.0d-10,jspin(k) - djspint(k))
*
* Ensure that the star does not spin up beyond break-up, and transfer
* the excess angular momentum back to the orbit.
*
            ospbru = twopi*SQRT(mass(k)*aursun**3/radx(k)**3)
            jspbru = (k2str(k)*(mass(k)-massc(k))*radx(k)*radx(k) +
     &                k3*massc(k)*radc(k)*radc(k))*ospbru
            if(jspin(k).gt.jspbru)then
               mew = 1.d0
               if(djtx(2).gt.0.d0)then
                  mew = MIN(mew,(jspin(k) - jspbru)/djtx(2))
               endif
               djorb = djorb - (jspin(k) - jspbru)
               jspin(k) = jspbru
* If excess material should not be accreted, activate next line.
*              dm22 = (1.d0 - mew)*dm22
            endif
*
 603     continue
*
* Update the masses.
*
         kstar(j2) = kst
         mass(j1) = mass(j1) - dm1 - dms(j1)
         if(kstar(j1).le.1.or.kstar(j1).eq.7) mass0(j1) = mass(j1)
         mass(j2) = mass(j2) + dm22 - dms(j2)
         if(kstar(j2).le.1.or.kstar(j2).eq.7) mass0(j2) = mass(j2)
*
         dm1 = dm1 + dms(j1)
         if(dtm.gt.tiny)then
            if(dm1.gt.tiny)then
               dtmx = MIN(dtmx,dmmax*mass(j1)*dtm/dm1)
            endif
            dmdt(j1) = -1.d0*dm1/dt
            dmdt(j2) = (dm22-dms(j2))/dt
         endif
*
* For a HG star check if the initial mass can be reduced. 
*
         if(kstar(j1).eq.2.and.mass0(j1).le.zpars(3))then
            m0 = mass0(j1)
            mass0(j1) = mass(j1)
            CALL star(kstar(j1),mass0(j1),mass(j1),tmsnew,tn,tscls,
     &                lums,GB,zpars)
            if(GB(9).lt.massc(j1))then
               mass0(j1) = m0
            endif
         endif
         if(kstar(j2).eq.2.and.mass0(j2).le.zpars(3))then
            m0 = mass0(j2)
            mass0(j2) = mass(j2)
            CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &                lums,GB,zpars)
            if(GB(9).lt.massc(j2))then
               mass0(j2) = m0
            endif
         endif
*
         ecc = ecc - delet
         ecc = MAX(ecc,0.d0)
         if(ecc.lt.1.0d-10) ecc = 0.d0
*
         if(ecc.ge.1.d0) goto 135
*
* Ensure that Jorb does not become negative which could happen if the 
* primary overfills its Roche lobe initially. In this case we simply 
* allow contact to occur.
*
         jorb = MAX(1.d0,jorb - djorb)
         sep = (mass(1) + mass(2))*jorb*jorb/
     &         ((mass(1)*mass(2)*twopi)**2*aursun**3*(1.d0-ecc*ecc))
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
*
      endif
*
* Always rejuvenate the secondary and age the primary if they are on
* the main sequence.
*
      if(kstar(j1).le.2.or.kstar(j1).eq.7)then
         CALL star(kstar(j1),mass0(j1),mass(j1),tmsnew,tn,tscls,
     &             lums,GB,zpars)
         if(kstar(j1).eq.2)then
            aj(j1) = tmsnew + (tscls(1) - tmsnew)*(aj(j1)-tms(j1))/
     &                        (tbgb(j1) - tms(j1))
         else
            aj(j1) = tmsnew/tms(j1)*aj(j1)
         endif
         epoch(j1) = tphys - aj(j1)
      endif
*
      if(kstar(j2).le.2.or.kstar(j2).eq.7)then
         CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &             lums,GB,zpars)
         if(kstar(j2).eq.2)then
            aj(j2) = tmsnew + (tscls(1) - tmsnew)*(aj(j2)-tms(j2))/
     &                        (tbgb(j2) - tms(j2))
         elseif((mass(j2).lt.0.35d0.or.mass(j2).gt.1.25d0).
     &           and.kstar(j2).ne.7)then
            aj(j2) = tmsnew/tms(j2)*aj(j2)*(mass(j2) - dm22)/mass(j2)
         else
            aj(j2) = tmsnew/tms(j2)*aj(j2)
         endif
         epoch(j2) = tphys - aj(j2)
      endif
*
* Obtain the stellar parameters for the next step.
*
 85   tphys = tphys + dtm
      do 90 , k = 1,2
         age = tphys - epoch(k)
         m0 = mass0(k)
         mt = mass(k)
         mc = massc(k)
*
         kw = kstar(k)
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(m0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kw,mc,rc,me,re,k2)
*
* Check for a supernova and correct the semi-major axis if so.
*
         if(kw.ne.kstar(k).and.kstar(k).le.12.and.
     &      (kw.eq.13.or.kw.eq.14))then
            dms(k) = mass(k) - mt
            print*,'kick 5a- kw,k,mf,m1,m2,ecc,sepf,jorb,vs=',
     &      kw,k,mass(k),mt,mass(3-k),ecc,sep,jorb,vs(1),vs(2),vs(3)
            CALL kickv(kw,mass(k),mt,mass(3-k),ecc,sep,jorb,vs)
            print*,'kick 5b- kw,k,mf,m1,m2,ecc,sepf,jorb,vs=',
     &      kw,k,mass(k),mt,mass(3-k),ecc,sep,jorb,vs(1),vs(2),vs(3)
            print*,'called kickv 2 in evolv2b, ecc',ecc
            do kk = 1,3
               vd(kk) = vd(kk) + vs(kk)
            enddo
            if(ecc.gt.1.d0)then
               ospin(k) = 2.0d+08
               jspin(k) = k3*rc*rc*mc*ospin(k)
               kstar(k) = kw
               mass(k) = mt
               epoch(k) = tphys - age
               goto 135
            endif
            tb = (sep/aursun)*SQRT(sep/(aursun*(mt+mass(3-k))))
            oorb = twopi/tb
         endif
         if(kw.ne.kstar(k))then
            change = .true.
            mass(k) = mt
            if(kw.eq.15)then
               kstar(k) = kw
               goto 135
            endif
            mass0(k) = m0
            epoch(k) = tphys - age
            if(kw.eq.13.or.kw.eq.14)then
               ospin(k) = 2.0d+08
               jspin(k) = k3*rc*rc*mc*ospin(k)
            endif
         endif
*
*     Determine stellar evolution timescale for nuclear burning types.
*
         if(kw.le.9)then
            CALL deltat(kw,age,tm,tn,tscls,dt,dtr)
            dtmi(k) = MIN(dt,dtr)
*           dtmi(k) = dtr
            dtmi(k) = MAX(1.0d-07,dtmi(k))
         else
            dtmi(k) = 1.0d+10
         endif
*        dtmi(k) = MAX((tn-age),1.0d-07)
*
* Save relevent solar quantities.
*
         aj(k) = age
         kstar(k) = kw
         rad(k) = rm
         radx(k) = rm
         lumin(k) = lum
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
*
 90   continue
*
      do 100 , k = 1,2
         q(k) = mass(k)/mass(3-k)
         rol(k) = rl(q(k))*sep
 100  continue
      if(rad(j1).gt.rol(j1)) radx(j1) = rol(j1)
      do 110 , k = 1,2
         ospin(k) = jspin(k)/(k2str(k)*(mass(k)-massc(k))*radx(k)*
     &              radx(k) + k3*massc(k)*radc(k)*radc(k))
 110  continue
*
      if(tphys.ge.tphysf) goto 140
*
* Test whether the primary still fills its Roche lobe.
*
      if(rad(j1).gt.rol(j1))then
*
* Test for a contact system
*
         if(rad(j2).gt.rol(j2)) goto 130
         if(check.and.change) goto 140
         iter = iter + 1
         goto 8
      else
         dtm = 0.d0
         change = .true.
         goto 4
      endif
*
 130  continue
*
* Contact system.
*
      coel = .true.
      change = .true.
*
* If *1 or *2 is giant-like this will be common-envelope evolution.
*
      m1ce = mass(j1)
      m2ce = mass(j2)
*
      if(kstar(j1).ge.2.and.kstar(j1).le.9.and.kstar(j1).ne.7)then
         CALL comenv(mass0(j1),mass(j1),massc(j1),aj(j1),jspin(j1),
     &               kstar(j1),mass0(j2),mass(j2),massc(j2),aj(j2),
     &               jspin(j2),kstar(j2),zpars,ecc,sep,jorb,coel)
         com = .true.
      elseif(kstar(j2).ge.2.and.kstar(j2).le.9.and.kstar(j2).ne.7)then
         CALL comenv(mass0(j2),mass(j2),massc(j2),aj(j2),jspin(j2),
     &               kstar(j2),mass0(j1),mass(j1),massc(j1),aj(j1),
     &               jspin(j1),kstar(j1),zpars,ecc,sep,jorb,coel)
         com = .true.
      else
         CALL mix(mass0,mass,aj,kstar,zpars)
      endif
      epoch(1) = tphys - aj(1)
      epoch(2) = tphys - aj(2)
      if(.not.coel)then
*
* Next step should be made without changing the time.
*
         if(ecc.gt.1.d0)then
            if(kstar(1).ge.13)then
               rc = corerd(kstar(1),mass(1),mass(1),zpars(2))
               ospin(1) = jspin(1)/(k3*rc*rc*mass(1))
            endif
            if(kstar(2).ge.13)then
               rc = corerd(kstar(2),mass(2),mass(2),zpars(2))
               ospin(2) = jspin(2)/(k3*rc*rc*mass(2))
            endif
            goto 135
         endif
         dtm = 0.d0
*
* Reset orbital parameters as separation may have changed.
*
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
         goto 4
      endif
*
 135  continue
*
      change = .true.
      sgl = .true.
      if(kstar(1).eq.15) rad(1) = 0.d0
      if(kstar(2).eq.15) rad(2) = 0.d0
      dmdt(1) = 0.d0
      dmdt(2) = 0.d0
*
      if(kstar(1).ne.15.or.kstar(2).ne.15)then
         if(com)then
            com = .false.
         endif
         if(kstar(2).eq.15)then
            kmax = 1
            dtmi(2) = tphysf
         elseif(kstar(1).eq.15)then
            kmin = 2
            dtmi(1) = tphysf
         endif
         ecc = -1.d0
         sep = 0.d0
         dtm = 0.d0
         coel = .false.
         goto 4
      endif
*
 140  continue
*
      tphys = MIN(tphys,tphysf)
      kw = MAX(kstar(1),kstar(2))
      if(sgl.or.ecc.lt.0.d0.or.ecc.ge.1.d0.or.kw.eq.15)then
         ecc = -1.d0
         tb = -1.d0
      else
         tb = tb*yeardy
      endif
*
* Set next timestep. 
*
      dtm = 1.0d+03
      if(kstar(1).ne.15) dtm = MIN(dtm,dtmi(1))
      if(kstar(2).ne.15) dtm = MIN(dtm,dtmi(2))
      if(tb.gt.0.d0.and.dtmx.gt.tiny) dtm = MIN(dtm,dtmx)
      dtm = MAX(dtm,1.0d-07*tphys)
      tphysf = tphys + dtm
*
      RETURN
      END
***
***
      SUBROUTINE gntage(mc,mt,kw,zpars,m0,aj)
*
* A routine to determine the age of a giant from its core mass and type.
*
*       Author : C. A. Tout
*       Date   : 24th September 1996
*       Revised: 21st February 1997 to include core-helium-burning stars
*
*       Rewritten: 2nd January 1998 by J. R. Hurley to be compatible with
*                  the new evolution routines and to include new stellar
*                  types.
*
      implicit none
*
      integer kw
      integer j,jmax
      parameter(jmax=30)
*
      real*8 mc,mt,m0,aj,tm,tn
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 mmin,mmax,mmid,dm,f,fmid,dell,derl,lum
      real*8 macc,lacc,tiny
      parameter(macc=0.00001d0,lacc=0.0001d0,tiny=1.0d-14)
      real*8 mcx,mcy
*
      real*8 mcheif,mcagbf,mheif,mbagbf,mcgbf,lmcgbf,lbgbf,lbgbdf
      external mcheif,mcagbf,mheif,mbagbf,mcgbf,lmcgbf,lbgbf,lbgbdf
*
* This should only be entered with KW = 3, 4, 5, 6 or 9
*
* First we check that we don't have a CheB star 
* with too small a core mass.
      if(kw.eq.4)then
* Set the minimum CHeB core mass using M = Mflash
         mcy = mcheif(zpars(2),zpars(2),zpars(10))
         if(mc.le.mcy) kw = 3
*        if(mc.le.mcy) WRITE(66,*)' GNTAGE4: changed to 3'
      endif
*
* Next we check that we don't have a GB star for M => Mfgb
      if(kw.eq.3)then
* Set the maximum GB core mass using M = Mfgb
         mcy = mcheif(zpars(3),zpars(2),zpars(9))
         if(mc.ge.mcy)then
            kw = 4
            aj = 0.d0
*           WRITE(66,*)' GNTAGE3: changed to 4'
         endif
      endif
*
      if(kw.eq.6)then
*
* We try to start the star from the start of the SAGB by
* setting Mc = Mc,TP.
*
         mcy = 0.44d0*2.25d0 + 0.448d0
         if(mc.gt.mcy)then
* A type 6 with this sized core mass cannot exist as it should
* already have become a NS or BH as a type 5. 
* We set it up so that it will.
            mcx = (mc + 0.35d0)/0.773d0
         elseif(mc.ge.0.8d0)then
            mcx = (mc - 0.448d0)/0.44d0
         else
            mcx = mc
         endif
         m0 = mbagbf(mcx)
         if(m0.lt.tiny)then
* Carbon core mass is less then the minimum for the start of SAGB.
* This must be the case of a low-mass C/O or O/Ne WD with only a
* very small envelope added or possibly the merger of a helium star
* with a main sequence star. We will set m0 = mt and then reset the
* core mass to allow for some helium to be added to the C/O core.
            kw = 14
*           WRITE(66,*)' GNTAGE6: changed to 4'
         else
            CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(13) + 2.d0*tiny
         endif
      endif
*
      if(kw.eq.5)then
*
* We fit a Helium core mass at the base of the AGB.
*
         m0 = mbagbf(mc)
         if(m0.lt.tiny)then
* Helium core mass is less then the BAGB minimum.
            kw = 14
*           WRITE(66,*)' GNTAGE5: changed to 4'
         else
            CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(2) + tscls(3) + 2.d0*tiny
         endif
      endif
*
*
      if(kw.eq.4)then
*
* The supplied age is actually the fractional age, fage, of CHeB lifetime
* that has been completed, ie. 0 <= aj <= 1.
*
         if(aj.lt.0.d0.or.aj.gt.1.d0)then
*           WRITE(99,*)' FATAL ERROR! GNTAGE4: fage out of bounds '
*           WRITE(99,*)' FAGE ',aj
*           WRITE(*,*)' STOP: FATAL ERROR '
*           CALL exit(0)
*           STOP
            aj = 0.d0
         endif
* Get the minimum, fage=1, and maximum, fage=0, allowable masses
         mcy = mcagbf(zpars(2))
         if(mc.ge.mcy)then
            mmin = mbagbf(mc)
         else
            mmin = zpars(2)
         endif
         mmax = mheif(mc,zpars(2),zpars(10))
         if(aj.lt.tiny)then
            m0 = mmax
            goto 20
         elseif(aj.ge.1.d0)then
            m0 = mmin
            goto 20
         endif
* Use the bisection method to find m0
         fmid = (1.d0-aj)*mcheif(mmax,zpars(2),zpars(10)) +
     &          aj*mcagbf(mmax) - mc
         f = (1.d0-aj)*mcheif(mmin,zpars(2),zpars(10)) +
     &       aj*mcagbf(mmin) - mc
         if(f*fmid.ge.0.d0)then
* This will probably occur if mc is just greater than the minimum
* allowed mass for a CHeB star and fage > 0.
            kw = 3
*           WRITE(66,*)' GNTAGE4: changed to 3'
            goto 90
         endif
         m0 = mmin
         dm = mmax - mmin
         do 10 , j = 1,jmax
            dm = 0.5d0*dm
            mmid = m0 + dm
            fmid = (1.d0-aj)*mcheif(mmid,zpars(2),zpars(10)) +
     &             aj*mcagbf(mmid) - mc
            if(fmid.lt.0.d0) m0 = mmid
            if(ABS(dm).lt.macc.or.ABS(fmid).lt.tiny) goto 20
            if(j.eq.jmax)then
*              WRITE(99,*)' FATAL ERROR! GNTAGE4: root not found '
*              WRITE(*,*)' STOP: FATAL ERROR '
*              CALL exit(0)
*              STOP
               m0 = mt
               aj = 0.d0
            endif
 10      continue
 20      continue
*
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         aj = tscls(2) + aj*tscls(3)
*
      endif
*
 90   continue
*
      if(kw.eq.3)then
*
* First we double check that we don't have a GB star for M => Mfgb
         mcy = mcheif(zpars(3),zpars(2),zpars(9))
         if(mc.ge.mcy)then
*           WRITE(99,*)' GNTAGE3: star too big for GB '
*           WRITE(*,*)' STOP: FATAL ERROR '
*           CALL exit(0)
*           STOP
            mc = 0.99d0*mcy
         endif
* Next we find an m0 so as to place the star at the BGB
         mcx = mcheif(zpars(2),zpars(2),zpars(9))
         if(mc.gt.mcx)then
            m0 = mheif(mc,zpars(2),zpars(9))
         else
* Use Newton-Raphson to find m0 from Lbgb
            m0 = zpars(2)
            CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
            lum = lmcgbf(mc,GB)
            j = 0
 30         continue
            dell = lbgbf(m0) - lum
            if(ABS(dell/lum).le.lacc) goto 40
            derl = lbgbdf(m0)
            m0 = m0 - dell/derl
            j = j + 1
            if(j.eq.jmax)then
*              WRITE(99,*)' FATAL ERROR! GNTAGE3: root not found '
*              WRITE(*,*)' STOP: FATAL ERROR '
*              CALL exit(0)
*              STOP
               m0 = zpars(2)
               m0 = MAX(m0,mt)
               goto 40
            endif
            goto 30
 40         continue
         endif
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         aj = tscls(1) + 1.0d-06*(tscls(2) - tscls(1))
*
      endif
*
      if(kw.eq.8.or.kw.eq.9)then
*
* We make a post-MS naked helium star.
* To make things easier we put the star at the TMS point
* so it actually begins as type 8.
*
         kw = 8
         mmin = mc
         CALL star(kw,mmin,mc,tm,tn,tscls,lums,GB,zpars)
         mcx = mcgbf(lums(2),GB,lums(6))
         if(mcx.ge.mc)then
*           WRITE(99,*)' FATAL ERROR! GNTAGE9: mmin too big '
*           WRITE(*,*)' STOP: FATAL ERROR '
*           CALL exit(0)
*           STOP
            m0 = mt
            goto 80
         endif
         f = mcx - mc
         mmax = mt
         do 50 , j = 1,jmax
            CALL star(kw,mmax,mc,tm,tn,tscls,lums,GB,zpars)
            mcy = mcgbf(lums(2),GB,lums(6))
            if(mcy.gt.mc) goto 60
            mmax = 2.d0*mmax
            if(j.eq.jmax)then
*              WRITE(99,*)' FATAL ERROR! GNTAGE9: mmax not found '
*              WRITE(*,*)' STOP: FATAL ERROR '
*              CALL exit(0)
*              STOP
               m0 = mt
               goto 80
            endif
 50      continue
 60      continue
         fmid = mcy - mc
* Use the bisection method to find m0
         if(f*fmid.ge.0.d0)then
*           WRITE(99,*)' FATAL ERROR! GNTAGE9: root not bracketed '
*           WRITE(*,*)' STOP: FATAL ERROR '
*           CALL exit(0)
*           STOP
            m0 = mt
            goto 80
         endif
         m0 = mmin
         dm = mmax - mmin
         do 70 , j = 1,jmax
            dm = 0.5d0*dm
            mmid = m0 + dm
            CALL star(kw,mmid,mc,tm,tn,tscls,lums,GB,zpars)
            mcy = mcgbf(lums(2),GB,lums(6))
            fmid = mcy - mc
            if(fmid.lt.0.d0) m0 = mmid
            if(ABS(dm).lt.macc.or.ABS(fmid).lt.tiny) goto 80
            if(j.eq.jmax)then
*              WRITE(99,*)' FATAL ERROR! GNTAGE9: root not found '
*              WRITE(*,*)' STOP: FATAL ERROR '
*              CALL exit(0)
*              STOP
               m0 = mt
               goto 80
            endif
 70      continue
 80      continue
*
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         aj = tm + 1.0d-10*tm
*
      endif
*
      if(kw.eq.14)then
*
         kw = 4
         m0 = mt
         mcy = mcagbf(m0)
         aj = mc/mcy
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         if(m0.le.zpars(2))then
            mcx = mcgbf(lums(4),GB,lums(6))
         else
            mcx = mcheif(m0,zpars(2),zpars(10))
         end if
         mc = mcx + (mcy - mcx)*aj
         aj = tscls(2) + aj*tscls(3)
      endif
*
      RETURN
      END
***
***
      SUBROUTINE hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2)
*
*
*       H-R diagram for population I stars.
*       -----------------------------------
*
*       Computes the new mass, luminosity, radius & stellar type.
*       Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR.
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout;
*       24th October 1995 to include metallicity;
*       14th November 1996 to include naked helium stars;
*       28th February 1997 to allow accretion induced supernovae.
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include Z=0.001 as well as Z=0.02, convective overshooting,
*       MS hook and more elaborate CHeB
*
      implicit none
*
      integer kw,kwp
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
*
      real*8 mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,rc,menv,renv,k2
      real*8 mch,mlp,tiny
      parameter(mch=1.44d0,mlp=12.d0,tiny=1.0d-14)
      real*8 mass0,mt0,mtc
      REAL*8 neta,bwind,mxns
      COMMON /VALUE1/ neta,bwind,mxns
* 
      real*8 thook,thg,tbagb,tau,tloop,taul,tauh,tau1,tau2,dtau,texp
      real*8 lx,ly,dell,alpha,beta,eta
      real*8 rx,ry,delr,rzams,rtms,gamma,rmin,taumin,rg
      parameter(taumin=5.0d-08)
      real*8 mcmax,mcx,mcy,mcbagb,lambda
      real*8 am,xx,fac,rdgen,mew,lum0,kap,zeta,ahe,aco
      parameter(lum0=7.0d+04,kap=-0.5d0,ahe=4.d0,aco=16.d0)
*
      real*8 thookf,tblf
      real*8 lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      real*8 rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      real*8 rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      real*8 mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
      external thookf,tblf
      external lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      external rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      external rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      external mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
*
*
*       ---------------------------------------------------------------------
*       MASS    Stellar mass in solar units (input: old; output: new value).
*       AJ      Current age in Myr.
*       MT      Current mass in solar units (used for R).
*       TM      Main sequence time.
*       TN      Nuclear burning time.
*       TSCLS   Time scale for different stages.
*       LUMS    Characteristic luminosity.
*       GB      Giant Branch parameters
*       ZPARS   Parameters for distinguishing various mass intervals.
*       R       Stellar radius in solar units.
*       TE      Effective temperature (suppressed).
*       KW      Classification type (0 - 15).
*       MC      Core mass.
*       ---------------------------------------------------------------------
*
*
* Make evolutionary changes to stars that have not reached KW > 5.
*
      mass0 = mass
      if(mass0.gt.100.d0) mass = 100.d0
      mt0 = mt
      if(mt0.gt.100.d0) mt = 100.d0
*
c      print*,'kw',kw
      if(kw.gt.6) goto 90
*
      tbagb = tscls(2) + tscls(3)
      thg = tscls(1) - tm
*
      rzams = rzamsf(mass)
      rtms = rtmsf(mass)
*
      if(aj.lt.tscls(1))then
*
*        Either on MS or HG
*
         rg = rgbf(mt,lums(3))
*
         if(aj.lt.tm)then
*
*           Main sequence star.
*
            mc = 0.d0
            tau = aj/tm
            thook = thookf(mass)*tscls(1)
            zeta = 0.01d0
            tau1 = MIN(1.d0,aj/thook)
            tau2 = MAX(0.d0,
     &             MIN(1.d0,(aj-(1.d0-zeta)*thook)/(zeta*thook)))
*
            dell = lhookf(mass,zpars(1))
            dtau = tau1**2 - tau2**2
            alpha = lalphf(mass)
            beta = lbetaf(mass)
            eta = lnetaf(mass)
            lx = LOG10(lums(2)/lums(1))
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**eta +
     &              (lx - alpha - beta)*tau**2 - dell*dtau
            else
               xx = alpha*tau + (lx - alpha)*tau**2 - dell*dtau
            endif
            lum = lums(1)*10.d0**xx
*
            delr = rhookf(mass,zpars(1))
            dtau = tau1**3 - tau2**3
            alpha = ralphf(mass)
            beta = rbetaf(mass)
            gamma = rgammf(mass)
            rx = LOG10(rtms/rzams)
* Note that the use of taumin is a slightly pedantic attempt to
* avoid floating point underflow. It IS overkill!
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**10 + gamma*tau**40 +
     &              (rx - alpha - beta - gamma)*tau**3 - delr*dtau
            else
               xx = alpha*tau + (rx - alpha)*tau**3 - delr*dtau
            endif
            r = rzams*10.d0**xx
*
            if(mass.lt.(zpars(1)-0.3d0))then
               kw = 0
* This following is given by Chris for low mass MS stars which will be 
* substantially degenerate. We need the Hydrogen abundance, X, which we 
* calculate from Z assuming that the helium abundance, Y, is calculated
* according to Y = 0.24 + 2*Z
               rdgen = 0.0258d0*((1.d0+zpars(11))**(5.d0/3.d0))*
     &                          (mass**(-1.d0/3.d0))
               r = MAX(rdgen,r)
            else
               kw = 1
            endif
*
         else 
*
*           Star is on the HG
*
            mcx = mc
            if(mass.le.zpars(2))then
               mc = mcgbf(lums(3),GB,lums(6))
            elseif(mass.le.zpars(3))then
               mc = mcheif(mass,zpars(2),zpars(9))
            else
               mc = mcheif(mass,zpars(2),zpars(10))
            endif
            eta = mctmsf(mass)
            tau = (aj - tm)/thg
            mc = ((1.d0 - tau)*eta + tau)*mc
            mc = MAX(mc,mcx)
*
* Test whether core mass has reached total mass.
*
            if(mc.ge.mt)then
               aj = 0.d0
               if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
                  mc = 0.d0
                  mass = mt
                  kw = 7
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               else
*
* Zero-age helium white dwarf.
*
                  mc = mt
                  mass = mt
                  kw = 10
               endif
            else
               lum = lums(2)*(lums(3)/lums(2))**tau
               if(mass.le.zpars(3))then
                  rx = rg
               else
* He-ignition and end of HG occur at Rmin
                  rmin = rminf(mass)
                  ry = ragbf(mt,lums(4),zpars(2))
                  rx = MIN(rmin,ry)
                  if(mass.le.mlp)then
                     texp = log(mass/mlp)/log(zpars(3)/mlp)
                     rx = rg
                     rx = rmin*(rx/rmin)**texp
                  endif
                  tau2 = tblf(mass,zpars(2),zpars(3))
                  if(tau2.lt.tiny) rx = ry
               endif
               r = rtms*(rx/rtms)**tau
               kw = 2
            endif
*
         endif
*
* Now the GB, CHeB and AGB evolution.
*
      elseif(aj.lt.tscls(2))then
*
*        Red Giant.
*
         kw = 3
         lum = lgbtf(aj,GB(1),GB,tscls(4),tscls(5),tscls(6))
         if(mass.le.zpars(2))then
* Star has a degenerate He core which grows on the GB
            mc = mcgbf(lum,GB,lums(6))
         else
* Star has a non-degenerate He core which may grow, but
* only slightly, on the GB
            tau = (aj - tscls(1))/(tscls(2) - tscls(1))
            mcx = mcheif(mass,zpars(2),zpars(9))
            mcy = mcheif(mass,zpars(2),zpars(10))
            mc = mcx + (mcy - mcx)*tau
         endif
         r = rgbf(mt,lum)
         rg = r
         if(mc.ge.mt)then
            aj = 0.d0
            if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
               mc = 0.d0
               mass = mt
               kw = 7
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            else
*
* Zero-age helium white dwarf.
*
               mc = mt
               mass = mt
               kw = 10
            endif
         endif
*
      elseif(aj.lt.tbagb)then
*
*       Core helium burning star.
*
         if(kw.eq.3.and.mass.le.zpars(2))then
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(2)
         endif
         if(mass.le.zpars(2))then
            mcx = mcgbf(lums(4),GB,lums(6))
         else
            mcx = mcheif(mass,zpars(2),zpars(10))
         end if
         tau = (aj - tscls(2))/tscls(3)
         mc = mcx + (mcagbf(mass) - mcx)*tau
*
         if(mass.le.zpars(2))then
            lx = lums(5)
            ly = lums(7)
            rx = rzahbf(mt,mc,zpars(2))
            rg = rgbf(mt,lx)
            rmin = rg*zpars(13)**(mass/zpars(2))
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            ry = ragbf(mt,ly,zpars(2))
            if(rmin.lt.rx)then
               taul = (log(rx/rmin))**(1.d0/3.d0)
            else
               rmin = rx
               taul = 0.d0
            endif
            tauh = (log(ry/rmin))**(1.d0/3.d0)
            tau2 = taul*(tau - 1.d0) + tauh*tau
            r = rmin*exp(abs(tau2)**3)
            rg = rg + tau*(ry - rg)
            lum = lx*(ly/lx)**(tau**texp)
         elseif(mass.gt.zpars(3))then
*
* For HM stars He-ignition takes place at Rmin in the HG, and CHeB
* consists of a blue phase (before tloop) and a RG phase (after tloop).
*
            tau2 = tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            rmin = rminf(mass)
            rg = rgbf(mt,lums(4))
            rx = ragbf(mt,lums(4),zpars(2))
            rmin = MIN(rmin, rx)
            if(mass.le.mlp) then
               texp = log(mass/mlp)/log(zpars(3)/mlp)
               rx = rg
               rx = rmin*(rx/rmin)**texp
            else
               rx = rmin
            end if
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            lum = lums(4)*(lums(7)/lums(4))**(tau**texp)
            if(aj.lt.tloop)then
               ly = lums(4)*(lums(7)/lums(4))**(tau2**texp)
               ry = ragbf(mt,ly,zpars(2))
               taul = 0.d0
               if(ABS(rmin-rx).gt.tiny)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               endif
               tauh = 0.d0
               if(ry.gt.rmin) tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tscls(2))/(tau2*tscls(3))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rg + tau*(ry - rg)
            else
               r = ragbf(mt,lum,zpars(2))
               rg = r
            end if
         else
*
* For IM stars CHeB consists of a RG phase (before tloop) and a blue
* loop (after tloop).
*
            tau2 = 1.d0 - tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            if(aj.lt.tloop)then
               tau = (tloop - aj)/(tau2*tscls(3))
               lum = lums(5)*(lums(4)/lums(5))**(tau**3)
               r = rgbf(mt,lum)
               rg = r
            else
               lx = lums(5)
               ly = lums(7)
               rx = rgbf(mt,lx)
               rmin = rminf(mt)
               texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
               ry = ragbf(mt,ly,zpars(2))
               if(rmin.lt.rx)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               else
                  rmin = rx
                  taul = 0.d0
               endif
               tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tloop)/(tscls(3) - (tloop - tscls(2)))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rx + tau*(ry - rx)
               lum = lx*(ly/lx)**(tau**texp)
            endif
         endif
* 
* Test whether core mass exceeds total mass.
*
         if(mc.ge.mt)then
*
* Evolved MS naked helium star.
*
            kw = 7
            xx = (aj - tscls(2))/tscls(3)
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = xx*tm
         else
            kw = 4
         endif
*
      else
*
*        Asymptotic Red Giant.
*
* On the AGB the He core mass remains constant until at Ltp it
* is caught by the C core mass and they grow together.
*
         mcbagb = mcagbf(mass)
         mcx = mcgbtf(tbagb,GB(8),GB,tscls(7),tscls(8),tscls(9))
         mcmax = MAX(MAX(mch,0.773d0*mcbagb-0.35d0),1.05d0*mcx)
*
         if(aj.lt.tscls(13))then
            mcx = mcgbtf(aj,GB(8),GB,tscls(7),tscls(8),tscls(9))
            mc = mcbagb
            lum = lmcgbf(mcx,GB)
            if(mt.le.mc)then
*
* Evolved naked helium star as the envelope is lost but the
* star has not completed its interior burning. The star becomes
* a post-HeMS star.
*
               kw = 9
               mt = mc
               mass = mt
               mc = mcx
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(mc.le.GB(7))then
                  aj = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                            (mc**(1.d0-GB(5)))
               else
                  aj = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                            (mc**(1.d0-GB(6)))
               endif
               aj = MAX(aj,tm)
               goto 90
            else
               kw = 5
            endif
         else
c            print*,'set kw = 6'
            kw = 6
            mc = mcgbtf(aj,GB(2),GB,tscls(10),tscls(11),tscls(12))
            lum = lmcgbf(mc,GB)
c            print*, lum,mc,GB
*
* Approximate 3rd Dredge-up on AGB by limiting Mc.
*
            lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
            tau = tscls(13)
            mcx = mcgbtf(tau,GB(2),GB,tscls(10),tscls(11),tscls(12))
            mcy = mc
            mc = mc - lambda*(mcy-mcx)
            mcx = mc
            mcmax = MIN(mt,mcmax)   
c            print*,'kw,mc,lum,lambda,tau,mcx,mcy,mc,mcx,mcmax',kw,mc,
c     &           lum,lambda,tau,mcx,mcy,mc,mcx,mcmax
         endif
         r = ragbf(mt,lum,zpars(2))
c         print*,r,mt,lum,zpars(2)
         rg = r
*
* Mc,x represents the C core mass and we now test whether it
* exceeds either the total mass or the maximum allowed core mass.
*
         if(mcx.ge.mcmax)then
            aj = 0.d0
            mc = mcmax
            if(mc.lt.mch)then
               if(ifflag.ge.1)then
*
* Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800. 
*
                  if(zpars(14).ge.1.0d-08)then
                     mc = MIN(0.36d0+0.104d0*mass,0.58d0+0.061d0*mass)
                     mc = MAX(0.54d0+0.042d0*mass,mc)
                     if(mass.lt.1.d0) mc = 0.46d0
                  else
                     mc = MIN(0.29d0+0.178d0*mass,0.65d0+0.062d0*mass)
                     mc = MAX(0.54d0+0.073d0*mass,mc)
                  endif
                  mc = MIN(mch,mc)
               endif
*
               mt = mc
               if(mcbagb.lt.1.6d0)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                  kw = 11
               else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                  kw = 12
               endif
               mass = mt
*
            else
               if(mcbagb.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                  kw = 15
                  aj = 0.d0
                  mt = 0.d0
                  lum = 1.0d-10
                  r = 1.0d-10
               else
                  if(nsflag.eq.0)then
                     mt = 1.17d0 + 0.09d0*mc
                  elseif(nsflag.ge.1)then
*
* Use NS/BH mass given by Belczynski et al. 2002, ApJ, 572, 407. 
*
                     if(mc.lt.2.5d0)then
                        mcx = 0.161767d0*mc + 1.067055d0
                     else
                        mcx = 0.314154d0*mc + 0.686088d0
                     endif
                     if(mc.le.5.d0)then
                        mt = mcx
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                     endif
                  endif
                  mc = mt
                  if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                     kw = 13
                  else
*
* Zero-age Black hole
*
                     kw = 14
                  endif  
               endif
            endif
         endif
*
      endif
*
 90   continue
*
      if(kw.ge.7.and.kw.le.9)then
*
* Naked Helium Star
*
         rzams = rzhef(mt)
         rx = rzams
         if(aj.lt.tm)then
*
* Main Sequence
*
            kw = 7
            tau = aj/tm
            am = MAX(0.d0,0.85d0-0.08d0*mass)
            lum = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
            am = MAX(0.d0,0.4d0-0.22d0*LOG10(mt))
            r = rx*(1.d0+am*(tau-tau**6))
            rg = rx
* Star has no core mass and hence no memory of its past
* which is why we subject mass and mt to mass loss for
* this phase.
            mc = 0.d0
            if(mt.lt.zpars(10)) kw = 10
         else
*
* Helium Shell Burning
*
            kw = 8
            lum = lgbtf(aj,GB(8),GB,tscls(4),tscls(5),tscls(6))
            r = rhehgf(mt,lum,rx,lums(2))
            rg = rhegbf(lum)
            if(r.ge.rg)then
               kw = 9
               r = rg
            endif
            mc = mcgbf(lum,GB,lums(6))
            mtc = MIN(mt,1.45d0*mt-0.31d0)
            mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
            if(mc.ge.mcmax)then
               aj = 0.d0
               mc = mcmax
               if(mc.lt.mch)then
                  if(mass.lt.1.6d0)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                     mt = MAX(mc,(mc+0.31d0)/1.45d0)
                     kw = 11
                  else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                     mt = mc
                     kw = 12
                  endif
                  mass = mt
               else
                  if(mass.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                     kw = 15
                     aj = 0.d0
                     mt = 0.d0
                     lum = 1.0d-10
                     r = 1.0d-10
                  else
                     if(nsflag.eq.0)then
                        mt = 1.17d0 + 0.09d0*mc
                     elseif(nsflag.ge.1)then
                        if(mc.lt.2.5d0)then
                           mcx = 0.161767d0*mc + 1.067055d0
                        else
                           mcx = 0.314154d0*mc + 0.686088d0
                        endif
                        if(mc.le.5.d0)then
                           mt = mcx
                        elseif(mc.lt.7.6d0)then
                           mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                        endif
                     endif
                     mc = mt
                     if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                        kw = 13
                     else
*
* Zero-age Black hole
*
                        kw = 14
                     endif
                  endif  
               endif
            endif
         endif
      endif
*
      if(kw.ge.10.and.kw.le.12)then
*
*        White dwarf.
*
         mc = mt
         if(mc.ge.mch)then
*
* Accretion induced supernova with no remnant
* unless WD is ONe in which case we assume a NS 
* of minimum mass is the remnant. 
*
            if(kw.eq.12)then
               kw = 13
               mt = 1.3d0
            else
               kw = 15
               aj = 0.d0
               mt = 0.d0
               lum = 1.0d-10
               r = 1.0d-10
            endif
         else
*
            if(kw.eq.10)then
               xx = ahe
            else
               xx = aco
            endif
*
            if(wdflag.eq.0)then
*
* Mestel cooling
*
               lum = 635.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.4d0
*
            elseif(wdflag.ge.1)then
*
* modified-Mestel cooling
*
               if(aj.lt.9000.0)then
                  lum = 300.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.18d0
               else
                  fac = (9000.1d0*xx)**5.3d0
                  lum = 300.d0*fac*mt*zpars(14)/(xx*(aj+0.1d0))**6.48d0
               endif
*
            endif
*
            r = 0.0115d0*SQRT(MAX(1.48204d-06,(mch/mt)**(2.d0/3.d0)
     &                                      - (mt/mch)**(2.d0/3.d0)))
            r = MIN(0.1d0,r)
            if(mt.lt.0.0005d0) r = 0.09d0
            if(mt.lt.0.000005d0) r = 0.009d0
*
         endif
      endif
*
      if(kw.eq.13)then
*
*        Neutron Star.
*
         mc = mt
         if(mc.gt.mxns)then
*
* Accretion induced Black Hole?
*
            kw = 14
            aj = 0.d0
         else
            lum = 0.02d0*(mt**0.67d0)/(MAX(aj,0.1d0))**2
            r = 1.4d-05
         endif
      endif
*
      if(kw.eq.14)then
*
*        Black hole
*
         mc = mt
         lum = 1.0d-10
         r = 4.24d-06*mt
      endif
*
* Calculate the core radius and the luminosity and radius of the
* remnant that the star will become.
*
      tau = 0.d0
      if(kw.le.1.or.kw.eq.7)then
         rc = 0.d0
      elseif(kw.le.3)then
         if(mass.gt.zpars(2))then
            lx = lzhef(mc)
            rx = rzhef(mc)
            rc = rx
         else
            if(wdflag.eq.0)then
               lx = 635.d0*mc*zpars(14)/((ahe*0.1d0)**1.4d0)
            elseif(wdflag.ge.1)then
               lx = 300.d0*mc*zpars(14)/((ahe*0.1d0)**1.18d0)
            endif
            rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &           (mch/mc)**(2.d0/3.d0)-(mc/mch)**(2.d0/3.d0)))
            rc = 5.d0*rx
         endif
      elseif(kw.eq.4)then
         tau = (aj - tscls(2))/tscls(3)
         kwp = 7
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         am = MAX(0.d0,0.85d0-0.08d0*mc)
         lx = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
         rx = rzhef(mc)
         am = MAX(0.d0,0.4d0-0.22d0*LOG10(mc))
         rx = rx*(1.d0+am*(tau-tau**6))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.eq.5)then
         kwp = 9
         if(tn.gt.tbagb) tau = 3.d0*(aj-tbagb)/(tn-tbagb)
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         lx = lmcgbf(mcx,GB)
         if(tau.lt.1.d0) lx = lums(2)*(lx/lums(2))**tau
         rx = rzhef(mc)
         rx = MIN(rhehgf(mc,lx,rx,lums(2)),rhegbf(lx))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.le.9)then
         if(wdflag.eq.0)then
            lx = 635.d0*mc*zpars(14)/((aco*0.1d0)**1.4d0)
         elseif(wdflag.ge.1)then
            lx = 300.d0*mc*zpars(14)/((aco*0.1d0)**1.18d0)
         endif
         rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &        (mch/mc)**(2.d0/3.d0) - (mc/mch)**(2.d0/3.d0)))
         rc = 5.d0*rx
      else
         rc = r
         menv = 1.0d-10
         renv = 1.0d-10
         k2 = 0.21d0
      endif
*
* Perturb the luminosity and radius due to small envelope mass.
*
      if(kw.ge.2.and.kw.le.9.and.kw.ne.7)then
         mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
         if(kw.ge.8) mew = ((mtc-mc)/mtc)*5.d0
         if(mew.lt.1.d0)then
            xx = lpertf(mt,mew)
            lum = lx*(lum/lx)**xx
            if(r.le.rx)then
               xx = 0.d0
            else
               xx = rpertf(mt,mew,r,rx)
            endif
            r = rx*(r/rx)**xx
         endif
         rc = MIN(rc,r)
      endif
*
* Calculate mass and radius of convective envelope, and envelope
* gyration radius.
*
      if(kw.lt.10)then
         CALL mrenv(kw,mass,mt,mc,lum,r,rc,aj,tm,lums(2),lums(3),
     &              lums(4),rzams,rtms,rg,menv,renv,k2)
      endif
*
      if(mass.gt.99.99d0)then
         mass = mass0
      endif
      if(mt.gt.99.99d0)then
         mt = mt0
      endif
*
      return
      end
***
***
      SUBROUTINE instar
*
*
*       Initialization of collision matrix.
*       ------------------------
*
      implicit none
      integer i,j,ktype(0:14,0:14)
      common /TYPES/ ktype
*
*       Initialize stellar collision matrix.
*
      ktype(0,0) = 1
      do 10 , j = 1,6
         ktype(0,j) = j
         ktype(1,j) = j
 10   continue
      ktype(0,7) = 4
      ktype(1,7) = 4
      do 15 , j = 8,12
         if(j.ne.10)then
            ktype(0,j) = 6
         else
            ktype(0,j) = 3
         endif
         ktype(1,j) = ktype(0,j)
 15   continue
      ktype(2,2) = 3
      do 20 , i = 3,14
         ktype(i,i) = i
 20   continue
      ktype(5,5) = 4
      ktype(7,7) = 1
      ktype(10,10) = 15
      ktype(13,13) = 14
      do 25 , i = 2,5
         do 30 j = i+1,12
            ktype(i,j) = 4
 30      continue
 25   continue
      ktype(2,3) = 3
      ktype(2,6) = 5
      ktype(2,10) = 3
      ktype(2,11) = 5
      ktype(2,12) = 5
      ktype(3,6) = 5
      ktype(3,10) = 3
      ktype(3,11) = 5
      ktype(3,12) = 5
      ktype(6,7) = 4
      ktype(6,8) = 6
      ktype(6,9) = 6
      ktype(6,10) = 5 
      ktype(6,11) = 6
      ktype(6,12) = 6
      ktype(7,8) = 8
      ktype(7,9) = 9
      ktype(7,10) = 7
      ktype(7,11) = 9
      ktype(7,12) = 9
      ktype(8,9) = 9
      ktype(8,10) = 7
      ktype(8,11) = 9
      ktype(8,12) = 9
      ktype(9,10) = 7
      ktype(9,11) = 9
      ktype(9,12) = 9
      ktype(10,11) = 9
      ktype(10,12) = 9
      ktype(11,12) = 12
      do 35 , i = 0,12
         ktype(i,13) = 13
         ktype(i,14) = 14
 35   continue
      ktype(13,14) = 14
*
* Increase common-envelope cases by 100.
      do 40 , i = 0,9
         do 45 , j = i,14
            if(i.le.1.or.i.eq.7)then
               if(j.ge.2.and.j.le.9.and.j.ne.7)then
                  ktype(i,j) = ktype(i,j) + 100
               endif
            else
               ktype(i,j) = ktype(i,j) + 100
            endif
 45      continue
 40   continue
*
*       Assign the remaining values by symmetry.
      do 50 , i = 1,14
         do 55 , j = 0,i-1
            ktype(i,j) = ktype(j,i)
 55      continue
 50   continue
*
      return
      end
***
***
!! Subroutine written by Arkadiusz Hypki (ahypki@camk.edu.pl) - 13.03.2009
!!
!! Calculates velocities of binary components/star after SN explosion
!! This subroutine presents approach from:
!! ApJ 650:303-325, 10 Oct 2006 Belczynski et al.
!!
!! Variables description:
!!  - em = Eccentric anomaly
!!  - mm = Mean anomaly
!!  - sep = semi-major axis (Rsun)
!!  - phi = is the angle from the Z-axis, called also as colatitude or
!!    zenith, smaller than PI
!!  - theta = is the angle from the X-axis (as in the polar coordinates)
!!    smaller than 2*PI
!!  - vk = kick velocity (km/s)
!!
*
      SUBROUTINE kickv(kw,m1,m1n,m2,ecc,sep,jorb,oVs)
      implicit none
      CHARACTER ftos*(18);
      CHARACTER vtos*(60);
*
      integer idum
      COMMON /VALUE3/ idum
*
!! this COMMON is used to store velocity of the second star (companion
!! of the Supernova)
      Real*8 oVs2(3)
      COMMON /COMPANION/ oVs2 
*
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      external ran3
*
      integer kw,k,i,bhflag
      real*8 mm,em,dif,r,u1,u2,vk,v(4),s,theta,phi,sphi,cphi,stheta,
     &       ctheta,vr2,vk2,vn2,oVs(3),sigma,m1,m2,m1n,redMass,ecc, 
     &       sep,jorb,pi,piHalf,twopi,gmrkm,yearsc,rsunkm,EPS
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      COMMON /VALUE4/ sigma,bhflag
      real ran3,xx
*
      real*8 f1,f2,f3   ! temporal variables
      real*8 vr ! relative velocity (scalar value)
      real*8 vrVec(3) ! cartesian coordinates of the relative velocity
      real*8 vkVec(3) ! cartesian coordinates of the kick velocity
!! coordinate system I variables
      real*8 rVec(3)
      real*8 v1IVec(3), v2IVec(3)
      real*8 v1iIVec(3), v2iIVec(3) ! velocities after the kick
!! coordinate system II variables
      real*8 vCMIIVec(3) ! CM velocity
      real*8 vIIVec(3) ! relative velocity 
      real*8 nIIVec(3)
      real*8 JIIVec(3) ! angular momentum 
!! coordinate system III variables
      real*8 nIIIVec(3)
      real*8 pIII, eccIII, EIII, alphaIII, vfIII
      real*8 vfIIIVec(3), nfIIIVec(3)
      real*8 cosPhii, cosPhif, sinPhii, sinPhif, sf, cf
!! final velocitites in coordinate system I - velocities at infinity
      real*8 v1fIVec(3), v2fIVec(3)
      real*8 rotMat(3,3), rotMatInv(3,3)
*
      pi = ACOS(-1.d0)
      piHalf = pi*0.5d0
      twopi = 2.d0*pi
      EPS = 0.0001d0
      oVs2(1) = 0.d0
      oVs2(2) = 0.d0
      oVs2(3) = 0.d0
*
!! Conversion factor to ensure velocities are in km/s using mass and
!! radius in solar units.
      gmrkm = 1.906125d+05
!! M = E - e * sin(E)  <=>  mm = em - ecc * sin(em)
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
         mm = xx*twopi
         em = mm
*
         do 2
            dif = em - ecc*SIN(em) - mm
            if(ABS(dif/mm).le.1.0d-04) goto 3
            em = em - (dif/(1.d0 - ecc*COS(em)))
 2       continue
*
 3       continue
*
         r = sep*(1.d0 - ecc*COS(em))
!! Find the initial relative velocity scalar and vector.
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2) 
      else
         vr = 0.d0
         vr2 = 0.d0
      endif
!! Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
!! Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      do 20 k = 1,2
         u1 = RAN3(idum)
         u2 = RAN3(idum)
!! Generate two velocities from polar coordinates S & THETA.
         s = sigma*SQRT(-2.d0*LOG(1.d0 - u1))
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
*
      vk2 = v(1)**2 + v(2)**2 + v(3)**2
      vk = SQRT(vk2)

      if(kw.eq.14.and.bhflag.eq.0)then
!!     in the case of black holes and bhflag == 0 there is no kick velocity
         call DEBUG("black hole - there is no kick velocity")
         vk2 = 0.d0
         vk = 0.d0
      endif
!! drawing phi and theta angles - random orientation of the kick velocity
      sphi = -1.d0 + 2.d0*u1 
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
!! convert kick velocity into cartesian coordinates
      vkVec(1) = vk * sphi * ctheta;
      vkVec(2) = vk * sphi * stheta;
      vkVec(3) = vk * cphi;
      call DEBUG("kick velocity vector: "//vtos(vkVec))
      call DEBUG("kick velocity (length): "//ftos(vk))
                                   
!! coordinate system I - coordinate system before the SN explosion.
!! Center is in the center of mass. X axis is pointed to the pericenter
!! coordinates of the relative position and relative velocity
      call carth_coord_of_rel_vel_pos(em,ecc,sep,m1,m2,rVec,vrVec)
      f1 = m2 / (m1 + m2)
      v1IVec(1) = f1 * vrVec(1)
      v1IVec(2) = f1 * vrVec(2)
      v1IVec(3) = f1 * vrVec(3)
      f1 = - m1 / (m1 + m2)
      v2IVec(1) = f1 * vrVec(1)
      v2IVec(2) = f1 * vrVec(2)
      v2IVec(3) = f1 * vrVec(3)
!! velocities after the SN explosion in the coordinate system I
      do i = 1, 3
        v1iIVec(i) = v1IVec(i) + vkVec(i)
        v2iIVec(i) = v2IVec(i) 
      enddo
!! velocity of the center of mass (CM) after the explosion in
!! coordinate system II - vCMIIVec
      f1 = 1.d0/(m1n + m2)
      do i = 1, 3
        vCMIIVec(i) = f1 * (m1n * v1iIVec(i) + m2 * v2iIVec(i)) !! (6)
      enddo
!! in the case of only one star (no companion) new velocity of the
!! center of mass (vCMIIVec) is what we are looking for
      if (sep.eq.-1.d0.and.ecc.eq.0.d0) then
          do i = 1, 3
            oVs(i) = vCMIIVec(i)
          enddo
          call DEBUG("Single star case");
          call DEBUG("Final velocity: "//vtos(vCMIIVec))
          call DEBUG("Final velocity length: "//ftos(sqrt(
     &    vCMIIVec(1)*vCMIIVec(1)+vCMIIVec(2)*vCMIIVec(2)+
     &    vCMIIVec(3)*vCMIIVec(3))))
          return
      endif
!! relative velocity in coord sys II - vIIVec
      do i = 1, 3 !! vIIVec = v1iIVec - v2iIVec
        vIIVec(i) = vrVec(i) + vkVec(i);
      enddo
!! compute reduced mass and then angular momentum
      redMass = (m1n * m2) / (m1n + m2)
      do i = 1, 3
        nIIVec(i) = rVec(i) / r
      enddo
      f1 = redMass * r
      JIIVec(1) = f1 * (nIIVec(2) * vIIVec(3) - nIIVec(3) * vIIVec(2));
      JIIVec(2) = f1 * (nIIVec(3) * vIIVec(1) - nIIVec(1) * vIIVec(3));
      JIIVec(3) = f1 * (nIIVec(1) * vIIVec(2) - nIIVec(2) * vIIVec(1));
!! Coordinate system III - angular momentum is aligned with z-axis -
!! tranformation from II to III is rotation by matrix R
      call buildRotMat(JIIVec, rotMat)
!! inverve for of the rotation matrix is the same rotation matrix
!! but transposed
      call mTranspose(rotMat, rotMatInv)
!! move vectors n to coord sys III
      do i = 1, 3
        nIIIVec(i)=rotMat(1, i)*nIIVec(1) + rotMat(2, i)*nIIVec(2) +
     &             rotMat(3, i)*nIIVec(3)
      enddo
!! set angular momentum of the binary stars
      jorb = sqrt(JIIVec(1)*JIIVec(1) + JIIVec(2)*JIIVec(2) +
     &       JIIVec(3)*JIIVec(3))
!! in the III coordinate system two stars move on a hyperbolic orbit,
!! binary star is disrupted
      alphaIII = gmrkm * m1n * m2
!! Determine the magnitude of the new relative velocity and orbit parameters.
      vn2 = vIIVec(1)*vIIVec(1) + vIIVec(2)*vIIVec(2) +
     &      vIIVec(3)*vIIVec(3)
      EIII = redMass * vn2 * 0.5d0 - alphaIII / r
      eccIII = sqrt(1.d0 + (2.0d0*EIII*jorb*jorb) /
     &         (alphaIII * alphaIII * redMass))
      pIII = (jorb*jorb) / (alphaIII * redMass)
!! Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
!! if Energy in III is positive then binary is disrupted
      if (EIII.gt.0.d0) then
          call DEBUG("Binary is disrupted")
          ecc = 1.1d0
!! final relative velocity in IIIr
          vfIII = sqrt(2.d0 * EIII / redMass) !! (11)
!! final relative velocity in III at infinity is parallel to the
!! separation vector nfIII. 
!! it remains to find the angle between nIIIVec and nfIIIVec
          cosPhif = -(1.d0 / eccIII)
          sinPhif = sqrt(1.0d0 - cosPhif * cosPhif)
          cosPhii = (pIII / r - 1.0d0) / eccIII !! (12)
          if (cosPhii > 0.999999d0) then
             cosPhii = 1.0d0
             sinPhii = 0.0d0
          else
             sinPhii = sqrt(1.0d0 - cosPhii * cosPhii)
          endif                                                              }
*
          f2 = vIIVec(1) * nIIVec(1) + vIIVec(2) * nIIVec(2) +
     &         vIIVec(3) * nIIVec(3)
          if(f2.lt.0.0d0) sinPhii = -sinPhii
*
          cf = cosPhif * cosPhii + sinPhif * sinPhii
          sf = sinPhif * cosPhii - cosPhif * sinPhii
!! velocity at infinity after - rotation of vector nIII to nfIII
          VfIIIVec(1) = vfIII * (nIIIVec(1) * cf - nIIIVec(2) * sf)
          VfIIIVec(2) = vfIII * (nIIIVec(1) * sf + nIIIVec(2) * cf)
          VfIIIVec(3) = 0.0;
!! finally one can calculate velocities of the individual stars
!! in I coordinate system
          do i = 1, 3
            f1 = 0.0d0
            f2 = 0.0d0
            do k = 1, 3
                f1 = f1 + rotMatInv(k, i) * (m2/(m1n+m2)) * vfIIIVec(k)
                f2 = f2 + rotMatInv(k, i)*(-m1n/(m1n + m2))*vfIIIVec(k)
            enddo
            v1fIVec(i) = vCMIIVec(i) + f1 
            v2fIVec(i) = vCMIIVec(i) + f2 
          enddo
!! printing final velocities of each star
          call DEBUG("final vel 1: "//vtos(v1fIVec))
          f3 = sqrt(v1fIVec(1)*v1fIVec(1) + v1fIVec(2)*v1fIVec(2) +
     &              v1fIVec(3)*v1fIVec(3))
          call DEBUG("final vel 1 length: "//ftos(f3))
          call DEBUG("final vel 2: "//vtos(v2fIVec))
          f3 = sqrt(v2fIVec(1)*v2fIVec(1) + v2fIVec(2)*v2fIVec(2) +
     &              v2fIVec(3)*v2fIVec(3))
          call DEBUG("final vel 2 length: "//ftos(f3))

          do i = 1, 3
            oVs(i) = v1fIVec(i)
            oVs2(i) = v2fIVec(i)
          enddo
      else
          call DEBUG("Binary survives")
! rewriting new eccentricity; new jorb and separation were calculated before
          ecc = eccIII
!! new binary velocity for the center of mass
          do i = 1, 3
            oVs(i) = vCMIIVec(i)
          enddo
          call DEBUG("Final CM vector: "//vtos(vCMIIVec))
          f1 = sqrt(vCMIIVec(1)*vCMIIVec(1) + vCMIIVec(2)*vCMIIVec(2) +
     &              vCMIIVec(3)*vCMIIVec(3))
          call DEBUG("Final CM vec (length): "//ftos(f1))
      endif

*
      RETURN
      END
*
*
!=============================================================================
!=============================================================================
        subroutine buildRotMat(JVec, rotMat)
        implicit none
        real*8 JVec(3),rotMat(3,3),EPS,length,xy,s1,s2,c1,c2
*
        EPS = 0.0000001d0
*
        xy = JVec(1) * JVec(1) + JVec(2) * JVec(2);
        length = sqrt(xy + JVec(3) * JVec(3));
        xy = sqrt(xy)
*
        if (xy.gt.EPS) then
            s1 = JVec(2) / xy;
            c1 = JVec(1) / xy;
        else
            s1 = 0.0;
            c1 = 1.0;
        endif
        if (length.gt.EPS) then
            s2 = xy / length;
            c2 = JVec(3) / length;
        else
            s2 = 0.0;
            c2 = 1.0;
        endif
!! Below rotation matrix (rotMat) is a result of a multiplication of two rotation
!! matrices: rotMat = Ry * Rz. Ry is clockwise rotation matrix around Y axis
!! and the second one is also clockwise rotation matrix but around Z axis.
!! Rotation around Z axis is first so matrix Rz is on the right side of Ry
        rotMat(1, 1) = c1 * c2;
        rotMat(2, 1) = s1 * c2;
        rotMat(3, 1) = -s2;
*
        rotMat(1, 2) = -s1;
        rotMat(2, 2) = c1;
        rotMat(3, 2) = 0.0;
*
        rotMat(1, 3) = c1 * s2;
        rotMat(2, 3) = s1 * s2;
        rotMat(3, 3) = c2;
*
        end
*
*
!=============================================================================
!=============================================================================
      subroutine carth_coord_of_rel_vel_pos(em,e,sep,m1,m2,rVec,vrVec)
      implicit none
      real*8 em,e,sep,m1,m2,rVec(3),vrVec(3),vc,G
*
!! in the case if kickv subroutine is called for one star
      if (sep.le.0.0d0) then
          rVec(1) = 0.0d0
          rVec(2) = 0.0d0
          rVec(3) = 0.0d0
          vrVec(1) = 0.0d0
          vrVec(2) = 0.0d0
          vrVec(3) = 0.0d0
	      return
      endif
*
      rVec(1) = sep * (cos(em) - e)
      rVec(2) = sep * sqrt(1 - e*e) * sin(em)
      rVec(3) = 0.d0
*                  
      G = 1.906125d+05 !! mass and radius in solar units, velocity km/s                 !! in km/s
      vc = sqrt(G * (m1 + m2) / sep)
*
      vrVec(1) = - vc * (sin(em)) / (1 - e * cos(em))
      vrVec(2) = vc * (sqrt(1-e*e) * cos(em)) / (1 - e*cos(em))
      vrVec(3) = 0.0d0
*
      end
*
*
!=============================================================================
!=============================================================================
        subroutine mTranspose(iMatrix, oMatrix)
        implicit none
        Real*8 iMatrix(3,3),oMatrix(3,3)
        Integer i, k

        do i = 1, 3
            do k = 1, 3
                oMatrix(i, k) = iMatrix(k, i)
            enddo
        enddo

        end
*
*
!=============================================================================
!=============================================================================
        subroutine debug(m1)
        implicit none
        CHARACTER m1*(*)

        write (*,*) "DEBUG: ", m1

        end
*
*
!=============================================================================
!=============================================================================
        FUNCTION ftos(f)
        implicit none
        Real*8 f
        CHARACTER ftos*(18);
!        WRITE(ftos,'(1pE16.8)') f !! scientific format
        WRITE(ftos,'(F16.4)') f !! 12 digits where 5 digits are after dot
        return;
        END
*
*
!=============================================================================
!=============================================================================
        FUNCTION vtos(iVec)
        implicit none
        Real*8 iVec(3)
        CHARACTER vtos*(60);
!        WRITE(vtos,'(A,E16.4,A,E16.4,A,E16.4,A)') "[",iVec(1),", ",
        WRITE(vtos,'(A,F16.4,A,F16.4,A,F16.4,A)') "[",iVec(1),", ",
     &             iVec(2),", ",iVec(3),"]"
        return;
        END
*
*
c      SUBROUTINE kickv(kw,m1,m1n,m2,ecc,sep,jorb,vs)
c      implicit none
*
c      integer kw,k
c      INTEGER idum
c      COMMON /VALUE3/ idum
c      INTEGER idum2,iy,ir(32)
c      COMMON /RAND3/ idum2,iy,ir
c      integer bhflag
c      real*8 m1,m2,m1n,ecc,sep,jorb,ecc2
c      real*8 pi,twopi,gmrkm,yearsc,rsunkm
c      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
c      real*8 mm,em,dif,der,del,r
c      real*8 u1,u2,vk,v(4),s,theta,phi
c      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
c     real*8 vr,vr2,vk2,vn2,hn2
c      real*8 mu,cmu,vs(3),v1,v2,mx1,mx2
c      real*8 sigma
c      COMMON /VALUE4/ sigma,bhflag
c      real ran3,xx
c      external ran3
*
c      do k = 1,3
c         vs(k) = 0.d0
c      enddo
*
c      pi = ACOS(-1.d0)
c      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
c      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
c      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
c         xx = RAN3(idum)
c         mm = xx*twopi
c         em = mm
c 2       dif = em - ecc*SIN(em) - mm
c         if(ABS(dif/mm).le.1.0d-04) goto 3
c         der = 1.d0 - ecc*COS(em)
c         del = dif/der
c         em = em - del
c         goto 2
c 3       continue
c         r = sep*(1.d0 - ecc*COS(em))
*
* Find the initial relative velocity vector.
c         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
c         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
c         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
c         vr = SQRT(vr2)
c      else
c         vr = 0.d0
c         vr2 = 0.d0
c         salpha = 0.d0
c         calpha = 0.d0
c      endif
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
c      do 20 k = 1,2
c         u1 = RAN3(idum)
c         u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
c         s = sigma*SQRT(-2.d0*LOG(1.d0 - u1))
c         theta = twopi*u2
c         v(2*k-1) = s*COS(theta)
c         v(2*k) = s*SIN(theta)
c 20   continue
c      vk2 = v(1)**2 + v(2)**2 + v(3)**2
c      vk = SQRT(vk2)
c      print*,'kickv: vk, sigma',vk,sigma
c      if(kw.eq.14.and.bhflag.eq.0)then
c         vk2 = 0.d0
c         vk = 0.d0
c      endif
c      sphi = -1.d0 + 2.d0*u1
c      phi = ASIN(sphi)
c      cphi = COS(phi)
c      stheta = SIN(theta)
c      ctheta = COS(theta)
c      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*
* Determine the magnitude of the new relative velocity.
c      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha-stheta*cphi*calpha)
* Calculate the new semi-major axis.
c      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
c      sep = 1.d0/sep
c      if(sep.le.0.d0)then
c         ecc = 1.1d0
c         goto 90
c      endif
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
c      v1 = vk2*sphi*sphi
c      v2 = (vk*ctheta*cphi-vr*salpha)**2
c      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
c      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
c      ecc2 = MAX(0.d0,ecc2)
c      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
c      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the angle between the new and old orbital angular
* momentum vectors.
c      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
c      mu = ACOS(cmu)
* Calculate the components of the velocity of the new centre-of-mass.
c 90   mx1 = vk*m1n/(m1n+m2)
c      mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
c      vs(1) = mx1*ctheta*cphi + mx2*salpha
c      vs(2) = mx1*stheta*cphi + mx2*calpha
c      vs(3) = mx1*sphi
*
c      RETURN
c      END
***
***
      SUBROUTINE MIX(M0,M,AJ,KS,ZPARS)
*
*     Author : J. R. Hurley
*     Date :   7th July 1998
*
*       Evolution parameters for mixed star.
*       ------------------------------------
*
      implicit none
*
      INTEGER KS(2),I1,I2,K1,K2,KW,ICASE
      INTEGER KTYPE(0:14,0:14)
      COMMON /TYPES/ KTYPE
      REAL*8 M0(2),M(2),AJ(2),ZPARS(20)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TMS1,TMS2,TMS3,TN
      REAL*8 M01,M02,M03,M1,M2,M3,AGE1,AGE2,AGE3,MC3,MCH
      PARAMETER(MCH=1.44D0)
      REAL*8 NETA,BWIND,MXNS
      COMMON /VALUE1/ NETA,BWIND,MXNS
*
*
*       Define global indices with body #I1 being most evolved.
      IF(KS(1).GE.KS(2))THEN
          I1 = 1
          I2 = 2
      ELSE
          I1 = 2
          I2 = 1
      END IF
*
*       Specify case index for collision treatment.
      K1 = KS(I1)
      K2 = KS(I2)
      ICASE = KTYPE(K1,K2)
*     if(icase.gt.100) WRITE(66,*)' MIX ERROR ICASE>100 ',icase,k1,k2
*
*       Determine evolution time scales for first star.
      M01 = M0(I1)
      M1 = M(I1)
      AGE1 = AJ(I1)
      CALL star(K1,M01,M1,TMS1,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain time scales for second star.
      M02 = M0(I2)
      M2 = M(I2)
      AGE2 = AJ(I2)
      CALL star(K2,M02,M2,TMS2,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Check for planetary systems - defined as HeWDs and low-mass WDs!
      IF(K1.EQ.10.AND.M1.LT.0.05)THEN
         ICASE = K2
         IF(K2.LE.1)THEN
            ICASE = 1
            AGE1 = 0.D0
         ENDIF
      ELSEIF(K1.GE.11.AND.M1.LT.0.5.AND.ICASE.EQ.6)THEN
         ICASE = 9
      ENDIF
      IF(K2.EQ.10.AND.M2.LT.0.05)THEN
         ICASE = K1
         IF(K1.LE.1)THEN
            ICASE = 1
            AGE2 = 0.D0
         ENDIF
      ENDIF
*
*       Specify total mass.
      M3 = M1 + M2
      M03 = M01 + M02
      KW = ICASE
      AGE3 = 0.d0
*
*       Restrict merged stars to masses less than 100 Msun. 
      IF(M3.GE.100.D0)THEN
         M3 = 99.D0
         M03 = MIN(M03,M3)
      ENDIF
*
*       Evaluate apparent age and other parameters.
*
      IF(ICASE.EQ.1)THEN
*       Specify new age based on complete mixing.
         IF(K1.EQ.7) KW = 7
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = 0.1d0*TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
      ELSEIF(ICASE.EQ.3.OR.ICASE.EQ.6.OR.ICASE.EQ.9)THEN
         MC3 = M1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.4)THEN
         MC3 = M1
         AGE3 = AGE1/TMS1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.7)THEN
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = TMS3*(AGE2*M2/TMS2)/M3
      ELSEIF(ICASE.LE.12)THEN
*       Ensure that a new WD has the initial mass set correctly.
         M03 = M3
         IF(ICASE.LT.12.AND.M3.GE.MCH)THEN
            M3 = 0.D0
            KW = 15
         ENDIF
      ELSEIF(ICASE.EQ.13.OR.ICASE.EQ.14)THEN
*       Set unstable Thorne-Zytkow object with fast mass loss of envelope 
*       unless the less evolved star is a WD, NS or BH. 
         IF(K2.LT.10)THEN
            M03 = M1
            M3 = M1
         ENDIF
         IF(ICASE.EQ.13.AND.M3.GT.MXNS) KW = 14
      ELSEIF(ICASE.EQ.15)THEN
         M3 = 0.D0
      ELSEIF(ICASE.GT.100)THEN
*       Common envelope case which should only be used after COMENV.
         KW = K1
         AGE3 = AGE1
         M3 = M1
         M03 = M01
      ELSE
*       This should not be reached.
        KW = 1
        M03 = M3
      ENDIF
*
* Put the result in *1.
*
      KS(1) = KW
      KS(2) = 15
      M(1) = M3
      M(2) = 0.D0
      M0(1) = M03
      AJ(1) = AGE3
*
      RETURN
      END
***
***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,rl,z
      real*8 dml,dms,dmt,p0,x,mew,lum0,kap,neta,bwind,mxns
      parameter(lum0=7.0d+04,kap=-0.5d0)
      common /value1/ neta,bwind,mxns
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
      dms = 0.d0
      if(lum.gt.4000.d0)then
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**(1.d0/2.d0)
      endif
      if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
         dml = neta*4.0d-13*r*lum/mt
         if(rl.gt.0.d0) dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641, 
* for high pulsation periods on AGB.
         if(kw.eq.5.or.kw.eq.6)then
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dml = MAX(dml,dmt)
         endif
         if(kw.gt.6)then
            dms = MAX(dml,1.0d-13*lum**(3.d0/2.d0))
         else
            dms = MAX(dml,dms)
            mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
            if(mew.lt.1.d0)then
               dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
               dms = MAX(dml,dms)
            end if
* LBV-like mass loss beyond the Humphreys-Davidson limit.
            x = 1.0d-5*r*sqrt(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0)then
               dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
               dms = dms + dml
            endif
         endif
      endif
*
      mlwind = dms
*
      return
      end
***
***
      SUBROUTINE mrenv(kw,mass,mt,mc,lum,rad,rc,aj,tm,ltms,lbgb,lhei,
     &                 rzams,rtms,rg,menv,renv,k2e)
      implicit none
      integer kw
      real*8 mass,mt,mc,lum,rad,rc,aj,tm
      real*8 k2e,menv,menvg,menvt,menvz,renv,renvg,renvt,renvz
      real*8 A,B,C,D,E,F,x,y
      real*8 k2bgb,k2g,k2z,logm,logmt,lbgb,ltms,lhei,rg,rtms,rzams
      real*8 teff,tebgb,tetms,tau,tauenv,tautms
*
* A function to estimate the mass and radius of the convective envelope,
* as well as the gyration radius of the envelope.
* N.B. Valid only for Z=0.02!
*
* The following input is needed from HRDIAG:
*   kw = stellar type
*   mass = zero-age stellar mass
*   mt = actual mass
*   mc = core mass (not really needed, can also be done outside subroutine)
*   lum = luminosity
*   rad = radius
*   rc = core radius (not really needed...)
*   aj = age
*   tm = main-sequence lifetime
*   ltms = luminosity at TMS, lums(2)
*   lbgb = luminosity at BGB, lums(3)
*   lhei = luminosity at He ignition, lums(4)
*   rzams = radius at ZAMS
*   rtms = radius at TMS
*   rg = giant branch or Hayashi track radius, approporaite for the type. 
*        For kw=1 or 2 this is radius at BGB, and for kw=4 either GB or 
*        AGB radius at present luminosity.
*
      logm = log10(mass)
      A = MIN(0.81d0,MAX(0.68d0,0.68d0+0.4d0*logm))
      C = MAX(-2.5d0,MIN(-1.5d0,-2.5d0+5.d0*logm))
      D = -0.1d0
      E = 0.025d0
*
* Zero-age and BGB values of k^2.
*
      k2z = MIN(0.21d0,MAX(0.09d0-0.27d0*logm,0.037d0+0.033d0*logm))
      if(logm.gt.1.3d0) k2z = k2z - 0.055d0*(logm-1.3d0)**2
      k2bgb = MIN(0.15d0,MIN(0.147d0+0.03d0*logm,0.162d0-0.04d0*logm))
*
      if(kw.ge.3.and.kw.le.6)then
*
* Envelope k^2 for giant-like stars; this will be modified for non-giant
* CHeB stars or small envelope mass below.
* Formula is fairly accurate for both FGB and AGB stars if M <= 10, and
* gives reasonable values for higher masses. Mass dependence is on actual
* rather than ZA mass, expected to work for mass-losing stars (but not
* tested!). The slightly complex appearance is to insure continuity at 
* the BGB, which depends on the ZA mass.
*
         logmt = log10(mt)
         F = 0.208d0 + 0.125d0*logmt - 0.035d0*logmt**2
         B = 1.0d+04*mt**(3.d0/2.d0)/(1.d0+0.1d0*mt**(3.d0/2.d0))
         x = ((lum-lbgb)/B)**2
         y = (F - 0.033d0*log10(lbgb))/k2bgb - 1.d0
         k2g = (F - 0.033d0*log10(lum) + 0.4d0*x)/(1.d0+y*(lbgb/lum)+x)
      elseif(kw.eq.9)then
*
* Rough fit for for HeGB stars...
*
         B = 3.0d+04*mt**(3.d0/2.d0)
         x = (MAX(0.d0,lum/B-0.5d0))**2
         k2g = (k2bgb + 0.4d0*x)/(1.d0 + 0.4d0*x)
      else
         k2g = k2bgb
      endif
*
      if(kw.le.2)then
         menvg = 0.5d0
         renvg = 0.65d0
      elseif(kw.eq.3.and.lum.lt.3.d0*lbgb)then
*
* FGB stars still close to the BGB do not yet have a fully developed CE.
*
         x = MIN(3.d0,lhei/lbgb)
         tau = MAX(0.d0,MIN(1.d0,(x-lum/lbgb)/(x-1.d0)))
         menvg = 1.d0 - 0.5d0*tau**2
         renvg = 1.d0 - 0.35d0*tau**2
      else
         menvg = 1.d0
         renvg = 1.d0
      endif
*
      if(rad.lt.rg)then
*
* Stars not on the Hayashi track: MS and HG stars, non-giant CHeB stars,
* HeMS and HeHG stars, as well as giants with very small envelope mass.
*
         
         if(kw.le.6)then
*
* Envelope k^2 fitted for MS and HG stars.
* Again, pretty accurate for M <= 10 but less so for larger masses.
* [Note that this represents the whole star on the MS, so there is a 
* discontinuity in stellar k^2 between MS and HG - okay for stars with a 
* MS hook but low-mass stars should preferably be continous...]
*
* For other types of star not on the Hayashi track we use the same fit as 
* for HG stars, this is not very accurate but has the correct qualitative 
* behaviour. For CheB stars this is an overestimate because they appear
* to have a more centrally concentrated envelope than HG stars.
*
            k2e = (k2z-E)*(rad/rzams)**C + E*(rad/rzams)**D
         elseif(kw.eq.7)then
* Rough fit for naked He MS stars.
            tau = aj/tm
            k2e = 0.08d0 - 0.03d0*tau
         elseif(kw.le.9)then
* Rough fit for HeHG stars.
            k2e = 0.08d0*rzams/rad
         endif
*
* tauenv measures proximity to the Hayashi track in terms of Teff.
* If tauenv>0 then an appreciable convective envelope is present, and
* k^2 needs to be modified.
*
         if(kw.le.2)then
            teff = sqrt(sqrt(lum)/rad)
            tebgb = sqrt(sqrt(lbgb)/rg)
            tauenv = MAX(0.d0,MIN(1.d0,(tebgb/teff-A)/(1.d0-A)))
         else
            tauenv = MAX(0.d0,MIN(1.d0,(sqrt(rad/rg)-A)/(1.d0-A)))
         endif
*
         if(tauenv.gt.0.d0)then
            menv = menvg*tauenv**5
            renv = renvg*tauenv**(5.d0/4.d0)
            if(kw.le.1)then
* Zero-age values for CE mass and radius.
               x = MAX(0.d0,MIN(1.d0,(0.1d0-logm)/0.55d0))
               menvz = 0.18d0*x + 0.82d0*x**5
               renvz = 0.4d0*x**(1.d0/4.d0) + 0.6d0*x**10
               y = 2.d0 + 8.d0*x
* Values for CE mass and radius at start of the HG.
               tetms = sqrt(sqrt(ltms)/rtms)
               tautms = MAX(0.d0,MIN(1.d0,(tebgb/tetms-A)/(1.d0-A)))
               menvt = menvg*tautms**5
               renvt = renvg*tautms**(5.d0/4.d0)
* Modified expressions during MS evolution.
               tau = aj/tm
               if(tautms.gt.0.d0)then
                  menv = menvz + tau**y*menv*(menvt - menvz)/menvt
                  renv = renvz + tau**y*renv*(renvt - renvz)/renvt
               else
                  menv = 0.d0
                  renv = 0.d0
               endif
               k2e = k2e + tau**y*tauenv**3*(k2g - k2e)
            else
               k2e = k2e + tauenv**3*(k2g - k2e)
            endif
         else
            menv = 0.d0
            renv = 0.d0
         endif
      else
*
* All other stars should be true giants.
*
         menv = menvg
         renv = renvg
         k2e = k2g
      endif
*
      menv = menv*(mt - mc)
      renv = renv*(rad - rc)
      menv = MAX(menv,1.0d-10)
      renv = MAX(renv,1.0d-10)
*
      return
      end
***
***
      REAL FUNCTION ran3(IDUM)
*
* Random number generator from Numerical Recipes, Press et al. pg 272.
*
      IMPLICIT NONE
      INTEGER j,k,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      PARAMETER(im1=2147483563,im2=2147483399,ia1=40014,ia2=40692)
      PARAMETER(iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32)
      INTEGER idum
      INTEGER idum2,iy,ir(ntab)
      COMMON /RAND3/ idum2,iy,ir
      DATA idum2/123456789/, iy/0/, ir/ntab*0/
      REAL am
*
      am = 1.0/float(im1)
      imm1 = im1 - 1
      ndiv = 1 + imm1/ntab
*
      if(idum.le.0)then
         idum = MAX(-idum,1)
         idum2 = idum
         do 11 , j = ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum = idum + im1
            if(j.le.ntab) ir(j) = idum
 11      continue
         iy = ir(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum = idum + im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2
      j = 1 + iy/ndiv
      iy = ir(j) - idum2
      ir(j) = idum
      if(iy.lt.1) iy = iy + imm1
      ran3 = am*iy
*
      RETURN
      END
***
***
      REAL*8 FUNCTION RL(Q)
      IMPLICIT NONE
      REAL*8 Q,P
*
* A function to evaluate R_L/a(q), Eggleton 1983.
*
      P = Q**(1.d0/3.d0)
      RL = 0.49d0*P*P/(0.6d0*P*P + LOG(1.d0+P))
*
      RETURN
      END
***
***
      SUBROUTINE star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
*
*       Stellar luminosity & evolution time. 
*       ------------------------------------
*
      implicit none
*
      integer kw
*
      real*8 mass,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 tgb,tbagb,mch,mcmax,mc1,mc2,mcbagb,dx,am
      real*8 lambda,tau,mtc,mass0
      parameter(mch=1.44d0)
*
      real*8 lzamsf,lzahbf,lzhef
      real*8 tbgbf,thookf,tHef,themsf,mcgbf,mcagbf,mcheif,mcgbtf
      real*8 ltmsf,lbgbf,lHeIf,lHef,lbagbf,lmcgbf
      external lzamsf,lzahbf,lzhef
      external tbgbf,thookf,tHef,themsf,mcgbf,mcagbf,mcheif,mcgbtf
      external ltmsf,lbgbf,lHeIf,lHef,lbagbf,lmcgbf
*
*       Computes the characteristic luminosities at different stages (LUMS),
*       and various timescales (TSCLS).
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout
*       and 24th October 1995 to include metallicity
*       and 13th December 1996 to include naked helium stars
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include Z=0.001 as well as Z=0.02, convective overshooting,
*       MS hook and more elaborate CHeB. It now also sets the Giant
*       Branch parameters relevant to the mass of the star.
*
*       ------------------------------------------------------------
*       Times: 1; BGB              2; He ignition   3; He burning
*              4; Giant t(inf1)    5; Giant t(inf2) 6; Giant t(Mx)
*              7; FAGB t(inf1)     8; FAGB t(inf2)  9; FAGB  t(Mx)
*             10; SAGB t(inf1)    11; SAGB t(inf2) 12; SAGB  t(Mx)
*             13; TP              14; t(Mcmax)     
*
*       LUMS:  1; ZAMS             2; End MS        3; BGB
*              4; He ignition      5; He burning    6; L(Mx)
*              7; BAGB             8; TP
*
*       GB:    1; effective A(H)   2; A(H,He)       3; B
*              4; D                5; p             6; q
*              7; Mx               8; A(He)         9; Mc,BGB
*
*       ------------------------------------------------------------
*
*
      mass0 = mass
      if(mass0.gt.100.d0) mass = 100.d0
*
      if(kw.ge.7.and.kw.le.9) goto 90
      if(kw.ge.10) goto 95
*
* MS and BGB times
*
      tscls(1) = tbgbf(mass)
      tm = MAX(zpars(8),thookf(mass))*tscls(1)
*
* Zero- and terminal age main sequence luminosity
*
      lums(1) = lzamsf(mass)
      lums(2) = ltmsf(mass)
*
* Set the GB parameters
*
      GB(1) = MAX(-4.8d0,MIN(-5.7d0+0.8d0*mass,-4.1d0+0.14d0*mass))
      GB(1) = 10.d0**GB(1)
      GB(2) = 1.27d-05
      GB(8) = 8.0d-05
      GB(3) = MAX(3.0d+04,500.d0 + 1.75d+04*mass**0.6d0)
      if(mass.le.2.0)then
         GB(4) = zpars(6)
         GB(5) = 6.d0
         GB(6) = 3.d0
      elseif(mass.lt.2.5)then
         dx = zpars(6) - (0.975d0*zpars(6) - 0.18d0*2.5d0)
         GB(4) = zpars(6) - dx*(mass - 2.d0)/(0.5d0)
         GB(5) = 6.d0 - (mass - 2.d0)/(0.5d0)
         GB(6) = 3.d0 - (mass - 2.d0)/(0.5d0)
      else
         GB(4) = MAX(-1.d0,0.5d0*zpars(6) - 0.06d0*mass)
         GB(4) = MAX(GB(4),0.975d0*zpars(6) - 0.18d0*mass)
         GB(5) = 5.d0
         GB(6) = 2.d0
      endif
      GB(4) = 10.d0**GB(4)
      GB(7) = (GB(3)/GB(4))**(1.d0/(GB(5)-GB(6)))
* Change in slope of giant L-Mc relation.
      lums(6) = GB(4)*GB(7)**GB(5)
*
      if(mass.lt.0.1d0.and.kw.le.1) goto 96
*
* HeI ignition luminosity
*
      lums(4) = lHeIf(mass,zpars(2)) 
      lums(7) = lbagbf(mass,zpars(2))
      if(mass.le.zpars(3))then
* Base of the giant branch luminosity
         lums(3) = lbgbf(mass) 
* Set GB timescales 
         tscls(4) = tscls(1) + (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &              ((GB(4)/lums(3))**((GB(5)-1.d0)/GB(5)))
         tscls(6) = tscls(4) - (tscls(4) - tscls(1))*((lums(3)/lums(6))
     &              **((GB(5)-1.d0)/GB(5)))
         tscls(5) = tscls(6) + (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &              ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
* Set Helium ignition time
         if(lums(4).le.lums(6))then
            tscls(2) = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &                      ((GB(4)/lums(4))**((GB(5)-1.d0)/GB(5)))
         else
            tscls(2) = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &                      ((GB(3)/lums(4))**((GB(6)-1.d0)/GB(6)))
         endif
         tgb = tscls(2) - tscls(1)
         if(mass.le.zpars(2))then
            mc1 = mcgbf(lums(4),GB,lums(6))
            mc2 = mcagbf(mass)
            lums(5) = lzahbf(mass,mc1,zpars(2))
            tscls(3) = tHef(mass,mc1,zpars(2))
         else
            lums(5) = lHef(mass)*lums(4)
            tscls(3) = tHef(mass,1.d0,zpars(2))*tscls(1)
         endif
      else
* Note that for M>zpars(3) there is no GB as the star goes from
* HG -> CHeB -> AGB. So in effect tscls(1) refers to the time of
* Helium ignition and not the BGB.
         tscls(2) = tscls(1)
         tscls(3) = tHef(mass,1.d0,zpars(2))*tscls(1)
* This now represents the luminosity at the end of CHeB, ie. BAGB
         lums(5) = lums(7)
* We set lums(3) to be the luminosity at the end of the HG
         lums(3) = lums(4)
      endif
*
* Set the core mass at the BGB.
*
      if(mass.le.zpars(2))then
         GB(9) = mcgbf(lums(3),GB,lums(6))
      elseif(mass.le.zpars(3))then
         GB(9) = mcheif(mass,zpars(2),zpars(9))
      else
         GB(9) = mcheif(mass,zpars(2),zpars(10))
      endif
*
* FAGB time parameters
*
      tbagb = tscls(2) + tscls(3)
      tscls(7) = tbagb + (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &               ((GB(4)/lums(7))**((GB(5)-1.d0)/GB(5)))
      tscls(9) = tscls(7) - (tscls(7) - tbagb)*((lums(7)/lums(6))
     &                                    **((GB(5)-1.d0)/GB(5)))
      tscls(8) = tscls(9) + (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &               ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
*
* Now to find Ltp and ttp using Mc,He,tp
*
      mcbagb = mcagbf(mass)
      mc1 = mcbagb
      if(mc1.ge.0.8d0.and.mc1.lt.2.25d0)then
* The star undergoes dredge-up at Ltp causing a decrease in Mc,He
         mc1 = 0.44d0*mc1 + 0.448d0
      endif
      lums(8) = lmcgbf(mc1,GB)
      if(mc1.le.GB(7))then
         tscls(13) = tscls(7) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                   (mc1**(1.d0-GB(5)))
      else
         tscls(13) = tscls(8) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                   (mc1**(1.d0-GB(6)))
      endif
*
* SAGB time parameters
*
      if(mc1.le.GB(7))then
         tscls(10) = tscls(13) + (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &               ((GB(4)/lums(8))**((GB(5)-1.d0)/GB(5)))
         tscls(12) = tscls(10) - (tscls(10) - tscls(13))*
     &               ((lums(8)/lums(6))**((GB(5)-1.d0)/GB(5)))
         tscls(11) = tscls(12) + (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &               ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
      else
         tscls(10) = tscls(7)
         tscls(12) = tscls(9)
         tscls(11) = tscls(13) + (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &               ((GB(3)/lums(8))**((GB(6)-1.d0)/GB(6)))
      endif
*
* Get an idea of when Mc,C = Mc,C,max on the AGB
      tau = tscls(2) + tscls(3)
      mc2 = mcgbtf(tau,GB(8),GB,tscls(7),tscls(8),tscls(9))
      mcmax = MAX(MAX(mch,0.773d0*mcbagb - 0.35d0),1.05d0*mc2)
*
      if(mcmax.le.mc1)then
         if(mcmax.le.GB(7))then
            tscls(14) = tscls(7) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                      (mcmax**(1.d0-GB(5)))
         else
            tscls(14) = tscls(8) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                      (mcmax**(1.d0-GB(6)))
         endif
      else
* Star is on SAGB and we need to increase mcmax if any 3rd
* dredge-up has occurred.
         lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
         mcmax = (mcmax - lambda*mc1)/(1.d0 - lambda)
         if(mcmax.le.GB(7))then
            tscls(14) = tscls(10) - (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &                      (mcmax**(1.d0-GB(5)))
         else
            tscls(14) = tscls(11) - (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &                      (mcmax**(1.d0-GB(6)))
         endif
      endif
      tscls(14) = MAX(tbagb,tscls(14))
      if(mass.ge.100.d0)then
         tn = tscls(2)
         goto 100
      endif
*
* Calculate the nuclear timescale - the time of exhausting
* nuclear fuel without further mass loss.
* This means we want to find when Mc = Mt which defines Tn and will
* be used in determining the timestep required. Note that after some 
* stars reach Mc = Mt there will be a Naked Helium Star lifetime
* which is also a nuclear burning period but is not included in Tn.
*
      if(ABS(mt-mcbagb).lt.1.0d-14.and.kw.lt.5)then
         tn = tbagb
      else
* Note that the only occurence of Mc being double-valued is for stars
* that have a dredge-up. If Mt = Mc where Mc could be the value taken
* from CHeB or from the AGB we need to check the current stellar type.
         if(mt.gt.mcbagb.or.(mt.ge.mc1.and.kw.gt.4))then
            if(kw.eq.6)then
               lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
               mc1 = (mt - lambda*mc1)/(1.d0 - lambda)
            else
               mc1 = mt
            endif
            if(mc1.le.GB(7))then
               tn = tscls(10) - (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &                         (mc1**(1.d0-GB(5)))
            else
               tn = tscls(11) - (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &                         (mc1**(1.d0-GB(6)))
            endif
         else
            if(mass.gt.zpars(3))then
               mc1 = mcheif(mass,zpars(2),zpars(10))
               if(mt.le.mc1)then
                  tn = tscls(2)
               else
                  tn = tscls(2) + tscls(3)*((mt - mc1)/(mcbagb - mc1))
               endif
            elseif(mass.le.zpars(2))then
               mc1 = mcgbf(lums(3),GB,lums(6))
               mc2 = mcgbf(lums(4),GB,lums(6))
               if(mt.le.mc1)then
                  tn = tscls(1)
               elseif(mt.le.mc2)then
                  if(mt.le.GB(7))then
                     tn = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &                               (mt**(1.d0-GB(5)))
                  else
                     tn = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &                               (mt**(1.d0-GB(6)))
                  endif
               else
                  tn = tscls(2) + tscls(3)*((mt - mc2)/(mcbagb - mc2))
               endif
            else
               mc1 = mcheif(mass,zpars(2),zpars(9))
               mc2 = mcheif(mass,zpars(2),zpars(10))
               if(mt.le.mc1)then
                  tn = tscls(1)
               elseif(mt.le.mc2)then
                  tn = tscls(1) + tgb*((mt - mc1)/(mc2 - mc1))
               else
                  tn = tscls(2) + tscls(3)*((mt - mc2)/(mcbagb - mc2))
               endif
            endif
         endif
      endif
      tn = MIN(tn,tscls(14))
*
      goto 100
*
 90   continue
*
* Calculate Helium star Main Sequence lifetime.
*
      tm = themsf(mass)
      tscls(1) = tm
*
* Zero- and terminal age Helium star main sequence luminosity
*
      lums(1) = lzhef(mass)
      am = MAX(0.d0,0.85d0-0.08d0*mass)
      lums(2) = lums(1)*(1.d0+0.45d0+am)
*
* Set the Helium star GB parameters
*
      GB(8) = 8.0d-05
      GB(3) = 4.1d+04
      GB(4) = 5.5d+04/(1.d0+0.4d0*mass**4)
      GB(5) = 5.d0
      GB(6) = 3.d0
      GB(7) = (GB(3)/GB(4))**(1.d0/(GB(5)-GB(6)))
* Change in slope of giant L-Mc relation.
      lums(6) = GB(4)*GB(7)**GB(5)
*
*** Set Helium star GB timescales 
*
      mc1 = mcgbf(lums(2),GB,lums(6))
      tscls(4) = tm + (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                                          mc1**(1.d0-GB(5))
      tscls(6) = tscls(4) - (tscls(4) - tm)*((GB(7)/mc1)
     &                                     **(1.d0-GB(5)))
      tscls(5) = tscls(6) + (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                                       GB(7)**(1.d0-GB(6))
*
* Get an idea of when Mc = MIN(Mt,Mc,C,max) on the GB
      mtc = MIN(mt,1.45d0*mt-0.31d0)
      if(mtc.le.0.d0) mtc = mt
      mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
      if(mcmax.le.GB(7))then
         tscls(14) = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                   (mcmax**(1.d0-GB(5)))
      else
         tscls(14) = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                   (mcmax**(1.d0-GB(6)))
      endif
      tscls(14) = MAX(tscls(14),tm)
      tn = tscls(14)
*
      goto 100
*
 95   continue
      tm = 1.0d+10
      tscls(1) = tm
 96   continue
      tn = 1.0d+10
*
 100  continue
      mass = mass0
*
      return
      end
***
***
      SUBROUTINE zcnsts(z,zpars)
* 
      implicit none
      integer kw
*
      real*8 z,zpars(20)
      real*8 tm,tn,tscls(20),lums(10),GB(10)
      real*8 lzs,dlzs,lz,lzd,dum1,m1,m2,rr,rb,mhefl,lhefl,thefl,lx
      real*8 tbgbf,thef,lbagbf,lheif,lhef,lzahbf
      real*8 rgbf,ragbf,rminf,mcgbf
      external tbgbf,thef,lbagbf,lheif,lhef,lzahbf
      external rgbf,ragbf,rminf,mcgbf
*
      include 'zdata.h'
      real*8 msp(200),gbp(200),c(5)
      common /MSCFF/ msp
      common /GBCFF/ gbp
      data c /3.040581d-01, 8.049509d-02, 8.967485d-02,
     &        8.780198d-02, 2.219170d-02/
*
*       ------------------------------------------------------------
*
*      zpars:  1; M below which hook doesn't appear on MS, Mhook.
*              2; M above which He ignition occurs non-degenerately, Mhef.
*              3; M above which He ignition occurs on the HG, Mfgb. 
*              4; M below which C/O ignition doesn't occur, Mup.
*              5; M above which C ignites in the centre, Mec.
*              6; value of log D for M<= zpars(3)
*              7; value of x for Rgb propto M^(-x)
*              8; value of x for tMS = MAX(tHOOK,x*tBGB)
*              9; constant for McHeIf when computing Mc,BGB, mchefl.
*             10; constant for McHeIf when computing Mc,HeI, mchefl.
*             11; hydrogen abundance.
*             12; helium abundance.
*             13; constant x in rmin = rgb*x**y used by LM CHeB.
*             14; z**0.4 to be used for WD L formula.
*
*       ------------------------------------------------------------
*
      lzs = log10(z/0.02d0)
      dlzs = 1.d0/(z*log(10.d0))
      lz = log10(z)
      lzd = lzs + 1.d0
*
      zpars(1) = 1.0185d0 + lzs*(0.16015d0 + lzs*0.0892d0)
      zpars(2) = 1.995d0 + lzs*(0.25d0 + lzs*0.087d0)
      zpars(3) = 16.5d0*z**0.06d0/(1.d0 + (1.0d-04/z)**1.27d0)
      zpars(4) = MAX(6.11044d0 + 1.02167d0*lzs, 5.d0)
      zpars(5) = zpars(4) + 1.8d0
      zpars(6) = 5.37d0 + lzs*0.135d0
      zpars(7) = c(1) + lzs*(c(2) + lzs*(c(3) + lzs*(c(4) + lzs*c(5))))
      zpars(8) = MAX(0.95d0,MAX(0.95d0-(10.d0/3.d0)*(z-0.01d0),
     &           MIN(0.99d0,0.98d0-(100.d0/7.d0)*(z-0.001d0))))
***
* Lzams
      msp(1) = xz(1)+lzs*(xz(2)+lzs*(xz(3)+lzs*(xz(4)+lzs*xz(5))))
      msp(2) = xz(6)+lzs*(xz(7)+lzs*(xz(8)+lzs*(xz(9)+lzs*xz(10))))
      msp(3) = xz(11)+lzs*(xz(12)+lzs*(xz(13)+lzs*(xz(14)+lzs*xz(15))))
      msp(4) = xz(16)+lzs*(xz(17)+lzs*(xz(18)+lzs*(xz(19)+lzs*xz(20))))
      msp(5) = xz(21)+lzs*(xz(22)+lzs*(xz(23)+lzs*(xz(24)+lzs*xz(25))))
      msp(6) = xz(26)+lzs*(xz(27)+lzs*(xz(28)+lzs*(xz(29)+lzs*xz(30))))
      msp(7) = xz(31)+lzs*(xz(32)+lzs*(xz(33)+lzs*(xz(34)+lzs*xz(35))))
* Rzams
      msp(8) = xz(36)+lzs*(xz(37)+lzs*(xz(38)+lzs*(xz(39)+lzs*xz(40))))
      msp(9) = xz(41)+lzs*(xz(42)+lzs*(xz(43)+lzs*(xz(44)+lzs*xz(45))))
      msp(10) = xz(46)+lzs*(xz(47)+lzs*(xz(48)+lzs*(xz(49)+lzs*xz(50))))
      msp(11) = xz(51)+lzs*(xz(52)+lzs*(xz(53)+lzs*(xz(54)+lzs*xz(55))))
      msp(12) = xz(56)+lzs*(xz(57)+lzs*(xz(58)+lzs*(xz(59)+lzs*xz(60))))
      msp(13) = xz(61)
      msp(14) = xz(62)+lzs*(xz(63)+lzs*(xz(64)+lzs*(xz(65)+lzs*xz(66))))
      msp(15) = xz(67)+lzs*(xz(68)+lzs*(xz(69)+lzs*(xz(70)+lzs*xz(71))))
      msp(16) = xz(72)+lzs*(xz(73)+lzs*(xz(74)+lzs*(xz(75)+lzs*xz(76))))
* Tbgb 
      msp(17) = xt(1)+lzs*(xt(2)+lzs*(xt(3)+lzs*xt(4)))
      msp(18) = xt(5)+lzs*(xt(6)+lzs*(xt(7)+lzs*xt(8)))
      msp(19) = xt(9)+lzs*(xt(10)+lzs*(xt(11)+lzs*xt(12)))
      msp(20) = xt(13)+lzs*(xt(14)+lzs*(xt(15)+lzs*xt(16)))
      msp(21) = xt(17)
* dTbgb/dz
      msp(117) = dlzs*(xt(2)+lzs*(2.d0*xt(3)+3.d0*lzs*xt(4)))
      msp(118) = dlzs*(xt(6)+lzs*(2.d0*xt(7)+3.d0*lzs*xt(8)))
      msp(119) = dlzs*(xt(10)+lzs*(2.d0*xt(11)+3.d0*lzs*xt(12)))
      msp(120) = dlzs*(xt(14)+lzs*(2.d0*xt(15)+3.d0*lzs*xt(16)))
* Thook
      msp(22) = xt(18)+lzs*(xt(19)+lzs*(xt(20)+lzs*xt(21)))
      msp(23) = xt(22)
      msp(24) = xt(23)+lzs*(xt(24)+lzs*(xt(25)+lzs*xt(26)))
      msp(25) = xt(27)+lzs*(xt(28)+lzs*(xt(29)+lzs*xt(30)))
      msp(26) = xt(31)
* Ltms 
      msp(27) = xl(1)+lzs*(xl(2)+lzs*(xl(3)+lzs*(xl(4)+lzs*xl(5))))
      msp(28) = xl(6)+lzs*(xl(7)+lzs*(xl(8)+lzs*(xl(9)+lzs*xl(10))))
      msp(29) = xl(11)+lzs*(xl(12)+lzs*(xl(13)+lzs*xl(14)))
      msp(30) = xl(15)+lzs*(xl(16)+lzs*(xl(17)+lzs*(xl(18)+lzs*xl(19))))
      msp(27) = msp(27)*msp(30)
      msp(28) = msp(28)*msp(30)
      msp(31) = xl(20)+lzs*(xl(21)+lzs*(xl(22)+lzs*xl(23)))
      msp(32) = xl(24)+lzs*(xl(25)+lzs*(xl(26)+lzs*xl(27)))
* Lalpha
      m2 = 2.d0
      msp(33) = xl(28)+lzs*(xl(29)+lzs*(xl(30)+lzs*xl(31)))
      msp(34) = xl(32)+lzs*(xl(33)+lzs*(xl(34)+lzs*xl(35)))
      msp(35) = xl(36)+lzs*(xl(37)+lzs*(xl(38)+lzs*xl(39)))
      msp(36) = xl(40)+lzs*(xl(41)+lzs*(xl(42)+lzs*xl(43)))
      msp(37) = MAX(0.9d0,1.1064d0+lzs*(0.415d0+0.18d0*lzs))
      msp(38) = MAX(1.d0,1.19d0+lzs*(0.377d0+0.176d0*lzs))
      if(z.gt.0.01d0)then
         msp(37) = MIN(msp(37),1.d0)
         msp(38) = MIN(msp(38),1.1d0)
      endif
      msp(39) = MAX(0.145d0,0.0977d0-lzs*(0.231d0+0.0753d0*lzs))
      msp(40) = MIN(0.24d0+lzs*(0.18d0+0.595d0*lzs),0.306d0+0.053d0*lzs)
      msp(41) = MIN(0.33d0+lzs*(0.132d0+0.218d0*lzs),
     &              0.3625d0+0.062d0*lzs)
      msp(42) = (msp(33)+msp(34)*m2**msp(36))/
     &          (m2**0.4d0+msp(35)*m2**1.9d0)
* Lbeta
      msp(43) = xl(44)+lzs*(xl(45)+lzs*(xl(46)+lzs*(xl(47)+lzs*xl(48))))
      msp(44) = xl(49)+lzs*(xl(50)+lzs*(xl(51)+lzs*(xl(52)+lzs*xl(53))))
      msp(45) = xl(54)+lzs*(xl(55)+lzs*xl(56))
      msp(46) = MIN(1.4d0,1.5135d0+0.3769d0*lzs)
      msp(46) = MAX(0.6355d0-0.4192d0*lzs,MAX(1.25d0,msp(46)))
* Lhook
      msp(47) = xl(57)+lzs*(xl(58)+lzs*(xl(59)+lzs*xl(60)))
      msp(48) = xl(61)+lzs*(xl(62)+lzs*(xl(63)+lzs*xl(64)))
      msp(49) = xl(65)+lzs*(xl(66)+lzs*(xl(67)+lzs*xl(68)))
      msp(50) = xl(69)+lzs*(xl(70)+lzs*(xl(71)+lzs*xl(72)))
      msp(51) = MIN(1.4d0,1.5135d0+0.3769d0*lzs)
      msp(51) = MAX(0.6355d0-0.4192d0*lzs,MAX(1.25d0,msp(51)))
* Rtms
      msp(52) = xr(1)+lzs*(xr(2)+lzs*(xr(3)+lzs*(xr(4)+lzs*xr(5))))
      msp(53) = xr(6)+lzs*(xr(7)+lzs*(xr(8)+lzs*(xr(9)+lzs*xr(10))))
      msp(54) = xr(11)+lzs*(xr(12)+lzs*(xr(13)+lzs*(xr(14)+lzs*xr(15))))
      msp(55) = xr(16)+lzs*(xr(17)+lzs*(xr(18)+lzs*xr(19)))
      msp(56) = xr(20)+lzs*(xr(21)+lzs*(xr(22)+lzs*xr(23)))
      msp(52) = msp(52)*msp(54)
      msp(53) = msp(53)*msp(54)
      msp(57) = xr(24)
      msp(58) = xr(25)+lzs*(xr(26)+lzs*(xr(27)+lzs*xr(28)))
      msp(59) = xr(29)+lzs*(xr(30)+lzs*(xr(31)+lzs*xr(32)))
      msp(60) = xr(33)+lzs*(xr(34)+lzs*(xr(35)+lzs*xr(36)))
      msp(61) = xr(37)+lzs*(xr(38)+lzs*(xr(39)+lzs*xr(40)))
*
      msp(62) = MAX(0.097d0-0.1072d0*(lz+3.d0),MAX(0.097d0,MIN(0.1461d0,
     &              0.1461d0+0.1237d0*(lz+2.d0))))
      msp(62) = 10.d0**msp(62)
      m2 = msp(62) + 0.1d0
      msp(63) = (msp(52)+msp(53)*msp(62)**msp(55))/
     &          (msp(54)+msp(62)**msp(56))
      msp(64) = (msp(57)*m2**3+msp(58)*m2**msp(61)+
     &           msp(59)*m2**(msp(61)+1.5d0))/(msp(60)+m2**5)
* Ralpha
      msp(65) = xr(41)+lzs*(xr(42)+lzs*(xr(43)+lzs*xr(44)))
      msp(66) = xr(45)+lzs*(xr(46)+lzs*(xr(47)+lzs*xr(48)))
      msp(67) = xr(49)+lzs*(xr(50)+lzs*(xr(51)+lzs*xr(52)))
      msp(68) = xr(53)+lzs*(xr(54)+lzs*(xr(55)+lzs*xr(56)))
      msp(69) = xr(57)+lzs*(xr(58)+lzs*(xr(59)+lzs*(xr(60)+lzs*xr(61))))
      msp(70) = MAX(0.9d0,MIN(1.d0,1.116d0+0.166d0*lzs))
      msp(71) = MAX(1.477d0+0.296d0*lzs,MIN(1.6d0,-0.308d0-1.046d0*lzs))
      msp(71) = MAX(0.8d0,MIN(0.8d0-2.d0*lzs,msp(71)))
      msp(72) = xr(62)+lzs*(xr(63)+lzs*xr(64))
      msp(73) = MAX(0.065d0,0.0843d0-lzs*(0.0475d0+0.0352d0*lzs))
      msp(74) = 0.0736d0+lzs*(0.0749d0+0.04426d0*lzs)
      if(z.lt.0.004d0) msp(74) = MIN(0.055d0,msp(74))
      msp(75) = MAX(0.091d0,MIN(0.121d0,0.136d0+0.0352d0*lzs))
      msp(76) = (msp(65)*msp(71)**msp(67))/(msp(66) + msp(71)**msp(68))
      if(msp(70).gt.msp(71))then
         msp(70) = msp(71)
         msp(75) = msp(76)
      endif
* Rbeta
      msp(77) = xr(65)+lzs*(xr(66)+lzs*(xr(67)+lzs*xr(68)))
      msp(78) = xr(69)+lzs*(xr(70)+lzs*(xr(71)+lzs*xr(72)))
      msp(79) = xr(73)+lzs*(xr(74)+lzs*(xr(75)+lzs*xr(76)))
      msp(80) = xr(77)+lzs*(xr(78)+lzs*(xr(79)+lzs*xr(80)))
      msp(81) = xr(81)+lzs*(xr(82)+lzs*lzs*xr(83))
      if(z.gt.0.01d0) msp(81) = MAX(msp(81),0.95d0)
      msp(82) = MAX(1.4d0,MIN(1.6d0,1.6d0+lzs*(0.764d0+0.3322d0*lzs)))
* Rgamma
      msp(83) = MAX(xr(84)+lzs*(xr(85)+lzs*(xr(86)+lzs*xr(87))),
     &              xr(96)+lzs*(xr(97)+lzs*xr(98)))
      msp(84) = MIN(0.d0,xr(88)+lzs*(xr(89)+lzs*(xr(90)+lzs*xr(91))))
      msp(84) = MAX(msp(84),xr(99)+lzs*(xr(100)+lzs*xr(101)))
      msp(85) = xr(92)+lzs*(xr(93)+lzs*(xr(94)+lzs*xr(95)))
      msp(85) = MAX(0.d0,MIN(msp(85),7.454d0+9.046d0*lzs))
      msp(86) = MIN(xr(102)+lzs*xr(103),MAX(2.d0,-13.3d0-18.6d0*lzs))
      msp(87) = MIN(1.5d0,MAX(0.4d0,2.493d0+1.1475d0*lzs))
      msp(88) = MAX(1.d0,MIN(1.27d0,0.8109d0-0.6282d0*lzs))
      msp(88) = MAX(msp(88),0.6355d0-0.4192d0*lzs)
      msp(89) = MAX(5.855420d-02,-0.2711d0-lzs*(0.5756d0+0.0838d0*lzs))
* Rhook
      msp(90) = xr(104)+lzs*(xr(105)+lzs*(xr(106)+lzs*xr(107)))
      msp(91) = xr(108)+lzs*(xr(109)+lzs*(xr(110)+lzs*xr(111)))
      msp(92) = xr(112)+lzs*(xr(113)+lzs*(xr(114)+lzs*xr(115)))
      msp(93) = xr(116)+lzs*(xr(117)+lzs*(xr(118)+lzs*xr(119)))
      msp(94) = MIN(1.25d0,
     &          MAX(1.1d0,1.9848d0+lzs*(1.1386d0+0.3564d0*lzs)))
      msp(95) = 0.063d0 + lzs*(0.0481d0 + 0.00984d0*lzs)
      msp(96) = MIN(1.3d0,MAX(0.45d0,1.2d0+2.45d0*lzs))
* Lneta
      if(z.gt.0.0009d0)then
         msp(97) = 10.d0
      else
         msp(97) = 20.d0
      endif
* Lbgb 
      gbp(1) = xg(1)+lzs*(xg(2)+lzs*(xg(3)+lzs*xg(4)))
      gbp(2) = xg(5)+lzs*(xg(6)+lzs*(xg(7)+lzs*xg(8)))
      gbp(3) = xg(9)+lzs*(xg(10)+lzs*(xg(11)+lzs*xg(12)))
      gbp(4) = xg(13)+lzs*(xg(14)+lzs*(xg(15)+lzs*xg(16)))
      gbp(5) = xg(17)+lzs*(xg(18)+lzs*xg(19))
      gbp(6) = xg(20)+lzs*(xg(21)+lzs*xg(22))
      gbp(3) = gbp(3)**gbp(6)
      gbp(7) = xg(23)
      gbp(8) = xg(24)
* Lbagb
* set gbp(16) = 1.d0 until it is reset later with an initial
* call to Lbagbf using mass = zpars(2) and mhefl = 0.0
      gbp(9) = xg(25) + lzs*(xg(26) + lzs*xg(27))
      gbp(10) = xg(28) + lzs*(xg(29) + lzs*xg(30))
      gbp(11) = 15.d0
      gbp(12) = xg(31)+lzs*(xg(32)+lzs*(xg(33)+lzs*xg(34)))
      gbp(13) = xg(35)+lzs*(xg(36)+lzs*(xg(37)+lzs*xg(38)))
      gbp(14) = xg(39)+lzs*(xg(40)+lzs*(xg(41)+lzs*xg(42)))
      gbp(15) = xg(43)+lzs*xg(44)
      gbp(12) = gbp(12)**gbp(15)
      gbp(14) = gbp(14)**gbp(15)
      gbp(16) = 1.d0
* Rgb
      gbp(17) = -4.6739d0-0.9394d0*lz
      gbp(17) = 10.d0**gbp(17)
      gbp(17) = MAX(gbp(17),-0.04167d0+55.67d0*z)
      gbp(17) = MIN(gbp(17),0.4771d0-9329.21d0*z**2.94d0)
      gbp(18) = MIN(0.54d0,0.397d0+lzs*(0.28826d0+0.5293d0*lzs))
      gbp(19) = MAX(-0.1451d0,-2.2794d0-lz*(1.5175d0+0.254d0*lz))
      gbp(19) = 10.d0**gbp(19)
      if(z.gt.0.004d0)then
         gbp(19) = MAX(gbp(19),0.7307d0+14265.1d0*z**3.395d0)
      endif
      gbp(20) = xg(45)+lzs*(xg(46)+lzs*(xg(47)+lzs*(xg(48)+
     &          lzs*(xg(49)+lzs*xg(50)))))
      gbp(21) = xg(51)+lzs*(xg(52)+lzs*(xg(53)+lzs*(xg(54)+lzs*xg(55))))
      gbp(22) = xg(56)+lzs*(xg(57)+lzs*(xg(58)+lzs*(xg(59)+
     &          lzs*(xg(60)+lzs*xg(61)))))
      gbp(23) = xg(62)+lzs*(xg(63)+lzs*(xg(64)+lzs*(xg(65)+lzs*xg(66))))
* Ragb
      gbp(24) = MIN(0.99164d0-743.123d0*z**2.83d0,
     &              1.0422d0+lzs*(0.13156d0+0.045d0*lzs))
      gbp(25) = xg(67)+lzs*(xg(68)+lzs*(xg(69)+lzs*(xg(70)+
     &          lzs*(xg(71)+lzs*xg(72)))))
      gbp(26) = xg(73)+lzs*(xg(74)+lzs*(xg(75)+lzs*(xg(76)+lzs*xg(77))))
      gbp(27) = xg(78)+lzs*(xg(79)+lzs*(xg(80)+lzs*(xg(81)+
     &          lzs*(xg(82)+lzs*xg(83)))))
      gbp(28) = xg(84)+lzs*(xg(85)+lzs*(xg(86)+lzs*(xg(87)+lzs*xg(88))))
      gbp(29) = xg(89)+lzs*(xg(90)+lzs*(xg(91)+lzs*(xg(92)+
     &          lzs*(xg(93)+lzs*xg(94)))))
      gbp(30) = xg(95)+lzs*(xg(96)+lzs*(xg(97)+lzs*(xg(98)+
     &          lzs*(xg(99)+lzs*xg(100)))))
      m1 = zpars(2) - 0.2d0
      gbp(31) = gbp(29) + gbp(30)*m1
      gbp(32) = MIN(gbp(25)/zpars(2)**gbp(26),gbp(27)/zpars(2)**gbp(28))
* Mchei
      gbp(33) = xg(101)**4
      gbp(34) = xg(102)*4.d0
* Mcagb
      gbp(35) = xg(103)+lzs*(xg(104)+lzs*(xg(105)+lzs*xg(106)))
      gbp(36) = xg(107)+lzs*(xg(108)+lzs*(xg(109)+lzs*xg(110)))
      gbp(37) = xg(111)+lzs*xg(112)
      gbp(35) = gbp(35)**4
      gbp(36) = gbp(36)*4.d0
      gbp(37) = gbp(37)**4
* Lhei
* set gbp(41) = -1.d0 until it is reset later with an initial
* call to Lheif using mass = zpars(2) and mhefl = 0.0
      gbp(38) = xh(1)+lzs*xh(2)
      gbp(39) = xh(3)+lzs*xh(4)
      gbp(40) = xh(5)
      gbp(41) = -1.d0
      gbp(42) = xh(6)+lzs*(xh(7)+lzs*xh(8))
      gbp(43) = xh(9)+lzs*(xh(10)+lzs*xh(11))
      gbp(44) = xh(12)+lzs*(xh(13)+lzs*xh(14))
      gbp(42) = gbp(42)**2
      gbp(44) = gbp(44)**2
* Lhe
      gbp(45) = xh(15)+lzs*(xh(16)+lzs*xh(17))
      if(lzs.gt.-1.d0)then
         gbp(46) = 1.d0 - xh(19)*(lzs+1.d0)**xh(18)
      else
         gbp(46) = 1.d0
      endif
      gbp(47) = xh(20)+lzs*(xh(21)+lzs*xh(22))
      gbp(48) = xh(23)+lzs*(xh(24)+lzs*xh(25))
      gbp(45) = gbp(45)**gbp(48)
      gbp(47) = gbp(47)**gbp(48)
      gbp(46) = gbp(46)/zpars(3)**0.1d0+(gbp(46)*gbp(47)-gbp(45))/
     &          zpars(3)**(gbp(48)+0.1d0)
* Rmin
      gbp(49) = xh(26)+lzs*(xh(27)+lzs*(xh(28)+lzs*xh(29)))
      gbp(50) = xh(30)+lzs*(xh(31)+lzs*(xh(32)+lzs*xh(33)))
      gbp(51) = xh(34)+lzs*(xh(35)+lzs*(xh(36)+lzs*xh(37)))
      gbp(52) = 5.d0+xh(38)*z**xh(39)
      gbp(53) = xh(40)+lzs*(xh(41)+lzs*(xh(42)+lzs*xh(43)))
      gbp(49) = gbp(49)**gbp(53)
      gbp(51) = gbp(51)**(2.d0*gbp(53))
* The
* set gbp(57) = -1.d0 until it is reset later with an initial
* call to Thef using mass = zpars(2), mc = 0.0  and mhefl = 0.0
      gbp(54) = xh(44)+lzs*(xh(45)+lzs*(xh(46)+lzs*xh(47)))
      gbp(55) = xh(48)+lzs*(xh(49)+lzs*xh(50))
      gbp(55) = MAX(gbp(55),1.d0)
      gbp(56) = xh(51)
      gbp(57) = -1.d0
      gbp(58) = xh(52)+lzs*(xh(53)+lzs*(xh(54)+lzs*xh(55)))
      gbp(59) = xh(56)+lzs*(xh(57)+lzs*(xh(58)+lzs*xh(59)))
      gbp(60) = xh(60)+lzs*(xh(61)+lzs*(xh(62)+lzs*xh(63)))
      gbp(61) = xh(64)+lzs*xh(65)
      gbp(58) = gbp(58)**gbp(61)
      gbp(60) = gbp(60)**5
* Tbl
      dum1 = zpars(2)/zpars(3)
      gbp(62) = xh(66)+lzs*xh(67)
      gbp(62) = -gbp(62)*log10(dum1)
      gbp(63) = xh(68)
      if(lzd.gt.0.d0) then
         gbp(64) = 1.d0-lzd*(xh(69)+lzd*(xh(70)+lzd*xh(71)))
      else
         gbp(64) = 1.d0
      end if
      gbp(65) = 1.d0-gbp(64)*dum1**gbp(63)
      gbp(66) = 1.d0 - lzd*(xh(77) + lzd*(xh(78) + lzd*xh(79)))
      gbp(67) = xh(72) + lzs*(xh(73) + lzs*(xh(74) + lzs*xh(75)))
      gbp(68) = xh(76)
* Lzahb
      gbp(69) = xh(80) + lzs*(xh(81) + lzs*xh(82))
      gbp(70) = xh(83) + lzs*(xh(84) + lzs*xh(85))
      gbp(71) = 15.d0
      gbp(72) = xh(86)
      gbp(73) = xh(87)
* Rzahb
      gbp(75) = xh(88) + lzs*(xh(89) + lzs*(xh(90) + lzs*xh(91)))
      gbp(76) = xh(92) + lzs*(xh(93) + lzs*(xh(94) + lzs*xh(95)))
      gbp(77) = xh(96) + lzs*(xh(97) + lzs*(xh(98) + lzs*xh(99)))
***
* finish Lbagb
      mhefl = 0.d0
      lx = lbagbf(zpars(2),mhefl)
      gbp(16) = lx
* finish LHeI
      dum1 = 0.d0
      lhefl = lheif(zpars(2),mhefl)
      gbp(41) = (gbp(38)*zpars(2)**gbp(39)-lhefl)/
     &          (EXP(zpars(2)*gbp(40))*lhefl)
* finish THe
      thefl = thef(zpars(2),dum1,mhefl)*tbgbf(zpars(2))
      gbp(57) = (thefl-gbp(54))/(gbp(54)*EXP(gbp(56)*zpars(2)))
* finish Tblf
      rb = ragbf(zpars(3),lheif(zpars(3),zpars(2)),mhefl)
      rr = 1.d0 - rminf(zpars(3))/rb
      rr = MAX(rr,1.0d-12)
      gbp(66) = gbp(66)/(zpars(3)**gbp(67)*rr**gbp(68))
* finish Lzahb
      gbp(74) = lhefl*lHef(zpars(2))
***
      kw = 0
      tm = 0.d0
      tn = 0.d0
      CALL star(kw,zpars(2),zpars(2),tm,tn,tscls,lums,GB,zpars)
      zpars(9) = mcgbf(lums(3),GB,lums(6))
      zpars(10) = mcgbf(lums(4),GB,lums(6))
* set the hydrogen and helium abundances
      zpars(11) = 0.76d0 - 3.d0*z
      zpars(12) = 0.24d0 + 2.d0*z
* set constant for low-mass CHeB stars
      zpars(13) = rminf(zpars(2))/
     &            rgbf(zpars(2),lzahbf(zpars(2),zpars(9),zpars(2)))
* 
      zpars(14) = z**0.4d0
*
      return
      end
***
***
      real*8 FUNCTION lzamsf(m)
      implicit none
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Lzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      lzamsf = (a(1)*m**5*mx + a(2)*m**11)/
     &         (a(3) + m**3 + a(4)*m**5 + a(5)*m**7 +
     &          a(6)*m**8 + a(7)*m**9*mx)
*
      return
      end
***
      real*8 FUNCTION rzamsf(m)
      implicit none
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Rzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      rzamsf = ((a(8)*m**2 + a(9)*m**6)*mx + a(10)*m**11 +
     &          (a(11) + a(12)*mx)*m**19)/
     &         (a(13) + a(14)*m**2 + 
     &          (a(15)*m**8 + m**18 + a(16)*m**19)*mx)
*
      return
      end
***
      real*8 FUNCTION tbgbf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the BGB or to
* Helium ignition if no FGB exists.
* (JH 24/11/97)
*
      tbgbf = (a(17) + a(18)*m**4 + a(19)*m**(11.d0/2.d0) + m**7)/
     &        (a(20)*m**2 + a(21)*m**7)
*
      return
      end
***
      real*8 FUNCTION tbgbdf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt mass.
* (JH 24/11/97)
*
      mx = SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*m**5*mx + m**7
      df = 4.d0*a(18)*m**3 + 5.5d0*a(19)*m**4*mx + 7.d0*m**6
      g = a(20)*m**2 + a(21)*m**7
      dg = 2.d0*a(20)*m + 7.d0*a(21)*m**6
      tbgbdf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION tbgdzf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt Z.
* (JH 14/12/98)
*
      mx = m**5*SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*mx + m**7
      df = a(117) + a(118)*m**4 + a(119)*mx
      g = a(20)*m**2 + a(21)*m**7
      dg = a(120)*m**2
      tbgdzf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION thookf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the end of the MS
* hook ( for those models that have one ) as a fraction of 
* the lifetime to the BGB
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      thookf = 1.d0 - 0.01d0*MAX(a(22)/m**a(23),a(24)+a(25)/m**a(26))
      thookf = MAX(thookf,0.5d0)
*
      return
      end
***
      real*8 FUNCTION ltmsf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the luminosity at the end of the MS
* (JH 24/11/97)
*
      ltmsf = (a(27)*m**3 + a(28)*m**4 + a(29)*m**(a(32)+1.8d0))/
     &        (a(30) + a(31)*m**5 + m**a(32))
* 
      return
      end
***
      real*8 FUNCTION lalphf(m)
      implicit none
      real*8 m,mcut,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity alpha coefficent.
* (JH 24/11/97)
*
      mcut = 2.d0
      if(m.ge.mcut)then
         lalphf = (a(33) + a(34)*m**a(36))/(m**0.4d0 + a(35)*m**1.9d0)
      else
         if(m.le.0.5d0)then
            lalphf = a(39)
         elseif(m.le.0.7d0)then
            lalphf = a(39) + ((0.3d0 - a(39))/0.2d0)*(m - 0.5d0)
         elseif(m.le.a(37))then
            lalphf = 0.3d0 + ((a(40)-0.3d0)/(a(37)-0.7d0))*(m - 0.7d0)
         elseif(m.le.a(38))then
            lalphf = a(40) + ((a(41)-a(40))/(a(38)-a(37)))*(m - a(37))
         else
            lalphf = a(41) + ((a(42)-a(41))/(mcut-a(38)))*(m - a(38))
         endif
      endif
*
      return
      end
***
      real*8 FUNCTION lbetaf(m)
      implicit none
      real*8 m,a1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity beta coefficent.
* (JH 24/11/97)
*
      lbetaf = a(43) - a(44)*m**a(45)
      lbetaf = MAX(lbetaf,0.d0)
      if(m.gt.a(46).and.lbetaf.gt.0.d0)then
         a1 = a(43) - a(44)*a(46)**a(45)
         lbetaf = a1 - 10.d0*a1*(m - a(46))
         lbetaf = MAX(lbetaf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION lnetaf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity neta exponent.
* (JH 24/11/97)
*
      if(m.le.1.d0)then
         lnetaf = 10.d0
      elseif(m.ge.1.1d0)then
         lnetaf = 20.d0
      else
         lnetaf = 10.d0 + 100.d0*(m - 1.d0)
      endif
      lnetaf = MIN(lnetaf,a(97))
*
      return
      end
***
      real*8 FUNCTION lhookf(m,mhook)
      implicit none
      real*8 m,mhook,a2,a(200)
      common /MSCFF/ a
*
* A function to evalute the luminosity at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         lhookf = 0.d0
      elseif(m.ge.a(51))then
         lhookf = MIN(a(47)/m**a(48),a(49)/m**a(50))
      else
         a2 = MIN(a(47)/a(51)**a(48),a(49)/a(51)**a(50))
         lhookf = a2*((m-mhook)/(a(51)-mhook))**0.4d0
      endif
*
      return
      end
***
      real*8 FUNCTION rtmsf(m)
      implicit none
      real*8 m,m2,rchk,a(200)
      common /MSCFF/ a
      real*8 rzamsf
      external rzamsf
*
* A function to evaluate the radius at the end of the MS
* Note that a safety check is added to ensure Rtms > Rzams
* when extrapolating the function to low masses. 
* (JH 24/11/97)
*
      m2 = a(62) + 0.1d0
      if(m.le.a(62))then
         rchk = 1.5d0*rzamsf(m)
         rtmsf = MAX(rchk,(a(52) + a(53)*m**a(55))/(a(54) + m**a(56)))
      elseif(m.ge.m2)then
         rtmsf = (a(57)*m**3+a(58)*m**a(61)+a(59)*m**(a(61)+1.5d0))/
     &           (a(60) + m**5)
      else
         rtmsf = a(63) + ((a(64) - a(63))/0.1d0)*(m - a(62))
      endif
* 
      return
      end
***
      real*8 FUNCTION ralphf(m)
      implicit none
      real*8 m,a5,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius alpha coefficent.
* (JH 24/11/97)
*
      if(m.le.0.5d0)then
         ralphf = a(73)
      elseif(m.le.0.65d0)then
         ralphf = a(73) + ((a(74) - a(73))/0.15d0)*(m - 0.5d0)
      elseif(m.le.a(70))then
         ralphf = a(74) + ((a(75)-a(74))/(a(70)-0.65d0))*(m - 0.65d0)
      elseif(m.le.a(71))then
         ralphf = a(75) + ((a(76) - a(75))/(a(71) - a(70)))*(m - a(70))
      elseif(m.le.a(72))then
         ralphf = (a(65)*m**a(67))/(a(66) + m**a(68))
      else
         a5 = (a(65)*a(72)**a(67))/(a(66) + a(72)**a(68))
         ralphf = a5 + a(69)*(m - a(72))
      endif
*
      return
      end
***
      real*8 FUNCTION rbetaf(m)
      implicit none
      real*8 m,m2,m3,b2,b3,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius beta coefficent.
* (JH 24/11/97)
*
      m2 = 2.d0
      m3 = 16.d0
      if(m.le.1.d0)then
         rbetaf = 1.06d0
      elseif(m.le.a(82))then
         rbetaf = 1.06d0 + ((a(81)-1.06d0)/(a(82)-1.d0))*(m-1.d0)
      elseif(m.le.m2)then
         b2 = (a(77)*m2**(7.d0/2.d0))/(a(78) + m2**a(79))
         rbetaf = a(81) + ((b2-a(81))/(m2-a(82)))*(m-a(82))
      elseif(m.le.m3)then
         rbetaf = (a(77)*m**(7.d0/2.d0))/(a(78) + m**a(79))
      else
         b3 = (a(77)*m3**(7.d0/2.d0))/(a(78) + m3**a(79))
         rbetaf = b3 + a(80)*(m - m3)
      endif
      rbetaf = rbetaf - 1.d0
*
      return
      end
***
      real*8 FUNCTION rgammf(m)
      implicit none
      real*8 m,m1,b1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius gamma coefficent.
* (JH 24/11/97)
*
      m1 = 1.d0
      if(m.gt.(a(88)+0.1d0))then
         rgammf = 0.d0
      else
         b1 = MAX(0.d0,a(83) + a(84)*(m1-a(85))**a(86))
         if(m.le.m1)then
            rgammf = a(83) + a(84)*ABS(m-a(85))**a(86)
         elseif(m.le.a(88))then
            rgammf = b1 + (a(89) - b1)*((m - m1)/(a(88) - m1))**a(87)
         else
            if(a(88).gt.m1) b1 = a(89)
            rgammf = b1 - 10.d0*b1*(m - a(88))
         endif
         rgammf = MAX(rgammf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION rhookf(m,mhook)
      implicit none
      real*8 m,mhook,m2,b2,a(200)
      common /MSCFF/ a
*
* A function to evalute the radius at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         rhookf = 0.d0
      elseif(m.le.a(94))then
         rhookf = a(95)*SQRT((m-mhook)/(a(94)-mhook))
      elseif(m.le.2.d0)then
         m2 = 2.d0
         b2 = (a(90) + a(91)*m2**(7.d0/2.d0))/
     &        (a(92)*m2**3 + m2**a(93)) - 1.d0
         rhookf = a(95) + (b2-a(95))*((m-a(94))/(m2-a(94)))**a(96)
      else
         rhookf = (a(90) + a(91)*m**(7.d0/2.d0))/
     &            (a(92)*m**3 + m**a(93)) - 1.d0
      endif
*
      return
      end
***
      real*8 FUNCTION lbgbf(m)
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the luminosity at the end of the 
* FGB ( for those models that have one )
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      lbgbf = (a(1)*m**a(5) + a(2)*m**a(8))/
     &        (a(3) + a(4)*m**a(7) + m**a(6))
* 
      return
      end
***
      real*8 FUNCTION lbgbdf(m)
      real*8 m,a(200)
      real*8 f,df,g,dg
      common /GBCFF/ a
*
* A function to evaluate the derivitive of the Lbgb function.
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      f = a(1)*m**a(5) + a(2)*m**a(8)
      df = a(5)*a(1)*m**(a(5)-1.d0) + a(8)*a(2)*m**(a(8)-1.d0)
      g = a(3) + a(4)*m**a(7) + m**a(6)
      dg = a(7)*a(4)*m**(a(7)-1.d0) + a(6)*m**(a(6)-1.d0)
*
      lbgbdf = (df*g - f*dg)/(g*g)
* 
      return
      end
***
      real*8 FUNCTION lbagbf(m,mhefl)
      implicit none
      real*8 m,mhefl,a4,a(200)
      common /GBCFF/ a
*
* A function to evaluate the BAGB luminosity. (OP 21/04/98)
* Continuity between LM and IM functions is ensured by setting
* gbp(16) = lbagbf(mhefl,0.0) with gbp(16) = 1.0.
*
      a4 = (a(9)*mhefl**a(10) - a(16))/(exp(mhefl*a(11))*a(16))
*
      if(m.lt.mhefl)then
         lbagbf = a(9)*m**a(10)/(1.d0 + a4*exp(m*a(11)))
      else
         lbagbf = (a(12) + a(13)*m**(a(15)+1.8d0))/(a(14) + m**a(15))
      endif
*
      return
      end
***
      real*8 FUNCTION rgbf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the GB.
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbf = a1*(lum**a(18) + a(17)*lum**a(19))
*
      return
      end
***
      real*8 FUNCTION rgbdf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the GB (as f(L)).
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbdf = a1*(a(18)*lum**(a(18)-1.d0) + 
     &            a(17)*a(19)*lum**(a(19)-1.d0))
*
      return
      end
***
      real*8 FUNCTION ragbf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the AGB.
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbf = a1*(lum**a(18) + a(17)*lum**a4)
*
      return
      end
***
      real*8 FUNCTION ragbdf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the AGB (as f(L)).
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbdf = a1*(a(18)*lum**(a(18)-1.d0) + 
     &             a(17)*a4*lum**(a4-1.d0))
*
      return
      end
***
      real*8 FUNCTION mctmsf(m)
      implicit none
      real*8 m,m525
*
* A function to evaluate core mass at the end of the MS as a 
* fraction of the BGB value, i.e. this must be multiplied by 
* the BGB value (see below) to give the actual core mass (JH 5/9/99)
*
      m525 = m**(21.d0/4.d0)
      mctmsf = (1.586d0 + m525)/(2.434d0 + 1.02d0*m525)
*
      return
      end
***
      real*8 FUNCTION mcheif(m,mhefl,mchefl)
      implicit none
      real*8 m,mhefl,mchefl,mcbagb,a3,a(200)
      common /GBCFF/ a
      real*8 mcagbf
      external mcagbf
*
* A function to evaluate core mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars  (OP 25/11/97)
*
      mcbagb = mcagbf(m)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      mcheif = MIN(0.95d0*mcbagb,(a3 + a(33)*m**a(34))**(1.d0/4.d0))
*
      return
      end
***
      real*8 FUNCTION mheif(mc,mhefl,mchefl)
      implicit none
      real*8 mc,mhefl,mchefl,m1,m2,a3,a(200)
      common /GBCFF/ a
      real*8 mbagbf
      external mbagbf
*
* A function to evaluate mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars by inverting
* mcheif
*
      m1 = mbagbf(mc/0.95d0)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      m2 = ((mc**4 - a3)/a(33))**(1.d0/a(34))
      mheif = MAX(m1,m2)
*
      return
      end
***
      real*8 FUNCTION mcagbf(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate core mass at the BAGB  (OP 25/11/97)
*
      mcagbf = (a(37) + a(35)*m**a(36))**(1.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION mbagbf(mc)
      implicit none
      real*8 mc,mc4,a(200)
      common /GBCFF/ a
*
* A function to evaluate mass at the BAGB by inverting mcagbf.
*
      mc4 = mc**4
      if(mc4.gt.a(37))then
         mbagbf = ((mc4 - a(37))/a(35))**(1.d0/a(36))
      else
         mbagbf = 0.d0
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate Mc given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         mcgbtf = ((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (1.d0/(1.d0-GB(5)))
      else
         mcgbtf = ((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (1.d0/(1.d0-GB(6)))
      endif
*
      return
      end
***
      real*8 FUNCTION lgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate L given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         lgbtf = GB(4)*(((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (GB(5)/(1.d0-GB(5))))
      else
         lgbtf = GB(3)*(((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (GB(6)/(1.d0-GB(6))))
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbf(lum,GB,lx)
      implicit none
      real*8 lum,GB(10),lx
*
* A function to evaluate Mc given L for GB, AGB and NHe stars
*
      if(lum.le.lx)then
         mcgbf = (lum/GB(4))**(1.d0/GB(5))
      else
         mcgbf = (lum/GB(3))**(1.d0/GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lmcgbf(mc,GB)
      implicit none
      real*8 mc,GB(10)
*
* A function to evaluate L given Mc for GB, AGB and NHe stars
*
      if(mc.le.GB(7))then
         lmcgbf = GB(4)*(mc**GB(5))
      else
         lmcgbf = GB(3)*(mc**GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lHeIf(m,mhefl)
      implicit none
      real*8 m,mhefl,a(200)
      common /GBCFF/ a
*
* A function to evaluate He-ignition luminosity  (OP 24/11/97)
* Continuity between the LM and IM functions is ensured with a first
* call setting lhefl = lHeIf(mhefl,0.0)
*
      if(m.lt.mhefl)then
         lHeIf = a(38)*m**a(39)/(1.d0 + a(41)*EXP(m*a(40)))
      else
         lHeIf = (a(42) + a(43)*m**3.8d0)/(a(44) + m**2)
      endif
*
      return
      end
***
      real*8 FUNCTION lHef(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the ratio LHe,min/LHeI  (OP 20/11/97)
* Note that this function is everywhere <= 1, and is only valid
* for IM stars
*
      lHef = (a(45) + a(46)*m**(a(48)+0.1d0))/(a(47) + m**a(48))
*
      return
      end
***
      real*8 FUNCTION rminf(m)
      implicit none
      real*8 m,mx,a(200)
      common /GBCFF/ a
*
* A function to evaluate the minimum radius during He-burning
* for IM & HM stars  (OP 20/11/97)
*
      mx = m**a(53)
      rminf = (a(49)*m + (a(50)*m)**a(52)*mx)/(a(51) + mx)
*
      return
      end
***
      real*8 FUNCTION tHef(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a(200)
      common /GBCFF/ a
      real*8 themsf
      external themsf
*
* A function to evaluate the He-burning lifetime.  (OP 26/11/97)
* For IM & HM stars, tHef is relative to tBGB.
* Continuity between LM and IM stars is ensured by setting
* thefl = tHef(mhefl,0.0,,0.0), and the call to themsf ensures
* continuity between HB and NHe stars as Menv -> 0.
*
      if(m.le.mhefl)then
         mm = MAX((mhefl - m)/(mhefl - mc),1.0d-12)
         tHef = (a(54) + (themsf(mc) - a(54))*mm**a(55))*
     &          (1.d0 + a(57)*EXP(m*a(56)))
      else
         mm = m**5
         tHef = (a(58)*m**a(61) + a(59)*mm)/(a(60) + mm)
      endif
*
      return
      end
***
      real*8 FUNCTION tblf(m,mhefl,mfgb)
      implicit none
      real*8 m,mhefl,mfgb,mr,m1,m2,r1,a(200)
      common /GBCFF/ a
      real*8 lheif,rminf,ragbf
      external lheif,rminf,ragbf
*
* A function to evaluate the blue-loop fraction of the He-burning
* lifetime for IM & HM stars  (OP 28/01/98)
*
      mr = mhefl/mfgb
      if(m.le.mfgb) then
         m1 = m/mfgb
         m2 = log10(m1)/log10(mr)
         m2 = max(m2,1.0d-12)
         tblf = a(64)*m1**a(63) + a(65)*m2**a(62)
      else
         r1 = 1.d0 - rminf(m)/ragbf(m,lheif(m,mhefl),mhefl)
         r1 = max(r1,1.0d-12)
         tblf = a(66)*m**a(67)*r1**a(68)
      end if
      tblf = MIN(1.d0,MAX(0.d0,tblf))
      if(tblf.lt.1.0d-10) tblf = 0.d0
*
      return
      end
***
      real*8 FUNCTION lzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a4,a5,a(200)
      common /GBCFF/ a
      real*8 lzhef
      external lzhef
*
* A function to evaluate the ZAHB luminosity for LM stars. (OP 28/01/98)
* Continuity with LHe,min for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to lzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      a5 = lzhef(mc)
      a4 = (a(69) + a5 - a(74))/((a(74) - a5)*exp(a(71)*mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      lzahbf = a5 + (1.d0 + a(72))*a(69)*mm**a(70)/
     &         ((1.d0 + a(72)*mm**a(73))*(1.d0 + a4*EXP(m*a(71))))
*
      return
      end
***
      real*8 FUNCTION rzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,rx,ry,mm,f,a(200)
      common /GBCFF/ a
      real*8 rzhef,rgbf,lzahbf
*
* A function to evaluate the ZAHB radius for LM stars. (OP 28/01/98)
* Continuity with R(LHe,min) for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to rzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      rx = rzhef(mc)
      ry = rgbf(m,lzahbf(m,mc,mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      f = (1.d0 + a(76))*mm**a(75)/(1.d0 + a(76)*mm**a(77))
      rzahbf = (1.d0 - f)*rx + f*ry
*
      return
      end
***
      real*8 FUNCTION lzhef(m)
      implicit none
      real*8 m,m15
*
* A function to evaluate Naked Helium star 'ZAMS' luminosity
*
      m15 = m*SQRT(m)
      lzhef = 1.5262d+04*m**(41.d0/4.d0)/
     &        (0.0469d0 + m**6*(31.18d0 + m15*(29.54d0 + m15)))
*
      return
      end
***
      real*8 FUNCTION rzhef(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star 'ZAMS' radius
*
      rzhef = 0.2391d0*m**4.6d0/(0.0065d0 + (0.162d0 + m)*m**3)
*
      return
      end
***
      real*8 FUNCTION themsf(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star main sequence lifetime
*
      themsf = (0.4129d0 + 18.81d0*m**4 + 1.853d0*m**6)/m**(13.d0/2.d0)
*
      return
      end
***
      real*8 FUNCTION rhehgf(m,lum,rx,lx)
      implicit none
      real*8 m,lum,rx,lx,cm
*
* A function to evaluate Helium star radius on the Hertzsprung gap 
* from its mass and luminosity. 
*
      cm = 2.0d-03*m**(5.d0/2.d0)/(2.d0 + m**5)
      rhehgf = rx*(lum/lx)**0.2d0 + 0.02d0*(EXP(cm*lum) - EXP(cm*lx))
*
      return
      end
***
      real*8 FUNCTION rhegbf(lum)
      implicit none
      real*8 lum
*
* A function to evaluate Helium star radius on the giant branch. 
*
      rhegbf = 0.08d0*lum**(3.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION lpertf(m,mew)
      implicit none
      real*8 m,mew
      real*8 b,c
*
* A function to obtain the exponent that perturbs luminosity.
*
      b = 0.002d0*MAX(1.d0,2.5d0/m)
      c = 3.d0
      lpertf = ((1.d0 + b**c)*((mew/b)**c))/(1.d0+(mew/b)**c)
*
      return
      end
***
      real*8 FUNCTION rpertf(m,mew,r,rc)
      implicit none
      real*8 m,mew,r,rc
      real*8 a,b,c,q
*
* A function to obtain the exponent that perturbs radius.
*
      a = 0.1d0
      b = 0.006d0*MAX(1.d0,2.5d0/m)
      c = 3.d0
      q = log(r/rc)
      rpertf = ((1.d0 + b**c)*((mew/b)**c)*(mew**(a/q)))/
     &         (1.d0+(mew/b)**c)
*
      return
      end
***
      real*8 FUNCTION vrotf(m)
      implicit none
      real*8 m
*
      vrotf = 330.d0*m**3.3d0/(15.d0 + m**3.45d0)
*
      return
      end
***
      real*8 FUNCTION celamf(kw,m,lum,rad,rzams,menv,fac)
      implicit none
      integer kw
      real*8 m,lum,rad,rzams,menv,fac
      real*8 lam1,lam2,m1,logm,logl
      real*8 aa,bb,cc,dd
*
* A function to estimate lambda for common-envelope.
*
      if(fac.ge.0.d0)then
*
* No fits yet for naked He stars...
*
         if(kw.gt.6)then
            celamf = 0.5d0
            goto 90
         endif
*
         if(menv.gt.0.d0)then
* Formulae for giant-like stars; also used for HG and CHeB stars close
* to the Hayashi track.
            logl = log10(lum)
            logm = log10(m)
            if(kw.le.5)then
               m1 = m
               if(kw.gt.3) m1 = 100.d0
               lam1 = 3.d0/(2.4d0 + 1.d0/m1**(3.d0/2.d0)) - 0.15d0*logl
               lam1 = MIN(lam1,0.8d0)
            else
               lam1 = -3.5d0 - 0.75d0*logm + logl
            endif
            if(kw.gt.3)then
               lam2 = MIN(0.9d0,0.58d0 + 0.75d0*logm) - 0.08d0*logl
               if(kw.lt.6)then
                  lam1 = MIN(lam2,lam1)
               else
                  lam1 = MAX(lam2,lam1)
                  lam1 = MIN(lam1,1.d0)
               endif
            endif
            lam1 = 2.d0*lam1
            if(fac.gt.0.d0)then
* Use a fraction FAC of the ionization energy in the energy balance.
               if(kw.le.3)then
                  aa = MIN(1.2d0*(logm - 0.25d0)**2 - 0.7d0,-0.5d0)
               else
                  aa = MAX(-0.2d0 - logm,-0.5d0)
               endif
               bb = MAX(3.d0 - 5.d0*logm,1.5d0)
               cc = MAX(3.7d0 + 1.6d0*logm,3.3d0 + 2.1d0*logm)
               lam2 = aa + ATAN(bb*(cc - logl))
               if(kw.le.3)then
                  dd = MAX(0.d0,MIN(0.15d0,0.15d0 - 0.25d0*logm))
                  lam2 = lam2 + dd*(logl - 2.d0)
               endif
               lam2 = MAX(lam2,1.d-2)
               lam2 = MAX(1.d0/lam2,lam1)
               if(fac.ge.1.d0)then
                  lam1 = lam2
               else
                  lam1 = lam1 + fac*(lam2 - lam1)
               endif
            endif
         endif
*
         if(menv.lt.1.d0)then
* Formula for HG stars; also reasonable for CHeB stars in blue loop.
            lam2 = 0.42d0*(rzams/rad)**0.4d0
* Alternatively:
*           lam2 = 0.3d0*(rtms/rad)**0.4d0
            lam2 = 2.d0*lam2
         endif
*
         if(menv.le.0.d0)then
            celamf = lam2
         elseif(menv.ge.1.d0)then
            celamf = lam1
         else
* Interpolate between HG and GB values depending on conv. envelope mass.
            celamf = lam2 + sqrt(menv)*(lam1 - lam2)
         endif
*
 90      continue
*
      else
         celamf = -1.d0*fac
      endif
*
      return
      end
***
C ***
C       PROGRAM colour
C c-------------------------------------------------------------c
C c
C c     Example program to show how to compute colours 
C c     using Kurucz and Bergeron (WDs) model atmosphere data. 
C c     Requires data files: 
C c         Kurucz.dat
C c         (Kurucz, 1992, Proc. IAU Symp. 149, p. 225)
C c         wdhyd.dat
C c         (Bergeron, Wesemael & Beauchamp, 1995, PASP, 107, 1047)
C c     Written by Jarrod Hurley 10/08/06. 
C c
C c-------------------------------------------------------------c
C       implicit none
C *
C       INTEGER kw
C *
C       REAL*8 z,massi,logl,logr,logte,logz,bolc
C       REAL*8 mv,mb,bminv,vmini,uminb
C       REAL*8 mv1,mv2,bminv1,bminv2
C *
C * Set metallicity. 
C *
C       z = 0.001d0
C       logz = log10(z/0.02)
C *
C * Read in the bolometric correction tables. 
C *
C       CALL ubvtab
C *
C * Mass, luminosity and radius of a star (as well as metallicity) 
C * are required to compute colours. Here are some test values: 
C * 2 Msun turn-off star, Z = 0.001
C       massi = 2.0
C       logl = 1.8355
C       logr = 0.3955
C * 1.3 Msun turn-off star, Z = 0.02
C *     massi = 1.3 
C *     logl = 0.7574    
C *     logr = 0.3369
C * 0.55 Msun cool WD 
C *     massi = 0.55
C *     logl = -3.1694   
C *     logr = -1.8712
C *
C       logte = 3.762d0 + 0.25d0*logl - 0.5d0*logr
C       if(kw.ge.10.and.kw.le.12)then
C          CALL wd2ubv(logl,logte,massi,bolc,mv,
C      &               uminb,bminv,vmini)
C       else
C          CALL lt2ubv(logl,logte,massi,logz,bolc,mv,
C      &               uminb,bminv,vmini)
C       endif
C       mb = mv + bminv
C *
C       WRITE(*,*)
C       WRITE(*,*)' Single Star:'
C       WRITE(*,*)'              Mv ',mv
C       WRITE(*,*)'              Mb ',mb
C       WRITE(*,*)'             U-B ',uminb
C       WRITE(*,*)'             B-V ',bminv
C       WRITE(*,*)'             V-I ',vmini
C *
C * Here is an example for an unresolved binary. 
C *
C       mv1 = mv
C       bminv1 = bminv
C *
C * Take a 0.55 Msun WD companion. 
C       massi = 0.55
C       logl = -3.1694   
C       logr = -1.8712
C *
C       logte = 3.762d0 + 0.25d0*logl - 0.5d0*logr
C       if(kw.ge.10.and.kw.le.12)then
C          CALL wd2ubv(logl,logte,massi,bolc,mv2,
C      &               uminb,bminv2,vmini)
C       else
C          CALL lt2ubv(logl,logte,massi,logz,bolc,mv2,
C      &               uminb,bminv2,vmini)
C       endif
C *
C * Add fluxes to get combined magnitudes. 
C       mv = -2.5d0*log10(10.d0**(-0.4d0*mv1) + 10.d0**(-0.4d0*mv2))
C       mb = -2.5d0*log10(10.d0**(-0.4d0*(mv1+bminv1)) +
C      &                  10.d0**(-0.4d0*(mv2+bminv2)))
C       bminv = mb - mv
C *
C       WRITE(*,*)
C       WRITE(*,*)' Binary Star:'
C       WRITE(*,*)'              Mv ',mv
C       WRITE(*,*)'              Mb ',mb
C       WRITE(*,*)'             B-V ',bminv
C *
C       STOP
C       END
***
      SUBROUTINE ubvtab
      implicit none
*
      INTEGER i,j,k,n
      INTEGER nzgr,ntgr,nggr
      PARAMETER(nzgr=8,ntgr=61,nggr= 11)
      INTEGER ntgr2,nggr2
      PARAMETER(ntgr2=91,nggr2=5)
*
      REAL feh
      REAL*8 zgr(nzgr),tgr(ntgr),ggr(nggr),ubv(nzgr,ntgr,nggr,5)
      COMMON /ubvdata/ zgr,tgr,ggr,ubv
      REAL*8 wtgr(ntgr2),wggr(nggr2),wubv(ntgr2,nggr2,5)
      COMMON /wubvdata/ wtgr,wggr,wubv
*
      OPEN(23,file='Kurucz.dat',
     &        form='formatted',status='old')
      do k = 1, nzgr
         do i = 1, ntgr
            do j = 1, nggr
               read(23,*)feh,tgr(i),ggr(j),(ubv(k,i,j,n),n=1,5)
            end do
            tgr(i) = log10(tgr(i))
         end do
c....... zgr=log(Z/0.02), assuming X=0.76-3*Z and Z(sun)=0.02
         zgr(k) = -log10((3.d0 + 37.425d0*10.d0**(-feh))/38.d0)
*        zgr(k) = -log10(0.07895 + 0.92105*10.0**(-feh))
      end do
      CLOSE(23)
*
      OPEN(24,file='wdhyd.dat',
     &        form='formatted',status='old')
      do j = 1,nggr2
         do i = 1,ntgr2
            read(24,*)wtgr(i),wggr(j),(wubv(i,j,k),k=1,5)
         enddo
      enddo
      do i = 1,ntgr2
         wtgr(i) = log10(wtgr(i))
      enddo
      CLOSE(24)
*
      RETURN
      END
***
      subroutine lt2ubv(logl,logt,mass,logz,bolc,mv,uminb,bminv,vmini)
      implicit none
c.... computes values of Mv, U-B, B-V and V-I for given log L, log T,
c.... mass and log(Z/0.02)
*
      integer k,ig,it,iz
      integer nzgr,ntgr,nggr
      parameter(nzgr=8,ntgr=61,nggr=11)
      integer indx
      external indx
      real*8 zgr(nzgr),tgr(ntgr),ggr(nggr),ubv(nzgr,ntgr,nggr,5)
      common /ubvdata/ zgr,tgr,ggr,ubv
      real*8 cm(5),gconst
      parameter(gconst=-10.6071d0)
      real*8 logl,logt,mass,logz,mv,uminb,bminv,vmini
      real*8 logm,logg,dg1,dg2,dt1,dt2,dz1,dz2,mbol,bolc
*
      logm = dlog10(mass)
      logg = logm + 4.d0*logt - logl + gconst
c.... find indices of log Z, log g and log T to interpolate between.
c.... don't allow extrapolation outside log Z and log g grid.
      ig = indx(logg,ggr,nggr)
      it = indx(logt,tgr,ntgr)
      iz = indx(logz,zgr,nzgr)
      dg1 = (logg - ggr(ig-1))/(ggr(ig) - ggr(ig-1))
      dg1 = max(0.d0, min(1.d0, dg1))
      dg2 = 1.d0 - dg1
      dt1 = (logt - tgr(it-1))/(tgr(it) - tgr(it-1))
      dt2 = 1.d0 - dt1
      dz1 = (logz - zgr(iz-1))/(zgr(iz) - zgr(iz-1))
      dz1 = max(0.d0, min(1.d0, dz1))
      dz2 = 1.d0 - dz1
      do k = 1, 5
         cm(k) = ((ubv(iz,it,ig,k)*dg1 + ubv(iz,it,ig-1,k)*dg2)*dt1
     &            + (ubv(iz,it-1,ig,k)*dg1 +
     &            ubv(iz,it-1,ig-1,k)*dg2)*dt2)*dz1 +
     &           ((ubv(iz-1,it,ig,k)*dg1 +
     &            ubv(iz-1,it,ig-1,k)*dg2)*dt1 +
     &           (ubv(iz-1,it-1,ig,k)*dg1 +
     &            ubv(iz-1,it-1,ig-1,k)*dg2)*dt2)*dz2
      enddo
      mbol = 4.75d0 - 2.5d0*logl
      bolc = cm(1)
      mv = mbol - bolc
      uminb = cm(2)
      bminv = cm(3)
      vmini = cm(4) + cm(5)
*
      return
      end
***
      integer function indx(ax,xx,nx)
c.....finds index of ax in monotonously increasing or decreasing array xx
      implicit none
      integer nx,j,jl,jh
      real*8 ax,xx(nx),sx
*
      sx = xx(nx) - xx(1)
      jl = 1
      jh = nx
 1    if (jh-jl.gt.1) then
         j = (jh + jl)/2
         if ((ax-xx(j))*sx.gt.0.d0) then
            jl = j
         else
            jh = j
         end if
         goto 1
      end if
      indx = jh
*
      return
      end
***
      subroutine wd2ubv(logl,logt,mass,bolc,mv,uminb,bminv,vmini)
c.... computes values of Mv, U-B, B-V and V-I for given log L, log T,
c.... mass and log(Z/0.02)
      implicit none
      integer k,ig,it
      integer ntgr,nggr
      parameter(ntgr=91,nggr=5)
      integer indx
      external indx
      real*8 tgr(ntgr),ggr(nggr),ubv(ntgr,nggr,5)
      common /wubvdata/ tgr,ggr,ubv
      real*8 cm(5),gconst
      parameter(gconst = -10.6071d0)
      real*8 mbol,bolc,mv,uminb,bminv,vmini
      real*8 logm,logg,dg1,dg2,dt1,dt2,mass,logt,logl
*
      logm = dlog10(mass)
      logg = logm + 4.d0*logt - logl + gconst
c.... find indices of log g and log T to interpolate between.
c.... don't allow extrapolation outside log g grid.
      ig = indx(logg,ggr,nggr)
      it = indx(logt,tgr,ntgr)
      dg1 = (logg - ggr(ig-1))/(ggr(ig) - ggr(ig-1))
      dg1 = MAX(0.d0,MIN(1.d0,dg1))
      dg2 = 1.d0 - dg1
      dt1 = (logt - tgr(it-1))/(tgr(it) - tgr(it-1))
      dt2 = 1.d0 - dt1
      do k = 1,5
         cm(k) = (ubv(it,ig,k)*dg1 + ubv(it,ig-1,k)*dg2)*dt1
     &         + (ubv(it-1,ig,k)*dg1 + ubv(it-1,ig-1,k)*dg2)*dt2
      enddo
      mbol = 4.75d0 - 2.5d0*logl
      bolc = cm(1)
      mv = mbol - bolc
      uminb = cm(2)
      bminv = cm(3)
      vmini = cm(4) + cm(5)
      return
      end
***
