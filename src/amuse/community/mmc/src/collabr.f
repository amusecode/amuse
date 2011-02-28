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

      write(6,*)   '  collabr alive'
      call flush(6)

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
        write(6,*)   '  collabr alive', r(im), ikind(i),rtid
        call flush(6)

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
