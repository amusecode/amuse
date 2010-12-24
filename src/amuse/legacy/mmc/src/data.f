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
      real*8 xtemp(3,nmax),vtemp(3,nmax)
      real*4 ran2
*
      integer i,k,n,ntt,nss,nbb,ibb,ii,ixb,ibx,im
*
      data g1,g2,g3,g4 /0.19d0,1.55d0,0.05d0,0.6d0/
      DATA g5,g6,g7,g8 /0.75,0.04,0.25,1.04/
      data lambda,chi,delta,eta,pminl/28.d0,0.75d0,45.d0,2.5d0,1.d0/
*
      print*,'Entered data'

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
                     write(6,*) 'bins:',i,sm1,sm2,sm2in,qbirth,qin,
     &                    ebirth,
     &                            ein,pbirth,pin,rperi1,rperi,rho,rhop
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
      print*,'cases 3,5,6'
      if (imodel.eq.3) then
         call king
      elseif (imodel.eq.5) then
         call poly
      elseif (imodel.eq.6) then
         if (w0.lt.0.d0.or.w0.gt.1.d0) then
            print*,'bad w0 in data.f, stopping',w0
            stop
         endif
         print*,'calling plumix',body(1),body(n),n
         call flush(6)
         call plumix(body,n,xtemp,vtemp,w0)
         print*,(xtemp(k,1),k=1,3),(vtemp(k,2),k=1,3)
cFill the required arrays
         do i = 1,n
            do k = 1,3
               x(i,k) = xtemp(k,i)
               xdot(i,k) = vtemp(k,i)
            enddo
         enddo
cFollowing ensures that rplum has same significance (when S=0) as for Plummer model
         b1 = 1.5d0*twopi/16.0d0
         rtidkg = rplum*b1/sqrt(2.d0**(2.d0/3.d0) - 1.d0)
         rtid = rtidkg
c         stop
      else
         print*,'imodel = ',imodel,' stopping'
         stop
      endif
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
