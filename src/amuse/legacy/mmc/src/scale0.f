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
*