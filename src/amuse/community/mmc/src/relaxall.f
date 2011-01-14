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
*cello            print*,'escape- nup,tescp(l) = ',l,nup,tescp(l)
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
