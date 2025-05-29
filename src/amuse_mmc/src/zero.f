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
