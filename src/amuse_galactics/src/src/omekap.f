      subroutine omekap(r,fom,fka)

      include 'parameters'

      dimension rr(nmax),om(nmax),om2(nmax),psi(nmax),
     +                   ak(nmax),ak2(nmax)
      dimension psirr(nmax)
      data ifirst /0/
      save ifirst, a, b
      save rr,om,om2,ak,ak2,nreq

c     read omega, potential and potential'' from the frequencies file,
c     spline-fit pot(r), and then use the 2nd derivatives of the
c     potential in the table of kappa. these can then in turn be
c     spline-fitted. should arrange for the boundary condition at r=0 to
c     be even for pot, omega and kappa.

      if (ifirst.eq.0) then
        open(13,file='freqdbh.dat',status='old',err=999)
        read(13,*)
        read(13,*)
        read(13,*)
        nreq=0
c read only every other step, since remainder was linearly 
c interpolated anyway
 77     read(13,*,end=99) rr(nreq+2),t,t,t,vc,t,t,psi(nreq+2)
     *       ,psirr(nreq+2)
        om(nreq+2)=vc/rr(nreq+2)
        nreq=nreq+1
        read(13,*,end=99) 
        goto 77
c extrapolate omega (linearly) and potential (quadratically) to zero
99      close(13)
        om(1)=2*om(2)-om(3)
        psi(1)=(4*psi(2)-psi(3))/3.
        psirr(1)=(4*psirr(2)-psirr(3))/3.
        rr(1)=0.
        call splined(rr,om,nreq,1.e32,1.e32,om2)
c calculate epicycle frequencies
        ak(1)=2*om(1)
        do ir=2,nreq
          ak(ir)=sqrt(max(0.,psirr(ir)+3*om(ir)**2))
        enddo
        call splined(rr,ak,nreq,0.,1.e32,ak2)
c a and b are coefficients of potential -a/r -b/r^3 +c 
c which give omeka and kappa as observed
c at the outer radius of the model---these are used to extrapolate if needed.
c in this potential om^2=a/r^3+3b/r^5, kap^2=a/r^3-3b/r^5 
        a=rr(nreq)**3/2.*(om(nreq)**2+ak(nreq)**2)
        b=rr(nreq)**5/6.*(om(nreq)**2-ak(nreq)**2)
        c=a/rr(nreq)**3+b/rr(nreq)**5+psi(nreq)
        open(11,file='omekap.dat',status='new',err=23)
        do ir=1,nreq
          write(11,'(7g16.8)') rr(ir),om(ir),om2(ir),ak(ir),ak2(ir),
     +                          psi(ir),psirr(ir)
       enddo
        do ir=1,10
           rx=rr(nreq)*(1+0.1*ir)
           r3=a/rx**3
           r5=3*b/rx**5
           fom=sqrt(r3+r5)
           fka=sqrt(r3-r5)
           write(11,'(7g16.8)') rx,fom,0.,fka,0.,
     +                          c-(r3+r5/3.)*rx*rx,-2*r3-4*r5
           enddo
        close(11)
        write(0,*) 'wrote a table of frequencies in file omekap.dat.'
 23     continue
        ifirst=1
      endif

      if (r.gt.rr(nreq)) then
         r3=a/r**3
         r5=3*b/r**5
         fom=sqrt(r3+r5)
         fka=sqrt(r3-r5)
      else
         call splintd(rr,om,om2,nreq,r,fom)
         call splintd(rr,ak,ak2,nreq,r,fka)
      endif
      return
999   write(0,*) 'no file with frequencies 
     +   corresponding to the potential!'
      stop
      end
