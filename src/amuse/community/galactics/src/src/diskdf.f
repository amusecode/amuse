      program diskdf

      include 'commonblocks'

      dimension drat(0:nrmax),
     +          dz2rat(0:nrmax),
     +          fzrat(0:nrmax)
      dimension d0rat(0:nrmax),d1rat(0:nrmax),d2rat(0:nrmax)
     +     ,d3rat(0:nrmax),
     +          d4rat(0:nrmax)
      dimension fplot(0:10*nrmax),rplot(0:10*nrmax)
      include 'equivalence'
      character*60 filename

      common /fileroot/ filename

      write(0,*) 'Central velocity dispersion, scale length of sigR**2?'
      read(*,*) sigr0,disksr
      write(0,*) 'number of radial steps for correction fns (min. 6)?'
      read(*,*) nrspl
      write(0,*) 'number of iterations?'
      read(*,*) niter
      filename='dbh.dat'
      call readharmfile(filename,ibuf,jbuf,kbuf)
      psi00=pot(0.0,0.0)
      rfid = 2.5*rdisk

      call omekap(rfid,fom,fka)
      sigr = sqrt(sigr2(rfid))
      sigden = diskdensf(rfid,0.0)*2.0*zdisk
c
c a close estimate of the disk surface density
c
      sigrcrit = 3.36*sigden/fka
      qtoomre = sigr/sigrcrit
      write(0,*) 'Toomre Q = ',qtoomre, ' at R = 2.5 R_d'
      open(file='toomre2.5',status='replace',unit=11)
      write(11,*) qtoomre
      close(11)

c first call initializes the spline fitting functions by the way
      call omekap(outdisk,fom,fka)
      if (sqrt(sigr2(outdisk)).gt.drtrunc*fka) then
         write(0,*) '** WARNING **: you have asked for an outer vel'
         write(0,*) '  dispersion which is larger than that used in the'
         write(0,*) '  disk truncation. Strange things may happen.'
         write(0,*) sqrt(sigr2(outdisk))/(drtrunc*fka),fka
      endif
      drspl = (outdisk + 2*drtrunc)/nrspl
      do ir=0,nrspl
         rr(ir)=drspl*ir
         fdrat(ir)=1
         fszrat(ir)=1
      enddo
      call splined(rr(0),fdrat(0),nrspl+1,1.e32,1.e32,drat2(0))
      call splined(rr(0),fszrat(0),nrspl+1,1.e32,1.e32,szrat2(0))

      open(12,file='cordbh.dat',status='unknown')

c Make niter iterations on the correction functions for midplane density
c and vertical velocity dispersion.
      do icor=1,niter
         do ir=0,nrspl
            r=rr(ir)
            call omekap(r,fom,fka)
            sigr=sqrt(sigr2(r))
            sigz=sqrt(sigz2(r))
            vc=r*fom
            d0=0
            dz=0
            dz2=0
            dz3=0
            dsz=0
c     integrate analytically over vr and vz, numerically over vphi
            dvr=0.1*sigr
            do ivt=1,101
               vt=vc+(ivt-51)*dvr
               df=diskdf5intez(vt,r,0.0)
               c=0.333333333333*(mod(ivt,2)*2+2)
               d0=d0+c*dvr*df
               df=diskdf5intez(vt,r,zdisk)
               dz2=dz2+c*dvr*df
            enddo
            drat(ir)=d0/diskdensf(r,0.)
            dz2rat(ir)=dz2/diskdensf(r,zdisk)
            fzrat(ir)=
     +           log(diskdensf(r,0.)/diskdensf(r,zdisk))/log(d0/dz2)
         enddo
         write(0,*)
         nrd=max(1,nrspl/10)
         write(0,'(''      r: '',11f6.3)') (rr(ir),ir=0,nrspl,nrd)
         write(0,'(f6.3,''  d'',11f6.3)') 0.,(drat(ir),ir=0,nrspl,nrd)
         write(0,'(f6.3,''  d'',11f6.3)') zdisk,
     +        (dz2rat(ir),ir=0,nrspl,nrd)

         do ir=0,nrspl
            fdrat(ir)=fdrat(ir)/drat(ir)
            fszrat(ir)=fszrat(ir)/fzrat(ir)
         enddo
         call splined(rr(0),fdrat(0),nrspl+1,1.e32,1.e32,drat2(0))
         call splined(rr(0),fszrat(0),nrspl+1,1.e32,1.e32,szrat2(0))
      enddo
c Finally, write out a more detailed table of residuals.
      do ir=0,nrspl
         r=rr(ir)
         call omekap(r,fom,fka)
         sigr=sqrt(sigr2(r))
         sigz=sqrt(sigz2(r))
         vc=r*fom
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
c     integrate analytically over vr and vz, numerically over vphi
         dvr=0.1*sigr
         do ivt=1,101
            vt=vc+(ivt-51)*dvr
            df=diskdf5intez(vt,r,0.)
            c=0.333333333333*(mod(ivt,2)*2+2)
            d0=d0+c*dvr*df
            df=diskdf5intez(vt,r,0.5*zdisk)
            d1=d1+c*dvr*df
            df=diskdf5intez(vt,r,zdisk)
            d2=d2+c*dvr*df
            df=diskdf5intez(vt,r,1.5*zdisk)
            d3=d3+c*dvr*df
            df=diskdf5intez(vt,r,2*zdisk)
            d4=d4+c*dvr*df
         enddo
         d0rat(ir)=d0/diskdensf(r,0.)
         d1rat(ir)=d1/diskdensf(r,0.5*zdisk)
         d2rat(ir)=d2/diskdensf(r,zdisk)
         d3rat(ir)=d3/diskdensf(r,1.5*zdisk)
         d4rat(ir)=d4/diskdensf(r,2*zdisk)
      enddo
      write(0,*)
      nrd=max(1,(nrspl-1)/10+1)
      write(0,'(''      r: '',11f6.3)') (rr(ir),ir=0,nrspl,nrd)
      write(0,'(f6.3,''  d'',11f6.3)') 0.,(d0rat(ir),ir=0,nrspl,nrd)
      write(0,'(f6.3,''  d'',11f6.3)') 
     +     0.5*zdisk,(d1rat(ir),ir=0,nrspl,nrd)
      write(0,'(f6.3,''  d'',11f6.3)') zdisk,(d2rat(ir),ir=0,nrspl,nrd)
      write(0,'(f6.3,''  d'',11f6.3)') 
     +     1.5*zdisk,(d3rat(ir),ir=0,nrspl,nrd)
      write(0,'(f6.3,''  d'',11f6.3)') 
     +     2*zdisk,(d4rat(ir),ir=0,nrspl,nrd)
      write(12,'(''#'', 2g17.7,1x,i4)') sigr0,disksr,nrspl
      do i=0,nrspl
         write(12,'(1x,3g17.7)') rr(i),fdrat(i),fszrat(i)
      enddo
      close(12)
      
      end
