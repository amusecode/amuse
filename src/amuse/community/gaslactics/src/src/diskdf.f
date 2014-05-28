c to fix: interpolation in rcirc does not seem to work well. after 
c that, hope that things come together!

! fip: changed central sig0 to be equal (except for a factor sigroz) to
! central z dispersion

      program diskdf
      parameter(nrmax=1000)
      parameter (pi=3.1415926535)
      dimension drat(0:nrmax),
     +          dz2rat(0:nrmax),
     +          fzrat(0:nrmax)
      dimension d0rat(0:nrmax),d1rat(0:nrmax),d2rat(0:nrmax),d3rat(0:nrmax),
     +          d4rat(0:nrmax)
      dimension fplot(0:10*nrmax),rplot(0:10*nrmax)
      dimension ibuf(100)
      character*60 filename

      common /fileroot/ filename
      common /potconstants/ apot(20,0:20000), frad(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /diskpars/ sigr0, disksr, nrdisk
      common /splines/ rr(0:nrmax),fdrat(0:nrmax),drat2(0:nrmax),
     +                 fszrat(0:nrmax),szrat2(0:nrmax), nrspl

c      write(0,*) 'Root name for potential?'
c      read(*,'(a60)') filename
      write(0,*) 'Central velocity dispersion, scale length of sigR**2?'
      read(*,*) sigroz,disksr
      write(0,*) 'number of radial steps for correction fns (min. 6)?'
      read(*,*) nrspl
      write(0,*) 'number of iterations?'
      read(*,*) niter


      filename='dbh.dat'
      call readharmfile(filename,ibuf)

      psi00=pot(0.0,0.0)
      rmax = nrdisk*dr

      sigr0=sigroz*sqrt(3.1415*zdisk*diskdensf(0.,0.)*2.0*zdisk)

      drspl = (outdisk + 2*drtrunc)/nrspl
      do ir=0,nrspl
        rr(ir)=drspl*ir
        fdrat(ir)=1
        fszrat(ir)=1
      enddo
c      drspl=(outdisk-2*drtrunc)/(nrspl-4)
c      do ir=0,nrspl-4
c        rr(ir)=drspl*ir
c        fdrat(ir)=1
c        fszrat(ir)=1
c        enddo
c      do ir=3,0,-1
c        rr(nrspl-ir)=outdisk+(2-ir)*drtrunc
c        fdrat(nrspl-ir)=1
c        fszrat(nrspl-ir)=1
c     enddo
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
c integrate analytically over vr and vz, numerically over vphi
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
        fzrat(ir)=log(diskdensf(r,0.)/diskdensf(r,zdisk))/log(d0/dz2)
c        write(*,*) 'density=',d0
        enddo
       write(0,*)
       nrd=max(1,nrspl/10)
       write(0,'(''      r: '',11f6.3)') (rr(ir),ir=0,nrspl,nrd)
       write(0,'(f6.3,''  d'',11f6.3)') 0.,(drat(ir),ir=0,nrspl,nrd)
       write(0,'(f6.3,''  d'',11f6.3)') zdisk,(dz2rat(ir),ir=0,nrspl,nrd)
       do ir=0,nrspl
        fdrat(ir)=fdrat(ir)/drat(ir)
        fszrat(ir)=fszrat(ir)/fzrat(ir)
        enddo
       call splined(rr(0),fdrat(0),nrspl+1,1.e32,1.e32,drat2(0))
       call splined(rr(0),fszrat(0),nrspl+1,1.e32,1.e32,szrat2(0))
c       do ir=0,nrspl
c        r = rr(ir)
c        call splintd(rr(0),fdrat(0),drat2(0),nrspl+1,r,fd)
c        call splintd(rr(0),fszrat(0),szrat2(0),nrspl+1,r,fz)
c        write(11,'(x,5g12.5)') rr(ir),fdrat(ir), fd, fszrat(ir), fz
c        enddo
c       do ir=0,nrspl*10-1
c        r=rr(ir/10)+(rr(ir/10+1)-rr(ir/10))*0.1*mod(ir,10)
c        call splintd(rr(0),fdrat(0),drat2(0),nrspl+1,r,fd)
c        call splintd(rr(0),fszrat(0),szrat2(0),nrspl+1,r,fz)
c        write(11,'(x,3g12.5)') r,fd,fz
c        enddo
       print*,'sigr0,sigz0:', sigr0,sqrt(sigz2(0.))
       sigr0=sigroz*sqrt(sigz2(0.))
       enddo

      rfid = 2.5*rdisk
      call omekap(rfid,fom,fka)
      sigr = sqrt(sigr2(rfid))
      sigden = diskdensf(rfid,0.0)*2.0*zdisk
      print*,'central surfdens:',diskdensf(0.,0.0)*2.0*zdisk
c
c a close estimate of the disk surface density
c
      sigrcrit = 3.36*sigden/fka
      qtoomre = sigr/sigrcrit
      write(0,*) 'Toomre Q = ',qtoomre, ' at R = 2.5 R_d'
c first call initializes the spline fitting functions by the way
      call omekap(outdisk,fom,fka)

      write(0,'('' ==> Outer dispersion, truncation dispersion:'',2g10.3)')
     +  sqrt(sigr2(outdisk)),drtrunc*fka
      if (sqrt(sigr2(outdisk)).gt.drtrunc*fka) then
        write(0,*) '** WARNING **: you have asked for an outer velocity'
        write(0,*) '  dispersion which is larger than that used in the'
        write(0,*) '  disk truncation. Strange things may happen.'
        endif



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
c integrate analytically over vr and vz, numerically over vphi
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
       write(0,'(f6.3,''  d'',11f6.3)') 0.5*zdisk,(d1rat(ir),ir=0,nrspl,nrd)
       write(0,'(f6.3,''  d'',11f6.3)') zdisk,(d2rat(ir),ir=0,nrspl,nrd)
       write(0,'(f6.3,''  d'',11f6.3)') 1.5*zdisk,(d3rat(ir),ir=0,nrspl,nrd)
       write(0,'(f6.3,''  d'',11f6.3)') 2*zdisk,(d4rat(ir),ir=0,nrspl,nrd)
       write(12,'(''#'', 2g17.7,x,i4)') sigr0,disksr,nrspl
       do i=0,nrspl
          write(12,'(x,3g17.7)') rr(i),fdrat(i),fszrat(i)
       enddo
c       close(11)
       close(12)



c make plots of the correction functions.
       
       call pgbeg(0,'/ps',2,2)
       call pgsch(2.)
       write(*,*)

       do ir=0,nrspl*10-1
          rplot(ir)=rr(ir/10)+(rr(ir/10+1)-rr(ir/10))*0.1*mod(ir,10)
       enddo
       do ir=0,nrspl*10-1
          call splintd(rr(0),fdrat(0),drat2(0),nrspl+1,rplot(ir),fplot(ir))
       enddo
       
       call pgenv(0.,rplot(nrspl*10-1),0.,2.,0,0)
       call modstamp
       call pgline(nrspl*10,rplot(0),fplot(0))
       call pglabel('R',' ','Surfden Correction')
       
       do ir=0,nrspl*10-1
          call splintd(rr(0),fszrat(0),szrat2(0),nrspl+1,rplot(ir),fplot(ir))
       enddo
       
       call pgenv(0.,rplot(nrspl*10-1),0.,2.,0,0)
       call pgline(nrspl*10,rplot(0),fplot(0))
       call pglabel('R',' ','sigz**2 Correction')

       call modstamp

       call pgend
       

      end
