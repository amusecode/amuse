      program plotforce

      real rplot(0:200),zplot(0:200),vc(0:200),
     :     vd(0:200),vh(0:200),vb(0:200),fz(0:200)
      character*60 toplbl,filename
      integer*4 ibuf1(15)
      character ans

      open(20,file='scales',status='old',err=5)
      read(20,*) rscale
      write(*,*) 'Will normalize model so that vcirc(rsun)=200.'
      write(*,*) 'and will write this to "scales"'
      close(20)
      goto 6

 5    write(*,*) 'No files "scales" found.'
      write(*,*) 'Will use default length scale [4kpc] instead.'
      write(*,*) 'Will normalize forces so that vcirc(8kpc)=200.'
      rscale=4.

 6    continue

      write(*,*) 'Will plot rotation curve and Kz.'
      write(*,*) 'For rotation curve: radius bin [kpc], no. bins?'
      read(*,*) drplot,nrplot
      write(*,*) 'For Kz:  z-bin [kpc], no. bins?'
      read(*,*) dzplot,nzplot

c inititially store the radial forces in the velocity arrays

      filename='h.dat'
      call readharmfile(filename,ibuf1)

      write(*,*) 'Calculating halo rotation curve...'
      do ir=0,nrplot
         rp=ir*drplot/rscale
         rplot(ir)=rp*rscale
         call force(rp,0.,vh(ir),fz(0),p)
!	 print*,rp,vh(ir),fz(0),p
      enddo

      filename='b.dat'
      call readharmfile(filename,ibuf1)
      write(*,*) 'Calculating bulge rotation curve...'
      do ir=0,nrplot
         rp=rplot(ir)/rscale
         call force(rp,0.,vb(ir),fz,p)
      enddo

      filename='dbh.dat'
      call readharmfile(filename,ibuf1)
c      call readdiskdf(filename,ibuf1)      
c (this also reads in the dbh.dat file)
      write(*,*) 'Calculating total rotation curve...'
      do ir=0,nrplot
         rp=rplot(ir)/rscale
         call force(rp,0.,vc(ir),fz,p)
      enddo

      rsun=8.
      call force(rsun/rscale,0.,fr,fz,p)
      vcsun=sqrt(max(0.,-rsun/rscale*fr))
      vscale=200./vcsun
      emscale=rscale*vscale**2/4.3e-6

      write(*,*) 'Scales for mass, distance, speed:'
      write(*,*) emscale,' Msun, ',rscale,' kpc, ',vscale,' km/s.'
      densscale=emscale/rscale**3
      dfscale=densscale/vscale**3
      open(16,file='scales',status='unknown')
      write(16,*) rscale,'    radius scale [kpc]'
      write(16,*) vscale,'    speed scale [km/s]'
      write(16,*) emscale,'    mass scale [Msun]'
      close(16)
      vsun=vcsun*vscale
      write(*,*) 'Circular speed at R=',rsun,' kpc: ',vsun,' km/s.'

c now turn forces into circular speeds.
      vmax=-1e31
      write(*,*) 'Calculating disk rotation curve...'
      do ir=0,nrplot
         vd(ir)=vc(ir)-vh(ir)-vb(ir)
         rp=rplot(ir)/rscale
         vb(ir)=sqrt(max(0.,-rp*vb(ir)))*vscale
         vh(ir)=sqrt(max(0.,-rp*vh(ir)))*vscale
         vd(ir)=sqrt(max(0.,-rp*vd(ir)))*vscale
         vc(ir)=sqrt(max(0.,-rp*vc(ir)))*vscale
         vmax=max(vc(ir),vmax)
      enddo

      fmax=-1e31
      do iz=0,nzplot
         zp=iz*dzplot
         zplot(iz)=zp
         call force(rsun/rscale,zp/rscale,fr,fz(iz),p)
         fz(iz)=-fz(iz)*vscale**2/rscale/1000.
         fmax=max(fmax,fz(iz))
      enddo

 1    call pgbeg(0,'?',1,2)
      write(*,*)
      call pgsch(2.0)
      call pgenv(0.,rplot(nrplot),0.,vmax*1.15,0,0)
      call pgline(nrplot+1,rplot(0),vb(0))
      call pgline(nrplot+1,rplot(0),vd(0))
      call pgline(nrplot+1,rplot(0),vh(0))
      call pgline(nrplot+1,rplot(0),vc(0))
      call pgpt(1,8.,200.,12)
      call pglabel('R [kpc]','Vcirc [km/s]','ROTATION CURVES')
      
      call pgenv(0.,zplot(nzplot),0.,fmax*1.15,0,0)
      call pgline(nzplot+1,zplot(0),fz(0))
      write(toplbl,'(''z-FORCE AT R='',f6.3,'' kpc'')') rsun
      call pglabel('z [kpc]','Kz [(km/s)**2 /pc]',toplbl)
      call pgpt(1,1.1,1.9,12)
      call modstamp
      call pgend

      write(*,*) 'Same plot, different device?'
      read(*,'(a)') ans
      if (ans.ne.'n' .and. ans.ne.'N') goto 1
      end

