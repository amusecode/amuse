
      program plothalo

      external halodens,diskdensf,bulgedens
      character*60 filename
      character ans

      integer*4 ibuf1(15)

      filename='dbh.dat'

c      call readdiskdf(filename,ibuf1)      
      call readharmfile(filename,ibuf1)      

      open(20,file='scales',status='old',err=5)
      read(20,*) rscale
      close(20)
      goto 6

 5    write(*,*) 'No files "scales" found.'
      write(*,*) 'Will use default length and speed scales instead.'
      rscale=1.
 6    continue
      write(*,*) 'Scale for distances: ',rscale,' kpc.'

      write(*,*) 'Will make panels for inner and outer halo density.'
      write(*,*) 'Inner panel: bin size [kpc], no. of bins (max 100)?'
      read(*,*) din,nin
      write(*,*) 'Outer panel: bin size [kpc], no. of bins (max 100)?'
      read(*,*) dout,nout
      write(*,*) 'Will also make a panel for bulge density.'
      write(*,*) 'Bulge panel: bin size [kpc], no. of bins (max 100)?'
      read(*,*) dbu,nbu

 1    call pgbeg(0,'?',1,1)
      write(*,*)
      call pgvport(0.05,0.32,0.6,0.93)
      call contourden(nin,din,halodens,1.,rscale)
      call pglabel('R [kpc]','z [kpc]','Inner Halo Density')
      call pgvport(0.37,0.64,0.6,0.93)
      call contourden(nbu,dbu,bulgedens,1.,rscale)
      call pglabel('R [kpc]','z [kpc]','Bulge Density')
      call pgvport(0.69,0.96,0.6,0.93)
      call contourden(nout,dout,halodens,1.,rscale)
      call pglabel('R [kpc]','z [kpc]','Outer Halo Density')
      call pgvport(0.05,0.96,0.1,0.42)
      call contourden(nin,din,diskdensf,0.3,rscale)
      call pglabel('R [kpc]','z [kpc]','Disk Density')
      call modstamp
      call pgend
      write(*,*) 'Density contour levels jump by factor 2.'
      write(*,*) 'Same plot, different device?'
      read(*,'(a)') ans
      if (ans.ne.'n' .and. ans.ne.'N') goto 1
      end
