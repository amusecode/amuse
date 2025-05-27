program read_snapshot

   !------------------------------------------------------------------
   ! This Fortran 90 routine shows how a snapshot of GADGET
   ! which is distributed into several files may be read in.
   ! In this example, coodinates of particles type 1 are read.
   !
   ! Each dump consists of several files, each with a variable 
   ! number of particles. However, these files all have the same format. 
   ! The sequence of information in the files is as follows:
   !
   ! Header
   ! Positions
   ! Velocities
   ! Particle ID's
   !  
   ! Masses (for particle types with vanishing mass entry in the header)
   ! Internal Energy (Gas)
   ! Density (Gas)
   !
   !-------------------------------------------------------------------


   include 'globals.h'
   integer :: NFILES        

   character*200 fileroot,filename, fnumber

   integer*4 npfile(1:6), nall(1:6)
   real*8    massarr(1:6)
   real*8    a
   real*8    z
   integer*4 unused(34)
   integer*4 fn,i,nstart,flag_sfr,flag_feedback
   integer*4 N,Ntot,Nwithmass,Ntotwithmass,nmass


   real*4,allocatable    :: rpos(:,:),rvel(:,:),rmass(:)
   integer*4,allocatable :: id(:)
   real :: mscale,vscale,tscale,lscale
   call initmem(nbodsmax,nsphmax,ncells)

   unitm_in_msun=1.e9
   unitl_in_kpc=1.
   timescale=1.488e19*unitl_in_kpc*sqrt(unitl_in_kpc/unitm_in_msun)
   lengthscale=3.086e21*unitl_in_kpc
   velscale=lengthscale/timescale

   print*,'single gadget file conversion'
   print*,'fileroot?'
   read*,fileroot
   print*,'mscale (in msun)?'
   read*,mscale
   mscale=mscale/unitm_in_msun
   print*,'vscale (in km/s)?'
   read*,vscale
   vscale=1.e5*vscale/velscale
   print*,'tscale (in Myr)?'
   read*,tscale
   tscale=tscale/timescale*year*1.e6
   print*,'m,v,t factor:',mscale,vscale,tscale
   lscale=1
   
   filename=trim(fileroot)

   print *,'opening...  '//trim(filename)

   ! now, read in the header

   open (1, file=filename, form='unformatted')
   read (1) npfile, massarr ,a, z, flag_sfr,flag_feedback, nall, unused

   do i=1,6
      print *,'Type,npfile,massarr',i, nall(i), massarr(i)
   end do
   
   Ntot= sum(nall)
   N=sum(npfile)

   Nwithmass=0
   do i=1,6
    if(massarr(i).EQ.0) Nwithmass=Nwithmass+npfile(i)  
   enddo

   allocate(rpos(1:3,1:N),rvel(1:3,1:N),id(1:N))
   if(Nwithmass.gt.0) allocate(rmass(1:Nwithmass))

   read (1) rpos
   read (1) rvel
   read (1) id
   if(Nwithmass.gt.0) read (1) rmass
   close (1)
 
   print *,'Done with reading.'
 
   nbodies=Ntot
   nsph=npfile(1)
   nstar=sum(npfile(3:6))
   nbh=npfile(6)
   totptag=nbodies
   tnow=a*tscale   
   massres=0.

   if(nbodies.GT.nbodsmax.OR.nsph.GT.nsphmax) then
    print*,"particle overflow:",nbodies,nsph
    stop
   endif

   outputfile=fileroot
   firstsnap=0
   uentropy=.FALSE.
   output=0
   output(1:3)=1
   
   do i=1,nbodies  
   pos(i,1:3)=rpos(1:3,i)*lscale
   vel(i,1:3)=rvel(1:3,i)*vscale
   enddo
   nmass=1
   nstart=1
   do i=1,6
    if(massarr(i).eq.0) then
      mass(nstart:nstart+npfile(i)-1)=rmass(nmass:nmass+npfile(i)-1)
      nstart=nstart+npfile(i)
      nmass=nmass+npfile(i)
    else
      mass(nstart:nstart+npfile(i)-1)=massarr(i)
      nstart=nstart+npfile(i)  
    endif
   enddo
   mass=mass*mscale
  
   if(Nwithmass.gt.0) deallocate(rmass)
   deallocate(rpos,rvel,id)

   print*,'fill ethermal (T= 10^4 K)'
   output(11)=0   
   if(nsph.gt.0) then
    output(11)=1
    ethermal(1:nsph) =1.e4/4.86e5 ! 10^4 K
   endif
   
   print*,'fill tform (age=12.e9 yr)'
   output(5)=1
   tform(1:nbodies)=0.
   allocate(rmass(nstar))
   call random_number(rmass)
   tform(nbodies-nstar+1:nbodies)=-rmass*800
deallocate(rmass)

   call outbods(trim(fileroot)//".new")
end program




