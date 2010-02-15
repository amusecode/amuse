
program read_snapshot

   include 'globals.h'
   character*200 :: filename
   character*200 :: outfilename
   integer :: np,i,ioer
   real :: pmass,x(3),v(3)
   real :: mscale,vscale,tscale,lscale
   real :: epsg,epsdm,epsstar,ether
   call initmem(nbodsmax,nsphmax,ncells)

!   unitm_in_msun=1.e9
!   unitl_in_kpc=1.
!   timescale=1.488e19*unitl_in_kpc*sqrt(unitl_in_kpc/unitm_in_msun)
!   lengthscale=3.086e21*unitl_in_kpc
!   velscale=lengthscale/timescale

   print*,'simple file conversion'
   print*,'filename?'
   read*,filename
!   print*,'mscale (in msun)?'
!   read*,mscale
!   mscale=mscale/unitm_in_msun
!   print*,'vscale (in km/s)?'
!   read*,vscale
!   vscale=1.e5*vscale/velscale
!   print*,'tscale (in Myr)?'
!   read*,tscale
!   tscale=tscale/timescale*year*1.e6
!   print*,'m,v,t factor:',mscale,vscale,tscale
!   lscale=1
   
   print *,'opening...  '//trim(filename)

   open (1, file=trim(filename), form='unformatted')
   read (1,iostat=ioer) np,nsph,nstar
   if(ioer.NE.0) then 
    rewind(1)
    read (1,iostat=ioer) np
    nsph=0
    nstar=0
   endif
   print*, 'nbodies:',np,nsph,nstar
   nbodies=np
!   nsph=0
!   nstar=0
   totptag=nbodies
   tnow=0.   
   massres=0.

   outputfile=trim(filename)
   firstsnap=0
   nsnap=0
   uentropy=.FALSE.
   
   output=0
   output(1:3)=1
   
   do i=1,nbodies  
    read(1) pmass,x,v
    mass(i)=pmass
    pos(i,1:3)=x
    vel(i,1:3)=v
   enddo
  
   pmass=0;x=0;v=0
   do i=1,nbodies  
    pmass=pmass+mass(i)
    x=x+mass(i)*pos(i,1:3)    
    v=v+mass(i)*vel(i,1:3)    
   enddo
   x=x/pmass
   v=v/pmass
   
   do i=1,nbodies
    pos(i,1:3)=pos(i,1:3)-x
    vel(i,1:3)=vel(i,1:3)-v
   enddo

   print*,'ethermal? (0 for T= 10^4 K)'
   read*, ether
   if(ether.LE.0) ether=1.e4/4.86e5
   output(11)=0   
   if(nsph.gt.0) then
    output(11)=1
    ethermal(1:nsph) = ether
   endif
   
   print*,'fill tform (age=12.e9 yr)'
   output(5)=1
   tform(1:nbodies)=0.
   call random_number(tform(nbodies-nstar+1:nbodies))
   tform(nbodies-nstar+1:nbodies)=-tform(nbodies-nstar+1:nbodies)*800

   print*,'provide epsgrav? (0=no, 1=yes)'
   read*,output(4)
   if(output(4).EQ.1) then
    print*,'eps_gas, eps_DM, eps_stars?'
    read*,epsg,epsdm,epsstar
    epsgrav(1:nsph)=epsg
    epsgrav(nsph+1:nbodies-nstar)=epsdm
    epsgrav(nbodies-nstar+1:nbodies)=epsstar
   endif
   outfilename=trim(filename)//'.new'
   call outbods(outfilename)
end program




