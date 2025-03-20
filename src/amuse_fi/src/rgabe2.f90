program read_snapshot
 use sphray_io_mod
   include 'globals.h'
   integer :: NFILES        
   type(sphray_snap_type) :: s

   character*200 fileroot,filename
   character*6 nostring

   real :: codemass,codevel,mscale,vscale,tscale,lscale
   integer no,k,i
   call initmem(nbodsmax,nsphmax,ncells)

   unitm_in_msun=1.e9
   unitl_in_kpc=1.
   timescale=1.488e19*unitl_in_kpc*sqrt(unitl_in_kpc/unitm_in_msun)
   lengthscale=3.086e21*unitl_in_kpc
   velscale=lengthscale/timescale

   print*,'gabe file conversion'
   print*,'filename?'
   read*,fileroot
   
   filename=trim(fileroot)

   k=index(filename,'.')
   if(k.le.5.or.k.ge.200) then 
    print*,'filename not of expected xxx_nnn.y form'
    stop
   endif
   print*,fileroot(k-3:k-1)

   read(fileroot(k-2:k),*) no
   print *,'opening...  '//trim(filename)
   print*,' detected No.', no


   call read_sphray_snapshot(s,filename)

   print *,'Done with reading.'

   codemass=s%CGSlum * s%CGStime**3/ s%CGSlen**2
   codevel=s%CGSlen/s%CGStime

   print*,'derived codemass, codevel:',codemass,codevel

   if(s%Nfiles.NE.1) then
    print*,'numfiles>1 not yet..'
    stop
   endif 
   
    print*,' assumed that all is gas'
   if(s%NparSnap.GT.nsphmax) then
    print*,'max number of particle exceed'
    stop
   endif 
   nsph=s%NparSnap
   nbodies=nsph

   if(s%Datapresent(1))  nbexist(1:nsph)=s%id(1:nsph)
   if(s%Datapresent(2))  pos(1:nsph,1)=s%pos(1:nsph,1)
   if(s%Datapresent(2))  pos(1:nsph,2)=s%pos(1:nsph,2)
   if(s%Datapresent(2))  pos(1:nsph,3)=s%pos(1:nsph,3)
   if(s%Datapresent(5))  vel(1:nsph,1)=s%vel(1:nsph,1)
   if(s%Datapresent(5))  vel(1:nsph,2)=s%vel(1:nsph,2)
   if(s%Datapresent(5))  vel(1:nsph,3)=s%vel(1:nsph,3)
   if(s%Datapresent(8))  hsmooth(1:nsph)=s%hsml(1:nsph) /2
   if(s%Datapresent(9))  rho(1:nsph)=s%rho(1:nsph)
   if(s%Datapresent(10))  mass(1:nsph)=s%mass(1:nsph)
   if(s%Datapresent(11))  temperat(1:nsph)=s%T(1:nsph)
   if(s%Datapresent(12))  elecfrac(1:nsph)=s%xHII(1:nsph)
   if(s%Datapresent(13))  elecfrac(1:nsph)=elecfrac(1:nsph)+0.1*s%xHeII(1:nsph)
   if(s%Datapresent(14))  elecfrac(1:nsph)=elecfrac(1:nsph)+0.1*2*s%xHeIII(1:nsph)
      
   mass(1:nsph)=mass(1:nsph)/unitm_in_msun/solarmass*codemass
   pos(1:nsph,1:3)=pos(1:nsph,1:3)/unitl_in_kpc*s%CGSlen/kpc
   hsmooth(1:nsph)=hsmooth(1:nsph)/unitl_in_kpc*s%CGSlen/kpc
   vel(1:nsph,1:3)=vel(1:nsph,1:3)/velscale*codevel
   
   totptag=nbodies
   tnow=0
!   tnow=s%Time/timescale*s%CGStime
   massres=0.

   firstsnap=0
   uentropy=.FALSE.
   output=0
   output(1:3)=1
   output(13)=1
   output(17:18)=1
   
   write(nostring,'(I6.6)') no   
   call outbods(trim(fileroot(1:k-5)//'.'//nostring))
end program




