program makecube
 use MakeMapMod
 use StarsMod
 include 'globals.h'
 integer nba,omap
 character*80 filenaam
 character(len=5) :: mapmode
 call initmem(nbodsmax,nsphmax,ncells)

 interface
   subroutine WriteMap(filename,tauname)
     character, optional :: filename*80,tauname*80
   end subroutine
 end interface
 
 CALL set_parameters(0)

 print*,'filenaam?'
 read*,filenaam
 if(filenaam.EQ."debug") goto 5
 CALL readbods(filenaam)

 CALL heattabel
        
 CALL initpars

 IF(periodic) CALL initbc

 CALL initnewstar

 if(usepm) then
  if(verbosity.GT.0) print*,' ...initPM...'
  call initpm
 endif

 CALL postprocessread
 
 nba=nbandsQ()

 if(input(42).EQ.0) then
  call nbexistinit
 endif

5 continue
  
 call InitMap

 print*,' nbands=', nba
 print*,'  HI         = HI,'
 print*,'  1...nbands = band i,'
 print*,'  H2         = H_2,'
 print*,'  Ha         = Halpha,'
 print*,'  FUV        = FUV,'
 print*,'  stars      = stellar mass,'
 print*,'  gas        = gas mass,'
 print*,'  dark       = DM mass,'
 print*,'  all        = all mass,'
 print*,'  new        = new star mass,'
 print*,'  unst       = unstable gas,'
 print*,'  C          = C,'
 print*,'  CO         = CO' 
 print*,'  Ccool      = C+ cooling' 
 print*,' map?'
 read*, mapmode

 call projectparticles(mapmode)
 call WriteMap
 call EndMap

print*,'other map?(1=yes)'
read*,omap
if(omap.eq.1) goto 5
 
end program

