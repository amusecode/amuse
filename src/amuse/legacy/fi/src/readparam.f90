subroutine set_parameters_to_defaults
  include 'globals.h'

  datadir='./datadir/'
  inputfile='input'
  outputfile='output'
  firstsnap=0
  nsteps=0
  stepout=5
  steplog=5
  dtime=1.
  bh_tol=0.5
  eps=0.
  adaptive_eps=.FALSE.
  targetnn=32
  nn_tol=0.1
  usequad=.FALSE.
  directsum=.FALSE.
  max_tbin=4096
  minppbin=1
  usesph=.TRUE.
  sphinit=.TRUE.
  uentropy=.TRUE.
  isotherm=.FALSE.
  sph_visc='sph '
  epsgas=0.005
  gamma=1.6666667
  alpha=0.5
  beta=1.
  epssph=0.01
  courant=0.3
  eps_is_h=.TRUE.
  consthsm=0.2
  nsmooth=64
  nsmtol=0.1
  symmetry='hk'
  radiate=.FALSE.
  graineff=0.05
  crionrate=3.6
  heat_par1=0.0
  heat_par2=0.0
  cool_par=1.
  unitl_in_kpc=1.0
  unitm_in_msun=1.e9
  halofile='none'
  fixthalo=.FALSE.
  selfgrav=.TRUE.
  starform=.FALSE.
  cosmo=.FALSE.
  comove=.FALSE.
  pboxsize=300.
  optdepth=0.
  tcollfac=1.
  masscrit=1.e5
  removgas=0.25
  sfeff=0.25
  tbubble=3.e7
  sne_eff=0.
  tsnbeg=3.e6
  feedback='fuv '
  verbosity=0
  hupdatemethod='mass'
  smoothinput=.FALSE.
  consph=.TRUE.
  sfmode='gerritsen'
  rhomax=100
  sqrttstp=.FALSE.
  acc_tstp=.TRUE.
  tstepcrit=1.
  tstpcr2=0.25
  freetstp=.FALSE.
  freev=0.5
  freea=.35
  freevexp=0.
  freeaexp=-1.
  gdgop=.TRUE.
  gdgtol=0.01

  output=0

end subroutine

subroutine write_parameters(restart)
  include 'globals.h'
  integer :: restart
  integer :: out(64),i,no,nout=64
  character(len=256) :: outfile=""
  character(len=256) :: outline=""
  character(len=20) :: keys(64)
  data keys/ "mass","pos","vel","eps","tform","acc", &
"","","","rho","ethermal","entropy","hsmooth","fuvheat","esnthdt", &
"tcollaps","temperat","elecfrac","csound","pressure","hsmdivv","mumaxdvh", &
"hsmcurlv","vdisp","h2frac","","","","","","","","","starfuv","snentropy", &
"itimestp","","","","phi","phiext","nbexist","","","","","","","","","","", &
"","","","","","","","","","","","" /

  if(restart.EQ.0) then
    outfile='parameters_used.txt'
  else
    outfile='restart_parameters_used.txt'
  endif    

  open(unit=upars,file=TRIM(outfile), status='unknown')
  
  write(upars, '(" datadir:   ",a)')  '"'//TRIM(datadir)//'"'
  write(upars, '(" inputfile:   ",a)')  '"'//TRIM(inputfile)//'"'
  write(upars, '(" outputfile:   ",a)')  '"'//TRIM(outputfile)//'"'
  write(upars,*) "firstsnap:   ",  firstsnap
  write(upars,*) "nsteps:   ",  nsteps
  write(upars,*) "stepout:   ",  stepout
  write(upars,*) "steplog:   ",  steplog
  write(upars,*) "dtime:   ",  dtime
  write(upars,*) "bh_tol:   ",  bh_tol
  write(upars,*) "eps:   ",  eps
  write(upars,*) "adaptive_eps:   ",  adaptive_eps
  write(upars,*) "targetnn:   ",  targetnn
  write(upars,*) "nn_tol:   ",  nn_tol
  write(upars,*) "usequad:   ",  usequad
  write(upars,*) "directsum:   ",  directsum
  write(upars,*) "periodic:   ",  periodic
  write(upars,*) "max_tbin:   ",  max_tbin
  write(upars,*) "minppbin:   ",  minppbin
  write(upars,*) "usesph:   ",  usesph
  write(upars,*) "sphinit:   ",  sphinit
  write(upars,*) "uentropy:   ",  uentropy
  write(upars,*) "isotherm:   ",  isotherm
  write(upars, '(" sph_visc:   ",a)')  '"'//TRIM(sph_visc)//'"'
  write(upars,*) "epsgas:   ",  epsgas
  write(upars,*) "gamma:   ",  gamma
  write(upars,*) "alpha:   ",  alpha
  write(upars,*) "beta:   ",  beta
  write(upars,*) "epssph:   ",  epssph
  write(upars,*) "courant:   ",  courant
  write(upars,*) "eps_is_h:   ",  eps_is_h 
  write(upars,*) "consthsm:   ",  consthsm
  write(upars,*) "nsmooth:   ",  nsmooth
  write(upars,*) "nsmtol:   ",  nsmtol
  write(upars, '(" symmetry:   ",a)')  '"'//TRIM(symmetry)//'"'
  write(upars,*) "radiate:   ",  radiate
  write(upars,*) "graineff:   ",  graineff
  write(upars,*) "crionrate:   ",  crionrate
  write(upars,*) "heat_par1:   ",  heat_par1
  write(upars,*) "heat_par2:   ",  heat_par2
  write(upars,*) "cool_par:   ",  cool_par
  write(upars,*) "unitl_in_kpc:   ",  unitl_in_kpc
  write(upars,*) "unitm_in_msun:   ",  unitm_in_msun
  write(upars, '(" halofile:   ",a)')  '"'//TRIM(halofile)//'"'
  write(upars,*) "fixthalo:   ",  fixthalo
  write(upars,*) "selfgrav:   ",  selfgrav
  write(upars,*) "starform:   ",  starform
  write(upars,*) "cosmo:   ",  cosmo
  write(upars,*) "comove:   ",  comove
  write(upars,*) "pboxsize:   ",  pboxsize
  write(upars,*) "optdepth:   ",  optdepth
  write(upars,*) "tcollfac:   ",  tcollfac
  write(upars,*) "masscrit:   ",  masscrit
  write(upars,*) "removgas:   ",  removgas
  write(upars,*) "sfeff:   ",  sfeff
  write(upars,*) "tbubble:   ",  tbubble
  write(upars,*) "sne_eff:   ",  sne_eff
  write(upars,*) "tsnbeg:   ",  tsnbeg
  write(upars, '(" feedback:   ",a)')  '"'//TRIM(feedback)//'"'
  write(upars,*) "verbosity:   ",  verbosity
  write(upars, '(" hupdatemethod:   ",a)')  '"'//TRIM(hupdatemethod)//'"'
  write(upars,*) "smoothinput:   ",  smoothinput
  write(upars,*) "consph:   ",  consph
  write(upars, '(" sfmode:   ",a)')  '"'//TRIM(sfmode)//'"'
  write(upars,*) "rhomax:   ",  rhomax
  write(upars,*) "sqrttstp:   ",  sqrttstp
  write(upars,*) "acc_tstp:   ",  acc_tstp
  write(upars,*) "tstepcrit:   ",  tstepcrit
  write(upars,*) "tstpcr2:   ",  tstpcr2
  write(upars,*) "freetstp:   ",  freetstp
  write(upars,*) "freev:   ",  freev
  write(upars,*) "freea:   ",  freea
  write(upars,*) "freevexp:   ",  freevexp
  write(upars,*) "freeaexp:   ",  freeaexp
  write(upars,*) "gdgop:   ",  gdgop
  write(upars,*) "gdgtol:   ",  gdgtol

  no=0
  outline=''
  do i=1,nout
    if(output(i).NE.0.AND.keys(i).NE."") then
      if(no.NE.0) outline=TRIM(outline)//','
      outline=TRIM(outline)//TRIM(keys(i))
      no=1
    endif  
  enddo
  write(upars,'("output:   ",a)') TRIM(outline)
  no=0
  outline=''
  do i=1,nout
    if(output(i).EQ.0.AND.keys(i).NE."") then
      if(no.NE.0) outline=TRIM(outline)//','
      outline=TRIM(outline)//TRIM(keys(i))
      no=1
    endif  
  enddo
  write(upars,'("non-output:   ",a)') TRIM(outline)
  
  close(upars)

end subroutine

subroutine get_output(line,unit,output,set)
  integer :: output(64),nout=64
  character(len=256) line
  character(len=20) key
  character(len=256) :: outline=""
  character(len=20) :: keys(64)
  data keys/"mass","pos","vel","eps","tform","acc",&
  "","","","rho","ethermal","entropy","hsmooth","fuvheat","esnthdt",      &
  "tcollaps","temperat","elecfrac","csound","pressure","hsmdivv","mumaxdvh",&
  "hsmcurlv","vdisp","h2frac","","","","","","","","","starfuv","snentropy",&
  "itimestp","","","","phi","phiext","nbexist","","","","","","","","","","", &
  "","","","","","","","","","","","" /
  integer :: i,j,unit,set
  i=0
  do
    line=line(i+1:)
    i=index(line,',')
    if(i.EQ.0) then 
      key=line
    else
      if(i.EQ.1) cycle
      key=line(1:i-1)
    endif
    key=adjustl(key) 
    if(key.EQ.'') then
      read(unit,'(a)',iostat=i) line
      line=adjustl(line)
      cycle
    endif
    j=0
    do while(j.LE.nout)
      if(keys(j).EQ.key) then
        output(j)=set
        exit
      endif  
      j=j+1
    enddo
    if(j.EQ.nout+1) then
       print*,'not found:', key 
       stop
    endif   
    if(i.eq.0) exit 
  enddo
end subroutine

subroutine read_parameters(restart)
  include 'globals.h'
  integer :: restart
  character(len=20) key, value
  character(len=256) line
  character(len=20) varc1,varc2
  character :: commentchar(6)=(/'!','#','$','%',';',':'/)
  integer ierr,i,isok
  
  print*,' ...reading run info...'
  open(unit=upars,file=parsfile,status='OLD',iostat=ierr)
  if(ierr.NE.0.AND.restart.EQ.0) then
    print*,' cannot read parameters from file: ',parsfile
    print*,' using defaults'
    return
  endif

  isok=0
  do 
    read(upars,'(a)',iostat=ierr) line
    if(ierr.NE.0) exit
    line=adjustl(line)
    i=index(line,':')
    if(any(line(1:1).EQ.commentchar).OR.line.EQ."") cycle
    if(i.LE.1) then
      print*,"in file ",TRIM(parsfile),", I don't understand:"
      print*, TRIM(line)
      print*,"(comments start with: %)"
      isok=1
      cycle
    endif
    key=line(1:i-1)
    line=line(i+1:)
    line=adjustl(line)
    i=index(line,' ')

    select case(key)
    
      case ("datadir")
        if(restart.EQ.0) read(line,*,iostat=ierr) datadir        
      case ("inputfile")
        if(restart.EQ.0) read(line,*,iostat=ierr) inputfile       
      case ("outputfile")
        read(line,*,iostat=ierr) outputfile       
      case ("firstsnap")
        if(restart.EQ.0) read(line,*,iostat=ierr) firstsnap      
      case ("nsteps")
        read(line,*,iostat=ierr) nsteps         
      case ("stepout")
        read(line,*,iostat=ierr) stepout        
      case ("steplog")
        read(line,*,iostat=ierr) steplog        
      case ("dtime")
        if(restart.EQ.0) read(line,*,iostat=ierr) dtime          
      case ("bh_tol")
        if(restart.EQ.0) read(line,*,iostat=ierr) bh_tol            
      case ("eps")
        if(restart.EQ.0) read(line,*,iostat=ierr) eps            
      case ("adaptive_eps")
        if(restart.EQ.0) read(line,*,iostat=ierr) adaptive_eps       
      case ("targetnn")
        if(restart.EQ.0) read(line,*,iostat=ierr) targetnn       
      case ("nn_tol")
        if(restart.EQ.0) read(line,*,iostat=ierr) nn_tol         
      case ("usequad")
        if(restart.EQ.0) read(line,*,iostat=ierr) usequad        
      case ("directsum")
        if(restart.EQ.0) read(line,*,iostat=ierr) directsum       
      case ("periodic")
      case ("max_tbin")
        read(line,*,iostat=ierr) max_tbin       
      case ("minppbin")
        read(line,*,iostat=ierr) minppbin       
      case ("usesph")
        if(restart.EQ.0) read(line,*,iostat=ierr) usesph         
      case ("sphinit")
        if(restart.EQ.0) read(line,*,iostat=ierr) sphinit        
      case ("uentropy")
        if(restart.EQ.0) read(line,*,iostat=ierr) uentropy       
      case ("isotherm")
        if(restart.EQ.0) read(line,*,iostat=ierr) isotherm       
      case ("sph_visc")
        if(restart.EQ.0) read(line,*,iostat=ierr) sph_visc       
      case ("epsgas")
        if(restart.EQ.0) read(line,*,iostat=ierr) epsgas         
      case ("gamma")
        if(restart.EQ.0) read(line,*,iostat=ierr) gamma          
      case ("alpha")
        if(restart.EQ.0) read(line,*,iostat=ierr) alpha          
      case ("beta")
        if(restart.EQ.0) read(line,*,iostat=ierr) beta           
      case ("epssph")
        if(restart.EQ.0) read(line,*,iostat=ierr) epssph         
      case ("courant")
        if(restart.EQ.0) read(line,*,iostat=ierr) courant        
      case ("eps_is_h")
        if(restart.EQ.0) read(line,*,iostat=ierr) eps_is_h       
      case ("consthsm")
        if(restart.EQ.0) read(line,*,iostat=ierr) consthsm       
      case ("nsmooth")
        if(restart.EQ.0) read(line,*,iostat=ierr) nsmooth        
      case ("nsmtol")
        if(restart.EQ.0) read(line,*,iostat=ierr) nsmtol         
      case ("symmetry")
        if(restart.EQ.0) read(line,*,iostat=ierr) symmetry       
      case ("radiate")
        if(restart.EQ.0) read(line,*,iostat=ierr) radiate        
      case ("graineff")
        if(restart.EQ.0) read(line,*,iostat=ierr) graineff       
      case ("crionrate")
        if(restart.EQ.0) read(line,*,iostat=ierr) crionrate      
      case ("heat_par1")
        if(restart.EQ.0) read(line,*,iostat=ierr) heat_par1        
      case ("heat_par2")
        if(restart.EQ.0) read(line,*,iostat=ierr) heat_par2       
      case ("cool_par")
        if(restart.EQ.0) read(line,*,iostat=ierr) cool_par       
      case ("unitl_in_kpc")
        if(restart.EQ.0) read(line,*,iostat=ierr) unitl_in_kpc        
      case ("unitm_in_msun")
        if(restart.EQ.0) read(line,*,iostat=ierr) unitm_in_msun       
      case ("halofile")
        if(restart.EQ.0) read(line,*,iostat=ierr) halofile       
      case ("fixthalo")
        if(restart.EQ.0) read(line,*,iostat=ierr) fixthalo       
      case ("selfgrav")
        if(restart.EQ.0) read(line,*,iostat=ierr) selfgrav       
      case ("starform")
        if(restart.EQ.0) read(line,*,iostat=ierr) starform       
      case ("cosmo")
        if(restart.EQ.0) read(line,*,iostat=ierr) cosmo          
      case ("comove")
        if(restart.EQ.0) read(line,*,iostat=ierr) comove         
      case ("pboxsize")
        if(restart.EQ.0) read(line,*,iostat=ierr) pboxsize       
      case ("optdepth")
        if(restart.EQ.0) read(line,*,iostat=ierr) optdepth       
      case ("tcollfac")
        if(restart.EQ.0) read(line,*,iostat=ierr) tcollfac       
      case ("masscrit")
        if(restart.EQ.0) read(line,*,iostat=ierr) masscrit       
      case ("removgas")
        if(restart.EQ.0) read(line,*,iostat=ierr) removgas       
      case ("sfeff")
        if(restart.EQ.0) read(line,*,iostat=ierr) sfeff          
      case ("tbubble")
        if(restart.EQ.0) read(line,*,iostat=ierr) tbubble        
      case ("sne_eff")
        if(restart.EQ.0) read(line,*,iostat=ierr) sne_eff        
      case ("tsnbeg")
        if(restart.EQ.0) read(line,*,iostat=ierr) tsnbeg         
      case ("feedback")
        if(restart.EQ.0) read(line,*,iostat=ierr) feedback       
      case ("verbosity")
        if(restart.EQ.0) read(line,*,iostat=ierr) verbosity           
      case ("hupdatemethod")
        if(restart.EQ.0) read(line,*,iostat=ierr) hupdatemethod  
      case ("smoothinput")
        if(restart.EQ.0) read(line,*,iostat=ierr) smoothinput    
      case ("consph")
        if(restart.EQ.0) read(line,*,iostat=ierr) consph         
      case ("sfmode")
        if(restart.EQ.0) read(line,*,iostat=ierr) sfmode         
      case ("rhomax")
        if(restart.EQ.0) read(line,*,iostat=ierr) rhomax         
      case ("sqrttstp")
        if(restart.EQ.0) read(line,*,iostat=ierr) sqrttstp        
      case ("acc_tstp")
        if(restart.EQ.0) read(line,*,iostat=ierr) acc_tstp        
      case ("tstepcrit")
        if(restart.EQ.0) read(line,*,iostat=ierr) tstepcrit      
      case ("tstpcr2")
        if(restart.EQ.0) read(line,*,iostat=ierr) tstpcr2        
      case ("freetstp")
        if(restart.EQ.0) read(line,*,iostat=ierr) freetstp       
      case ("freev")
        if(restart.EQ.0) read(line,*,iostat=ierr) freev          
      case ("freea")
        if(restart.EQ.0) read(line,*,iostat=ierr) freea          
      case ("freevexp")
        if(restart.EQ.0) read(line,*,iostat=ierr) freevexp       
      case ("freeaexp")
        if(restart.EQ.0) read(line,*,iostat=ierr) freeaexp       
      case ("gdgop")
        if(restart.EQ.0) read(line,*,iostat=ierr) gdgop          
      case ("gdgtol")
        if(restart.EQ.0) read(line,*,iostat=ierr) gdgtol         
      case ("output")
        if(restart.EQ.0) call get_output(line,upars,output,1)
      case ("non-output")
        if(restart.EQ.0) call get_output(line,upars,output,0)
      case default
        print*,'unknown key: ', key
        isok=1
        cycle
    end select
    if(ierr.NE.0) then
      print*,'invalid value for: ', key
      isok=1
      cycle
    endif  
  enddo
   
  close(upars)
  
  if(isok.NE.0) then
    call terror("parameterfile error")
  endif
end subroutine

subroutine check_parameters
  include 'globals.h'
  integer nerror,i
  
  if(firstsnap.LT.0.OR.firstsnap.GT.999999) call terror(' parameter error: firstsnap ')
  if(nsteps.LT.0.OR.nsteps.GT.10000000) call terror(' parameter error: nsteps ')
  if(stepout.LT.0) call terror(' parameter error: stepout ')
  if(steplog.LT.0) call terror(' parameter error: steplog ')
  if(dtime.LE.0.0.OR.dtime.GT.1.e20) call terror(' parameter error: dtime ')
  if((.NOT.directsum).AND.(bh_tol.LT.0.0.OR.bh_tol.GT.1.5)) call terror(' parameter error: bh_tol ')
  if(eps.LT.0.0.OR.eps.GT.1.e20) call terror(' parameter error: eps ')
  if(adaptive_eps.AND.(targetnn.LE.0.OR.targetnn.GT.nbodsmax)) call terror(' parameter error: targetnn ')
  if(max_tbin.LT.1.OR.max_tbin.GT.2**20) call terror(' parameter error: max_tbin ')
  if(minppbin.LE.0.OR.minppbin.GT.nbodsmax) call terror(' parameter error: minppbin ')
  if(usesph) THEN
    if(sph_visc.NE.'bulk'.AND.sph_visc.NE.'sphv'.AND. &
       sph_visc.NE.'sph ') call terror(' parameter error: sph_visc ')
    if(sph_visc.EQ.'bulk') call terror(' check bulk viscosity!')
    if(epsgas.LE.0.0) call terror(' parameter error: epsgas ')
    if(gamma.LE.0.0.OR.gamma.GT.20.) call terror(' parameter error: gamma ')
    if(alpha.LT.0.0.OR.alpha.GT.10) call terror(' parameter error: alpha ')
    if(beta.LT.0.0.OR.beta.GT.10) call terror(' parameter error: beta ')
    if(epssph.LE.0.) call terror(' parameter error: epssph ')
    if(courant.LE.0.0.OR.courant.GT.1.) call terror(' parameter error: courant ')
    if(consthsm.LT.0.AND.hupdatemethod.EQ.'none') call terror(' parameter error: consthsm ')
    if(nsmooth.LE.0.OR.nsmooth.GT.nbodsmax.AND.hupdatemethod.EQ.'none') call terror(' parameter error: nsmooth ')
    if(symmetry.NE.'hk'.AND.symmetry.NE.'be') call terror(' parameter error: symmetry ')
    if(graineff.LT.0.0.OR.graineff.GT.1) call terror(' parameter error: graineff ')
    if(crionrate.LT.0.0) call terror(' parameter error: crionrate ')
    if(heat_par1.LT.0.0) call terror(' parameter error: heat_par1 ')
    if(heat_par2.LT.0.0) call terror(' parameter error: heat_par2 ')
    if(cool_par.LT.0.0) call terror(' parameter error: cool_par ')
    if(isotherm.AND.ABS(gamma-1.0).GT.1.e-6) call terror(' gamma must be 1.0 for isothermal gas')
    if(uentropy.AND.isotherm) call terror(' uentropy must be .FALSE. for isothermal gas')
  endif

  if(pboxsize.LE.0.0) call terror(' parameter error: pboxsize ')
  if(cosmo.OR.comove) call terror(' cosmo/comove not implemented yet')
  if((.NOT.cosmo).AND.comove) call terror(' comove allowed only if cosmo is true ')
  if(usequad.AND.periodic) call terror(' quad. option allowed only with vacuum b.c. ')
  if(usequad.AND.usepm) call terror(' quad. option not allowed for PM gravity ')
  if(starform.AND.(.NOT.radiate.OR.isotherm)) call terror(' radiate=.TRUE., isotherm=.FALSE. when starform=.TRUE.')
  if(.NOT.isotherm.AND.gamma.EQ.1.0) call terror(' gamma may not be 1.0 for non isothermal')
  if(starform.OR.radiate) then
    if(optdepth.LT.0) call terror('parameter error: optdepth ')
    if(tcollfac.LT.0) call terror('parameter error: tcollfac ')
    if(masscrit.LT.0) call terror('parameter error: masscrit ')
    if(removgas.LT.0.0.OR.removgas.GT.1.0) call terror('parameter error: removgas ')
    if(sfeff.LE.0.0.OR.sfeff.GT.1.0) call terror(' parameter error: sfeff ')
    if(tbubble.LT.0.0) call terror(' parameter error: tbubble ')
    if(sne_eff.GT.100.0.OR.sne_eff.LT.0.0) call terror(' parameter error: sne_eff ')
    if(sfmode.NE.'nieuw'.AND.sfmode.NE.'gerritsen'.AND. &
       sfmode.NE.'molecular'.AND.sfmode.NE.'molecular2') call terror(' parameter error: sfmode ')	      
    if(tsnbeg.LT.0.0) call terror(' parameter error: tsnbeg ')
    if(feedback.NE.'fuv '.AND.feedback.NE.'pres'.AND.feedback.NE.'kine'.AND. &
       feedback.NE.'solo'.AND.feedback.NE.'solh') call terror(' parameter error: feedback ')
  endif
        
end subroutine



subroutine set_parameters(restart)
  include 'globals.h'
  integer restart,ierr,ierr2
  character(len=256) line


  open(unit=upars,file=parsfile,status='OLD',iostat=ierr)
  read(upars,*,iostat=ierr2) line
  close(upars)
  
  if(line(1:16).EQ."C**********Basic".OR. &
    ierr.NE.0) then
    print*,"old style input"
    if(restart.EQ.0) call inparams
    if(restart.EQ.1) call someparams
  else
    if(restart.EQ.0)  call set_parameters_to_defaults
    call read_parameters(restart)
  endif
  
  call write_parameters(restart)
  call check_parameters

end subroutine

