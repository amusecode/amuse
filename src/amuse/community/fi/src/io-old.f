C this file contains the input/output routines
C ol -style
C***********************************************************************
C
C
                          SUBROUTINE defaultparams
C
C
C=======================================================================
 
	INCLUDE 'globals.h'

        INTEGER i,kk

         datadir='./datadir/'
         inputfile='input'
	 outputfile='output'
	 firstsnap=0
         nsteps=0
         stepout=5
         steplog=5
         dtime=1.
         bh_tol=0.5
         eps=1.
         adaptive_eps=.FALSE.
         targetnn=32
         nn_tol=0.1
         usequad=.FALSE.
         directsum=.FALSE.
!         sboundary=boundary

         max_tbin=4096
         minppbin=1

         usesph=.TRUE.
         sphinit=.TRUE.
         uentropy=.TRUE.
         isotherm=.FALSE.
         sph_visc='sph '
         epsgas=0.005
         gamma=1.66667
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
         output(1:3)=1

         RETURN
        END



C***********************************************************************
C
C
                          SUBROUTINE inparams
C
C
C=======================================================================
 
	INCLUDE 'globals.h'

        CHARACTER*9 sboundary,boundary
        CHARACTER *1 pcomment
        INTEGER i,ioerror,kk
        LOGICAL bdummy

        if(.not.periodic) boundary='vacuum'
        if(periodic) boundary='periodic'

        print*,' ...reading run info...'
        OPEN(UNIT=upars,FILE=parsfile,STATUS='OLD',IOSTAT=ioerror)
        IF(ioerror.NE.0) THEN
	 PRINT*,' cannot read parameters from file: ',parsfile
         PRINT*,' using defaults'

	 CALL defaultparams
	
	 RETURN
	 	 
	ENDIF

        READ(upars,'(a)') pcomment
        READ(upars,'(a)') pcomment
        READ(upars,'(a)') datadir
	datadir=adjustl(datadir)
        i=index(datadir,' ')
        if(i.ne.0) datadir=datadir(1:i)

        READ(upars,'(a)') inputfile
	READ(upars,'(a)') outputfile
	READ(upars,*) firstsnap
        READ(upars,*) nsteps
        READ(upars,*) stepout
        READ(upars,*) steplog
        READ(upars,*) dtime
        READ(upars,*) bh_tol
        READ(upars,*) eps
        READ(upars,*) adaptive_eps
        READ(upars,*) targetnn
        READ(upars,*) nn_tol
        READ(upars,*) usequad
        READ(upars,*) directsum
        READ(upars,'(a)') sboundary
	IF(sboundary.NE.boundary) 
     &	  print*,'warning: compiled with: boundary=',boundary

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) max_tbin
        READ(upars,*) minppbin

        READ(upars,'(a)') pcomment

        READ(upars,*) usesph
        READ(upars,*) sphinit
        READ(upars,*) uentropy
        READ(upars,*) isotherm
        READ(upars,'(a)') sph_visc
        READ(upars,*) epsgas
        READ(upars,*) gamma
        READ(upars,*) alpha
        READ(upars,*) beta
        READ(upars,*) epssph
        READ(upars,*) courant
        READ(upars,*) 
        READ(upars,*) eps_is_h
        READ(upars,*) consthsm
        READ(upars,*) nsmooth
        READ(upars,*) nsmtol

        READ(upars,'(a)') symmetry
        READ(upars,*) 
        READ(upars,*)  
        READ(upars,*)  

        READ(upars,'(a)') pcomment

        READ(upars,*) radiate
        READ(upars,*) graineff
	READ(upars,*) crionrate
        READ(upars,*) heat_par1
        READ(upars,*) heat_par2
        READ(upars,*) cool_par
        READ(upars,*) 
        READ(upars,*)  ! mh/kboltz
        READ(upars,*) unitl_in_kpc
        READ(upars,*) unitm_in_msun
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*)
        READ(upars,*)
        READ(upars,*)

        READ(upars,'(a)') pcomment

        READ(upars,'(a)') halofile
	halofile=adjustl(halofile)
        i=index(halofile,' ')
        if(i.ne.0) halofile=halofile(1:i)

        READ(upars,*) fixthalo
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*)
        READ(upars,*)
        READ(upars,*)
        READ(upars,*)
        READ(upars,*)
        READ(upars,*)

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*)

        READ(upars,'(a)') pcomment

        READ(upars,*) selfgrav
        READ(upars,*) starform
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*)

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) cosmo
        READ(upars,*) comove
        READ(upars,*) 
        READ(upars,*) pboxsize
        READ(upars,*) 

        READ(upars,'(a)') pcomment
        READ(upars,*) optdepth
        READ(upars,*) tcollfac
        READ(upars,*) masscrit
        READ(upars,*) removgas
        READ(upars,*) sfeff
        READ(upars,*) 
	READ(upars,*) tbubble
        READ(upars,*) sne_eff
        READ(upars,*)
        READ(upars,*)
        READ(upars,*) tsnbeg
        READ(upars,*) feedback
        READ(upars,*) bdummy
        if(bdummy) verbosity=99
        
	READ(upars,'(a)') pcomment
        READ(upars,'(a)') hupdatemethod
        READ(upars,*)  smoothinput
        READ(upars,*) consph
	READ(upars,'(a)') sfmode
        READ(upars,'(a)') 
        READ(upars,*) rhomax
	READ(upars,*) 
	READ(upars,*) sqrttstp
	READ(upars,*) acc_tstp
	READ(upars,*) tstepcrit
	READ(upars,*) tstpcr2
	READ(upars,*) freetstp
	READ(upars,*) freev
	READ(upars,*) freea
	READ(upars,*) freevexp
	READ(upars,*) freeaexp
	READ(upars,*) gdgop
	READ(upars,*) gdgtol

        READ(upars,'(a)') pcomment

        output=0
	kk=0
	DO  i=1,64
	 IF(kk+1.EQ.i) THEN
          READ(upars,*) kk,pcomment, output(kk)
	 ENDIF
	ENDDO
		
        CLOSE(UNIT=upars)

        IF(bdummy) THEN
	 print*,'  > detailed info: loud =.TRUE.'   
         print*,'  >  (< a > t: min average max total)'
	ELSE
	 print*,'  > sparse info: loud=.FALSE.' 
        ENDIF
	
	RETURN
        END
	

C***********************************************************************
C
C
                          SUBROUTINE someparams
C
C
C=======================================================================
 
	INCLUDE 'globals.h'

        CHARACTER*9 sboundary,boundary
        CHARACTER *1 pcomment,pdummy
        INTEGER i,ioerror,kk,idummy
	REAL rdummy
        LOGICAL bdummy

        if(.not.periodic) boundary='vacuum'
        if(periodic) boundary='periodic'

        print*,' ...scanning run info...'
        OPEN(UNIT=upars,FILE=parsfile,STATUS='OLD',IOSTAT=ioerror)
        IF(ioerror.NE.0) THEN
	 PRINT*,' cannot read parameters from file: ',parsfile
	 PRINT*,' thus no changes'
	 RETURN	 	 
	ENDIF

        READ(upars,'(a)') pcomment
        READ(upars,'(a)') pcomment
        READ(upars,'(a)') datadir
        READ(upars,'(a)') pdummy !inputfile
	READ(upars,'(a)') outputfile
	READ(upars,*) idummy !firstsnap
        READ(upars,*) nsteps
        READ(upars,*) stepout
        READ(upars,*) steplog
        READ(upars,*) rdummy !dtime
        READ(upars,*) rdummy !bh_tol
        READ(upars,*) rdummy !eps
        READ(upars,*) bdummy !adaptive_eps
        READ(upars,*) idummy !targetnn
        READ(upars,*) rdummy !nn_tol
        READ(upars,*) bdummy !usequad
        READ(upars,*) bdummy !directsum
        READ(upars,'(a)') pdummy !sboundary

        READ(upars,'(a)') pcomment

        READ(upars,*)
        READ(upars,*) max_tbin
        READ(upars,*) minppbin

        CLOSE(UNIT=upars)
        PRINT*,'...So far so good...'
	RETURN
        
        READ(upars,'(a)') pcomment

        READ(upars,*) usesph
        READ(upars,*) sphinit
        READ(upars,*) uentropy
        READ(upars,*) isotherm
        READ(upars,'(a)') sph_visc
        READ(upars,*) epsgas
        READ(upars,*) gamma
        READ(upars,*) alpha
        READ(upars,*) beta
        READ(upars,*) epssph
        READ(upars,*) courant
        READ(upars,*) 
        READ(upars,*) eps_is_h
        READ(upars,*) consthsm
        READ(upars,*) nsmooth
        READ(upars,*) nsmtol

        READ(upars,'(a)') symmetry
        READ(upars,*)
        READ(upars,*) 
        READ(upars,*) 

        READ(upars,'(a)') pcomment

        READ(upars,*) radiate
        READ(upars,*) graineff
	READ(upars,*) crionrate
        READ(upars,*) heat_par1
        READ(upars,*) heat_par2
        READ(upars,*) cool_par
        READ(upars,*) 
        READ(upars,*)
        READ(upars,*) unitl_in_kpc
        READ(upars,*) unitm_in_msun
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*)
        READ(upars,*)

        READ(upars,'(a)') pcomment

        READ(upars,'(a)') halofile
	halofile=adjustl(halofile)
        i=index(halofile,' ')
        if(i.ne.0) halofile=halofile(1:i)

        READ(upars,*) fixthalo
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) 

        READ(upars,'(a)') pcomment

        READ(upars,*) selfgrav
        READ(upars,*) starform
        READ(upars,*) 
        READ(upars,*)
        READ(upars,*)

        READ(upars,'(a)') pcomment

        READ(upars,*) 
        READ(upars,*) 
        READ(upars,*) cosmo
        READ(upars,*) comove
        READ(upars,*) 
        READ(upars,*) pboxsize
        READ(upars,*)

        READ(upars,'(a)') pcomment
        READ(upars,*) optdepth
        READ(upars,*) tcollfac
        READ(upars,*) masscrit
        READ(upars,*) removgas
        READ(upars,*) sfeff
        READ(upars,*) 
	READ(upars,*) tbubble
        READ(upars,*) sne_eff
        READ(upars,*)
        READ(upars,*)
        READ(upars,*) tsnbeg
        READ(upars,*) feedback
        READ(upars,*) bdummy
        if(bdummy) verbosity=99
        
	READ(upars,'(a)') pcomment
        READ(upars,'(a)') hupdatemethod
        READ(upars,*)  smoothinput
        READ(upars,*) consph
	READ(upars,'(a)') sfmode
        READ(upars,'(a)') 
        READ(upars,*) rhomax
	READ(upars,*) 
	READ(upars,*) sqrttstp
	READ(upars,*) acc_tstp
	READ(upars,*) tstepcrit
	READ(upars,*) tstpcr2
	READ(upars,*) freetstp
	READ(upars,*) freev
	READ(upars,*) freea
	READ(upars,*) freevexp
	READ(upars,*) freeaexp
	READ(upars,*) gdgop
	READ(upars,*) gdgtol

        READ(upars,'(a)') pcomment

        output=0
	kk=0
	DO  i=1,64
	 IF(kk+1.EQ.i) THEN
          READ(upars,*) kk,pcomment, output(kk)
	 ENDIF
	ENDDO
		
        IF(bdummy) THEN
	 print*,'  > detailed info: loud =.TRUE.'   
         print*,'  >  (< a > t: min average max total)'
	ELSE
	 print*,'  > sparse info: loud=.FALSE.' 
        ENDIF
	
	RETURN
        END

