; IDL procedure to plot decay of linear wave amplitude with finite
; resistivity and/or viscosity
; waveflag=[0,1,2]=fast,Alfven,slow
;
; usage: plotdecay,'LinWave.0000.bin',10,0.02,0.01,0
;

PRO plotdecay,filename,nfiles,nu,eta,wave_flag
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz

nstart = STRPOS(filename,'.')
vymax = fltarr(nfiles)
time  = fltarr(nfiles)

; read first file, compute vymax

readbin,filename
vymax[0] = max(vy)
time[0]=0.0

; read rest of files

FOR i=1,nfiles-1 DO BEGIN
  number=STRMID(filename,nstart+1,4)
  inumber=FIX(number)+1
  IF inumber LT 10 THEN BEGIN
    number=STRING(inumber,FORMAT='(I1)')
    STRPUT,filename,number,nstart+4
  ENDIF
  IF inumber GE 10 AND inumber LT 100 THEN BEGIN
    number=STRING(inumber,FORMAT='(I2)')
    STRPUT,filename,number,nstart+3
  ENDIF
  IF inumber GE 100 AND inumber LT 1000 THEN BEGIN
    number=STRING(inumber,FORMAT='(I3)')
    STRPUT,filename,number,nstart+2
  ENDIF
  IF inumber GE 1000 AND inumber LT 10000 THEN BEGIN
    number=STRING(inumber,FORMAT='(I3)')
    STRPUT,filename,number,nstart+1
  ENDIF
  IF inumber ge 10000 THEN BEGIN
    print,'ERROR in ANIM_PLOT, Invalid Filename'
    RETURN
  ENDIF
  readbin,filename

;
; Compute maximum Vy
; Time assumes tlim=5.0
;
  vymax[i] = max(vy)
  time[i]=(i*5.0)/nfiles
ENDFOR

plot,time,(vymax/vymax[0]),YTITLE='MAX(Vy)',XTITLE='t/t_s'

; overplot analytic solution, assuming k=2\pi, Cf=2,CA=1,B^2/d=13/4,Cs=1/2

vylin = fltarr(nfiles)
vylin[0]=1.0

if wave_flag EQ 0 THEN BEGIN
  gamma = 5.264*(4.0*((7.0/3.0)*nu + eta) - 3.25*nu - nu/3.0 - (nu+eta))
ENDIF
if wave_flag EQ 1 THEN gamma = 19.739*(eta + nu)
if wave_flag EQ 2 THEN BEGIN
  gamma = 5.264*(3.25*nu + nu/3.0 + (nu+eta) - 0.25*((7.0/3.0)*nu + eta))
ENDIF

FOR i=1,nfiles-1 DO vylin[i] = exp(-gamma*time[i])

oplot,time,vylin,linestyle=1

END
