;+
; NAME:
;       SEDOV1D
;
; PURPOSE:
;       This procedure reads 1D profiles from a RAMSES LOG ASCII file
;       for the Sedov planar blast wave test and plot the results versus
;       the analytical solution.
;
; CATEGORY:
;       Plotting hydro test case.
;
; CALLING SEQUENCE:
;       SEDOV1D,e0=e0,d0=d0,file=file,xr=xr,symsize=symsize
;
; OPTIONAL INPUTS:
;       e0:     the total energy in the half-plane (x>0). Default: 0.5
;
;       d0:     the initial mass density. Default: 1.
;
;       file:   if set, input the scalar string containing the name of
;       the file to be read. Otherwise, a PICKFILE widget is launched. 
;
;       xr:     the X axis range. Default: [0,1]
;
;       symsize: the size of the square symbol to plot the numerical
;       solution. Default=1.
;
; OUTPUTS:
;       None.
;;
; EXAMPLE:
;       To plot the Sedov planar blast wave profiles from a RAMSES LOG
;       ASCII file, type: 
;
;               SEDOV1D,file='sedov1d.log'
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro sedov1d,e0=e0,d0=d0,file=file,xr=xr,symsize=symsize

if not keyword_set(file)    then file=pickfile(/READ)
if not keyword_set(symsize) then symsize=1.
if not keyword_set(xr)      then xr=[0,0.5]
if not keyword_set(d0)      then d0=1.
if not keyword_set(e0)      then e0=0.5

yrl=[3,10]
yrd=[0.,6.]
yru=[0.,3.]
yrp=[0.,7.]
ss=symsize
tek_color

; Read the numerical solution
rd_1d,output,file=file
ntime=n_elements(output)

; Compute the analytical solution
sedovana,ra,da,ua,pa,gamma=1.4,dim=1.

scalera=dblarr(ntime)
scaleua=dblarr(ntime)
scaleda=dblarr(ntime)
scalepa=dblarr(ntime)
for i=0,ntime-1 do begin
    t=(*output[i]).t
    scalera(i)=(e0/d0)^(1./3.)*t^(2./3.)
    scaleda(i)=d0
    scalepa(i)=(e0/d0)^(2./3.)*t^(-2./3.)
    scaleua(i)=d0*(e0/d0)^(1./3.)*t^(-1./3.)
endfor

!p.multi=[0,2,2]
plot,[0.,0.],color=1,xr=xr,yr=yru,/ys,/xs $
  ,symsize=ss,xtitle='!17x',ytitle='u',/nodata,xtickformat='(F3.1)'
for i=0,ntime-1 do begin
    oplot,(*output[i]).x,(*output[i]).u,psym=6,color=1,symsize=ss
    oplot,ra*scalera(i),ua*scaleua(i),color=2
endfor

!p.multi=[1,2,2]
plot,[0.,0.],color=1,xr=xr,yr=yrd,/ys,/xs $
  ,symsize=ss,xtitle='x',ytitle='!7q!17',/nodata,xtickformat='(F3.1)'
for i=0,ntime-1 do begin
    oplot,(*output[i]).x,(*output[i]).d,psym=6,color=1,symsize=ss
    oplot,ra*scalera(i),da*scaleda(i),color=2
endfor

!p.multi=[2,2,2]
plot,[0.,0.],psym=6,color=1,xr=xr,yr=yrp,/ys,/xs $
  ,symsize=ss,xtitle='x',ytitle='P',/nodata,xtickformat='(F3.1)'
for i=0,ntime-1 do begin
    oplot,(*output[i]).x,(*output[i]).p,psym=6,color=1,symsize=ss
    oplot,ra*scalera(i),pa*scalepa(i),color=2
endfor

!p.multi=[3,2,2]
plot,[0.,0.],xr=xr,yr=yrl,/ys,/xs,color=1 $
  ,xtitle='x',ytitle='Level',/nodata,xtickformat='(F3.1)'
for i=0,ntime-1 do begin
    oplot,(*output[i]).x,(*output[i]).l,psym=10,color=1,symsize=ss
endfor

end
;###################################################
;###################################################
;###################################################
