;+
; NAME:
;       PROF1D
;
; PURPOSE:
;       This procedure reads 1D profiles from a RAMSES LOG ASCII file
;       and plot the results in 4 windows (rho, u, P, levels).
;
; CATEGORY:
;       Plotting hydro 1D profiles.
;
; CALLING SEQUENCE:
;       PROF1D,file=file,xr=xr,symsize=symsize,itime=itime
;
; OPTIONAL INPUTS:
;       file:   if set, input the scalar string containing the name of
;       the file to be read. Otherwise, a PICKFILE widget is launched. 
;
;       xr:     the X axis range. Default: [0,1]
;
;       symsize: the size of the square symbol to plot the numerical
;       solution. Default=1.
;
;       itime:  if set, plot profiles for output time number
;       itime. Otherwise, plot all output times.
;
; OUTPUTS:
;       None.
;;
; EXAMPLE:
;       To plot profiles from a RAMSES LOG ASCII file, type:
;
;               PROF1D,file='hydro.log',itime=2
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro prof1d,file=file,xr=xr,dr=dr,ur=ur,pr=pr,lr=lr $
           ,symsize=symsize,step=step,itime=itime,noerase=noerase

if not keyword_set(file) then file=pickfile(/READ)
if not keyword_set(symsize) then symsize=1.
if not keyword_set(step) then step=1

ss=symsize
tek_color
rd_1d,output,file=file
ntime=n_elements(output)
if keyword_set(itime) then begin
    itimemax=MIN([itime,ntime])
    itimemax=MAX([itime,1])
    itimemin=itimemax
endif else begin
    itimemin=1
    itimemax=ntime
endelse
for i=itimemin,itimemax,step do begin
    ttt=(*output[i-1])
    xmin=min(ttt.x) & xmax=max(ttt.x)
    dmin=min(ttt.d) & dmax=max(ttt.d)
    umin=min(ttt.u) & umax=max(ttt.u)
    pmin=min(ttt.p) & pmax=max(ttt.p)
    lmin=min(ttt.l) & lmax=max(ttt.l)
    if(i eq itimemin) then begin
        xmin0=xmin & xmax0=xmax
        dmin0=dmin & dmax0=dmax
        umin0=umin & umax0=umax
        pmin0=pmin & pmax0=pmax
        lmin0=lmin & lmax0=lmax
    endif else begin
        xmin0=min([xmin0,xmin]) & xmax=max([xmax0,xmax])
        dmin0=min([dmin0,dmin]) & dmax=max([dmax0,dmax])
        umin0=min([umin0,umin]) & umax=max([umax0,umax])
        pmin0=min([pmin0,pmin]) & pmax=max([pmax0,pmax])
        lmin0=min([lmin0,lmin]) & lmax=max([lmax0,lmax])
    endelse
endfor
if not keyword_set(xr) then xr=[xmin0,xmax0]
if not keyword_set(dr) then dr=[dmin0,dmax0]
if not keyword_set(pr) then pr=[pmin0,pmax0]
if not keyword_set(ur) then ur=[umin0,umax0]
if not keyword_set(lr) then lr=[lmin0,lmax0]

!p.multi=[0,2,2]
plot,(*output[itimemax-1]).x,(*output[itimemax-1]).u $
  ,color=1,xr=xr,yr=ur,/xs,/ys,symsize=ss,xtitle='!17x',ytitle='u' $
  ,noerase=noerase
for i=itimemin,itimemax,step do begin
    oplot,(*output[i-1]).x,(*output[i-1]).u,color=1,symsize=ss
endfor

!p.multi=[3,2,2]
plot,(*output[itimemax-1]).x,(*output[itimemax-1]).d $
  ,color=1,xr=xr,yr=dr,/xs,/ys,symsize=ss,xtitle='x',ytitle='!7q!17' $
  ,noerase=noerase
for i=itimemin,itimemax,step do begin
    oplot,(*output[i-1]).x,(*output[i-1]).d,color=1,symsize=ss
endfor

!p.multi=[2,2,2]
plot,(*output[itimemax-1]).x,(*output[itimemax-1]).p $
  ,color=1,xr=xr,yr=pr,/xs,/ys,symsize=ss,xtitle='x',ytitle='P' $
  ,noerase=noerase
for i=itimemin,itimemax,step do begin
    oplot,(*output[i-1]).x,(*output[i-1]).p,color=1,symsize=ss
endfor

!p.multi=[1,2,2]
plot,(*output[itimemax-1]).x,(*output[itimemax-1]).l $
  ,xr=xr,yr=lr,/xs,/ys,psym=10,color=1,symsize=ss,xtitle='x',ytitle='Level' $
  ,noerase=noerase
for i=itimemin,itimemax,step do begin
    oplot,(*output[i-1]).x,(*output[i-1]).l,color=1,psym=10,symsize=ss
endfor

end
;###################################################
;###################################################
;###################################################
