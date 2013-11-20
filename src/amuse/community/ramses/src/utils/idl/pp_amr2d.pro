;+
; NAME:
;	PP_AMR2D
;
; PURPOSE:
;	This procedure plots to the current device the mesh structure
;	stored in memory. Valid for 2D RAMSES simulations. Only octs
;	boundaries are shown.
;
; CATEGORY:
;	Plotting routines.
;
; CALLING SEQUENCE:
;       PP_AMR2D, Grid, SCALE=scale, NOERASE=noerase, NODATA=nodata,
;                       XR=xr, YR=yr, X0=x0, Y0=y0, LEVEL=level,
;                       COLOR=color 
;
; INPUTS:
;       Grid: structure defining the AMR grid.
;
; OPTIONAL INPUTS:
;       SCALE:   if set, defines the scale of the
;       plot. Scale=1. corresponds to the full box size, while
;       scale=0.5 corresponds to half the box size and so on. 
;       Default: 1.
;
;       NOERASE: standard IDL option.
;
;       NODATA:  standard IDL option.
;
;       XR:      if set, the routine plots only octs that lie within
;       these boundaries, defined along the X axis. Default: the
;       whole box ! 
;
;       YR:      Same for the Y axis.
;
;       X0:      Move the center to this new coordinate. Default: the
;       center of the box.
;       
;       Y0:      Move the center to this new coordinate. Default: the
;       center of the box.
;      
;       LEVEL:   If set, the routine plots only octs belonging to
;       level = level. Otherwise, all levels are shown (default).
;       
;       COLOR:   If set, each level is shown with a different
;       color. Otherwise, all levels are shown in white (default).
;
; OUTPUTS:
;       None.
;       
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To plot the mesh structure currently stored in memory, type:
;
;	        PP_AMR2D, Grid, /COLOR
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro plot_square,x,y,color=color,regular=regular
ss=size(x)
if(ss(0) eq 1)then begin
    if keyword_set(regular) then begin
        oplot,[x(0),x(1),x(3),x(2),x(0)],[y(0),y(1),y(3),y(2),y(0)],color=color
    endif else begin
        plots,x(0),y(0),color=color
        plots,x(1),y(1),/continue,color=color
        plots,x(3),y(3),/continue,color=color
        plots,x(2),y(2),/continue,color=color
        plots,x(0),y(0),/continue,color=color
    endelse
endif
if (ss(0) eq 2) then begin
    n=ss(1)
    for i=0L,n-1L do begin
        if keyword_set(regular) then begin
            oplot,[x(i,0),x(i,1),x(i,3),x(i,2),x(i,0)],$
              [y(i,0),y(i,1),y(i,3),y(i,2),y(i,0)],color=color
        endif else begin
            plots,[x(i,0),x(i,1),x(i,3),x(i,2),x(i,0)],$
              [y(i,0),y(i,1),y(i,3),y(i,2),y(i,0)],color=color
        endelse
    endfor
endif
end

;###################################################
;###################################################
;###################################################
pro pp_amr2d,grid,scale=scale,noerase=noerase $
           ,xr=xr,yr=yr,nodata=nodata,cpu=cpu,cmin=cmin,cmax=cmax $
           ,lmin=lmin,lmax=lmax $
           ,x0=x0,y0=y0,level=level,color=color,regular=regular,expand=expand

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'pp_amr2d'
    RETURN
ENDIF

ndim=grid.ndim
nlevelmax=grid.nlevelmax
ngrid=grid.ngrid
ncpu=grid.ncpu

if ndim lt 2 then begin
    print,'Mesh should have 2 or 3 dimensions'
    print,'but ndim=',ndim
    return
endif

if not(keyword_set(scale)) then scale = 1.0
scalein=scale/sqrt(3.0d0)
if keyword_set(color) then tek_color
if not keyword_set(x0) then x0=0.5d0
if not keyword_set(y0) then y0=0.5d0
if not keyword_set(xr)then xr=[0.0d0,1.0d0]
if not keyword_set(yr)then yr=[0.0d0,1.0d0]

xz=[x0,x0]+scalein*[-0.5,0.5]
yz=[y0,y0]+scalein*[-0.5,0.5]*!d.y_size/!d.x_size
if not keyword_set(regular) then begin
    !p.t3d=1
    scale3,xrange=xz,yrange=yz,ax=90,az=0
    if not keyword_set(noerase) then erase
endif
if not keyword_set(lmin) then lmin=1
if not keyword_set(lmax) then lmax=nlevelmax
if not keyword_set(cmin) then cmin=1
if not keyword_set(cmax) then cmax=ncpu

leveldown=lmin
levelup=lmax
if keyword_set(level)then begin
    leveldown=level
    levelup=level
endif
for ilevel=leveldown,levelup do begin
    for icpu=cmin-1,cmax-1 do begin
        if not keyword_set(color) then begin
            color1=!d.table_size-1
        endif else begin
            if keyword_set(cpu) then color1=icpu+1 else $
              color1=ilevel
        endelse
        dx=0.5d0^(ilevel-1)
        if(ngrid(ilevel-1,icpu) gt 0)then begin
            mesh2=(*grid.level[ilevel-1,icpu])
            x=mesh2.xg(*,0)
            y=mesh2.xg(*,1)
            if not keyword_set(expand) then begin
                ind=where(x gt xr[0] and x lt xr[1] and $
                          y gt yr[0] and y lt yr[1], nok)
            endif else begin
                ind=where(x gt xr[0]-dx and x lt xr[1]+dx and $
                          y gt yr[0]-dx and y lt yr[1]+dx, nok)
            endelse
            if not keyword_set(nodata) and nok gt 0 then begin
                xf=dblarr(ngrid(ilevel-1,icpu),4)
                yf=dblarr(ngrid(ilevel-1,icpu),4)
                for i=0,1 do begin
                    for j=0,1 do begin
                        ind1=i+2*j
                        xf(ind,ind1)=x(ind)+(double(i)-0.5d0)*dx
                        yf(ind,ind1)=y(ind)+(double(j)-0.5d0)*dx
                    endfor
                endfor
                plot_square,xf(ind,*),yf(ind,*),color=color1,regular=regular
            endif
        endif
    endfor
endfor

if not keyword_set(regular) then !p.t3d=0
end
;###################################################
;###################################################
;###################################################
