;+
; NAME:
;       TV2D
;
; PURPOSE:
;       This procedure plots to the current device a color map
;       of a given hydro variable for the AMR structure stored in
;       memory. Valid for 2D RAMSES simulations only.
;
; CATEGORY:
;       Plotting routines.
;
; CALLING SEQUENCE:
;       TV2D, Grid, Hydro, TYPE = type, CLT = clt, LOG = log, 
;                          XR = xr, YR = yr, 
;                          LEVELMAX = levelmax, SAVE = save, BYTE = byte,
;                          VMIN = vmin, VMAX = vmax
;
; INPUTS
;       Grid: structure containing the AMR mesh.
;
;       Hydro: structure containing the hydro variables.
;
; OPTIONAL INPUTS:
;       TYPE:   if set, defines the variable to map. If not set (TYPE
;       = 0), the routine plots the level for each cell. TYPE = 1, 2,
;       3, 4... corresponds to hydro variable number 1, 2, 3, 4... 
;
;       CLT: color table number.
;
;       LOG:  if set, the routine uses a logarithmic scaling.
;
;       XR:      if set, defines the map boundaries for the X
;       axis. Default: the whole box.
;
;       YR:      same for the Y axis.
;
;       LEVELMAX: if set, specifies the maximum level to be shown.
;       Default: 0.
;       
;       SAVE:     if set, the image is outputted to the
;       two-dimensional array save. If not, the image is outputted to
;       the current plotting device.
;      
;       BYTE:     if set, the image is converted to bytes (valid if
;       the SAVE keyword is used only).
;       
;       VMIN:     if set, the minimum value to be plotted.
;
;       VMAX:     if set, the maximum value to be plotted.
;
; OUTPUTS:
;       None.
;       
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To plot the mass density as a color map with 5 levels of refinement:
;
;               TV2D, grid, hydro, type=1, levelmax=5, clt=33
;
;       To save the same image into a 2D floating point array:
;
;               d = fltarr(512,512)
;               TV2D, grid, hydro, type=1, levelmax=5, save=d
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro tv2d,grid,hydro,type=type,clt=clt $
                   ,log=log $
                   ,xr=xr,yr=yr $
                   ,lmin=lmin,lmax=lmax $
                   ,save=save,byte=byte $
                   ,vmin=vmin,vmax=vmax,dummy=dummy $
                   ,verbose=verbose,showgrid=showgrid

IF N_PARAMS() NE 2 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'tv2d'
    RETURN
ENDIF

ncpu=grid.ncpu
ndim=grid.ndim
nlevelmax=grid.nlevelmax
ngrid=grid.ngrid

if ndim ne 2 then begin
    print,'Mesh should have 2 dimensions'
    print,'but ndim=',ndim
    return
endif

if keyword_set(clt) then loadct,clt
if not keyword_set(dummy) then dummy=-1.d32
if keyword_set(type) then begin
    if type gt 0 then begin
        if type gt hydro.nvar then begin
            print,'Invalid variable type=',type
            print,'for nvar=',hydro.nvar
            return
        endif
    endif
endif

if not keyword_set(lmin) then leveldown=1 else leveldown=lmin
if not keyword_set(lmax) then levelup=nlevelmax else levelup=lmax
if levelup   lt 1 then levelup   = 1
if leveldown lt 1 then leveldown = 1
if levelup   gt nlevelmax then levelup   = nlevelmax
if leveldown gt nlevelmax then leveldown = nlevelmax
if levelup   lt leveldown then levelup   = leveldown

if not keyword_set(xr) then xr=[0.,1.0]
if not keyword_set(yr) then yr=[0.,1.0]

xmin=MAX([-1,xr(0)]) & xmax=MIN([2,xr(1)])
ymin=MAX([-1,yr(0)]) & ymax=MIN([2,yr(1)])

maxl = max([xmax-xmin,ymax-ymin])
lminn = fix(-alog(maxl)/alog(2.)+0.0001)+1L
leveldown=max([lminn,leveldown])

if not keyword_set(save) then begin
    ximmax=!d.x_size
    yimmax=!d.y_size
endif else begin
    ssz=size(save)
    ximmax=ssz(1)
    yimmax=ssz(2)
endelse

maxpix=max([ximmax,yimmax])
lmaxn = fix(alog(maxpix)/alog(2.))+lminn;-1L
levelup=min([lmaxn,levelup])

print,leveldown,levelup,format='("lmin=",I2," lmax=",I2)'

dxmax=0.5^(leveldown-1)
scale2=2^(levelup-leveldown+1)

xmin=xmin/dxmax & xmax=xmax/dxmax
ymin=ymin/dxmax & ymax=ymax/dxmax
imin=floor(xmin) & imax=ceil(xmax) 
jmin=floor(ymin) & jmax=ceil(ymax) 
nimx=(imax-imin)*scale2 & nimy=(jmax-jmin)*scale2
iimin=fix((xmin-imin)*scale2) & iimax=fix((xmax-xmin)*scale2)+iimin-1
jjmin=fix((ymin-jmin)*scale2) & jjmax=fix((ymax-ymin)*scale2)+jjmin-1
xeff=(iimax-iimin+1) & yeff=(jjmax-jjmin+1)
zxmin=(double(iimin  )/scale2+imin)*dxmax
zxmax=(double(iimax+1)/scale2+imin)*dxmax
zymin=(double(jjmin  )/scale2+jmin)*dxmax
zymax=(double(jjmax+1)/scale2+jmin)*dxmax
xrn=[zxmin,zxmax]
yrn=[zymin,zymax]
print,xrn,yrn,format='("xr=[",E10.3,",",E10.3,"], yr=[",E10.3,",",E10.3,"]")'
print,xeff,yeff,format='("Computing image of size=",i4,"x",i4)'

im=FLTARR(nimx,nimy)+dummy

for icpu=0,ncpu-1 do begin
for ilevel=leveldown,levelup do begin

    dx=0.5d0^(ilevel-leveldown+1)
    scale_loc=fix(scale2*dx)
    nimx_loc=2^(ilevel-leveldown+1)*(imax-imin)
    nimy_loc=2^(ilevel-leveldown+1)*(jmax-jmin)
    if keyword_set(verbose) then $
      print,ilevel,ngrid(ilevel-1,icpu),format='("Level=",i2," ngrid=",i6)'

    if scale_loc gt 0 and ngrid(ilevel-1,icpu) gt 0 then begin

        mesh2=(*grid.level[ilevel-1,icpu])
        x=mesh2.xg(*,0)/dxmax
        y=mesh2.xg(*,1)/dxmax
        ind=where( x+dx gt xmin and x-dx lt xmax and $
                   y+dx gt ymin and y-dx lt ymax)
        if(ind(0) ne -1) then begin
            x=mesh2.xg(ind,0)/dxmax-dx-imin
            y=mesh2.xg(ind,1)/dxmax-dx-jmin
            subim=fltarr(nimx_loc,nimy_loc)+dummy
            x=x*2^(ilevel-leveldown+1)
            y=y*2^(ilevel-leveldown+1)
            ind_x=fix(x)
            ind_y=fix(y)

            for i=0,1 do begin
                for j=0,1 do begin
                    ind_cell=i+2*j
                    active=mesh2.son(ind,ind_cell)
                    if(ilevel lt levelup) then begin
                        ind2=where(active eq 0)
                    endif else begin
                        ind2=where(active gt -1)
                    endelse

                    if not keyword_set(type)then begin
                        q=ilevel+0.d0*mesh2.son(ind,ind_cell)
                    endif else if type lt 0 then begin
                        q=icpu+0.d0*mesh2.son(ind,ind_cell)
                    endif else begin
                        mesh=(*hydro.levelh[ilevel-1,icpu])
                        q=mesh.u(ind,ind_cell,type-1)
                    endelse
                    
                    if (ind2(0) ne -1) then begin
                        ind_xx=ind_x(ind2)+i
                        ind_yy=ind_y(ind2)+j
                        for ii2=0L,n_elements(ind2)-1L do begin
                            subim(ind_xx(ii2),ind_yy(ii2))= q(ind2(ii2))
                        endfor
                    endif

                endfor
             endfor
            im=im > ( REBIN(subim,nimx,nimy,/sample) > im)
            
       endif
    endif
endfor
endfor
im=im(iimin:iimax,jjmin:jjmax)

xscale=double(ximmax)/double(xeff)
yscale=double(yimmax)/double(yeff)
lscale=min([xscale,yscale])
;imtv=CONGRID(im,lscale*xeff,lscale*yeff)
imtv=REBIN(im,lscale*xeff,lscale*yeff,/sample)


inddummy=where(imtv le dummy,ndummy)
notdummy=where(imtv gt dummy,nnotdummy)
if(ndummy gt 0)then begin
    TVLCT,r,g,b,/get
;    r(0) = 0   & g(0) = 0   & b(0) = 0 ; reserve for black
;    r(1) = 255 & g(1) = 255 & b(1) = 255 ; reserve for white
    r(2) = 175 & g(2) = 175 & b(2) = 175 ; reserve for neutral grey
    TVLCT,r,g,b
endif

vmin0=0. & vmax0=0.
if nnotdummy gt 0 then vmax0=max(imtv(notdummy))
if nnotdummy gt 0 then vmin0=min(imtv(notdummy))
print,vmin0,vmax0,format='("Min=",E10.3," Max=",E10.3)'

if not keyword_set(vmax) then vmax=vmax0
if not keyword_set(vmin) then vmin=vmin0

vmax1=vmax
vmin1=vmin

if keyword_set(log) then begin
    if nnotdummy gt 0 then begin
        imtv(notdummy)=alog10(imtv(notdummy))
        vmax1=alog10(vmax)
        vmin1=alog10(vmin)
    endif
endif

if not keyword_set(save) then begin
    image=bytscl(imtv,top=MIN([!d.table_size,256])-35 $
                 ,min=vmin1,max=vmax1)+byte(33)
    if(ndummy gt 0)then image(inddummy)=2
    tv,image
    if keyword_set(showgrid) then begin
        TVLCT,r,g,b,/get
        mmm = !d.table_size - 1
        r(mmm) = 0 & g(mmm) = 0 & b(mmm) = 0 ; reserve for white
        TVLCT,r,g,b
        zoom=max([(xrn(1)-xrn(0)),(yrn(1)-yrn(0))])
        pp_amr2d,grid,xr=xr,yr=yr,x0=xrn(0)+0.5*zoom,y0=yrn(0)+0.5*zoom $
          ,scale=zoom,/noer,lmin=leveldown,lmax=levelup,/expand
    endif
endif
if keyword_set(save) then begin
    if keyword_set(byte)then begin
        image=bytscl(imtv,top=MIN([!d.table_size,256])-35,min=vmin1,max=vmax1)+byte(33)
        if(ndummy gt 0)then image(inddummy)=2
        save=image
    endif else begin
        save=imtv
    endelse
endif
imtv=1.
im=1.
subim=1.

return
end
;###################################################
;###################################################
;###################################################
