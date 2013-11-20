;+
; NAME:
;       RAY3D
;
; PURPOSE:
;       This procedure performs an orthogonal ray-tracing across the
;       AMR mesh for a given emissivity (specified by the user) or for
;       one of the hydro variable, and output a color map to the
;       current graphic device. Valid for 3D RAMSES simulations only. 
;
; CATEGORY:
;       Plotting routines.
;
; CALLING SEQUENCE:
;       RAY3D, Grid, Hydro, 
;          TYPE = type, CLT = clt, LOG = log,
;          XR = xr, YR = yr, ZR = zr, AVERAGE = average,
;          XPROJ = xproj, YPROJ = YPROJ, ZPROJ = zproj, 
;          LMIN = lmin, LMAX = lmax, SAVE = save,
;          VMIN = vmin, VMAX = vmax
;
; INPUTS:
;       Grid: structure containing the AMR mesh.
;
;       Hydro: structure containing the hydro variables.
;       
; OPTIONAL INPUTS:
;       TYPE:   character string. If set, defines the type of hydro
;       variable to map, with TYPE = '1', '2', '3', etc. If not set,
;       the routine plots the level of refinement (Default). If TYPE =
;       'userdef', the routine calls procedure emissivity.pro that
;       must be provided and compiled by the user.
;
;       CLT: color table number.
;
;       LOG:  if set, the routine uses a logarithmic scaling.
;
;       AVERAGE:  if set, an average projection method is used. If not
;       (default), the maximum of the chosen variable along each ray
;       is projected.
;
;       XR:      if set, defines the volume boundaries for the X
;       axis. Default: the 2 central cells only.
;
;       YR:      same for the Y axis.
;
;       ZR:      same for the Y axis.
;
;       XPROJ:   set the projection axis to be the X axis.
;
;       YPROJ:   set the projection axis to be the Y axis.
;
;       ZPROJ:   set the projection axis to be the Z axis.
;
;       LMIN:    if set, specifies the minimum level to be shown.
;       Default: 1.
;       
;       LMAX:    if set, specifies the maximum level to be shown.
;       Default: 1.
;       
;       SAVE:     if set, the image is outputted to the
;       two-dimensional array save. If not, the image is outputted to
;       the current plotting device.
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
;       To project the mass density along the Z axis for the whole box
;       and output the color map to the current graphic device in a
;       logarithmic scaling:
;
;               RAY3D, Grid, Hydro, xr=[0,1], yr=[0,1],
;               zr=[0,1], /zproj, 
;               type=1, lmax=3, /log
;
;       To project the mass density along the Z axis for a thin slice
;       of 1 coarse cell thickness and output the color map to a 2D
;       array:
;
;               d = fltarr(512,512)
;               RAY3D, Grid, Hydro, xr=[0,1], yr=[0,1],
;               zr=[0.4,0.6], /zproj, 
;               type=1, lmax=3, save=d
;
;       To define an external emissivity, load first a procedure
;       called "emissivity" that contains the desired emission
;       law. Here is an example of such a routine:
;
;               PRO EMISSIVITY, Var, Emi
;               d = Var(*,0)
;               Emi = d^2
;               end
;
;       Here, on input, Var is a 2D array with the index of the hydro
;       variables as 2nd dimension. Emi is a 1D array containing the
;       specified emission law. Then type:
;
;
;               RAY3D, Grid, Hydro, xr=[0,1], yr=[0,1],
;               zr=[0,1], /zproj, 
;               type='userdef', lmax=3, /log
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro ray3d,grid,hydro $
        ,type=type,clt=clt $
        ,log=log,average=average $
        ,xr=xr,yr=yr,zr=zr $
        ,xproj=xproj,yproj=yproj,zproj=zproj $
        ,lmin=lmin,lmax=lmax,save=save $
        ,vmin=vmin,vmax=vmax,verbose=verbose,dummy=dummy

IF N_PARAMS() NE 2 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'ray3d'
    RETURN
ENDIF

ndim=grid.ndim
nlevelmax=grid.nlevelmax
ngrid=grid.ngrid
nvar=hydro.nvar
ncpu=grid.ncpu

if ndim ne 3 then begin
    print,'Mesh should have 3 dimensions'
    print,'but ndim=',ndim
    return
endif

; Check that type is character
if keyword_set(type) then type=string(type)

; Set up color table
if keyword_set(clt) then begin
    loadct,clt
endif

; Set up image size
if not keyword_set(save)then begin
    ximax=!d.x_size
    yimax=!d.y_size
endif else begin
    ss=size(save)
    ximax=ss(1)
    yimax=ss(2)
endelse

; Set up projection geometry and scale
if not keyword_set(lmin) then leveldown=1 else leveldown=lmin
if not keyword_set(lmax) then levelup  =1 else levelup  =lmax
if leveldown lt 1 then leveldown = 1
if leveldown gt nlevelmax then leveldown = nlevelmax
if levelup   lt leveldown then levelup   = leveldown
if levelup   gt nlevelmax then levelup   = nlevelmax

if not keyword_set(xr) then xr=[0.,1.]
if not keyword_set(yr) then yr=[0.,1.]
if not keyword_set(zr) then zr=[0.,1.]

range=fltarr(2,3)
range(0:1,0)=xr
range(0:1,1)=yr
range(0:1,2)=zr

if keyword_set(zproj) then begin
    idir=0 & jdir=1 & kdir=2
endif else if keyword_set(yproj) then begin
    idir=0 & jdir=2 & kdir=1
endif else if keyword_set(xproj) then begin
    idir=1 & jdir=2 & kdir=0
endif else begin
    print,'Choose projection axis'
    return
endelse

; rotate axis
xmin=MAX([-1,range(0,idir)]) & xmax=MIN([2.0,range(1,idir)])
ymin=MAX([-1,range(0,jdir)]) & ymax=MIN([2.0,range(1,jdir)])
zmin=MAX([-1,range(0,kdir)]) & zmax=MIN([2.0,range(1,kdir)])

; cell size
dxmax=0.5^leveldown
dxmin=0.5^levelup
scale2=dxmax/dxmin

; geometry
xmin=xmin/dxmax  & xmax=xmax/dxmax
ymin=ymin/dxmax  & ymax=ymax/dxmax
zmin=zmin/dxmax  & zmax=zmax/dxmax

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
zzmin=zmin*dxmax
zzmax=zmax*dxmax

; new range
rangen=range
rangen(0:1,idir)=[zxmin,zxmax]
rangen(0:1,jdir)=[zymin,zymax]
rangen(0:1,kdir)=[zzmin,zzmax]
print,rangen,format='("xr=[",E10.3,",",E10.3,"], yr=[",E10.3,",",E10.3,"], zr=[",E10.3,",",E10.3,"]")'

; image size
nxmem=nimx  & nymem=nimy
immin=iimin & immax=iimax
jmmin=jjmin & jmmax=jjmax
xsize=xeff  & ysize=yeff
hhh=(zmax-zmin)
if not keyword_set(average) then hhh=1.0
print,xsize,ysize,format='("Projection on image of size=",i4,"x",i4)'
imtv=fltarr(nxmem,nymem)
if not keyword_set(average) then imtv=imtv-1.d15

; project AMR levels
for icpu=0,ncpu-1 do begin
for ilevel=leveldown,levelup do begin

    dx=0.5d0^ilevel
    dx_loc=dx/dxmax
    scale_loc=fix(dxmax/dx)
    hhh_loc=dx_loc/hhh
    nimx_loc=scale_loc*(imax-imin)
    nimy_loc=scale_loc*(jmax-jmin)
    if keyword_set(verbose) then print,ilevel,icpu,ngrid(ilevel-1,icpu)$
      ,format='("Level=",i2," cpu=",i2," ngrid=",i6)'

    if scale_loc gt 0 and ngrid(ilevel-1,icpu) gt 0 then begin

        x=(*grid.level[ilevel-1,icpu]).xg(*,idir)/dxmax
        y=(*grid.level[ilevel-1,icpu]).xg(*,jdir)/dxmax
        z=(*grid.level[ilevel-1,icpu]).xg(*,kdir)/dxmax
        ind=where( x ge xmin-dx_loc and x le xmax+dx_loc and $
                   y ge ymin-dx_loc and y le ymax+dx_loc and $
                   z ge zmin-dx_loc and z le zmax+dx_loc, nok)

        if(nok gt 0) then begin

            x=(*grid.level[ilevel-1,icpu]).xg(ind,idir)/dxmax-dx_loc-imin
            y=(*grid.level[ilevel-1,icpu]).xg(ind,jdir)/dxmax-dx_loc-jmin
            z=(*grid.level[ilevel-1,icpu]).xg(ind,kdir)/dxmax-dx_loc
            subimtv=fltarr(nimx_loc,nimy_loc)
            if not keyword_set(average) then subimtv=subimtv-1.d15
            x=x*scale_loc
            y=y*scale_loc
            ind_x=fix(x)
            ind_y=fix(y)

            for k=0,1 do begin
            for j=0,1 do begin
            for i=0,1 do begin

                ind_cell=i*2^idir+j*2^jdir+k*2^kdir
                active=(*grid.level[ilevel-1,icpu]).son(ind,ind_cell)
                if(ilevel lt levelup) then begin
                    ind2=where(active eq  0, nleaf)
                endif else begin
                    ind2=where(active gt -1, nleaf)
                endelse
                
                if not keyword_set(type)then begin
                    q=ilevel+0.d0*active
                endif else if type eq 'userdef' then begin
                    mesh=(*hydro.levelh[ilevel-1,icpu])
                    variables=mesh.u(ind,ind_cell,0:nvar-1)
                    variables=reform(variables,n_elements(ind),nvar)
                    emissivity,variables,q
                endif else if type lt 0 then begin
                    q=icpu+0.d0*active
                endif else begin
                    mesh=(*hydro.levelh[ilevel-1,icpu])
                    q=mesh.u(ind,ind_cell,type-1)
                endelse
                
                if (nleaf gt 0) then begin
                    ind_xx=ind_x(ind2)+i
                    ind_yy=ind_y(ind2)+j
                    zc1=z(ind2)+(k  )*dx_loc
                    zc2=z(ind2)+(k+1)*dx_loc
                    indc=where(zc2 ge zmax,nc)
                    if nc gt 0 then zc2(indc)=zmax
                    indc=where(zc1 le zmin,nc)
                    if nc gt 0 then zc1(indc)=zmin
                    weight=(zc2-zc1)/dx_loc
                    indc=where(weight le 0.,nc)
                    if nc gt 0 then weight(indc)=0.
                    okcell=ind_xx ge 0        and $
                      ind_xx lt nimx_loc and $
                      ind_yy ge 0        and $
                      ind_yy lt nimy_loc
                    if not keyword_set(average) then begin
                        for ii2=0L,n_elements(ind2)-1L do begin
                            if(okcell(ii2)) then begin
                                subimtv(ind_xx(ii2),ind_yy(ii2))= $
                                  subimtv(ind_xx(ii2),ind_yy(ii2)) $
                                  > ( q(ind2(ii2) ) > $
                                      subimtv(ind_xx(ii2),ind_yy(ii2)) )
                            endif
                        endfor
                    endif else begin
                        for ii2=0L,n_elements(ind2)-1L do begin
                            if(okcell(ii2)) then begin
                                subimtv(ind_xx(ii2),ind_yy(ii2))= $
                                  subimtv(ind_xx(ii2),ind_yy(ii2)) $
                                  +q(ind2(ii2))*hhh_loc*weight(ii2)
                            endif
                        endfor
                    endelse
                endif
                
            endfor
            endfor
            endfor
            if keyword_set(average) then begin
                imtv=imtv+REBIN(subimtv,nxmem,nymem,/sample)
            endif else begin
                imtv=imtv > ( REBIN(subimtv,nxmem,nymem,/sample) > imtv )
            endelse
        endif

    endif

endfor
endfor
imtv=imtv(immin:immax,jmmin:jmmax)

ndummy=0
if keyword_set(dummy) then begin
    inddummy=where(imtv le dummy,ndummy)
    notdummy=where(imtv gt dummy,nnotdummy)
    if(ndummy gt 0)then begin
        TVLCT,r,g,b,/get
        r(2) = 175 & g(2) = 175 & b(2) = 175 ; reserve for neutral grey
        TVLCT,r,g,b
    endif
    vmax0=max(imtv(notdummy))
    vmin0=min(imtv(notdummy))
endif else begin
    vmax0=max(imtv)
    vmin0=min(imtv)
endelse
print,vmin0,vmax0,format='("Min=",E10.3," Max=",E10.3)'

xscale=float(ximax)/xsize
yscale=float(yimax)/ysize
lscale=min([xscale,yscale])
imtv=CONGRID(imtv,lscale*xsize,lscale*ysize)

if keyword_set(save)then begin
    if keyword_set(verbose) then $
      print,'Outputing image to variable'
    save=imtv
    return
endif else begin
    print,'Outputing image to screen'
    if (keyword_set(vmin)) then $
      print,vmin,format='("New min=",E10.3)'
    if (keyword_set(vmax)) then $
      print,vmax,format='("New max=",E10.3)'
    if not keyword_set(vmin)then vmin=vmin0
    if not keyword_set(vmax)then vmax=vmax0
    if keyword_set(log) then begin
        indok=where(imtv gt 0.,nn1)
        notok=where(imtv le 0.,nn2)
        vmin=MAX([vmin,min(imtv(indok))])
        imtv(indok)=alog10(imtv(indok))
        if(nn2 gt 0) then imtv(notok)=alog10(vmin)
        vmin=alog10(vmin)
        vmax=alog10(vmax)
    endif
    image=bytscl(imtv,top=MIN([!d.n_colors,256])-35,min=vmin,max=vmax)+byte(33)
    if(ndummy gt 0)then image(inddummy)=2
    tv,image
endelse

; Free memory
x=0
y=0
z=0
q=0
imtv=0
image=0
subim=0
subimtv=0
mesh=0
mesh2=0
ind=0
ind2=0
ind_x=0
ind_y=0
ind_xx=0
ind_yy=0
active=0
variables=0


return
end
;###################################################
;###################################################
;###################################################
