;+
; NAME:
;       AMR2CELL
;
; PURPOSE:
;       This procedure extract all data from HYDRO and AMR structures
;       as a list of "leaf" cells (refined cells are discarded).
;
; CATEGORY:
;       Data analysis.
;
; CALLING SEQUENCE:
;       AMR2CELL, Grid, Hydro, Cell,
;                          XR = xr, YR = yr, ZR = zr, 
;                          LMIN=lmin, LMAX = lmax, VERBOSE = verbose
;
; INPUTS
;       Grid: structure containing the AMR mesh.
;
;       Hydro: structure containing the hydro variables.
;
; OPTIONAL INPUTS:
;       XR:      if set, defines the map boundaries for the X
;       axis. Default: the whole box.
;
;       YR:      same for the Y axis.
;
;       ZR:      same for the Z axis.
;
;       LMIN:    if set, specifies the minimum level to be shown.
;       Default: 1.
;
;       LMAX:    if set, specifies the maximum level to be shown.
;       Default: 1.
;
;       VERBOSE: if set, output infos to screen. Default: not set.
;
; OUTPUTS:
;       Cell: structure containing the cell variables. The
;             variables are the followings: 
;             Point.n    = (Integer scalar) number of cells.
;             Point.dx   = (1D Float array) size of each cell.
;             Point.x    = (1D Float array) x position of each cell.
;             Point.y    = (1D Float array) y position of each cell.
;             Point.z    = (1D Float array) z position of each cell.
;             Point.var  = (2D Float array) hydro variables of each cell.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To extract AMR cells up to level 5, type:
;
;               AMR2CELL, grid, hydro, cell, lmax=5
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
;###################################################
;###################################################
;###################################################
pro amr2cell,grid,hydro,amr $
         ,xr=xr,yr=yr,zr=zr $
         ,lmin=lmin,lmax=lmax, verbose=verbose

IF N_PARAMS() NE 3 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'amrcontour'
    RETURN
ENDIF

ncpu=grid.ncpu
ndim=grid.ndim
nlevelmax=grid.nlevelmax
ngrid=grid.ngrid
nvar=hydro.nvar

if ndim ne 3 and ndim ne 2then begin
    print,'Mesh should have 2 or 3 dimensions'
    print,'but ndim=',ndim
    return
endif

; Set up geometry and scale
if not keyword_set(lmin) then leveldown=1 else leveldown=lmin
if not keyword_set(lmax) then levelup=nlevelmax else levelup=lmax
if levelup   lt 1 then levelup   = 1
if leveldown lt 1 then leveldown = 1
if levelup   gt nlevelmax then levelup   = nlevelmax
if leveldown gt nlevelmax then leveldown = nlevelmax
if not keyword_set(xr) then xr=[0.,1.]
if not keyword_set(yr) then yr=[0.,1.]
if not keyword_set(zr) then zr=[0.,1.]
xmin=MAX([0.,xr(0)]) & xmax=MIN([1.,xr(1)])
ymin=MAX([0.,yr(0)]) & ymax=MIN([1.,yr(1)])
zmin=MAX([0.,zr(0)]) & zmax=MIN([1.,zr(1)])

;===========================================================
if keyword_set(verbose) then print,'Counting leaf cells...'
;===========================================================

; Loop over levels
nleaf=0L

if (ndim eq 3) then begin

for ilevel=leveldown,levelup do begin
    for icpu=0,ncpu-1 do begin
        dx=0.5d0^ilevel
        if(ngrid(ilevel-1,icpu) gt 0)then begin

        mesh2=(*grid.level[ilevel-1,icpu])
        x=mesh2.xg(*,0)
        y=mesh2.xg(*,1)
        z=mesh2.xg(*,2)
        ind=where( x+dx gt xmin and x-dx lt xmax and $
                   y+dx gt ymin and y-dx lt ymax and $
                   z+dx gt zmin and z-dx lt zmax, nkeep)
        if(nkeep gt 0) then begin
            for i=0,1 do begin
            for j=0,1 do begin
            for k=0,1 do begin
                ind_cell=i+2*j+4*k
                active=mesh2.son(ind,ind_cell)
                nok1=n_elements(active)
                if ilevel lt levelup then indok=where(active eq 0, nok1)
                nleaf=nleaf+nok1
            endfor
            endfor
            endfor
        endif

        if keyword_set(verbose) then $
          print,ilevel,nleaf,format='("Level=",i2," nleaf=",i8)'
        endif
    endfor
endfor

; Define AMR structure
amr={n:nleaf,dx:fltarr(nleaf),x:fltarr(nleaf),y:fltarr(nleaf),z:fltarr(nleaf),var:fltarr(nleaf,nvar)}

endif else begin

for ilevel=leveldown,levelup do begin
    for icpu=0,ncpu-1 do begin
        dx=0.5d0^ilevel
        if(ngrid(ilevel-1,icpu) gt 0)then begin

        mesh2=(*grid.level[ilevel-1,icpu])
        x=mesh2.xg(*,0)
        y=mesh2.xg(*,1)
        ind=where( x+dx gt xmin and x-dx lt xmax and $
                   y+dx gt ymin and y-dx lt ymax, nkeep)
        if(nkeep gt 0) then begin
            for i=0,1 do begin
            for j=0,1 do begin
                ind_cell=i+2*j
                active=mesh2.son(ind,ind_cell)
                nok1=n_elements(active)
                if ilevel lt levelup then indok=where(active eq 0, nok1)
                nleaf=nleaf+nok1
            endfor
            endfor
        endif

        if keyword_set(verbose) then $
          print,ilevel,nleaf,format='("Level=",i2," nleaf=",i8)'
        endif
    endfor
endfor

; Define AMR structure
amr={n:nleaf,dx:fltarr(nleaf),x:fltarr(nleaf),y:fltarr(nleaf),var:fltarr(nleaf,nvar)}


endelse




;=============================
if keyword_set(verbose) then print,'Building point list...'
;=============================

; Loop over levels
nstart=0L

if (ndim eq 3) then begin

for ilevel=leveldown,levelup do begin
    for icpu=0,ncpu-1 do begin
    dx=0.5d0^ilevel
    if ngrid(ilevel-1,icpu) gt 0 then begin

    if keyword_set(verbose) then print,ilevel,ngrid(ilevel-1,icpu),format='("Level=",i2," ngrid=",i6)'
    mesh =(*hydro.levelh[ilevel-1,icpu])
    mesh2=(*grid.level  [ilevel-1,icpu])
    x=mesh2.xg(*,0)
    y=mesh2.xg(*,1)
    z=mesh2.xg(*,2)
    ind=where( x+dx gt xmin and x-dx lt xmax and $
               y+dx gt ymin and y-dx lt ymax and $
               z+dx gt zmin and z-dx lt zmax, nkeep)
    if(nkeep gt 0) then begin
        x=mesh2.xg(ind,0)
        y=mesh2.xg(ind,1)
        z=mesh2.xg(ind,2)
        for i=0,1 do begin
        for j=0,1 do begin
        for k=0,1 do begin
            ind_cell=i+2*j+4*k
            active=mesh2.son(ind,ind_cell)
            xc=x+(i-0.5)*dx
            yc=y+(j-0.5)*dx
            zc=z+(k-0.5)*dx
            nok1=n_elements(active)
            if levelup gt ilevel then $
              indok=where(active eq 0, nok1) else $
              indok=LINDGEN(n_elements(active))
            if nok1 gt 0L then begin
                amr.dx(nstart:nstart+nok1-1L)=dx
                amr.x (nstart:nstart+nok1-1L)=xc(indok)
                amr.y (nstart:nstart+nok1-1L)=yc(indok)
                amr.z (nstart:nstart+nok1-1L)=zc(indok)
                for ivar=0,nvar-1 do begin
                    tmp=mesh.u(ind(indok),ind_cell,ivar)
                    amr.var(nstart:nstart+nok1-1L,ivar)=tmp
                endfor
                tmp=0
            endif
            nstart=nstart+nok1
        endfor          
        endfor
        endfor
    endif

    endif
    endfor
endfor

endif else begin


for ilevel=leveldown,levelup do begin
    for icpu=0,ncpu-1 do begin
    dx=0.5d0^ilevel
    if ngrid(ilevel-1,icpu) gt 0 then begin
    if keyword_set(verbose) then print,ilevel,ngrid(ilevel-1,icpu),format='("Level=",i2," ngrid=",i6)'
    mesh =(*hydro.levelh[ilevel-1,icpu])
    mesh2=(*grid.level  [ilevel-1,icpu])
    x=mesh2.xg(*,0)
    y=mesh2.xg(*,1)
    ind=where( x+dx gt xmin and x-dx lt xmax and $
               y+dx gt ymin and y-dx lt ymax, nkeep)
    if(nkeep gt 0) then begin
        x=mesh2.xg(ind,0)
        y=mesh2.xg(ind,1)
        for i=0,1 do begin
        for j=0,1 do begin
            ind_cell=i+2*j
            active=mesh2.son(ind,ind_cell)
            xc=x+(i-0.5)*dx
            yc=y+(j-0.5)*dx
            nok1=n_elements(active)
            if levelup gt ilevel then $
              indok=where(active eq 0, nok1) else $
              indok=LINDGEN(n_elements(active))
            if nok1 gt 0L then begin
                amr.dx(nstart:nstart+nok1-1L)=dx
                amr.x (nstart:nstart+nok1-1L)=xc(indok)
                amr.y (nstart:nstart+nok1-1L)=yc(indok)
                for ivar=0,nvar-1 do begin
                    tmp=mesh.u(ind(indok),ind_cell,ivar)
                    amr.var(nstart:nstart+nok1-1L,ivar)=tmp
                endfor
                tmp=0
            endif
            nstart=nstart+nok1
        endfor
        endfor
    endif

    endif
    endfor
endfor

endelse


; Free memory
mesh=0
mesh2=0
active=0
subim=0
indok=0
var=0
x=0
y=0
z=0
ddd=0
xxx=0
yyy=0
zzz=0
vvv=0

return
end
;###################################################
;###################################################
;###################################################
