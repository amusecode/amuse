;+
; NAME:
;       MESHIJK
;
; PURPOSE:
;       This functions returns a regular 3D Cartesian array
;       from a 3D AMR data set and for the chosen variable. The array
;       size is computed automatically to contains all the desired
;       region (XR, YR, ZR) at the desired resolution (LEVELMAX).
;
; CATEGORY:
;       Data analysis.
;
; CALLING SEQUENCE:
;       Cube = MESHIJK(Grid3d, Hydro3d, XR = xr, YR = yr, ZR = zr,
;               TYPE=type, LMIN=lmin, LMAX=lmax)
;
; INPUTS
;       Grid3d: structure containing the 3D AMR mesh.
;
;       Hydro3d: structure containing the 3D hydro variables.
;
; OPTIONAL INPUTS:
;       XR:      if set, defines the map boundaries for the X
;       axis. Default: the whole box. 
;
;       YR:      same for the Y axis.
;
;       ZR:      same for the Z axis.
;
;       TYPE:    The variable to extract (1 for density, 2 for vx,
;       etc...). For zero or a negative value, the routine extract an
;       integer corresponding to the AMR level. Default: TYPE = 0
;
;       LMIN:    The minimum resolution (coarser level). Default: 1.
;
;       LMAX:    The maximum resolution (finer level). Default: 1.
;
; OUTPUTS:
;       Cube: a structure containing the 3D regular cartesian array of
;       the extracted variable, the type of the extracted variable,
;       the x, y and z range, and the array dimensions.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To extract an array in the center of the box and for the
;       density, type: 
;
;               Cube = MESHIJK(Grid3d, Hydro3d, xr=[0.4,0.6],
;               yr=[0.4,0.6], zr=[0.4,0.6], type=1, lmin=2, lmax=6)
;
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
function meshijk, grid, hydro, type=type, xr=xr, yr=yr, zr=zr $
                  , lmin=lmin, lmax=lmax

IF N_PARAMS() NE 2 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'meshijk'
    RETURN,0.
ENDIF

ncpu=grid.ncpu
ndim=grid.ndim
nlevelmax=grid.nlevelmax
ngrid=grid.ngrid
nvar=hydro.nvar

if ndim ne 3 then begin
    print,'Mesh should have 3 dimensions'
    print,'but ndim=',ndim
    return,0.
endif

; Set variable to extract
if keyword_set(type) then type=string(type)

; Set up levels to consider
if not keyword_set(lmin) then leveldown=1 else leveldown=lmin
if not keyword_set(lmax) then levelup=1 else levelup=lmax
if levelup lt 1 then levelup = 1
if leveldown lt 1 then leveldown = 1
if levelup   gt nlevelmax then levelup   = nlevelmax
if leveldown gt nlevelmax then leveldown = nlevelmax

; Set up cube size
if not keyword_set(xr) then xr=[0.,1.0]
if not keyword_set(yr) then yr=[0.,1.0]
if not keyword_set(zr) then zr=[0.,1.0]

xmin=MAX([0.,xr(0)]) & xmax=MIN([1.0,xr(1)])
ymin=MAX([0.,yr(0)]) & ymax=MIN([1.0,yr(1)])
zmin=MAX([0.,zr(0)]) & zmax=MIN([1.0,zr(1)])

; Set resolutions
dxmax=0.5^leveldown
dxmin=0.5^levelup
scale2=dxmax/dxmin

xmin=xmin/dxmax & xmax=xmax/dxmax
ymin=ymin/dxmax & ymax=ymax/dxmax
zmin=zmin/dxmax & zmax=zmax/dxmax
imin=floor(xmin) & imax=ceil(xmax) 
jmin=floor(ymin) & jmax=ceil(ymax) 
kmin=floor(zmin) & kmax=ceil(zmax) 
nimx=(imax-imin)*scale2
nimy=(jmax-jmin)*scale2
nimz=(kmax-kmin)*scale2 
xmin=imin & xmax=imax
ymin=jmin & ymax=jmax
zmin=kmin & zmax=kmax
xrnew=[imin*dxmax,imax*dxmax]
yrnew=[jmin*dxmax,jmax*dxmax]
zrnew=[kmin*dxmax,kmax*dxmax]

print,nimx,nimy,nimz,format='("Computing array of size=",i3,"x",i3,"x",i3)'

save=DBLARR(nimx,nimy,nimz)

;======================
; Project AMR levels
;======================
for icpu=0,ncpu-1 do begin
    for ilevel=leveldown,levelup do begin

        dx=0.5d0^(ilevel)/dxmax
        nimx_loc=(imax-imin)/dx
        nimy_loc=(jmax-jmin)/dx
        nimz_loc=(kmax-kmin)/dx
        ;print,ilevel,ngrid(ilevel-1,icpu),format='("Level=",i2," ngrid=",i6)'
        print,(icpu+1),ilevel,ngrid(ilevel-1,icpu),format='("cpu=",i3," Level=",i2," ngrid=",i6)'
        if ngrid(ilevel-1,icpu) gt 0 then begin
            mesh=(*hydro.levelh[ilevel-1,icpu])
            x=(*grid.level[ilevel-1,icpu]).xg(*,0)/dxmax
            y=(*grid.level[ilevel-1,icpu]).xg(*,1)/dxmax
            z=(*grid.level[ilevel-1,icpu]).xg(*,2)/dxmax
            ind=where( x+dx gt xmin and x-dx lt xmax and $
                       y+dx gt ymin and y-dx lt ymax and $
                       z+dx gt zmin and z-dx lt zmax, nok)
            
            if(nok gt 0) then begin
                    
                x=(x(ind)-dx-imin)/dx
                y=(y(ind)-dx-jmin)/dx
                z=(z(ind)-dx-kmin)/dx
                
                subcube=dblarr(nimx_loc,nimy_loc,nimz_loc)
                
                ind_x=fix(x)
                ind_y=fix(y)
                ind_z=fix(z)

                for k=0,1 do begin
                for j=0,1 do begin
                for i=0,1 do begin
                    ind_cell=i+2*j+4*k
                    active=(*grid.level[ilevel-1,icpu]).son(ind,ind_cell)

                    if not keyword_set(type) then begin
                        q=ilevel+0.*mesh.u(ind,ind_cell,1)
                    endif else if type eq 'userdef' then begin
                        variables=mesh.u(ind,ind_cell,0:nvar-1)
                        variables=reform(variables,n_elements(ind),nvar)
                        emissivity,variables,q
                    endif else if(type gt 0 and type le nvar)then begin
                        q=mesh.u(ind,ind_cell,type-1)
                    endif else begin
                        q=ilevel+0.*mesh.u(ind,ind_cell,1)
                    endelse
                
                    if(ilevel lt levelup) then begin
                        ind2=where(active eq 0, nok2)
                    endif else begin
                        ind2=where(active gt -1, nok2)
                    endelse
                    
                    if (nok2 gt 0) then begin

                    ind_xx=ind_x(ind2)+i
                    ind_yy=ind_y(ind2)+j
                    ind_zz=ind_z(ind2)+k

                    for ii2=0L,n_elements(ind2)-1L do begin
                    if(    ind_xx(ii2) ge 0 and ind_xx(ii2) lt nimx_loc $
                       and ind_yy(ii2) ge 0 and ind_yy(ii2) lt nimy_loc $
                       and ind_zz(ii2) ge 0 and ind_zz(ii2) lt nimz_loc )then $
                      subcube(ind_xx(ii2),ind_yy(ii2),ind_zz(ii2))= $
                      subcube(ind_xx(ii2),ind_yy(ii2),ind_zz(ii2)) $
                      +q(ind2(ii2))
                    endfor

                    endif
                endfor
                endfor
                endfor
                save=TEMPORARY(save)+REBIN(subcube,nimx,nimy,nimz,/sample)
            endif
        endif
    endfor
endfor

; Free memory
subcube=0
active=0
q=0
mesh=0
x=0
y=0
z=0
ind=0
ind2=0
ind_x=0
ind_xx=0
ind_y=0
ind_yy=0
ind_z=0
ind_zz=0

; Store in sructure
if not keyword_set(type) then type=0
cube={nx:fix(nimx), ny:fix(nimy), nz:fix(nimz), xr:xrnew, yr:yrnew, zr:zrnew, $
      type:type, data:TEMPORARY(save)}

return,cube
end
;###################################################
;###################################################
;###################################################
