;+
; NAME:
;       BFIELDIJK
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
function bfieldijk, grid, hydro, xr=xr, yr=yr, zr=zr $
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

imin=floor(xmin/dxmax) & imax=ceil(xmax/dxmax) 
jmin=floor(ymin/dxmax) & jmax=ceil(ymax/dxmax) 
kmin=floor(zmin/dxmax) & kmax=ceil(zmax/dxmax) 
nimx=(imax-imin)*scale2
nimy=(jmax-jmin)*scale2
nimz=(kmax-kmin)*scale2 
xmin=imin*dxmax & xmax=imax*dxmax
ymin=jmin*dxmax & ymax=jmax*dxmax
zmin=kmin*dxmax & zmax=kmax*dxmax
xrnew=[xmin,xmax]
yrnew=[ymin,ymax]
zrnew=[zmin,zmax]

print,nimx,nimy,nimz,format='("Computing array of size=",i3,"x",i3,"x",i3)'

bxl=FLTARR(nimx,nimy,nimz)
bxr=FLTARR(nimx,nimy,nimz)
byl=FLTARR(nimx,nimy,nimz)
byr=FLTARR(nimx,nimy,nimz)
bzl=FLTARR(nimx,nimy,nimz)
bzr=FLTARR(nimx,nimy,nimz)

;======================
; Project AMR levels
;======================
for icpu=0,ncpu-1 do begin
    for ilevel=leveldown,levelup do begin

        dx=0.5d0^(ilevel)

        print,ilevel,ngrid(ilevel-1,icpu),format='("Level=",i2," ngrid=",i6)'

        if ngrid(ilevel-1,icpu) gt 0 then begin

            x=(*grid.level[ilevel-1,icpu]).xg(*,0)
            y=(*grid.level[ilevel-1,icpu]).xg(*,1)
            z=(*grid.level[ilevel-1,icpu]).xg(*,2)

            ind=where( x+dx gt xmin and x-dx lt xmax and $
                       y+dx gt ymin and y-dx lt ymax and $
                       z+dx gt zmin and z-dx lt zmax, nok)
            
            if(nok gt 0) then begin
                    
                x=x(ind)
                y=y(ind)
                z=z(ind)
                                
                for k=0,1 do begin
                    for j=0,1 do begin
                        for i=0,1 do begin

                            ind_cell=i+2*j+4*k
                            active=(*grid.level[ilevel-1,icpu]).son(ind,ind_cell)
                            if(ilevel lt levelup) then begin
                                ind2=where(active eq 0, nok2)
                            endif else begin
                                ind2=where(active gt -1, nok2)
                            endelse
                            
                            mesh=(*hydro.levelh[ilevel-1,icpu])
                            cxl=mesh.u(ind,ind_cell,4)
                            cyl=mesh.u(ind,ind_cell,5)
                            czl=mesh.u(ind,ind_cell,6)
                            cxr=mesh.u(ind,ind_cell,7)
                            cyr=mesh.u(ind,ind_cell,8)
                            czr=mesh.u(ind,ind_cell,9)
                            
                            if (nok2 gt 0) then begin
                            
                                nxc=2L^(levelup-ilevel)                                
                                for ii2=0L,n_elements(ind2)-1L do begin

                                    xl=x(ind2(ii2))+(double(i)-1.0)*dx
                                    yl=y(ind2(ii2))+(double(j)-1.0)*dx
                                    zl=z(ind2(ii2))+(double(k)-1.0)*dx

                                    il=fix((xl-xmin)/dxmin)
                                    jl=fix((yl-ymin)/dxmin)
                                    kl=fix((zl-zmin)/dxmin)

                                    if(il ge 0 and il lt nimx and $
                                       jl ge 0 and jl lt nimy and $
                                       kl ge 0 and kl lt nimz) then begin

                                        for ii=0L,nxc-1L do begin
                                            for jj=0L,nxc-1L do begin
                                                for kk=0L,nxc-1L do begin
                                                    xxl=double(ii)/double(nxc)
                                                    xxr=(double(ii)+1.0d0)/double(nxc)
                                                    yyl=double(jj)/double(nxc)
                                                    yyr=(double(jj)+1.0d0)/double(nxc)
                                                    zzl=double(kk)/double(nxc)
                                                    zzr=(double(kk)+1.0d0)/double(nxc)
                                                    
                                                    bxl(il+ii,jl+jj,kl+kk)=cxr(ind2(ii2))*xxl+cxl(ind2(ii2))*(1.0-xxl)
                                                    bxr(il+ii,jl+jj,kl+kk)=cxr(ind2(ii2))*xxr+cxl(ind2(ii2))*(1.0-xxr)
                                                    byl(il+ii,jl+jj,kl+kk)=cyr(ind2(ii2))*yyl+cyl(ind2(ii2))*(1.0-yyl)
                                                    byr(il+ii,jl+jj,kl+kk)=cyr(ind2(ii2))*yyr+cyl(ind2(ii2))*(1.0-yyr)
                                                    bzl(il+ii,jl+jj,kl+kk)=czr(ind2(ii2))*zzl+czl(ind2(ii2))*(1.0-zzl)
                                                    bzr(il+ii,jl+jj,kl+kk)=czr(ind2(ii2))*zzr+czl(ind2(ii2))*(1.0-zzr)
                                                    
                                                endfor
                                            endfor
                                        endfor

                                    endif

                                endfor
                                
                            endif
                        endfor
                    endfor
                endfor
            endif
        endif
    endfor
endfor

; Store in sructure
if not keyword_set(type) then type=0
cube={nx:fix(nimx), ny:fix(nimy), nz:fix(nimz), $
      xr:xrnew, yr:yrnew, zr:zrnew, $
      bxl:bxl, byl:byl, bzl:bzl, $
      bxr:bxr, byr:byr, bzr:bzr }


; Free memory
active=0
q=0
mesh=0
mesh2=0
save=0.
subcube=0.
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

return,cube
end
;###################################################
;###################################################
;###################################################
