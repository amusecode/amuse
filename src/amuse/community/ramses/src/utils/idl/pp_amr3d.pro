;+
; NAME:
;	PP_AMR3D
;
; PURPOSE:
;	This procedure plots to the current device the mesh structure
;	stored in memory. Valid for 3D RAMSES simulations. Only octs
;	boundaries are shown.
;
; CATEGORY:
;	Plotting routines.
;
; CALLING SEQUENCE:
;       PP_AMR3D, Grid, SCALE=scale, AX=ax, AY=ay, AZ=az, 
;                       NOERASE=noerase, NODATA=nodata, 
;                       XR=xr, YR=yr, ZR=zr,
;                       X0=x0, Y0=y0, Z0=z0, 
;                       LEVEL=level, COLOR=color 
;
; INPUTS:
;       Grid = Structure defining the AMR grid.
;
; OPTIONAL INPUTS:
;       SCALE:   if set, defines the scale of the
;       plot. SCALE = 1. corresponds to the full box size, while
;       SCALE = 0.5 corresponds to half the box size and so on. 
;       Default: 1.
;
;       AX: Angle of rotation about the X axis in degrees. Default is
;       0. 
;
;       AY: Angle of rotation about the Y axis in degrees. Default is
;       0. 
;
;       AZ: Angle of rotation about the Z axis in degrees. Default is
;       0. 
;
;       PERSP: Position of the eye along the Z-axis for perspective
;       view. Default: no perspective.
;
;       NOERASE: standard IDL option.
;
;       NODATA:  standard IDL option.
;
;       NOAXIS:  If set, X, Y and Z axis are not shown. Otherwise, X
;       axis is shown in red, Y axis in blue and Z axis in green.
;
;       XR:      if set, the routine plots only octs that lie within
;       these boundaries, defined along the X axis. Default: the
;       whole box ! 
;
;       YR:      Same for the Y axis.
;
;       ZR:      Same for the Z axis.
;
;       X0:      Move the center to this new coordinate. Default: the
;       center of the chosen X range.
;       
;       Y0:      Move the center to this new coordinate. Default: the
;       center of the chosen Y range.
;      
;       Z0:      Move the center to this new coordinate. Default: the
;       center of the chosen Z range.
;      
;       LEVEL:     If set, the routine plots only octs belonging to
;       level LEVEL. Otherwise, all levels are shown (default).
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
;       To plot the mesh structure currently stored in variable GRID,
;       with the computational box head on, type:
;
;	        PP_AMR3D, GRID, /COLOR, ax=30, ay=30
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro plot_cube,x,y,z,color=color
ss=size(x)
if(ss(0) eq 1)then begin
    plots,x(0),y(0),z(0),/t3d          ,color=color
    plots,x(1),y(1),z(1),/t3d,/continue,color=color
    plots,x(3),y(3),z(3),/t3d,/continue,color=color
    plots,x(2),y(2),z(2),/t3d,/continue,color=color
    plots,x(0),y(0),z(0),/t3d,/continue,color=color
    plots,x(0),y(0),z(0),/t3d          ,color=color
    plots,x(4),y(4),z(4),/t3d,/continue,color=color
    plots,x(1),y(1),z(1),/t3d          ,color=color
    plots,x(5),y(5),z(5),/t3d,/continue,color=color
    plots,x(2),y(2),z(2),/t3d          ,color=color
    plots,x(6),y(6),z(6),/t3d,/continue,color=color
    plots,x(3),y(3),z(3),/t3d          ,color=color
    plots,x(7),y(7),z(7),/t3d,/continue,color=color
    plots,x(4),y(4),z(4),/t3d          ,color=color
    plots,x(5),y(5),z(5),/t3d,/continue,color=color
    plots,x(7),y(7),z(7),/t3d,/continue,color=color
    plots,x(6),y(6),z(6),/t3d,/continue,color=color
    plots,x(4),y(4),z(4),/t3d,/continue,color=color
endif

if (ss(0) eq 2L) then begin
    n=ss(2)
    for i=0L,n-1L do begin
        plots,[x(0,i),x(1,i),x(3,i),x(2,i),x(0,i)],$
              [y(0,i),y(1,i),y(3,i),y(2,i),y(0,i)],$
              [z(0,i),z(1,i),z(3,i),z(2,i),z(0,i)],/t3d,color=color
        plots,[x(0,i),x(4,i)],[y(0,i),y(4,i)],[z(0,i),z(4,i)],/t3d,color=color
        plots,[x(1,i),x(5,i)],[y(1,i),y(5,i)],[z(1,i),z(5,i)],/t3d,color=color
        plots,[x(2,i),x(6,i)],[y(2,i),y(6,i)],[z(2,i),z(6,i)],/t3d,color=color
        plots,[x(3,i),x(7,i)],[y(3,i),y(7,i)],[z(3,i),z(7,i)],/t3d,color=color
        plots,[x(4,i),x(5,i),x(7,i),x(6,i),x(4,i)],$
              [y(4,i),y(5,i),y(7,i),y(6,i),y(4,i)],$
              [z(4,i),z(5,i),z(7,i),z(6,i),z(4,i)],/t3d,color=color
    endfor
endif
end
;###################################################
;###################################################
;###################################################
pro pp_amr3d,grid,scale=scale,ax=ax,ay=ay,az=az,persp=persp $
                 ,noerase=noerase,nodata=nodata $
                 ,xr=xr,yr=yr,zr=zr $
                 ,x0=x0,y0=y0,z0=z0 $
                 ,cmin=cmin,cmax=cmax $
                 ,lmin=lmin,lmax=lmax,level=level $
                 ,color=color,noaxis=noaxis,cpu=cpu

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'pp_amr3d'
    RETURN
ENDIF

if grid.ndim ne 3 then begin
    print,'Mesh should have 3 dimensions'
    print,'but ndim=',grid.ndim
    return
endif

ncpu=grid.ncpu

tek_color

if not keyword_set(scale) then scale = 1.0
if not keyword_set(xr)    then xr=[0,1.]
if not keyword_set(yr)    then yr=[0,1.]
if not keyword_set(zr)    then zr=[0,1.]
if not keyword_set(x0)    then x0=0.5*(xr[0]+xr[1])
if not keyword_set(y0)    then y0=0.5*(yr[0]+yr[1])
if not keyword_set(z0)    then z0=0.5*(zr[0]+zr[1])

xz=[x0,x0]+scale*[-0.5,0.5]
yz=[y0,y0]+scale*[-0.5,0.5]
zz=[z0,z0]+scale*[-0.5,0.5]

if not keyword_set(noerase) then erase

!p.t3d=1
CREATE_VIEW $
  ,XMIN=xz[0],XMAX=xz[1] $
  ,YMIN=yz[0],YMAX=yz[1] $
  ,ZMIN=zz[0],ZMAX=zz[1] $
  ,AX=ax,AY=ay,AZ=az $
  ,WINX=!d.x_size,WINY=!d.y_size,PERSP=PERSP

leveldown=1
levelup=grid.nlevelmax
if keyword_set(lmin) then leveldown=lmin
if keyword_set(lmax) then levelup  =lmax
if not keyword_set(cmin) then cmin=1
if not keyword_set(cmax) then cmax=ncpu

if keyword_set(nodata) then goto,skip

; Loop over levels
for ilevel=leveldown,levelup do begin
; Mesh size at current level
    dx=0.5d0^(ilevel-1)
    dy=dx
    dz=dx

    for icpu=cmin-1,cmax-1 do begin
    if(grid.ngrid(ilevel-1,icpu) gt 0)then begin
        
        mesh2=(*grid.level[ilevel-1,icpu])
        x=mesh2.xg(*,0)
        y=mesh2.xg(*,1)
        z=mesh2.xg(*,2)
        ind=where(x-0.5*dx ge xr[0] and x+0.5*dx le xr[1] and $
                  y-0.5*dy ge yr[0] and y+0.5*dy le yr[1] and $
                  z-0.5*dz ge zr[0] and z+0.5*dz le zr[1], nplot)        

        if nplot gt 0 then begin
            xf=dblarr(8,grid.ngrid(ilevel-1,icpu))
            yf=dblarr(8,grid.ngrid(ilevel-1,icpu))
            zf=dblarr(8,grid.ngrid(ilevel-1,icpu))
            for i=0,1 do begin
                for j=0,1 do begin
                    for k=0,1 do begin
                        ind1=i+2*j+4*k
                        xf(ind1,ind)=x(ind)+(i-0.5)*dx
                        yf(ind1,ind)=y(ind)+(j-0.5)*dy
                        zf(ind1,ind)=z(ind)+(k-0.5)*dz
                    endfor
                endfor
            endfor
            if keyword_set(color) then begin
                if keyword_set(cpu) then begin
                    plot_cube,xf(*,ind),yf(*,ind),zf(*,ind),color=icpu+1
                endif else begin
                    if keyword_set(lmin)then begin
                        plot_cube,xf(*,ind),yf(*,ind),zf(*,ind),color=ilevel-lmin+1
                    endif else begin
                        plot_cube,xf(*,ind),yf(*,ind),zf(*,ind),color=ilevel
                    endelse
                endelse
            endif else begin
                plot_cube,xf(*,ind),yf(*,ind),zf(*,ind),color=1
            endelse
        endif

    endif
    endfor
endfor

skip:
; Draw x, y and z-range cube boundaries
plots,xr[0],yr[0],zr[0],/t3d          ,color=1
plots,xr[1],yr[0],zr[0],/t3d,/continue,color=1
plots,xr[1],yr[1],zr[0],/t3d,/continue,color=1
plots,xr[0],yr[1],zr[0],/t3d,/continue,color=1
plots,xr[0],yr[0],zr[0],/t3d,/continue,color=1
plots,xr[0],yr[0],zr[1],/t3d          ,color=1
plots,xr[1],yr[0],zr[1],/t3d,/continue,color=1
plots,xr[1],yr[1],zr[1],/t3d,/continue,color=1
plots,xr[0],yr[1],zr[1],/t3d,/continue,color=1
plots,xr[0],yr[0],zr[1],/t3d,/continue,color=1
plots,xr[0],yr[0],zr[0],/t3d          ,color=1
plots,xr[0],yr[0],zr[1],/t3d,/continue,color=1
plots,xr[1],yr[0],zr[0],/t3d          ,color=1
plots,xr[1],yr[0],zr[1],/t3d,/continue,color=1
plots,xr[0],yr[1],zr[0],/t3d          ,color=1
plots,xr[0],yr[1],zr[1],/t3d,/continue,color=1
plots,xr[1],yr[1],zr[0],/t3d          ,color=1
plots,xr[1],yr[1],zr[1],/t3d,/continue,color=1

; Draw 3 axis
if not keyword_set(noaxis) then begin
plots,xr[0],yr[0],zr[0],/t3d
plots,xr[1],yr[0],zr[0],/t3d,/continue,color=2,thick=2
plots,xr[0],yr[0],zr[0],/t3d
plots,xr[0],yr[1],zr[0],/t3d,/continue,color=3,thick=2
plots,xr[0],yr[0],zr[0],/t3d
plots,xr[0],yr[0],zr[1],/t3d,/continue,color=4,thick=2
endif

; Free memory
ind=0
mesh2=0
x=0
xf=0
y=0
yf=0
z=0
zf=0

end
;###################################################
;###################################################
;###################################################

