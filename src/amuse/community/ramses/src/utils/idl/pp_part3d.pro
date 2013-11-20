;+
; NAME:
;	PP_PART3D
;
; PURPOSE:
;	This procedure plots particles to the current device.
;	It uses the 3D view defined by PP_AMR3D. Particle colors 
;       are coded using particle level (default) or particle density.
;	The last option must be used after routine MAKE_DENS has been
;	used, and particle density array has been stored to a file.
;
; CATEGORY:
;	Plotting routines.
;
; CALLING SEQUENCE:
;       PP_PART3D, Part, XR=xr, YR=yr, ZR=zr, DENSITY=density
;
; INPUTS:
;       Part = Structure defining particles.
;
; OPTIONAL INPUTS:
;       XR:      if set, the routine plots only particle that lie
;       within these boundaries, defined along the X axis. Default:
;       the whole box ! 
;
;       YR:      Same for the Y axis.
;
;       ZR:      Same for the Z axis.
;
;       DENSITY: If set, the color coding corresponds to particle
;       local density. In this case, DENSITY should be set to a
;       filename corresponding to the particle density array created
;       by routine MAKE_DENS. If not set, the color coding corresponds
;       to the particle level (default).
;
;       MASS: If set, the color coding corresponds to particle
;       mass. If not set, the color coding corresponds
;       to the particle level (default).
;
; OUTPUTS:
;       None.
;       
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To plot particles with a color coding corresponding to the
;       local density, type:
;
;               ; Set the 3D view
;	        PP_AMR3D, GRID, /COLOR, ax=90, az=0, /NODATA
;
;               ; Plot particles
;               PP_PART3D, PART, DENSITY='density_file.dat' 
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro pp_part3d, part, xr=xr, yr=yr, zr=zr $
               ,density=density, dmin=dmin, dmax=dmax, swap=swap $
               ,mass=mass, time=time, ital=ital

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'pp_part3d'
    RETURN
ENDIF

if not keyword_set(xr) then xr=[min(part.xp(*,0)),max(part.xp(*,0))]
if not keyword_set(yr) then yr=[min(part.xp(*,1)),max(part.xp(*,1))]
if not keyword_set(zr) then zr=[min(part.xp(*,2)),max(part.xp(*,2))]

if !p.t3d eq 0 then begin
    !p.t3d=1
    CREATE_VIEW $
      ,XMIN=xr[0],XMAX=xr[1] $
      ,YMIN=yr[0],YMAX=yr[1] $
      ,ZMIN=zr[0],ZMAX=zr[1] $
      ,AX=0,AY=0,AZ=0 $
      ,WINX=!d.x_size,WINY=!d.y_size
endif

ind=where(   part.xp(*,0) ge min(xr) and part.xp(*,0) le max(xr) $
         and part.xp(*,1) ge min(yr) and part.xp(*,1) le max(yr) $
         and part.xp(*,2) ge min(zr) and part.xp(*,2) le max(zr) )
npart=n_elements(ind)
print,'npart=',npart
if npart eq 0L then return

x=part.xp(ind,0)
y=part.xp(ind,1)
z=part.xp(ind,2)
m=part.mp(ind)

if keyword_set(time) then begin
    d=part.ap(ind)
    ind=where(d gt 0.,nstar)
    print,'nstar=',nstar
    if nstar le 0 then return
    d=d(ind)
    x=x(ind)
    y=y(ind)
    z=z(ind)
    npart=nstar
    print,'Min and max age=',min(d),max(d)
    print,'sorting...'
    if keyword_set(dmin) then begin
        ind2=where(d ge dmin)
        d=d(ind2)
        x=x(ind2)
        y=y(ind2)
        z=z(ind2)
        npart=n_elements(ind2)
        ind2=0
    endif
    ind_col=sort(d)
    if not keyword_set(dmin) then dmin = min(d)
    if not keyword_set(dmax) then dmax = max(d)
    print,'New min and new max=',dmin,dmax
    color=15+BYTSCL(alog10(d),MIN=alog10(dmin),MAX=alog10(dmax) $
                    ,TOP=!d.TABLE_SIZE-16)
    psym_col=3
endif else if keyword_set(density) then begin
    d=part.mp(ind)
    print,'Min and max density=',min(d),max(d)
    print,'sorting...'
    if keyword_set(dmin) then begin
        ind2=where(d ge dmin)
        d=d(ind2)
        x=x(ind2)
        y=y(ind2)
        z=z(ind2)
        npart=n_elements(ind2)
        ind2=0
    endif
    ind_col=sort(d)
    if not keyword_set(dmin) then dmin = min(d)
    if not keyword_set(dmax) then dmax = max(d)
    print,'New min and new max=',dmin,dmax
    color=15+BYTSCL(alog10(d),MIN=alog10(dmin),MAX=alog10(dmax) $
                    ,TOP=!d.TABLE_SIZE-16)
    psym_col=3
endif else if keyword_set(mass) then begin
    m_max=max(part.mp)
    m_min=min(part.mp)
    print,m_min,m_max
    color=3+0.*m
    ind=where(m eq m_min)
    if ind(0) ne -1 then color(ind)=2
    ind=where(m eq m_max)
    if ind(0) ne -1 then color(ind)=4
    ind_col=sort(color)
    print,minmax(color)
    psym_col=6
endif else if keyword_set(ital) then begin
    color=fix(0.*m)
    ind=where(m gt max(m)/10.)
    color(ind)=2
    ind=where(m lt max(m)/10.)
    color(ind)=1
    ind_col=sort(color)
    ind_col=reverse(ind_col)
    psym_col=3
endif else begin
    psym_col=3
    print,'sorting...'
    ind_col=LINDGEN(npart)
    color=INTARR(npart)
    ind_zero=where(color eq 0, n_zero)
    if n_zero gt 0 then color(ind_zero)=1
endelse

for i=0L,npart-1L do begin
    plots,x(ind_col(i)),y(ind_col(i)),z(ind_col(i)) $
;      ,/t3d,color=1,psym=3;color(ind_col(i)),psym=psym_col
      ,/t3d,color=color(ind_col(i)),psym=psym_col
endfor

; Free memory
x=0.
y=0.
z=0.
l=0
d=0.
ind_col=0
color=0
ind_lev=0

end
