;+
; NAME:
;       SUB3D
;
; PURPOSE:
;       This procedure extract all data from AMR, HYDRO and PART structures
;       as new AMR, HYDRO and PART structures.
;
; CATEGORY:
;       Data analysis.
;
; CALLING SEQUENCE:
;       SUB3D, Grid, Hydro, Part, Grid_new, Hydro_new, Part_new, 
;                          XR = xr, YR = yr, ZR = zr, 
;                          VERBOSE = verbose
;
; INPUTS
;       Grid: structure containing the AMR mesh.
;
;       Hydro: structure containing the hydro variables.
;
;       Part: structure containing the particle variables.
;
; OPTIONAL INPUTS:
;       XR:      if set, defines the map boundaries for the X
;       axis. Default: the whole box.
;
;       YR:      same for the Y axis.
;
;       ZR:      same for the Z axis.
;
;       VERBOSE: if set, output infos to screen. Default: not set.
;
; OUTPUTS:
;       Grid_new: structure containing the new AMR mesh.
;
;       Hydro_new: structure containing the new hydro variables.
;
;       Part_new: structure containing the new particle variables.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To extract a sub-cube, type:
;
;               SUB3D, grid, hydro, part, grid1, hydro1, part1, 
;               xr=[0,0.2],yr=[0,0.2], zr=[0,0.3]
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
pro sub3d,grid,hydro,part,grid1,hydro1,part1 $
         ,xr=xr,yr=yr,zr=zr,verbose=verbose

IF N_PARAMS() NE 6 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'subcube3d'
    RETURN
ENDIF

ncpu=grid.ncpu
ndim=grid.ndim
nlevelmax=grid.nlevelmax
ngrid=grid.ngrid
nvar=hydro.nvar

if ndim ne 3 then begin
    print,'Mesh should have 3 dimensions'
    print,'but ndim=',ndim
    return
endif

; Free memory associated to grid1, part1 and hydro1
del_amr,grid1
del_hydro,hydro1
part1=0.

nncpu=ncpu
nndim=3
nnlevelmax=0
nngridtot=0L
nngrid=lonarr(nlevelmax,ncpu)
nnlevel=PTRARR(nlevelmax,ncpu)
nnlevelh=PTRARR(nlevelmax,ncpu)
nnvar=nvar

; Set up geometry and scale
if not keyword_set(xr) then xr=[0.,1.]
if not keyword_set(yr) then yr=[0.,1.]
if not keyword_set(zr) then zr=[0.,1.]
xmin=MAX([0.,xr(0)]) & xmax=MIN([1.,xr(1)])
ymin=MAX([0.,yr(0)]) & ymax=MIN([1.,yr(1)])
zmin=MAX([0.,zr(0)]) & zmax=MIN([1.,zr(1)])

;=====================
; Extract fine levels
;=====================
for ilevel=1,nlevelmax do begin
    bool=0
    for icpu=0,ncpu-1 do begin
        
        dx=0.5d0^ilevel
        if keyword_set(verbose) then $
          print,ilevel,ngrid(ilevel-1,icpu),format='("Level=",i2," ngrid=",i6)'
        
        if ngrid(ilevel-1,icpu) gt 0 then begin
            
            xg=(*grid.level[ilevel-1,icpu]).xg(*,0)
            yg=(*grid.level[ilevel-1,icpu]).xg(*,1)
            zg=(*grid.level[ilevel-1,icpu]).xg(*,2)
            
            ind=where( xg gt xmin-dx and xg lt xmax+dx and $
                       yg gt ymin-dx and yg lt ymax+dx and $
                       zg gt zmin-dx and zg lt zmax+dx, nok)
            
            if(nok gt 0) then begin
; Define structure at fine level
                nngridtot=nngridtot+nok
                xg=xg(ind)
                yg=yg(ind)
                zg=zg(ind)
                if(bool eq 0) then begin
                    nnlevelmax=nnlevelmax+1
                    bool=1
                endif
                nngrid(ilevel-1,icpu)=nok
                mesh3 ={ilevel:ilevel,nc:nok,xg:fltarr(nok,nndim) $
                        ,son:lonarr(nok,8)}
                mesh3h={ilevel:ilevel,nc:nok,u:fltarr(nok,8,nvar)}
                mesh3.xg(*,0)=xg
                mesh3.xg(*,1)=yg
                mesh3.xg(*,2)=zg
                for i=0,1 do begin
                for j=0,1 do begin
                for k=0,1 do begin
                    ind_cell=i+2*j+4*k
                    mesh3.son(*,ind_cell) $
                      =(*grid.level[ilevel-1,icpu]).son(ind,ind_cell)
                    for ivar=0,nvar-1 do begin
                        mesh3h.u(*,ind_cell,ivar) $
                          =(*hydro.levelh[ilevel-1,icpu]).u(ind,ind_cell,ivar)
                    endfor
                endfor
                endfor
                endfor
                pc=ptr_new(mesh3)
                nnlevel(ilevel-1,icpu)=pc
                pc=ptr_new(mesh3h)
                nnlevelh(ilevel-1,icpu)=pc
            endif
        endif
    endfor
endfor
                        
grid1={ncpu:nncpu,ndim:nndim,aexp:grid.aexp,nlevelmax:nnlevelmax $
       ,boxlen:grid.boxlen,ngridtot:nngridtot $
       ,ngrid:nngrid(0:nnlevelmax-1,0:ncpu-1) $
       ,level:nnlevel(0:nnlevelmax-1,0:ncpu-1)}

hydro1={gamma:hydro.gamma,ndim:nndim,nlevemax:nnlevelmax $
        ,nvar:nvar,levelh:nnlevelh(0:nnlevelmax-1,0:ncpu-1)}

;==================
; Extract particles
;==================
ind=where(   part.xp(*,0) ge min(xr) and part.xp(*,0) le max(xr) $
         and part.xp(*,1) ge min(yr) and part.xp(*,1) le max(yr) $
         and part.xp(*,2) ge min(zr) and part.xp(*,2) le max(zr),npart1 )
if keyword_set(verbose) then print,'npart=',npart1
if npart1 eq 0L then begin
    ndim1=part.ndim
    part1={ndim:ndim1,npart:0L}
endif else begin
    ndim1=part.ndim
    xp=fltarr(npart1,ndim1)
    vp=fltarr(npart1,ndim1)
    mp=fltarr(npart1)
    lp=lonarr(npart1)
    xp(*,0)=part.xp(ind,0)
    xp(*,1)=part.xp(ind,1)
    xp(*,2)=part.xp(ind,2)
    vp(*,0)=part.vp(ind,0)
    vp(*,1)=part.vp(ind,1)
    vp(*,2)=part.vp(ind,2)
    lp=part.lp(ind)
    mp=part.mp(ind)
    part1={ndim:ndim1,npart:npart1,xp:xp,vp:vp,mp:mp,lp:lp}
    ok_dens=strpos(strjoin(tag_names(part),' '),'DP') ge 0
    if ok_dens then begin
        dp=part.dp(ind)
        part1=create_struct(part1,'dp',dp)
    endif
    ok_time=strpos(strjoin(tag_names(part),' '),'AP') ge 0
    if ok_time then begin
        ap=part.ap(ind)
        part1=create_struct(part1,'ap',ap)
    endif
endelse

return
end
;###################################################
;###################################################
;###################################################
