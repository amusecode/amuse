;+
; NAME:
;       CUT3D
;
; PURPOSE:
;       This procedure extracts a slice from a 3D AMR data set and
;       converts it to a 2D AMR data set.
;
; CATEGORY:
;       Data analysis.
;
; CALLING SEQUENCE:
;       CUT3D, Grid3d, Hydro3d, Grid2d, Hydro2d, 
;                          XCUT = xcut, YCUT = ycut, ZCUT = zcut
;
; INPUTS
;       Grid3d: structure containing the 3D AMR mesh.
;
;       Hydro3d: structure containing the 3D hydro variables.
;
; OPTIONAL INPUTS:
;       XCUT:    if set, defines the X coordinate of the cutting
;       plane. 
;
;       YCUT:    if set, defines the Y coordinate of the cutting
;       plane. 
;
;       ZCUT:    if set, defines the Z coordinate of the cutting
;       plane. 
;
; OUTPUTS:
;       Grid2d: structure containing the 2D AMR mesh defining the cut.
;
;       Hydro2d: structure containing the 2D hydro variables in the cut.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To extract a slice at coordinate z = 16 in the 3D AMR mesh,
;       type: 
;
;               CUT3D, Grid3d, Hydro3d, Grid2d, Hydro2d, zcut=16.
;
;       To plot the density in the resulting plane as a color map:
;
;               TV2D, Grid2d, Hydro2d, type=1, levelmax=3
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro cut3d,grid ,hydro,grid1,hydro1 $
              ,xcut=xcut,ycut=ycut,zcut=zcut

IF N_PARAMS() NE 4 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'cutslice3d'
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

del_amr,grid1
del_hydro,hydro1

nncpu=ncpu
nndim=2
nnlevelmax=0
nngrid=lonarr(nlevelmax,ncpu)
nnlevel=PTRARR(nlevelmax,ncpu)
nnlevelh=PTRARR(nlevelmax,ncpu)
nnvar=nvar

; Set up maximum levels
leveldown=1
levelup=nlevelmax

; Check that a cut plane has been properly chosen
if not keyword_set(xcut) and $
   not keyword_set(ycut) and $
   not keyword_set(zcut) then begin
    print,'Choose a cut coordinate'
    return
endif

; Set geometry
if keyword_set(xcut) then begin
    xcut=MAX([xcut,1.d-8])
    xcut=MIN([xcut,1d0-1.d-8])
endif
if keyword_set(ycut) then begin
    ycut=MAX([ycut,1.d-8])
    ycut=MIN([ycut,1d0-1.d-8])
endif
if keyword_set(zcut) then begin
    zcut=MAX([zcut,1.d-8])
    zcut=MIN([zcut,1d0-1.d-8])
endif

; Define new 2D AMR and HYDRO structures

; Project fine levels
for ilevel=leveldown,levelup do begin
bool=0
for icpu=0,ncpu-1 do begin

    dx=0.5d0^ilevel
    if keyword_set(verbose) then print,ilevel,ngrid(ilevel-1,icpu)$
      ,format='("Level=",i2," ngrid=",i6)'

    if ngrid(ilevel-1,icpu) gt 0 then begin

    x=(*grid.level[ilevel-1,icpu]).xg(*,0)
    y=(*grid.level[ilevel-1,icpu]).xg(*,1)
    z=(*grid.level[ilevel-1,icpu]).xg(*,2)

    if keyword_set(xcut) then begin
        ind=where( x ge xcut-dx and x le xcut+dx, ngridok)
        if(ngridok gt 0) then begin
            mesh1=(*hydro.levelh[ilevel-1,icpu])
            mesh2=(*grid.level[ilevel-1,icpu])
            xg=x(ind)
            yg=y(ind)
            zg=z(ind)
            for i=0,1 do begin
                xc=xg+(double(i)-0.5d0)*dx
                indok=where(xc ge xcut-dx/2. and xc lt xcut+dx/2,nok)
                if(nok gt 0)then begin
                    if(bool eq 0) then begin
                        nnlevelmax=nnlevelmax+1
                        bool=1
                    endif
                    nngrid(ilevel-1,icpu)=nok
                    ; Define structure at fine level
                    mesh3={ilevel:ilevel,nc:nok $
                           ,xg:dblarr(nok,nndim), son:lonarr(nok,4)}
                    mesh3h={ilevel:ilevel,nc:nok,u:dblarr(nok,4,nvar)}
                    mesh3.xg(*,0)=yg(indok)
                    mesh3.xg(*,1)=zg(indok)
                    for j=0,1 do begin
                        for k=0,1 do begin
                            ind_cell=i+2*j+4*k
                            ind_cell2=j+2*k
                            xxx=mesh2.son(*,ind_cell)
                            yyy=xxx(ind)
                            mesh3.son(*,ind_cell2)=yyy(indok)
                            for ivar=0,nvar-1 do begin
                                xxx=mesh1.u(*,ind_cell,ivar)
                                yyy=xxx(ind)
                                mesh3h.u(*,ind_cell2,ivar)=yyy(indok)
                            endfor
                        endfor
                    endfor
                    pc=ptr_new(mesh3)
                    nnlevel(ilevel-1,icpu)=pc
                    pc=ptr_new(mesh3h)
                    nnlevelh(ilevel-1,icpu)=pc
                endif
            endfor
        endif
    endif

    if keyword_set(ycut) then begin
        ind=where( y ge ycut-dx and y le ycut+dx, ngridok)
        if(ngridok gt 0) then begin
            mesh1=(*hydro.levelh[ilevel-1,icpu])
            mesh2=(*grid.level[ilevel-1,icpu])
            xg=x(ind)
            yg=y(ind)
            zg=z(ind)            
            for j=0,1 do begin
                yc=yg+(double(j)-0.5d0)*dx
                indok=where(yc ge ycut-dx/2. and yc lt ycut+dx/2,nok)
                if(nok gt 0)then begin
                    if(bool eq 0) then begin
                        nnlevelmax=nnlevelmax+1
                        bool=1
                    endif
                    nngrid(ilevel-1,icpu)=nok
                    ; Define structure at fine level
                    mesh3={ilevel:ilevel,nc:nok $
                           ,xg:dblarr(nok,nndim), son:lonarr(nok,4)}
                    mesh3h={ilevel:ilevel,nc:nok,u:dblarr(nok,4,nvar)}
                    mesh3.xg(*,0)=xg(indok)
                    mesh3.xg(*,1)=zg(indok)
                    for i=0,1 do begin
                        for k=0,1 do begin
                            ind_cell=i+2*j+4*k
                            ind_cell2=i+2*k
                            xxx=mesh2.son(*,ind_cell)
                            yyy=xxx(ind)
                            mesh3.son(*,ind_cell2)=yyy(indok)
                            for ivar=0,nvar-1 do begin
                                xxx=mesh1.u(*,ind_cell,ivar)
                                yyy=xxx(ind)
                                mesh3h.u(*,ind_cell2,ivar)=yyy(indok)
                            endfor
                        endfor
                    endfor
                    pc=ptr_new(mesh3)
                    nnlevel(ilevel-1,icpu)=pc
                    pc=ptr_new(mesh3h)
                    nnlevelh(ilevel-1,icpu)=pc
                endif
            endfor
        endif
    endif
    
    if keyword_set(zcut) then begin
        ind=where( z ge zcut-dx and z le zcut+dx, ngridok)
        if(ngridok gt 0) then begin
            mesh1=(*hydro.levelh[ilevel-1,icpu])
            mesh2=(*grid.level[ilevel-1,icpu])
            xg=x(ind)
            yg=y(ind)
            zg=z(ind)            
            for k=0,1 do begin
                zc=zg+(double(k)-0.5d0)*dx
                indok=where(zc ge zcut-dx/2. and zc lt zcut+dx/2.,nok)
                if(nok gt 0)then begin
                    if(bool eq 0) then begin
                        nnlevelmax=nnlevelmax+1
                        bool=1
                    endif
                    nngrid(ilevel-1,icpu)=nok
                    ; Define structure at fine level
                    mesh3={ilevel:ilevel,nc:nok $
                           ,xg:dblarr(nok,nndim), son:lonarr(nok,4)}
                    mesh3h={ilevel:ilevel,nc:nok,u:dblarr(nok,4,nvar)}
                    mesh3.xg(*,0)=xg(indok)
                    mesh3.xg(*,1)=yg(indok)
                    for i=0,1 do begin
                        for j=0,1 do begin
                            ind_cell=i+2*j+4*k
                            ind_cell2=i+2*j
                            xxx=mesh2.son(*,ind_cell)
                            yyy=xxx(ind)
                            mesh3.son(*,ind_cell2)=yyy(indok)
                            for ivar=0,nvar-1 do begin
                                xxx=mesh1.u(*,ind_cell,ivar)
                                yyy=xxx(ind)
                                mesh3h.u(*,ind_cell2,ivar)=yyy(indok)
                            endfor
                        endfor
                    endfor
                    pc=ptr_new(mesh3)
                    nnlevel(ilevel-1,icpu)=pc
                    pc=ptr_new(mesh3h)
                    nnlevelh(ilevel-1,icpu)=pc
                endif
            endfor
        endif
    endif

    endif

endfor
endfor

grid1={ ncpu:nncpu,ndim:nndim,time:grid.time,aexp:grid.aexp $
       ,nlevelmax:nnlevelmax $
       ,boxlen:grid.boxlen,ngridtot:grid.ngridtot $
       ,ngrid:nngrid(0:nnlevelmax-1,0:ncpu-1) $
       ,level:nnlevel(0:nnlevelmax-1,0:ncpu-1) };$
;       ,unit_l:grid.unit_l $
;       ,unit_d:grid.unit_d $
;       ,unit_t:grid.unit_t }

hydro1={ncpu:nncpu,nlevelmax:nnlevelmax,gamma:hydro.gamma,ndim:nndim,nlevemax:nnlevelmax $
        ,nvar:nvar,levelh:nnlevelh(0:nnlevelmax-1,0:ncpu-1)}
return
end
;###################################################
;###################################################
;###################################################
