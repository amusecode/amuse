;+
; NAME:
;	RD_GRAV
;
; PURPOSE:
;	This procedure reads particles from a RAMSES GRAV file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_GRAV, Grav, FILE=file, SWAP=swap, NVAR=nvar, NCPU=ncpu,
;	ICPU=icpu, VERBOSE=verbose
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;	FILE:   if set, input the scalar string containing the name of
;	        the file to be read. Otherwise, a PICKFILE widget is
;	        launched.  
;
;       SWAP:   if set, reverse the bit ordering (Little Endian versus
;               Big Endian)
;
;       IVAR:   first grav variable to be read. Default: 1.
;
;       NVAR:   if set, load in memory only NVAR variables. Default:
;               all variables from ivar to nvar_max. 
;	
;       ICPU:   first cpu file to be read. Default: 1.
;
;       NCPU:   number of cpu files to read, starting from
;               icpu. Default: all files from icpu to ncpu_max.
;
; OUTPUTS:
;	Grav:  Store gravity vector in structure Grav.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; EXAMPLE:
;       To read on a SGI architecture a RAMSES GRAV file created on a
;       COMPAQ Workstation, type:
;
;	        RD_GRAV,h,file='grav_00001.out',/swap
;
;       If the file was generated on the same IEEE system, just type:
;
;               RD_GRAV,h,file='grav_00001.out'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
pro rd_grav, grav, file=file, swap=swap, ivar=ivar, icpu=icpu, verbose=verbose, nout=nout

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_grav'
    RETURN
ENDIF

if not keyword_set(icpu) then icpu=0
if not keyword_set(ivar) then ivar=0
if icpu eq 0 then jcpu=1 else jcpu=icpu

suffix=getcarnum(jcpu)
if not keyword_set(file) and not keyword_set(nout) then begin
    key='*grav*.out'+suffix(jcpu-1)
    file=DIALOG_PICKFILE(/READ,filter=key)
endif
if keyword_set(nout) then begin
    suffnout=getcarnum(nout)
    file='output_'+suffnout(nout-1)+'/grav_'+suffnout(nout-1)+'.out'
endif
if not keyword_set(file) then return
base_offset=strpos(file,'.out')+4
file_base=strmid(file,0,base_offset)

; Freeing memory associated to grav structure
del_hydro,grav

; Initialize header variables
ncpu_run=0L & ndim=0L & nlevelmax=0L & ng=0L & nvar_run=0L & ilevel=0L & nboundary=0L
gamma=0d0

; Read first file to get header
print,'Reading file ',trim(file_base)
file=trim(file_base+suffix(jcpu-1))
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ncpu_run
readu,1,ndim
readu,1,nlevelmax
readu,1,nboundary
close,1

; Write header to screen
print,'ncpu      =',ncpu_run
print,'ndim      =',ndim
print,'nlevelmax =',nlevelmax
if nboundary eq 0 then begin
    print,"Periodic boundary conditions"
endif

; Allocate arrays
ncpu=ncpu_run & nvar=ndim
icpumin=1L & icpumax=ncpu & listmax=ncpu & vartot=nvar
if icpu gt 0 then begin
    icpumin=icpu & icpumax=icpu & listmax=ncpu+nboundary
endif
if ivar gt 0 then begin
    vartot=1
endif

suffix=getcarnum(ncpu_run)
ncell=2L^ndim
levelh=PTRARR(nlevelmax,listmax)

; Loop over cpu files
for jcpu=icpumin,icpumax do begin
    file=trim(file_base+suffix(jcpu-1))
    if keyword_set(verbose) then print,'Reading file ',file
    openr,1,file,/f77_unformatted,swap_endian=swap

; Read header
    readu,1,cpu_run
    readu,1,ndim
    readu,1,nlevelmax
    readu,1,nboundary

; Read fine levels
    nlevel=0L & ilevel=0L & ng=0L
    kcpumin=1L & kcpumax=nboundary+ncpu
    for i=0L,nlevelmax-1L do begin
        for kcpu=kcpumin,kcpumax do begin
            readu,1,ilevel
            readu,1,ng
            if(ng gt 0L)then begin
                if keyword_set(verbose) then print,ilevel,ng $
                  ,format='("Level ",i2," has ",i6," grids")'
                nlevel=nlevel+1L
                mesh={ilevel:ilevel,nc:ng,u:dblarr(ng,ncell,vartot) }
                xx=dblarr(ng)
                for icell=0,ncell-1 do begin
                    for jvar=1L,ndim do begin
                        readu,1,xx
                        if ivar eq 0 then begin
                            mesh.u(*,icell,jvar-1)=xx
                        endif else begin
                            if jvar eq ivar then begin
                                mesh.u(*,icell,0)=xx
                            endif
                        endelse
                    endfor
                endfor
                if icpu eq 0 then begin
                    if(kcpu eq jcpu)then begin
                        pc=ptr_new(mesh)
                        levelh(i,jcpu-1)=pc
                    endif
                endif else begin
                    pc=ptr_new(mesh)
                    levelh(i,kcpu-1)=pc
                endelse
                xx=0.
            endif
        endfor
    endfor
    close,1
endfor

;levelh=levelh[0:nlevelmax-1,0:ncpu-1]

grav={ncpu:listmax,gamma:1.0,ndim:ndim,nlevelmax:nlevelmax,nvar:vartot,levelh:levelh}

end
;###################################################
;###################################################
;###################################################
