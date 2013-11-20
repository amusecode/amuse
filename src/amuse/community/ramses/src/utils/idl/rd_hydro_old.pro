;+
; NAME:
;	RD_HYDRO
;
; PURPOSE:
;	This procedure reads particles from a RAMSES HYDRO file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_HYDRO, Hydro, FILE=file, SWAP=swap, NVAR=nvar, NCPU=ncpu,
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
;       IVAR:   first hydro variable to be read. Default: 1.
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
;	Hydro:  Store hydrodynamics variables in structure Hydro.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; EXAMPLE:
;       To read on a SGI architecture a RAMSES HYDRO file created on a
;       COMPAQ Workstation, type:
;
;	        RD_HYDRO,h,file='hydro_00001.out',/swap
;
;       If the file was generated on the same IEEE system, just type:
;
;               RD_HYDRO,h,file='hydro_00001.out'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
pro rd_hydro, hydro, file=file, swap=swap, ivar=ivar, nvar=nvar $
              , icpu=icpu, ncpu=ncpu, verbose=verbose, nout=nout

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_hydro'
    RETURN
ENDIF

if not keyword_set(icpu) then icpu=1
if not keyword_set(ivar) then ivar=1

suffix=getcarnum(icpu)
if not keyword_set(file) and not keyword_set(nout) then begin
    key='*hydro*.out'+suffix(icpu-1)
    file=DIALOG_PICKFILE(/READ,filter=key)
endif
if keyword_set(nout) then begin
    suffnout=getcarnum(nout)
    file='output_'+suffnout(nout-1)+'/hydro_'+suffnout(nout-1)+'.out'
endif
if not keyword_set(file) then return
base_offset=strpos(file,'.out')+4
file_base=strmid(file,0,base_offset)

; Freeing memory associated to hydro structure
del_hydro,hydro

ncpu_run=0L & ndim=0L & nlevelmax=0L & ng=0L & nvar_run=0L & ilevel=0L

print,'Reading file ',trim(file_base)
file=trim(file_base+suffix(icpu-1))
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ncpu_run
readu,1,nvar_run
readu,1,ndim
readu,1,nlevelmax
readu,1,gamma
close,1

print,'ncpu      =',ncpu_run
print,'ndim      =',ndim
print,'nlevelmax =',nlevelmax
print,'nvar      =',nvar_run
print,'gamma     =',gamma


if not keyword_set(ncpu) then ncpu=ncpu_run-icpu+1
if not keyword_set(nvar) then nvar=nvar_run-ivar+1 

suffix=getcarnum(ncpu_run)

ncell=2L^ndim
levelh=PTRARR(nlevelmax,ncpu)

;Read levels
for jcpu=0,ncpu-1 do begin
    file=trim(file_base+suffix(jcpu+icpu-1))
    if keyword_set(verbose) then print,'Reading file ',file
    openr,1,file,/f77_unformatted,swap_endian=swap
    readu,1,cpu_run
    readu,1,nvar_run
    readu,1,ndim
    readu,1,nlevelmax
    readu,1,gamma
    nlevel=0L & ilevel=0L & ng=0L
    for i=0L,nlevelmax-1L do begin
        readu,1,ilevel
        readu,1,ng
        if(ng gt 0L)then begin
            if keyword_set(verbose) then print,ilevel,ng $
              ,format='("Level ",i2," has ",i6," grids")'
            nlevel=nlevel+1L
            mesh={ilevel:ilevel,nc:ng,u:fltarr(ng,ncell,nvar) }
            xx=fltarr(ng)
            for icell=0,ncell-1 do begin
                for jvar=0L,nvar_run-1L do begin
                    readu,1,xx
                    if jvar ge (ivar-1) and $
                       jvar lt (ivar+nvar-1) then begin
                        mesh.u(*,icell,jvar-ivar+1)=xx
                    endif
                endfor
            endfor
            pc=ptr_new(mesh)
            levelh(i,jcpu)=pc
            xx=0.
        endif
    endfor
    close,1
endfor

;levelh=levelh[0:nlevelmax-1,0:ncpu-1]

hydro={ncpu:ncpu,gamma:gamma,ndim:ndim,nlevelmax:nlevelmax,nvar:nvar,levelh:levelh}

end
;###################################################
;###################################################
;###################################################
