;+
; NAME:
;	RAMSES2TIPSY
;
; PURPOSE:
;	This procedure reads particles from a RAMSES PART file and
;	create a TIPSY ASCII particles file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RAMSES2TIPSY, FILE=file, SWAP=swap, TIPSY=tipsyfile
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;	FILE:    if set, input the scalar string containing the name of
;         	 the file to be read. Otherwise, a PICKFILE widget is
;         	 launched.  
;
;       SWAP:    if set, reverse the bit ordering (Little Endian versus
;                BigEndian)
;
;	TIPSYFILE:    if set, input the scalar string containing the name of
;         	 the tipsy file to be written. Otherwise, the name of
;         	 the tipsy file will be tipsy.ascii.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To read on a SGI architecture a RAMSES PART file created on a
;       COMPAQ Workstation, type:
;
;	        RAMSES2TIPSY, file='part_00001.out',/swap
;
;       If the file was generated on the same IEEE system, just type:
;
;               RAMSES2TPISY, file='part_00001.out'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 19/09/2006.
;                       e-mail: Romain.Teyssier@cea.fr
;-
pro RAMSES2TIPSY, part, file=file, swap=swap, tipsy=tipsy $
             ,ncpu=ncpu, icpu=icpu , verbose=verbose, nout=nout  $
             ,star=star, time=time

if not keyword_set(time) then time=1.0d0

if not keyword_set(icpu) then icpu=1

suffix=getcarnum(icpu)

if not keyword_set(file) and not keyword_set(nout) then begin
    key='*.out'+suffix(icpu-1)
    file=DIALOG_PICKFILE(/READ,filter=key)
endif

if not keyword_set(tipsy) then tipsy='tipsy.ascii'

if keyword_set(nout) then begin
    suffnout=getcarnum(nout)
    if keyword_set(star) then begin
        file='output_'+suffnout(nout-1)+'/star_'+suffnout(nout-1)+'.out'
    endif else begin
        file='output_'+suffnout(nout-1)+'/part_'+suffnout(nout-1)+'.out'
    endelse
endif
if not keyword_set(file) then return
base_offset=strpos(file,'.out')+4
file_base=strmid(file,0,base_offset)

ncpu_run=0L & ndim=0L

print,'Reading file ',trim(file_base)
file=trim(file_base+suffix(icpu-1))
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ncpu_run
readu,1,ndim
close,1

print,'ncpu      =',ncpu_run
print,'ndim      =',ndim

if not keyword_set(ncpu) then ncpu=ncpu_run-icpu+1

suffix=getcarnum(icpu+ncpu-1)

npartp=0L & npart=0L & ndum=0L

; Compute total number of particle
for jcpu=0,ncpu-1 do begin
    file=trim(file_base+suffix(jcpu+icpu-1))
    if keyword_set(verbose) then print,'Reading file ',trim(file)
    openr,1,file,/f77_unformatted,swap_endian=swap
    readu,1,ncpu_run
    readu,1,ndim
    readu,1,npartp
    close,1
    if keyword_set(verbose) then print,ndim,npartp $
      ,format='("ndim=",I1," npart=",I8)'
    npart=npart+npartp
endfor
print,'npart     =',npart

; Allocate memory
part={  xp:fltarr(npart,ndim) $
       ,vp:fltarr(npart,ndim) $
       ,id:lonarr(npart) $
       ,mp:fltarr(npart)}

iskip=0L
for jcpu=0,ncpu-1 do begin
    file=trim(file_base+suffix(jcpu+icpu-1))
    openr,1,file,/f77_unformatted,swap_endian=swap
    readu,1,ncpu_run
    readu,1,ndim
    readu,1,npartp
    if(npartp gt 0)then begin
    xx=fltarr(npartp)
    i1=iskip & i2=iskip+npartp-1L
    readu,1,xx
    part.xp(i1:i2,0)=xx
    readu,1,xx
    part.xp(i1:i2,1)=xx
    readu,1,xx
    part.xp(i1:i2,2)=xx
    readu,1,xx
    part.vp(i1:i2,0)=xx
    readu,1,xx
    part.vp(i1:i2,1)=xx
    readu,1,xx
    part.vp(i1:i2,2)=xx
    readu,1,xx
    part.mp(i1:i2)=xx
    xx=0.
    id=lonarr(npartp)
    readu,1,id
    part.id(i1:i2)=id
    id=0L
    endif
    close,1
    iskip=iskip+npartp
endfor
close,1

dummy=1d-30
print,'Writing file ',trim(tipsy)
openw,2,tipsy
printf,2,npart,ndum,ndum
printf,2,ndim
printf,2,time
for i=0L,npart-1L do printf,2,part.mp(i)
for i=0L,npart-1L do printf,2,part.xp(i,0)
for i=0L,npart-1L do printf,2,part.xp(i,1)
for i=0L,npart-1L do printf,2,part.xp(i,2)
for i=0L,npart-1L do printf,2,part.vp(i,0)
for i=0L,npart-1L do printf,2,part.vp(i,1)
for i=0L,npart-1L do printf,2,part.vp(i,2)
for i=0L,npart-1L do printf,2,dummy
for i=0L,npart-1L do printf,2,dummy
close,2

end


