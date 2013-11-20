;+
; NAME:
;	RD_PART
;
; PURPOSE:
;	This procedure reads particles from a RAMSES PART file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_PART,Part, FILE=file, SWAP=swap, ICPU=icpu, NCPU=ncpu,
;	DENSITY=density, TIME=time
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
;       ICPU:    first cpu file to be read. Default: 1.
;
;       NCPU:    number of cpu files to read, starting from
;                icpu. Default: all files from icpu to ncpu_max.  
;
;       DENSITY: if set, read a file containing the SPH density for
;                each particle. Can also be set to the density file name.
;	
;       TIME:    if set, read in FILE the birth eopch of each
;                particle. Usefull for star formation runs.
;	
; OUTPUTS:
;	Part:   structure containing particles positions, velocities,
;	masses and levels. If DENSITY is set, it also contains SPH
;	densities. If TIME is set, it also contains the birth epoch.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To read on a SGI architecture a RAMSES PART file created on a
;       COMPAQ Workstation, type:
;
;	        RD_PART, part, file='part_00001.out',/swap
;
;       If the file was generated on the same IEEE system, just type:
;
;               RD_PART, part, file='part_00001.out'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
pro rd_part, part, file=file, swap=swap, density=density, velocity=velocity $
             ,ncpu=ncpu, icpu=icpu ,time=time, metal=metal, verbose=verbose, nout=nout

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_part'
    RETURN
ENDIF

if not keyword_set(icpu) then icpu=1

suffix=getcarnum(icpu)

if not keyword_set(file) and not keyword_set(nout) then begin
    key='*.out'+suffix(icpu-1)
    file=DIALOG_PICKFILE(/READ,filter=key)
endif
if keyword_set(nout) then begin
    suffnout=getcarnum(nout)
    file='output_'+suffnout(nout-1)+'/part_'+suffnout(nout-1)+'.out'
endif
if not keyword_set(file) then return
base_offset=strpos(file,'.out')+4
file_base=strmid(file,0,base_offset)

; Free memory associated to particles
part=0.

ncpu_run=0L & ndim=0L
nstar=0L & nsink=0L
mstar=0d0

print,'Reading file ',trim(file_base)
file=trim(file_base+suffix(icpu-1))
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ncpu_run
readu,1,ndim
readu,1
readu,1
readu,1,nstar
readu,1
readu,1
readu,1,nsink
close,1

print,'ncpu      =',ncpu_run
print,'ndim      =',ndim
if nsink gt 0 then print,'nsink     =',nsink
if nstar gt 0 then print,'nstar     =',nstar

if not keyword_set(ncpu) then ncpu=ncpu_run-icpu+1

suffix=getcarnum(icpu+ncpu-1)

npartp=0L & npart=0L

; Compute total number of particle
mstar_tot=0d0
for jcpu=0,ncpu-1 do begin
    file=trim(file_base+suffix(jcpu+icpu-1))
    if keyword_set(verbose) then print,'Reading file ',trim(file)
    openr,1,file,/f77_unformatted,swap_endian=swap
    readu,1,ncpu_run
    readu,1,ndim
    readu,1,npartp
    readu,1
    readu,1
    readu,1,mstar
    readu,1
    readu,1
    close,1
    if keyword_set(verbose) then print,ndim,npartp $
      ,format='("ndim=",I1," npart=",I8)'
    npart=npart+npartp
    mstar_tot=mstar_tot+mstar
endfor
if nstar gt 0 then print,'mstar     =',mstar_tot
print,'npart     =',npart

; Allocate memory
if not keyword_set(density) and not keyword_set(time) then begin
    if not keyword_set(velocity) then begin
        part={ ndim:ndim $
               ,npart:npart $
               ,xp:fltarr(npart,ndim) $
               ,id:lonarr(npart) $
               ,mp:fltarr(npart)}
    endif else begin
        part={ ndim:ndim $
               ,npart:npart $
               ,xp:fltarr(npart,ndim) $
               ,vp:fltarr(npart,ndim) $
               ,id:lonarr(npart) $
               ,mp:fltarr(npart)}
    endelse
endif else if not keyword_set(time) then begin
    if not keyword_set(velocity) then begin
    part={ ndim:ndim $
          ,npart:npart $
          ,xp:fltarr(npart,ndim) $
          ,id:lonarr(npart) $
          ,mp:fltarr(npart) $
          ,dp:fltarr(npart) }
    endif else begin
        part={ ndim:ndim $
               ,npart:npart $
               ,xp:fltarr(npart,ndim) $
               ,vp:fltarr(npart,ndim) $
               ,id:lonarr(npart) $
               ,mp:fltarr(npart) $
               ,dp:fltarr(npart) }
    endelse
endif else if not keyword_set(density) then begin
    if not keyword_set(velocity) then begin
        part={ ndim:ndim $
               ,npart:npart $
               ,xp:fltarr(npart,ndim) $
               ,id:lonarr(npart) $
               ,mp:fltarr(npart) $
               ,ap:fltarr(npart) }
    endif else begin
        part={ ndim:ndim $
               ,npart:npart $
               ,xp:fltarr(npart,ndim) $
               ,vp:fltarr(npart,ndim) $
               ,id:lonarr(npart) $
               ,mp:fltarr(npart) $
               ,ap:fltarr(npart) }
    endelse
endif else begin
    if not keyword_set(velocity) then begin
        part={ ndim:ndim $
               ,npart:npart $
               ,xp:fltarr(npart,ndim) $
               ,id:lonarr(npart) $
               ,mp:fltarr(npart) $
               ,dp:fltarr(npart) $
               ,ap:fltarr(npart) }
    endif else begin
        part={ ndim:ndim $
               ,npart:npart $
               ,xp:fltarr(npart,ndim) $
               ,vp:fltarr(npart,ndim) $
               ,id:lonarr(npart) $
               ,mp:fltarr(npart) $
               ,dp:fltarr(npart) $
               ,ap:fltarr(npart) }
    endelse
endelse

iskip=0L
for jcpu=0,ncpu-1 do begin
    file=trim(file_base+suffix(jcpu+icpu-1))
    openr,1,file,/f77_unformatted,swap_endian=swap
    readu,1,ncpu_run
    readu,1,ndim
    readu,1,npartp
    readu,1
    readu,1
    readu,1
    readu,1
    readu,1

    if(npartp gt 0)then begin
        xx=dblarr(npartp)
        i1=iskip & i2=iskip+npartp-1L
        readu,1,xx
        part.xp(i1:i2,0)=xx
        if (ndim gt 1)then begin
            readu,1,xx
            part.xp(i1:i2,1)=xx
        endif
        if (ndim gt 2)then begin
            readu,1,xx
            part.xp(i1:i2,2)=xx
        endif
        readu,1,xx
        if keyword_set(velocity) then part.vp(i1:i2,0)=xx
        if (ndim gt 1)then begin
            readu,1,xx
            if keyword_set(velocity) then part.vp(i1:i2,1)=xx
        endif
        if (ndim gt 2)then begin
            readu,1,xx
            if keyword_set(velocity) then part.vp(i1:i2,2)=xx
        endif
        readu,1,xx
        part.mp(i1:i2)=xx
        id=lonarr(npartp)
        readu,1,id
        part.id(i1:i2)=id
        readu,1,id
        if(nstar gt 0)then begin
            if keyword_set(time) then begin 
                readu,1,xx
                part.ap(i1:i2)=xx
            endif
            if keyword_set(metal) then begin 
                readu,1,xx
                part.zp(i1:i2)=xx
            endif
        endif
        xx=0d0
        id=0L
    endif
    close,1
    iskip=iskip+npartp
endfor

if keyword_set(density) then begin    
    if density eq 1 then begin
        file=DIALOG_PICKFILE(/READ,filter='*.dis')        
    endif else begin
        file=density
    endelse
    openr,1,file,swap_endian=swap,/f77_unf
;    openr,1,file
    nflag=0L & ntot=0L & nparmi=0L & nparma=0L & nparbuffer=0L
    readu,1,nflag,ntot,nparmi,nparma,nparbuffer
;    readf,1,ntot
    if not (ntot eq npart) then begin
        print,'file '+file+' not compatible'
        print,'ntot=',ntot
        close,1
        return
    endif
    print,ntot

    nblocs=ntot/nparbuffer
    res = npart - nblocs*nparbuffer
    print,nblocs,nparbuffer
    for ibloc=0,nblocs-1 do begin
        nread=nparbuffer
        i1=ibloc*nparbuffer
        i2=i1+nread-1
        print,ibloc,i1,i2,nread
        xx=dblarr(nread)
        readu,1,xx
        part.mp(i1:i2)=part.mp(i1:i2)/xx^3
    endfor
    nread=res
    i1=nblocs*nparbuffer
    i2=i1+nread-1
    xx=dblarr(nread)
    readu,1,xx
    part.mp(i1:i2)=part.mp(i1:i2)/xx^3
    close,1

endif

end


