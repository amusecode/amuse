;+
; NAME:
;	RD_CONE
;
; PURPOSE:
;	This procedure reads particles from a RAMSES PART file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_PART,Part, FILE=file, SWAP=swap,
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
pro rd_cone, part, file=file, swap=swap, density=density, velocity=velocity $
             ,ncpu=ncpu, time=time, verbose=verbose, nout=nout  $
             ,star=star, noalloc=noalloc

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_cone'
    RETURN
ENDIF

suffix=getcarnum(1024)

if not keyword_set(file) and not keyword_set(nout) then begin
    key='cone*.out*'
    file=DIALOG_PICKFILE(/READ,filter=key)
endif
if keyword_set(nout) then begin
    suffnout=getcarnum(nout)
    file='cone_'+suffnout(nout-1)+'/cone_'+suffnout(nout-1)+'#'
endif
if not keyword_set(file) then return
base_offset=strpos(file,'#')
file_base=strmid(file,0,base_offset)+'.out'

; Free memory associated to particles
part=0.

ncpu_run=0L & ndim=3L

ok=0
icpu=1
while (ok eq 0) do begin
   file=trim(file_base+suffix(icpu-1)+'.txt')
   if file_test(file) then begin
      print,'Reading file ',trim(file)
      openr,1,file
      readf,1,ncpu_run
      close,1
      ok=1
   endif
   icpu=icpu+1
endwhile

print,'ncpu      =',ncpu_run

if not keyword_set(ncpu) then ncpu=ncpu_run

suffix=getcarnum(ncpu)

npartp=0L & npart=0L & nstride=0L

dnpart=0d0
; Compute total number of particle
for jcpu=0,ncpu-1 do begin
    file=trim(file_base+suffix(jcpu))+'.txt'
    if keyword_set(verbose) then print,'Reading file ',trim(file)
    if file_test(file) then begin
       openr,1,file
       readf,1,ncpu_run
       readf,1,nstride
       readf,1,npartp
       close,1
       if keyword_set(verbose) then print,ndim,npartp $
               ,format='("ndim=",I1," npart=",I8)'
       npart=npart+npartp
       dnpart=dnpart+double(npartp)
    endif
endfor
print,'npart     =',npart

; Allocate memory
if not keyword_set(noalloc)then begin
    part={ ndim:ndim $
           ,npart:npart $
           ,xp:fltarr(npart,ndim) $
;       ,vp:fltarr(npart,ndim) $
    ,zp:fltarr(npart)}
endif  else begin
    part={ ndim:ndim $
           ,npart:npart}
endelse

iskip=0L
xmin=1d10
xmax=-1d10
rmin=1d10
rmax=-1d10
for jcpu=0,ncpu-1 do begin
    file=trim(file_base+suffix(jcpu))+'.txt'
    if file_test(file) then begin
       openr,1,file
       readf,1,ncpu_run
       readf,1,nstride
       readf,1,npartp
       close,1

       if npartp gt 0 then begin
          file=trim(file_base+suffix(jcpu))
          openr,1,file,/f77_unformatted,swap_endian=swap
          readu,1
          readu,1
          readu,1
          print,jcpu+1,npartp,nstride,sqrt(abs(rmin)),sqrt(abs(rmax))
          nblocs = npartp/nstride
          res = npartp - nblocs*nstride
          
;    print,nblocs,res
          for ibloc=0L,nblocs-1L do begin
             nread=nstride
             i1=iskip+ibloc*nstride
             i2=i1+nread-1L
;        print,iskip,ibloc,i1,i2,nread
             rr=fltarr(nread) 
             xx=fltarr(nread) 
             for idim=0L,ndim-1L do begin
                readu,1,xx
                xmax=max(xx,xmax)
                xmin=min(xx,xmin)
                rr=rr+xx^2
                if not keyword_set(noalloc)then begin
                   part.xp(i1:i2,idim)=xx
                endif
                readu,1,xx
;            part.vp(i1:i2,idim)=xx
             endfor        
             rmin=min([rmin,min(rr)])
             rmax=max([rmax,max(rr)])
             readu,1,xx
             if not keyword_set(noalloc)then begin
                part.zp(i1:i2)=xx
             endif
          endfor
          nread=res
          i1=iskip+nblocs*nstride
          i2=i1+nread-1
;    print,iskip,ibloc,i1,i2,nread
          xx=fltarr(nread)
          rr=fltarr(nread)
          for idim=0L,ndim-1L do begin
             readu,1,xx
             xmax=max(xx,xmax)
             xmin=min(xx,xmin)
             rr=rr+xx^2
             if not keyword_set(noalloc)then begin
                part.xp(i1:i2,idim)=xx
             endif
             readu,1,xx
;        part.vp(i1:i2,idim)=xx
          endfor
          rmin=min([rmin,min(rr)])
          rmax=max([rmax,max(rr)])
          readu,1,xx
          if not keyword_set(noalloc)then begin
             part.zp(i1:i2)=xx
          endif
          close,1
          iskip=iskip+npartp
       endif
    endif
endfor

print,sqrt(rmin),sqrt(rmax)
print,npart
end


