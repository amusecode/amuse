;+
; NAME:
;	MK_TIPSY
;
; PURPOSE:
;	This procedure
;	create a TIPSY ASCII particles file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RAMSES2TIPSY, FILE=file, SWAP=swap, TIPSY=tipsyfile
;
; INPUTS:
;       x, y, z, vx, vy, vz, m
;
; OPTIONAL INPUTS:
;	FILE:    if set, input the scalar string containing the name of
;         	 the file to be read. Otherwise, a PICKFILE widget is
;         	 launched.  
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 19/09/2006.
;                       e-mail: Romain.Teyssier@cea.fr
;-
pro MK_TIPSY, xp, yp, zp, up, vp, wp, mp, file=file, time=time

if not keyword_set(time) then time=1.0d0

if not keyword_set(file) then file='part.ascii'

npart=n_elements(mp)
dummy=1d-30
ndum=0
ndim=3

print,'Writing file ',trim(file)
openw,2,file
printf,2,npart,ndum,ndum
printf,2,ndim
printf,2,time
for i=0L,npart-1L do printf,2,mp(i)
for i=0L,npart-1L do printf,2,xp(i)
for i=0L,npart-1L do printf,2,yp(i)
for i=0L,npart-1L do printf,2,zp(i)
for i=0L,npart-1L do printf,2,up(i)
for i=0L,npart-1L do printf,2,vp(i)
for i=0L,npart-1L do printf,2,wp(i)
for i=0L,npart-1L do printf,2,dummy
for i=0L,npart-1L do printf,2,dummy
close,2

end


