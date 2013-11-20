;+
; NAME:
;	PP_HISTO
;
; PURPOSE:
;	This procedure reads an 2D histogram and makes the plot.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	PP_HISTO, FILE=file, SWAP=swap, VERBOSE=verbose 
;
; OPTIONAL INPUTS:
;	FILE:   if set, input the scalar string containing the name of
;	        the file to be read. Otherwise, a PICKFILE widget is
;	        launched.  
;
;       SWAP:   if set, reverse the bit ordering (Little Endian versus
;               Big Endian)
;
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;
;	        PP_HISTO, file='toto.dat',/swap
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2003.
;                       e-mail: Romain.Teyssier@cea.fr
;-
pro pp_histo, file=file, swap=swap, verbose=verbose

if not keyword_set(file) then begin
    file=DIALOG_PICKFILE(/READ)
endif
if not keyword_set(file) then return

if keyword_set(verbose) then print,'Reading file ',trim(file)
nx=0L & ny=0L
dmin=0d0 & dmax=0d0
tmin=0d0 & tmax=0d0

openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,nx,ny
image=fltarr(nx,ny)
readu,1,image
readu,1,dmin,dmax
readu,1,tmin,tmax
close,1

print,dmin,dmax
print,tmin,tmax

tt=tmin+FINDGEN(ny)/ny*(tmax-tmin)
dd=dmin+FINDGEN(nx)/nx*(dmax-dmin)

mycontour,image,dd,tt,ncont=300,/log,min=1d-10,/tab
return

bad_luck:  print,'I/O Error, exiting...'
           close,1
           return



end
;###################################################
;###################################################
;###################################################
