;+
; NAME:
;	RD_FRAME
;
; PURPOSE:
;	This procedure reads an 2D movie image
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_FRAME, Image, FILE=file, SWAP=swap, VERBOSE=verbose 
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
;	Image:   store the image in array Image
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;
;	        RD_FRAME, Image, file='toto.dat',/swap
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2003.
;                       e-mail: Romain.Teyssier@cea.fr
;-
pro rd_frame, image, file=file, swap=swap, verbose=verbose

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_frame'
    RETURN
ENDIF

if not keyword_set(file) then begin
    file=DIALOG_PICKFILE(/READ)
endif
if not keyword_set(file) then return

if keyword_set(verbose) then print,'Reading file ',trim(file)
nx=0L & ny=0L
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,t,dx,dy,dz
readu,1,nx,ny
image=fltarr(nx,ny)
readu,1,image
close,1

return

bad_luck:  print,'I/O Error, exiting...'
           close,1
           return

end
;###################################################
;###################################################
;###################################################
