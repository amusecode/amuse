;-------------------------------------------------------------
;+
; NAME:
;       REPCHR
; PURPOSE:
;       Replace all occurrences of one character with another in a text string.
; CATEGORY:
; CALLING SEQUENCE:
;       new = repchr(old, c1, [c2])
; INPUTS:
;       old = original text string.          in
;       c1 = character to replace.           in
;       c2 = character to replace it with.   in
;            default is space.
; KEYWORD PARAMETERS:
; OUTPUTS:
;       new = edited string.                 out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner.  28 Oct, 1986.
;       Johns Hopkins Applied Physics Lab.
;       RES 1 Sep, 1989 --- converted to SUN.
;       R. Sterner, 27 Jan, 1993 --- dropped reference to array.
;
; Copyright (C) 1986, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------
 
	FUNCTION REPCHR, OLD, C1, C2, help=hlp
 
	if (n_params(0) lt 2) or keyword_set(help) then begin
	  print,' Replace all occurrences of one character with another '+$
	    'in a text string.'
	  print,' new = repchr(old, c1, [c2])'
	  print,'   old = original text string.          in'
	  print,'   c1 = character to replace.           in'
	  print,'   c2 = character to replace it with.   in'
	  print,'        default is space.'
	  print,'   new = edited string.                 out'
	  return, -1
	endif
 
	B = BYTE(OLD)			   ; convert string to a byte array.
	CB1 = BYTE(C1)			   ; convert char 1 to byte.
	W = WHERE(B EQ CB1[0])		   ; find occurrences of char 1.
	IF W[0] EQ -1 THEN RETURN, OLD	   ; if none, return old string.
	IF N_PARAMS(0) LT 3 THEN C2 = ' '  ; default char 2 is space.
	CB2 = BYTE(C2)			   ; convert char 2 to byte.
	B[W] = CB2[0]			   ; replace char 1 by char 2.
	RETURN, STRING(B)		   ; return new string.
	END
