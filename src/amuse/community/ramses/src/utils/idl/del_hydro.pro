;+
; NAME:
;	DEL_HYDRO
;
; PURPOSE:
;	This procedure destroys the hydro structure in a clean way.
;
; CATEGORY:
;	Memory management.
;
; CALLING SEQUENCE:
;	DEL_HYDRO, Hydro
;
; OPTIONAL INPUTS:
;	None.
;	
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To destroy properly the HYDRO structure and free the allocated
;       memory, type:
;
;	        DEL_HYDRO, Hydro
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
pro del_hydro, hydro, verbose=verbose

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'del_hydro'
    RETURN
ENDIF

if n_tags(hydro) ne 0 then begin
    if keyword_set(verbose) then $
      print,'Freeing memory for hydro structure...'
    for icpu=0,hydro.ncpu-1 do begin
        for i=0,hydro.nlevelmax-1 do begin
            ptr_free,hydro.levelh[i,icpu]
        endfor
    endfor
    hydro=0
endif

end
;###################################################
;###################################################
;###################################################
