;+
; NAME:
;	DEL_AMR
;
; PURPOSE:
;	This procedure destroys the mesh structure in a clean way.
;
; CATEGORY:
;	Memory management.
;
; CALLING SEQUENCE:
;	DEL_AMR, Grid
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
;       To destroy properly the AMR structure and free the allocated
;       memory, type:
;
;	        DEL_AMR, Grid
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
pro del_amr, grid, verbose=verbose

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'del_amr'
    RETURN
ENDIF

; Free memory associated to grid
if n_tags(grid) ne 0 then begin
    if keyword_set(verbose) then $
      print,'Freeing memory for grid structure...'
    for icpu=0,grid.ncpu-1 do begin
        for i=0,grid.nlevelmax-1 do begin
            ptr_free,grid.level[i,icpu]
        endfor
    endfor
    grid=0
endif

end
;###################################################
;###################################################
;###################################################
