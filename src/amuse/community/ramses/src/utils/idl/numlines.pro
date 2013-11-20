function numlines,file
;+
; NAME:
;     NUMLINES() 
; PURPOSE:
;     Return the number of lines in a file
;
; CALLING SEQUENCE:
;     nl = NUMLINES( filename )
; INPUT:
;     filename = name of file, scalar string
; OUTPUT:
;     nl = number of lines in the file, scalar longword
;          Set to -1 if the number of lines could not be determined
; METHOD:
;     If Unix then spawn to wc; otherwise read 1 line at a time and count
;
; PROCEDURE CALLS:
;     EXPAND_TILDE(), SPEC_DIR()
; MODIFICATION HISTORY:
;     W. Landsman                              February 1996
;     Use /bin/sh shell with wc under Unix     March 1997
;     Use EXPAND_TILDE() under Unix         September 1997
;     Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2

 if N_params() EQ 0 then begin
        print,'Syntax - nl = NUMLINES( file)
        return,-1
 endif

 nl = -1L
 openr,lun,file,/get_lun, ERROR = err
 if err NE 0 then begin
        if !VERSION.OS eq "vms" then file = spec_dir(file,'DAT') else $
        file = spec_dir(file)
        message,'ERROR - Unable to open file '+ file,/CON
        return,-1
 endif

 if !VERSION.OS_FAMILY EQ 'unix' then begin
         free_lun,lun
         if strpos(file,'~') GE 0 then file = expand_tilde(file)
         spawn,'wc -l < '+file, result, /sh    
         return,long(result[0])
 endif else begin                 ;=====>> Loop through file counting lines  
        On_ioerror,NOASCII
        nl = 0l
        tmp = ' '
         while not eof(lun) do begin
          readf,lun,tmp
          nl = nl + 1
          endwhile
         free_lun,lun
         return,nl
 endelse

NOASCII:
  message,'Error reading file ' + string(file),/CON
  return,-1
 end
