;
;+
; NAME:
;       STR_TOKEN
; PURPOSE:
;       Retrieve portion of string up to token.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = str_token( old, token )
; INPUTS:
;       old   -- String to be split.  Contains text after    in, out
;                token on output.
;       token -- Token to use in splitting old.              in
; KEYWORD PARAMETERS:
;       /TRIM -- set to remove leading blanks from old 
;                before returning.
;       /HELP -- print useful message and exit.
; OUTPUTS:
;       new   -- portion of string up to token.              out
;       old   -- portion of old after token.                 out, in
; COMMON BLOCKS:
; SIDE EFFECTS:
;       Input parameter old is modified.
; NOTES:
;       Token may be one or more characters.
;       If token is not found, returns old and sets old to ''.
; EXAMPLE:
;       If old is 'foo44 bar', then str_token( old, '44' ) would return
;       'foo', and upon return, old will be left with ' bar'.  If /TRIM
;       were set, old would be 'bar' on return.
;
;       If old='xyz', then new=str_token(old,'a') would return with
;       new='xyz' and old=''.
; THANKS:
;       To D. Linder who wrote GETTOK, part of the goddard library,
;       upon which this is based.
; MODIFICATION HISTORY:
;       $Id: str_token.pro,v 1.1 2000/06/14 19:09:22 mcraig Exp $
;       $Log: str_token.pro,v $
;       Revision 1.1  2000/06/14 19:09:22  mcraig
;       Changed name of strtok str_token to avoid conflict in IDL 5.3.
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Added built in help.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: Rel_2_1_2 $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
FUNCTION Str_token, string, token, $
                 TRIM=trim, HELP=Help

; Back to the caller if error occurs.
    On_error, 2

    IF (n_params() NE 2) OR keyword_set(Help) THEN BEGIN 
        offset = '   '
        print, offset+'Retrieve portion of string up to token.'
        print, offset+'new = str_token( old, token )'
        print, offset+'Inputs:'
        print, offset+offset+'old   -- String to be split.  Contains text after    in, out'
        print, offset+offset+'         token on output.'
        print, offset+offset+'token -- Token to use in splitting old.              in'
        print, offset+'Keywords:'
        print, offset+offset+'/TRIM -- set to remove leading blanks from old '
        print, offset+offset+'         before returning.'
        print, offset+offset+'/HELP -- print useful message and exit.'
        print, offset+'Outputs:'
        print, offset+offset+'new   -- portion of string up to token.              out'
        print, offset+offset+'old   -- portion of old after token.                 out, in'
        print, offset+'Side effects:'
        print, offset+offset+'Input parameter old is modified.'
        print, offset+'Notes:'
        print, offset+offset+'Token may be one or more characters.'
        print, offset+offset+"If token is not found, returns old and sets old to ''."
        print, offset+'Examples:'
        print, offset+offset+"If old is 'foo44 bar', then str_token( old, '44' ) would return'"
        print, offset+offset+"  'foo', and upon return, old will be left with ' bar'.  If /TRIM"
        print, offset+offset+"  were set, old would be 'bar' on return."
;'
        print, offset+offset+"If old='xyz', then new=str_token(old,'a') would return with"
        print, offset+offset+"  new='xyz' and old=''."
        return, -1
    ENDIF 

    pos = strpos(string, token)

    IF (pos GE 0) THEN BEGIN
        front = strmid(string, 0, pos) 
        string = strmid(string, pos + strlen(token), strlen(string))
        IF keyword_set(trim) THEN string = strtrim(string, 1)
        return, front
    ENDIF
    
    front = string
    string = ''
    return, front
    
END
