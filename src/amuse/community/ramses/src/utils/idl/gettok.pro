function gettok,st,char
;+
; NAME:
;	GETTOK                                    
; PURPOSE:
;	Retrieve the first part of the string up to a specified character
; EXPLANATION:
;	GET TOKen - Retrieve first part of string until the character char 
;	is encountered.
;
; CALLING SEQUENCE:
;	token = gettok( st, char )
;
; INPUT:
;	char - character separating tokens, scalar string
;
; INPUT-OUTPUT:
;	st - (scalar) string to get token from (on output token is removed)
;
; OUTPUT:
;	token - scalar string value is returned 
;
; EXAMPLE:
;	If ST is 'abc=999' then gettok(ST,'=') would return
;	'abc' and ST would be left as '999' 
;
; HISTORY
;	version 1  by D. Lindler APR,86
;	Remove leading blanks    W. Landsman (from JKF)    Aug. 1991
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------
  On_error,2                           ;Return to caller

; if char is a blank treat tabs as blanks

  tab = string(9b)
  while strpos(st,tab) GE 0 do begin    ;Search for tabs
	pos = strpos(st,tab)
	strput,st,' ',pos
  endwhile

  st = strtrim(st,1)              ;Remove leading blanks

; find character in string

  pos = strpos(st,char)
  if pos EQ -1 then begin         ;char not found?
	token = st
 	st = ''
 	return, token
  endif

; extract token

 token = strmid(st,0,pos)
 len = strlen(st)
 if pos EQ (len-1) then st = '' else st = strmid(st,pos+1,len-pos-1)

;  Return the result.

 return,token
 end
