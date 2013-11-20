function getcarnum,nframes
ndigitmx = floor( alog10( 99999 ) ) + 1
car = strarr(nframes)
for i = 1L, nframes do begin
       a = string(i)
       if (i gt 0) then $
         ndigit =  floor( alog10( float(i) ) ) + 1 $
       else $
          ndigit = 1 	
       for j = 1L, ndigitmx - ndigit do begin
           a = '0' + a
       endfor         
       car(i-1) = strcompress(a, /remove_all)
;       print,car(i-1)
endfor
return,car
end
