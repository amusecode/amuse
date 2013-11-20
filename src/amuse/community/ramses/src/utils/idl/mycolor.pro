pro mycolor

if !d.name eq 'PS' then begin
;              0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
    red   = [  0,  0,255,  0,  0,255,140,255,140,200,255,200,  0,150,255]
    green = [  0,  0,  0,  0,255,255,  0,140,140,  0,  0,255,140,150,255]
    blue  = [  0,  0,  0,255,  0,  0,140,  0,255,150,255,120,  0,150,255]
;=> 0=black, 1=black, 2=red, 3=blue, 4=green, 5=yellow, 6=purple, 7=brown, 8=light blue
endif else begin
;              0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
    red   = [  0,255,255,  0,  0,255,140,255,140,200,255,200,  0,150,255]
    green = [  0,255,  0,  0,255,255,  0,140,140,  0,  0,255,140,150,255]
    blue  = [  0,255,  0,255,  0,  0,140,  0,255,150,255,120,  0,150,255]
;=> 0=black, 1=white, 2=red, 3=blue, 4=green, 5=yellow, 6=purple, 7=brown, 8=light blue
endelse

TVLCT, red, green, blue

end
