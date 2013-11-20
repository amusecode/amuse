function mystring,n

n=long(n)
if n lt 0L or n gt 99999L then begin
    print,'input value too large...'
    return,' '
endif

if n le 9L then begin
    return,'0000'+string(n,format='(I1)')
endif
if n le 99L then begin
    return,'000'+string(n,format='(I2)')
endif
if n le 999L then begin
    return,'00'+string(n,format='(I3)')
endif
if n le 9999L then begin
    return,'0'+string(n,format='(I4)')
endif
if n le 99999L then begin
    return,string(n,format='(I5)')
endif


end


