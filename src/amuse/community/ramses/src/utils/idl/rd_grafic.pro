pro rd_grafic,d,h,file=file,verbose=verbose,header=header,swap_endian=swap_endian


n1=0L
n2=n1
n3=n1
dxini=0.0
xoff1=0.0
xoff2=xoff1
xoff3=xoff1
astart=0.0
omega_m=0.0
omega_l=0.0
h0=0.0

openr, 1,file,/f77_unf,swap_endian=swap_endian
readu,1,n1,n2,n3,dxini,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
if keyword_set(verbose) then print,n1,n2,n3,dxini
din_plane=fltarr(n1,n2)
d=fltarr(n1,n2,n3)
for i3=0,n3-1 do begin
    if keyword_set(verbose) then print,i3
    readu,1,din_plane
    d(*,*,i3)=din_plane
endfor
close,1

if keyword_set(header) then begin
    h={n1:n1, n2:n2, n3:n3, dx: dxini, xoff1: xoff1, xoff2: xoff2, xoff3: xoff3, astart:astart, omega_m:omega_m, omega_l:omega_l, h0:h0}
endif

end
