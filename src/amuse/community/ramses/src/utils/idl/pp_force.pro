function fana,r,rsoft
;return,r^3/(rsoft^3+r^3)
return,( (r/rsoft)^3 < 1. )
end

pro ppforce

mycolor
t0='!17PM 32!u3!n'
t1='!17AMR 32!u3!n'
erase
!p.multi=[0,2,2]
force,file='force_ndim3_nexpand2_lmax0_iseed123.dat',lbox=32,title=t0
for iseed=124,133 do begin
file='force_ndim3_nexpand2_lmax0_iseed'+strcompress(string(iseed,format='(i3)'))+'.dat'
force,file=file,lbox=32,/noerase
endfor
soft=1.
r=10.^(3.*(FINDGEN(100)/100.)-2.)
oplot,r,fana(r,soft),color=2,thick=2

!p.multi=[1,2,2]
force,file='force_ndim3_nexpand2_lmax6_iseed123.dat',lbox=32,title=t1+' + 6 levels'
for iseed=124,133 do begin
file='force_ndim3_nexpand2_lmax6_iseed'+strcompress(string(iseed,format='(i3)'))+'.dat'
force,file=file,lbox=32,/noerase
endfor
soft=1./64
r=10.^(3.*(FINDGEN(100)/100.)-2.)
oplot,r,fana(r,soft),color=2,thick=2

!p.multi=[2,2,2]
force,file='force_ndim3_nexpand2_lmax4_iseed123.dat',lbox=32,title=t1+' + 4 levels'
for iseed=124,133 do begin
file='force_ndim3_nexpand2_lmax4_iseed'+strcompress(string(iseed,format='(i3)'))+'.dat'
force,file=file,lbox=32,/noerase
endfor
soft=1./16
r=10.^(3.*(FINDGEN(100)/100.)-2.)
oplot,r,fana(r,soft),color=2,thick=2

!p.multi=[3,2,2]
force,file='force_ndim3_nexpand2_lmax2_iseed123.dat',lbox=32,title=t1+' + 2 levels'
for iseed=124,133 do begin
file='force_ndim3_nexpand2_lmax2_iseed'+strcompress(string(iseed,format='(i3)'))+'.dat'
force,file=file,lbox=32,/noerase
endfor
soft=1./4
r=10.^(3.*(FINDGEN(100)/100.)-2.)
oplot,r,fana(r,soft),color=2,thick=2


end
