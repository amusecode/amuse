pro pp_galic


readcol,'ic_part',x,y,z,u,v,w,m

massmin=min(m)
massmax=max(m)
massint=max(m(where(m ne massmax)))

print,massmin,massint,massmax


tek_color
indmin=where(m eq massmin)
indint=where(m eq massint)
indmax=where(m eq massmax)

xr=[max([min(x),min(y),min(z)]),min([max(x),max(y),max(z)])]
xr=[-20,20]

window,0,xs=512,ys=512
plot,x,y,/xs,/ys,/nodata,xr=xr,yr=xr,xtitle='!17x (kpc)',ytitle='!17y (kpc)'
oplot,x(indmax),y(indmax),color=4,psym=3
oplot,x(indint),y(indint),color=3,psym=3
oplot,x(indmin),y(indmin),color=2,psym=3

window,1,xs=512,ys=512
plot,x,z,/xs,/ys,/nodata,xr=xr,yr=xr,xtitle='!17x (kpc)',ytitle='!17z (kpc)'
oplot,x(indmax),z(indmax),color=4,psym=3
oplot,x(indint),z(indint),color=3,psym=3
oplot,x(indmin),z(indmin),color=2,psym=3

window,2,xs=512,ys=512
plot,y,z,/xs,/ys,/nodata,xr=xr,yr=xr,xtitle='!17y (kpc)',ytitle='!17z (kpc)'
oplot,y(indmax),z(indmax),color=4,psym=3
oplot,y(indint),z(indint),color=3,psym=3
oplot,y(indmin),z(indmin),color=2,psym=3


window,3,xs=512,ys=512
readcol,'Vcirc.dat',r,v
plot,r,v,xtitl='!17 r (kpc)',ytit='!17V!dcirv!n (km/s)',xr=[0,max(xr)]






end
