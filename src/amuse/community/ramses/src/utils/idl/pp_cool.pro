; function Teq,nH,cool=cool

; ix=alog10(nH/min(cool.n))
; ix=ix/(alog10(max(cool.n)/min(cool.n)))
; ix=ix*n_elements(cool.n)

; T2eq=INTERPOLATE(cool.teq,ix)   
; T2=1.+0.*T2eq
; iy=alog10(T2/min(cool.t))   
; iy=iy/(alog10(max(cool.t)/min(cool.t)))      
; iy=iy*n_elements(cool.t)   
; mu=10.^INTERPOLATE(alog10(cool.mu),ix,iy)   

; return,t2eq*mu
; end

function cmp_cool,cool,nH,T2

n1=cool.n1
n2=cool.n2
dlog_nH=double(n1-1)/(cool.n(n1-1)-cool.n(0))
dlog_T2=double(n2-1)/(cool.T(n2-1)-cool.T(0))
h=1d0/dlog_T2
h2=h*h
h3=h2*h

facH=alog10(nH)
ix=facH-min(cool.n)
ix=ix*dlog_nH
ixp1=ix+1

facT=alog10(T2)
iy=facT-min(cool.t)
iy=iy*dlog_T2
iyp1=iy+1

w1H=(cool.n(ixp1)-facH)*dlog_nH
w2H=(facH-cool.n(ix))*dlog_nH


yy=facT-cool.t(iy)
yy2=yy*yy
yy3=yy2*yy

fa=cool.cool(ix,iy)*w1H+cool.cool(ixp1,iy)*w2H
fb=cool.cool(ix,iyp1)*w1H+cool.cool(ixp1,iyp1)*w2H
fprimea=cool.cool_prime(ix,iy)*w1H+cool.cool_prime(ixp1,iy)*w2H
fprimeb=cool.cool_prime(ix,iyp1)*w1H+cool.cool_prime(ixp1,iyp1)*w2H
alpha=fprimea
beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
cmp_cool=10d0^(fa+alpha*yy+beta*yy2+gamma*yy3)

return,cmp_cool

end

function cmp_heat,cool,nH,T2

n1=cool.n1
n2=cool.n2
dlog_nH=double(n1-1)/(cool.n(n1-1)-cool.n(0))
dlog_T2=double(n2-1)/(cool.T(n2-1)-cool.T(0))
h=1d0/dlog_T2
h2=h*h
h3=h2*h

facH=alog10(nH)
ix=facH-min(cool.n)
ix=ix*dlog_nH
ixp1=ix+1

facT=alog10(T2)
iy=facT-min(cool.t)
iy=iy*dlog_T2
iyp1=iy+1

w1H=(cool.n(ixp1)-facH)*dlog_nH
w2H=(facH-cool.n(ix))*dlog_nH


yy=facT-cool.t(iy)
yy2=yy*yy
yy3=yy2*yy

fa=cool.heat(ix,iy)*w1H+cool.heat(ixp1,iy)*w2H
fb=cool.heat(ix,iyp1)*w1H+cool.heat(ixp1,iyp1)*w2H
fprimea=cool.heat_prime(ix,iy)*w1H+cool.heat_prime(ixp1,iy)*w2H
fprimeb=cool.heat_prime(ix,iyp1)*w1H+cool.heat_prime(ixp1,iyp1)*w2H
alpha=fprimea
beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
cmp_heat=10d0^(fa+alpha*yy+beta*yy2+gamma*yy3)

return,cmp_heat

end

function cmp_metal,cool,nH,T2

n1=cool.n1
n2=cool.n2
dlog_nH=double(n1-1)/(cool.n(n1-1)-cool.n(0))
dlog_T2=double(n2-1)/(cool.T(n2-1)-cool.T(0))
h=1d0/dlog_T2
h2=h*h
h3=h2*h

facH=alog10(nH)
ix=facH-min(cool.n)
ix=ix*dlog_nH
ixp1=ix+1

facT=alog10(T2)
iy=facT-min(cool.t)
iy=iy*dlog_T2
iyp1=iy+1

w1H=(cool.n(ixp1)-facH)*dlog_nH
w2H=(facH-cool.n(ix))*dlog_nH


yy=facT-cool.t(iy)
yy2=yy*yy
yy3=yy2*yy

fa=cool.metal(ix,iy)*w1H+cool.metal(ixp1,iy)*w2H
fb=cool.metal(ix,iyp1)*w1H+cool.metal(ixp1,iyp1)*w2H
fprimea=cool.metal_prime(ix,iy)*w1H+cool.metal_prime(ixp1,iy)*w2H
fprimeb=cool.metal_prime(ix,iyp1)*w1H+cool.metal_prime(ixp1,iyp1)*w2H
alpha=fprimea
beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
cmp_metal=10d0^(fa+alpha*yy+beta*yy2+gamma*yy3)

return,cmp_metal

end



pro pp_cool,cin,z=z

if not keyword_set(z) then z=0.0

n1=1000
n2=1000

nH=10d0^(FINDGEN(n1)/double(n1)*14.-8.)
T2=10d0^(FINDGEN(n2)/double(n2)*10.-2.)
net=dblarr(n1,n2)

for j=0,n1-1 do begin
    nn=nH(j)
    cool=cmp_cool(cin,nn,T2)
    heat=cmp_heat(cin,nn,T2)
    metal=cmp_metal(cin,nn,T2)
    
    net(j,*)=ABS(heat-cool-z*metal)
endfor
hhh=alog10(net)       
mycontour,hhh,alog10(nH),alog10(t2),ncont=50,max=-20,min=-28.,xtitle='!17log n!dH!n (cm!u-3!n)',ytitle='!17log T (K)',title='!17Net cooling rate (erg cm!u3!n)'

end

