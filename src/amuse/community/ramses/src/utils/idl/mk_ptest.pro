pro mk_ptest,filename=filename,npart=npart,cartesian=cartesian,iseed=iseed


if not keyword_set(filename) then filename='part_test.ascii'
if not keyword_set(npart) then npart=32768L

x0=RANDOMU(iseed,1)
y0=RANDOMU(iseed,1)
z0=RANDOMU(iseed,1)
x0=x0(0)
y0=y0(0)
z0=z0(0)
if keyword_set(cartesian) then begin
    x=RANDOMU(iseed,npart)
    y=RANDOMU(iseed,npart)
    z=RANDOMU(iseed,npart)
endif else begin
    r=0.5*10.d0^(3.*RANDOMU(iseed,npart)-3.)
    costheta=2.d0*RANDOMU(iseed,npart)-1.d0
    phi=RANDOMU(iseed,npart)*2.d0*!DPI
    sintheta=sqrt(1.0d0-costheta^2)
    cosphi=cos(phi)
    sinphi=sin(phi)
    x=0.5d0+r*cosphi*sintheta
    y=0.5d0+r*sinphi*sintheta
    z=0.5d0+r*costheta
    x(0)=0.5d0
    y(0)=0.5d0
    z(0)=0.5d0
endelse

openw,1,filename
for i=0L,npart-1L do begin
    xx=x(i)
    yy=y(i)
    zz=z(i)
    vx=0.
    vy=0.
    vz=0.
    mm=0.0
    if i eq 0 then mm=1.0
    printf,1,xx,yy,zz,vx,vy,vz,mm,format='(7(F17.14,1x))'
endfor
close,1








end
