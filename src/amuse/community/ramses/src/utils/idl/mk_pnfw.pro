pro mk_pnfw,filename=filename,npart=npart,iseed=iseed


if not keyword_set(filename) then filename='part_nfw.ascii'
if not keyword_set(npart) then npart=32768L


rs=0.1
rtot=0.5
mtot=1.0

c=rtot/rs
ms=mtot/3./(alog(1.+c)-c/(1.+c))
print,ms


x0=0.5
y0=0.5
z0=0.5

distrib=RANDOMU(iseed,npart)
r=distrib*0.0

for i=0L,npart-1L do begin
    x_left=1d-6
    x_right=c
    err=1.0
    iter=1
    while(err gt 1d-6 and iter lt 100) do begin
        x_mil=0.5*(x_left+x_right)
        p_mil=3.*ms*(alog(1.0+x_mil)-x_mil/(1.0+x_mil))
        if (p_mil lt distrib(i)) then begin
            x_left=x_mil 
        endif else begin
            x_right=x_mil
        endelse
        err=ABS(x_left-x_right)
        iter=iter+1
    endwhile
    r(i)=x_mil*rs
endfor
costheta=2.d0*RANDOMU(iseed,npart)-1.d0
phi=RANDOMU(iseed,npart)*2.d0*!DPI
sintheta=sqrt(1.0d0-costheta^2)
cosphi=cos(phi)
sinphi=sin(phi)
x=0.5d0+r*cosphi*sintheta
y=0.5d0+r*sinphi*sintheta
z=0.5d0+r*costheta

openw,1,filename
for i=0L,npart-1L do begin
    xx=x(i)
    yy=y(i)
    zz=z(i)
    vx=0.
    vy=0.
    vz=0.
    mm=1.0
    printf,1,xx,yy,zz,vx,vy,vz,mm,format='(7(F17.14,1x))'
endfor
close,1








end
