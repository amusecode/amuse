pro ewald,x,y,z,fx,fy,fz,potential=potential
; Assumes L=1
print,'Computing alias using Ewald summation method...'

alpha = 4.0d0

if keyword_set(potential) then begin
    phi=double(!DPI/alpha^2+x-x)
endif else begin
    fx=double(0.0+x-x)
    fy=double(0.0+y-y)
    fz=double(0.0+z-z)
endelse

x0=double(x)
y0=double(y)
z0=double(z)

for i=-4,4 do begin
    for j=-4,4 do begin
        for k=-4,4 do begin
            x1=x0-double(i)
            y1=y0-double(j)
            z1=z0-double(k)
            r=sqrt(x1*x1+y1*y1+z1*z1+1d-15)
            erfc=double(1.d0-ERRORF(double(alpha*r)))
            fact=(erfc+2.d0*alpha/SQRT(!DPI)*r*double(EXP(-alpha^2*r^2)))
            if keyword_set(potential) then begin
                phi=phi-erfc/r
            endif else begin
                fx=fx-x1/r^3*fact
                fy=fy-y1/r^3*fact
                fz=fz-z1/r^3*fact
            endelse
        endfor
    endfor
endfor

u=!DPI/alpha^2+2.d0*alpha/sqrt(!DPI)
for i=-4,4 do begin
    for j=-4,4 do begin
        for k=-4,4 do begin
            h2=double(i*i+j*j+k*k)
            if (h2 ne 0.)then begin
                arg=double(2.d0*!DPI*(i*x0+j*y0+k*z0))
                if keyword_set(potential) then begin
                    fact=1.d0/!DPI/h2*EXP(-!DPI^2*h2/alpha^2)*cos(arg)
                    phi=phi-fact
                endif else begin
                    fact=2.d0/h2*EXP(-!DPI^2*h2/alpha^2)*sin(arg)
                    fx=fx-double(i)*fact
                    fy=fy-double(j)*fact
                    fz=fz-double(k)*fact
                endelse
                u=u-double(1.d0-ERRORF(alpha*sqrt(h2)))/sqrt(h2) $
                  -1.d0/!DPI/h2*EXP(-!DPI^2*h2/alpha^2)
            endif
        endfor
    endfor
endfor

r=sqrt(x0*x0+y0*y0+z0*z0+1d-15)
; Remove central particle contribution
if keyword_set(potential) then begin
    phi=phi+1.0d0/r
    fx=phi
endif else begin
    fx=fx+x0/r^3
    fy=fy+y0/r^3
    fz=fz+z0/r^3
endelse
print,u

end
       
