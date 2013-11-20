pro make_histo,a,h,histogram,n=n,xr=xr,yr=yr,partition=partition,cool=cool


if not keyword_set(partition) then partition=1L
if not keyword_set(n) then n=200L
if not keyword_set(xr) then xr=[1d-3,1d09]
if not keyword_set(yr) then yr=[1d+3,1d08]
partition=long(partition)

histo=dblarr(n+2,n+2)
aexp=a.aexp
dx=1./partition
dy=1./partition
dz=1./partition
xz=[0.,dx]
yz=[0.,dy]
zz=[0.,dz]

;xz=0.5+[-1.,1.]/16.
;yz=0.5+[-1.,1.]/16.
;zz=0.5+[-1.,1.]/16.

for ii=0L,partition-1L do begin
    for jj=0L,partition-1L do begin
        for kk=0L,partition-1L do begin
            amr2cell,a,h,c,xr=ii*dx+xz,yr=jj*dy+yz,zr=kk*dz+zz
            d=c.var(*,0)
            p=c.var(*,4)
            vol=(c.dx)^3
            p=p/d*(a.unit_l/a.unit_t)^2*1.66d-24/1.38d-16
            d=d/0.045*0.3

            print,'rho =',minmax(d)

            if keyword_set(cool) then begin
                nH=c.var(*,0)*0.3*1.88d-29*0.7^2/1.66d-24/aexp^3*0.76
                T2=c.var(*,4)/c.var(*,0)*(1.d5/aexp)^2*1.66d-24/1.38d-16
                
                ix=alog10(nH/min(cool.n))
                ix=ix/(alog10(max(cool.n)/min(cool.n)))
                ix=ix*n_elements(cool.n)
                
                T2eq=INTERPOLATE(cool.teq,ix)
                
                T2=T2/T2eq
                iy=alog10(T2/min(cool.t))
                iy=iy/(alog10(max(cool.t)/min(cool.t)))
                iy=iy*n_elements(cool.t)
                
                mu=10.^INTERPOLATE(alog10(cool.mu),ix,iy)
                p=p*mu

                print,'T   =',minmax(p)
            endif else begin
                print,'T/mu=',minmax(p)
            endelse

            p=ABS(p)

            i0=(alog10(d)-alog10(xr[0]))/(alog10(xr[1])-alog10(xr[0]))*n+1
            indok=where(i0 ge (n+1),nok)
            if nok gt 0 then i0(indok)=n+1
            indok=where(i0 le 0,nok)
            if nok gt 0 then i0(indok)=0
            j0=(alog10(p)-alog10(yr[0]))/(alog10(yr[1])-alog10(yr[0]))*n+1
            indok=where(j0 ge (n+1),nok)
            if nok gt 0 then j0(indok)=n+1
            indok=where(j0 le 0,nok)
            if nok gt 0 then j0(indok)=0
            for i=0L,n_elements(d)-1L do begin
                i1=i0(i)
                j1=j0(i)
                histo(i1,j1)=histo(i1,j1)+d(i)*vol(i)
            endfor
        endfor
    endfor
endfor
histo=histo(1:n,1:n)
print,'Mass fraction=',total(histo,/double)
x=xr[0]*10.^(FINDGEN(n)/double(n-1)*alog10(xr[1]/xr[0]))
y=yr[0]*10.^(FINDGEN(n)/double(n-1)*alog10(yr[1]/yr[0]))

histogram={h:histo,d:x,t:y}

end
