pro im,nout,type=type

xr=[0,1]
;xr=[0.25,0.75]
yr=xr


rd_amr,a,nout=nout
rd_hydro,h,nout=nout

n=2L^a.nlevelmax
current=fltarr(n,n)

if keyword_set(type) then  begin
    tv2d,a,h,type=fix(type),save=current,xr=xr,yr=yr
endif else begin
    bx=fltarr(n,n)
    by=bx
    bxl=bx & bxr=bx
    byl=by & byr=by
    tv2d,a,h,type= 6,save=bxl,xr=xr,yr=yr
    tv2d,a,h,type= 7,save=byl,xr=xr,yr=yr
    tv2d,a,h,type= 9,save=bxr,xr=xr,yr=yr
    tv2d,a,h,type=10,save=byr,xr=xr,yr=yr
    bx=0.5*(bxl+bxr)
    by=0.5*(byl+byr)
    current=bx
    goto,norm
    for i=0L,n-1L do begin
        im1=max([0,i-1])
        ip1=min([n-1L,i+1])
        for j=0,n-1L do begin
            jm1=max([0,j-1])
            jp1=min([n-1L,j+1])
            current(i,j)=by(ip1,j)-by(im1,j)-(bx(i,jp1)-bx(i,jm1))
        endfor
    endfor
    goto,skip
    norm:
    current=bx^2+by^2
    skip:
endelse

if keyword_set(type) then begin
    mytv,bytscl(current)
endif else begin
    mytv,bytscl(current),/noadj;*1d3*double(n)/(xr(1)-xr(0)),min=-0.04,max=0.08),/noadj
endelse
print,min(current)
print,max(current)

; if keyword_set(type) then begin
;     im=bytscl(current,min=-3d-7,max=8d-7)
; endif else begin
;     im=bytscl(current,min=-3d-7,max=8d-7)
; endelse
; write_gif,'toto.gif',im,/multiple


del_amr,a
del_hydro,h



end
