pro mov,nout=nout,iout=iout,type=type,filename=filename,xr=xr,yr=yr,vmin=vmin,vmax=vmax,lmin=lmin,lmax=lmax,log=log,showgrid=showgrid,dummy=dummy,xs=xs,ys=ys

if not keyword_set(iout) then iout=1
if not keyword_set(nout) then nout=1
if not keyword_set(filename) then filename='movie_default.gif'
if not keyword_set(xs) then xs=256
if not keyword_set(ys) then ys=256

window,0,xs=xs,ys=ys

; rd_amr,a,nout=iout
; rd_hydro,h,nout=iout
; tv2d,a,h,type=type,vmin=vmin,vmax=vmax,lmin=lmin,lmax=lmax,xr=xr,yr=yr,log=log,showgrid=showgrid,dummy=dummy 
mkall,iout

scr=tvrd(true=3)
im=color_quan(scr,3,r,g,b,cube=6)
write_gif,filename,im,r,g,b

for i=iout+1,nout do begin

;     rd_amr,a,nout=i
;     rd_hydro,h,nout=i
;     tv2d,a,h,type=type,vmin=vmin,vmax=vmax,lmin=lmin,lmax=lmax,xr=xr,yr=yr,log=log,showgrid=showgrid,dummy=dummy 
    mkall,i

    scr=tvrd(true=3)
    im=color_quan(scr,3,r,g,b,cube=6)
    write_gif,filename,im,r,g,b,/multiple
endfor

 write_gif,filename,/close

end


    
