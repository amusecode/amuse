pro visu_grid,a,h,cpu=cpu,type=type,clt=clt $
                   ,log=log $
                   ,xr=xr,yr=yr $
                   ,lmin=lmin,lmax=lmax $
                   ,save=save,byte=byte $
                   ,vmin=vmin,vmax=vmax,dummy=dummy $
                   ,verbose=verbose,showgrid=showgrid


make_hydro=0
if N_PARAMS() eq 2 THEN make_hydro=1

print,'Left button: pick coordinates and zoom in'
print,'Right button: pick coordinates and zoom out'
print,'Middle button: exit'

x0=0.5d0
y0=0.5d0
scale=1.0d0

start:

;print,scale
xr=x0+[-scale/2.,scale/2.]
yr=y0+[-scale/2.,scale/2.]
if make_hydro eq 1 then begin
    tv2d,a,h,xr=xr,yr=yr,type=type,clt=clt $
                   ,log=log $
                   ,lmin=lmin,lmax=lmax $
                   ,save=save,byte=byte $
                   ,vmin=vmin,vmax=vmax,dummy=dummy $
                   ,verbose=verbose,showgrid=showgrid
endif else begin
    pp_amr2d,a,x0=x0,y0=y0,scale=scale,/col,cpu=cpu
endelse

cursor,x,y,/dev
choice=!mouse.button
if choice eq 1 then begin
    xmin=x0+scale*(-0.5d0)
    ymin=y0+scale*(-0.5d0)*!d.y_size/!d.x_size
    x1=double(x)/(!d.x_size)
    y1=double(y)/(!d.x_size)
    x0=x1*scale+xmin
    y0=y1*scale+ymin
    scale=scale/2.0d0
    x0=max([x0,scale/2.])
    x0=min([x0,1.0d0-scale/2.])
    y0=max([y0,scale/2.])
    y0=min([y0,1.0d0-scale/2.])
    goto,start
endif else if choice eq 4 then begin
    x1=double(x)/(!d.x_size)
    y1=double(y)/(!d.y_size)
    x0=x0-0.5d0*scale+x1*scale
    y0=y0-0.5d0*scale+y1*scale
    scale=scale*2.0d0
    x0=max([x0,scale/2.])
    x0=min([x0,1.0d0-scale/2.])
    y0=max([y0,scale/2.])
    y0=min([y0,1.0d0-scale/2.])
    if scale ge 1.0 then begin
        scale=1.0
        x0=0.5
        y0=0.5
    endif
    goto,start
endif else begin
    goto,final
endelse
    
final:


end

