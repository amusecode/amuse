;##################################################
;##               ROUTINE force.pro              ##
;## For problems, contact Romain.Teyssier@cea.fr ##
;##################################################
pro force,file=file,lbox=lbox,soft=soft,title=title,noerase=noerase,isolated=isolated,psym=psym

if not keyword_set(psym) then psym=3

if not keyword_set(lbox)then begin
    print,'You must specified the box size'
    return
endif

rd_part,p,file=file,/vel
xp=p.xp
vp=p.vp
ndimp=p.ndim
npart=p.npart

ind0=where(p.id eq 1)
x0=xp(ind0,*)

print,x0,format='("Centre=",3(1x,e10.3))'

if ndimp gt 0 then begin
    x=xp(0:npart-2,0)-x0(0)
    fx=vp(0:npart-2,0)
endif
if ndimp gt 1 then begin
    y=xp(0:npart-2,1)-x0(1)
    fy=vp(0:npart-2,1)
endif
if ndimp gt 2 then begin
    z=xp(0:npart-2,2)-x0(2)
    fz=vp(0:npart-2,2)
endif

case 1 of
    ndimp eq 1: begin
        r=x
        fr=ABS(fx)
        fe=2.d0*!DPI*(1.0d0-ABS(r/lbox*2.0d0))
        ft=ABS(fe-fr)
    end 
    ndimp eq 2: begin
        r=sqrt(x^2+y^2)
        fr=(fx*x+fy*y)/r
        f=sqrt(fx*fx+fy*fy)
        ft=sqrt(f^2-fr^2)
        fe=2.0d0/r
    end
    ndimp eq 3: begin
        r=sqrt(x^2+y^2+z^2+1d-15)

        if not keyword_set(isolated) then begin
                                ; Substract alias
            ewald,x/lbox,y/lbox,z/lbox,fex,fey,fez
            fx=fx-fex/lbox^2
            fy=fy-fey/lbox^2
            fz=fz-fez/lbox^2
        endif

        fr=(fx*x+fy*y+fz*z)/r
        f=sqrt(fx^2+fy^2+fz^2)
        ft=sqrt(f^2-fr^2)

        fe=1.d0/r/r

    end
    else: begin
        print,'ndim invalid'
        print,ndimp
        return
    end
endcase
 
if not keyword_set(title) then title=' '

;if not keyword_set(noerase) then begin
;    plot_oi,r,abs(fr)/fe-1. $
;      ,xr=[0.01,10.],yr=[-1,1] $
;      ,/xs,/ys,psym=3 $
;      ,xtitle='!17r',ytitle='!7D!17F/F' $
;      ,charsize=1. $
;      ,title=title  
;endif else begin
;    oplot,r,abs(fr)/fe-1.,psym=3
;endelse

if not keyword_set(noerase) then begin
    plot_oo,r,abs(fr)/fe $
      ,xr=[0.0001,10.],yr=[1.1d-5,9] $
      ,/xs,/ys,psym=psym $
      ,xtitle='!17r',ytitle='!17F/F!dtrue!n' $
      ,charsize=1. $
      ,title=title,/nodata
    oplot,r,abs(fr)/fe,psym=psym,color=180
endif else begin
    oplot,r,abs(fr)/fe,psym=psym,color=180
endelse
oplot,r,abs(ft)/fe,psym=psym,color=100

end

pro errp
!p.charsize=1.5
mycolor
x3d=[0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003]
e3d=[1.197e-2,1.189e-2,8.600e-3,6.960e-3,6.646e-3,6.583e-3,6.578e-3,6.576e-3]

x2d=[0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003]
e2d=[8.491e-3,5.600e-3,3.590e-3,3.628e-3,3.548e-3,3.523e-3,3.535e-3,3.534e-3]

x1d=[0.1     ,0.03    ,0.01    ,0.003   ,0.001   ,0.0003  ,0.0001  ,0.00003]
e1d=[7.902e-4,1.266e-3,1.130e-3,3.301e-4,1.312e-4,1.299e-5,3.215e-6,1.224e-6]

plot_oi,x3d,e3d*100,/xs,yr=[0.,1.5],psym=-6,xtitle='!7e!17',ytitle='Error (%)'
oplot,x2d,e2d*100,psym=-6,color=2
oplot,x1d,e1d*100,psym=-6,color=3
end
;###################################################
;###################################################
;###################################################
