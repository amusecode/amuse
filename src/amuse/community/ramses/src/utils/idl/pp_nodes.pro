
pro pp_nodes,tree,xr=xr,yr=yr,zr=zr,xproj=xproj,yproj=yproj,zproj=zproj,verbose=verbose,color=color

if not keyword_set(xr) then xr=[min(tree.x),max(tree.x)]
if not keyword_set(yr) then yr=[min(tree.y),max(tree.y)]
if not keyword_set(zr) then zr=[min(tree.z),max(tree.z)]
if not keyword_set(color) then color=255

xp=0
yp=0
zp=0
if not keyword_set(xproj) and  not keyword_set(yproj) and  not keyword_set(zproj) then zp=1
if keyword_set(xproj) then xp=1
if keyword_set(yproj) then yp=1
if keyword_set(zproj) then zp=1

t=FINDGEN(49)/49.*2.*!PI
a=cos(t)
b=sin(t)

jj=0L
for i=0,tree.n-1 do begin
    if (tree.x(i) ge xr(0) and tree.x(i) le xr(1) and $
        tree.y(i) ge yr(0) and tree.y(i) le yr(1) and $
        tree.z(i) ge zr(0) and tree.z(i) le zr(1) and $
        tree.idchild(i) ge 0) then begin
;        radius=tree.radius(i)
        radius=(tree.truemass(i)/(4.*!PI/3.*200.))^(1./3.)
        jj=jj+1L
        if keyword_set(verbose) then begin
            print,'Node #',i
            print,tree.x(i),tree.y(i),tree.z(i)
            print,tree.truemass(i),tree.mass(i)
        endif
        if zp eq 1 then begin
            xdum=tree.x(i)+a*radius
            ydum=tree.y(i)+b*radius
            zdum=tree.z(i)
        endif else if yp eq 1 then begin
            xdum=tree.x(i)+a*radius
            ydum=tree.y(i)
            zdum=tree.z(i)+b*radius
        endif else if xp eq 1 then begin
            xdum=tree.x(i)
            ydum=tree.y(i)+a*radius
            zdum=tree.z(i)+b*radius
        endif
        plots,xdum,ydum,zdum,color=color
    endif
endfor
print,'Number of active nodes=',jj
end

