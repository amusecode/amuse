pro rd_nodes,tree,leaf=leaf,verbose=verbose,filename=filename,fof=fof,isolated=isolated,multiple=multiple

if not keyword_set(filename) then begin
    file=DIALOG_PICKFILE(/READ)
endif else begin
    file=filename
endelse

openr,1,file
nnodes=0L
readf,1,nnodes,ndum
if keyword_set(verbose) then print,nnodes
level=lonarr(nnodes)
idmother=lonarr(nnodes)
idchild=lonarr(nnodes)
nsister=lonarr(nnodes)
idsister=lonarr(nnodes)
rhosad=fltarr(nnodes)
rhoave=fltarr(nnodes)
rhomax=fltarr(nnodes)
radius=fltarr(nnodes)
mass=fltarr(nnodes)
truemass=fltarr(nnodes)
x=fltarr(nnodes)
y=fltarr(nnodes)
z=fltarr(nnodes)

for i=0L,nnodes-1L do begin
    readf,1,idl,levell,idmotherl,idchildl,nsisterl,idsisterl,rhosadl,rhoavel,rhomaxl,radiusl,massl,truemassl,xl,yl,zl
    readf,1
    if keyword_set(verbose) then begin
        print,'===================================='
        print,idl,levell,idmotherl,idchildl,nsisterl,idsisterl
        print,rhosadl,rhoavel,rhomaxl,radiusl,massl
        print,truemassl,xl,yl,zl
        print
        read,itest
    endif
    level(i)=levell
    idmother(i)=idmotherl
    idchild(i)=idchildl
    nsister(i)=nsisterl
    idsister(i)=idsisterl
    rhosad(i)=rhosadl
    rhoave(i)=rhoavel
    rhomax(i)=rhomaxl
    radius(i)=radiusl
    mass(i)=massl
    truemass(i)=truemassl
    x(i)=xl
    y(i)=yl
    z(i)=zl
endfor

close,1

if keyword_set(leaf) then begin
; select leaf nodes ONLY
    new_idchild=lonarr(nnodes)
    new_radius=fltarr(nnodes)
    new_mass=fltarr(nnodes)
    new_truemass=fltarr(nnodes)
    new_x=fltarr(nnodes)
    new_y=fltarr(nnodes)
    new_z=fltarr(nnodes)
    j=0L
    for i=0L,nnodes-1L do begin
        if(idchild(i) eq 0)then begin
            new_idchild (j)=0L
            new_radius  (j)=radius(i)
            new_mass    (j)=mass(i)
            new_truemass(j)=truemass(i)
            new_x       (j)=x(i)
            new_y       (j)=y(i)
            new_z       (j)=z(i)
            j=j+1L
        endif
    endfor
    tree={n:j, idchild:new_idchild(0L:j-1L), radius:new_radius(0L:j-1L), mass:new_mass(0L:j-1L), truemass:new_truemass(0L:j-1L),x:new_x(0L:j-1L), y:new_y(0L:j-1L), z:new_z(0L:j-1L)}
endif else if keyword_set(fof) then begin
; select fof halos ONLY
    new_idchild=lonarr(nnodes)
    new_radius=fltarr(nnodes)
    new_mass=fltarr(nnodes)
    new_truemass=fltarr(nnodes)
    new_x=fltarr(nnodes)
    new_y=fltarr(nnodes)
    new_z=fltarr(nnodes)
    j=0L
    for i=0L,nnodes-1L do begin
        if(level(i) eq 1)then begin
            if(idchild(i) gt 0)then begin
                new_idchild(j)=nsister(idchild(i)-1)
            endif else begin
                new_idchild(j)=0L
            endelse
            if keyword_set(isolated) then begin
                test=( new_idchild(j) eq 0L )
            endif else if keyword_set(multiple) then begin
                test=( new_idchild(j) gt 0L )
            endif else begin
                test=1
            endelse
            if (test) then begin
                new_radius  (j)=radius(i)
                new_mass    (j)=mass(i)
                new_truemass(j)=truemass(i)
                new_x       (j)=x(i)
                new_y       (j)=y(i)
                new_z       (j)=z(i)
                j=j+1L
            endif
        endif
    endfor
    tree={n:j, idchild:new_idchild(0L:j-1L), radius:new_radius(0L:j-1L), mass:new_mass(0L:j-1L), truemass:new_truemass(0L:j-1L),x:new_x(0L:j-1L), y:new_y(0L:j-1L), z:new_z(0L:j-1L)}
endif else begin
    tree={n:nnodes, level:level, idmother:idmother, idchild:idchild, nsister:nsister, idsister:idsister, rhosad:rhosad, rhoave:rhoave, rhomax:rhomax, radius:radius, mass:mass, truemass:truemass,x:x, y:y, z:z}
    
endelse

;.r
ind=sort(tree.mass)
ind=reverse(ind)
for i=1L,tree.n-1L do begin
    xmax=tree.x(ind(i))
    ymax=tree.y(ind(i))
    zmax=tree.z(ind(i))
    x0=xmax
    y0=ymax
    z0=zmax
    r0=tree.radius(ind(i))
    x0=max([x0,2*r0])
    x0=min([x0,(1.0d0-2*r0)])
    y0=max([y0,2*r0])
    y0=min([y0,(1.0d0-2*r0)])
    z0=max([z0,2*r0])
    z0=min([z0,(1.0d0-2*r0)])
    xr=x0+[-2*r0,+2*r0]
    yr=y0+[-2*r0,+2*r0]
    zr=z0+[-2*r0,+2*r0]
;    ray3d,a,h,lmin=7,lmax=11,type=1,/log,/ave,/zproj,xr=xr,yr=yr,zr=zr
;    wait,1
    print,i,tree.mass(ind(i)),tree.x(ind(i)),tree.y(ind(i)),tree.z(ind(i)),$
     xr,yr,zr, $
     format='(I6,1x,I8,1x,"   [",e10.4,",",e10.4,",",e10.4,"]","   xr=[",e10.4,",",e10.4,"]",",yr=[",e10.4,",",e10.4,"]",",zr=[",e10.4,",",e10.4,"]")'
endfor

end

