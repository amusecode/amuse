pro plot_group

xz=[1,3]
yz=[1,3]
zz=[1,3]

!p.t3d=1
CREATE_VIEW $
  ,XMIN=xz[0],XMAX=xz[1] $
  ,YMIN=yz[0],YMAX=yz[1] $
  ,ZMIN=zz[0],ZMAX=zz[1] $
  ,AX=0,AY=0,AZ=0 $
  ,WINX=!d.x_size,WINY=!d.y_size

end

pro group,p,group,prefix=prefix,xr=xr,yr=yr,zr=zr

if not keyword_set(xr) then xr=[1,3]
if not keyword_set(yr) then yr=[1,3]
if not keyword_set(zr) then zr=[1,3]
if not keyword_set(prefix) then prefix='amas9_hires'

print,'Reading DENMAX file'
openr,1,'amas3_regroup.tag',/f77_unf
np=0L
ng=0L
readu,1,np,ng
print,np,ng
if not (np eq p.npart) then begin
    print,'file '+prefix+'.tag not compatible'
    return
endif
all_group=lonarr(np)
readu,1,all_group
close,1
  
print,minmax(all_group)

print,'Reading DENSITY file'
openr,1,prefix+'.den'
readf,1,ntot
if not (ntot eq p.npart) then begin
    print,'file '+prefix+'.den not compatible'
    return
endif
d=fltarr(ntot)
readf,1,d
close,1
dmin = min(d)
dmax = max(d)
color=15+BYTSCL(alog10(d),MIN=alog10(dmin),MAX=alog10(dmax) $
                ,TOP=!d.TABLE_SIZE-16)

ind=where(   p.xp(*,0) gt min(xr) and p.xp(*,0) lt max(xr) $
         and p.xp(*,1) gt min(yr) and p.xp(*,1) lt max(yr) $
         and p.xp(*,2) gt min(zr) and p.xp(*,2) lt max(zr) )
ind1=where(d(ind) eq max(d(ind)))
print,'Maximum density at ',p.xp(ind(ind1),0),p.xp(ind(ind1),1),p.xp(ind(ind1),2)
ind=where(   p.xp(*,0) gt min(xr) and p.xp(*,0) lt max(xr) $
         and p.xp(*,1) gt min(yr) and p.xp(*,1) lt max(yr) $
         and p.xp(*,2) gt min(zr) and p.xp(*,2) lt max(zr) $
         and d(*) gt 100.)

group=all_group(ind)
x=p.xp(ind,0)
y=p.xp(ind,1)
z=p.xp(ind,2)
u=p.vp(ind,0)
v=p.vp(ind,1)
w=p.vp(ind,2)
m=p.mp(ind)
d=d(ind)
col=color(ind)
mpart=min(m)
print,'particle mass=',mpart,max(m)

jmax=long(n_elements(ind)/100)
i=long(randomu(123,jmax)*n_elements(ind))
eps=0.01
for j=0L,jmax-1L do begin
    iter=0L
    x0=x(i(j))
    y0=y(i(j))
    z0=z(i(j))
    print,x0,y0,z0,i(j),j
    error=10.
    while (error gt 0.0000001) do begin
        r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)
        ind2=where(r lt eps,n_ok)
        x1=x0
        y1=y0
        z1=z0
        x0=total(x(ind2),/double)/double(n_ok)
        y0=total(y(ind2),/double)/double(n_ok)
        z0=total(z(ind2),/double)/double(n_ok)
        iter=iter+1
        error=sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2)
        print,iter,n_ok,x0,y0,z0,error
    endwhile

    r=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)
    ind2=where(r lt eps,n_ok)

    plots,x(ind2),y(ind2),z(ind2),color=2,/t3d,psym=3


endfor








nmax=LONG(max(group))
print,nmax

xg=DBLARR(nmax+1)
yg=DBLARR(nmax+1)
zg=DBLARR(nmax+1)
ug=DBLARR(nmax+1)
vg=DBLARR(nmax+1)
wg=DBLARR(nmax+1)
mg=DBLARR(nmax+1)

print,'Computing group properties'
print,'Please wait...'
t=FINDGEN(49)/49.*2.*!PI
a=cos(t)
b=sin(t)
usersym,a,b
test=0

for i=0,nmax do begin
    ind_group=where(group eq i, n_ok)
    if n_ok gt 0 then begin
        print,i+1
        imax=ind_group(where(d(ind_group) eq max(d(ind_group))))
        xg(i)=x(imax)
        yg(i)=y(imax)
        zg(i)=z(imax)
        print,xg(i),yg(i),zg(i)
        xg(i)=total(x(ind_group)*d(ind_group),/double)/total(d(ind_group))
        yg(i)=total(y(ind_group)*d(ind_group),/double)/total(d(ind_group))
        zg(i)=total(z(ind_group)*d(ind_group),/double)/total(d(ind_group))
        print,xg(i),yg(i),zg(i)
        ug(i)=total(u(ind_group))/double(n_ok)
        vg(i)=total(v(ind_group))/double(n_ok)
        wg(i)=total(w(ind_group))/double(n_ok)
        mg(i)=total(m(ind_group),/double)
        print,mg(i),mg(i)/mpart,n_ok
;         plots,x(ind_group),y(ind_group),z(ind_group),psym=3,color=4,/t3d
;         read,test
;         erase
;         c_loc=col(ind_group)
;         x_loc=x  (ind_group)
;         y_loc=y  (ind_group)
;         z_loc=z  (ind_group)
;         ind_col=sort(c_loc)
;         for j=0L,n_ok-1L do begin
;             plots,x_loc(ind_col(j)),y_loc(ind_col(j)),z_loc(ind_col(j)) $
;               ,/t3d,color=c_loc(ind_col(j)),psym=3
;         endfor
;         read,test
;         for j=0L,n_ok-1L do begin
;             plots,x_loc(ind_col(j)),y_loc(ind_col(j)),z_loc(ind_col(j)) $
;               ,/t3d,color=0,psym=3
;         endfor
    endif
endfor

ind_final=where(mg gt 0)
ng=n_elements(ind_final)
xg=xg(ind_final)
yg=yg(ind_final)
zg=zg(ind_final)
ug=ug(ind_final)
vg=vg(ind_final)
wg=wg(ind_final)
mg=mg(ind_final)

group={n:ng, x:xg, y:yg, z:zg, u:ug, v:vg, w:wg, m:mg}

return
end

pro get_group,p,list,g_number,prefix=prefix,xr=xr,yr=yr,zr=zr

if not keyword_set(xr) then xr=[1,3]
if not keyword_set(yr) then yr=[1,3]
if not keyword_set(zr) then zr=[1,3]
if not keyword_set(prefix) then prefix='amas9_hires'

print,'Reading DENMAX file'
openr,1,prefix+'.tag',/f77_unf
np=0L
ng=0L
readu,1,np,ng
if not (nread eq p.npart) then begin
    print,'file '+prefix+'.tag not compatible'
    return
endif
all_group=fltarr(np)
readf,1,all_group
close,1
  
print,'Reading DENSITY file'
openr,1,prefix+'.den'
readf,1,ntot
if not (ntot eq p.npart) then begin
    print,'file '+prefix+'.den not compatible'
    return
endif
d=fltarr(ntot)
readf,1,d
close,1

ind=where(   p.xp(*,0) gt min(xr) and p.xp(*,0) lt max(xr) $
         and p.xp(*,1) gt min(yr) and p.xp(*,1) lt max(yr) $
         and p.xp(*,2) gt min(zr) and p.xp(*,2) lt max(zr) )
ind1=where(d(ind) eq max(d(ind)))
print,'Maximum density at ',p.xp(ind(ind1),0),p.xp(ind(ind1),1),p.xp(ind(ind1),2)
ind=where(   p.xp(*,0) gt min(xr) and p.xp(*,0) lt max(xr) $
         and p.xp(*,1) gt min(yr) and p.xp(*,1) lt max(yr) $
         and p.xp(*,2) gt min(zr) and p.xp(*,2) lt max(zr) $
         and all_group gt 0.)

group=all_group(ind)
x=p.xp(ind,0)
y=p.xp(ind,1)
z=p.xp(ind,2)
u=p.vp(ind,0)
v=p.vp(ind,1)
w=p.vp(ind,2)
m=p.mp(ind)

print,'Extracting group',g_number
print,'Please wait...'
i=g_number
ind_group=where(group eq i, n_ok)

list={n:n_elements(ind_group), x:x(ind_group), y:y(ind_group), z:z(ind_group), u:u(ind_group), v:v(ind_group), w:w(ind_group), m:m(ind_group)}

return
end

pro p_group,g,mc=mc

if not keyword_set(mc) then mc=min(g.m)

t=FINDGEN(49)/49.*2.*!PI
a=cos(t)
b=sin(t)
usersym,a,b

for i=0,g.n-1 do if g.m(i) ge mc then plots,g.x(i),g.y(i),g.z(i),psym=8,symsize=2.*(g.m(i))^(1./3.),color=4,/t3d

end




