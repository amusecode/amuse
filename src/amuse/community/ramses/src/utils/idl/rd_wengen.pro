pro rd_wengen,d,u,v,w,p,dir=dir,boxlen=boxlen,nx=nx,ny=ny,nz=nz

if not keyword_set(dir) then dir='test1h-256'
if not keyword_set(boxlen) then boxlen=2500.
if not keyword_set(nx) then nx=256L
if not keyword_set(ny) then ny=nx
if not keyword_set(nz) then nz=nx

gamma=1.6666667

d=fltarr(nx,ny,nz)
openr,1,dir+'/'+dir+'.density'
for k=0L,nz-1L do begin
    for j=0L,ny-1L do begin
        for i=0L,nx-1L do begin
            readf,1,dd
            d(i,j,k)=dd
        endfor
    endfor
endfor
close,1
u=fltarr(nx,ny,nz)
openr,1,dir+'/'+dir+'.xvelocity'
for k=0L,nz-1L do begin
    for j=0L,ny-1L do begin
        for i=0L,nx-1L do begin
            readf,1,uu
            u(i,j,k)=uu
        endfor
    endfor
endfor
close,1
v=fltarr(nx,ny,nz)
openr,1,dir+'/'+dir+'.yvelocity'
for k=0L,nz-1L do begin
    for j=0L,ny-1L do begin
        for i=0L,nx-1L do begin
            readf,1,vv
            v(i,j,k)=vv
        endfor
    endfor
endfor
close,1
w=fltarr(nx,ny,nz)
openr,1,dir+'/'+dir+'.zvelocity'
for k=0L,nz-1L do begin
    for j=0L,ny-1L do begin
        for i=0L,nx-1L do begin
            readf,1,ww
            w(i,j,k)=ww
        endfor
    endfor
endfor
close,1
p=fltarr(nx,ny,nz)
openr,1,dir+'/'+dir+'.ethermal'
for k=0L,nz-1L do begin
    for j=0L,ny-1L do begin
        for i=0L,nx-1L do begin
            readf,1,pp
            p(i,j,k)=pp
        endfor
    endfor
endfor
close,1
u=u/boxlen
v=v/boxlen
w=w/boxlen
p=(gamma-1.0)*p*d/boxlen^2


end

