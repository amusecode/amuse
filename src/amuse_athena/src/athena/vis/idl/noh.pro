; makes 1-D scatter plot of d vs r for spherical Noh shock test (noh.c)
; plots EVERY grid point.  May create too large a plot for 3D grids (bigger
; that 100^3)
;
PRO noh_plot,filename
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi

readbin,filename
 
denr=fltarr(nx*ny*nz)
r=fltarr(nx*ny*nz)

for i=0,nx-1 DO BEGIN
  for j=0,ny-1 DO BEGIN
  for k=0,nz-1 DO BEGIN
    x2 = (float(i) + 0.5)*(float(i) + 0.5)
    y2 = (float(j) + 0.5)*(float(j) + 0.5)
    z2 = (float(k) + 0.5)*(float(k) + 0.5)
    r[i*ny*nz + j*nz + k] = sqrt(x2 + y2 + z2)
    denr[i*ny*nz + j*nz + k] = d[i,j,k]
  ENDFOR
  ENDFOR
ENDFOR
r=r/nx

plot,r,denr,XTITLE='R',YTITLE='D',psym=3

END
