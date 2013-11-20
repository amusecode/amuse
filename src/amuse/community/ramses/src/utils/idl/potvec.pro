;###################################################
;###################################################
;###################################################
function potvec,grid,hydro

IF N_PARAMS() NE 2 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'potvec'
    RETURN,0
ENDIF

ncpu=grid.ncpu
ndim=grid.ndim
nlevelmax=grid.nlevelmax
ngrid=grid.ngrid

if ndim ne 2 then begin
    print,'Mesh should have 2 dimensions'
    print,'but ndim=',ndim
    return,0
endif

nx=2L^nlevelmax
ny=nx
print,nx,ny,format='("Computing image of size=",i4,"x",i4)'

bxl=FLTARR(nx,ny)
bxr=FLTARR(nx,ny)
byl=FLTARR(nx,ny)
byr=FLTARR(nx,ny)

for icpu=0,ncpu-1 do begin

    for ilevel=1,nlevelmax do begin

        dx=0.5d0^ilevel
        dxmin=0.5d0^nlevelmax

        if keyword_set(verbose) then $
          print,ilevel,ngrid(ilevel-1,icpu),format='("Level=",i2," ngrid=",i6)'

        if ngrid(ilevel-1,icpu) gt 0 then begin

            mesh2=(*grid.level[ilevel-1,icpu])

            x=mesh2.xg(*,0)
            y=mesh2.xg(*,1)

            for i=0,1 do begin
                for j=0,1 do begin

                    ind_cell=i+2*j
                    active=mesh2.son(*,ind_cell)
                    ind2=where(active eq 0)
                    
                    mesh=(*hydro.levelh[ilevel-1,icpu])

                    cxl=mesh.u(*,ind_cell,5)
                    cyl=mesh.u(*,ind_cell,6)
                    cxr=mesh.u(*,ind_cell,8)
                    cyr=mesh.u(*,ind_cell,9)
                    
                    if (ind2(0) ne -1) then begin

                        nxc=2L^(nlevelmax-ilevel)
                        for ii2=0L,n_elements(ind2)-1L do begin
                            
                            xl=x(ind2(ii2))+(double(i)-1.0d0)*dx
                            xr=x(ind2(ii2))+(double(i))*dx
                            yl=y(ind2(ii2))+(double(j)-1.0d0)*dx
                            yr=y(ind2(ii2))+(double(j))*dx
                            
                            il=fix(xl/dxmin)
                            jl=fix(yl/dxmin)

                            for ii=0L,nxc-1L do begin
                                for jj=0L,nxc-1L do begin
                                    xxl=double(ii)/double(nxc)
                                    xxr=(double(ii)+1.0d0)/double(nxc)
                                    yyl=double(jj)/double(nxc)
                                    yyr=(double(jj)+1.0d0)/double(nxc)

                                    bxl(il+ii,jl+jj)=cxr(ind2(ii2))*xxl+cxl(ind2(ii2))*(1.0-xxl)
                                    bxr(il+ii,jl+jj)=cxr(ind2(ii2))*xxr+cxl(ind2(ii2))*(1.0-xxr)
                                    byl(il+ii,jl+jj)=cyr(ind2(ii2))*yyl+cyl(ind2(ii2))*(1.0-yyl)
                                    byr(il+ii,jl+jj)=cyr(ind2(ii2))*yyr+cyl(ind2(ii2))*(1.0-yyr)
                                    
                                endfor
                            endfor
                        endfor
                        
                    endif
                    
                endfor
            endfor
            
        endif
    endfor
endfor


print,'Check divergence free...'
print,minmax((bxr-bxl)+(byr-byl))

az=fltarr(nx+1,ny+1)
az(0,0)=0.
for j=0,ny-1 do begin
    for i=1,nx do begin
        az(i,j)=az(i-1,j)-byl(i-1,j)
    endfor
    az(0,j+1)=az(0,j)+bxl(0,j)
endfor
for i=1,nx-1 do begin
    az(i,ny)=az(i,ny-1)+bxl(i,ny-1)
endfor
az(nx,ny)=az(nx,ny-1)+bxr(nx-1,ny-1)

return,az
end
;###################################################
;###################################################
;###################################################
