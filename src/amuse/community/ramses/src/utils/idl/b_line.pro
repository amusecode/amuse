;+       
;       
; NAME :  B_LINE       
; PURPOSE :       
;         Calculates a field line given a 3-dimensional B field       
; CALLING SEQUENCE :       
;         B_LINE, Bx, By, Bz, R0, R, Dx=Dx, Dy=Dy, Dz=Dz,Ds=Ds       
; INPUTS:       
;         Bx, By, Bz : 3-dimensional arrays of B components       
;         R0 : starting point (1-d array with 3 elements)        
; KEYWORD OPTIONAL INPUTS:       
;         Dx, Dy, Dz: increment values in x, y and z coordiantes       
;                  (Default=1)       
;                     
;         Ds : increment value of the arc length along the field       
;              lines( default=min(Dx, Dy, Dz))       
;              if Ds is negative, then the integration is done       
;              in the direction opposite to magnetic field.       
; OUPUTS :       
;         R : 2-d array. A series of spatial points defining       
;         the field line. 3 x N elements. N : number of       
;         elements.       
; RESTRICTION:       
;         R0 should be inside the spatial domain defined by       
;                          0 =< x =< (Nx-1)*Dx       
;                          0 =< y =< (Ny-1)*Dy       
;                          0 =< z =< (Nz-1)*Dz
;        and integration will be done while        
;        while the spatial points are inside that domain.       
; Modfication History          
;  June 1997, Jongchul Chae       
;-       
pro b_line, bx, by, bz,r0, r, dx=dx, dy=dy,dz=dz, ds=ds
       
nx = n_elements(bx(*,0,0))       
ny = n_elements(bx(0,*,0))       
nz = n_elements(bz(0,0,*))       
   
if n_elements(dx) eq 0 then dx=1.       
if n_elements(dy) eq 0 then dy=1.       
if n_elements(dz) eq 0 then dz=1.
if n_elements(ds) eq 0 then ds=dx<dy<dz
       
xmax=(nx-1)*dx       
ymax=(ny-1)*dy 
zmax=(nz-1)*dy 

r2=r0       
inside = r2(0)*(r2(0)-xmax) le 0. and  $       
         r2(1)*(r2(1)-ymax) le 0. and  $       
         r2(2)*(r2(2)-zmax) le 0.
      
if not inside then begin                
    tmp=widget_message('Starting point is outside the spatial domain.',/error)       
    return       
endif else begin       
    r=[transpose(r2)]   
    b2 = [interpolate(bx, r2(0)/dx, r2(1)/dy, r2(2)/dz), $       
          interpolate(by, r2(0)/dx, r2(1)/dy, r2(2)/dz), $       
          interpolate(bz, r2(0)/dx, r2(1)/dy, r2(2)/dz)]       
    b2=b2/norm(b2)       
endelse       
       
while inside do begin       
    r1=r2       
    b1=b2     
    iter2=0     
    repeat begin       
        r2_0 =r2       
        r2=r1+(b1+b2)*0.5*ds
        b2 = [interpolate(bx, r2(0)/dx, r2(1)/dy, r2(2)/dz), $       
              interpolate(by, r2(0)/dx, r2(1)/dy, r2(2)/dz), $       
              interpolate(bz, r2(0)/dx, r2(1)/dy, r2(2)/dz)]   
        b2=b2/norm(b2)             
        iter2=iter2+1       
    endrep until (norm(r2_0-r2) le ds*0.2) or iter2 ge 10       

    r=[r, transpose(r2)]       
    inside = r2(0)*(r2(0)-xmax) le 0. and  $       
             r2(1)*(r2(1)-ymax) le 0. and  $       
             r2(2)*(r2(2)-zmax) le 0. 
      
endwhile        

end       

  
                         
                       
       
