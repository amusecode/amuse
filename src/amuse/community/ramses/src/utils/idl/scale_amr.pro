pro scale_amr,a,h,scale_factor

print,'Rescaling AMR, HYDRO and PART variables to new box size'
a.boxlen=a.boxlen*scale_factor
print,'New box size (h-1 Mpc)=',a.boxlen
print,'Scaling hydro variables '
(*h.levelh[0]).u(*,*,*,1)=(*h.levelh[0]).u(*,*,*,1)*scale_factor
(*h.levelh[0]).u(*,*,*,2)=(*h.levelh[0]).u(*,*,*,2)*scale_factor
(*h.levelh[0]).u(*,*,*,3)=(*h.levelh[0]).u(*,*,*,3)*scale_factor
(*h.levelh[0]).u(*,*,*,4)=(*h.levelh[0]).u(*,*,*,4)*scale_factor^2
(*h.levelh[0]).u(*,*,*,5)=(*h.levelh[0]).u(*,*,*,5)*scale_factor^2
for i=1,a.nlevelmax-1 do begin
    (*h.levelh[i]).u(*,*,1)=(*h.levelh[i]).u(*,*,1)*scale_factor
    (*h.levelh[i]).u(*,*,2)=(*h.levelh[i]).u(*,*,2)*scale_factor
    (*h.levelh[i]).u(*,*,3)=(*h.levelh[i]).u(*,*,3)*scale_factor
    (*h.levelh[i]).u(*,*,4)=(*h.levelh[i]).u(*,*,4)*scale_factor^2
    (*h.levelh[i]).u(*,*,5)=(*h.levelh[i]).u(*,*,5)*scale_factor^2
endfor

end
