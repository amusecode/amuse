pro mk_grafic,vect,boxlen=boxlen,dir=dir,offset=offset

if not keyword_set(dir) then dir='test'

vect=float(vect)

if ((size(vect))(0) eq 3) then begin
    print,'added dimension for mk_grafic'
    vec=fltarr(1,(size(vect))(1),(size(vect))(2),(size(vect))(3))
    vec(0,*,*,*)=vect
    vect=vec
endif


nvar=(size(vect))(1)

if not keyword_set(boxlen) then boxlen=1.0
if not keyword_set(offset) then offset=0.
n1=long((size(vect))(2))
n2=long((size(vect))(3))
n3=long((size(vect))(4))
print,'boxlen=',boxlen
dxini=float(boxlen)/float(n1)
print,'dxini',dxini
xoff1=offset
xoff2=offset
xoff3=offset
astart=1.0
omega_m=0.0
omega_l=0.0
h0=0.0

if (dir eq 'here') then cd,current=dir
if (file_test(dir)eq 0) then file_mkdir,dir


fileout=['ic_d','ic_u','ic_v','ic_w','ic_p','ic_bxleft','ic_byleft','ic_bzleft','ic_bxright','ic_byright','ic_bzright']

for ivar=0,nvar-1 do begin    
    openw,1,dir+'/'+fileout(ivar),/f77_unf
    writeu,1,n1,n2,n3,dxini,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
    for i3=0,n3-1 do begin
        init_plane=reform(vect(ivar,*,*,i3),n1,n2)
        writeu,1,init_plane
    endfor
    close,1
endfor

end
