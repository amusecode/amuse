pro extract_part,file=file,swap=swap

if not keyword_set(file) then file=DIALOG_PICKFILE(/READ,filter='*part*')
print,'Opening file ',file
ndimp=0L
npart=0L

openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ndimp
readu,1,npart
print,ndimp,npart,format='("ndim=",I1," npart=",I8)'    

xp=dblarr(npart,ndimp)
readu,1,xp

xp(*,0)=xp(*,0)-32.+0.2
xp(*,1)=xp(*,1)-32.-0.1
xp(*,2)=xp(*,2)-32.-0.3

r2=xp(*,0)^2+xp(*,1)^2+xp(*,2)^2
ind=where(r2 lt 0.25,nok)
nok=long(nok)
print,'npart selected=',nok
if nok eq 0 then begin
    close,1
    return
endif

fileout=file+'.ext'
print,'Opening file ',fileout
openw,2,fileout,/f77_unformatted
writeu,2,ndimp
writeu,2,nok

x=dblarr(nok,3)
x(*,0)=xp(ind,0)
x(*,1)=xp(ind,1)
x(*,2)=xp(ind,2)
writeu,2,x
xp=0.
x=0.

vp=dblarr(npart,ndimp)
readu,1,vp
v=dblarr(nok,3)
v(*,0)=vp(ind,0)
v(*,1)=vp(ind,1)
v(*,2)=vp(ind,2)
writeu,2,v
vp=0.
v=0.

mp=dblarr(npart)
readu,1,mp
m=dblarr(nok)
m=mp(ind)
writeu,2,m
mp=0.
m=0.

lp=lonarr(npart)
readu,1,lp
close,1
l=dblarr(nok)
l=lp(ind)
writeu,2,l
close,1
close,2
lp=0.
l=0.

end
